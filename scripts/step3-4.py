#!/usr/bin/env python3
"""
Step 3.4: Extract edges from raw datasets
Streams through all datasets and outputs edges to CSV file in format:
source_raw_id | relation | target_raw_id | dataset_id
"""

import os
import json
import pandas as pd
from pathlib import Path
from xml.etree.ElementTree import iterparse
import csv
from datetime import datetime

# ===================== CONFIG =====================
RAW_DATA_DIR = Path("data/raw")
RELATION_VOCAB_FILE = Path("data/processed/relationship_vocab.csv")
OUTPUT_CSV = Path(f"data/processed/extracted_edges_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv")

# ===================== LOAD RELATIONSHIP VOCAB =====================
relation_df = pd.read_csv(RELATION_VOCAB_FILE)
RELATION_SET = set()

for r in relation_df['relationship']:
    if pd.isna(r):
        continue
    r = str(r).strip()
    # Filter out junk
    if r in ('-', '=', '>'):
        continue
    if r.upper().startswith('INFORES:'):
        continue
    # Normalize
    RELATION_SET.add(r.upper().replace(' ', '_'))

print(f"Loaded {len(RELATION_SET)} valid relationships from {RELATION_VOCAB_FILE}")

# ===================== HELPER FUNCTIONS =====================

def normalize_predicate(predicate: str) -> str:
    """Return a cleaned uppercase predicate or RELATED_TO if unknown."""
    if not predicate or pd.isna(predicate):
        return "RELATED_TO"
    pred = str(predicate).upper().replace(' ', '_')
    if pred in RELATION_SET:
        return pred
    return "RELATED_TO"

def stream_json_objects(file_handle):
    """Yield JSON objects from file."""
    content = file_handle.read().strip()
    try:
        data = json.loads(content)
        if isinstance(data, list):
            yield from data
        elif isinstance(data, dict):
            yield data
        return
    except json.JSONDecodeError:
        file_handle.seek(0)
        for line in file_handle:
            line = line.strip()
            if line:
                try:
                    yield json.loads(line)
                except:
                    continue

# ===================== EXTRACTION FUNCTIONS =====================

def extract_from_csv_tsv(file_path: Path, dataset_id: str, csv_writer):
    """Extract edges from CSV/TSV files."""
    delimiter = '\t' if file_path.suffix == '.tsv' else ','

    for chunk in pd.read_csv(
        file_path,
        chunksize=50000,
        delimiter=delimiter,
        low_memory=False,
        on_bad_lines='skip'
    ):
        cols = {c.lower(): c for c in chunk.columns}

        subj_col = next(
            (cols[c] for c in ['subject', 'subject_id', 'entity_1', 'source', 'source_id'] if c in cols),
            None
        )
        obj_col = next(
            (cols[c] for c in ['object', 'object_id', 'entity_2', 'target', 'target_id'] if c in cols),
            None
        )
        pred_col = next(
            (cols[c] for c in ['predicate', 'relation', 'relationship', 'edge_label'] if c in cols),
            None
        )

        if not subj_col or not obj_col:
            continue

        for _, row in chunk.iterrows():
            src = str(row[subj_col]).strip()
            tgt = str(row[obj_col]).strip()

            if not src or not tgt or src.lower() == 'nan' or tgt.lower() == 'nan':
                continue

            raw_pred = str(row[pred_col]).strip() if pred_col else ""
            relation = normalize_predicate(raw_pred)

            csv_writer.writerow([src, relation, tgt, dataset_id])


def extract_from_json(file_path: Path, dataset_id: str, csv_writer):
    """Extract edges from JSON files."""
    with open(file_path, 'r', encoding='utf-8') as f:
        for obj in stream_json_objects(f):
            src = (
                obj.get('subject')
                or obj.get('subject_id')
                or obj.get('entity_1')
                or obj.get('source')
                or obj.get('source_id')
            )
            tgt = (
                obj.get('object')
                or obj.get('object_id')
                or obj.get('entity_2')
                or obj.get('target')
                or obj.get('target_id')
            )

            if not src or not tgt:
                continue

            src = str(src).strip()
            tgt = str(tgt).strip()

            raw_pred = (
                obj.get('predicate')
                or obj.get('relation')
                or obj.get('relationship')
                or obj.get('edge_label')
                or ""
            )

            relation = normalize_predicate(raw_pred)
            csv_writer.writerow([src, relation, tgt, dataset_id])


# You can also similarly adapt extract_from_xml and extract_from_fasta as before, using normalize_predicate()

# ===================== MAIN =====================

with open(OUTPUT_CSV, 'w', newline='', encoding='utf-8') as csvfile:
    writer = csv.writer(csvfile, delimiter='|')
    writer.writerow(['source_raw_id', 'relation', 'target_raw_id', 'dataset_id'])

    for dataset_dir in sorted(RAW_DATA_DIR.iterdir()):
        if not dataset_dir.is_dir():
            continue
        dataset_id = dataset_dir.name
        print(f"Processing dataset: {dataset_id}")
        for file in dataset_dir.iterdir():
            if not file.is_file():
                continue
            if file.suffix.lower() in ['.csv', '.tsv']:
                extract_from_csv_tsv(file, dataset_id, writer)
            elif file.suffix.lower() == '.json':
                extract_from_json(file, dataset_id, writer)
            else:
                continue  # Ignore unsupported files

print(f"✅ 3.4 COMPLETE — Edges saved to {OUTPUT_CSV}")
