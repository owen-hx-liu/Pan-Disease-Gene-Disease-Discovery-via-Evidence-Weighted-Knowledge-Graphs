#!/usr/bin/env python3
"""
Step 3.5 — Dataset audit & relationship vocabulary discovery

• Streams files safely
• Never fabricates relationships
• Handles CSV / TSV / JSON / XML / FASTA
• Skips Orphadata-style files gracefully
"""

import os
import json
import pandas as pd
import xml.etree.ElementTree as ET
import re
from Bio import SeqIO

# ================= USER SETTINGS =================

DATA_DIR = "data/raw"
SAMPLE_SIZE = 2000

RELATIONSHIP_HINTS = [
    "relation",
    "relationship",
    "predicate",
    "edge",
    "association",
    "interaction",
    "type"
]

IDENTIFIER_HINTS = [
    "subject", "object", "source", "target",
    "entity", "gene", "disease", "protein", "id"
]

SUPPORTED_EXTENSIONS = (".csv", ".tsv", ".json", ".xml", ".fasta", ".fa")

# =================================================


# ---------- Normalization ----------

def normalize_relationship(rel):
    if not rel:
        return None
    rel = str(rel).strip()
    if not rel:
        return None
    rel = re.sub(r"[^A-Za-z0-9]+", "_", rel)
    rel = re.sub(r"_+", "_", rel)
    return rel.upper().strip("_")


# ---------- Column inference ----------

def infer_relationship_column(columns):
    for c in columns:
        if not isinstance(c, str):
            continue
        c_low = c.lower()
        for hint in RELATIONSHIP_HINTS:
            if hint in c_low:
                return c
    return None


def has_identifier_columns(columns):
    for c in columns:
        if not isinstance(c, str):
            continue
        c_low = c.lower()
        for h in IDENTIFIER_HINTS:
            if h in c_low:
                return True
    return False


# ---------- Relationship extraction ----------

def extract_relationships_df(df, global_relationships):
    rel_col = infer_relationship_column(df.columns)
    if not rel_col:
        return 0

    count = 0
    for v in df[rel_col].dropna():
        norm = normalize_relationship(v)
        if norm:
            global_relationships.add(norm)
            count += 1
    return count


# ---------- CSV / TSV ----------

def sample_csv(filepath, global_relationships):
    sep = "\t" if filepath.endswith(".tsv") else ","
    try:
        reader = pd.read_csv(filepath, sep=sep, chunksize=SAMPLE_SIZE, low_memory=False)
        chunk = next(reader)
    except Exception as e:
        return False, f"CSV read failed: {e}"

    if not has_identifier_columns(chunk.columns):
        return False, "No identifier-like columns detected"

    rels = extract_relationships_df(chunk, global_relationships)
    return True, f"CSV OK — sampled {len(chunk)} rows — relationships found: {rels}"


# ---------- JSON ----------

def sample_json(filepath, global_relationships):
    try:
        with open(filepath, "r", encoding="utf-8") as f:
            data = json.load(f)
    except Exception as e:
        return False, f"JSON read failed: {e}"

    if not isinstance(data, list):
        return False, "JSON is not a list of records"

    df = pd.DataFrame(data[:SAMPLE_SIZE])

    if df.empty:
        return False, "JSON empty after sampling"

    if not has_identifier_columns(df.columns):
        return False, "No identifier-like columns detected"

    rels = extract_relationships_df(df, global_relationships)
    return True, f"JSON OK — sampled {len(df)} rows — relationships found: {rels}"


# ---------- XML ----------

def sample_xml(filepath, global_relationships):
    rows = []
    try:
        for _, elem in ET.iterparse(filepath, events=("end",)):
            row = {}
            for child in elem:
                if isinstance(child.tag, str):
                    row[child.tag.lower()] = child.text
            if row:
                rows.append(row)
            elem.clear()
            if len(rows) >= SAMPLE_SIZE:
                break
    except Exception as e:
        return False, f"XML parse failed: {e}"

    if not rows:
        return False, "No records extracted from XML"

    df = pd.DataFrame(rows)

    if not has_identifier_columns(df.columns):
        return False, "No identifier-like columns detected"

    rels = extract_relationships_df(df, global_relationships)
    return True, f"XML OK — streamed {len(df)} records — relationships found: {rels}"


# ---------- FASTA ----------

def sample_fasta(filepath, global_relationships):
    count = 0
    try:
        for _ in SeqIO.parse(filepath, "fasta"):
            count += 1
            if count >= SAMPLE_SIZE:
                break
    except Exception as e:
        return False, f"FASTA read failed: {e}"

    return True, f"FASTA OK — sampled {count} sequences (no relationships)"


# ---------- File discovery ----------

def discover_files(root_dir):
    collected = []
    for root, _, files in os.walk(root_dir):
        for f in files:
            if f.lower().endswith(SUPPORTED_EXTENSIONS):
                collected.append(os.path.join(root, f))
    return collected


# ---------- Main ----------

def validate_datasets(file_list):
    global_relationships = set()

    print("\n=== DATASET VALIDATION REPORT ===\n")

    for file in file_list:
        print(f"Checking: {file}")

        if file.endswith((".csv", ".tsv")):
            ok, msg = sample_csv(file, global_relationships)
        elif file.endswith(".json"):
            ok, msg = sample_json(file, global_relationships)
        elif file.endswith(".xml"):
            ok, msg = sample_xml(file, global_relationships)
        elif file.endswith((".fasta", ".fa")):
            ok, msg = sample_fasta(file, global_relationships)
        else:
            ok, msg = False, "Unsupported format"

        icon = "✅" if ok else "❌"
        print(f"{icon} {msg}\n")

    print("\n=== UNIQUE NORMALIZED RELATIONSHIPS ===\n")
    for r in sorted(global_relationships):
        print(r)

    print(f"\nTotal unique relationships: {len(global_relationships)}\n")


if __name__ == "__main__":
    files = discover_files(DATA_DIR)
    print(f"\nDiscovered {len(files)} dataset files\n")
    validate_datasets(files)
