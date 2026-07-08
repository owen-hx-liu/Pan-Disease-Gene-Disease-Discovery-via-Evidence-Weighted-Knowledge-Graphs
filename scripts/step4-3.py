#!/usr/bin/env python3
"""
Step 4.3 — Canonicalize edges using global ID mapping

Input:
- extracted_edges_*.csv
- entity_id_map.csv

Output:
- canonical_edges.csv
"""

import pandas as pd
from pathlib import Path
import csv
from collections import defaultdict

# ===================== CONFIG =====================

EDGE_FILE = Path("data/processed/extracted_edges_20260205_213914.csv")
ID_MAP_FILE = Path("data/processed/entity_id_map.csv")
OUTPUT_FILE = Path("data/processed/canonical_edges.csv")

CHUNK_SIZE = 200_000

# ===================== LOAD ID MAP =====================

print("Loading entity ID map...")

id_map_df = pd.read_csv(ID_MAP_FILE)

RAW_TO_CANON = dict(
    zip(
        id_map_df["raw_id"].astype(str),
        id_map_df["canonical_id"].astype(str)
    )
)

print(f"Loaded {len(RAW_TO_CANON):,} ID mappings")

# ===================== PROCESS EDGES =====================

dropped = defaultdict(int)
written = 0

with open(OUTPUT_FILE, "w", newline="", encoding="utf-8") as out_f:
    writer = csv.writer(out_f, delimiter="|")
    writer.writerow(["source_id", "relation", "target_id", "dataset_id"])

    for chunk in pd.read_csv(
        EDGE_FILE,
        delimiter="|",
        chunksize=CHUNK_SIZE,
        dtype=str
    ):
        for _, row in chunk.iterrows():

            src_raw = row["source_raw_id"]
            tgt_raw = row["target_raw_id"]
            rel = row["relation"]
            dataset = row["dataset_id"]

            src = RAW_TO_CANON.get(src_raw)
            tgt = RAW_TO_CANON.get(tgt_raw)

            if not src:
                dropped["missing_source"] += 1
                continue

            if not tgt:
                dropped["missing_target"] += 1
                continue

            writer.writerow([src, rel, tgt, dataset])
            written += 1

print("\n✅ 4.3 COMPLETE")
print(f"Canonical edges written: {written:,}")

if dropped:
    print("\nDropped edges:")
    for k, v in dropped.items():
        print(f"  {k}: {v:,}")
