import csv
from collections import defaultdict

# ===================== CONFIG =====================

EDGES_INPUT = "data\processed\edges_deduplicated.csv"
ENTITY_MAP = "data\processed\entity_id_map.csv"
OUTPUT = "canonical_edges.csv"

# Relations that must NEVER self-loop
INVALID_SELF_LOOP_RELATIONS = {
    "biolink:orthologous_to",
    "biolink:interacts_with"
}

# ===================== LOAD ENTITY MAP =====================

print("Loading entity ID map...")

id_map = {}

with open(ENTITY_MAP, "r", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    for row in reader:
        raw = row["raw_id"].strip()
        canon = row["canonical_id"].strip()
        id_map[raw] = canon

print(f"Loaded {len(id_map):,} ID mappings")

# ===================== PROCESS EDGES =====================

print("Canonicalizing edges...")

total_edges = 0
canonicalized = 0
raw_leaks = 0
self_loops_removed = 0

written = 0

with open(EDGES_INPUT, "r", encoding="utf-8") as fin, \
     open(OUTPUT, "w", encoding="utf-8", newline="") as fout:

    reader = csv.DictReader(fin)
    fieldnames = reader.fieldnames
    writer = csv.DictWriter(fout, fieldnames=fieldnames)
    writer.writeheader()

    for row in reader:
        total_edges += 1

        src = row["source_id"].strip()
        tgt = row["target_id"].strip()
        rel = row["relation"].strip().lower()

        # Canonicalize
        src_c = id_map.get(src)
        tgt_c = id_map.get(tgt)

        if not src_c or not tgt_c:
            raw_leaks += 1
            continue

        row["source_id"] = src_c
        row["target_id"] = tgt_c
        canonicalized += 1

        # Remove invalid self-loops
        if src_c == tgt_c and rel in INVALID_SELF_LOOP_RELATIONS:
            self_loops_removed += 1
            continue

        writer.writerow(row)
        written += 1

# ===================== REPORT =====================

print("\n================ CANONICALIZATION REPORT ================")
print(f"Input edges:                {total_edges:,}")
print(f"Canonicalized edges:        {canonicalized:,}")
print(f"Raw ID leaks removed:       {raw_leaks:,}")
print(f"Invalid self-loops removed: {self_loops_removed:,}")
print(f"Final edges written:        {written:,}")
print("==========================================================")
print(f"Output file: {OUTPUT}")
