import csv
import sys
from pathlib import Path

INPUT_PATH = Path("data/processed/edges_clean.csv")
OUTPUT_PATH = Path("data/processed/relationships.csv")

if not INPUT_PATH.exists():
    sys.exit(f"Missing input file: {INPUT_PATH}")

print("Loading cleaned edges...")

edges_written = 0

OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)

with open(INPUT_PATH, newline="", encoding="utf-8") as infile, \
     open(OUTPUT_PATH, "w", newline="", encoding="utf-8") as outfile:

    reader = csv.DictReader(infile)

    required = {
        "source_id",
        "relation",
        "target_id",
        "weight",
        "dataset_sources"
    }

    missing = required - set(reader.fieldnames)
    if missing:
        sys.exit(f"Missing required columns: {missing}")

    writer = csv.writer(outfile)

    writer.writerow([
        ":START_ID",
        ":END_ID",
        ":TYPE",
        "weight:int",
        "dataset_sources"
    ])

    for row in reader:
        src = row["source_id"].strip()
        tgt = row["target_id"].strip()

        if not src or not tgt:
            continue

        rel = row["relation"].strip().upper()
        if rel.startswith("BIOLINK:"):
            rel = rel.replace("BIOLINK:", "")
        rel = rel.replace(" ", "_")

        weight = row.get("weight", "1")
        sources = row.get("dataset_sources", "")

        writer.writerow([src, tgt, rel, weight, sources])
        edges_written += 1

print(f"Relationships written: {edges_written:,}")
print(f"relationships.csv written to {OUTPUT_PATH}")
print("Neo4j relationships file ready")
