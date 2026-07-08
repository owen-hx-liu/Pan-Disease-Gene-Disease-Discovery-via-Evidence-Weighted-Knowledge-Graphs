import pandas as pd
from pathlib import Path
import json
import csv

RAW_DATA_DIR = Path("data/raw")
OUTPUT_CSV = Path("data/processed/relationship_vocab.csv")

all_relationships = set()

def scan_csv_tsv(file_path: Path):
    try:
        delimiter = "\t" if file_path.suffix.lower() == ".tsv" else ","
        df = pd.read_csv(file_path, sep=delimiter, low_memory=False, on_bad_lines="skip")
        for col in df.columns:
            if "relation" in col.lower() or "predicate" in col.lower() or "edge" in col.lower():
                all_relationships.update(df[col].dropna().astype(str).str.strip().unique())
    except Exception as e:
        print(f"Warning: Failed to process {file_path}: {e}")

def scan_json(file_path: Path):
    import json
    try:
        data = json.load(open(file_path, encoding="utf-8"))
        if isinstance(data, dict):
            data = [data]
        for obj in data:
            for key in obj:
                if "relation" in key.lower() or "predicate" in key.lower() or "edge" in key.lower():
                    all_relationships.add(str(obj[key]).strip())
    except Exception as e:
        print(f"Warning: Failed to process {file_path}: {e}")

for dataset_folder in RAW_DATA_DIR.iterdir():
    if not dataset_folder.is_dir():
        continue
    for file_path in dataset_folder.rglob("*"):
        if not file_path.is_file():
            continue
        if file_path.suffix.lower() in [".csv", ".tsv"]:
            scan_csv_tsv(file_path)
        elif file_path.suffix.lower() == ".json":
            scan_json(file_path)

relationship_vocab = [r.upper().replace(" ", "_") for r in all_relationships if r.strip()]

OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["relationship"])
    for r in sorted(relationship_vocab):
        writer.writerow([r])

print(f"3.2 COMPLETE — Found {len(relationship_vocab)} unique relationships")
print(f"Saved to: {OUTPUT_CSV}")
