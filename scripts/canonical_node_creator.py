import csv
from pathlib import Path

edges_path = "data/processed/edges_clean.csv"
nodes_path = "data/processed/canonical_nodes.csv"

nodes_set = set()

with open(edges_path, newline='', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        nodes_set.add(row['source_id'])
        nodes_set.add(row['target_id'])

Path(nodes_path).parent.mkdir(parents=True, exist_ok=True)
with open(nodes_path, 'w', newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=['canonical_id'])
    writer.writeheader()
    for nid in sorted(nodes_set):
        writer.writerow({'canonical_id': nid})

print(f"Created canonical_nodes.csv with {len(nodes_set):,} nodes")
