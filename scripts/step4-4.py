import pandas as pd
from collections import defaultdict

# ================= CONFIG =================

INPUT_EDGES = "data/processed/canonical_edges.csv"
OUTPUT_EDGES = "data/processed/edges_deduplicated.csv"

REQUIRED_COLUMNS = {
    "source_id",
    "relation",
    "target_id",
    "dataset_id"
}

# ================= LOAD =================

print("Loading canonical edges...")
edges = pd.read_csv(INPUT_EDGES, sep="|")


missing = REQUIRED_COLUMNS - set(edges.columns)
if missing:
    raise ValueError(f"Missing required columns: {missing}")

# ================= GROUP =================

edge_groups = defaultdict(set)

for _, row in edges.iterrows():
    key = (
        str(row["source_id"]).strip(),
        str(row["relation"]).strip(),
        str(row["target_id"]).strip()
    )
    edge_groups[key].add(str(row["dataset_id"]).strip())

# ================= OUTPUT =================

rows = []

for (src, rel, tgt), datasets in edge_groups.items():
    rows.append({
        "source_id": src,
        "relation": rel,
        "target_id": tgt,
        "weight": len(datasets),
        "dataset_sources": ";".join(sorted(datasets))
    })

out_df = pd.DataFrame(rows)
out_df.to_csv(OUTPUT_EDGES, index=False)

print("✅ STEP 4.4 COMPLETE")
print(f"Original edges: {len(edges)}")
print(f"Deduplicated edges: {len(out_df)}")
print(f"Saved to: {OUTPUT_EDGES}")
