#!/usr/bin/env python3
"""
merge_integrated_edges.py -- Merge Monarch edges_clean.csv with the newly
extracted non-Monarch edges into one integrated, deduplicated edge list.

Dedup key = (source_id, relation, target_id). For duplicate edges the
dataset_sources are unioned (so an edge supported by both Monarch and, say,
Gene2Phenotype is correctly attributed to both -- this also makes genuine
multi-source support visible for consensus weighting). Self-loops are dropped.

Output: data/processed/edges_clean_integrated.csv  (same schema as edges_clean)
"""
import pandas as pd
from config import EDGES_CLEAN_CSV, PROCESSED_DIR

ADDITIONAL = PROCESSED_DIR / "additional_edges.csv"
OUT = PROCESSED_DIR / "edges_clean_integrated.csv"


def norm_sources(series):
    """Union ';'-separated source tokens within a duplicate group."""
    toks = set()
    for v in series:
        if isinstance(v, str):
            toks.update(t for t in v.split(";") if t)
    return ";".join(sorted(toks))


def main():
    base = pd.read_csv(EDGES_CLEAN_CSV, dtype=str)
    add = pd.read_csv(ADDITIONAL, dtype=str)
    print(f"Monarch edges:    {len(base):,}")
    print(f"Additional edges: {len(add):,}")

    both = pd.concat([base, add], ignore_index=True)
    both = both[both["source_id"] != both["target_id"]]          # drop self-loops
    before = len(both)

    grp = (both.groupby(["source_id", "relation", "target_id"], sort=False)
                .agg(weight=("weight", "size"),
                     dataset_sources=("dataset_sources", norm_sources))
                .reset_index())
    grp["weight"] = grp["weight"].astype(int)

    grp.to_csv(OUT, index=False)
    n_multi = (grp["dataset_sources"].str.contains(";")).sum()
    n_nodes = len(set(grp["source_id"]) | set(grp["target_id"]))
    print(f"After concat (pre-dedup, no self-loops): {before:,}")
    print(f"Integrated unique edges: {len(grp):,}")
    print(f"  edges with >1 source token: {n_multi:,}")
    print(f"  nodes: {n_nodes:,}")
    print(f"  relations: {grp['relation'].nunique()}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
