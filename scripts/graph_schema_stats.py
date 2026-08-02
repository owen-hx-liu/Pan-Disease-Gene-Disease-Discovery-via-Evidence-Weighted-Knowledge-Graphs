"""Precompute the schema/degree data that Figures S1 and S2 read.

Reads ``data/processed/edges_clean_integrated.csv`` once and writes a single
frozen JSON, ``data/processed/graph_schema.json``, containing everything the two
supplementary figures need -- so ``make_figures.py`` never re-touches the 360 MB
edge list and, like the other figures, only reads from a result JSON.

  * ``degree_counts``   -- {degree: n_nodes}; the full (integer) degree
                           distribution for Figure S2 (log-log).
  * ``metagraph``       -- directed category->category edges aggregated over all
                           relations, with the dominant relation per pair, for the
                           Figure S1 schema/metagraph.
  * ``category_nodes``  -- node count per category (sizes the S1 nodes).

Pandas/numpy only (runs in the restricted sandbox), fully deterministic.

Usage:  python scripts/graph_schema_stats.py
"""
from __future__ import annotations

import json
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kg_categories import category_of  # noqa: E402

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
INPUT_CSV = os.path.join(REPO_ROOT, "data", "processed", "edges_clean_integrated.csv")
OUTPUT_JSON = os.path.join(REPO_ROOT, "data", "processed", "graph_schema.json")

# Show the relation name without the BIOLINK: prefix in the schema figure.
def _short_relation(rel: str) -> str:
    return str(rel).split(":", 1)[-1]


def main() -> None:
    if not os.path.exists(INPUT_CSV):
        sys.exit(f"Input not found: {INPUT_CSV}")

    df = pd.read_csv(
        INPUT_CSV,
        usecols=["source_id", "relation", "target_id"],
        dtype={"source_id": "string", "relation": "string", "target_id": "string"},
    ).dropna(subset=["source_id", "target_id"]).reset_index(drop=True)

    src, tgt = df["source_id"], df["target_id"]

    # --- Degree distribution (endpoint multiset; self-loops count twice) ------
    degree = pd.concat([src, tgt], ignore_index=True).value_counts()
    nodes = degree.index
    deg_vc = degree.value_counts()  # degree value -> number of nodes with it
    degree_counts = {str(int(k)): int(v) for k, v in sorted(deg_vc.items())}

    # --- Node category (for S1 node sizing) -----------------------------------
    node_category = pd.Series(nodes, dtype="string").map(category_of)
    cat_vc = node_category.value_counts()
    category_nodes = {str(k): int(v) for k, v in sorted(cat_vc.items(), key=lambda kv: -int(kv[1]))}

    # --- Metagraph: category -> category edges over all relations -------------
    src_cat = src.map(category_of)
    tgt_cat = tgt.map(category_of)
    meta = pd.DataFrame({
        "s": src_cat.to_numpy(),
        "t": tgt_cat.to_numpy(),
        "rel": df["relation"].map(_short_relation).to_numpy(),
    })
    # Total directed edges per (source category, target category).
    pair_tot = meta.groupby(["s", "t"]).size()
    # Dominant relation per pair (deterministic: max count, then name).
    per_rel = meta.groupby(["s", "t", "rel"]).size().reset_index(name="n")
    dominant = (
        per_rel.sort_values(["s", "t", "n", "rel"], ascending=[True, True, False, True])
        .groupby(["s", "t"], as_index=False)
        .first()
    )
    dom_map = {(r.s, r.t): (r.rel, int(r.n)) for r in dominant.itertuples()}

    metagraph = []
    for (s, t), n in sorted(pair_tot.items(), key=lambda kv: -int(kv[1])):
        rel, rel_n = dom_map[(s, t)]
        metagraph.append({
            "source": str(s),
            "target": str(t),
            "n_edges": int(n),
            "dominant_relation": rel,
            "dominant_relation_edges": rel_n,
        })

    out = {
        "input_csv": os.path.relpath(INPUT_CSV, REPO_ROOT).replace("\\", "/"),
        "n_nodes": int(len(nodes)),
        "n_edges_directed": int(len(df)),
        "degree_counts": degree_counts,
        "category_nodes": category_nodes,
        "metagraph": metagraph,
    }
    os.makedirs(os.path.dirname(OUTPUT_JSON), exist_ok=True)
    with open(OUTPUT_JSON, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2, sort_keys=False)

    print(f"Nodes {out['n_nodes']:,} | directed edges {out['n_edges_directed']:,}")
    print(f"Distinct degrees: {len(degree_counts)} | categories: {len(category_nodes)}"
          f" | metagraph edges: {len(metagraph)}")
    print(f"Wrote {os.path.relpath(OUTPUT_JSON, REPO_ROOT)}")


if __name__ == "__main__":
    main()
