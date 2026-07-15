"""Deterministic graph statistics for the integrated KG edge list (pandas only).

Reads ``data/processed/edges_clean_integrated.csv`` (columns
``source_id,relation,target_id,weight,dataset_sources``) and writes a summary to
``data/processed/graph_stats.json``, also printing a human-readable report.

No network/graph libraries (networkx/igraph/scipy) are used -- everything is
computed with pandas/numpy so it runs in the restricted sandbox. Results are
fully deterministic (no sampling, sorted outputs).

Usage:  python scripts/graph_stats.py
"""
from __future__ import annotations

import json
import os
import sys

import pandas as pd

# Import the shared category lookup (same scripts/ directory).
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kg_categories import category_of  # noqa: E402

# --- Paths (relative to repo root; no absolute paths in deliverables) --------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
INPUT_CSV = os.path.join(REPO_ROOT, "data", "processed", "edges_clean_integrated.csv")
OUTPUT_JSON = os.path.join(REPO_ROOT, "data", "processed", "graph_stats.json")


def _prefix(node_id: str) -> str:
    """Namespace prefix = text before the first ':' (upper-cased)."""
    return node_id.split(":", 1)[0].upper()


def _counts_dict(series: pd.Series) -> dict:
    """value_counts as a plain dict, sorted by count desc then key asc (deterministic)."""
    vc = series.value_counts()
    items = sorted(vc.items(), key=lambda kv: (-int(kv[1]), str(kv[0])))
    return {str(k): int(v) for k, v in items}


def main() -> None:
    if not os.path.exists(INPUT_CSV):
        sys.exit(f"Input not found: {INPUT_CSV}")

    # Only the columns we need; strings for IDs/relation.
    df = pd.read_csv(
        INPUT_CSV,
        usecols=["source_id", "relation", "target_id"],
        dtype={"source_id": "string", "relation": "string", "target_id": "string"},
    )
    # Drop rows with missing endpoints (shouldn't happen in a clean file).
    df = df.dropna(subset=["source_id", "target_id"]).reset_index(drop=True)

    src = df["source_id"]
    tgt = df["target_id"]

    # --- Edge counts ---------------------------------------------------------
    n_edges_directed = int(len(df))

    self_loop_mask = (src == tgt).to_numpy()
    n_self_loops = int(self_loop_mask.sum())

    # Unique undirected edges = distinct unordered node pairs, excluding self-loops.
    non_loop = df.loc[~self_loop_mask, ["source_id", "target_id"]]
    lo = non_loop[["source_id", "target_id"]].min(axis=1)
    hi = non_loop[["source_id", "target_id"]].max(axis=1)
    n_edges_unique_undirected = int(
        pd.DataFrame({"lo": lo, "hi": hi}).drop_duplicates().shape[0]
    )

    # --- Nodes ---------------------------------------------------------------
    # Endpoint multiset: node degree = number of incident edge-endpoints
    # (in + out; a self-loop contributes 2 to that node, as it should).
    endpoints = pd.concat([src, tgt], ignore_index=True)
    degree = endpoints.value_counts()  # index=node, value=total degree
    nodes = degree.index

    n_nodes = int(len(nodes))
    degree_mean = float(degree.mean())
    degree_median = float(degree.median())
    degree_max = int(degree.max())

    # --- Per-namespace node counts (prefix before ':') -----------------------
    node_prefix = pd.Series(nodes, dtype="string").map(_prefix)
    per_namespace = _counts_dict(node_prefix)

    # --- Per-category node counts (via kg_categories.category_of) ------------
    node_category = pd.Series(nodes, dtype="string").map(category_of)
    per_category = _counts_dict(node_category)

    # --- Per-relation edge counts --------------------------------------------
    per_relation = _counts_dict(df["relation"].fillna("<NA>"))

    stats = {
        "input_csv": os.path.relpath(INPUT_CSV, REPO_ROOT).replace("\\", "/"),
        "n_nodes": n_nodes,
        "n_edges_directed": n_edges_directed,
        "n_edges_unique_undirected": n_edges_unique_undirected,
        "n_self_loops": n_self_loops,
        "degree": {
            "mean": degree_mean,
            "median": degree_median,
            "max": degree_max,
        },
        "n_namespaces": len(per_namespace),
        "n_relations": len(per_relation),
        "n_categories": len(per_category),
        "per_namespace_node_counts": per_namespace,
        "per_category_node_counts": per_category,
        "per_relation_edge_counts": per_relation,
    }

    os.makedirs(os.path.dirname(OUTPUT_JSON), exist_ok=True)
    with open(OUTPUT_JSON, "w", encoding="utf-8") as fh:
        json.dump(stats, fh, indent=2, sort_keys=False)

    # --- Human-readable summary ---------------------------------------------
    print("=" * 64)
    print("Graph statistics -- edges_clean_integrated.csv")
    print("=" * 64)
    print(f"Nodes (unique)              : {n_nodes:,}")
    print(f"Edges (directed rows)       : {n_edges_directed:,}")
    print(f"Edges (unique undirected)   : {n_edges_unique_undirected:,}")
    print(f"Self-loops                  : {n_self_loops:,}")
    print(f"Degree mean / median / max  : {degree_mean:.2f} / {degree_median:.1f} / {degree_max:,}")

    def _print_block(title: str, d: dict, top: int | None = None) -> None:
        print("\n" + title)
        print("-" * len(title))
        items = list(d.items())
        shown = items if top is None else items[:top]
        width = max((len(k) for k, _ in shown), default=0)
        for k, v in shown:
            print(f"  {k:<{width}}  {v:>12,}")
        if top is not None and len(items) > top:
            print(f"  ... ({len(items) - top} more)")

    _print_block(f"Per-category node counts ({len(per_category)})", per_category)
    _print_block(f"Per-namespace node counts ({len(per_namespace)})", per_namespace)
    _print_block(f"Per-relation edge counts ({len(per_relation)})", per_relation, top=25)

    print(f"\nWrote {os.path.relpath(OUTPUT_JSON, REPO_ROOT)}")


if __name__ == "__main__":
    main()
