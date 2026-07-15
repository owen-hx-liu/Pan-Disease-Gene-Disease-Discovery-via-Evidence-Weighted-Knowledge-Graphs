"""Deterministic graph statistics for the Hetionet robustness graph (pandas only).

Reads ``data/external/hetionet/hetionet-v1.0-edges.sif`` (tab-separated columns
``source,metaedge,target`` where every ID is ``Type::id``) and writes a summary
to ``data/processed/graph_stats_hetionet.json``, also printing a report.

Hetionet is the second, independent graph used for the cross-graph robustness
check (C4). It is NEVER merged with the Monarch graph. Node "category" here is
simply the metanode string before ``::`` (NOT kg_categories.category_of, which
is Monarch-specific). Output structure mirrors ``graph_stats.py`` so the tables
builder can consume Monarch and Hetionet uniformly.

No network/graph libraries are used; fully deterministic (sorted outputs).

Usage:  python scripts/graph_stats_hetionet.py
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd

# --- Paths (relative to repo root; no absolute paths in deliverables) --------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
INPUT_SIF = os.path.join(
    REPO_ROOT, "data", "external", "hetionet", "hetionet-v1.0-edges.sif"
)
NODES_TSV = os.path.join(
    REPO_ROOT, "data", "external", "hetionet", "hetionet-v1.0-nodes.tsv"
)
OUTPUT_JSON = os.path.join(REPO_ROOT, "data", "processed", "graph_stats_hetionet.json")

# The gene->disease target metaedge (Disease-associates-Gene) for the audit task.
TARGET_METAEDGE = "DaG"


def _metanode(node_id: str) -> str:
    """Hetionet metanode = text before the '::' (e.g. 'Gene', 'Disease')."""
    return node_id.split("::", 1)[0]


def _counts_dict(series: pd.Series) -> dict:
    """value_counts as a plain dict, sorted by count desc then key asc (deterministic)."""
    vc = series.value_counts()
    items = sorted(vc.items(), key=lambda kv: (-int(kv[1]), str(kv[0])))
    return {str(k): int(v) for k, v in items}


def main() -> None:
    if not os.path.exists(INPUT_SIF):
        sys.exit(f"Input not found: {INPUT_SIF}")

    df = pd.read_csv(
        INPUT_SIF,
        sep="\t",
        dtype={"source": "string", "metaedge": "string", "target": "string"},
    )
    df = df.dropna(subset=["source", "target"]).reset_index(drop=True)

    src = df["source"]
    tgt = df["target"]

    # --- Edge counts ---------------------------------------------------------
    n_edges_directed = int(len(df))

    self_loop_mask = (src == tgt).to_numpy()
    n_self_loops = int(self_loop_mask.sum())

    # Unique undirected edges: order each pair (lo,hi) then dedup. Vectorized
    # with numpy -- pandas row-wise min/max(axis=1) over millions of strings is
    # prohibitively slow.
    a = src[~self_loop_mask].to_numpy()
    b = tgt[~self_loop_mask].to_numpy()
    swap = a > b
    lo = np.where(swap, b, a)
    hi = np.where(swap, a, b)
    n_edges_unique_undirected = int(
        len(set(zip(lo.tolist(), hi.tolist())))
    )

    # --- Nodes (degree = incident endpoints, in+out) -------------------------
    endpoints = pd.concat([src, tgt], ignore_index=True)
    degree = endpoints.value_counts()
    nodes = degree.index

    n_nodes = int(len(nodes))
    # Full node catalog (includes isolated nodes not present in any edge); this
    # is the canonical published Hetionet figure (~47,031).
    n_nodes_catalog = None
    if os.path.exists(NODES_TSV):
        n_nodes_catalog = int(
            pd.read_csv(NODES_TSV, sep="\t", usecols=["id"]).shape[0]
        )
    degree_mean = float(degree.mean())
    degree_median = float(degree.median())
    degree_max = int(degree.max())

    # --- Per-metanode node counts (category = string before '::') ------------
    node_category = pd.Series(nodes, dtype="string").map(_metanode)
    per_category = _counts_dict(node_category)

    # --- Per-metaedge edge counts --------------------------------------------
    per_relation = _counts_dict(df["metaedge"].fillna("<NA>"))
    n_target_edges = int(per_relation.get(TARGET_METAEDGE, 0))

    stats = {
        "graph": "hetionet-v1.0",
        "input_sif": os.path.relpath(INPUT_SIF, REPO_ROOT).replace("\\", "/"),
        "never_merged_with": "monarch (edges_clean_integrated.csv)",
        "n_nodes": n_nodes,
        "n_nodes_catalog": n_nodes_catalog,
        "n_edges_directed": n_edges_directed,
        "n_edges_unique_undirected": n_edges_unique_undirected,
        "n_self_loops": n_self_loops,
        "degree": {
            "mean": degree_mean,
            "median": degree_median,
            "max": degree_max,
        },
        "target_metaedge": TARGET_METAEDGE,
        "n_target_edges": n_target_edges,
        "n_categories": len(per_category),
        "n_relations": len(per_relation),
        "per_category_node_counts": per_category,
        "per_relation_edge_counts": per_relation,
    }

    os.makedirs(os.path.dirname(OUTPUT_JSON), exist_ok=True)
    with open(OUTPUT_JSON, "w", encoding="utf-8") as fh:
        json.dump(stats, fh, indent=2, sort_keys=False)

    # --- Human-readable summary ---------------------------------------------
    print("=" * 64)
    print("Graph statistics -- hetionet-v1.0-edges.sif")
    print("=" * 64)
    print(f"Nodes (in >=1 edge)         : {n_nodes:,}")
    if n_nodes_catalog is not None:
        print(f"Nodes (full catalog)        : {n_nodes_catalog:,}")
    print(f"Edges (directed rows)       : {n_edges_directed:,}")
    print(f"Edges (unique undirected)   : {n_edges_unique_undirected:,}")
    print(f"Self-loops                  : {n_self_loops:,}")
    print(f"Degree mean / median / max  : {degree_mean:.2f} / {degree_median:.1f} / {degree_max:,}")
    print(f"Metanodes / metaedges       : {len(per_category)} / {len(per_relation)}")
    print(f"Target metaedge {TARGET_METAEDGE!r:>5} edges : {n_target_edges:,}")

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

    _print_block(f"Per-metanode node counts ({len(per_category)})", per_category)
    _print_block(f"Per-metaedge edge counts ({len(per_relation)})", per_relation, top=30)

    print(f"\nWrote {os.path.relpath(OUTPUT_JSON, REPO_ROOT)}")


if __name__ == "__main__":
    main()
