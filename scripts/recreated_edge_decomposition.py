#!/usr/bin/env python
"""
How much of the surviving R2 score is the ~3% of held-out test triples that the
degree-preserving rewiring recreates by chance, and how much is the entity-
frequency prior over the remaining 97%?

The rewiring is degree-biased by construction, so a swap can land a held-out
(gene, relation, disease) triple back in the training graph. Those edges are
genuinely leaked and they score near-perfectly. This script splits every stored
per-edge reciprocal rank on the canonical R2 graph into the recreated subset and
the rest, so the manuscript can state the two contributions separately instead of
leaving a reader to wonder whether the residual MRR *is* the leak.

Reads only frozen result files and the frozen splits; writes one JSON.

    python scripts/recreated_edge_decomposition.py
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SPLITS = ROOT / "data" / "processed" / "splits"
RESULTS = ROOT / "data" / "processed" / "results"
OUT = RESULTS / "null" / "recreated_edge_decomposition.json"

TEST = SPLITS / "test.csv"
R2_TRAIN = SPLITS / "train_R2_degree_null_seed42.csv"

KGE_MODELS = ["TransE", "RotatE", "DistMult", "ComplEx"]
BASELINES = ["Random", "CommonNeighbors", "AdamicAdar", "Jaccard",
             "PreferentialAttachment"]
SEEDS = [42, 1, 7]


def read_triples(path: Path) -> list[tuple[str, str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return [(r["source_id"], r["relation"], r["target_id"])
                for r in csv.DictReader(fh)]


def recreated_mask() -> np.ndarray:
    """True where the canonical R2 rewiring put the held-out triple back."""
    test = read_triples(TEST)
    train = set(read_triples(R2_TRAIN))
    return np.array([t in train for t in test], dtype=bool)


def split_mrr(rr: np.ndarray, mask: np.ndarray) -> dict:
    """MRR overall, on the recreated edges, and on everything else."""
    overall = float(rr.mean())
    kept = float(rr[~mask].mean())
    return {
        "MRR_all": overall,
        "MRR_recreated": float(rr[mask].mean()),
        "MRR_excluding_recreated": kept,
        # How much of the reported R2 MRR the recreated edges account for.
        "contribution_of_recreated": overall - kept,
        "contribution_pct_of_R2": 100.0 * (overall - kept) / overall if overall else 0.0,
    }


def main() -> None:
    mask = recreated_mask()
    n_rec = int(mask.sum())
    print(f"canonical R2 recreates {n_rec}/{mask.size} test triples "
          f"({100 * n_rec / mask.size:.2f}%)")

    out: dict = {
        "description": "R2 MRR split into the recreated held-out triples and the rest.",
        "canonical_R2_file": R2_TRAIN.name,
        "n_test_edges": int(mask.size),
        "n_recreated": n_rec,
        "frac_recreated": n_rec / mask.size,
        "protocol": "sampled-50neg (the main protocol)",
        "methods": {},
    }

    # --- KGE arm: per-edge ranks live in one file per (model, regime, seed) ---
    for model in KGE_MODELS:
        r0, r2 = [], []
        for seed in SEEDS:
            d0 = json.loads((RESULTS / "kge" / f"kge_{model}_R0_seed{seed}.json").read_text())
            d2 = json.loads((RESULTS / "kge" / f"kge_{model}_R2_seed{seed}.json").read_text())
            r0.append(np.asarray(d0["reciprocal_ranks"], dtype=float))
            r2.append(np.asarray(d2["reciprocal_ranks"], dtype=float))
        rec = split_mrr(np.mean(r2, axis=0), mask)
        rec["MRR_R0"] = float(np.mean(r0))
        rec["collapse_pct_reported"] = 100.0 * (rec["MRR_all"] / rec["MRR_R0"] - 1.0)
        rec["collapse_pct_excluding_recreated"] = 100.0 * (
            rec["MRR_excluding_recreated"] / rec["MRR_R0"] - 1.0)
        out["methods"][model] = rec

    # --- topological arm: one file per regime, all methods and seeds inside ---
    bl0 = json.loads((RESULTS / "baselines" / "baselines_R0.json").read_text())
    bl2 = json.loads((RESULTS / "baselines" / "baselines_R2.json").read_text())
    for method in BASELINES:
        rr2 = np.mean(np.asarray(bl2["per_edge_reciprocal_ranks"][method], dtype=float), axis=0)
        rr0 = np.mean(np.asarray(bl0["per_edge_reciprocal_ranks"][method], dtype=float), axis=0)
        rec = split_mrr(rr2, mask)
        rec["MRR_R0"] = float(rr0.mean())
        rec["collapse_pct_reported"] = 100.0 * (rec["MRR_all"] / rec["MRR_R0"] - 1.0)
        rec["collapse_pct_excluding_recreated"] = 100.0 * (
            rec["MRR_excluding_recreated"] / rec["MRR_R0"] - 1.0)
        out["methods"][method] = rec

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2))
    print(f"wrote {OUT}")

    hdr = f"{'method':<24}{'R0':>8}{'R2':>8}{'R2 excl':>9}{'recreated':>11}{'drop':>8}{'drop excl':>11}"
    print(hdr)
    for name, m in out["methods"].items():
        print(f"{name:<24}{m['MRR_R0']:>8.3f}{m['MRR_all']:>8.3f}"
              f"{m['MRR_excluding_recreated']:>9.3f}{m['MRR_recreated']:>11.3f}"
              f"{m['collapse_pct_reported']:>8.1f}{m['collapse_pct_excluding_recreated']:>11.1f}")


if __name__ == "__main__":
    main()
