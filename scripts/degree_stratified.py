"""Degree-stratified re-analysis of the held-out gene-disease test set (Unit: reviewer #5).

Reviewer request: "Run a low-degree vs high-degree analysis, since it directly tests
your central claim." The central claim is that node degree is the universal leak.

This needs NO retraining. Every KGE result JSON already stores ``reciprocal_ranks`` --
one tie-averaged reciprocal rank per test edge, in the fixed row order of
``splits/test.csv`` (verified: len == 4228, mean == the file's reported MRR). We simply:

  1. Look up each test edge's gene and disease degree in the R0 TRAINING graph
     (splits/node_degree.csv -- the exact popularity the model saw; entities absent from
     training are the cold-start cases and get degree 0).
  2. Bin the 4,228 test edges by that degree into equal-count quartiles Q1..Q4 (lowest ->
     highest), plus a headline bottom-half vs top-half split.
  3. Recompute MRR and Hits@10 WITHIN each stratum, per KGE model and regime, per seed,
     then report mean +/- SD across the three seeds (42, 1, 7) -- the same uncertainty
     convention the rest of the paper uses.

We stratify along TWO axes and this distinction is the finding:
  * gene (query, source) degree -- essentially FLAT (hub genes are not easier), a clean
    negative control.
  * disease (ranked target) degree -- a clear monotonic gradient. Under the sampled-
    negative protocol the true disease is ranked against type-matched negative diseases,
    so it is the disease's popularity that the ranking can exploit. Crossed with the R2
    degree-null this localises the leak: high-degree (popular) diseases barely drop under
    R2 (their score is almost pure degree), whereas the rare diseases that motivate gene
    prioritisation collapse ~60% (their score was almost all genuine structure).

Only the two headline KGE models (TransE, RotatE) are covered here, because only the KGE
runs persist per-edge ranks. The topological baselines store aggregate MRR only; adding
them needs a one-line per-edge dump in run_baselines.py + a fast re-run, done separately.

Output: ``data/processed/results/degree_stratified/degree_stratified.json`` -- consumed by
make_tables.py (Table 5) and make_figures.py (Fig 4). No number is hand-typed downstream.

Reusing this on another method arm
----------------------------------
The stratification only needs per-edge reciprocal ranks in test.csv row order, so any arm
that writes the harness's run-JSON schema can be re-analysed with NO retraining. The
optional flags below point the loader at a different results directory / model list; every
default reproduces the frozen production behaviour byte-for-byte, so plain
``python scripts/degree_stratified.py`` is unchanged. The message-passing arm uses:

    python scripts/degree_stratified.py --kge-dir data/processed/results/gnn \\
        --models RGCN --regimes R0 R2 --no-baselines \\
        --out data/processed/results/gnn/degree_stratified.json

Usage:  python scripts/degree_stratified.py
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
RESULTS = os.path.join(REPO_ROOT, "data", "processed", "results")
SPLITS = os.path.join(REPO_ROOT, "data", "processed", "splits")
OUT_DIR = os.path.join(RESULTS, "degree_stratified")
# Where the per-run JSONs live, and where the report goes. Overridable via --kge-dir /
# --out; the defaults are the frozen production locations.
KGE_DIR = os.path.join(RESULTS, "kge")
OUT_PATH = os.path.join(OUT_DIR, "degree_stratified.json")
WITH_BASELINES = True

# Canonical presentation order (topological baselines first, then KGE), so the emitted
# tables/figures group methods by class.
BASELINES = ["Random", "CommonNeighbors", "AdamicAdar", "Jaccard", "PreferentialAttachment"]
KGE_MODELS = ["TransE", "RotatE"]
METHOD_ORDER = BASELINES + KGE_MODELS
# Class every KGE the repo trains, not just the two the frozen Table 5 stratifies. A model
# absent here falls through to parse_args' "message_passing" default, which would mislabel a
# matched KGE reference control passed via --models as message passing -- inverting the one
# contrast the GNN arm exists to make. Extra keys are inert: only `covered` models are emitted.
ALL_KGE_MODELS = KGE_MODELS + ["DistMult", "ComplEx"]
METHOD_CLASS = {**{m: "baseline" for m in BASELINES}, **{m: "kge" for m in ALL_KGE_MODELS}}
REGIMES = ["R0", "R1", "R2", "R3"]
SEEDS = [42, 1, 7]
N_QUANTILES = 4  # quartile gradient Q1 (low degree) .. Q4 (high degree)


# --------------------------------------------------------------------------- IO
def _die(msg: str) -> None:
    sys.exit(f"[degree_stratified] {msg}")


def load_test_edges() -> tuple[list[str], list[str]]:
    """(genes, diseases) for each test edge, in the row order of test.csv -- the order the
    per-edge reciprocal_ranks arrays are aligned to. source_id = gene (the query),
    target_id = disease (the entity actually ranked under the sampled-negative protocol)."""
    path = os.path.join(SPLITS, "test.csv")
    if not os.path.exists(path):
        _die(f"missing {path}")
    genes, diseases = [], []
    with open(path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            genes.append(row["source_id"])
            diseases.append(row["target_id"])
    return genes, diseases


def load_degrees() -> dict[str, int]:
    """{node_id: degree} over the R0 training graph (splits/node_degree.csv)."""
    path = os.path.join(SPLITS, "node_degree.csv")
    if not os.path.exists(path):
        _die(f"missing {path}")
    deg: dict[str, int] = {}
    with open(path, encoding="utf-8") as fh:
        reader = csv.reader(fh)
        next(reader)  # header: ",degree"
        for row in reader:
            if len(row) >= 2 and row[0]:
                deg[row[0]] = int(row[1])
    return deg


def load_reciprocal_ranks(model: str, regime: str) -> np.ndarray | None:
    """(n_seeds, n_test) matrix of per-edge reciprocal ranks; None if any seed is missing."""
    arrs = []
    for seed in SEEDS:
        path = os.path.join(KGE_DIR, f"kge_{model}_{regime}_seed{seed}.json")
        if not os.path.exists(path):
            return None
        with open(path, encoding="utf-8") as fh:
            rr = json.load(fh).get("reciprocal_ranks")
        if rr is None:
            return None
        arrs.append(np.asarray(rr, dtype=float))
    lens = {a.shape[0] for a in arrs}
    if len(lens) != 1:
        _die(f"{model}/{regime}: seed files disagree on length {lens}")
    return np.vstack(arrs)


def load_baseline_reciprocal_ranks(regime: str) -> dict[str, np.ndarray]:
    """{method: (n_seeds, n_test) matrix} for the topological baselines of one regime.

    Reads the ``per_edge_reciprocal_ranks`` block that run_baselines.py now persists (one
    list per seed, aligned to test.csv order). Returns {} if the file predates that block,
    so an old baselines_<regime>.json degrades gracefully (baselines simply omitted)."""
    path = os.path.join(RESULTS, "baselines", f"baselines_{regime}.json")
    if not os.path.exists(path):
        return {}
    with open(path, encoding="utf-8") as fh:
        payload = json.load(fh)
    per = payload.get("per_edge_reciprocal_ranks")
    if not per:
        return {}
    out: dict[str, np.ndarray] = {}
    for method, seed_lists in per.items():
        out[method] = np.asarray(seed_lists, dtype=float)  # (n_seeds, n_test)
    return out


# --------------------------------------------------------------- stratification
def build_strata(gene_degree: np.ndarray) -> dict:
    """Assign each test edge to an equal-count quartile by gene degree, and to a
    bottom-half / top-half headline split. Returns index masks and human-readable ranges."""
    n = gene_degree.shape[0]
    order = np.argsort(gene_degree, kind="stable")
    # Equal-count quartiles: split the degree-sorted edge list into N_QUANTILES chunks.
    chunks = np.array_split(order, N_QUANTILES)
    strata: dict = {"quartiles": {}, "halves": {}}
    for i, idx in enumerate(chunks, start=1):
        d = gene_degree[idx]
        strata["quartiles"][f"Q{i}"] = {
            "idx": idx,
            "degree_min": int(d.min()),
            "degree_max": int(d.max()),
            "n_edges": int(idx.shape[0]),
        }
    half = n // 2
    lo_idx, hi_idx = order[:half], order[half:]
    strata["halves"]["low"] = {
        "idx": lo_idx, "degree_min": int(gene_degree[lo_idx].min()),
        "degree_max": int(gene_degree[lo_idx].max()), "n_edges": int(lo_idx.shape[0])}
    strata["halves"]["high"] = {
        "idx": hi_idx, "degree_min": int(gene_degree[hi_idx].min()),
        "degree_max": int(gene_degree[hi_idx].max()), "n_edges": int(hi_idx.shape[0])}
    return strata


def stratum_metrics(rr: np.ndarray, idx: np.ndarray) -> dict:
    """MRR and Hits@10 within one stratum, mean +/- SD across seeds.

    rr is (n_seeds, n_test). rank = 1 / reciprocal_rank (the harness stores tie-averaged
    reciprocal ranks), so Hits@10 = fraction of edges with rank <= 10."""
    sub = rr[:, idx]                       # (n_seeds, n_edges_in_stratum)
    mrr_per_seed = sub.mean(axis=1)        # one MRR per seed
    rank = 1.0 / sub
    hits10_per_seed = (rank <= 10).mean(axis=1)
    return {
        "MRR_mean": float(mrr_per_seed.mean()),
        "MRR_sd": float(mrr_per_seed.std(ddof=0)),
        "Hits@10_mean": float(hits10_per_seed.mean()),
        "Hits@10_sd": float(hits10_per_seed.std(ddof=0)),
        "n_edges": int(idx.shape[0]),
    }


def stratify_axis(axis_name: str, node_degree: np.ndarray, cells: dict) -> dict:
    """Run the quartile + half stratification along one degree axis (gene or disease).

    node_degree : per-test-edge degree of the entity we are stratifying on.
    cells       : {"Model|Regime": rr matrix (n_seeds, n_test)} already loaded.
    """
    strata = build_strata(node_degree)
    q = strata["quartiles"]
    block: dict = {
        "axis": axis_name,
        "strata": {
            "quartiles": {k: {kk: vv for kk, vv in v.items() if kk != "idx"}
                          for k, v in q.items()},
            "halves": {k: {kk: vv for kk, vv in v.items() if kk != "idx"}
                       for k, v in strata["halves"].items()},
        },
        "results": {},
    }
    for key, rr in cells.items():
        cell = {"overall": stratum_metrics(rr, np.arange(rr.shape[1])),
                "quartiles": {}, "halves": {}}
        for name, s in q.items():
            cell["quartiles"][name] = stratum_metrics(rr, s["idx"])
        for name, s in strata["halves"].items():
            cell["halves"][name] = stratum_metrics(rr, s["idx"])
        block["results"][key] = cell
    return block


def parse_args() -> None:
    """Optional overrides for re-analysing a different method arm. Every default equals the
    frozen production configuration, so calling the script with no arguments is unchanged."""
    global KGE_DIR, OUT_PATH, KGE_MODELS, REGIMES, SEEDS, METHOD_ORDER, METHOD_CLASS
    global WITH_BASELINES
    ap = argparse.ArgumentParser(description="Degree-stratified re-analysis of stored "
                                             "per-edge ranks (no retraining).")
    ap.add_argument("--kge-dir", default=KGE_DIR,
                    help="directory holding kge_<model>_<regime>_seed<k>.json run files")
    ap.add_argument("--models", nargs="+", default=KGE_MODELS,
                    help="models to stratify (must store per-edge reciprocal_ranks)")
    ap.add_argument("--regimes", nargs="+", default=REGIMES)
    ap.add_argument("--seeds", nargs="+", type=int, default=SEEDS)
    ap.add_argument("--out", default=OUT_PATH, help="output JSON path")
    ap.add_argument("--no-baselines", action="store_true",
                    help="skip the topological-baseline per-edge dumps (they belong to the "
                         "full-graph arm and are not matched to a subgraph run)")
    a = ap.parse_args()
    KGE_DIR, OUT_PATH, KGE_MODELS = a.kge_dir, a.out, list(a.models)
    REGIMES, SEEDS = list(a.regimes), list(a.seeds)
    WITH_BASELINES = not a.no_baselines
    # A model this script has never heard of is a non-baseline method arm; class it as
    # message passing (the reason these flags exist) rather than crashing on a KeyError.
    for m in KGE_MODELS:
        METHOD_CLASS.setdefault(m, "message_passing")
        if m not in METHOD_ORDER:
            METHOD_ORDER.append(m)


def main() -> None:
    parse_args()
    print("Degree-stratified re-analysis (no retraining; reads frozen per-edge ranks)...")
    genes, diseases = load_test_edges()
    degrees = load_degrees()
    gene_deg = np.array([degrees.get(g, 0) for g in genes], dtype=float)
    dis_deg = np.array([degrees.get(d, 0) for d in diseases], dtype=float)
    n_gene_missing = sum(1 for g in genes if g not in degrees)
    n_dis_missing = sum(1 for d in diseases if d not in degrees)
    print(f"  test edges: {len(genes)}")
    print(f"  genes  absent from training (degree 0): {n_gene_missing}")
    print(f"  diseases absent from training (degree 0): {n_dis_missing}")

    # Load per-edge ranks once for every available cell (KGE seed files + baseline dumps).
    cells: dict = {}
    present = set()
    for model in KGE_MODELS:
        for regime in REGIMES:
            rr = load_reciprocal_ranks(model, regime)
            if rr is None:
                print(f"  skip {model}/{regime}: per-edge ranks not found for all seeds")
                continue
            cells[f"{model}|{regime}"] = rr
            present.add(model)
    for regime in (REGIMES if WITH_BASELINES else []):
        bl = load_baseline_reciprocal_ranks(regime)
        if not bl:
            print(f"  note {regime}: no per-edge baseline ranks (run_baselines.py not re-run?)")
        for method, rr in bl.items():
            if rr.shape[1] != len(genes):
                _die(f"{method}/{regime}: {rr.shape[1]} ranks != {len(genes)} test edges")
            cells[f"{method}|{regime}"] = rr
            present.add(method)
    if not cells:
        _die("no per-edge ranks found -- nothing to stratify")
    covered = [m for m in METHOD_ORDER if m in present]

    out: dict = {
        "meta": {
            "description": "Held-out gene-disease test MRR/Hits@10 stratified by node degree "
                           "in the R0 training graph. No retraining: recomputed from the "
                           "per-edge reciprocal ranks stored by each run. Stratified along "
                           "two axes -- the query gene (source) and the ranked disease "
                           "(target). Under the sampled-negative protocol the true disease is "
                           "ranked against type-matched negative diseases, so DISEASE degree "
                           "is the axis that exposes popularity bias.",
            "degree_source": "splits/node_degree.csv (R0 training-graph degree)",
            "n_test_edges": len(genes),
            "n_genes_absent_from_training": n_gene_missing,
            "n_diseases_absent_from_training": n_dis_missing,
            "n_quantiles": N_QUANTILES,
            "seeds": SEEDS,
            "methods_covered": covered,
            "method_class": {m: METHOD_CLASS[m] for m in covered},
            "models_covered": [m for m in covered if METHOD_CLASS[m] != "baseline"],
        },
        "by_gene_degree": stratify_axis("gene (query)", gene_deg, cells),
        "by_disease_degree": stratify_axis("disease (ranked target)", dis_deg, cells),
    }

    dest = OUT_PATH
    os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)
    with open(dest, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
    print(f"  wrote {os.path.relpath(dest, REPO_ROOT)}")

    for axis_key, axis_block in (("by_gene_degree", out["by_gene_degree"]),
                                 ("by_disease_degree", out["by_disease_degree"])):
        q = axis_block["strata"]["quartiles"]
        print(f"\n  === {axis_block['axis']} degree :: R0 quartiles (Q1 low -> Q4 high) ===")
        rng = "  ".join(f"{k}[{v['degree_min']}-{v['degree_max']}]" for k, v in q.items())
        print(f"    ranges: {rng}")
        for model in covered:
            r0 = axis_block["results"].get(f"{model}|R0")
            if not r0:
                continue
            qs = "  ".join(f"{k}={r0['quartiles'][k]['MRR_mean']:.3f}" for k in q)
            lo = r0["halves"]["low"]["MRR_mean"]
            hi = r0["halves"]["high"]["MRR_mean"]
            ratio = hi / lo if lo > 0 else float("inf")
            print(f"    {model:23s} MRR by quartile: {qs}   | high/low half = {ratio:.2f}x")


if __name__ == "__main__":
    main()
