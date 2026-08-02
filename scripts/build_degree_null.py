#!/usr/bin/env python3
"""
build_degree_null.py -- R2 degree-preserving, TYPE-preserving permutation null of
the R0 training graph, for the leakage audit (Unit 7, Machine B).

WHY
---
The leakage audit decomposes each method's R0 ranking performance into "how much is
pure node degree" vs "how much is real topological structure". This script builds the
degree null: it randomizes the training graph's wiring while holding every node's degree
exactly fixed, then re-scores the topological baselines on it under the SAME evaluation
protocol used for R0. A method whose R0 signal was really just degree collapses toward
its null here; a method that used genuine structure keeps a gap above the null.

WHAT MAKES THIS CORRECT (the gotchas)
-------------------------------------
* TYPE-preserving swaps, NOT a global rewire. A naive rewire of the whole graph would
  connect a gene to a gene where a gene->disease edge used to be -- nonsense, and it
  would change the degree *composition* seen by the relation-agnostic scorers. We rewire
  WITHIN each relation type separately (Hetionet-XSwap / Gu-et-al. style). Each relation
  subgraph is built DIRECTED (source->target), so directed double-edge swaps preserve
  every node's per-relation out- and in-degree and can never cross a relation boundary or
  flip a source/target (gene/disease) role. Total per-node degree is therefore identical.
* Preserve exactly the right thing (assert it). After swapping we assert the per-node
  degree sequence over the whole recombined graph is byte-identical to the original and
  that per-relation edge counts are unchanged. A degree null with a changed degree
  sequence is a bug (CLAUDE.md).
* Hold the evaluation FIXED for a clean permutation test. Negatives, category pools and
  the known-tail exclusions are built ONCE from R0 (lib_eval._category_pools /
  _known_tails_by_head) and the R0 test edges + seed=42 are reused for every replicate.
  The ONLY thing that changes per replicate is the scoring adjacency, so any change in MRR
  is attributable to the wiring alone.
* Sanity anchor. PreferentialAttachment reads only degree, which the null preserves, so
  its null MRR must ~= its R0 MRR. If PA moves materially, the swap changed degrees (bug);
  we assert this.

INPUT
-----
    data/processed/splits/train.csv   (source_id,relation,target_id)  -- the R0 train graph
    data/processed/splits/test.csv    -- ranking targets (via lib_eval.load_regime("R0"))

OUTPUT
------
    data/processed/splits/train_R2_degree_null_seed42.csv    first replicate; fixed name so
                                                             lib_eval's R2 regime + `run_baselines
                                                             --regimes R2` resolve to it.
    data/processed/splits/null/degree_null_rep<k>.csv        every replicate (k=0..N-1), unless
                                                             --no-write-replicates is given.
    data/processed/results/null/degree_null.json             per method: real R0 MRR, mean null
                                                             MRR, per-replicate null MRRs, and a
                                                             permutation p-value.

A replicate CSV is ~291 MB, so a large permutation distribution is disk-bound rather than
compute-bound: 200 replicates would need ~58 GB. --no-write-replicates keeps each replicate
in memory, scores it, records its per-method MRR, and discards it. The canonical rewiring is
still written, so the R2 regime file and every result derived from it are unchanged.

RUNTIME PREREQUISITES (this is Machine B's heavy job)
-----------------------------------------------------
Needs python-igraph + numpy + pandas, plus scripts/lib_eval.py and scripts/run_baselines.py
on the path and the R0 split files on disk. It CANNOT run in the no-package analysis sandbox
(no igraph); run it on the box that holds the splits. `--self-test` validates the swap
invariants on a tiny synthetic graph and needs no split files.

Usage:
    python scripts/build_degree_null.py                     # 10 replicates, seed 42
    python scripts/build_degree_null.py --replicates 20
    python scripts/build_degree_null.py --replicates 200 --no-write-replicates
    python scripts/build_degree_null.py --self-test         # fast, file-free invariant check
"""
from __future__ import annotations

import argparse
import json
import random
import time
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import igraph as ig
except ImportError as _e:  # pragma: no cover - environment guard
    raise SystemExit(
        "build_degree_null.py requires python-igraph (`pip install python-igraph`). "
        f"Import failed: {_e}")

try:
    import lib_eval
    import run_baselines
except ImportError:  # allow running from the repo root
    from scripts import lib_eval, run_baselines

try:
    from config import PROCESSED_DIR
    DEFAULT_SPLITS = Path(PROCESSED_DIR) / "splits"
    DEFAULT_RESULTS = Path(PROCESSED_DIR) / "results" / "null"
except Exception:
    DEFAULT_SPLITS = Path("data/processed/splits")
    DEFAULT_RESULTS = Path("data/processed/results/null")

# lib_eval.load_regime("R2") and `run_baselines --regimes R2` resolve to this exact
# filename via split_manifest.json, so the canonical (first) replicate is written here
# regardless of --seed.
R2_CANONICAL_NAME = "train_R2_degree_null_seed42.csv"
COLS = ["source_id", "relation", "target_id"]


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


# --------------------------------------------------------------------------- #
# Degree-preserving, type-preserving swap
# --------------------------------------------------------------------------- #
def seed_igraph(seed: int) -> None:
    """Seed igraph's RNG deterministically across python-igraph versions.

    Newer python-igraph uses its own Mersenne Twister; older versions delegate to
    Python's ``random``. We seed both routes. Bit-for-bit reproducibility across
    machines still requires a pinned igraph version (record it in requirements).
    """
    random.seed(seed)
    setter = getattr(ig, "set_random_number_generator", None)
    if setter is not None:
        try:
            setter(random)
        except Exception:
            pass
    set_seed = getattr(ig, "set_random_seed", None)
    if set_seed is not None:
        try:
            set_seed(seed)
        except Exception:
            pass


def swap_one_relation(src: np.ndarray, tgt: np.ndarray, swap_mult: int, seed: int):
    """Degree-preserving double-edge swaps within a single relation.

    ``src``/``tgt`` are aligned object arrays of node IDs for one relation (each row an
    edge source->target). Builds a DIRECTED simple igraph, runs ~``swap_mult`` * |E|
    double-edge swaps (mode="simple": no self-loops or parallel edges introduced), and
    returns the rewired (src, tgt) arrays. Asserts per-node out- and in-degree are
    identical before/after, which guarantees total degree is preserved and that no edge
    crossed the source/target (e.g. gene/disease) role -- the type-preservation guarantee.
    """
    m = len(src)
    codes, uniques = pd.factorize(np.concatenate([src, tgt]))
    src_codes = codes[:m]
    tgt_codes = codes[m:]

    if m < 2:  # nothing to swap; degree trivially preserved
        return src.copy(), tgt.copy()

    g = ig.Graph(n=len(uniques),
                 edges=list(zip(src_codes.tolist(), tgt_codes.tolist())),
                 directed=True)
    seed_igraph(seed)
    # "simple" = never introduce self-loops or parallel edges. The keyword was renamed
    # mode -> allowed_edge_types in igraph 1.0; try the new name, fall back for older ones.
    try:
        g.rewire(n=swap_mult * m, allowed_edge_types="simple")
    except TypeError:
        g.rewire(n=swap_mult * m, mode="simple")

    new = np.asarray(g.get_edgelist(), dtype=np.int64)
    assert len(new) == m, f"rewire changed edge count: {len(new)} != {m}"
    new_src_codes, new_tgt_codes = new[:, 0], new[:, 1]

    nu = len(uniques)
    assert np.array_equal(np.bincount(src_codes, minlength=nu),
                          np.bincount(new_src_codes, minlength=nu)), \
        "out-degree sequence changed by rewire"
    assert np.array_equal(np.bincount(tgt_codes, minlength=nu),
                          np.bincount(new_tgt_codes, minlength=nu)), \
        "in-degree sequence changed by rewire"

    return uniques[new_src_codes], uniques[new_tgt_codes]


def build_replicate(df: pd.DataFrame, swap_mult: int, seed: int) -> pd.DataFrame:
    """Rewire every relation of ``df`` independently and recombine into one edge list.

    Relations are processed in sorted order for determinism. Asserts, on the recombined
    graph, that (a) the per-node total-incidence degree sequence equals the original and
    (b) per-relation edge counts are unchanged -- the correctness guarantees for the null.
    """
    parts = []
    # Deterministic per-relation seed: sorted-relation index folded with the base
    # seed (NOT Python's built-in hash(), which is per-process randomized). Distinct
    # relations get independent swaps; the same (seed, relation) always reproduces.
    for idx, rel in enumerate(sorted(df["relation"].unique())):
        sub = df[df["relation"] == rel]
        rel_seed = (seed * 100_003 + idx) % (2**31 - 1)
        s, t = swap_one_relation(sub["source_id"].to_numpy(),
                                 sub["target_id"].to_numpy(),
                                 swap_mult, seed=rel_seed)
        parts.append(pd.DataFrame({"source_id": s, "relation": rel, "target_id": t}))
    null_df = pd.concat(parts, ignore_index=True)[COLS]

    orig_deg = pd.concat([df["source_id"], df["target_id"]]).value_counts().sort_index()
    null_deg = pd.concat([null_df["source_id"], null_df["target_id"]]).value_counts().sort_index()
    assert orig_deg.equals(null_deg), "degree sequence NOT preserved -- swap bug"
    assert (df["relation"].value_counts().sort_index()
            .equals(null_df["relation"].value_counts().sort_index())), \
        "per-relation edge counts changed -- swap bug"
    return null_df


# --------------------------------------------------------------------------- #
# Scoring (evaluation held fixed to R0; only the adjacency changes)
# --------------------------------------------------------------------------- #
def count_recreated_test_edges(rep_edges: np.ndarray, test_lookup: set) -> int:
    """How many held-out R0 test triples the rewiring recreates exactly.

    A degree-preserving swap is free to land a (head, relation, tail) that happens to be
    one of the held-out test triples, which puts a true positive back into the adjacency
    the null is scored on. That can only help the null, so the permutation test stays
    conservative, but the rate belongs in the methods disclosure. ``test_lookup`` is the
    set of R0 test triples, built once and reused across replicates.
    """
    rel_mask = pd.Index(rep_edges[:, 1]).isin({t[1] for t in test_lookup})
    sub = rep_edges[rel_mask]
    return sum(1 for e in zip(sub[:, 0], sub[:, 1], sub[:, 2]) if e in test_lookup)


def score_methods(rep_edges, train0, test0, hubs0, pools, known,
                  methods, n_neg, hub_cap, eval_seed):
    """Score ``methods`` on a replicate's adjacency, ranking the fixed R0 test edges.

    ``rep_edges`` is the replicate's (N,3) [source,relation,target] array (its wiring is
    the only thing that varies). pools/known/test0/eval_seed are the R0 evaluation, held
    fixed and passed straight through to rank_test_edges. Returns method -> ranking_metrics.
    """
    adj, deg = run_baselines.build_adjacency(rep_edges)
    scorers = run_baselines.make_scorers(adj, deg, hub_cap, seed=eval_seed)
    out = {}
    for name in methods:
        fn = scorers[name]
        ranks = lib_eval.rank_test_edges(
            fn, train0, test0, hubs0, False,
            n_neg=n_neg, seed=eval_seed, pools=pools, known=known)
        out[name] = lib_eval.ranking_metrics(ranks)
    return out


# --------------------------------------------------------------------------- #
# Self-test: swap invariants on a tiny synthetic multi-relation graph (no files)
# --------------------------------------------------------------------------- #
def _self_test(swap_mult=10, seed=42):
    rng = np.random.default_rng(0)
    genes = [f"HGNC:{i}" for i in range(60)]
    dis = [f"MONDO:{i:07d}" for i in range(120)]

    rows, seen = [], set()
    while len(seen) < 400:  # bipartite gene -> disease
        k = (genes[rng.integers(0, 60)], "BIOLINK:GENE_ASSOCIATED_WITH_CONDITION",
             dis[rng.integers(0, 120)])
        if k not in seen:
            seen.add(k); rows.append(k)
    while len(seen) < 600:  # symmetric-ish gene <-> gene interactions
        a, b = genes[rng.integers(0, 60)], genes[rng.integers(0, 60)]
        if a == b:
            continue
        k = (a, "BIOLINK:INTERACTS_WITH", b)
        if k not in seen:
            seen.add(k); rows.append(k)
    df = pd.DataFrame(rows, columns=COLS)

    null_df = build_replicate(df, swap_mult, seed)   # asserts degree + count invariants

    # Type preservation: the bipartite relation stays HGNC -> MONDO after swapping.
    gd = null_df[null_df["relation"].str.endswith("GENE_ASSOCIATED_WITH_CONDITION")]
    assert gd["source_id"].str.startswith("HGNC:").all(), "gene->disease source role broken"
    assert gd["target_id"].str.startswith("MONDO:").all(), "gene->disease target role broken"

    # PA proxy: undirected distinct-neighbour degree (what PreferentialAttachment reads)
    # must be ~unchanged, so PA's null MRR will track its R0 MRR.
    _, deg_o = run_baselines.build_adjacency(df[COLS].to_numpy(dtype=object))
    _, deg_n = run_baselines.build_adjacency(null_df[COLS].to_numpy(dtype=object))
    drift = max(abs(deg_o[n] - deg_n.get(n, 0)) for n in deg_o)
    assert drift <= 2, f"distinct-neighbour degree drifted by {drift} (>2) -- unexpected"

    # The wiring must actually change (this is a permutation, not a copy).
    same = df.merge(null_df, on=COLS, how="inner")
    log(f"self-test: {len(df):,} edges, {df['relation'].nunique()} relations; "
        f"{len(same):,}/{len(df):,} edges unchanged; max degree drift={drift}")
    assert len(same) < len(df), "swap did not change any edge -- rewire is a no-op"
    log("SELF-TEST PASSED: degree + per-relation counts + type roles preserved; wiring randomized.")


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--splits-dir", default=str(DEFAULT_SPLITS))
    ap.add_argument("--out-results", default=str(DEFAULT_RESULTS))
    ap.add_argument("--replicates", type=int, default=10,
                    help="number of null replicates for the permutation distribution")
    ap.add_argument("--seed", type=int, default=42, help="base seed for the swaps")
    ap.add_argument("--swap-mult", type=int, default=10,
                    help="double-edge swaps per relation = swap-mult * |E_relation|")
    ap.add_argument("--n-neg", type=int, default=50, help="negatives per test edge (match R0)")
    ap.add_argument("--hub-cap", type=int, default=2000, help="CN/AA hub cap (match R0)")
    ap.add_argument("--methods", default=",".join(run_baselines.METHOD_ORDER))
    ap.add_argument("--eval-seed", type=int, default=42,
                    help="FIXED negative-sampling seed reused for every replicate")
    ap.add_argument("--pa-tol", type=float, default=0.05,
                    help="max |PA null MRR - PA R0 MRR| before the sanity assert fails")
    ap.add_argument("--resume", action="store_true",
                    help="reuse existing replicate CSVs instead of re-running their swaps")
    ap.add_argument("--force-canonical", action="store_true",
                    help="overwrite an existing train_R2_degree_null_seed42.csv. This is a "
                         "frozen, hash-registered artifact that this script does not "
                         "reproduce bit-for-bit; regenerating it invalidates every published "
                         "R2 result and requires re-hashing split_manifest.json")
    ap.add_argument("--no-write-replicates", action="store_true",
                    help="score each replicate in memory and discard it instead of writing a "
                         "~291 MB CSV per replicate; the canonical R2 file is still written")
    ap.add_argument("--self-test", action="store_true",
                    help="validate swap invariants on a synthetic graph; needs no split files")
    args = ap.parse_args()

    if args.self_test:
        _self_test(args.swap_mult, args.seed)
        return

    methods = [m.strip() for m in args.methods.split(",") if m.strip()]
    splits_dir = Path(args.splits_dir)
    null_dir = splits_dir / "null"
    if not args.no_write_replicates:
        null_dir.mkdir(parents=True, exist_ok=True)
    out_results = Path(args.out_results)
    out_results.mkdir(parents=True, exist_ok=True)

    # -- R0: load once; build the fixed evaluation (pools/known) once -------- #
    reg0 = lib_eval.load_regime("R0", splits_dir)
    train0, test0, hubs0 = reg0
    log(f"R0 loaded: train={len(train0):,} test={len(test0):,} "
        f"(train_file={reg0.train_file}, test_file={reg0.test_file})")
    t0 = time.time()
    pools = lib_eval._category_pools(train0, hubs0)
    known = lib_eval._known_tails_by_head(train0, test0, test0[:, 0])
    log(f"built fixed category pools + known-tails from R0 ({time.time() - t0:.1f}s)")

    df = pd.DataFrame(train0, columns=COLS)  # exact R0 train graph, for swapping

    # -- real R0 baseline MRRs (same harness/seed as the null) --------------- #
    log("scoring baselines on the REAL R0 graph ...")
    real = score_methods(train0, train0, test0, hubs0, pools, known,
                          methods, args.n_neg, args.hub_cap, args.eval_seed)
    for name in methods:
        log(f"  R0 {name:23s} MRR={real[name]['MRR']:.4f} "
            f"H@10={real[name]['Hits@10']:.4f}")

    # -- replicates: swap -> (optional checkpoint CSV) -> score -------------- #
    # Built once: the R0 test triples, for the per-replicate recreation count.
    test_lookup = set(zip(test0[:, 0], test0[:, 1], test0[:, 2]))
    null_mrr = {name: [] for name in methods}          # method -> [per-rep MRR]
    null_metrics = {name: [] for name in methods}       # method -> [per-rep full metrics]
    recreated = []                                      # per-rep recreated test triples
    rep_files = []
    for k in range(args.replicates):
        rep_seed = args.seed + k
        rep_path = null_dir / f"degree_null_rep{k}.csv"
        if args.resume and rep_path.exists():
            log(f"[rep {k}] resume: reading {rep_path.name}")
            rep_edges = lib_eval._load_edges(rep_path)
            rep_files.append(str(rep_path))
        else:
            t0 = time.time()
            null_df = build_replicate(df, args.swap_mult, rep_seed)
            t_swap = time.time() - t0
            if args.no_write_replicates:
                log(f"[rep {k}] swapped in memory ({len(null_df):,} edges, {t_swap:.1f}s)")
            else:
                null_df.to_csv(rep_path, index=False)
                rep_files.append(str(rep_path))
                log(f"[rep {k}] swapped + wrote {rep_path.name} "
                    f"({len(null_df):,} edges, {t_swap:.1f}s swap, "
                    f"{time.time() - t0 - t_swap:.1f}s write)")
            rep_edges = null_df[COLS].to_numpy(dtype=object)
            del null_df
        if k == 0:  # canonical R2 file so lib_eval / run_baselines resolve R2
            canon = splits_dir / R2_CANONICAL_NAME
            # The canonical file is a FROZEN artifact: its sha256 is registered in
            # split_manifest.json and every published R2 number (Table 2's KGE column,
            # the degree-stratified analysis, the full-ranking robustness run) was
            # computed on it. This script does NOT reproduce it bit-for-bit -- igraph
            # seeds its own RNG, so rep 0 here is a statistically equivalent but
            # different rewiring. Overwriting it would silently invalidate those
            # results, so an existing file is left alone unless --force-canonical.
            if canon.exists() and not args.force_canonical:
                log(f"[rep 0] canonical {canon.name} already exists; leaving it untouched "
                    f"(pass --force-canonical to regenerate, which invalidates every "
                    f"published R2 number and the split_manifest.json hash)")
            else:
                pd.DataFrame(rep_edges, columns=COLS).to_csv(canon, index=False)
                log(f"[rep 0] wrote canonical {canon.name}")

        t0 = time.time()
        n_rec = count_recreated_test_edges(rep_edges, test_lookup)
        recreated.append(n_rec)
        log(f"[rep {k}] recreates {n_rec}/{len(test0):,} test triples "
            f"({100 * n_rec / len(test0):.2f}%, {time.time() - t0:.1f}s)")

        t0 = time.time()
        m = score_methods(rep_edges, train0, test0, hubs0, pools, known,
                          methods, args.n_neg, args.hub_cap, args.eval_seed)
        for name in methods:
            null_mrr[name].append(m[name]["MRR"])
            null_metrics[name].append(m[name])
        log(f"[rep {k}] scored ({time.time() - t0:.1f}s): "
            + "  ".join(f"{n}={m[n]['MRR']:.4f}" for n in methods))
        del rep_edges  # the replicate itself is not needed past this point

    # -- the frozen canonical graph's own recreation count -------------------- #
    # Measured on the file, not on replicate 0, because the two are different draws.
    # This is the number the manuscript quotes for the canonical rewiring.
    canon_path = splits_dir / R2_CANONICAL_NAME
    canon_recreated = None
    if canon_path.exists():
        t0 = time.time()
        canon_recreated = count_recreated_test_edges(
            lib_eval._load_edges(canon_path), test_lookup)
        log(f"canonical {R2_CANONICAL_NAME} recreates {canon_recreated}/{len(test0):,} "
            f"test triples ({100 * canon_recreated / len(test0):.2f}%, "
            f"{time.time() - t0:.1f}s)")

    # -- aggregate: permutation p-value per method --------------------------- #
    N = args.replicates
    methods_out = {}
    for name in methods:
        real_mrr = real[name]["MRR"]
        nm = np.asarray(null_mrr[name], dtype=float)
        count_ge = int(np.sum(nm >= real_mrr))
        methods_out[name] = {
            "real_R0_MRR": real_mrr,
            "real_R0_metrics": {k: real[name][k] for k in real[name]},
            "null_MRR_mean": float(nm.mean()),
            "null_MRR_sd": float(nm.std(ddof=1)) if N > 1 else 0.0,
            "null_MRR_per_rep": [float(x) for x in nm],
            "real_minus_null_mean": float(real_mrr - nm.mean()),
            # p = fraction of null replicates as good as / better than the real graph.
            "p_value": count_ge / N,
            # add-one variant (never reports exactly 0; standard permutation correction).
            "p_value_plus_one": (1 + count_ge) / (N + 1),
            "n_null_ge_real": count_ge,
            "structure_beyond_degree": bool(count_ge / N < 0.05),
        }

    # -- PA sanity anchor: null MRR must track the R0 MRR -------------------- #
    pa = methods_out.get("PreferentialAttachment")
    pa_sanity = None
    if pa is not None:
        diff = abs(pa["null_MRR_mean"] - pa["real_R0_MRR"])
        pa_sanity = {"real": pa["real_R0_MRR"], "null_mean": pa["null_MRR_mean"],
                     "abs_diff": diff, "tol": args.pa_tol, "pass": bool(diff <= args.pa_tol)}
        log(f"PA sanity: real={pa['real_R0_MRR']:.4f} null={pa['null_MRR_mean']:.4f} "
            f"|diff|={diff:.4f} (tol {args.pa_tol})")
        assert diff <= args.pa_tol, (
            f"PreferentialAttachment null MRR ({pa['null_MRR_mean']:.4f}) drifted from its "
            f"R0 MRR ({pa['real_R0_MRR']:.4f}) by {diff:.4f} > tol {args.pa_tol}. PA reads "
            "only degree, which the null preserves, so this means the swap changed degrees "
            "-- a bug. Investigate before trusting the other methods' nulls.")

    payload = {
        "created": time.strftime("%Y-%m-%d %H:%M:%S"),
        "description": "R2 degree-preserving, type-preserving permutation null of the R0 "
                       "training graph; baselines re-scored with the R0 evaluation held fixed.",
        "seed": args.seed,
        "eval_seed": args.eval_seed,
        "n_replicates": N,
        "swap_multiplier": args.swap_mult,
        "n_neg": args.n_neg,
        "hub_cap": args.hub_cap,
        "splits_dir": str(splits_dir),
        "r0_train_file": reg0.train_file,
        "r0_test_file": reg0.test_file,
        "n_train_edges": int(len(train0)),
        "n_test_edges": int(len(test0)),
        "canonical_R2_file": R2_CANONICAL_NAME,
        "replicate_files": rep_files,
        "replicates_written": not args.no_write_replicates,
        # How often a degree-preserving swap happens to put a held-out test triple back
        # into the scoring adjacency. Only ever helps the null, so the test stays
        # conservative; the manuscript discloses the distribution.
        "recreated_test_edges": {
            "per_rep": [int(x) for x in recreated],
            "mean": float(np.mean(recreated)) if recreated else 0.0,
            "min": int(np.min(recreated)) if recreated else 0,
            "max": int(np.max(recreated)) if recreated else 0,
            "mean_frac": float(np.mean(recreated) / len(test0)) if recreated else 0.0,
            "min_frac": float(np.min(recreated) / len(test0)) if recreated else 0.0,
            "max_frac": float(np.max(recreated) / len(test0)) if recreated else 0.0,
            # The frozen canonical R2 graph, measured on the file itself. It is a
            # separate draw from replicate 0, so it gets its own count.
            "canonical_file": canon_recreated,
            "n_test_edges": int(len(test0)),
        },
        "eval_protocol": "pools/known/negatives(seed) built once from R0; only the scoring "
                         "adjacency changes per replicate.",
        "degree_sequence_preserved": True,  # asserted per replicate in build_replicate
        "methods": methods_out,
        "pa_sanity": pa_sanity,
    }
    out_json = out_results / "degree_null.json"
    out_json.write_text(json.dumps(payload, indent=2))
    log(f"wrote {out_json}")

    # -- print the reportable table ----------------------------------------- #
    print("\n=== R2 degree-preserving null: real R0 vs null MRR "
          f"(N={N} replicates) ===")
    print(f"{'method':24s} {'real_R0':>8s} {'null_mean':>10s} {'null_sd':>8s} "
          f"{'real-null':>10s} {'p(null>=real)':>14s}")
    for name in methods:
        mo = methods_out[name]
        print(f"{name:24s} {mo['real_R0_MRR']:8.4f} {mo['null_MRR_mean']:10.4f} "
              f"{mo['null_MRR_sd']:8.4f} {mo['real_minus_null_mean']:10.4f} "
              f"{mo['p_value']:14.3f}")
    rec = payload["recreated_test_edges"]
    print(f"\nTest triples recreated by the rewiring: mean {rec['mean']:.1f}/{len(test0):,} "
          f"({100 * rec['mean_frac']:.2f}%), range {rec['min']}-{rec['max']} "
          f"({100 * rec['min_frac']:.2f}-{100 * rec['max_frac']:.2f}%)"
          + (f"; canonical R2 file {rec['canonical_file']} "
             f"({100 * rec['canonical_file'] / len(test0):.2f}%)."
             if rec["canonical_file"] is not None else "."))
    print("\nInterpretation: a large real-null gap (small p) = genuine structure beyond "
          "degree;\na small gap (large p) = the R0 performance was mostly degree.")
    log("DONE.")


if __name__ == "__main__":
    main()
