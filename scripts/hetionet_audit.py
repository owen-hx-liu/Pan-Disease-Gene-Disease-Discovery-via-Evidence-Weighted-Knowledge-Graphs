#!/usr/bin/env python3
"""
hetionet_audit.py -- independent-graph robustness replication of the leakage audit.

The Monarch benchmark (build_deleaked_splits.py + run_baselines.py + build_degree_null.py)
showed that topological link-prediction scores on the gene->disease task are partly real
structure and partly a pure-degree artifact: under the R2 degree-preserving null, the
neighbourhood scorers (CommonNeighbors / AdamicAdar / Jaccard) collapse a lot, while
PreferentialAttachment (which reads only degree) is unchanged. This script asks whether
that pattern is a Monarch quirk or a property of biomedical KGs by re-running the SAME
harness on a completely independent graph -- Hetionet v1.0 -- and NEVER merging the two.

Why the ID remap is the trick
-----------------------------
lib_eval calls scripts/kg_categories.category_of internally (to pick type-matched disease
negatives). Hetionet node IDs are "Type::id"; category_of understands biomedical prefixes.
Rather than fork the harness, we remap Hetionet IDs to prefixes category_of already knows:

    Gene::9021              -> NCBIGENE:9021        (category 'gene')
    Disease::DOID:263       -> DOID:263             (category 'disease')
    Compound::DB00014       -> DRUGBANK:DB00014     (category 'compound')
    <OtherKind>::<id>       -> <OtherKind>:<id>     (category 'unknown', pass-through)

The pass-through rule is a *correctness* choice, not laziness: Hetionet Symptom IDs are
MeSH ("D000006") and Side Effect IDs are UMLS ("C0000727"). Remapping those to MESH:/UMLS:
would make category_of read them as *diseases* and poison the disease negative pool. Keeping
the Hetionet kind as the prefix ("Symptom:...", "SideEffect:...") -> 'unknown' avoids that.

The task and the regimes
------------------------
Target metaedge = DaG (Disease-associates-Gene), oriented gene->disease. The SIF stores it
Disease->Gene; we flip every DaG edge to (gene, DaG, disease) so 'gene' is the head, exactly
like Monarch's HGNC->disease target. Then, reusing the SAME split logic as
build_deleaked_splits.py and the SAME harness (lib_eval.py) + scorers (run_baselines.py):

    R0 Standard      train = full remapped graph minus the held-out DaG edges
    R1 Redundancy    R0 minus every train edge that DIRECTLY connects a held-out
                     gene-disease pair in either direction (inverse/duplicate/symmetric;
                     in Hetionet these are the DdG / DuG edges on a held pair)
    R2 Degree null   degree-preserving, type-preserving permutation of the R0 graph, via
                     build_degree_null.build_replicate (Unit 7). Evaluation is held FIXED to
                     R0 (pools / known-tails / negatives from R0); only the scoring adjacency
                     changes -- so any MRR change is attributable to wiring alone.
    R3               SKIPPED. Hetionet has no cross-species orthology (no ORTHOLOGOUS_TO /
                     MODEL_OF metaedge), so the orthology-blocked regime is not defined here.

Never merged with Monarch. The claim is qualitative: does the R0->R2 pattern reproduce?

Outputs
-------
    data/processed/results/hetionet/baselines_<regime>.json   (mean +/- 95% CI over 3 seeds)
    data/processed/results/hetionet/hetionet_audit_manifest.json  (provenance + counts + hashes)
    data/external/hetionet/                                   (consumed SIF + node table; gitignored)

Usage:
    .venv/bin/python scripts/hetionet_audit.py                       # R0,R1,R2 x seeds 42,1,7
    .venv/bin/python scripts/hetionet_audit.py --regimes R0 --seeds 42   # quick smoke test
    .venv/bin/python scripts/hetionet_audit.py --self-test           # remap unit check, no data
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import time
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import lib_eval
    import run_baselines
    from build_degree_null import build_replicate
    from kg_categories import category_of
except ImportError:  # allow running from the repo root
    from scripts import lib_eval, run_baselines
    from scripts.build_degree_null import build_replicate
    from scripts.kg_categories import category_of

try:
    import igraph as _ig
    IGRAPH_VERSION = _ig.__version__
except Exception:  # pragma: no cover
    IGRAPH_VERSION = "unknown"

try:
    from config import PROCESSED_DIR, REPO_ROOT
    RESULTS_DIR = Path(PROCESSED_DIR) / "results" / "hetionet"
    EXTERNAL_DIR = Path(REPO_ROOT) / "data" / "external" / "hetionet"
except Exception:
    RESULTS_DIR = Path("data/processed/results/hetionet")
    EXTERNAL_DIR = Path("data/external/hetionet")

DATA_CHANNEL = Path("~/Downloads/sfy2-data").expanduser()
# The data channel is Windows-origin: filenames contain literal backslashes on this POSIX box.
CHANNEL_EDGES = DATA_CHANNEL / "hetionet\\hetionet-v1.0-edges.sif"
CHANNEL_NODES = DATA_CHANNEL / "hetionet\\hetionet-v1.0-nodes.tsv"

TARGET_METAEDGE = "DaG"        # Disease-associates-Gene, oriented gene->disease here
COLS = ["source_id", "relation", "target_id"]
METHOD_ORDER = run_baselines.METHOD_ORDER
METRICS = ("MRR", "Hits@1", "Hits@3", "Hits@10", "AUROC", "AUPRC")

# Monarch R0->null MRR for side-by-side eyeballing (from build_degree_null.py, Unit 7).
MONARCH_REF = {
    "Random": (0.094, 0.094),
    "CommonNeighbors": (0.441, 0.249),
    "AdamicAdar": (0.458, 0.264),
    "Jaccard": (0.398, 0.202),
    "PreferentialAttachment": (0.589, 0.591),
}


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# --------------------------------------------------------------------------- #
# ID remap  (Type::id  ->  a prefix category_of understands)
# --------------------------------------------------------------------------- #
def remap_id(node: str) -> str:
    """Remap one Hetionet 'Kind::id' node ID. See the module docstring for the rules."""
    kind, _, rest = node.partition("::")
    if kind == "Gene":
        return "NCBIGENE:" + rest
    if kind == "Disease":
        return rest                      # already 'DOID:xxxx'
    if kind == "Compound":
        return "DRUGBANK:" + rest
    # Pass-through: keep the Hetionet kind as the prefix so category_of returns 'unknown'
    # (crucial for Symptom=MeSH and SideEffect=UMLS, which must NOT read as diseases).
    return kind.replace(" ", "") + ":" + rest


# --------------------------------------------------------------------------- #
# Consume Hetionet from the data channel (copy + hash), then load + remap
# --------------------------------------------------------------------------- #
def consume(refresh: bool) -> dict:
    """Copy the Hetionet SIF/nodes out of the data channel into data/external (gitignored)."""
    EXTERNAL_DIR.mkdir(parents=True, exist_ok=True)
    dest_edges = EXTERNAL_DIR / "hetionet-v1.0-edges.sif"
    dest_nodes = EXTERNAL_DIR / "hetionet-v1.0-nodes.tsv"
    for src, dst in ((CHANNEL_EDGES, dest_edges), (CHANNEL_NODES, dest_nodes)):
        if refresh or not dst.exists():
            if not src.exists():
                raise SystemExit(f"data-channel file not found: {src}\n"
                                 "Expected the Hetionet dump in ~/Downloads/sfy2-data/.")
            log(f"consume: {src.name} -> {dst}")
            shutil.copyfile(src, dst)
    prov = {
        "edges_sif": str(dest_edges),
        "edges_sha256": sha256(dest_edges),
        "nodes_tsv": str(dest_nodes),
        "nodes_sha256": sha256(dest_nodes),
        "source_channel": str(CHANNEL_EDGES),
        "note": "No upstream Hetionet hash existed in ARTIFACT_HASHES.txt / the channel "
                "manifest (those cover only the Monarch artifacts); these sha256 are recorded "
                "here to freeze the exact Hetionet bytes this audit consumed.",
    }
    log(f"edges sha256 = {prov['edges_sha256']}")
    return prov, dest_edges


def load_remapped_edges(sif_path) -> pd.DataFrame:
    """Load the SIF, remap IDs, orient DaG gene->disease, drop self-loops + dups.

    Returns a DataFrame with columns [source_id, relation, target_id]; every DaG row is
    (gene, 'DaG', disease); every other metaedge keeps its native source->target order.
    """
    log(f"loading {sif_path}")
    raw = pd.read_csv(sif_path, sep="\t", dtype=str,
                      names=["source", "metaedge", "target"], header=0)
    n0 = len(raw)

    # Remap unique node strings once (nodes repeat heavily), then map back -- fast.
    nodes = pd.unique(pd.concat([raw["source"], raw["target"]], ignore_index=True))
    remap = {n: remap_id(n) for n in nodes}
    src = raw["source"].map(remap).to_numpy()
    tgt = raw["target"].map(remap).to_numpy()
    rel = raw["metaedge"].to_numpy()

    # Orient the DaG target metaedge gene->disease (SIF stores it disease->gene): flip.
    is_dag = rel == TARGET_METAEDGE
    source_id = np.where(is_dag, tgt, src)   # gene on the source side for DaG
    target_id = np.where(is_dag, src, tgt)   # disease on the target side for DaG

    df = pd.DataFrame({"source_id": source_id, "relation": rel, "target_id": target_id})
    df = df[df["source_id"] != df["target_id"]]                       # drop self-loops
    df = df.drop_duplicates(subset=COLS).reset_index(drop=True)
    log(f"  {n0:,} SIF rows -> {len(df):,} edges after remap/self-loop/dedup "
        f"({df['relation'].nunique()} metaedges)")

    # Guard: the flip must make DaG a gene(head)->disease(tail) relation.
    dag = df[df["relation"] == TARGET_METAEDGE]
    assert (dag["source_id"].str.startswith("NCBIGENE:")).all(), "DaG head is not a gene"
    assert dag["target_id"].map(lambda t: category_of(t) == "disease").all(), \
        "DaG tail is not category 'disease'"
    log(f"  DaG target edges (gene->disease): {len(dag):,}")
    return df


# --------------------------------------------------------------------------- #
# Splits: 80/10/10 on DaG (same logic as build_deleaked_splits.py), then R1, R2
# --------------------------------------------------------------------------- #
def _triple_key(a, r, b):
    return a + "\x1f" + r + "\x1f" + b


def _pair_key(a, b):
    return (a + "\x1f" + b) if a <= b else (b + "\x1f" + a)   # unordered {a,b}


def build_splits(df: pd.DataFrame, seed: int, valid_frac: float, test_frac: float):
    """Replicate build_deleaked_splits' split logic on Hetionet's DaG task."""
    targets = df[df["relation"] == TARGET_METAEDGE].reset_index(drop=True)
    rng = np.random.default_rng(seed)
    idx = rng.permutation(len(targets))
    n_test = int(len(targets) * test_frac)
    n_val = int(len(targets) * valid_frac)
    test = targets.iloc[idx[:n_test]].reset_index(drop=True)
    valid = targets.iloc[idx[n_test:n_test + n_val]].reset_index(drop=True)
    log(f"split: targets={len(targets):,} -> train_targets={len(targets) - n_test - n_val:,} "
        f"valid={len(valid):,} test={len(test):,} (seed={seed})")

    held = pd.concat([test, valid], ignore_index=True)
    held_keys = set(_triple_key(a, r, b) for a, r, b in
                    held[COLS].itertuples(index=False, name=None))
    # R0 train = every edge that is NOT a held-out target edge.
    all_keys = (df["source_id"] + "\x1f" + df["relation"] + "\x1f" + df["target_id"])
    keep = ~all_keys.isin(held_keys)
    train0 = df[keep].reset_index(drop=True)
    assert len(train0) == len(df) - len(held_keys), "held edge count mismatch in R0 train"
    log(f"R0 train: {len(train0):,} edges (removed {len(held_keys):,} held-out DaG edges)")

    # R1 redundancy: drop every R0-train edge whose endpoints are a held pair (either
    # direction, any relation). The held-out DaG edges are already gone, so this removes
    # only *other* metaedges on a held gene-disease pair (Hetionet: DdG / DuG).
    held_pairs = set(_pair_key(a, b) for a, _r, b in
                     held[["source_id", "relation", "target_id"]].itertuples(index=False, name=None))
    # Vectorized unordered-pair key over the whole R0 train graph (avoids a row-wise apply).
    a = train0["source_id"].to_numpy()
    b = train0["target_id"].to_numpy()
    lo = np.where(a <= b, a, b)
    hi = np.where(a <= b, b, a)
    tp = pd.Series(lo, dtype=object) + "\x1f" + pd.Series(hi, dtype=object)
    r1_keep = ~tp.isin(held_pairs)
    train1 = train0[r1_keep].reset_index(drop=True)
    log(f"R1 train: {len(train1):,} edges (removed {len(train0) - len(train1):,} direct "
        f"held-pair edges)")

    test_arr = test[COLS].to_numpy(dtype=object)
    return train0, train1, test, valid, test_arr, len(targets), len(held_keys)


# --------------------------------------------------------------------------- #
# Evaluation (one regime): 5 scorers x N seeds, mean + bootstrap 95% CI
# --------------------------------------------------------------------------- #
def eval_regime(regime, adj_edges, test_arr, pools, known, seeds, n_neg, hub_cap):
    """Score the 5 baselines on ``adj_edges`` against the fixed ``test_arr`` ranking targets.

    ``pools`` / ``known`` are the negative-sampling context. For R0/R1 they are built from
    that regime's own train graph; for R2 they are held fixed to R0 (degree-null protocol),
    so only the scoring adjacency differs from R0.
    """
    t0 = time.time()
    adj, deg = run_baselines.build_adjacency(adj_edges)
    log(f"  [{regime}] adjacency: {len(adj):,} nodes ({time.time() - t0:.1f}s)")

    results = {m: {k: [] for k in METRICS} for m in METHOD_ORDER}
    per_edge_rr = {}
    for seed in seeds:
        scorers = run_baselines.make_scorers(adj, deg, hub_cap, seed)
        for name in METHOD_ORDER:
            t0 = time.time()
            ranks, pos, neg = lib_eval.rank_test_edges(
                scorers[name], adj_edges, test_arr, frozenset(), False,
                n_neg=n_neg, seed=seed, return_scores=True, pools=pools, known=known)
            rm = lib_eval.ranking_metrics(ranks)
            cm = lib_eval.classification_metrics(pos, neg)
            for k in ("MRR", "Hits@1", "Hits@3", "Hits@10"):
                results[name][k].append(rm[k])
            results[name]["AUROC"].append(cm["AUROC"])
            results[name]["AUPRC"].append(cm["AUPRC"])
            if seed == seeds[0]:
                per_edge_rr[name] = 1.0 / ranks
            log(f"  [{regime}] seed={seed} {name:23s} MRR={rm['MRR']:.4f} "
                f"H@10={rm['Hits@10']:.4f} AUROC={cm['AUROC']:.4f} ({time.time() - t0:.1f}s)")

    summary = {}
    for name in METHOD_ORDER:
        summary[name] = {}
        for metric, vals in results[name].items():
            v = np.asarray(vals, dtype=float)
            summary[name][metric] = {
                "mean": float(np.mean(v)),
                "sd": float(np.std(v, ddof=1)) if len(v) > 1 else 0.0,
                "per_seed": [float(x) for x in v],
            }
        # 95% CI on per-edge reciprocal ranks at the first seed (same as run_baselines).
        summary[name][f"MRR_ci_seed{seeds[0]}"] = lib_eval.bootstrap_ci(per_edge_rr[name])

    paired = {}
    if "AdamicAdar" in per_edge_rr:
        base = per_edge_rr["AdamicAdar"]
        for name, rr in per_edge_rr.items():
            if name != "AdamicAdar":
                paired[f"AdamicAdar_vs_{name}"] = lib_eval.paired_bootstrap_pvalue(base, rr)
    return summary, paired


def write_regime_json(out_dir, regime, summary, paired, meta):
    payload = {
        "graph": "hetionet-v1.0",
        "regime": regime,
        "task": "DaG (Disease-associates-Gene), oriented gene->disease",
        "never_merged_with": "monarch (edges_clean_integrated.csv)",
        **meta,
        "methods": summary,
        f"paired_vs_AdamicAdar_seed{meta['seeds'][0]}": paired,
        "created": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    p = out_dir / f"baselines_{regime}.json"
    p.write_text(json.dumps(payload, indent=2))
    log(f"wrote {p}")


# --------------------------------------------------------------------------- #
# Self-test: remap rules (no data files needed)
# --------------------------------------------------------------------------- #
def _self_test():
    cases = {
        "Gene::9021": ("NCBIGENE:9021", "gene"),
        "Disease::DOID:263": ("DOID:263", "disease"),
        "Compound::DB00014": ("DRUGBANK:DB00014", "compound"),
        "Anatomy::UBERON:0000002": ("Anatomy:UBERON:0000002", "unknown"),
        "Symptom::D000006": ("Symptom:D000006", "unknown"),      # MeSH must NOT read as disease
        "Side Effect::C0000727": ("SideEffect:C0000727", "unknown"),  # UMLS must NOT read as disease
        "Biological Process::GO:0071357": ("BiologicalProcess:GO:0071357", "unknown"),
    }
    for src, (want_id, want_cat) in cases.items():
        got_id = remap_id(src)
        got_cat = category_of(got_id)
        assert got_id == want_id, f"remap {src!r} -> {got_id!r} != {want_id!r}"
        assert got_cat == want_cat, f"category_of({got_id!r}) = {got_cat!r} != {want_cat!r}"
        log(f"  OK  {src:34s} -> {got_id:32s} [{got_cat}]")
    # The critical anti-pollution guarantees, stated as asserts:
    assert category_of(remap_id("Symptom::D000006")) != "disease", "Symptom leaked into disease pool"
    assert category_of(remap_id("Side Effect::C0000727")) != "disease", "SideEffect leaked into disease pool"
    log("SELF-TEST PASSED: remap keeps genes/diseases/compounds mapped and non-tail kinds 'unknown'.")


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regimes", default="R0,R1,R2",
                    help="comma-separated regimes to run (R3 is intentionally undefined)")
    ap.add_argument("--seeds", default="42,1,7", help="negative-sampling seeds")
    ap.add_argument("--n-neg", type=int, default=50)
    ap.add_argument("--hub-cap", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=42, help="split + degree-null base seed")
    ap.add_argument("--valid-frac", type=float, default=0.10)
    ap.add_argument("--test-frac", type=float, default=0.10)
    ap.add_argument("--swap-mult", type=int, default=10, help="R2 double-edge swaps per |E_rel|")
    ap.add_argument("--pa-tol", type=float, default=0.05,
                    help="max |PA R2 MRR - PA R0 MRR| before the degree-null sanity assert fails")
    ap.add_argument("--refresh", action="store_true", help="re-copy Hetionet from the data channel")
    ap.add_argument("--write-splits", action="store_true",
                    help="also dump the (large) R0/R1/R2 train CSVs for reproducibility")
    ap.add_argument("--out", default=str(RESULTS_DIR))
    ap.add_argument("--self-test", action="store_true", help="remap unit check; no data needed")
    args = ap.parse_args()

    if args.self_test:
        _self_test()
        return

    regimes = [r.strip().upper() for r in args.regimes.split(",") if r.strip()]
    seeds = [int(s) for s in args.seeds.split(",") if s.strip()]
    if "R3" in regimes:
        log("NOTE: R3 (orthology-blocked) is undefined for Hetionet (no cross-species "
            "orthology metaedge); dropping it.")
        regimes = [r for r in regimes if r != "R3"]
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    # -- consume + load + remap --------------------------------------------- #
    prov, sif_path = consume(args.refresh)
    df = load_remapped_edges(sif_path)

    # -- splits -------------------------------------------------------------- #
    train0, train1, test, valid, test_arr, n_targets, n_held = build_splits(
        df, args.seed, args.valid_frac, args.test_frac)
    train0_arr = train0[COLS].to_numpy(dtype=object)
    train1_arr = train1[COLS].to_numpy(dtype=object)

    # -- fixed R0 evaluation context (also reused, unchanged, for the R2 null) - #
    pools0 = lib_eval._category_pools(train0_arr, frozenset())
    known0 = lib_eval._known_tails_by_head(train0_arr, test_arr, test_arr[:, 0])
    n_disease = len(pools0.get("disease", []))
    log(f"disease negative pool: {n_disease:,} nodes (n_neg={args.n_neg})")

    # small provenance splits are always written; big train graphs only with --write-splits
    splits_dir = EXTERNAL_DIR / "splits"
    splits_dir.mkdir(parents=True, exist_ok=True)
    test[COLS].to_csv(splits_dir / "test.csv", index=False)
    valid[COLS].to_csv(splits_dir / "valid.csv", index=False)
    if args.write_splits:
        train0[COLS].to_csv(splits_dir / "train.csv", index=False)
        train1[COLS].to_csv(splits_dir / "train_R1_redundancy.csv", index=False)

    # -- build R2 degree null once (if requested) --------------------------- #
    train2_arr = None
    if "R2" in regimes:
        t0 = time.time()
        log(f"building R2 degree null (build_replicate, swap_mult={args.swap_mult}, "
            f"seed={args.seed}) over {len(train0):,} edges ...")
        null_df = build_replicate(train0[COLS], swap_mult=args.swap_mult, seed=args.seed)
        train2_arr = null_df[COLS].to_numpy(dtype=object)
        log(f"  R2 null built ({time.time() - t0:.1f}s); degree + per-relation counts asserted")
        if args.write_splits:
            null_df.to_csv(splits_dir / "train_R2_degree_null_seed42.csv", index=False)

    # -- run each regime ----------------------------------------------------- #
    common_meta = {
        "seeds": seeds, "n_neg": args.n_neg, "hub_cap": args.hub_cap,
        "split_seed": args.seed, "valid_frac": args.valid_frac, "test_frac": args.test_frac,
        "n_dag_target_edges": int(n_targets), "n_test_edges": int(len(test)),
        "n_valid_edges": int(len(valid)), "n_disease_pool": int(n_disease),
    }
    mrr_table = {}   # regime -> {method -> mean MRR}
    summaries = {}
    for regime in regimes:
        log(f"================= regime {regime} =================")
        if regime == "R0":
            adj_edges, pools, known = train0_arr, pools0, known0
            meta = {**common_meta, "n_train_edges": int(len(train0)),
                    "regime_desc": "standard: full remapped graph minus held-out DaG edges"}
        elif regime == "R1":
            pools = lib_eval._category_pools(train1_arr, frozenset())
            known = lib_eval._known_tails_by_head(train1_arr, test_arr, test_arr[:, 0])
            adj_edges = train1_arr
            meta = {**common_meta, "n_train_edges": int(len(train1)),
                    "removed_for_R1": int(len(train0) - len(train1)),
                    "regime_desc": "redundancy-controlled: R0 minus direct held-pair edges"}
        elif regime == "R2":
            # Degree-null protocol: adjacency from the null, evaluation FIXED to R0.
            adj_edges, pools, known = train2_arr, pools0, known0
            meta = {**common_meta, "n_train_edges": int(len(train2_arr)),
                    "swap_mult": args.swap_mult, "null_seed": args.seed,
                    "igraph_version": IGRAPH_VERSION,
                    "eval_fixed_to": "R0 (pools/known/negatives); only adjacency changes",
                    "regime_desc": "degree-preserving, type-preserving permutation null"}
        else:
            log(f"SKIP unknown regime {regime}")
            continue

        summary, paired = eval_regime(regime, adj_edges, test_arr, pools, known,
                                      seeds, args.n_neg, args.hub_cap)
        write_regime_json(out_dir, regime, summary, paired, meta)
        summaries[regime] = summary
        mrr_table[regime] = {m: summary[m]["MRR"]["mean"] for m in METHOD_ORDER}

    # -- degree-null sanity anchor: PA (degree only) must not move R0 -> R2 --- #
    pa_sanity = None
    if "R0" in mrr_table and "R2" in mrr_table:
        pa0 = mrr_table["R0"]["PreferentialAttachment"]
        pa2 = mrr_table["R2"]["PreferentialAttachment"]
        diff = abs(pa2 - pa0)
        pa_sanity = {"R0": pa0, "R2": pa2, "abs_diff": diff, "tol": args.pa_tol,
                     "pass": bool(diff <= args.pa_tol)}
        log(f"PA sanity: R0={pa0:.4f} R2={pa2:.4f} |diff|={diff:.4f} (tol {args.pa_tol})")
        assert diff <= args.pa_tol, (
            f"PreferentialAttachment MRR moved {diff:.4f} from R0 to R2 (> {args.pa_tol}). "
            "PA reads only degree, which the null preserves, so this means the swap changed "
            "degrees -- a bug. Investigate before trusting the other nulls.")

    # -- manifest ------------------------------------------------------------ #
    manifest = {
        "created": time.strftime("%Y-%m-%d %H:%M:%S"),
        "script": "scripts/hetionet_audit.py",
        "graph": "hetionet-v1.0",
        "purpose": "independent-graph replication of the Monarch leakage audit; NEVER merged "
                   "with Monarch",
        "provenance": prov,
        "remap_rules": {"Gene::N": "NCBIGENE:N", "Disease::DOID:x": "DOID:x",
                        "Compound::DB": "DRUGBANK:DB", "<Kind>::id": "<Kind>:id (-> 'unknown')"},
        "target_metaedge": TARGET_METAEDGE,
        "orientation": "gene->disease (SIF stores DaG disease->gene; flipped)",
        "regimes_run": regimes,
        "R3": "SKIPPED -- Hetionet has no cross-species orthology metaedge",
        "counts": common_meta,
        "pa_sanity": pa_sanity,
        "results": {r: f"baselines_{r}.json" for r in mrr_table},
    }
    (out_dir / "hetionet_audit_manifest.json").write_text(json.dumps(manifest, indent=2))
    log(f"wrote {out_dir / 'hetionet_audit_manifest.json'}")

    # -- printed MRR table: Hetionet regimes + Monarch reference for eyeballing - #
    print("\n=== Hetionet DaG (gene->disease) -- mean MRR over "
          f"{len(seeds)} seeds {tuple(seeds)} ===")
    header = f"{'method':24s}" + "".join(f"{r:>10s}" for r in regimes)
    if "R0" in mrr_table and "R2" in mrr_table:
        header += f"{'R0-R2':>10s}"
    header += f"{'  | Monarch R0->R2':>22s}"
    print(header)
    for m in METHOD_ORDER:
        row = f"{m:24s}" + "".join(f"{mrr_table[r][m]:10.4f}" for r in regimes)
        if "R0" in mrr_table and "R2" in mrr_table:
            row += f"{mrr_table['R0'][m] - mrr_table['R2'][m]:10.4f}"
        mr0, mr2 = MONARCH_REF.get(m, (float('nan'), float('nan')))
        row += f"   {mr0:6.3f}->{mr2:.3f}"
        print(row)
    print("\nClaim: the R0->R2 pattern reproduces qualitatively if the neighbourhood scorers "
          "(CommonNeighbors/AdamicAdar/Jaccard) drop from R0 to R2 while PreferentialAttachment "
          "stays put (degree-only) and Random sits at chance.")
    log("DONE.")


if __name__ == "__main__":
    main()
