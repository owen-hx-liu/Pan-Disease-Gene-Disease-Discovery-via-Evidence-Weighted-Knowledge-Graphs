#!/usr/bin/env python3
"""
build_deleaked_splits.py -- THE BACKBONE of the leakage-audit benchmark.

Builds the fixed train/valid/test splits and the leakage-controlled graph variants that
EVERY downstream method (baselines on Machine B, KGE on Machine A) must consume, so
all methods are evaluated on byte-identical held-out edges. See plan_audit/PAPER_SCOPE.md.

Prediction task: human gene -> disease
    head  = HGNC:*                     (human gene)
    rel   = BIOLINK:GENE_ASSOCIATED_WITH_CONDITION   (configurable)
    tail  = a disease node (category_of == 'disease', e.g. MONDO:*)

We hold out a random 10%/10% of these target edges as valid/test; everything else
(all non-target edges + the remaining 80% of target edges) is the training graph.

Four AUDIT REGIMES are supported by the files this script emits (the regimes themselves
are applied by the eval harness, lib_eval.py). Under the leakage-audit plan the regimes
are graded controls -- each removes one leakage source:

    R0 Standard              train = train.csv                     test = test.csv
    R1 Redundancy-controlled train = train_R1_redundancy.csv       test = test.csv
    R2 Degree-controlled     train = train_R2_degree_null_seed42.csv test = test.csv   <- built by Unit 7
    R3 Orthology-blocked     train = train_R3_orthology_blocked.csv test = test.csv

    R1 removes any *direct* edge in train between a held-out gene and its held-out
       disease (either direction, any relation): the exact reverse of the target, a
       duplicate association under another relation name, or a symmetric restatement --
       all trivially leak the held-out pair (OpenBioLink-style redundancy control).
    R2 is a degree-preserving permutation null; that file is NOT written here -- it is
       produced by Unit 7 (build_degree_null.py) from train.csv + node_degree.csv.
    R3 removes the cross-species bridges (ORTHOLOGOUS_TO, MODEL_OF by default) so a human
       gene->disease edge can no longer be recovered through a mouse/fish ortholog's
       phenotype profile. --strict-ortho also drops non-human gene HAS_PHENOTYPE edges.

Optional test cleaning (not a regime): if gene/disease labels are supplied, we also emit
test_tautology_filtered.csv = test.csv minus eponymous/obsolete pairs, usable as a
robustness check across all regimes.

Outputs (data/processed/splits/ by default):
    train.csv                        R0 training graph  (source_id,relation,target_id)
    valid.csv                        held-out target edges for KGE early stopping
    test.csv                         held-out target edges  = the ranking targets
    train_R1_redundancy.csv          R1 training graph (direct held-pair edges removed)
    train_R3_orthology_blocked.csv   R3 training graph (orthology/model-organism removed)
    test_tautology_filtered.csv      optional cleaned test (needs labels; else == test.csv)
    hub_nodes.txt                    nodes with degree > --hub-degree (graph stats / diagnostics)
    node_degree.csv                  node,degree over the R0 training graph (feeds the R2 null)
    split_manifest.json              seed, counts, regime defs, sha256 of every output

Deterministic: fixed --seed (default 42). Pure pandas/numpy -- no networkx/torch, so
Machine C (the coordinator laptop) can run it without the heavy stack.

Usage:
    python scripts/build_deleaked_splits.py            # full run, writes everything
    python scripts/build_deleaked_splits.py --dry-run  # compute + print counts, write nothing
    python scripts/build_deleaked_splits.py \
        --gene-symbols data/registry/hgnc_complete_set.txt \
        --disease-labels data/registry/mondo_labels.tsv    # enables tautology test cleaning
"""
import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from kg_categories import category_of
except ImportError:  # allow running from repo root
    from scripts.kg_categories import category_of

try:
    from config import PROCESSED_DIR
    DEFAULT_OUT = str(Path(PROCESSED_DIR) / "splits")
except Exception:
    DEFAULT_OUT = "data/processed/splits"


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def rel_suffix(r):
    """'BIOLINK:GENE_ASSOCIATED_WITH_CONDITION' -> 'GENE_ASSOCIATED_WITH_CONDITION'."""
    return r.split(":", 1)[1] if ":" in r else r


def _pair_keys(a, b):
    """Canonical, order-independent node-pair keys for two ID columns.

    _pair_keys(['g','d'], ['d','g']) -> both map to the same 'd\\x01g' key, so an edge
    and its reverse collide. Used to detect direct edges between held-out endpoints.
    """
    a = np.asarray(a, dtype=str)
    b = np.asarray(b, dtype=str)
    lo = np.where(a <= b, a, b)
    hi = np.where(a <= b, b, a)
    return np.char.add(np.char.add(lo, "\x01"), hi)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# --------------------------------------------------------------------------- #
# Load
# --------------------------------------------------------------------------- #
def load_edges(path):
    log(f"loading {path}")
    df = pd.read_csv(path, usecols=["source_id", "relation", "target_id"], dtype=str)
    n0 = len(df)
    df = df.dropna(subset=["source_id", "relation", "target_id"])
    df = df[df["source_id"] != df["target_id"]]                 # drop self-loops
    df = df.drop_duplicates(subset=["source_id", "relation", "target_id"])
    log(f"  {n0:,} rows -> {len(df):,} after dropna/self-loop/dedup")
    return df.reset_index(drop=True)


# --------------------------------------------------------------------------- #
# Target-edge identification + split
# --------------------------------------------------------------------------- #
def target_mask(df, target_relations, head_prefix):
    rels = {r.upper() for r in target_relations}
    rel_ok = df["relation"].map(lambda r: rel_suffix(r).upper() in rels)
    head_ok = df["source_id"].str.split(":", n=1).str[0].str.upper() == head_prefix.upper()
    tail_ok = df["target_id"].map(lambda t: category_of(t) == "disease")
    return rel_ok & head_ok & tail_ok


def eponymous(gene_symbol, disease_label):
    """True if the disease is (near-)defined by the gene -- tautological to 'predict'."""
    if not gene_symbol or not disease_label:
        return False
    g = gene_symbol.strip().upper()
    d = disease_label.strip().upper()
    if len(g) < 3:                       # too-short symbols cause false positives
        return False
    return g in d or f"{g}-RELATED" in d or f"{g} RELATED" in d


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--edges", default="data/processed/edges_clean_integrated.csv")
    ap.add_argument("--out", default=DEFAULT_OUT)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--valid-frac", type=float, default=0.10)
    ap.add_argument("--test-frac", type=float, default=0.10)
    ap.add_argument("--target-relations", default="GENE_ASSOCIATED_WITH_CONDITION",
                    help="comma-separated relation suffixes that define the task")
    ap.add_argument("--head-prefix", default="HGNC",
                    help="restrict the gene side to this (human) namespace")
    ap.add_argument("--block-relations", default="ORTHOLOGOUS_TO,MODEL_OF",
                    help="relation suffixes removed in the R3 orthology-blocked graph")
    ap.add_argument("--strict-ortho", action="store_true",
                    help="also drop non-human gene HAS_PHENOTYPE edges in R3")
    ap.add_argument("--hub-degree", type=int, default=2000,
                    help="nodes with train-graph degree above this are hubs (diagnostics)")
    ap.add_argument("--gene-symbols", default=None,
                    help="HGNC complete-set TSV (hgnc_id, symbol) -> enables tautology filter")
    ap.add_argument("--disease-labels", default=None,
                    help="TSV 'mondo_id<TAB>label' -> enables tautology/obsolete filter")
    ap.add_argument("--dry-run", action="store_true",
                    help="compute and print counts but write no files")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    target_relations = [s.strip() for s in args.target_relations.split(",") if s.strip()]
    block_relations = {s.strip().upper() for s in args.block_relations.split(",") if s.strip()}
    out = Path(args.out)

    df = load_edges(args.edges)

    # -- identify target (human gene -> disease) edges ----------------------- #
    tmask = target_mask(df, target_relations, args.head_prefix)
    targets = df[tmask].reset_index(drop=True)
    log(f"target edges ({'+'.join(target_relations)}, {args.head_prefix}->disease): "
        f"{len(targets):,}")
    if len(targets) < 100:
        raise SystemExit("ERROR: <100 target edges found -- check --target-relations / "
                         "--head-prefix against the actual relation names in the graph.")

    # -- split target edges 80/10/10 ---------------------------------------- #
    idx = rng.permutation(len(targets))
    n_test = int(len(targets) * args.test_frac)
    n_val = int(len(targets) * args.valid_frac)
    test_i = idx[:n_test]
    val_i = idx[n_test:n_test + n_val]
    test = targets.iloc[test_i].reset_index(drop=True)
    valid = targets.iloc[val_i].reset_index(drop=True)
    held_keys = set(map(tuple, pd.concat([test, valid])[["source_id", "relation", "target_id"]]
                        .itertuples(index=False, name=None)))

    # train = every edge that is NOT a held-out target edge
    all_keys = df[["source_id", "relation", "target_id"]].itertuples(index=False, name=None)
    keep = np.fromiter((k not in held_keys for k in all_keys), dtype=bool, count=len(df))
    train = df[keep].reset_index(drop=True)
    log(f"split: train={len(train):,}  valid={len(valid):,}  test={len(test):,}")

    # sanity: held-out edges must not appear in train
    assert len(set(map(tuple, train[["source_id", "relation", "target_id"]]
                       .itertuples(index=False, name=None))) & held_keys) == 0, \
        "LEAK: held-out edge found in train"

    # -- R1 redundancy-controlled training graph ---------------------------- #
    # Drop any DIRECT train edge between a held-out gene and its held-out disease
    # (either direction, any relation): reverse-of-target, duplicate association, or
    # symmetric restatement -- all trivially leak the held-out pair.
    held_pairs = pd.concat([test, valid], ignore_index=True)
    held_pair_keys = set(_pair_keys(held_pairs["source_id"], held_pairs["target_id"]).tolist())
    train_pair_keys = _pair_keys(train["source_id"], train["target_id"])
    redundant = np.isin(train_pair_keys, list(held_pair_keys))
    train_r1 = train[~redundant].reset_index(drop=True)
    log(f"R1 redundancy-controlled graph: {len(train_r1):,} edges "
        f"(removed {int(redundant.sum()):,} direct held-pair edges, any relation/direction)")

    # -- R3 orthology-blocked training graph (Monarch-only leakage) --------- #
    rel_suf = train["relation"].map(lambda r: rel_suffix(r).upper())
    block = rel_suf.isin(block_relations)
    if args.strict_ortho:
        nonhuman_gene = (train["source_id"].str.split(":", n=1).str[0].str.upper() != "HGNC") \
            & train["source_id"].map(lambda s: category_of(s) == "gene")
        block = block | (rel_suf.eq("HAS_PHENOTYPE") & nonhuman_gene)
    train_r3 = train[~block].reset_index(drop=True)
    log(f"R3 orthology-blocked graph: {len(train_r3):,} edges "
        f"(removed {int(block.sum()):,}: {sorted(block_relations)}"
        f"{' + non-human HAS_PHENOTYPE' if args.strict_ortho else ''})")

    # -- degrees + hubs over the R0 training graph (feeds the R2 null) ------- #
    deg = (pd.concat([train["source_id"], train["target_id"]])
           .value_counts())
    hubs = deg[deg > args.hub_degree]
    log(f"degree: {len(deg):,} nodes; {len(hubs):,} hubs (deg > {args.hub_degree})")

    # -- optional tautology-cleaned test set (needs labels; not a regime) ---- #
    taut_note = "labels not provided -- test_tautology_filtered == test (filter SKIPPED)"
    test_taut = test.copy()
    n_epo = n_obs = 0
    if args.gene_symbols and args.disease_labels:
        sym = pd.read_csv(args.gene_symbols, sep="\t", dtype=str,
                          usecols=lambda c: c in ("hgnc_id", "symbol"))
        sym_map = dict(zip(sym["hgnc_id"], sym["symbol"]))
        lab = pd.read_csv(args.disease_labels, sep="\t", dtype=str, header=None,
                          names=["mondo_id", "label"])
        lab_map = dict(zip(lab["mondo_id"], lab["label"]))
        g_sym = test["source_id"].map(sym_map.get)
        d_lab = test["target_id"].map(lab_map.get)
        obs = d_lab.fillna("").str.upper().str.startswith("OBSOLETE")
        epo = pd.Series([eponymous(a, b) for a, b in zip(g_sym, d_lab)], index=test.index)
        n_epo, n_obs = int(epo.sum()), int(obs.sum())
        test_taut = test[~(epo | obs)].reset_index(drop=True)
        taut_note = f"removed {n_epo} eponymous + {n_obs} obsolete -> {len(test_taut)} kept"
    log(f"tautology-cleaned test: {taut_note}")

    # -- manifest ----------------------------------------------------------- #
    manifest = {
        "created": time.strftime("%Y-%m-%d %H:%M:%S"),
        "seed": args.seed,
        "plan": "leakage-audit (R0 standard / R1 redundancy / R2 degree-null / R3 orthology)",
        "edges_source": args.edges,
        "task": {
            "head_prefix": args.head_prefix,
            "target_relations": target_relations,
            "tail_category": "disease",
            "n_target_edges": len(targets),
        },
        "counts": {
            "train": len(train), "valid": len(valid), "test": len(test),
            "train_R1_redundancy": len(train_r1),
            "removed_for_R1": int(redundant.sum()),
            "train_R3_orthology_blocked": len(train_r3),
            "removed_for_R3": int(block.sum()),
            "hub_nodes": len(hubs),
            "test_tautology_filtered": len(test_taut),
            "removed_eponymous": n_epo, "removed_obsolete": n_obs,
        },
        "regimes": {
            "R0_standard": {
                "train": "train.csv", "test": "test.csv", "hub_filter": False,
                "desc": "random split, full training graph (inflated reference)",
            },
            "R1_redundancy_controlled": {
                "train": "train_R1_redundancy.csv", "test": "test.csv", "hub_filter": False,
                "desc": "remove direct held-pair edges (inverse/duplicate/symmetric)",
            },
            "R2_degree_controlled": {
                "train": "train_R2_degree_null_seed42.csv", "test": "test.csv",
                "hub_filter": False,
                "desc": "degree-preserving permutation null (attributes perf to degree)",
                "produced_by": "Unit 7 build_degree_null.py -- NOT written by this script",
            },
            "R3_orthology_blocked": {
                "train": "train_R3_orthology_blocked.csv", "test": "test.csv",
                "hub_filter": False,
                "desc": "remove ORTHOLOGOUS_TO + MODEL_OF (Monarch-only orthology leakage)",
            },
        },
        "params": {
            "block_relations": sorted(block_relations),
            "strict_ortho": args.strict_ortho,
            "hub_degree": args.hub_degree,
            "valid_frac": args.valid_frac, "test_frac": args.test_frac,
        },
    }

    if args.dry_run:
        log("DRY RUN -- writing nothing. Manifest preview:")
        print(json.dumps(manifest, indent=2))
        return

    # -- write -------------------------------------------------------------- #
    out.mkdir(parents=True, exist_ok=True)
    writes = {
        "train.csv": train,
        "valid.csv": valid,
        "test.csv": test,
        "train_R1_redundancy.csv": train_r1,
        "train_R3_orthology_blocked.csv": train_r3,
        "test_tautology_filtered.csv": test_taut,
    }
    hashes = {}
    for name, frame in writes.items():
        p = out / name
        frame.to_csv(p, index=False)
        hashes[name] = sha256(p)
        log(f"wrote {name}: {len(frame):,} rows")
    hubs.rename("degree").to_csv(out / "hub_nodes.txt", header=["degree"])
    hashes["hub_nodes.txt"] = sha256(out / "hub_nodes.txt")
    deg.rename("degree").to_csv(out / "node_degree.csv", header=["degree"])
    hashes["node_degree.csv"] = sha256(out / "node_degree.csv")

    manifest["sha256"] = hashes
    with open(out / "split_manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)
    log("wrote split_manifest.json")
    log("DONE. R2 degree-null graph is built later by Unit 7. Add these sha256 to "
        "ARTIFACT_HASHES.txt and upload the split files to the shared data channel.")
    print(json.dumps(manifest["counts"], indent=2))


if __name__ == "__main__":
    main()
