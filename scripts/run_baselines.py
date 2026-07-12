#!/usr/bin/env python3
"""
run_baselines.py -- topological link-prediction baselines across the audit regimes.

The primary evidence for the leakage audit: five cheap topological scorers evaluated
under the SAME ranking harness (scripts/lib_eval.py) that the KGE models use, on the
byte-identical held-out gene->disease edges, under each de-leaking regime. The R0->R3
drop in MRR is the paper's headline (Figure 2).

Scorers (score_fn(gene, disease) -> float, higher = more likely a true edge):
    Random                  seeded uniform noise (negative control; must sit at chance)
    CommonNeighbors         |N(g) ∩ N(d)|, generic hubs (deg > --hub-cap) skipped
    AdamicAdar              sum_{w in N(g)∩N(d), 1<deg(w)<=hub-cap} 1/log(deg w)
    Jaccard                 |N(g) ∩ N(d)| / |N(g) ∪ N(d)|
    PreferentialAttachment  deg(g) * deg(d)

Adjacency is undirected and relation-agnostic (matches predict_adamic_adar.py). The
hub-cap (default 2000) reproduces the existing benchmark's AdamicAdar definition so R0
AdamicAdar MRR should land near the previously reported ~0.69.

Regimes come from split_manifest.json via lib_eval.load_regime. R2 (degree null) is
skipped automatically until Unit 7 builds train_R2_degree_null_seed42.csv.

Out: data/processed/results/baselines/baselines_<regime>.json  (one file per regime,
written as each regime finishes -> safe to interrupt/resume regime-by-regime).

Usage:
    python scripts/run_baselines.py                          # R0,R1,R3 x seeds 42,1,7
    python scripts/run_baselines.py --regimes R0 --seeds 42  # quick single-regime look
    python scripts/run_baselines.py --regimes R0,R1,R2,R3    # once Unit 7 has built R2
"""
import argparse
import json
import math
import time
from pathlib import Path

import numpy as np

try:
    import lib_eval
except ImportError:  # allow running from repo root
    from scripts import lib_eval

try:
    from config import PROCESSED_DIR
    DEFAULT_OUT = Path(PROCESSED_DIR) / "results" / "baselines"
except Exception:
    DEFAULT_OUT = Path("data/processed/results/baselines")

METHOD_ORDER = ["Random", "CommonNeighbors", "AdamicAdar", "Jaccard",
                "PreferentialAttachment"]


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def build_adjacency(train_edges):
    """Undirected neighbor sets + degree from a (N,3) edge array (relation ignored)."""
    adj = {}
    src = train_edges[:, 0]
    tgt = train_edges[:, 2]
    for a, b in zip(src.tolist(), tgt.tolist()):
        if a == b:
            continue
        s = adj.get(a)
        if s is None:
            adj[a] = s = set()
        s.add(b)
        s = adj.get(b)
        if s is None:
            adj[b] = s = set()
        s.add(a)
    deg = {n: len(nb) for n, nb in adj.items()}
    return adj, deg


def make_scorers(adj, deg, hub_cap, seed):
    """Build the five score_fns. Only Random depends on the seed."""
    log_ = math.log

    def common_neighbors(g, d):
        ng = adj.get(g)
        nd = adj.get(d)
        if not ng or not nd:
            return 0.0
        if len(ng) > len(nd):
            ng, nd = nd, ng
        c = 0
        for w in ng:
            if w in nd and deg[w] <= hub_cap:
                c += 1
        return float(c)

    def adamic_adar(g, d):
        ng = adj.get(g)
        nd = adj.get(d)
        if not ng or not nd:
            return 0.0
        if len(ng) > len(nd):
            ng, nd = nd, ng
        s = 0.0
        for w in ng:
            if w in nd:
                dw = deg[w]
                if 1 < dw <= hub_cap:
                    s += 1.0 / log_(dw)
        return s

    def jaccard(g, d):
        ng = adj.get(g)
        nd = adj.get(d)
        if not ng or not nd:
            return 0.0
        inter = len(ng & nd)
        if inter == 0:
            return 0.0
        return inter / (len(ng) + len(nd) - inter)

    def pref_attachment(g, d):
        return float(deg.get(g, 0) * deg.get(d, 0))

    rng = np.random.default_rng(seed + 10_000)

    def random_scorer(g, d):
        return float(rng.random())

    return {
        "Random": random_scorer,
        "CommonNeighbors": common_neighbors,
        "AdamicAdar": adamic_adar,
        "Jaccard": jaccard,
        "PreferentialAttachment": pref_attachment,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regimes", default="R0,R1,R3",
                    help="comma-separated regimes (R2 auto-skips until Unit 7 builds it)")
    ap.add_argument("--seeds", default="42,1,7",
                    help="comma-separated negative-sampling seeds")
    ap.add_argument("--methods", default=",".join(METHOD_ORDER),
                    help="comma-separated subset of methods to run")
    ap.add_argument("--n-neg", type=int, default=50)
    ap.add_argument("--hub-cap", type=int, default=2000,
                    help="skip shared neighbours with degree above this (CN/AA)")
    ap.add_argument("--splits-dir", default=None)
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    args = ap.parse_args()

    regimes = [r.strip() for r in args.regimes.split(",") if r.strip()]
    seeds = [int(s) for s in args.seeds.split(",") if s.strip()]
    want_methods = [m.strip() for m in args.methods.split(",") if m.strip()]
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    for regime in regimes:
        log(f"================= regime {regime} =================")
        try:
            reg = lib_eval.load_regime(regime, args.splits_dir)
        except FileNotFoundError as e:
            log(f"SKIP {regime}: {e}")
            continue
        train, test, hubs = reg
        log(f"building adjacency over {len(train):,} train edges ...")
        t0 = time.time()
        adj, deg = build_adjacency(train)
        log(f"  {len(adj):,} nodes in adjacency ({time.time() - t0:.1f}s)")
        t0 = time.time()
        pools = lib_eval._category_pools(train, hubs)      # once per regime
        known = lib_eval._known_tails_by_head(train, test, test[:, 0])
        log(f"  built category pools + known-tails ({time.time() - t0:.1f}s)")

        results = {}          # method -> metric -> [per-seed values]
        per_edge_rr = {}      # method -> per-edge reciprocal ranks at the FIRST seed
        for seed in seeds:
            scorers = make_scorers(adj, deg, args.hub_cap, seed)
            for name in want_methods:
                fn = scorers[name]
                t0 = time.time()
                ranks, pos, neg = lib_eval.rank_test_edges(
                    fn, train, test, hubs, reg.hub_filter,
                    n_neg=args.n_neg, seed=seed, return_scores=True,
                    pools=pools, known=known)
                rm = lib_eval.ranking_metrics(ranks)
                cm = lib_eval.classification_metrics(pos, neg)
                rec = results.setdefault(name, {k: [] for k in
                                                ("MRR", "Hits@1", "Hits@3", "Hits@10",
                                                 "AUROC", "AUPRC")})
                for k in ("MRR", "Hits@1", "Hits@3", "Hits@10"):
                    rec[k].append(rm[k])
                rec["AUROC"].append(cm["AUROC"])
                rec["AUPRC"].append(cm["AUPRC"])
                if seed == seeds[0]:
                    per_edge_rr[name] = 1.0 / ranks
                log(f"  {regime} seed={seed} {name:23s} "
                    f"MRR={rm['MRR']:.4f} H@10={rm['Hits@10']:.4f} "
                    f"AUROC={cm['AUROC']:.4f} ({time.time() - t0:.1f}s)")

        # -- summarize: mean +/- sd over seeds; within-seed bootstrap CI on RR ---
        summary = {}
        for name, rec in results.items():
            summary[name] = {}
            for metric, vals in rec.items():
                v = np.asarray(vals, dtype=float)
                summary[name][metric] = {
                    "mean": float(np.mean(v)),
                    "sd": float(np.std(v, ddof=1)) if len(v) > 1 else 0.0,
                    "per_seed": [float(x) for x in v],
                }
            summary[name]["MRR_ci_seed%d" % seeds[0]] = lib_eval.bootstrap_ci(
                per_edge_rr[name])

        # -- paired bootstrap: AdamicAdar vs each other method (first seed) ------
        paired = {}
        if "AdamicAdar" in per_edge_rr:
            base = per_edge_rr["AdamicAdar"]
            for name, rr in per_edge_rr.items():
                if name != "AdamicAdar":
                    paired[f"AdamicAdar_vs_{name}"] = lib_eval.paired_bootstrap_pvalue(
                        base, rr)

        payload = {
            "regime": regime,
            "regime_key": reg.name,
            "train_file": reg.train_file,
            "test_file": reg.test_file,
            "n_train_edges": int(len(train)),
            "n_test_edges": int(len(test)),
            "n_neg": args.n_neg,
            "hub_cap": args.hub_cap,
            "seeds": seeds,
            "methods": summary,
            f"paired_vs_AdamicAdar_seed{seeds[0]}": paired,
            "created": time.strftime("%Y-%m-%d %H:%M:%S"),
        }
        p = out / f"baselines_{regime}.json"
        p.write_text(json.dumps(payload, indent=2))
        log(f"wrote {p}")
        log(f"--- {regime} summary (mean MRR / Hits@10 / AUROC over {len(seeds)} seed(s)) ---")
        for name in METHOD_ORDER:
            if name in summary:
                s = summary[name]
                log(f"    {name:23s} MRR={s['MRR']['mean']:.4f}  "
                    f"H@10={s['Hits@10']['mean']:.4f}  AUROC={s['AUROC']['mean']:.4f}")
    log("DONE.")


if __name__ == "__main__":
    main()
