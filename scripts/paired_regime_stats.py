#!/usr/bin/env python3
"""paired_regime_stats.py -- cross-regime paired bootstrap for the KGE arm.

The manuscript reports each regime effect as a PAIRED difference over the 4,228 held-out
edges (for example "TransE -0.189, 95% CI -0.201 to -0.176, p < 0.001") rather than as a
bare percentage. Those numbers were previously produced ad hoc, drifted out of sync with
the frozen result JSONs when the KGE grid was re-run, and cannot be regenerated from the
repository. This script recomputes them all from the stored per-edge reciprocal ranks so
the prose cites reproducible values.

Convention (matches lib_eval and the rest of the benchmark):
  * Per-edge reciprocal ranks are averaged over the training seeds first, giving one
    aligned vector per (model, regime). The ranking negative-sampling seed is fixed at 42
    for every run, so element i is the same held-out edge in every regime and the pairing
    is valid.
  * delta = MRR(regime) - MRR(R0), so a drop is negative.
  * 1,000 bootstrap resamples, seed 42, two-sided p-value (lib_eval.paired_bootstrap_pvalue).

Usage:
    python scripts/paired_regime_stats.py
    python scripts/paired_regime_stats.py --results-dir data/processed/results/kge_fullrank
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib_eval as le  # noqa: E402

try:
    from config import PROCESSED_DIR
    RESULTS = Path(PROCESSED_DIR) / "results"
except Exception:
    RESULTS = Path("data/processed/results")

REGIMES = ["R1", "R2", "R3"]


def load_rr(results_dir):
    """(model, regime) -> {seed: per-edge reciprocal-rank array}."""
    rr = defaultdict(dict)
    for p in sorted(glob.glob(os.path.join(results_dir, "kge_*_seed*.json"))):
        try:
            d = json.loads(Path(p).read_text())
        except Exception:
            continue
        if "model" not in d or "reciprocal_ranks" not in d:
            continue
        rr[(d["model"], d["regime"])][d["seed"]] = np.asarray(d["reciprocal_ranks"], float)
    return rr


def mean_rr(rr, model, regime):
    per = rr.get((model, regime))
    if not per:
        return None
    mats = [per[s] for s in sorted(per)]
    n = {len(m) for m in mats}
    if len(n) != 1:
        raise SystemExit(f"[paired_regime_stats] {model}/{regime}: seeds disagree on edge "
                         f"count {n}; the runs are not aligned and cannot be paired.")
    return np.mean(np.vstack(mats), axis=0), sorted(per)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", default=str(RESULTS / "kge"))
    ap.add_argument("--out", default=None, help="output JSON (default: <results-dir>/paired_regime_stats.json)")
    args = ap.parse_args()
    out_path = Path(args.out) if args.out else Path(args.results_dir) / "paired_regime_stats.json"

    rr = load_rr(args.results_dir)
    models = sorted({m for m, _r in rr})
    if not models:
        raise SystemExit(f"[paired_regime_stats] no run JSONs with reciprocal_ranks in {args.results_dir}")

    out = {}
    lines = ["| model | regime | MRR R0 | MRR regime | delta | 95% CI | p | n edges | seeds |",
             "|---|---|---|---|---|---|---|---|---|"]
    for model in models:
        base = mean_rr(rr, model, "R0")
        if base is None:
            print(f"  skip {model}: no R0 run")
            continue
        a, seeds0 = base
        for regime in REGIMES:
            got = mean_rr(rr, model, regime)
            if got is None:
                continue
            b, seeds = got
            s = le.paired_bootstrap_pvalue(b, a)
            key = f"{model}|R0->{regime}"
            out[key] = dict(model=model, regime=regime, n_edges=s["n"], seeds=seeds,
                            mrr_R0=float(a.mean()), mrr_regime=float(b.mean()),
                            delta=s["delta"], ci_low=s["ci_low"], ci_high=s["ci_high"],
                            p_value=s["p_value"],
                            pct=100.0 * (float(b.mean()) - float(a.mean())) / float(a.mean()))
            pstr = "<0.001" if s["p_value"] < 0.001 else f"{s['p_value']:.3f}"
            lines.append(f"| {model} | R0->{regime} | {a.mean():.4f} | {b.mean():.4f} | "
                         f"{s['delta']:+.4f} | [{s['ci_low']:+.4f}, {s['ci_high']:+.4f}] | "
                         f"{pstr} | {s['n']} | {','.join(map(str, seeds))} |")
            print(f"  {model:9s} R0->{regime}: delta={s['delta']:+.4f} "
                  f"CI [{s['ci_low']:+.4f}, {s['ci_high']:+.4f}] p={pstr} "
                  f"({out[key]['pct']:+.1f}%)")

    payload = {"comparisons": out,
               "method": ("per-edge reciprocal ranks averaged over training seeds, then "
                          "paired bootstrap (1000 resamples, seed 42, two-sided); "
                          "delta = MRR(regime) - MRR(R0)"),
               "results_dir": str(args.results_dir)}
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2))
    out_path.with_suffix(".md").write_text("\n".join(lines) + "\n")
    print(f"wrote {out_path}")
    print(f"wrote {out_path.with_suffix('.md')}")


if __name__ == "__main__":
    main()
