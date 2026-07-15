"""Build the manuscript figures from the frozen result JSONs (no hand-typed numbers).

Matplotlib only (no seaborn), 300 dpi, deterministic, colorblind-safe. Every
value is read from a result JSON; the script exits non-zero if an input is
missing. One fixed method->color map (Okabe-Ito, a published CVD-safe palette)
is shared across all figures so a method keeps its colour everywhere.

  Fig 1 -- AUROC-vs-MRR dissociation (R0): AUROC compresses the method
           differences that MRR reveals; KGE saturates AUROC near 1.
  Fig 2 -- The leakage audit: MRR per method across R0->R1->R2->R3 with +/-SD
           bands (baseline and KGE panels). R1/R3 flat; the R2 degree-null is the
           visible drop.
  Fig 3 -- Degree-vs-structure decomposition: each method's R0 MRR split into a
           degree component (mean degree-null MRR) and a structure residual
           (R0 - null), with the permutation p-value annotated.

Usage:  python scripts/make_figures.py
"""
from __future__ import annotations

import json
import os
import sys

import matplotlib

matplotlib.use("Agg")  # headless / deterministic
import matplotlib.pyplot as plt  # noqa: E402

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
RESULTS = os.path.join(REPO_ROOT, "data", "processed", "results")
FIG_DIR = os.path.join(REPO_ROOT, "figures")

BASELINES = ["Random", "CommonNeighbors", "AdamicAdar", "Jaccard", "PreferentialAttachment"]
KGE_MODELS = ["TransE", "RotatE"]
REGIMES = ["R0", "R1", "R2", "R3"]

# Okabe-Ito CVD-safe palette; one fixed hue per method, used in every figure.
COLOR = {
    "Random": "#999999",
    "CommonNeighbors": "#56B4E9",
    "AdamicAdar": "#0072B2",
    "Jaccard": "#009E73",
    "PreferentialAttachment": "#E69F00",
    "TransE": "#D55E00",
    "RotatE": "#CC79A7",
}
MARKER = {
    "Random": "X",
    "CommonNeighbors": "o",
    "AdamicAdar": "s",
    "Jaccard": "^",
    "PreferentialAttachment": "D",
    "TransE": "P",
    "RotatE": "v",
}

plt.rcParams.update({
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "font.size": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": "#DDDDDD",
    "grid.linewidth": 0.6,
    "axes.axisbelow": True,
    "legend.frameon": False,
})


# --------------------------------------------------------------------------- IO
def _load(path: str) -> dict:
    if not os.path.exists(path):
        sys.exit(f"[make_figures] required input missing: {os.path.relpath(path, REPO_ROOT)}")
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def _baselines() -> dict[str, dict]:
    return {r: _load(os.path.join(RESULTS, "baselines", f"baselines_{r}.json")) for r in REGIMES}


def _kge() -> dict:
    return _load(os.path.join(RESULTS, "kge", "kge_summary.json"))["aggregate"]


def _degree_null() -> dict:
    return _load(os.path.join(RESULTS, "null", "degree_null.json"))


def _bl(bl_regime: dict, method: str, metric: str) -> tuple[float, float]:
    m = bl_regime["methods"][method][metric]
    return float(m["mean"]), float(m["sd"])


def _kge_val(kge: dict, model: str, regime: str, field: str) -> tuple[float, float]:
    rec = kge[f"{model}|{regime}|d64|e300"]
    return float(rec[f"{field}_mean"]), float(rec[f"{field}_sd"])


def _save(fig, name: str) -> None:
    os.makedirs(FIG_DIR, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG_DIR, f"{name}.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote figures/{name}.png and figures/{name}.pdf")


# ------------------------------------------------------------------- Figure 1
def fig1_dissociation(baselines: dict, kge: dict) -> None:
    fig, ax = plt.subplots(figsize=(5.0, 4.0))
    for m in BASELINES:
        mrr, _ = _bl(baselines["R0"], m, "MRR")
        auroc, _ = _bl(baselines["R0"], m, "AUROC")
        ax.scatter(mrr, auroc, s=70, color=COLOR[m], marker=MARKER[m],
                   edgecolor="white", linewidth=0.7, zorder=3, label=m)
    for m in KGE_MODELS:
        mrr, _ = _kge_val(kge, m, "R0", "MRR")
        auroc, _ = _kge_val(kge, m, "R0", "AUROC_type")
        ax.scatter(mrr, auroc, s=70, color=COLOR[m], marker=MARKER[m],
                   edgecolor="white", linewidth=0.7, zorder=3, label=m)

    ax.axhline(0.5, color="#BBBBBB", lw=1.0, ls=":", zorder=1)
    ax.set_xlabel("MRR (ranking quality)")
    ax.set_ylabel("AUROC (type-matched negatives)")
    ax.set_xlim(0, 1.0)
    ax.set_ylim(0.45, 1.02)
    ax.set_title("AUROC compresses the differences MRR reveals (R0)")
    ax.legend(loc="lower right", fontsize=7, ncol=1)
    _save(fig, "fig1_auroc_mrr_dissociation")


# ------------------------------------------------------------------- Figure 2
def fig2_audit(baselines: dict, kge: dict) -> None:
    x = list(range(len(REGIMES)))
    fig, (axb, axk) = plt.subplots(1, 2, figsize=(9.0, 4.0), sharey=True)

    for m in BASELINES:
        means, sds = [], []
        for r in REGIMES:
            mu, sd = _bl(baselines[r], m, "MRR")
            means.append(mu); sds.append(sd)
        lo = [mu - sd for mu, sd in zip(means, sds)]
        hi = [mu + sd for mu, sd in zip(means, sds)]
        axb.plot(x, means, color=COLOR[m], marker=MARKER[m], ms=6, lw=2, label=m, zorder=3)
        axb.fill_between(x, lo, hi, color=COLOR[m], alpha=0.15, zorder=2)

    for m in KGE_MODELS:
        means, sds = [], []
        for r in REGIMES:
            mu, sd = _kge_val(kge, m, r, "MRR")
            means.append(mu); sds.append(sd)
        lo = [mu - sd for mu, sd in zip(means, sds)]
        hi = [mu + sd for mu, sd in zip(means, sds)]
        axk.plot(x, means, color=COLOR[m], marker=MARKER[m], ms=6, lw=2, label=m, zorder=3)
        axk.fill_between(x, lo, hi, color=COLOR[m], alpha=0.15, zorder=2)

    for ax, title in ((axb, "Topological baselines"), (axk, "KGE models")):
        ax.set_xticks(x)
        ax.set_xticklabels(REGIMES)
        ax.set_xlabel("Leakage regime")
        ax.set_title(title)
        # Highlight the degree-null column (R2 = index 2).
        ax.axvspan(1.6, 2.4, color="#F0E442", alpha=0.18, zorder=0)
        ax.legend(fontsize=7)
    axb.set_ylabel("MRR (mean ± SD, 3 seeds)")
    fig.suptitle("Leakage audit across regimes R0–R3 (degree-null R2 shaded)", y=1.02)
    _save(fig, "fig2_leakage_audit")


# ------------------------------------------------------------------- Figure 3
def fig3_decomposition(dn: dict) -> None:
    methods = dn["methods"]
    order = [m for m in ["CommonNeighbors", "AdamicAdar", "Jaccard",
                         "PreferentialAttachment"] if m in methods]
    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    x = list(range(len(order)))
    degree_vals = [float(methods[m]["null_MRR_mean"]) for m in order]
    struct_vals = [float(methods[m]["real_minus_null_mean"]) for m in order]

    ax.bar(x, degree_vals, width=0.62, color="#BBBBBB", label="degree component (null MRR)",
           zorder=3)
    # 2px surface gap between stacked segments (thin white separator).
    bars = ax.bar(x, struct_vals, width=0.62, bottom=degree_vals,
                  color=[COLOR[m] for m in order],
                  label="structure residual (R0 − null)", zorder=3,
                  edgecolor="white", linewidth=1.2)

    for xi, m in zip(x, order):
        real = float(methods[m]["real_R0_MRR"])
        p = float(methods[m]["p_value"])
        p_str = "p<0.001" if p == 0.0 else f"p={p:.3f}"
        ax.text(xi, real + 0.012, p_str, ha="center", va="bottom", fontsize=7,
                color="#333333")

    ax.set_xticks(x)
    ax.set_xticklabels(order, rotation=20, ha="right")
    ax.set_ylabel("R0 MRR")
    ax.set_ylim(0, max(float(methods[m]["real_R0_MRR"]) for m in order) * 1.18)
    ax.set_title("How much of R0 MRR is pure degree vs genuine structure")
    ax.legend(fontsize=7, loc="upper left")
    _save(fig, "fig3_degree_vs_structure")


def main() -> None:
    print("Building manuscript figures from result JSONs...")
    baselines = _baselines()
    kge = _kge()
    dn = _degree_null()
    fig1_dissociation(baselines, kge)
    fig2_audit(baselines, kge)
    fig3_decomposition(dn)
    print(f"\nAll figures written to {os.path.relpath(FIG_DIR, REPO_ROOT)}/")


if __name__ == "__main__":
    main()
