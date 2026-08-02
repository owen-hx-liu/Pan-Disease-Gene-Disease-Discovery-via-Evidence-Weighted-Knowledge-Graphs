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
KGE_MODELS = ["TransE", "RotatE", "DistMult", "ComplEx"]
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
    # The seven Okabe-Ito hues are taken above; black is the remaining CVD-safe
    # option and stays legible in the greyscale print edition. ComplEx takes the
    # IBM colorblind-safe purple, which is distinct from both the pink (#CC79A7)
    # and the blue (#0072B2) already in use.
    "DistMult": "#000000",
    "ComplEx": "#785EF0",
}
MARKER = {
    "Random": "X",
    "CommonNeighbors": "o",
    "AdamicAdar": "s",
    "Jaccard": "^",
    "PreferentialAttachment": "D",
    "TransE": "P",
    "RotatE": "v",
    "DistMult": "*",
    "ComplEx": "h",
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


def _degree_stratified() -> dict:
    return _load(os.path.join(RESULTS, "degree_stratified", "degree_stratified.json"))


def _gnn() -> tuple[dict, dict]:
    """(degree stratification, regime delta) for the message-passing arm."""
    strat = _load(os.path.join(RESULTS, "gnn", "degree_stratified.json"))
    delta = _load(os.path.join(RESULTS, "gnn", "gnn_regime_delta.json"))
    return strat, delta


def _gnn_refcal() -> list[dict]:
    """Reference-decoder calibration points, one per subgraph fraction, frac ascending."""
    import glob
    out = []
    for path in sorted(glob.glob(os.path.join(RESULTS, "gnn_refcal", "*", "kge_summary.json"))):
        payload = _load(path)
        runs = payload.get("runs", [])
        if not runs:
            continue
        out.append({"frac": float(runs[0]["subgraph_frac"]),
                    "R0": float(payload["aggregate"]["DistMult|R0|d64|e300"]["MRR_mean"]),
                    "R2": float(payload["aggregate"]["DistMult|R2|d64|e300"]["MRR_mean"])})
    if not out:
        sys.exit("[make_figures] no gnn_refcal/*/kge_summary.json found")
    return sorted(out, key=lambda r: r["frac"])


def _schema() -> dict:
    # Supplementary S1/S2 read a single precomputed schema JSON
    # (scripts/graph_schema_stats.py), so this script never touches the raw edges.
    return _load(os.path.join(REPO_ROOT, "data", "processed", "graph_schema.json"))


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
    fig, ax = plt.subplots(figsize=(4.6, 4.6))
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
    ax.set_xlim(0, 1.02)
    ax.set_ylim(0, 1.02)
    ax.set_aspect("equal")
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
        # Add-one estimator (b+1)/(B+1): with B=200 replicates this floors at
        # ~0.005, so structural methods sit at the floor rather than "p<0.001".
        # Matches the Perm. p column of Table 3.
        p = float(methods[m]["p_value_plus_one"])
        p_str = f"p={p:.3f}"
        ax.text(xi, real + 0.012, p_str, ha="center", va="bottom", fontsize=7,
                color="#333333")

    ax.set_xticks(x)
    ax.set_xticklabels(order, rotation=20, ha="right")
    ax.set_ylabel("R0 MRR")
    ax.set_ylim(0, max(float(methods[m]["real_R0_MRR"]) for m in order) * 1.18)
    ax.set_title("How much of R0 MRR is pure degree vs genuine structure")
    ax.legend(fontsize=7, loc="upper left")
    _save(fig, "fig3_degree_vs_structure")


# ------------------------------------------------------------------- Figure 4
def fig4_degree_stratified(ds: dict) -> None:
    """Held-out R0 MRR by degree quartile along two axes -- a double dissociation. Left:
    query-gene degree; the overlap heuristics climb steeply (they need the gene to have
    neighbours) while KGE is flat. Right: ranked-disease degree; KGE and the pure-degree
    Preferential-Attachment anchor climb steeply while the overlap heuristics are flat.
    Different method classes lean on different entities' degree, but every method leans on
    degree; Random is flat on both."""
    methods = ds["meta"].get("methods_covered") or ds["meta"]["models_covered"]
    gene, dis = ds["by_gene_degree"], ds["by_disease_degree"]
    qnames = list(dis["strata"]["quartiles"].keys())
    x = list(range(len(qnames)))

    def series(block: dict, method: str):
        recs = block["results"][f"{method}|R0"]["quartiles"]
        mu = [recs[n]["MRR_mean"] for n in qnames]
        sd = [recs[n]["MRR_sd"] for n in qnames]
        return mu, sd

    fig, (axg, axd) = plt.subplots(1, 2, figsize=(9.5, 4.6), sharey=True)
    handles: list = []
    for ax, block, title in ((axg, gene, "Query-gene degree"),
                             (axd, dis, "Ranked-disease degree")):
        q = block["strata"]["quartiles"]
        for m in methods:
            mu, sd = series(block, m)
            lo = [a - s for a, s in zip(mu, sd)]
            hi = [a + s for a, s in zip(mu, sd)]
            (line,) = ax.plot(x, mu, color=COLOR[m], marker=MARKER[m], ms=6, lw=2,
                              label=m, zorder=3)
            ax.fill_between(x, lo, hi, color=COLOR[m], alpha=0.12, zorder=2)
            if ax is axg:
                handles.append(line)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{n}\n{q[n]['degree_min']}–{q[n]['degree_max']}" for n in qnames])
        ax.set_xlabel("Degree quartile (low → high)")
        ax.set_title(title)
    axg.set_ylabel("R0 MRR (mean ± SD, 3 seeds)")
    # Nine methods no longer fit inside an axes without covering the curves, so the legend
    # goes under the panels as one shared strip.
    fig.legend(handles=handles, labels=[h.get_label() for h in handles],
               fontsize=7.5, ncol=5, loc="lower center", bbox_to_anchor=(0.5, -0.10),
               columnspacing=1.4, handletextpad=0.5)
    fig.suptitle("A double dissociation: overlap heuristics ride query-gene degree; "
                 "KGE and Pref.-Attachment ride target-disease degree", y=1.00, fontsize=10)
    _save(fig, "fig4_degree_stratified")


# ------------------------------------------------------------------ Figure S1
# One fixed hue per node category (Okabe-Ito), used only by the schema figure.
CAT_COLOR = {
    "gene": "#0072B2",
    "disease": "#D55E00",
    "phenotype": "#CC79A7",
    "pathway": "#009E73",
    "variant": "#56B4E9",
    "compound": "#E69F00",
    "anatomy": "#F0E442",
    "therapy": "#999999",
}


def figS1_schema(schema: dict) -> None:
    """Schema / metagraph: the 8 node categories on a ring, sized by node count,
    joined by directed category->category edges whose width is the (log) total
    edge count. Each drawn edge (pairs >=10,000 edges) carries a small numbered
    marker at its midpoint; the numbered legend on the right names the dominant
    relation and edge count, so no long relation text clutters the graph. Smaller
    pairs are omitted for legibility (full list in graph_schema.json)."""
    import math

    from matplotlib.gridspec import GridSpec
    from matplotlib.patches import FancyArrowPatch

    cats = list(schema["category_nodes"].keys())          # already count-desc
    counts = schema["category_nodes"]
    n = len(cats)
    # Fixed clockwise ring; gene at top.
    ang = {c: math.pi / 2 - 2 * math.pi * i / n for i, c in enumerate(cats)}
    R = 1.0
    pos = {c: (R * math.cos(ang[c]), R * math.sin(ang[c])) for c in cats}

    max_cnt = max(counts.values())
    def node_r(c: str) -> float:                          # area ~ sqrt(count)
        return 0.055 + 0.135 * math.sqrt(counts[c] / max_cnt)

    LABEL_MIN = 10_000                                    # only draw major pairs
    edges = [e for e in schema["metagraph"] if e["n_edges"] >= LABEL_MIN]
    emax = math.log10(max(e["n_edges"] for e in edges))
    emin = math.log10(min(e["n_edges"] for e in edges))
    def lw(nn: int) -> float:
        return 1.2 + 6.0 * (math.log10(nn) - emin) / (emax - emin)

    fig = plt.figure(figsize=(12.5, 8.5))
    gs = GridSpec(1, 2, width_ratios=[1.0, 0.72], wspace=0.06)
    ax = fig.add_subplot(gs[0]); axl = fig.add_subplot(gs[1])
    # Extra right margin so the right-most node label (pathway) clears the legend.
    ax.set_xlim(-1.5, 1.95); ax.set_ylim(-1.5, 1.5)
    ax.set_aspect("equal"); ax.axis("off")
    axl.axis("off"); axl.set_xlim(0, 1); axl.set_ylim(0, 1)

    placed: list[tuple[float, float]] = []                # marker centres already drawn

    def draw_marker(x: float, y: float, num: int) -> None:
        ax.text(x, y, str(num), fontsize=8, fontweight="bold", color="#222222",
                ha="center", va="center", zorder=6,
                bbox=dict(boxstyle="circle,pad=0.28", fc="white", ec="#888888", lw=1.0))
        placed.append((x, y))

    def bezier(p0, c, p1, t):                             # quadratic Bezier point
        mt = 1.0 - t
        return (mt * mt * p0[0] + 2 * mt * t * c[0] + t * t * p1[0],
                mt * mt * p0[1] + 2 * mt * t * c[1] + t * t * p1[1])

    RAD = 0.18                                            # arc curvature
    arcs = []                                             # (num, P0, ctrl, P1) for pass 2
    for i, e in enumerate(edges, start=1):
        s, t = e["source"], e["target"]
        w = lw(e["n_edges"])
        (x0, y0), (x1, y1) = pos[s], pos[t]
        if s == t:                                        # self-loop, swung off the node label
            rr = node_r(s)
            ldir = ang[s] + math.radians(55)              # into the empty upper-left
            cx, cy = x0 + 1.5 * rr * math.cos(ldir), y0 + 1.5 * rr * math.sin(ldir)
            ax.add_patch(plt.Circle((cx, cy), rr * 0.62, fill=False, lw=w,
                                    color="#9A9A9A", zorder=1))
            draw_marker(cx + rr * 0.62 * math.cos(ldir), cy + rr * 0.62 * math.sin(ldir), i)
            continue
        dx, dy = x1 - x0, y1 - y0
        d = math.hypot(dx, dy) or 1.0
        ux, uy = dx / d, dy / d
        sx, sy = x0 + ux * node_r(s), y0 + uy * node_r(s)   # trim to node borders
        tx, ty = x1 - ux * node_r(t), y1 - uy * node_r(t)
        curve = RAD if cats.index(t) > cats.index(s) else -RAD
        ax.add_patch(FancyArrowPatch(
            (sx, sy), (tx, ty), connectionstyle=f"arc3,rad={curve}",
            arrowstyle="-|>", mutation_scale=13, lw=w, color="#9A9A9A",
            alpha=0.9, zorder=1, shrinkA=0, shrinkB=0))
        # Control point of matplotlib's Arc3(rad): midpoint + rad*(dy, -dx), where
        # (dx, dy) is the chord. Using the exact form keeps our markers on the
        # drawn curve (an earlier sign error put them on the mirror-image side).
        mx, my = (sx + tx) / 2, (sy + ty) / 2
        ctrl = (mx + curve * (ty - sy), my - curve * (tx - sx))
        arcs.append((i, (sx, sy), ctrl, (tx, ty)))

    # Pass 2: place each number ON its own arc, sliding along the curve to the
    # first position that clears every marker already placed (so none overlap and
    # each stays unambiguously on its arrow). Arrowhead/source ends are avoided.
    MIN_DIST = 0.17
    T_CANDS = [0.50, 0.40, 0.60, 0.33, 0.67, 0.27, 0.73]
    for i, p0, ctrl, p1 in arcs:
        best, best_gap = None, -1.0
        for tt in T_CANDS:
            x, y = bezier(p0, ctrl, p1, tt)
            gap = min((math.hypot(x - a, y - b) for a, b in placed), default=1e9)
            if gap >= MIN_DIST:
                best = (x, y)
                break
            if gap > best_gap:
                best, best_gap = (x, y), gap
        draw_marker(best[0], best[1], i)

    for c in cats:
        x, y = pos[c]
        ax.add_patch(plt.Circle((x, y), node_r(c), color=CAT_COLOR[c],
                                ec="white", lw=1.5, zorder=3))
        # Label just outside the node, radially outward, so text never overflows.
        lx = x + (node_r(c) + 0.11) * math.cos(ang[c])
        ly = y + (node_r(c) + 0.11) * math.sin(ang[c])
        ha = "center" if abs(math.cos(ang[c])) < 0.35 else ("left" if math.cos(ang[c]) > 0 else "right")
        va = "center" if abs(math.sin(ang[c])) < 0.35 else ("bottom" if math.sin(ang[c]) > 0 else "top")
        ax.text(lx, ly, f"{c}\n{counts[c]:,}", ha=ha, va=va, fontsize=10,
                fontweight="bold", color=CAT_COLOR[c], zorder=4)

    # ---- Numbered legend (right panel): number . src -> tgt : relation (count).
    # Rows are top-aligned right under the header (no gap), fixed row height.
    axl.text(0.0, 0.985, "Relations (edge width ∝ log edge count)", fontsize=9.5,
             fontweight="bold", va="top")
    top, row_h = 0.905, 0.086
    for i, e in enumerate(edges, start=1):
        yy = top - (i - 1) * row_h
        axl.text(0.0, yy, str(i), fontsize=8.5, fontweight="bold", color="#222222",
                 va="center", ha="center",
                 bbox=dict(boxstyle="circle,pad=0.28", fc="white", ec="#888888", lw=1.0))
        rel = e["dominant_relation"].replace("_", " ").title()
        axl.text(0.06, yy + 0.020,
                 f"{e['source']} → {e['target']}", fontsize=9, fontweight="bold",
                 va="center", color="#222222")
        axl.text(0.06, yy - 0.022,
                 f"{rel} · {e['n_edges']:,} edges", fontsize=8, va="center",
                 color="#666666")

    fig.suptitle("Graph schema (metagraph): node categories and their relations — "
                 f"{schema['n_nodes']:,} nodes, {schema['n_edges_directed']:,} directed edges",
                 fontsize=11, y=0.97)
    _save(fig, "figS1_schema")


# ------------------------------------------------------------------ Figure S2
def figS2_degree_distribution(schema: dict) -> None:
    """Node-degree distribution on log-log axes -- the heavy tail that motivates
    the degree-null. Points are the empirical count of nodes at each integer
    degree (number of nodes on the x-axis, node degree on the y-axis);
    median/mean/max degree are marked."""
    dc = schema["degree_counts"]
    degs = sorted(int(k) for k in dc)
    cnts = [dc[str(d)] for d in degs]

    total = sum(cnts)
    mean = sum(d * c for d, c in zip(degs, cnts)) / total
    # Median over the node population.
    half, run, median = total / 2, 0, degs[-1]
    for d, c in zip(degs, cnts):
        run += c
        if run >= half:
            median = d
            break
    dmax = degs[-1]

    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    ax.scatter(cnts, degs, s=14, color="#0072B2", edgecolor="none",
               alpha=0.7, zorder=3)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Number of nodes")
    ax.set_ylabel("Node degree (in + out)")

    xmax = max(cnts)
    for val, label, color in ((median, f"median {median}", "#009E73"),
                              (mean, f"mean {mean:.1f}", "#E69F00"),
                              (dmax, f"max {dmax:,}", "#D55E00")):
        ax.axhline(val, color=color, ls="--", lw=1.2, zorder=2)
        ax.text(xmax * 1.25, val, label, va="bottom", ha="right",
                fontsize=7, color=color)

    ax.set_xlim(right=xmax * 3)
    ax.set_title(f"Degree distribution is heavy-tailed "
                 f"({schema['n_nodes']:,} nodes)")
    _save(fig, "figS2_degree_distribution")


# ------------------------------------------------------------------ Figure S3
# The message-passing arm carries only two curves, and neither method is in the
# shared COLOR map (the R-GCN is not one of the audited methods and its DistMult
# is a *different*, subgraph-trained model from the one in Figures 1-4). Giving
# them their own two hues keeps a reader from reading the DistMult curve here as
# the DistMult of the main figures.
GNN_COLOR = {"RGCN": "#0072B2", "DistMult": "#D55E00"}
GNN_LABEL = {"RGCN": "R-GCN (message-passing encoder)",
             "DistMult": "DistMult (matched decoder control)"}


def figS3_gnn(strat: dict, delta: dict, refcal: list, kge: dict) -> None:
    """The message-passing probe: what survives, and why the rest does not.

    Left -- disease-degree stratification, the one measurement in this arm that does
    not route through R2 and is therefore untouched by the failed subgraph-validity
    gate. Same decoder, same graph, same negatives, same seeds; the encoder is the
    only difference, and the message-passing curve is the flatter of the two.

    Right -- the calibration sweep, run with the reference decoder alone, showing why
    the R0->R2 delta on this subgraph cannot be read as degree leakage: the leak
    returns monotonically with subgraph size and only clears the validity threshold
    at a fraction the encoder cannot be trained on in budget.
    """
    dis = strat["by_disease_degree"]
    q = dis["strata"]["quartiles"]
    qnames = list(q.keys())
    x = list(range(len(qnames)))

    full_r0 = float(kge["DistMult|R0|d64|e300"]["MRR_mean"])
    full_r2 = float(kge["DistMult|R2|d64|e300"]["MRR_mean"])
    full_drop = 100.0 * (full_r2 - full_r0) / full_r0
    min_ret = float(delta["subgraph_validity"]["min_retained_fraction"])

    fig, (axs, axc) = plt.subplots(1, 2, figsize=(9.5, 4.2))

    for model in ["RGCN", "DistMult"]:
        rec = dis["results"][f"{model}|R0"]["quartiles"]
        mu = [float(rec[n]["MRR_mean"]) for n in qnames]
        sd = [float(rec[n]["MRR_sd"]) for n in qnames]
        hi_lo = (float(dis["results"][f"{model}|R0"]["halves"]["high"]["MRR_mean"])
                 / float(dis["results"][f"{model}|R0"]["halves"]["low"]["MRR_mean"]))
        axs.plot(x, mu, color=GNN_COLOR[model], marker="o" if model == "RGCN" else "s",
                 ms=6, lw=2, zorder=3, label=f"{GNN_LABEL[model]} — high/low {hi_lo:.2f}×")
        axs.fill_between(x, [a - s for a, s in zip(mu, sd)], [a + s for a, s in zip(mu, sd)],
                         color=GNN_COLOR[model], alpha=0.14, zorder=2)
    axs.set_xticks(x)
    axs.set_xticklabels([f"{n}\n{q[n]['degree_min']}–{q[n]['degree_max']}" for n in qnames])
    axs.set_xlabel("Ranked-disease degree quartile (full-graph R0 degree)")
    axs.set_ylabel("R0 MRR (mean ± SD, 3 seeds)")
    axs.set_ylim(0, 1.0)
    axs.set_title("Message passing rides target popularity less\nthan its own decoder")
    axs.legend(fontsize=7, loc="upper left")

    # ---- right: retained share of the full-graph collapse vs subgraph fraction
    sub = float(delta["subgraph_validity"]["reference_models"]["DistMult"]["subgraph_drop"]) * 100
    pts = ([{"frac": 0.05, "drop": sub}]
           + [{"frac": r["frac"], "drop": 100.0 * (r["R2"] - r["R0"]) / r["R0"]} for r in refcal]
           + [{"frac": 1.0, "drop": full_drop}])
    pts = sorted(pts, key=lambda r: r["frac"])
    fx = [p["frac"] for p in pts]
    fy = [p["drop"] / full_drop for p in pts]
    ok = [v >= min_ret for v in fy]

    axc.plot(fx, fy, color="#666666", lw=1.6, zorder=2)
    axc.scatter([f for f, g in zip(fx, ok) if not g], [v for v, g in zip(fy, ok) if not g],
                s=70, color="#D55E00", marker="X", zorder=4, label="validity gate FAIL")
    axc.scatter([f for f, g in zip(fx, ok) if g], [v for v, g in zip(fy, ok) if g],
                s=70, color="#009E73", marker="o", zorder=4, label="validity gate PASS")
    axc.axhline(min_ret, color="#009E73", ls="--", lw=1.2, zorder=1)
    axc.text(0.052, min_ret + 0.02, f"validity threshold {min_ret:.0%}", fontsize=7,
             color="#009E73")
    # The production arm had to run at frac 0.05; everything trainable in budget
    # sits left of this line and everything valid sits right of it.
    axc.axvspan(0.11, 0.25, color="#999999", alpha=0.12, zorder=0)
    axc.text(0.165, 0.90, "trainable\n|  valid", fontsize=7, ha="center", color="#555555")
    axc.set_xscale("log")
    axc.set_xticks(fx)
    axc.set_xticklabels([f"{f:g}" for f in fx])
    axc.set_xlabel("Training-subgraph fraction (log scale)")
    axc.set_ylabel(f"Share of the full-graph R0→R2 collapse retained\n"
                   f"(1.0 = the full-graph {full_drop:.1f}%)")
    axc.set_ylim(0, 1.08)
    axc.set_title("The leak dissolves with the subsample —\nand returns only out of budget")
    axc.legend(fontsize=7, loc="upper left")

    fig.suptitle("The message-passing probe: the surviving measurement (left) and "
                 "the failed validity gate (right)", y=1.04, fontsize=10)
    _save(fig, "figS3_gnn_probe")


def main() -> None:
    print("Building manuscript figures from result JSONs...")
    baselines = _baselines()
    kge = _kge()
    dn = _degree_null()
    ds = _degree_stratified()
    fig1_dissociation(baselines, kge)
    fig2_audit(baselines, kge)
    fig3_decomposition(dn)
    fig4_degree_stratified(ds)
    schema = _schema()
    figS1_schema(schema)
    figS2_degree_distribution(schema)
    gnn_strat, gnn_delta = _gnn()
    figS3_gnn(gnn_strat, gnn_delta, _gnn_refcal(), kge)
    print(f"\nAll figures written to {os.path.relpath(FIG_DIR, REPO_ROOT)}/")


if __name__ == "__main__":
    main()
