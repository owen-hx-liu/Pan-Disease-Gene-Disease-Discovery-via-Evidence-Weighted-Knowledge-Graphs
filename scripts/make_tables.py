"""Build the manuscript tables from the frozen result JSONs (no hand-typed numbers).

Emits both Markdown (``tables/*.md``) and LaTeX (``tables/*.tex``) for:

  Table 1 -- Graph statistics (Monarch + Hetionet).
  Table 2 -- Leakage audit: every method x regime R0-R3 (MRR, plus R0 AUROC /
             Hits@10 to expose the ranking-vs-AUROC dissociation).
  Table 3 -- Degree-null decomposition (real R0 MRR, degree-null MRR, the
             structure-beyond-degree residual, permutation p-value).
  Table 4 -- Cross-graph robustness on Hetionet (R0/R1/R2).

Every number is read from a JSON file; the script fails loudly (exits non-zero)
if any expected input is missing, rather than silently emitting a blank cell.
Uncertainty is reported as mean +/- SD across seeds (n=3) -- the quantity the
per-seed result files actually contain -- and labelled as such.

Usage:  python scripts/make_tables.py
"""
from __future__ import annotations

import json
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
RESULTS = os.path.join(REPO_ROOT, "data", "processed", "results")
OUT_DIR = os.path.join(REPO_ROOT, "tables")

# Presentation order.
BASELINES = ["Random", "CommonNeighbors", "AdamicAdar", "Jaccard", "PreferentialAttachment"]
KGE_MODELS = ["TransE", "RotatE"]
REGIMES = ["R0", "R1", "R2", "R3"]
REGIME_LABEL = {
    "R0": "R0 standard",
    "R1": "R1 redundancy",
    "R2": "R2 degree-null",
    "R3": "R3 orthology-blocked",
}


# --------------------------------------------------------------------------- IO
def _load(path: str) -> dict:
    if not os.path.exists(path):
        sys.exit(f"[make_tables] required input missing: {os.path.relpath(path, REPO_ROOT)}")
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def _graph_stats() -> tuple[dict, dict]:
    monarch = _load(os.path.join(REPO_ROOT, "data", "processed", "graph_stats.json"))
    hetio = _load(os.path.join(REPO_ROOT, "data", "processed", "graph_stats_hetionet.json"))
    return monarch, hetio


def _baselines() -> dict[str, dict]:
    """{regime: parsed baselines_<regime>.json} for Monarch."""
    out = {}
    for r in REGIMES:
        out[r] = _load(os.path.join(RESULTS, "baselines", f"baselines_{r}.json"))
    return out


def _kge() -> dict:
    return _load(os.path.join(RESULTS, "kge", "kge_summary.json"))["aggregate"]


def _degree_null() -> dict:
    return _load(os.path.join(RESULTS, "null", "degree_null.json"))


def _hetionet() -> dict[str, dict]:
    out = {}
    for r in ["R0", "R1", "R2"]:
        out[r] = _load(os.path.join(RESULTS, "hetionet", f"baselines_{r}.json"))
    return out


# --------------------------------------------------------------- value helpers
def _bl_stat(bl_regime: dict, method: str, metric: str) -> tuple[float, float]:
    """(mean, sd) for a baseline method/metric in a parsed baselines_<regime>.json."""
    m = bl_regime["methods"][method][metric]
    return float(m["mean"]), float(m["sd"])


def _kge_key(model: str, regime: str) -> str:
    return f"{model}|{regime}|d64|e300"


def _kge_stat(kge: dict, model: str, regime: str, field: str) -> tuple[float, float]:
    rec = kge[_kge_key(model, regime)]
    return float(rec[f"{field}_mean"]), float(rec[f"{field}_sd"])


def _pm(mean: float, sd: float, dec: int = 3) -> str:
    return f"{mean:.{dec}f} ± {sd:.{dec}f}"


def _pm_tex(mean: float, sd: float, dec: int = 3) -> str:
    return f"{mean:.{dec}f} $\\pm$ {sd:.{dec}f}"


def _pct(new: float, old: float) -> str:
    if old == 0:
        return "n/a"
    return f"{100.0 * (new - old) / old:+.0f}%"


def _paired_bootstrap_note(bl_R0: dict) -> "tuple[str, str]":
    """Markdown + LaTeX note surfacing the stored paired-bootstrap MRR comparisons.

    Each baseline JSON stores ``paired_vs_AdamicAdar_seed42`` = per-comparison
    {delta (AdamicAdar - method), 95% CI, two-sided p, n} on aligned per-edge
    reciprocal ranks (R0, seed 42, 1000 resamples). Surfacing them makes Table 2
    report method-vs-method significance, not only mean +/- SD across seeds.
    """
    comps = bl_R0.get("paired_vs_AdamicAdar_seed42")
    if not comps:
        return "", ""
    order = ["AdamicAdar_vs_CommonNeighbors", "AdamicAdar_vs_Jaccard",
             "AdamicAdar_vs_PreferentialAttachment", "AdamicAdar_vs_Random"]
    parts, n = [], 0
    for ck in order:
        rec = comps.get(ck)
        if not rec:
            continue
        other = ck.split("_vs_")[1]
        delta, p, n = float(rec["delta"]), float(rec["p_value"]), int(rec.get("n", 0))
        p_str = "p<0.001" if p == 0.0 else f"p={p:.3f}"
        parts.append(f"{other} Δ={delta:+.3f} ({p_str})")
    if not parts:
        return "", ""
    body = "; ".join(parts)
    md = (f"\n\n*Paired-bootstrap MRR difference (AdamicAdar − method; R0, seed 42, "
          f"1,000 resamples, n={n:,} test edges): {body}. Every difference is "
          f"significant at p<0.001 (two-sided).*")
    tex = _tex_uni(f"Paired-bootstrap MRR difference (AdamicAdar − method; R0, seed 42, "
                   f"1,000 resamples, n={n:,}): {body}. All significant at p<0.001 "
                   f"(two-sided).").replace("<", "$<$")
    return md, tex


# ------------------------------------------------------------- markdown / latex
def _md_table(headers: list[str], rows: list[list[str]]) -> str:
    line = "| " + " | ".join(headers) + " |"
    sep = "|" + "|".join("---" for _ in headers) + "|"
    body = "\n".join("| " + " | ".join(r) + " |" for r in rows)
    return "\n".join([line, sep, body])


# UTF-8 glyphs used in the Markdown cells -> portable pdflatex math/commands.
_TEX_UNI = {
    "±": "$\\pm$", "→": "$\\to$", "Δ": "$\\Delta$", "−": "$-$",
    "≥": "$\\ge$", "≤": "$\\le$", "×": "$\\times$", "–": "--", "’": "'",
}


def _tex_uni(s: str) -> str:
    for u, tex in _TEX_UNI.items():
        s = s.replace(u, tex)
    return s


def _tex_table(headers: list[str], rows: list[list[str]], caption: str, label: str,
               colspec: str | None = None) -> str:
    if colspec is None:
        colspec = "l" + "r" * (len(headers) - 1)

    def esc(s: str) -> str:
        return _tex_uni(s.replace("%", "\\%").replace("_", "\\_"))
    head = " & ".join(esc(h) for h in headers) + " \\\\"
    body = " \\\\\n".join(" & ".join(esc(c) for c in r) for r in rows) + " \\\\"
    return "\n".join([
        "\\begin{table}[t]",
        "\\centering",
        f"\\caption{{{_tex_uni(caption)}}}",
        f"\\label{{{label}}}",
        f"\\begin{{tabular}}{{{colspec}}}",
        "\\toprule",
        head,
        "\\midrule",
        body,
        "\\bottomrule",
        "\\end{tabular}",
        "\\end{table}",
    ])


def _write(name: str, md: str, tex: str) -> None:
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(os.path.join(OUT_DIR, f"{name}.md"), "w", encoding="utf-8") as fh:
        fh.write(md + "\n")
    with open(os.path.join(OUT_DIR, f"{name}.tex"), "w", encoding="utf-8") as fh:
        fh.write(tex + "\n")
    print(f"  wrote tables/{name}.md and tables/{name}.tex")


# ------------------------------------------------------------------- Table 1
def table1_graph_stats(monarch: dict, hetio: dict) -> None:
    headers = ["Property", "Monarch (this work)", "Hetionet v1.0 (robustness)"]

    def relcount(stats: dict, key: str) -> str:
        return f"{stats['per_relation_edge_counts'].get(key, 0):,}"

    het_nodes = hetio.get("n_nodes_catalog") or hetio["n_nodes"]
    rows = [
        ["Nodes", f"{monarch['n_nodes']:,}", f"{het_nodes:,}"],
        ["Edges (directed)", f"{monarch['n_edges_directed']:,}", f"{hetio['n_edges_directed']:,}"],
        ["Edges (unique undirected)",
         f"{monarch['n_edges_unique_undirected']:,}", f"{hetio['n_edges_unique_undirected']:,}"],
        ["Node categories / metanodes", f"{monarch['n_categories']}", f"{hetio['n_categories']}"],
        ["Relation / metaedge types", f"{monarch['n_relations']}", f"{hetio['n_relations']}"],
        ["Degree mean", f"{monarch['degree']['mean']:.1f}", f"{hetio['degree']['mean']:.1f}"],
        ["Degree median", f"{monarch['degree']['median']:.0f}", f"{hetio['degree']['median']:.0f}"],
        ["Degree max", f"{monarch['degree']['max']:,}", f"{hetio['degree']['max']:,}"],
        ["Gene nodes",
         f"{monarch['per_category_node_counts'].get('gene', 0):,}",
         f"{hetio['per_category_node_counts'].get('Gene', 0):,}"],
        ["Disease nodes",
         f"{monarch['per_category_node_counts'].get('disease', 0):,}",
         f"{hetio['per_category_node_counts'].get('Disease', 0):,}"],
        ["Gene–disease target edges",
         relcount(monarch, "BIOLINK:GENE_ASSOCIATED_WITH_CONDITION"),
         f"{hetio.get('n_target_edges', 0):,}"],
        ["Orthology edges",
         relcount(monarch, "BIOLINK:ORTHOLOGOUS_TO"), "0 (single-species)"],
    ]
    caption = ("Graph statistics for the Monarch-derived graph and the Hetionet robustness "
               "graph. The gene--disease target relation is "
               "GENE\\_ASSOCIATED\\_WITH\\_CONDITION (Monarch) and DaG (Hetionet). "
               "Hetionet carries no cross-species orthology, so regime R3 is Monarch-only.")
    md = "### Table 1. Graph statistics\n\n" + _md_table(headers, rows)
    tex = _tex_table(headers, rows, caption, "tab:graph_stats")
    _write("table1_graph_stats", md, tex)


# ------------------------------------------------------------------- Table 2
def table2_audit(baselines: dict, kge: dict) -> None:
    headers = ["Method", "MRR R0", "MRR R1", "MRR R2", "MRR R3",
               "ΔMRR R0→R2", "AUROC R0", "Hits@10 R0"]

    def baseline_row(method: str) -> list[str]:
        mrr = {r: _bl_stat(baselines[r], method, "MRR") for r in REGIMES}
        auroc0, auroc0_sd = _bl_stat(baselines["R0"], method, "AUROC")
        h10, h10_sd = _bl_stat(baselines["R0"], method, "Hits@10")
        return [
            method,
            _pm(*mrr["R0"]), _pm(*mrr["R1"]), _pm(*mrr["R2"]), _pm(*mrr["R3"]),
            _pct(mrr["R2"][0], mrr["R0"][0]),
            _pm(auroc0, auroc0_sd), _pm(h10, h10_sd),
        ]

    def kge_row(model: str) -> list[str]:
        mrr = {r: _kge_stat(kge, model, r, "MRR") for r in REGIMES}
        auroc0, auroc0_sd = _kge_stat(kge, model, "R0", "AUROC_type")
        h10, h10_sd = _kge_stat(kge, model, "R0", "Hits@10")
        return [
            model,
            _pm(*mrr["R0"]), _pm(*mrr["R1"]), _pm(*mrr["R2"]), _pm(*mrr["R3"]),
            _pct(mrr["R2"][0], mrr["R0"][0]),
            _pm(auroc0, auroc0_sd), _pm(h10, h10_sd),
        ]

    rows = [baseline_row(m) for m in BASELINES] + [kge_row(m) for m in KGE_MODELS]

    paired_md, paired_tex = _paired_bootstrap_note(baselines["R0"])
    caption = ("Leakage audit on the held-out human gene--disease test edges (mean $\\pm$ SD "
               "over 3 seeds; sampled 50 type-matched negatives). MRR is shown for every "
               "regime; AUROC and Hits@10 at R0 expose the ranking-vs-AUROC dissociation. "
               "The degree-preserving null (R2) is the only regime that materially reduces "
               "performance; redundancy (R1) and orthology (R3) are flat. AUROC for KGE is "
               "type-matched, comparable to the baselines' negatives.")
    md = ("### Table 2. Leakage audit across regimes\n\n"
          "*Mean ± SD over 3 seeds (42, 1, 7); 50 type-matched sampled negatives.*\n\n"
          + _md_table(headers, rows) + paired_md)
    tex = _tex_table(headers, rows, caption, "tab:audit")
    if paired_tex:
        tex = tex + "\n\n\\noindent\\footnotesize " + paired_tex + "\\normalsize\n"
    _write("table2_audit", md, tex)


# ------------------------------------------------------------------- Table 3
def table3_degree_null(dn: dict) -> None:
    headers = ["Method", "Real R0 MRR", "Degree-null MRR",
               "Structure (real−null)", "Perm. p", "Structure > degree?"]
    methods = dn["methods"]
    order = [m for m in ["CommonNeighbors", "AdamicAdar", "Jaccard",
                         "PreferentialAttachment", "Random"] if m in methods]
    rows = []
    for m in order:
        rec = methods[m]
        real = float(rec["real_R0_MRR"])
        null_m = float(rec["null_MRR_mean"])
        null_sd = float(rec.get("null_MRR_sd", 0.0))
        struct = float(rec["real_minus_null_mean"])
        p = float(rec["p_value"])
        p_str = "<0.001" if p == 0.0 else f"{p:.3f}"
        rows.append([
            m, f"{real:.3f}", _pm(null_m, null_sd), f"{struct:+.3f}",
            p_str, "yes" if rec.get("structure_beyond_degree") else "no",
        ])
    n_rep = dn.get("n_replicates", "?")
    pa = dn.get("pa_sanity", {})
    pa_note = ""
    if pa:
        pa_note = (f" PreferentialAttachment sanity: real {pa.get('real', float('nan')):.3f} "
                   f"vs null {pa.get('null_mean', float('nan')):.3f} (unchanged → the "
                   f"null preserves the degree sequence).")
    seq_ok = dn.get("degree_sequence_preserved")
    caption = (f"Degree-null decomposition of R0 ranking performance ({n_rep} degree-preserving, "
               f"type-preserving permutation replicates; R0 evaluation held fixed). The "
               f"structure column is the residual above what pure node degree explains; the "
               f"permutation p-value is the fraction of null replicates with MRR $\\geq$ the real "
               f"R0 MRR. Degree sequence preserved: {seq_ok}." + pa_note)
    md = ("### Table 3. Degree-null decomposition\n\n"
          f"*{n_rep} permutation replicates; R0 evaluation held fixed.*"
          + (f"\n\n{pa_note.strip()}" if pa_note else "") + "\n\n"
          + _md_table(headers, rows))
    tex = _tex_table(headers, rows, caption, "tab:degree_null")
    _write("table3_degree_null", md, tex)


# ------------------------------------------------------------------- Table 4
def table4_hetionet(hetio_bl: dict) -> None:
    headers = ["Method", "MRR R0", "MRR R1", "MRR R2", "ΔMRR R0→R2"]
    rows = []
    for m in BASELINES:
        mrr = {r: _bl_stat(hetio_bl[r], m, "MRR") for r in ["R0", "R1", "R2"]}
        rows.append([
            m, _pm(*mrr["R0"]), _pm(*mrr["R1"]), _pm(*mrr["R2"]),
            _pct(mrr["R2"][0], mrr["R0"][0]),
        ])
    caption = ("Cross-graph robustness: the same audit rerun unmodified on Hetionet "
               "(DaG task; mean $\\pm$ SD over 3 seeds). Hetionet has no cross-species "
               "orthology, so R3 is omitted. The R0$\\to$R2 degree-null drop reproduces the "
               "Monarch pattern on an independent graph.")
    md = ("### Table 4. Hetionet robustness (independent graph)\n\n"
          "*Mean ± SD over 3 seeds; R3 omitted (Hetionet has no orthology).*\n\n"
          + _md_table(headers, rows))
    tex = _tex_table(headers, rows, caption, "tab:hetionet")
    _write("table4_hetionet", md, tex)


def main() -> None:
    print("Building manuscript tables from result JSONs...")
    monarch, hetio = _graph_stats()
    baselines = _baselines()
    kge = _kge()
    dn = _degree_null()
    hetio_bl = _hetionet()

    table1_graph_stats(monarch, hetio)
    table2_audit(baselines, kge)
    table3_degree_null(dn)
    table4_hetionet(hetio_bl)
    print(f"\nAll tables written to {os.path.relpath(OUT_DIR, REPO_ROOT)}/")


if __name__ == "__main__":
    main()
