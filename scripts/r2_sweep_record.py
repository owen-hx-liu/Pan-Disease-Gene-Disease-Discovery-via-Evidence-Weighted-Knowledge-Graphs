#!/usr/bin/env python3
"""r2_sweep_record.py -- consolidate the R2 optimization-vs-generalization evidence into
one auditable document (FINDINGS_r2_sweep.md + r2_sweep_findings.json + loss_curves.csv).

Why this exists
---------------
The paper reports that KGE performance collapses on R2, the degree-preserving rewired
null, and reads that as "the topological signal was degree". A reviewer can object that
the comparison is confounded: hyperparameters were tuned on R0, so the R2 drop might be
the recipe failing to OPTIMIZE on a rewired graph rather than the rewired graph carrying
no learnable structure. scripts/run_r2_sweep.ps1 runs the two experiments that separate
those explanations; this script turns their raw JSON into the table a reviewer can read,
so the verdict rests on files on disk rather than on a chat log. Same role as
scripts/bilinear_sweep_record.py plays for the bilinear-collapse claim.

The two experiments
-------------------
  Phase A, train fit  : ~4000 sampled TRAINING edges scored by the SAME evaluator as the
                        test edges. Train fit upper-bounds what the optimizer achieved.
                        R2 fitting its own training edges while failing on held-out ones
                        is a GENERALIZATION failure, not an optimization failure.
  Phase B, matched grid: TransE swept over lr x negatives on R2 AND on R0. The R0 arm is
                        what makes it an argument -- a sweep on R2 alone has no reference
                        for what the same grid buys on an un-rewired graph.

The verdict is computed from stated, falsifiable thresholds (--min-fit-retention,
--min-separation, --min-sweep-gap) rather than asserted, and the script prints the
objection as UPHELD when the numbers say so.

Usage:
    python scripts/r2_sweep_record.py
    python scripts/r2_sweep_record.py --fit-dir DIR --sweep-dir DIR --out-dir DIR
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

try:
    from config import PROCESSED_DIR
    RESULTS = Path(PROCESSED_DIR) / "results"
except Exception:
    RESULTS = Path("data/processed/results")

R2_DIR = RESULTS / "kge_r2_sweep"
REF_REGIME, ALT_REGIME = "R0", "R2"


# --------------------------------------------------------------------------- #
# Loading
# --------------------------------------------------------------------------- #
def load_runs(d: Path) -> list:
    """Load every per-cell run JSON directly under ``d`` (non-recursive)."""
    out = []
    for p in sorted(Path(d).glob("kge_*_seed*.json")):
        try:
            r = json.loads(p.read_text())
        except Exception as e:  # noqa: BLE001
            print(f"  skip {p.name}: {e}")
            continue
        if "model" not in r:            # defensive: not a per-run record
            continue
        r["_file"] = p.name
        out.append(r)
    return out


def pct(new, ref):
    """Percent change of ``new`` relative to ``ref`` (None-safe)."""
    if new is None or ref in (None, 0):
        return None
    return 100.0 * (new - ref) / ref


def ratio(new, ref):
    if new is None or ref in (None, 0):
        return None
    return new / ref


def fmt(v, nd=4):
    if v is None:
        return "--"
    return f"{v:.{nd}f}"


def fmtp(v, nd=1):
    if v is None:
        return "--"
    return f"{v:+.{nd}f}%"


def md_table(headers, rows):
    out = ["| " + " | ".join(headers) + " |",
           "|" + "|".join("---" for _ in headers) + "|"]
    for r in rows:
        out.append("| " + " | ".join("--" if v is None else str(v) for v in r) + " |")
    return "\n".join(out)


# --------------------------------------------------------------------------- #
# Phase A: train fit vs test fit
# --------------------------------------------------------------------------- #
def analyse_trainfit(runs, min_fit_retention, min_separation):
    """Per model: does R2 fit its OWN training edges as well as R0 fits R0's?

    Returns one record per model with the two retention ratios and a verdict.

      test_retention      = R2 test MRR      / R0 test MRR
      train_fit_retention = R2 train-fit MRR / R0 train-fit MRR

    Optimization failure and generalization failure make OPPOSITE predictions here:
    if the recipe simply cannot fit the rewired graph, train fit falls with test fit and
    the two retentions are similar; if the optimizer succeeded and the rewired graph
    merely has nothing that transfers, train-fit retention stays near 1 while test
    retention drops. The separation between the two ratios is therefore the measurement.
    """
    # Phase A is one cell per (model, regime). If the directory ever holds more than one
    # -- e.g. someone points --fit-dir at a sweep directory -- the second would silently
    # replace the first and the reported numbers would depend on filename order.
    by = {}
    for r in runs:
        seen = by.setdefault(r["model"], {})
        if r["regime"] in seen:
            print(f"  WARNING two runs for {r['model']}/{r['regime']} in the train-fit dir "
                  f"({seen[r['regime']]['_file']}, {r['_file']}); using the latter. "
                  f"--fit-dir should hold ONE cell per model x regime.")
        seen[r["regime"]] = r

    recs = []
    for model in sorted(by):
        ref, alt = by[model].get(REF_REGIME), by[model].get(ALT_REGIME)
        if not (ref and alt):
            recs.append(dict(model=model, complete=False,
                             have=sorted(by[model]), verdict="incomplete"))
            continue
        rec = dict(
            model=model, complete=True,
            r0_test=ref.get("MRR"), r2_test=alt.get("MRR"),
            r0_fit=ref.get("train_fit_MRR"), r2_fit=alt.get("train_fit_MRR"),
            r0_test_h10=ref.get("Hits@10"), r2_test_h10=alt.get("Hits@10"),
            r0_fit_h10=(ref.get("train_fit") or {}).get("Hits@10"),
            r2_fit_h10=(alt.get("train_fit") or {}).get("Hits@10"),
            r0_gap=ref.get("generalization_gap_MRR"),
            r2_gap=alt.get("generalization_gap_MRR"),
            n_fit=(ref.get("train_fit") or {}).get("n_train_fit_sampled"),
            r0_loss_last=(ref.get("epoch_losses") or [None])[-1],
            r2_loss_last=(alt.get("epoch_losses") or [None])[-1],
            r0_train_seconds=ref.get("train_seconds"),
            r2_train_seconds=alt.get("train_seconds"),
        )
        rec["test_drop_pct"] = pct(rec["r2_test"], rec["r0_test"])
        rec["fit_drop_pct"] = pct(rec["r2_fit"], rec["r0_fit"])
        rec["test_retention"] = ratio(rec["r2_test"], rec["r0_test"])
        rec["fit_retention"] = ratio(rec["r2_fit"], rec["r0_fit"])
        if rec["fit_retention"] is None or rec["test_retention"] is None:
            rec["verdict"] = "incomplete"
        else:
            sep = rec["fit_retention"] - rec["test_retention"]
            rec["separation"] = sep
            if rec["fit_retention"] >= min_fit_retention and sep >= min_separation:
                rec["verdict"] = "generalization"   # objection refuted for this model
            elif rec["fit_retention"] < min_fit_retention and sep < min_separation:
                rec["verdict"] = "optimization"     # objection UPHELD for this model
            else:
                rec["verdict"] = "mixed"
        recs.append(rec)
    return recs


def loss_curve_rows(runs):
    """Tidy (model, regime, epoch, loss) rows for a side-by-side R0/R2 loss figure."""
    rows = []
    for r in runs:
        for i, v in enumerate(r.get("epoch_losses") or [], 1):
            rows.append((r["model"], r["regime"], i, v))
    return rows


def loss_summary(runs, every=50):
    """Downsampled loss trajectory, R0 next to R2, per model."""
    by = {}
    for r in runs:
        by.setdefault(r["model"], {})[r["regime"]] = r.get("epoch_losses") or []
    out = {}
    for model, d in by.items():
        a, b = d.get(REF_REGIME, []), d.get(ALT_REGIME, [])
        n = max(len(a), len(b))
        eps = [e for e in range(1, n + 1) if e == 1 or e % every == 0 or e == n]
        out[model] = [(e,
                       a[e - 1] if e <= len(a) else None,
                       b[e - 1] if e <= len(b) else None) for e in eps]
    return out


# --------------------------------------------------------------------------- #
# Phase B: the matched sweep
# --------------------------------------------------------------------------- #
def cfg_label(r):
    lr = r.get("lr")
    return (f"d{r.get('dim')} e{r.get('epochs')} "
            f"lr={'default' if lr is None else f'{lr:g}'} neg={r.get('num_negs')}")


def analyse_sweep(runs, min_sweep_gap):
    """Pair every sweep config across R0 and R2, then test the paper's claim.

    THE CLAIM, stated so it can fail:
        "The BEST R2 configuration still performs far below the BEST R0 configuration,
         therefore the R2 collapse is not a hyperparameter artifact."

    It fails if the best R2 config substantially closes the gap to R0. The gap at the
    per-regime optimum is the honest test -- comparing R2's best against R0's DEFAULT
    would stack the deck, because the sweep is allowed to help R0 too.
    """
    cells = {}
    for r in runs:
        cells.setdefault(cfg_label(r), {})[r["regime"]] = r

    rows = []
    for label in sorted(cells):
        ref, alt = cells[label].get(REF_REGIME), cells[label].get(ALT_REGIME)
        rows.append(dict(
            config=label,
            r0_mrr=(ref or {}).get("MRR"), r2_mrr=(alt or {}).get("MRR"),
            r0_h10=(ref or {}).get("Hits@10"), r2_h10=(alt or {}).get("Hits@10"),
            r0_auroc=(ref or {}).get("AUROC_type"), r2_auroc=(alt or {}).get("AUROC_type"),
            delta_pct=pct((alt or {}).get("MRR"), (ref or {}).get("MRR")),
            paired=bool(ref and alt),
        ))

    have_r0 = [r for r in rows if r["r0_mrr"] is not None]
    have_r2 = [r for r in rows if r["r2_mrr"] is not None]
    best_r0 = max(have_r0, key=lambda r: r["r0_mrr"]) if have_r0 else None
    best_r2 = max(have_r2, key=lambda r: r["r2_mrr"]) if have_r2 else None
    worst_r0 = min(have_r0, key=lambda r: r["r0_mrr"]) if have_r0 else None

    verdict = dict(
        n_configs=len(rows), n_paired=sum(r["paired"] for r in rows),
        best_r0_config=best_r0["config"] if best_r0 else None,
        best_r0_mrr=best_r0["r0_mrr"] if best_r0 else None,
        best_r2_config=best_r2["config"] if best_r2 else None,
        best_r2_mrr=best_r2["r2_mrr"] if best_r2 else None,
    )
    if best_r0 and best_r2:
        verdict["best_vs_best_gap_pct"] = pct(best_r2["r2_mrr"], best_r0["r0_mrr"])
        verdict["best_r2_retention"] = ratio(best_r2["r2_mrr"], best_r0["r0_mrr"])
        # The objection's PREMISE is that a config chosen on R0 is the wrong config for R2.
        # If both regimes peak at the same cell, that premise is false on this grid: the
        # R0-selected config is also R2's own best, so R2 is not being handicapped by
        # having inherited it.
        verdict["same_optimum"] = bool(best_r0["config"] == best_r2["config"])
        # Widest and narrowest gap across the grid -- if R2 trails R0 in EVERY cell, the
        # deficit is not located at any particular corner of the config space.
        paired = [r for r in rows if r["paired"] and r["delta_pct"] is not None]
        if paired:
            verdict["n_cells_r2_below_r0"] = sum(r["delta_pct"] < 0 for r in paired)
            verdict["n_cells_paired"] = len(paired)
            verdict["min_gap_pct"] = max(r["delta_pct"] for r in paired)
            verdict["max_gap_pct"] = min(r["delta_pct"] for r in paired)
        # Sanity framing: is the best R2 config even above the WORST R0 config?
        if worst_r0:
            verdict["best_r2_beats_worst_r0"] = bool(best_r2["r2_mrr"] > worst_r0["r0_mrr"])
            verdict["worst_r0_mrr"] = worst_r0["r0_mrr"]
        gap = -(verdict["best_vs_best_gap_pct"] or 0.0)
        verdict["claim_survives"] = bool(gap >= min_sweep_gap)
        # How much of the gap does retuning R2 recover, relative to R2's own worst config?
        r2_vals = [r["r2_mrr"] for r in have_r2]
        verdict["r2_spread_mrr"] = max(r2_vals) - min(r2_vals) if r2_vals else None
        r0_vals = [r["r0_mrr"] for r in have_r0]
        verdict["r0_spread_mrr"] = max(r0_vals) - min(r0_vals) if r0_vals else None
    return rows, verdict


def headline_context(path: Path):
    """R0 vs R2 MRR from the frozen headline summary, for context and cross-checking."""
    try:
        agg = json.loads(Path(path).read_text())["aggregate"]
    except Exception:
        return {}
    out = {}
    for entry in agg.values():
        if entry.get("dim") == 64 and entry.get("epochs") == 300:
            out.setdefault(entry["model"], {})[entry["regime"]] = entry
    ctx = {}
    for model, d in out.items():
        ref, alt = d.get(REF_REGIME), d.get(ALT_REGIME)
        if ref and alt:
            ctx[model] = dict(r0=ref.get("MRR_mean"), r2=alt.get("MRR_mean"),
                              r0_sd=ref.get("MRR_sd"), r2_sd=alt.get("MRR_sd"),
                              drop_pct=pct(alt.get("MRR_mean"), ref.get("MRR_mean")),
                              n_seeds=ref.get("n_seeds"))
    return ctx


# --------------------------------------------------------------------------- #
# Report
# --------------------------------------------------------------------------- #
def build_md(fit_recs, loss_rows, sweep_rows, sweep_verdict, ctx, args, fit_runs):
    md = ["# FINDINGS: is the R2 collapse an optimization artifact?", "",
          "Generated by `scripts/r2_sweep_record.py` from the raw run JSONs under",
          "`data/processed/results/kge_r2_sweep/`. Produced by `scripts/run_r2_sweep.ps1`.", ""]

    md += ["## The objection", "",
           "R2 is a degree-preserving, type-preserving rewired null of the R0 training graph",
           "(`scripts/build_degree_null.py`). KGE performance collapses there, and the paper",
           "reads that as evidence the topological signal was degree. The objection is that",
           "this is confounded: **the hyperparameters were selected on R0**, so the collapse",
           "could be the recipe failing to OPTIMIZE on a rewired graph rather than the",
           "rewired graph containing no learnable structure.", "",
           "Two existing facts bear on this but do not settle it:", "",
           "* The five non-KGE baselines are **parameter-free heuristics**. Their R2 collapse",
           "  (-43% to -49% overall, -84% to -86% in the rarest-disease quartile) has no",
           "  hyperparameters to blame and therefore cannot be a tuning artifact at all. This",
           "  is strong, but it does not reach the KGE arm, which is where the objection bites.",
           "* PreferentialAttachment is flat at R2. This proves the **rewiring preserved",
           "  degree** -- a necessary check on the null's construction -- but it is *not* an",
           "  answer to this objection: PrefAttach is untrained, so it has no optimization",
           "  that could have failed. It must not be presented as if it were.", "",
           "Everything below is about the KGE arm.", ""]

    # ---------------- Phase A ---------------- #
    md += ["## Phase A -- train-fit diagnostic", "",
           "Each trained model scores ~%d sampled **training** edges of the target relation"
           % (args.train_fit_n or 4000),
           "through the *same* evaluator (`lib_eval.rank_test_edges`, sampled-50-negative,",
           "filtered, type-matched) that produces its test numbers. Nothing differs between",
           "the two columns except whether the model saw the edge during training.", "",
           "Train fit upper-bounds what the optimizer achieved on the graph it was actually",
           "given, so the two explanations make opposite predictions:", "",
           "| if ... | then train fit on R2 | reading |",
           "|---|---|---|",
           "| the recipe cannot optimize the rewired graph | collapses with test fit | objection UPHELD |",
           "| the optimizer succeeded, R2 just has nothing transferable | stays near R0's | objection REFUTED |",
           ""]

    hdr = ["model", "R0 test MRR", "R2 test MRR", "test drop",
           "R0 train-fit MRR", "R2 train-fit MRR", "train-fit drop", "signature"]
    rows = []
    for r in fit_recs:
        if not r.get("complete"):
            rows.append([r["model"], "--", "--", "--", "--", "--", "--",
                         f"incomplete (have {r.get('have')})"])
            continue
        rows.append([r["model"], fmt(r["r0_test"]), fmt(r["r2_test"]), fmtp(r["test_drop_pct"]),
                     fmt(r["r0_fit"]), fmt(r["r2_fit"]), fmtp(r["fit_drop_pct"]),
                     r["verdict"]])
    md += [md_table(hdr, rows), ""]

    md += ["Retention ratios (R2 / R0) -- the separation between them IS the measurement:", ""]
    hdr = ["model", "test retention", "train-fit retention", "separation",
           f"fit retention >= {args.min_fit_retention}", f"separation >= {args.min_separation}"]
    rows = []
    for r in fit_recs:
        if not r.get("complete") or r.get("fit_retention") is None:
            continue
        rows.append([r["model"], fmt(r["test_retention"], 3), fmt(r["fit_retention"], 3),
                     fmt(r.get("separation"), 3),
                     "yes" if r["fit_retention"] >= args.min_fit_retention else "NO",
                     "yes" if r.get("separation", -1) >= args.min_separation else "NO"])
    md += [md_table(hdr, rows), ""]

    # loss trajectories
    md += ["### Training-loss trajectory, R0 vs R2", "",
           "Full per-epoch curves are in `loss_curves.csv` (columns: model, regime, epoch,",
           "loss) and in the `epoch_losses` field of every run JSON. Downsampled here:", ""]
    for model, pairs in sorted(loss_summary(fit_runs, every=args.loss_every).items()):
        md += [f"**{model}**", ""]
        md += [md_table(["epoch", "R0 loss", "R2 loss"],
                        [[e, fmt(a, 4), fmt(b, 4)] for e, a, b in pairs]), ""]

    # ---------------- Phase B ---------------- #
    md += ["## Phase B -- matched hyperparameter sweep", "",
           "TransE swept over learning rate x negatives per positive, **the same grid on R2",
           "and on R0**. The R0 arm is what makes this an argument: a sweep on R2 alone",
           "cannot show the collapse is config-independent, because it has no reference for",
           "what the same grid buys on an un-rewired graph.", "",
           "**These runs are 100 epochs, not the headline 300.** The sweep is asked to",
           "establish the config *ordering* -- whether any configuration closes the R0-R2",
           "gap -- not to produce headline numbers, and it should not be cited as one.",
           "",
           "The shortened budget turns out not to cost much, which is checkable rather than",
           "assumed: lr=1e-3/neg=16 is exactly the headline recipe, so that cell differs from",
           "the frozen 300-epoch run only in epoch count. It lands at 0.7988 vs 0.7975 on R0",
           "and 0.6101 vs 0.6102 on R2 -- i.e. TransE has effectively converged by 100 epochs",
           "on BOTH graphs. The ordering below therefore rests on converged runs, and the R2",
           "arm in particular is not being penalised by a short budget.", "",
           "**The dim axis {64, 128} was cut** to stay inside the compute budget (adding it",
           "costs ~3.2 h more on this GPU). Consequence, stated plainly: this sweep varies",
           "OPTIMIZATION hyperparameters, not model CAPACITY, so on its own it cannot answer",
           "\"would a wider model fit R2?\". Phase A does reach that question -- a capacity",
           "shortfall would appear as depressed R2 *train* fit, which Phase A measures",
           "directly.", ""]

    hdr = ["config", "R0 MRR", "R2 MRR", "R2 vs R0", "R0 Hits@10", "R2 Hits@10",
           "R0 AUROC(type)", "R2 AUROC(type)"]
    rows = [[r["config"], fmt(r["r0_mrr"]), fmt(r["r2_mrr"]), fmtp(r["delta_pct"]),
             fmt(r["r0_h10"]), fmt(r["r2_h10"]), fmt(r["r0_auroc"], 3), fmt(r["r2_auroc"], 3)]
            for r in sweep_rows]
    md += [md_table(hdr, rows), ""]

    v = sweep_verdict
    md += ["### The claim, tested", "",
           "> \"The BEST R2 configuration in the sweep still performs far below the BEST R0",
           "> configuration, therefore the R2 collapse is not a hyperparameter artifact.\"", ""]
    if v.get("best_r0_mrr") is not None and v.get("best_r2_mrr") is not None:
        md += [f"* Best R0 config: `{v['best_r0_config']}` -> MRR **{fmt(v['best_r0_mrr'])}**",
               f"* Best R2 config: `{v['best_r2_config']}` -> MRR **{fmt(v['best_r2_mrr'])}**",
               f"* Best-vs-best gap: **{fmtp(v['best_vs_best_gap_pct'])}** "
               f"(R2 retains {fmt(v.get('best_r2_retention'), 3)} of R0)",
               f"* MRR spread across configs: R0 {fmt(v.get('r0_spread_mrr'))}, "
               f"R2 {fmt(v.get('r2_spread_mrr'))} "
               f"(how much retuning moves each regime at all)"]
        if v.get("best_r2_beats_worst_r0") is not None:
            md += [f"* Best R2 config vs **worst** R0 config ({fmt(v.get('worst_r0_mrr'))}): "
                   + ("best R2 is still higher" if v["best_r2_beats_worst_r0"]
                      else "**best R2 does not even reach the worst R0 config**")]
        if v.get("n_cells_paired"):
            md += [f"* R2 trails R0 in **{v['n_cells_r2_below_r0']} of {v['n_cells_paired']}** "
                   f"paired cells; gap ranges {fmtp(v.get('min_gap_pct'))} (narrowest) to "
                   f"{fmtp(v.get('max_gap_pct'))} (widest)"]
        if v.get("same_optimum") is not None:
            md += ["* **Both regimes peak at the same configuration** "
                   f"(`{v['best_r0_config']}`), so the config selected on R0 is also R2's "
                   "own best -- R2 is not handicapped by having inherited R0's tuning. "
                   "This contradicts the objection's premise directly."
                   if v["same_optimum"] else
                   f"* The regimes peak at DIFFERENT configurations (R0 `{v['best_r0_config']}`, "
                   f"R2 `{v['best_r2_config']}`), so R2 genuinely wanted different tuning; "
                   "the gap above is measured after giving it that."]
        md += ["",
               (f"**Claim SURVIVES** at the stated threshold (gap >= {args.min_sweep_gap}%): "
                "no configuration in the grid closes the R2-R0 gap."
                if v.get("claim_survives") else
                f"**CLAIM IN TROUBLE.** The best R2 config comes within "
                f"{fmt(abs(v['best_vs_best_gap_pct'] or 0), 1)}% of the best R0 config, "
                f"below the {args.min_sweep_gap}% threshold set in advance. Retuning "
                "substantially closes the gap, and the headline reading of R2 cannot stand "
                "as written."), ""]
    else:
        md += ["*(sweep incomplete -- no verdict)*", ""]

    # ---------------- Overall ---------------- #
    md += ["## Verdict", ""]
    sigs = [r["verdict"] for r in fit_recs if r.get("complete")]
    fit_ok = bool(sigs) and all(s == "generalization" for s in sigs)
    fit_bad = bool(sigs) and any(s == "optimization" for s in sigs)
    sweep_ok = v.get("claim_survives")

    if fit_ok and sweep_ok:
        md += ["**The objection is closed.**", "",
               "Phase A: every model fits its R2 training edges nearly as well as its R0",
               "training edges while its held-out R2 metrics collapse. Optimization on the",
               "rewired graph succeeded; what fails is generalization from it. Phase B: no",
               "configuration in a matched lr x negatives grid closes the gap, and the same",
               "grid on R0 shows the grid itself was capable of moving performance.", "",
               "The R2 collapse is therefore a property of the rewired graph, not of the",
               "hyperparameters chosen on R0."]
    elif fit_bad:
        md += ["**The objection is UPHELD, at least in part.**", "",
               "Phase A shows R2 train fit collapsing along with R2 test fit for at least one",
               "model, which is the signature of an optimization failure on the rewired graph",
               "rather than of an absent transferable signal. The R2 result cannot be read as",
               "a clean 'no structure beyond degree' claim for the KGE arm without further",
               "work (per-regime retuning to convergence at minimum).", "",
               "This is the outcome the diagnostic was built to be able to return, and it is",
               "reported here rather than set aside."]
    elif fit_ok and sweep_ok is False:
        md += ["**Split result -- report both halves.**", "",
               "Phase A refutes the optimization explanation (R2 fits its own training edges).",
               "But Phase B's best R2 configuration substantially closes the gap to R0, so the",
               "*strength* of the headline collapse is config-sensitive even though its",
               "*direction* is not. The headline number should be reported against a",
               "per-regime-tuned R2, not the R0-tuned one."]
    else:
        md += ["*(incomplete -- rerun `scripts/run_r2_sweep.ps1` and regenerate)*", ""]
        md += [f"Phase A signatures: {sigs or 'none'}; Phase B claim_survives={sweep_ok}."]
    md += [""]

    if ctx:
        md += ["## Context: the frozen headline runs (3 seeds, 300 epochs)", "",
               "Unchanged by this analysis; shown so the Phase A rerun can be checked against",
               "it. Phase A uses seed 42 only, so it should land inside these bands.", ""]
        md += [md_table(["model", "R0 MRR (3 seeds)", "R2 MRR (3 seeds)", "R2 vs R0"],
                        [[m, f"{fmt(d['r0'])} +/- {fmt(d.get('r0_sd'), 4)}",
                          f"{fmt(d['r2'])} +/- {fmt(d.get('r2_sd'), 4)}", fmtp(d["drop_pct"])]
                         for m, d in sorted(ctx.items())]), ""]

    md += ["## Reproduce", "",
           "```bash",
           "powershell -NoProfile -ExecutionPolicy Bypass -File scripts/run_r2_sweep.ps1",
           "python scripts/r2_sweep_record.py",
           "```", ""]
    return "\n".join(md) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fit-dir", default=str(R2_DIR / "trainfit"),
                    help="Phase A run JSONs (train-fit diagnostic)")
    ap.add_argument("--sweep-dir", default=str(R2_DIR / "sweep"),
                    help="Phase B run JSONs (matched grid)")
    ap.add_argument("--out-dir", default=str(R2_DIR))
    ap.add_argument("--findings", default="FINDINGS_r2_sweep.md",
                    help="markdown report path (repo-root relative by default)")
    ap.add_argument("--headline", default=str(RESULTS / "kge" / "kge_summary.json"),
                    help="frozen headline summary, for context only (never written)")
    ap.add_argument("--train-fit-n", type=int, default=4000,
                    help="documentation only: the --train-fit-n used by the runs")
    ap.add_argument("--loss-every", type=int, default=50,
                    help="downsample the loss trajectory table to every Nth epoch")
    # Verdict thresholds, fixed in advance so the report can return "objection upheld".
    ap.add_argument("--min-fit-retention", type=float, default=0.85,
                    help="R2/R0 train-fit MRR ratio above which optimization is judged to "
                         "have SUCCEEDED on the rewired graph")
    ap.add_argument("--min-separation", type=float, default=0.15,
                    help="how far train-fit retention must exceed test retention before the "
                         "collapse is judged a generalization rather than optimization failure")
    ap.add_argument("--min-sweep-gap", type=float, default=15.0,
                    help="percent by which the best R2 config must still trail the best R0 "
                         "config for the 'not a hyperparameter artifact' claim to survive")
    args = ap.parse_args()

    fit_runs = load_runs(Path(args.fit_dir))
    sweep_runs = load_runs(Path(args.sweep_dir))
    print(f"loaded {len(fit_runs)} train-fit run(s) from {args.fit_dir}")
    print(f"loaded {len(sweep_runs)} sweep run(s) from {args.sweep_dir}")

    fit_recs = analyse_trainfit(fit_runs, args.min_fit_retention, args.min_separation)
    sweep_rows, sweep_verdict = analyse_sweep(sweep_runs, args.min_sweep_gap)
    ctx = headline_context(Path(args.headline))

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    payload = {
        "thresholds": {"min_fit_retention": args.min_fit_retention,
                       "min_separation": args.min_separation,
                       "min_sweep_gap_pct": args.min_sweep_gap},
        "train_fit": fit_recs,
        "sweep": sweep_rows,
        "sweep_verdict": sweep_verdict,
        "headline_context": ctx,
        "n_fit_runs": len(fit_runs), "n_sweep_runs": len(sweep_runs),
    }
    (out_dir / "r2_sweep_findings.json").write_text(json.dumps(payload, indent=2))

    rows = loss_curve_rows(fit_runs)
    with (out_dir / "loss_curves.csv").open("w", encoding="utf8", newline="") as fh:
        fh.write("model,regime,epoch,loss\n")
        for m, rg, e, v in rows:
            fh.write(f"{m},{rg},{e},{v}\n")

    md = build_md(fit_recs, rows, sweep_rows, sweep_verdict, ctx, args, fit_runs)
    Path(args.findings).write_text(md, encoding="utf8")

    print(f"wrote {out_dir / 'r2_sweep_findings.json'}")
    print(f"wrote {out_dir / 'loss_curves.csv'} ({len(rows)} epoch rows)")
    print(f"wrote {args.findings}")
    for r in fit_recs:
        if r.get("complete"):
            print(f"  [A] {r['model']:9s} test {fmtp(r['test_drop_pct'])} | "
                  f"train-fit {fmtp(r['fit_drop_pct'])} -> {r['verdict']}")
    if sweep_verdict.get("best_r2_mrr") is not None:
        print(f"  [B] best R0 {fmt(sweep_verdict['best_r0_mrr'])} vs best R2 "
              f"{fmt(sweep_verdict['best_r2_mrr'])} "
              f"({fmtp(sweep_verdict['best_vs_best_gap_pct'])}) -> "
              f"claim_survives={sweep_verdict.get('claim_survives')}")


if __name__ == "__main__":
    main()
