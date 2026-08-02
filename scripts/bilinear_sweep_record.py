#!/usr/bin/env python3
"""bilinear_sweep_record.py -- consolidate the DistMult/ComplEx training sweep into
one auditable table (data/processed/results/kge/bilinear_sweep.{json,md}).

Why this exists
---------------
The manuscript previously reported DistMult and ComplEx as configuration-specific
failures. That claim is only defensible if the configurations actually tried are
enumerated. This script collects the single-variable sweep that identified the cause
(PyKEEN's per-model default regularizer) and writes it out in a form a reviewer can
read, so that whichever model ends up trained and whichever does not, the evidence
behind the statement is on disk rather than in a chat log.

Sources
-------
  * the diagnostic JSONs produced by the sweep runner (loss trajectory, pos-vs-neg
    score gap, embedding norms), passed via --diag,
  * any completed run JSONs under data/processed/results/kge_smoke*/ (MRR at the
    smoke epoch budget).

Usage:
    python scripts/bilinear_sweep_record.py --diag a.json b.json [--out-dir DIR]
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


def load_diag(paths):
    rows = []
    for p in paths:
        try:
            recs = json.loads(Path(p).read_text())
        except Exception as e:  # noqa: BLE001
            print(f"  skip {p}: {e}")
            continue
        for r in recs:
            if "error" in r:
                rows.append(dict(config=r["name"], outcome="error", detail=r["error"]))
                continue
            b, a = r["before"], r["after"]
            ep = r["per_epoch"]
            # "trained" = the loss left its initialization value AND the positive/negative
            # score gap opened up. Both are required: a loss that drifts while the gap
            # stays at zero is the degenerate solution reparameterized, not learning.
            moved = len(ep) > 1 and abs(ep[-1] - ep[0]) > 1e-3
            opened = abs(a["gap"]) > 1e-2
            rows.append(dict(
                config=r["name"],
                loss_first=ep[0] if ep else None,
                loss_last=ep[-1] if ep else None,
                gap_before=round(b["gap"], 4), gap_after=round(a["gap"], 4),
                rel_norm_before=round(b["rel_target_norm"], 3),
                rel_norm_after=round(a["rel_target_norm"], 3),
                ent_norm_after=round(a["ent_norm"], 3),
                epochs=len(ep),
                outcome="trains" if (moved and opened) else "collapses",
            ))
    return rows


def load_smokes(results_root: Path):
    out = []
    for d in sorted(results_root.glob("kge_smoke*")):
        for p in sorted(d.glob("kge_*_seed*.json")):
            try:
                r = json.loads(p.read_text())
            except Exception:
                continue
            if "model" not in r:
                continue
            # Runs predating the --regularizer flag have no such key and used PyKEEN's
            # per-model default; runs that set it explicitly store None meaning DISABLED.
            # Those two are opposite conditions, so absence must not render as "none".
            # The default differs per model (Lp for the bilinear models, none for TransE
            # and RotatE), so label it generically rather than asserting "Lp".
            reg = ("PyKEEN per-model default" if "regularizer" not in r
                   else ("none (disabled)" if r["regularizer"] is None else r["regularizer"]))
            out.append(dict(
                source=d.name, model=r["model"], regime=r["regime"], epochs=r["epochs"],
                loss=r.get("loss"), num_negs=r.get("num_negs"), loop=r.get("loop"),
                batch=r.get("train_batch"), regularizer=reg,
                lr=r.get("lr"), MRR=round(r["MRR"], 4),
                Hits10=round(r["Hits@10"], 4), AUROC_type=round(r["AUROC_type"], 4),
                train_seconds=r.get("train_seconds"),
            ))
    return out


def md_table(headers, rows):
    out = ["| " + " | ".join(headers) + " |",
           "|" + "|".join("---" for _ in headers) + "|"]
    for r in rows:
        out.append("| " + " | ".join("" if v is None else str(v) for v in r) + " |")
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--diag", nargs="*", default=[],
                    help="diagnostic result JSONs (loss/gap/norm sweep)")
    ap.add_argument("--out-dir", default=str(RESULTS / "kge"))
    ap.add_argument("--results-root", default=str(RESULTS))
    args = ap.parse_args()

    diag = load_diag(args.diag)
    smokes = load_smokes(Path(args.results_root))

    payload = {"sweep": diag, "smoke_runs": smokes,
               "criterion": "trains = training loss leaves its initialization value AND "
                            "the mean positive-minus-negative score gap exceeds 0.01; "
                            "collapses = loss pinned at its zero-gap value "
                            "(softplus -> ln 2 = 0.6931; NSSA margin 9 -> 4.5; "
                            "margin-ranking 1.0 -> 1.0) with the gap at zero."}
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "bilinear_sweep.json").write_text(json.dumps(payload, indent=2))

    md = ["# Bilinear KGE training sweep (DistMult, ComplEx)", "",
          payload["criterion"], "",
          "## Single-variable sweep (short runs; loss trajectory and score gap)", ""]
    if diag:
        # No pipe characters in headers: they would split the markdown columns.
        hdr = ["config", "epochs", "loss first", "loss last", "gap after",
               "rel norm before", "rel norm after", "outcome"]
        md.append(md_table(hdr, [[r.get("config"), r.get("epochs"), r.get("loss_first"),
                                  r.get("loss_last"), r.get("gap_after"),
                                  r.get("rel_norm_before"), r.get("rel_norm_after"),
                                  r.get("outcome")] for r in diag]))
    md += ["", "## Smoke runs scored by the shared harness (lib_eval)", ""]
    if smokes:
        hdr = ["source", "model", "regime", "epochs", "loss", "negs", "batch",
               "regularizer", "MRR", "Hits@10", "AUROC_type"]
        md.append(md_table(hdr, [[r["source"], r["model"], r["regime"], r["epochs"],
                                  r["loss"], r["num_negs"], r["batch"], r["regularizer"],
                                  r["MRR"], r["Hits10"], r["AUROC_type"]] for r in smokes]))
    (out_dir / "bilinear_sweep.md").write_text("\n".join(md) + "\n")
    print(f"wrote {out_dir/'bilinear_sweep.json'} ({len(diag)} sweep rows, {len(smokes)} smoke runs)")
    print(f"wrote {out_dir/'bilinear_sweep.md'}")


if __name__ == "__main__":
    main()
