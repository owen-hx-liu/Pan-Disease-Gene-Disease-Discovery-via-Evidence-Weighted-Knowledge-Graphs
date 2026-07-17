"""Worked biological case study: how degree leakage would mislead gene prioritization
(Unit: reviewer #6).

Reviewer request: "Include a short biological case study showing how benchmark leakage
could affect gene prioritization." We build one entirely from the frozen per-edge ranks
(no retraining, no hand-picked numbers): for every held-out gene the benchmark ranks the
true disease against 50 type-matched candidate diseases, and we already store that rank
per edge for RotatE/TransE (R0 and the R2 degree-null) and for the topological baselines.

Three pieces, all recomputed here and written to case_study.json:

  1. disease_degree_buckets -- RotatE MRR at R0 vs the R2 degree-null, and pure-degree
     Preferential-Attachment MRR, partitioned by the *ranked disease's* training degree.
     The pattern: on popular diseases RotatE's score is degree (barely drops under R2, and
     Preferential-Attachment alone matches or beats it); on rare diseases RotatE's score is
     genuine structure (collapses under R2) and degree is useless.

  2. hub_example -- one popular disease (insomnia, degree 801) that is the held-out target
     for many genes. A pure-degree predictor ranks it #1 more often than the KGE does, and
     the rank survives the degree-null: the "prediction" encodes that the disease is
     popular, not gene-specific biology.

  3. worked_examples -- a curated, biologically legible handful of real held-out edges
     (their numbers computed here, not typed): a degree-1 antisense RNA that the model
     confidently maps to the insomnia hub, versus canonical disease-defining associations
     (e.g. RS1->retinoschisis) that the model ranks near-last because the specific disease
     node is rare. The curation is transparent -- fixed IDs below -- and the buckets show
     the examples are representative, not cherry-picked.

Output: data/processed/results/case_study/case_study.json -- consumed by make_tables.py
(Table 6). Usage:  python scripts/case_study.py
"""
from __future__ import annotations

import csv
import json
import os
import sys

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
RESULTS = os.path.join(REPO_ROOT, "data", "processed", "results")
SPLITS = os.path.join(REPO_ROOT, "data", "processed", "splits")
OUT_DIR = os.path.join(RESULTS, "case_study")
SEEDS = [42, 1, 7]

HUB_DISEASE = "MONDO:0013600"  # insomnia (degree 801): a held-out target for many genes

# Contiguous partition of the ranked disease's training degree.
BUCKETS = [(0, 0, "cold-start (0)"), (1, 5, "rare (1–5)"), (6, 50, "uncommon (6–50)"),
           (51, 300, "common (51–300)"), (301, 10**9, "hub (301+)")]

# Curated, biologically legible held-out edges. IDs are fixed; every number is computed.
# "false_confidence": the model maps a barely-connected gene to the popular hub disease.
# "missed_rare": a canonical, disease-defining association the model ranks near-last
# because the specific disease node is sparse in the graph.
CURATED = [
    ("false_confidence", "HGNC:21553", "LIN28B-AS1", "MONDO:0013600", "insomnia"),
    ("missed_rare", "HGNC:10457", "RS1", "MONDO:0004579", "retinoschisis"),
    ("missed_rare", "HGNC:7765", "NF1", "MONDO:0003304", "plexiform neurofibroma"),
    ("missed_rare", "HGNC:6342", "KIT", "MONDO:0000544", "mucosal melanoma"),
    ("missed_rare", "HGNC:8138", "ONECUT1", "MONDO:0016391", "neonatal diabetes mellitus"),
]


def _die(msg: str) -> None:
    sys.exit(f"[case_study] {msg}")


def load_labels() -> dict[str, str]:
    path = os.path.join(REPO_ROOT, "data", "processed", "entity_vocab.csv")
    lab = {}
    with open(path, encoding="utf-8", errors="ignore") as fh:
        for row in csv.DictReader(fh):
            lab[row["id"]] = row["name"]
    return lab


def load_degrees() -> dict[str, int]:
    deg = {}
    with open(os.path.join(SPLITS, "node_degree.csv"), encoding="utf-8") as fh:
        next(fh)
        for r in csv.reader(fh):
            if len(r) >= 2 and r[0]:
                deg[r[0]] = int(r[1])
    return deg


def load_test_edges() -> tuple[list[str], list[str]]:
    genes, dis = [], []
    with open(os.path.join(SPLITS, "test.csv"), encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            genes.append(row["source_id"])
            dis.append(row["target_id"])
    return genes, dis


def kge_rr(model: str, regime: str) -> np.ndarray:
    arrs = []
    for s in SEEDS:
        p = os.path.join(RESULTS, "kge", f"kge_{model}_{regime}_seed{s}.json")
        if not os.path.exists(p):
            _die(f"missing {p}")
        arrs.append(json.load(open(p, encoding="utf-8"))["reciprocal_ranks"])
    return np.asarray(arrs, dtype=float).mean(0)


def baseline_rr(method: str, regime: str) -> np.ndarray:
    p = os.path.join(RESULTS, "baselines", f"baselines_{regime}.json")
    per = json.load(open(p, encoding="utf-8")).get("per_edge_reciprocal_ranks")
    if not per or method not in per:
        _die(f"{method}/{regime}: per-edge ranks missing (re-run run_baselines.py)")
    return np.asarray(per[method], dtype=float).mean(0)


def rank_of(rr: float) -> int:
    """Tie-averaged rank from a reciprocal rank (1 = top; 51 = below all 50 negatives)."""
    return int(round(1.0 / rr)) if rr > 0 else 51


def main() -> None:
    print("Building worked case study from frozen per-edge ranks (no retraining)...")
    lab, deg = load_labels(), load_degrees()
    genes, dis = load_test_edges()
    n = len(genes)
    dd = np.array([deg.get(d, 0) for d in dis])

    rot0, rot2 = kge_rr("RotatE", "R0"), kge_rr("RotatE", "R2")
    tr0 = kge_rr("TransE", "R0")
    pa0 = baseline_rr("PreferentialAttachment", "R0")
    aa0 = baseline_rr("AdamicAdar", "R0")

    # 1) disease-degree buckets ------------------------------------------------
    buckets = []
    for lo, hi, label in BUCKETS:
        m = np.where((dd >= lo) & (dd <= hi))[0]
        if len(m) == 0:
            continue
        buckets.append({
            "label": label, "deg_min": lo, "deg_max": (None if hi >= 10**9 else hi),
            "n_edges": int(len(m)),
            "RotatE_R0_MRR": round(float(rot0[m].mean()), 3),
            "RotatE_R2_MRR": round(float(rot2[m].mean()), 3),
            "PrefAttach_R0_MRR": round(float(pa0[m].mean()), 3),
        })

    # 2) hub disease profile ---------------------------------------------------
    hub = np.where(np.array(dis) == HUB_DISEASE)[0]
    if len(hub) == 0:
        _die(f"hub disease {HUB_DISEASE} not in test set")
    pct1 = lambda arr: round(100.0 * float(np.mean(arr[hub] >= 0.999)), 0)
    hub_example = {
        "disease_id": HUB_DISEASE, "disease_name": lab.get(HUB_DISEASE, HUB_DISEASE),
        "disease_degree": deg.get(HUB_DISEASE, 0), "n_test_genes": int(len(hub)),
        "pct_rank1_RotatE": pct1(rot0), "pct_rank1_PrefAttach": pct1(pa0),
        "pct_rank1_AdamicAdar": pct1(aa0),
        "RotatE_R0_MRR": round(float(rot0[hub].mean()), 3),
        "RotatE_R2_MRR": round(float(rot2[hub].mean()), 3),
    }

    # 3) curated worked examples ----------------------------------------------
    pos = {(genes[i], dis[i]): i for i in range(n)}
    worked = []
    for kind, gid, gdisp, did, ddisp in CURATED:
        i = pos.get((gid, did))
        if i is None:
            print(f"  WARN curated edge not in test set: {gdisp} -> {ddisp} ({gid},{did})")
            continue
        worked.append({
            "kind": kind, "gene_id": gid, "gene": gdisp, "gene_degree": int(deg.get(gid, 0)),
            "disease_id": did, "disease": ddisp, "disease_degree": int(deg.get(did, 0)),
            "RotatE_rank_R0": rank_of(rot0[i]), "RotatE_rank_R2": rank_of(rot2[i]),
            "TransE_rank_R0": rank_of(tr0[i]),
            "PrefAttach_rank_R0": rank_of(pa0[i]), "AdamicAdar_rank_R0": rank_of(aa0[i]),
        })

    out = {
        "meta": {
            "description": "Worked case study: degree-driven false confidence vs missed rare "
                           "associations, from stored per-edge ranks (sampled 50 type-matched "
                           "negatives; ranks are 1=top .. 51=below all negatives; MRR mean "
                           "over 3 seeds). No retraining.",
            "n_test_edges": n, "seeds": SEEDS, "hub_disease": HUB_DISEASE,
        },
        "disease_degree_buckets": buckets,
        "hub_example": hub_example,
        "worked_examples": worked,
    }
    os.makedirs(OUT_DIR, exist_ok=True)
    dest = os.path.join(OUT_DIR, "case_study.json")
    with open(dest, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
    print(f"  wrote {os.path.relpath(dest, REPO_ROOT)}")

    # console summary
    print("\n  disease-degree buckets (RotatE R0 -> R2 ; PrefAttach R0):")
    for b in buckets:
        print(f"    {b['label']:16s} n={b['n_edges']:4d}  RotatE {b['RotatE_R0_MRR']:.3f}"
              f" -> {b['RotatE_R2_MRR']:.3f}   PrefAttach {b['PrefAttach_R0_MRR']:.3f}")
    h = hub_example
    print(f"\n  hub: {h['disease_name']} (deg {h['disease_degree']}, {h['n_test_genes']} genes)"
          f" ranked #1 by RotatE {h['pct_rank1_RotatE']:.0f}% / PrefAttach "
          f"{h['pct_rank1_PrefAttach']:.0f}% / AdamicAdar {h['pct_rank1_AdamicAdar']:.0f}%")
    print("\n  worked examples (RotatE R0/R2, PrefAttach, AdamicAdar ranks):")
    for w in worked:
        print(f"    [{w['kind']:16s}] {w['gene']:11s}(deg {w['gene_degree']:4d}) -> "
              f"{w['disease']:26s}(deg {w['disease_degree']:3d}) | "
              f"RotatE {w['RotatE_rank_R0']:>2}/{w['RotatE_rank_R2']:>2}  "
              f"PA {w['PrefAttach_rank_R0']:>2}  AA {w['AdamicAdar_rank_R0']:>2}")


if __name__ == "__main__":
    main()
