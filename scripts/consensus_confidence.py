#!/usr/bin/env python3
"""
consensus_confidence.py -- Principled, independence-aware edge-confidence score.

Fixes PROJECT_REPORT.md §6.5 / §6.6: the old `consensus similar table.csv` was
degenerate (34% of rows shared one score), had malformed double-prefix IDs
(HGNC:HGNC:...), collapsed everything to a generic CONSENSUS_SIMILAR relation, and
treated non-independent sources as independent votes.

New model (noisy-OR over independent sources):
  Each source s has a reliability r_s = P(edge true | s asserts it).
  confidence(edge) = 1 - Π_{s in srcs'} (1 - r_s)
  where srcs' collapses aggregators: Monarch / DisGeNET / DGIdb aggregate other
  sources, so they are counted as a SINGLE aggregator vote rather than treated as
  independent of the primary sources they redistribute (see provenance table).

This yields a continuous score that rewards genuine corroboration by independent
*primary* sources, and is reported with a sensitivity analysis over r_s so the
ranking is shown to be robust to the exact weights.

Outputs (data/processed/consensus/):
  edge_confidence_summary.json   distribution + support-level counts
  top_corroborated_edges.csv     multi-source edges, highest confidence first
"""
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from config import PROCESSED_DIR
from provenance_table import db_of, AGGREGATORS   # reuse name canonicalization

EDGES = PROCESSED_DIR / "edges_clean_integrated.csv"
OUT_DIR = PROCESSED_DIR / "consensus"

# per-source reliability priors (defensible defaults; sensitivity-tested below)
RELIABILITY = {
    "Monarch": 0.80,          # large curated aggregator
    "Gene2Phenotype": 0.85,   # expert-curated gene-disease
    "Orphadata": 0.80,        # expert-curated rare disease
    "CiVIC": 0.85,            # expert-curated clinical cancer evidence
    "DGIdb": 0.70,            # aggregated drug-gene
    "DrugCentral": 0.75,      # curated drug-target bioactivity
    "GWAS": 0.65,             # statistical SNP->gene/trait associations
    "DIDA": 0.75,             # curated digenic interactions
    "_default": 0.70,
}
AGG = set(AGGREGATORS)        # aggregator DB names


def sources_of(token_str):
    return {db_of(t) for t in str(token_str).split(";") if t.strip()}


def independent_sources(srcs):
    """Collapse aggregators to one vote: keep all primaries + at most one aggregator."""
    primaries = {s for s in srcs if s not in AGG}
    aggs = {s for s in srcs if s in AGG}
    out = set(primaries)
    if aggs:
        out.add("__AGG__")     # single combined aggregator vote
    return out, primaries, aggs


def confidence(srcs, rel):
    indep, _p, aggs = independent_sources(srcs)
    p_true = 1.0
    for s in indep:
        if s == "__AGG__":
            # combined aggregator reliability = max of present aggregators
            r = max(rel.get(a, rel["_default"]) for a in aggs)
        else:
            r = rel.get(s, rel["_default"])
        p_true *= (1.0 - r)
    return 1.0 - p_true


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    conf_hist = Counter()                 # rounded confidence -> count
    support_hist = Counter()              # # independent sources -> count
    multi = []                            # rows with >1 source token
    n = 0
    # sensitivity: alternative reliability settings
    rel_lo = {k: max(0.5, v - 0.1) for k, v in RELIABILITY.items()}
    rel_hi = {k: min(0.95, v + 0.1) for k, v in RELIABILITY.items()}
    sens_base, sens_lo, sens_hi = [], [], []

    for chunk in pd.read_csv(EDGES, dtype=str, chunksize=500_000):
        for s, r, t, src in zip(chunk["source_id"], chunk["relation"],
                                chunk["target_id"], chunk["dataset_sources"]):
            n += 1
            srcs = sources_of(src)
            indep, primaries, aggs = independent_sources(srcs)
            c = confidence(srcs, RELIABILITY)
            conf_hist[round(c, 4)] += 1
            support_hist[len(indep)] += 1
            if len(srcs) > 1:
                row = (s, r, t, ";".join(sorted(srcs)), len(indep), round(c, 4))
                multi.append(row)
                sens_base.append(c)
                sens_lo.append(confidence(srcs, rel_lo))
                sens_hi.append(confidence(srcs, rel_hi))

    # sensitivity: rank stability of multi-source edges under different weights
    spear_lo = float(spearmanr(sens_base, sens_lo).correlation) if len(sens_base) > 2 else None
    spear_hi = float(spearmanr(sens_base, sens_hi).correlation) if len(sens_base) > 2 else None

    md = pd.DataFrame(multi, columns=["source_id", "relation", "target_id",
                                      "dataset_sources", "indep_support", "confidence"])
    md = md.sort_values(["confidence", "indep_support"], ascending=False)
    md.to_csv(OUT_DIR / "top_corroborated_edges.csv", index=False)

    summary = {
        "total_edges": n,
        "distinct_confidence_values": len(conf_hist),
        "max_single_value_share": max(conf_hist.values()) / n,
        "support_level_counts": dict(sorted(support_hist.items())),
        "multi_source_edges": len(multi),
        "reliability_priors": RELIABILITY,
        "sensitivity_spearman_low_weights": spear_lo,
        "sensitivity_spearman_high_weights": spear_hi,
        "confidence_quantiles": {
            q: float(np.quantile(np.repeat(list(conf_hist.keys()),
                                           list(conf_hist.values())), q))
            for q in (0.5, 0.9, 0.99)
        },
    }
    json.dump(summary, open(OUT_DIR / "edge_confidence_summary.json", "w"), indent=2)
    print(json.dumps(summary, indent=2))
    print(f"\nwrote {OUT_DIR/'edge_confidence_summary.json'} and top_corroborated_edges.csv")
    print(f"top corroborated edges: {len(md)}")


if __name__ == "__main__":
    main()
