# Consensus / edge-confidence score: rebuild & calibration

**Date:** 2026-06-18 · Script: `scripts/consensus_confidence.py` ·
Input: `data/processed/edges_clean_integrated.csv`

## What was wrong (PROJECT_REPORT §6.5–6.6)

The old `consensus similar table.csv`:
- was **degenerate** — 34% of its 1,299 rows shared the single score 3.2317 (86
  distinct values total), because every edge was Monarch-sourced with weight 1, so
  the `edge_support` term was constant and the score collapsed;
- had **malformed identifiers** with duplicated prefixes (`HGNC:HGNC:17640`,
  `MONDO:MONDO:0006882`);
- discarded all Biolink semantics into one generic `CONSENSUS_SIMILAR` relation;
- treated non-independent sources as independent votes.

## New model: independence-aware noisy-OR

For an edge asserted by source set *S* with per-source reliability
*r_s = P(edge true | s asserts it)*:

> confidence = 1 − Π_{s ∈ S′} (1 − r_s)

where **S′ collapses aggregators**: Monarch / DisGeNET / DGIdb redistribute other
databases, so they contribute a **single combined aggregator vote** rather than
counting as independent of the primary sources they carry (provenance table).
Priors: Monarch 0.80, Gene2Phenotype 0.85, Orphadata 0.80, DGIdb 0.70.

## Results (integrated graph, 5,814,437 edges)

| Independent-source support | # edges |
|---|---|
| 1 source | 5,809,622 (99.15%) |
| 2 sources | 4,779 |
| 3 sources | 36 |

- **Sensitivity:** perturbing every reliability prior by ±0.10 leaves the
  multi-source edge ranking **identical** (Spearman ρ = 1.00 for both the lower
  and higher weight settings) — the score is robust to the exact priors.
- **Top corroborated edges** (`data/processed/consensus/top_corroborated_edges.csv`):
  36 gene–disease links supported by **three independent sources**
  (Gene2Phenotype + Monarch + Orphadata, confidence 0.994), all
  `GENE_ASSOCIATED_WITH_CONDITION` / `CAUSES`. These are exactly the
  high-confidence hypotheses a consensus score should surface, and they are now
  real rather than artifacts.

## Honest limitation (ties to CRITICAL_FINDINGS Finding 1)

Because **99.15% of edges are still single-source (Monarch)**, the confidence
score only *discriminates* the ~4,851 multi-source edges; for the rest it is the
constant Monarch prior (0.80). This is no longer a bug — it is an accurate
reflection of current source coverage. The consensus score becomes broadly useful
only as more of the standalone databases are ingested (CiVIC, DisGeNET,
DrugCentral, GWAS, ClinGen, …); the three integrated so far (DGIdb, Orphadata,
Gene2Phenotype) already create the first genuine multi-source corroboration the
project has had. Treat consensus confidence as a **high-precision filter for the
corroborated subset**, not a graph-wide ranking, until coverage broadens.
