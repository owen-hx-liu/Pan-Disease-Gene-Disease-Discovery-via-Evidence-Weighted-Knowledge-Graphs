# Biological validation of novel gene–disease predictions

**Date:** 2026-06-19 · Method: Adamic-Adar (the benchmark-winning predictor) ·
Input: `edges_clean_integrated.csv` · Output:
`data/processed/predictions/adamic_adar_gene_disease.csv` (1,000 ranked,
filtered novel HGNC-gene → MONDO-disease hypotheses) · Script:
`scripts/predict_adamic_adar.py`

Predictions are filtered: no self-loops, no already-known edges, gene→disease type
constraint, and generic hubs (degree > 2,000) excluded as shared-neighbor evidence.

## Two layers of validation

### 1. Quantitative — held-out recovery (primary)
Under the proper held-out, filtered protocol (`LINK_PREDICTION_BENCHMARK.md`),
Adamic-Adar recovers held-out true edges with **MRR 0.689, Hits@1 0.59,
Hits@10 0.82** — i.e. for a held-out gene–disease edge it ranks the correct
partner first ~59% of the time out of 448,931 candidates. This is the rigorous,
quantitative evidence that the method recovers real associations.

### 2. Face validity — top novel predictions are established biology
The top-ranked *novel* (not-yet-an-edge) predictions are bona-fide gene–disease
associations confirmable from the medical literature:

| Rank | Gene | Disease | Known? |
|---|---|---|---|
| 1 | GCK | monogenic diabetes | ✔ GCK = MODY2 |
| 2 | HNF1A | monogenic diabetes | ✔ HNF1A = MODY3 |
| 3 | HNF4A | monogenic diabetes | ✔ HNF4A = MODY1 |
| 4 | RPGR | RPGR-related retinopathy | ✔ |
| 5 | PTEN | PTEN hamartoma tumor syndrome | ✔ |
| 6 | RPE65 | RPE65-related retinopathy | ✔ |
| 7 | DYSF | limb-girdle muscular dystrophy | ✔ dysferlinopathy |
| 8 | IDUA | mucopolysaccharidosis type I | ✔ Hurler |
| 10 | ATM | ATM-related cancer predisposition | ✔ |
| 12 | MYH7 | hypertrophic cardiomyopathy | ✔ |

The three top hits forming a clean **MODY (monogenic diabetes) cluster**
(GCK/HNF1A/HNF4A) is a strong signal: the method independently re-discovers a
coherent, correct gene family. That these are flagged "novel" means the *specific*
`gene –GENE_ASSOCIATED_WITH_CONDITION→ disease` triple was absent from the graph —
so Adamic-Adar is filling genuine, verifiable gaps in the KG.

## Honest caveats (state in the paper)
- **Some top hits are near-tautological / eponymous** (e.g. RPGR → "RPGR-related
  retinopathy"): the disease is defined by the gene, so recovery is easy. These
  inflate face validity; report them but don't over-claim.
- A few targets are **obsolete MONDO terms** (e.g. "obsolete Glanzmann's
  thrombasthenia") — the integrated graph carries some deprecated disease nodes;
  MONDO-obsolete filtering is a useful cleanup step before a final candidate list.
- **This is retrospective face validity, not prospective validation.** The
  strongest test is to hold out a *newer database release* (or a source withheld
  from training) and measure precision@k on edges that appear only there. The
  ingestion of DGIdb/Orphadata/Gene2Phenotype makes a leave-one-source-out test
  feasible as the next step.
- Predictions come from **Adamic-Adar, not the KGE models**, because KGE ranking
  was an order of magnitude worse (MRR ≤ 0.09). Basing hypotheses on the weaker
  method would be indefensible.
