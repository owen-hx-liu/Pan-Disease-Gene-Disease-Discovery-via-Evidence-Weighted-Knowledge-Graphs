# Link-Prediction Benchmark: KGE models vs. topological baselines

**Date:** 2026-06-19 · **Graph:** integrated KG `edges_clean_integrated.csv`
(448,931 nodes / 5,814,437 unique non-self-loop edges; Monarch + DGIdb + Orphadata +
Gene2Phenotype) · **Split:** random 80/10/10 (train 4,651,551 / valid 581,443 /
test 581,443), seed 42 · **KGE:** dim 128, 100 epochs, GPU (RTX 5070) ·
Script: `scripts/kge_benchmark.py`

This is the real, executed benchmark (PyKEEN + PyTorch on GPU) that the earlier
sandbox run could not perform. It fixes the original `TransE.py` (which evaluated on
its own training data, `training=tf, testing=tf`) and reports standard metrics WITH
baselines under one protocol.

## Results

Ranking = filtered, both-sides, all-entity (PyKEEN gold standard), evaluated on a
5,000-triple sample of the held-out test set. Classification = held-out positives
vs. negatives, **shared across all methods** so baselines and KGE are directly
comparable; `_type` columns use **type-matched hard negatives** (corrupt the tail
only with an entity of the same category).

| Method | MRR | Hits@1 | Hits@3 | Hits@10 | AUROC (rand) | AUPRC (rand) | AUROC (hard) | AUPRC (hard) |
|---|---|---|---|---|---|---|---|---|
| Random           | 0.092 | 0.022 | 0.063 | 0.207 | 0.500 | 0.501 | 0.498 | 0.503 |
| Common-Neighbors | 0.654 | 0.535 | 0.725 | 0.821 | 0.872 | 0.962 | 0.866 | 0.952 |
| **Adamic-Adar**  | **0.689** | **0.593** | **0.764** | **0.822** | 0.877 | 0.964 | 0.871 | 0.953 |
| TransE           | 0.044 | 0.012 | 0.050 | 0.097 | 0.957 | 0.962 | 0.890 | 0.880 |
| DistMult         | 0.002 | 0.002 | 0.002 | 0.002 | 0.526 | 0.534 | 0.521 | 0.526 |
| ComplEx          | 0.006 | 0.002 | 0.004 | 0.011 | 0.917 | 0.921 | 0.881 | 0.864 |
| RotatE           | 0.091 | 0.051 | 0.096 | 0.166 | **0.989** | **0.988** | **0.965** | **0.962** |

## The headline finding (two parts)

1. **On ranking, simple topology decisively beats all KGE models.** Adamic-Adar
   reaches MRR 0.689 / Hits@10 0.822; the best KGE (RotatE) reaches MRR 0.091 /
   Hits@10 0.166 — an order of magnitude worse at placing the true tail among all
   448,931 candidates. *Any KGE result on this graph must be reported against
   Adamic-Adar, and as configured none of them clear that bar on ranking.*

2. **Yet KGE wins the easy classification task.** RotatE/TransE/ComplEx score
   AUROC 0.88–0.99 separating true edges from random/type-matched negatives —
   *above* the baselines. This is the classic dichotomy: distinguishing a true edge
   from a random pair is easy (degree/structure alone nearly suffice), but
   **ranking** the correct partner against every candidate is the task that matters
   for discovery, and there topology dominates. Reporting only AUROC (as much of
   the KGE literature does) would have hidden this.

Random behaves exactly as theory predicts (AUROC ≈ 0.50, MRR ≈ H₅₁/51 ≈ 0.09),
validating the harness.

## Why KGE underperforms here (caveats — state these)

- **Untuned, single run.** dim=128, 100 epochs, default optimizer/loss/negative
  sampling, one seed. DistMult collapsed to ~0 (likely needs regularization /
  different loss); this is a tuning failure, not evidence DistMult is useless. A
  hyperparameter sweep with variance over ≥3 seeds is the obvious next step.
- **Strongly multi-relational + hub-dominated.** 27 relations and generic
  annotation hubs (GO:0005515 etc.) make all-entity ranking very hard for
  translational/bilinear models without filtering or relation-specific tuning.
- **AUPRC with random negatives is optimistic** (most node pairs are trivial
  non-edges); the `_type` (hard-negative) columns are the credible classification
  numbers and are uniformly lower, as expected.

## Predictions

Because KGE ranking is so weak, **novel-link hypotheses should be generated from
the winning method (Adamic-Adar), not KGE** — see
`data/processed/predictions/adamic_adar_gene_disease.csv` (filtered: no self-loops,
no known edges, gene→disease type-constrained). The KGE prediction path in
`kge_benchmark.py` (`--no-predict` was used here) remains available but is not the
basis for reported hypotheses.

## No temporal split (unchanged)

`temporalchart.py` simulates discovery years with `random.choices`; there is no
genuine edge date, so a time-based split is impossible. The split is random.
`kge_benchmark.py --year-col` enables a time split the moment real dates exist.

## Reproduce

```bash
python scripts/kge_benchmark.py \
    --edges data/processed/edges_clean_integrated.csv \
    --out data/processed/benchmark_integrated --epochs 100 --dim 128
# resumable (skips checkpointed models); --skip-kge for baselines only
```
