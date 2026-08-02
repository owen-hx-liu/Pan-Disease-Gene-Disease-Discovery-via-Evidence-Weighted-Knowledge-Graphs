### Table S2. Hyperparameter validation of the R2 degree null

*Sampled 50 type-matched negatives, seed 42. Panel A: train fit (held-in edges) against test fit (held-out edges), 300 epochs at the headline recipe. Panel B: one matched lr × negatives grid run on both regimes, 100 epochs. Retention is R2/R0; separation is train-fit retention − test retention.*

| Model / config | R0 MRR | R2 MRR | ΔMRR R0→R2 | R0 train-fit MRR | R2 train-fit MRR | Train-fit retention | Test retention | Separation |
|---|---|---|---|---|---|---|---|---|
| **Panel A. Train-fit diagnostic (300 epochs, seed 42, headline recipe)** |  |  |  |  |  |  |  |  |
| DistMult | 0.6665 | 0.4208 | -36.9% | 0.8842 | 0.7720 | 0.873 | 0.631 | 0.242 |
| TransE | 0.7975 | 0.6090 | -23.6% | 0.9120 | 0.8565 | 0.939 | 0.764 | 0.175 |
| **Panel B. Matched hyperparameter grid, same grid on both regimes (100 epochs, seed 42)** |  |  |  |  |  |  |  |  |
| d64 e100 lr=0.001 neg=16 | 0.7988 | 0.6101 | -23.6% | — | — | — | 0.764 | — |
| d64 e100 lr=0.001 neg=64 | 0.8167 | 0.6121 | -25.0% | — | — | — | 0.750 | — |
| d64 e100 lr=0.003 neg=16 | 0.7517 | 0.5470 | -27.2% | — | — | — | 0.728 | — |
| d64 e100 lr=0.003 neg=64 | 0.7798 | 0.5647 | -27.6% | — | — | — | 0.724 | — |
| d64 e100 lr=0.01 neg=16 | 0.5840 | 0.4015 | -31.2% | — | — | — | 0.688 | — |
| d64 e100 lr=0.01 neg=64 | 0.6428 | 0.4130 | -35.8% | — | — | — | 0.642 | — |

*Panel B verdict: the best R0 configuration (d64 e100 lr=0.001 neg=64, MRR 0.8167) and the best R2 configuration (d64 e100 lr=0.001 neg=64, MRR 0.6121) are the SAME configuration, so R2 is not handicapped by having inherited tuning chosen on R0 -- which contradicts the objection's premise directly. Best-vs-best gap -25.0% (R2 retains 0.750); R2 trails R0 in 6 of 6 paired cells, with the narrowest gap anywhere in the grid at -23.6%. The grid was capable of moving performance (R0 spans 0.2327 MRR, R2 0.2106), so its failure to close the gap is informative rather than a dead grid.*

*Panel B runs 100 epochs rather than the headline 300 because it is asked for configuration ORDERING, not headline numbers, and its cells should not be quoted as results; the cost of that is measurable rather than assumed, since the lr 1e-3 / neg 16 cell is exactly the headline recipe and differs from the frozen 300-epoch seed-42 run only in epoch count -- it lands +0.2% from it on R0 and +0.2% on R2, so TransE has effectively converged by 100 epochs on both graphs and the R2 arm is not penalised by the short budget.*

*Thresholds were fixed before the runs: train-fit retention ≥ 0.85, fit-vs-test separation ≥ 0.15, and a surviving best-vs-best gap ≥ 15%.*
