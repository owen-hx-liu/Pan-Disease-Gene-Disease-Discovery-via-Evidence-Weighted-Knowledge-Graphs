### Table S1. Full-ranking robustness of the regime effects

*Filtered full ranking against the whole disease pool (D per regime); mean ± SD over 3 seed(s) (1, 7, 42). Stricter analogue of the sampled-negative audit in Table 2: the R2 degree-null collapse and the smaller R1 drop both reproduce (absolute MRR is lower against the ~14k-disease pool).*

| Model | Regime | Pool D | MRR | Hits@1 | Hits@3 | Hits@10 | AUROC (type) | ΔMRR vs R0 |
|---|---|---|---|---|---|---|---|---|
| TransE | R0 standard | 14,307 | 0.079 ± 0.001 | 0.036 ± 0.001 | 0.072 ± 0.002 | 0.150 ± 0.003 | 0.962 ± 0.002 | ref |
| TransE | R1 redundancy | 14,252 | 0.071 ± 0.003 | 0.031 ± 0.002 | 0.066 ± 0.002 | 0.138 ± 0.008 | 0.954 ± 0.001 | -10% |
| TransE | R2 degree-null | 14,307 | 0.023 ± 0.002 | 0.005 ± 0.001 | 0.014 ± 0.003 | 0.042 ± 0.004 | 0.879 ± 0.001 | -71% |

*30 cold-start test genes (absent from the training vocabulary) receive the worst possible rank.*
