### Table S1. Full-ranking robustness of the regime effects

*Filtered full ranking against the whole disease pool (D per regime). Panel A: KGE arm, mean ± SD over 3 seeds (1, 7, 42). Panel B: topological arm, one deterministic run per regime. Stricter analogue of the sampled-negative audit in Table 2: the R2 degree-null collapse reproduces and deepens, R3 lowers no method's score, and the arm-specific R1 effect is larger here than under sampled negatives (absolute MRR is lower against the ~14k-disease pool).*

| Model | Regime | Pool D | MRR | Hits@1 | Hits@3 | Hits@10 | AUROC (type) | ΔMRR vs R0 |
|---|---|---|---|---|---|---|---|---|
| **Panel A. Knowledge-graph embedding** |  |  |  |  |  |  |  |  |
| TransE | R0 standard | 14,307 | 0.079 ± 0.001 | 0.036 ± 0.001 | 0.072 ± 0.002 | 0.150 ± 0.003 | 0.962 ± 0.002 | ref |
| TransE | R1 redundancy | 14,252 | 0.071 ± 0.003 | 0.031 ± 0.002 | 0.066 ± 0.002 | 0.138 ± 0.008 | 0.954 ± 0.001 | -10% |
| TransE | R2 degree-null | 14,307 | 0.023 ± 0.002 | 0.005 ± 0.001 | 0.014 ± 0.003 | 0.042 ± 0.004 | 0.879 ± 0.001 | -71% |
| TransE | R3 orthology-blocked | 14,137 | 0.080 ± 0.004 | 0.034 ± 0.003 | 0.077 ± 0.005 | 0.155 ± 0.008 | 0.963 ± 0.000 | +2% |
| RotatE | R0 standard | 14,307 | 0.108 ± 0.004 | 0.053 ± 0.002 | 0.107 ± 0.003 | 0.202 ± 0.003 | 0.967 ± 0.001 | ref |
| RotatE | R1 redundancy | 14,252 | 0.088 ± 0.002 | 0.037 ± 0.002 | 0.086 ± 0.002 | 0.178 ± 0.002 | 0.955 ± 0.001 | -19% |
| RotatE | R2 degree-null | 14,307 | 0.032 ± 0.002 | 0.006 ± 0.001 | 0.021 ± 0.001 | 0.067 ± 0.007 | 0.889 ± 0.000 | -70% |
| RotatE | R3 orthology-blocked | 14,137 | 0.111 ± 0.003 | 0.055 ± 0.003 | 0.110 ± 0.003 | 0.211 ± 0.004 | 0.969 ± 0.001 | +3% |
| DistMult | R0 standard | 14,307 | 0.050 ± 0.003 | 0.020 ± 0.002 | 0.045 ± 0.004 | 0.098 ± 0.005 | 0.937 ± 0.002 | ref |
| DistMult | R1 redundancy | 14,252 | 0.039 ± 0.002 | 0.014 ± 0.002 | 0.032 ± 0.004 | 0.076 ± 0.006 | 0.897 ± 0.002 | -22% |
| DistMult | R2 degree-null | 14,307 | 0.011 ± 0.001 | 0.001 ± 0.001 | 0.006 ± 0.001 | 0.020 ± 0.002 | 0.832 ± 0.002 | -78% |
| DistMult | R3 orthology-blocked | 14,137 | 0.058 ± 0.003 | 0.026 ± 0.003 | 0.050 ± 0.003 | 0.108 ± 0.002 | 0.943 ± 0.001 | +16% |
| ComplEx | R0 standard | 14,307 | 0.055 ± 0.006 | 0.020 ± 0.005 | 0.046 ± 0.008 | 0.114 ± 0.008 | 0.918 ± 0.004 | ref |
| ComplEx | R1 redundancy | 14,252 | 0.037 ± 0.002 | 0.012 ± 0.001 | 0.028 ± 0.001 | 0.072 ± 0.005 | 0.891 ± 0.004 | -32% |
| ComplEx | R2 degree-null | 14,307 | 0.010 ± 0.001 | 0.001 ± 0.000 | 0.004 ± 0.000 | 0.017 ± 0.003 | 0.799 ± 0.006 | -82% |
| ComplEx | R3 orthology-blocked | 14,137 | 0.059 ± 0.006 | 0.022 ± 0.004 | 0.052 ± 0.006 | 0.118 ± 0.014 | 0.917 ± 0.009 | +8% |
| **Panel B. Topological baselines** |  |  |  |  |  |  |  |  |
| Random | R0 standard | 14,307 | 0.001 | 0.000 | 0.000 | 0.001 | — | ref |
| Random | R1 redundancy | 14,252 | 0.001 | 0.000 | 0.000 | 0.001 | — | +0% |
| Random | R2 degree-null | 14,307 | 0.001 | 0.000 | 0.000 | 0.001 | — | +0% |
| Random | R3 orthology-blocked | 14,137 | 0.001 | 0.000 | 0.000 | 0.001 | — | +1% |
| Common-Neighbors | R0 standard | 14,307 | 0.118 | 0.076 | 0.125 | 0.192 | — | ref |
| Common-Neighbors | R1 redundancy | 14,252 | 0.116 | 0.074 | 0.123 | 0.191 | — | -2% |
| Common-Neighbors | R2 degree-null | 14,307 | 0.014 | 0.001 | 0.007 | 0.049 | — | -88% |
| Common-Neighbors | R3 orthology-blocked | 14,137 | 0.117 | 0.075 | 0.123 | 0.191 | — | -1% |
| Adamic-Adar | R0 standard | 14,307 | 0.133 | 0.091 | 0.143 | 0.214 | — | ref |
| Adamic-Adar | R1 redundancy | 14,252 | 0.130 | 0.088 | 0.140 | 0.211 | — | -2% |
| Adamic-Adar | R2 degree-null | 14,307 | 0.014 | 0.000 | 0.008 | 0.048 | — | -89% |
| Adamic-Adar | R3 orthology-blocked | 14,137 | 0.132 | 0.090 | 0.143 | 0.213 | — | -1% |
| Jaccard | R0 standard | 14,307 | 0.106 | 0.079 | 0.113 | 0.150 | — | ref |
| Jaccard | R1 redundancy | 14,252 | 0.104 | 0.077 | 0.112 | 0.150 | — | -2% |
| Jaccard | R2 degree-null | 14,307 | 0.009 | 0.003 | 0.005 | 0.020 | — | -91% |
| Jaccard | R3 orthology-blocked | 14,137 | 0.106 | 0.079 | 0.113 | 0.151 | — | -0% |
| Preferential-Attachment | R0 standard | 14,307 | 0.025 | 0.000 | 0.000 | 0.081 | — | ref |
| Preferential-Attachment | R1 redundancy | 14,252 | 0.025 | 0.000 | 0.000 | 0.081 | — | -0% |
| Preferential-Attachment | R2 degree-null | 14,307 | 0.025 | 0.000 | 0.000 | 0.081 | — | +1% |
| Preferential-Attachment | R3 orthology-blocked | 14,137 | 0.025 | 0.000 | 0.000 | 0.080 | — | -1% |

*Panel A ΔMRR vs R0: R1 -32% (ComplEx), -22% (DistMult), -19% (RotatE), -10% (TransE); R2 -82% (ComplEx), -78% (DistMult), -70% (RotatE), -71% (TransE); R3 +8% (ComplEx), +16% (DistMult), +3% (RotatE), +2% (TransE). R2 costs every embedding model 70-82% of its MRR, and R1 (redundancy) is the regime that separates the two arms: it costs the KGE models 10-32% of MRR against only 2% for the overlap heuristics, which is the arm-specific redundancy leak of the sampled-negative audit reappearing at larger magnitude under the stricter protocol. R3 lowers no model's MRR under this protocol either, but the size of the gain is not uniform: the translational pair stays within +3% while the bilinear pair rises 8-16%, so blocking the cross-species bridges is mildly helpful to the bilinear models rather than merely uninformative -- a stronger refutation of the orthology-leak hypothesis than a flat null, and in the opposite direction from leakage.*

*Panel B is a single deterministic run per regime (the topological scorers carry no seed; the Random control is the analytic chance floor), so no SD is reported, and AUROC is not defined for the ranking-only protocol. n=4,228 test edges; ranks were verified against a brute-force loop of the same scorers over the full pool on 150 sampled edges per regime (exact agreement). Under full ranking the R0→R2 collapse of the overlap heuristics deepens to 88-91% (it is 43-49% under 50 sampled negatives), while Preferential-Attachment moves by +1%: the degree-only control is unmoved by the degree-preserving null under the strictest protocol, so the sampled-negative effect sizes are lower bounds.*

*30 cold-start test genes (absent from the training vocabulary) receive the worst possible rank.*
