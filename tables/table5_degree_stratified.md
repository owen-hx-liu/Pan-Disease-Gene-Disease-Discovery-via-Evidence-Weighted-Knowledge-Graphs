### Table 5. Degree-stratified performance across method classes

*MRR mean over 3 seeds (42, 1, 7); 50 type-matched sampled negatives; recomputed from stored per-edge ranks (no retraining). Δ disease / Δ gene = Q4−Q1 lift along each degree axis; R2 drop = loss under the degree-null within that disease quartile.*

*Disease-degree quartiles: Q1 0–40; Q2 40–134; Q3 134–374; Q4 374–1248. Gene-degree quartiles: Q1 0–97; Q2 97–234; Q3 235–425; Q4 425–4913.*

| Method | R0 Q1 (rare) | R0 Q4 (popular) | Δ disease (Q4−Q1) | Δ gene (Q4−Q1) | R2 drop Q1 | R2 drop Q4 |
|---|---|---|---|---|---|---|
| Random | 0.090 | 0.090 | +0.001 | +0.003 | +0% | +0% |
| CommonNeighbors | 0.442 | 0.533 | +0.091 | +0.454 | -86% | -10% |
| AdamicAdar | 0.465 | 0.546 | +0.082 | +0.476 | -85% | -8% |
| Jaccard | 0.443 | 0.424 | -0.019 | +0.397 | -84% | -17% |
| PreferentialAttachment | 0.061 | 0.961 | +0.900 | -0.154 | +1% | +0% |
| TransE | 0.683 | 0.864 | +0.180 | -0.044 | -62% | -7% |
| RotatE | 0.682 | 0.916 | +0.234 | -0.067 | -62% | -3% |
