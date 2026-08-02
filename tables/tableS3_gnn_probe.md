### Table S3. The message-passing (R-GCN) probe

*Matched 5% subgraph of the R0 and R2 training graphs (292,604 of 5,852,083 edges, all 33,832 target edges retained in both), d64, 300 epochs, seeds 42/1/7, sampled 50 type-matched negatives. Absolute values are NOT comparable to Table 2 — different graph, different negative pool, 8.4% of test edges touching an out-of-vocabulary entity against 30 cold-start genes on the full graph.*

| Row | R0 MRR | R2 MRR | ΔMRR R0→R2 | Δ% | R0 Hits@10 | R2 Hits@10 | Retained of full-graph collapse | Validity gate |
|---|---|---|---|---|---|---|---|---|
| **Panel A. Matched 5% subgraph, d64, 300 epochs, 3 seeds (mean ± SD)** |  |  |  |  |  |  |  |  |
| R-GCN (message-passing encoder) | 0.459 ± 0.007 | 0.341 ± 0.009 | -0.1179 | -25.7% | 0.707 ± 0.007 | 0.655 ± 0.012 | no full-graph counterpart | inherits FAIL |
| DistMult (matched decoder control) | 0.600 ± 0.003 | 0.573 ± 0.008 | -0.0263 | -4.4% | 0.804 ± 0.009 | 0.793 ± 0.001 | 12% | FAIL |
| **Panel B. Subgraph-validity calibration — the reference decoder alone, R0 vs R2 by subgraph fraction (seed 42)** |  |  |  |  |  |  |  |  |
| frac 0.05 (292,604 train edges) | 0.5997 | 0.5734 | -0.0263 | -4.4% | — | — | 12% | FAIL |
| frac 0.11 (643,729 train edges) | 0.6333 | 0.5735 | -0.0598 | -9.4% | — | — | 26% | FAIL |
| frac 0.25 (1,463,021 train edges) | 0.6411 | 0.5087 | -0.1324 | -20.7% | — | — | 57% | PASS |
| frac 0.50 (2,926,042 train edges) | 0.6454 | 0.4624 | -0.1830 | -28.4% | — | — | 78% | PASS |
| frac 1.00 (5,852,083 train edges) | 0.6689 | 0.4259 | -0.2430 | -36.3% | — | — | 100% | definitional |

*Convergence gate: PASS on all 3 R0 seeds (training loss fell, MRR ≥ 0.35 against a chance floor of 0.089, Hits@10 ≥ 0.5, and MRR ≥ 0.6× the matched reference), so Panel A's delta does not come from an undertrained model.*

*Subgraph-validity gate: FAIL, and this governs how Panel A may be read. On the frac-0.05 subgraph the reference decoder retains only 12% of its known full-graph R0→R2 collapse (-4.4% here against -36.3% on the full graph), below the threshold of 50% fixed in advance. Subsampling dissolves the very leak under study, so no R0→R2 delta measured on this graph -- the R-GCN's included -- transfers to the full-graph claim. Panel A's -25.7% must not be read as degree leakage, nor placed beside the frozen KGE band of Table 2, which it happens to fall inside.*

*Panel B asks whether any fraction is both trainable and valid, and the answer is no. The retained collapse rises monotonically with subgraph size while R0 is nearly flat from frac 0.11 upward, so the subsample was never damaging the task, only the leak -- which also shows what the R-GCN needed was more graph, not more epochs. The gate first clears at frac 0.25, and a direct probe measured 25.2 minutes per R-GCN epoch there (about 31 days for the six production cells on the GPU used here), while every fraction the R-GCN can train on within budget fails the gate. The wall is memory-pressure thrashing expressed as time rather than an allocation failure, so there is no batch size to lower.*
