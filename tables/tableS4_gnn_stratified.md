### Table S4. Degree stratification of the message-passing arm

*MRR mean over 3 seeds; sampled 50 type-matched negatives; equal-count quartiles of full-graph R0 training degree (1,057 test edges per quartile), the same bins as Table 5; recomputed from stored per-edge ranks with no retraining.*

*Disease-degree quartiles: Q1 0–40; Q2 40–134; Q3 134–374; Q4 374–1248. Gene-degree quartiles: Q1 0–97; Q2 97–234; Q3 235–425; Q4 425–4913.*

| Cell | Q1 | Q2 | Q3 | Q4 | low half | high half | high/low |
|---|---|---|---|---|---|---|---|
| **Panel A. ranked-disease (target) degree axis — MRR by quartile of full-graph R0 degree** |  |  |  |  |  |  |  |
| R-GCN (message-passing encoder) — R0 | 0.256 | 0.459 | 0.571 | 0.548 | 0.357 | 0.560 | 1.57× |
| R-GCN (message-passing encoder) — R2 | 0.151 | 0.334 | 0.443 | 0.436 | 0.242 | 0.439 | 1.81× |
| DistMult (matched decoder control) — R0 | 0.233 | 0.599 | 0.784 | 0.783 | 0.416 | 0.784 | 1.88× |
| DistMult (matched decoder control) — R2 | 0.209 | 0.556 | 0.770 | 0.758 | 0.383 | 0.764 | 2.00× |
| **Panel B. query-gene (source) degree axis — MRR by quartile of full-graph R0 degree** |  |  |  |  |  |  |  |
| R-GCN (message-passing encoder) — R0 | 0.530 | 0.465 | 0.410 | 0.429 | 0.498 | 0.420 | 0.84× |
| R-GCN (message-passing encoder) — R2 | 0.356 | 0.348 | 0.323 | 0.336 | 0.352 | 0.330 | 0.94× |
| DistMult (matched decoder control) — R0 | 0.642 | 0.641 | 0.571 | 0.544 | 0.642 | 0.558 | 0.87× |
| DistMult (matched decoder control) — R2 | 0.576 | 0.614 | 0.557 | 0.547 | 0.595 | 0.552 | 0.93× |

*Under the sampled-negative protocol the true disease is ranked against type-matched negative diseases, so Panel A is the axis that exposes popularity bias, and it is the one message-passing measurement that does not route through R2: it is computed on full-graph degree bins from stored per-edge ranks, so the failed subgraph-validity gate (S3 Table) does not touch it. It shows the message-passing encoder is less disease-popularity-reliant than its own decoder (1.57× against 1.88× high/low, and 0.548 against 0.783 in the top quartile) on the identical graph, negatives, and seeds. Two qualifications belong with it. The R-GCN is the weaker model overall here (R0 MRR 0.459 against 0.600), and a weaker model has less performance to lose to popularity, so part of the gap may be level rather than architecture. And on the gene axis (Panel B) the two are nearly identical and both inverted, which retires an apparent encoder-versus-decoder contrast seen in a 40-epoch smoke run against an undertrained control.*
