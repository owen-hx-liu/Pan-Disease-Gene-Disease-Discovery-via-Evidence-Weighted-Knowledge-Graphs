### Table 6. Worked examples of degree-driven misranking

*Rank 1 = true disease ranked above all 50 type-matched negatives; 51 = below all of them. RotatE mean over 3 seeds; recomputed from stored per-edge ranks.*

| Case | Held-out edge (gene → disease) | Gene deg | Disease deg | RotatE R0 | RotatE R2 | Pref.-Att. | Adamic-Adar |
|---|---|---|---|---|---|---|---|
| degree false-positive | LIN28B-AS1 → insomnia | 1 | 801 | 1 | 1 | 1 | 26 |
| missed real association | RS1 → retinoschisis | 350 | 1 | 51 | 48 | 48 | 28 |
| missed real association | NF1 → plexiform neurofibroma | 1,065 | 1 | 43 | 49 | 48 | 36 |
| missed real association | KIT → mucosal melanoma | 1,387 | 1 | 49 | 45 | 48 | 30 |
| missed real association | ONECUT1 → neonatal diabetes mellitus | 277 | 1 | 50 | 50 | 46 | 26 |

*For the insomnia hub (degree 801; the held-out target of 86 genes) the true disease is ranked #1 for 97% of those genes by pure-degree Preferential-Attachment, 85% by RotatE, and only 43% by shared-neighbour Adamic-Adar — the degree predictor matches or beats the KGE, and RotatE's mean MRR here barely moves under the degree-null (0.932 to 0.882).*
