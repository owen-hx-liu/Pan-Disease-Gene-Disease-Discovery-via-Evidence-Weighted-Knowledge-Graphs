# Network structure vs. a degree-preserving null model

**Date:** 2026-06-18 · **Graph:** integrated KG (`edges_clean_integrated.csv`),
448,931 nodes / 5,537,413 undirected simple edges · **Seed:** 42 ·
**Null:** 5 configuration-model realizations (same degree sequence, edges rewired
at random) · Script: `scripts/null_model_tests.py`

This addresses PROJECT_REPORT.md §6.9 / §7.9: the network claims ("resilient",
hubs, clustering, communities) were previously reported as bare numbers with **no
null model**. A configuration model preserves every node's degree but destroys all
other structure, so any metric that differs from the null is genuine structure —
not an artifact of the degree distribution.

| Metric | Observed | Null mean ± sd | z-score | Interpretation |
|---|---|---|---|---|
| Avg. local clustering (sampled, n=3000) | **0.169** | 0.0169 ± 0.0007 | **+229.6** | ~10× the null — strong, highly significant community/triadic structure |
| Degree assortativity | **−0.074** | −0.0404 ± 0.0002 | **−153.8** | significantly *more* disassortative than null — hubs preferentially attach to low-degree nodes |
| Giant-component fraction | **0.975** | 0.9922 ± 0.0002 | **−92.3** | slightly *less* connected than null — real graph retains more periphery |

**Takeaways for the paper:**

- The high clustering (z ≈ +230) is the headline: the KG has far more local
  cohesion than its degree sequence alone would produce, which is the structural
  precondition for meaningful community detection and for neighborhood-based link
  prediction (and explains why Common-Neighbors/Adamic-Adar are such strong
  baselines in the link-prediction benchmark).
- The disassortativity (z ≈ −154) is consistent with a biomedical hub-and-spoke
  topology (generic annotation hubs like GO:0005515 linking to many specific
  genes), and is significantly stronger than chance.
- The giant component is marginally *smaller* than the null expectation, so the
  earlier "resilient / highly connected" framing should be stated carefully: the
  graph is dominated by one component, but **no more so than a random graph with
  the same degrees — in fact slightly less**. Drop or qualify any claim that the
  connectivity itself is exceptional.

**Caveats.** Clustering is estimated on a 3,000-node random sample per graph
(exact transitivity at this scale is dominated by ultra-high-degree hubs and is
expensive); the z-scores are large enough that sampling error does not change the
conclusions. Modularity-based community significance (`--modularity`) was not run
here because greedy modularity on 5.5M edges across multiple realizations is slow;
the clustering contrast already establishes non-random meso-scale structure.
