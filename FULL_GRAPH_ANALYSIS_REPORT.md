# Full-Graph Network Analysis — Results

*Tier-1, Task 1 of the publication roadmap: re-run all network analyses on the complete integrated knowledge graph (440,171 nodes / 5,767,668 edges) instead of the 606-node subset.*
*Generated June 18, 2026.*

---

## ⏩ 2026-06-19 addendum — re-run on the INTEGRATED graph with networkx + scipy

The analysis below was computed on the Monarch-only graph with hand-written numpy.
It has now been re-run on the **integrated graph** (Monarch + DGIdb + Orphadata +
Gene2Phenotype: 448,931 nodes / 5,537,413 undirected edges) using networkx + scipy
(`scripts/full_graph_analysis.py`), which also yields the **exact Fiedler value**
the sandbox could only bound. Outputs: `full_graph_analysis_integrated/`.

| Metric | Monarch-only (numpy, below) | **Integrated (networkx+scipy)** |
|---|---|---|
| Nodes / edges | 440,171 / 5.77M | 448,931 / 5.54M undirected |
| Top hub | GO:0005515 (20,641) | GO:0005515 (20,641) — unchanged |
| Giant component | 97.35% | 97.48% |
| Mean local clustering (sampled) | 0.280 | 0.171 |
| Largest Laplacian eigenvalue | 20,642 | 20,642.006 |
| Fiedler (algebraic connectivity) | bounded only | **0.037904 (exact)** |
| Communities (Louvain) | 11,471 (Q=0.457) | 5,162 (Q=0.487) |

Notes: (1) the lower sampled clustering reflects the newly added drug nodes
(CHEMBL/RXCUI/DRUGBANK), which are low-degree leaves attached to few genes;
(2) the community count differs partly because of the algorithm change
(networkx Louvain vs. the hand-written numpy version) — both show strong modular
structure (Q ≈ 0.46–0.49), which the null-model test confirms is far above chance
(`NETWORK_NULL_MODEL_REPORT.md`, clustering z ≈ +230). The exact Fiedler value
(0.038) supersedes the old subset estimate (≈0.014) and the bounded full-graph note.

---

## 1. What was done

The original community / centrality results (`science_fair_data_report.txt`) were computed on `consensus similar table.csv` — only **606 nodes and 1,299 edges (0.02 % of the graph)**. This pass re-runs every network-science analysis on the **entire** graph read directly from `data/processed/relationships.csv`.

Because the work ran in a sandbox without `scipy`, `networkx`, or `igraph` (and no package installation), every algorithm was re-implemented in `numpy`: CSR adjacency, PageRank, eigenvector centrality, connected components (union-find), Louvain modularity optimization, clustering, and sampled betweenness. All result tables and figures are in the `full_graph_analysis/` folder and the four `fullgraph_*.png` files.

A note on integrity: an initial run used a vectorized scatter to build the adjacency, which silently dropped edges; a clustering sanity check (a coefficient > 1, which is impossible) caught it, and **all adjacency-dependent metrics were recomputed** with a correct sort-based CSR. This is exactly the kind of self-check the paper will need.

---

## 2. Graph scale and structure (full graph)

| Metric | Subset (old report) | **Full graph (this run)** |
|---|---|---|
| Nodes | 606 | **440,171** |
| Edges | 1,299 | **5,494,276** unique undirected (5,767,668 directed; 2,417 self-loops removed) |
| Mean degree | ~4 | **24.96** |
| Max degree | — | **20,641** |
| Mean local clustering | — | **0.280** (sampled, n=3,000) |
| Connected components | — | **5,220** |
| Largest component (LCC) | — | **428,520 nodes = 97.4 %** |
| 2nd-largest component | — | 25 nodes |
| Isolated nodes | — | 288 |

The single most important structural fact: **97.4 % of all nodes lie in one connected component**, and the next-largest component has only 25 nodes. The integrated graph is one giant, densely interconnected core — not a collection of separate disease modules. Mean local clustering of 0.28 indicates substantial triadic closure (genes/phenotypes that share neighbors tend to be linked), typical of a real biological network.

### Node composition (identifier namespaces)
Mouse (MGI, 130,621) and zebrafish (ZFIN, 74,507) genes dominate, followed by human genes (HGNC 33,971; NCBIGene 16,625), Xenopus (28,529), Gene Ontology terms (22,586), worm (WB 16,840), ClinVar variants (15,150) and MONDO diseases (13,942). The graph is heavily **multi-organism** — a strength for orthology-based inference but something the paper must control for when making human-disease claims.

---

## 3. Centrality — who the hubs really are

On the full graph the top hubs are **generic annotation terms**, not the MONDO disease IDs the subset surfaced. This is both a sanity check (these are exactly the terms you'd expect to dominate) and a finding.

**Top hubs by degree**

| Rank | Node | Meaning | Degree |
|---|---|---|---|
| 1 | GO:0005515 | "protein binding" (most generic GO term) | 20,641 |
| 2 | UBERON:0000473 | testis | 18,172 |
| 3 | GO:0003674 | molecular_function (GO root) | 14,510 |
| 4 | GO:0046872 | metal ion binding | 13,553 |
| 5 | GO:0008150 | biological_process (GO root) | 13,074 |
| 6 | HGNC:10856 | (gene) | 8,775 |
| 7 | GO:0007165 | signal transduction | 8,192 |
| 8 | ZFA:0001094 | zebrafish anatomy term | 8,081 |
| 9 | UBERON:0002107 | liver | 7,943 |
| 10 | HP:0000007 | autosomal recessive inheritance | 7,557 |

PageRank gives a nearly identical ranking (testis and the generic GO terms on top). **Eigenvector ("spectral") centrality reproduces the original `spectral_integrity_report.txt` top-10 targets exactly** — GO:0005515, EMAPA:16894, EMAPA:17373, UBERON:0000992, UBERON:0000473, HGNC:31859, EMAPA:17544, EMAPA:16728, EMAPA:16846, EMAPA:17577 — which independently validates that this re-implementation matches the project's earlier `scipy`-based spectral computation.

**Bridge nodes (approximate betweenness).** Betweenness was estimated by sampled Brandes (k = 21 source nodes; exact betweenness on 5.5 M edges is not feasible in pure Python). The top bridge node is **HGNC:1101 = BRCA2**, followed by zebrafish/anatomy terms and several HGNC disease genes — biologically sensible bridges. Treat this ranking as **indicative only**; for the paper, recompute exact (or properly sampled) betweenness with a compiled library.

---

## 4. Community structure — the headline correction

This is where the full-graph result differs most from the subset.

- The subset reported **57 communities**.
- On the full graph, Louvain modularity optimization (resolution = 1, single level, converged in 14 local-moving passes) finds **11,471 communities** with **modularity Q = 0.457** — strong, genuine modular structure.
- Community sizes are highly skewed: the largest holds 94,009 nodes (21 %), the next 82,421; **201 communities have ≥ 100 nodes and 1,892 have ≥ 10 nodes**; 290 are singletons.

So the integrated graph is modular (Q ≈ 0.46 is a healthy value for a biological network), but with a few very large communities and a long tail of small ones. The subset's "57 clean communities" was an artifact of analyzing a tiny, near-bipartite slice; the real community landscape is far richer and should be characterized (e.g. annotate the top communities by disease area / organism) in the manuscript.

*Methodological note:* a true multi-level Louvain/Leiden run (with aggregation between levels) would likely push Q higher still and merge the long tail. The single-level result here is a lower bound on the achievable modularity and was chosen to fit the compute environment.

---

## 5. Spectral / connectivity analysis

- Largest Laplacian eigenvalue λ_max = **20,642** (equal to the max degree, as expected — a correctness check).
- The LCC is a single connected component, so the algebraic connectivity (Fiedler value) is strictly positive. The project's earlier full-graph `scipy` run reported **Fiedler = 0.0137** ("resilient / high redundancy"), consistent with the one-giant-core picture found here.
- An exact Fiedler value from this pure-`numpy` run would require a shift-invert sparse eigensolver (`scipy.sparse.linalg.eigsh(sigma=0)`); the Lanczos bottom-spectrum did not converge in the available iterations, so the connectivity conclusion is stated from the exact component structure plus the prior `scipy` value rather than re-estimated here.

---

## 6. Deliverables in the folder

Result tables — `full_graph_analysis/`:
`full_degree_top.csv`, `full_pagerank_top.csv`, `full_eigenvector_top.csv` (top-200 each), `full_betweenness_approx.csv`, `full_node_prefix_composition.csv`, `full_components.json`, `full_communities.json`, `full_spectral.json`, `full_spectral_report.txt`, `full_clustering.json`.

Figures (project root):
`fullgraph_degree_distribution.png` (heavy-tailed/scale-free-like), `fullgraph_top_hubs.png`, `fullgraph_community_sizes.png`, `fullgraph_node_composition.png`.

---

## 7. What this changes for the paper

1. **The community claim is now defensible.** Report Q = 0.457 over 11,471 communities on the full graph, not "57 communities" on 0.02 % of the data.
2. **Hubs are generic-annotation and anatomy terms** (protein binding, GO roots, testis/liver, inheritance modes). This matters: link-prediction models can score these trivially, so they should likely be down-weighted or excluded when ranking novel gene–disease predictions.
3. **One connected 97 % core + clustering 0.28 + Fiedler ≈ 0.014** together support a coherent "resilient, densely integrated" description of the resource.
4. **Multi-organism composition** (mouse/zebrafish/Xenopus dominate) must be handled explicitly before making human-disease discovery claims.

### Recommended next steps (to make these results publication-grade)
- Re-run in a full scientific stack (`igraph`/`scipy`/`networkx`) to get exact betweenness, multi-level Leiden communities, and the exact Fiedler value with confidence.
- Compare community structure against a degree-preserving null model (configuration model) to show Q = 0.457 is significant.
- Annotate the largest communities by predominant disease/organism/namespace and report what they represent biologically.
- Strip the ultra-generic hubs (GO:0005515, GO roots) for a "biologically informative" variant of the centrality and prediction analyses.

The substantive conclusion of Tier-1 Task 1: **all network analyses now run on the complete 440K-node / 5.5M-edge graph**, the results are internally consistent and corroborate the project's earlier full-graph spectral computation, and the community structure — the one analysis that had been run on the tiny subset — is now properly characterized on the real graph.
