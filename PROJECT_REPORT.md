# Biomedical Knowledge Graph & Link-Prediction Project — Comprehensive Report

*Prepared June 18, 2026 · Status: pre-publication review*

---

## ⏩ 2026-06-18 status update (addendum — original assessment below is preserved)

Substantial Tier-1/Tier-2 work was completed this day (full environment with
GPU + PyKEEN, not the restricted sandbox). The original sections below are left
unchanged for the record; current status:

| Roadmap item | Status |
|---|---|
| 6.1 Analysis on full graph, not 606-node subset | **Done** — all analyses run on the full/integrated graph |
| 6.2 Real train/test split (no train==test) | **Done** — `kge_benchmark.py`, 80/10/10 held-out, `TransE.py` evaluation removed |
| 6.3 Metrics + baselines | **Done** — TransE/DistMult/ComplEx/RotatE vs Random/CN/Adamic-Adar, MRR/Hits/AUROC/AUPRC, random + type-matched negatives |
| 6.4 Filter invalid predictions | **Done** — `novel_predictions_filtered.csv` (no self-loops/known/type-mismatch) |
| 6.5 Consensus score degeneracy | **Done** — rebuilt as independence-aware noisy-OR, sensitivity-tested (`CONSENSUS_CALIBRATION_REPORT.md`) |
| 6.6 Source non-independence | **Done** — aggregators collapsed in consensus; provenance table |
| 6.7 Secrets / absolute paths / repro | **Done** — `scripts/config.py`, `requirements.txt`, `README.md`, `.gitignore` (password was never committed) |
| 6.9 Null models for network claims | **Done** — `NETWORK_NULL_MODEL_REPORT.md` (z-scores vs configuration model) |
| **NEW blocker** — graph was 100% Monarch | **Partially fixed** — DGIdb + Orphadata + Gene2Phenotype ingested (+49K novel edges, first multi-source corroboration); remaining DBs pending. See `CRITICAL_FINDINGS.md` |

Still open: broaden ingestion to the remaining standalone DBs (CiVIC, DisGeNET,
DrugCentral, GWAS, ClinGen, LncRNADisease); biological/literature validation of
top predictions; embedding hyperparameter sweep with variance over ≥3 seeds.

---

## 1. Executive summary

This project builds a large, integrated **biomedical knowledge graph** from ~13 public data sources, resolves the entities to canonical identifiers, embeds the graph with two different representation-learning methods (node2vec and TransE), and uses those embeddings plus a "consensus" scoring scheme to **predict new gene–disease and disease–disease links** — framed as a drug-repurposing / disease-gene-discovery tool. On top of the graph it runs a battery of network-science analyses (community detection, spectral/algebraic-connectivity analysis, centrality, temporal growth) and ships several interactive HTML explorers and a Streamlit prediction interface.

The data-engineering half of the project is genuinely substantial: a clean, multi-format ingestion pipeline that produces a **440,171-node / 5,767,668-edge** graph with passing entity-resolution validation. The scientific/analysis half is currently the weak point: the headline predictive and network results are computed on a **606-node / 1,299-edge subset**, the link-prediction models have **methodological flaws that would not survive peer review** (notably no train/test split for TransE), and there is **no quantitative evaluation, no baseline comparison, and no biological validation** of the predictions. The gap between the scale of the data and the rigor of the analysis is the single biggest thing holding back publication.

This report documents the project in detail and then lays out exactly what needs to change to reach a peer-reviewed journal.

---

## 2. What the project does (scientific framing)

The core idea is a standard but powerful one in computational biology: many separate databases each describe pieces of the gene → variant → disease → phenotype → drug landscape. Individually each is incomplete and noisy. If you merge them into one graph where nodes are biological entities (genes, diseases, variants, phenotypes, anatomical structures, drugs) and edges are curated relationships (`interacts_with`, `has_phenotype`, `causes`, `treats`, etc.), you can then ask a machine-learning model to **predict edges that are not yet in the graph but are strongly implied by its structure**. Those predicted edges are hypotheses — e.g. "gene X is likely associated with disease Y" or "drug D may apply to disease Y" — that could guide research or drug repurposing.

The project implements this end to end: ingest → canonicalize → validate → embed → predict → analyze → visualize.

---

## 3. Data sources

The raw corpus lives under `data/raw/` as **90 source folders** spanning **13 distinct databases**, in CSV, TSV, JSON, XML, and FASTA formats:

| Database | What it contributes |
|---|---|
| **Monarch** (18 dumps) | Core multi-species gene–phenotype–disease associations (Biolink model) |
| **Orphadata** (26 XML/JSON dumps) | Rare-disease catalog, gene associations, inheritance |
| **CiVIC** (6) | Clinical interpretation of cancer variants |
| **DGIdb** (4) | Drug–gene interactions |
| **Gene2Phenotype** (7) | Curated gene–disease with confidence & allelic requirement |
| **DIDA** (4) | Digenic disease variants |
| **ClinGen / ClinPGx** (5) | Clinical gene–disease validity, pharmacogenomics |
| **DisGeNET** | Gene–disease associations with scores |
| **DrugCentral** | Drug–target interactions (`drug.target.interaction.tsv`) |
| **GWAS** | Genome-wide association study hits |
| **LncRNADisease / LNCipedia** | Long non-coding RNA–disease links, lncRNA sequences |

This is a strong, credible source list — comparable in spirit to established integrative resources (Monarch/KG, Hetionet, PrimeKG, the OpenBioLink benchmark). The breadth across human and model organisms (mouse/MGI, zebrafish/ZFIN, rat/RGD, fly/FB, worm/WB, Xenopus) is a real asset.

---

## 4. The pipeline, stage by stage

The work is organized as numbered scripts under `scripts/` that correspond to project "phases." Reconstructed flow:

**Phase 3 — Ingestion & extraction.** `dataset_scanner.py`, `identifier.py`, and `pipeline.py` stream every raw file regardless of format, detect which columns hold entity IDs, and emit standardized `(source, relation, target)` edges with provenance. `auto_relationship_mapper.py` / `relationship_finder_3-2.py` infer relationship types from the entity-type pair (e.g. Gene–Disease → `ASSOCIATED_WITH`). Output: timestamped `extracted_edges_*.csv`.

**Phase 3.4–3.5 — Canonicalization.** `canonical_node_creator.py` and `canonical_edge_creator.py` collapse the **695,917 raw IDs down to 458,475 canonical IDs** using identifier-priority rules (`id_priority.json`, `mondo_xrefs.tsv`) so that, e.g., the same gene referenced by HGNC, NCBIGene, and Ensembl becomes one node. `entity_vocab.py` / `entity_vocab_audit.py` build and audit the vocabulary.

**Phase 4 — Cleaning & validation.** `remove_self_loops.py`, dedup steps, and a validation harness (`step4-*.py`) produce `edges_clean.csv` and run an **entity-resolution validation report**. An early run was `FATAL` (9.4M raw-ID-leakage violations, 1,946 invalid self-loops); the later **`step_4_5_validation_report.txt` is a clean PASS** — 0 leakage, 0 self-loops, dedup and evidence-integrity passing. This iteration from FATAL → PASS is good engineering hygiene and is worth showing in the paper.

**Phase 5 — Graph assembly.** `step5-1.py` assigns every node a category from its ID prefix (gene / disease / variant / phenotype / anatomy / drug); `step5-2.py` writes the Neo4j-style `relationships.csv` (`:START_ID, :END_ID, :TYPE, weight, dataset_sources`). Final graph: **440,171 nodes, 5,767,668 edges**. The most common edge types are `INTERACTS_WITH` (1.85M), `EXPRESSED_IN` (1.26M), `HAS_PHENOTYPE` (1.07M), `ORTHOLOGOUS_TO` (511K).

**Phase 6 — Consensus scoring.** `step6-1.py` defines a confidence score combining edge support, disease-level, and gene-level evidence (weights 0.5 / 0.3 / 0.2, with thresholds like ≥2 datasets per edge). This feeds the `consensus similar table.csv`.

**Embedding & prediction.**
- `n2vextractor.py` → `node2vec_embeddings.csv` (structural embeddings).
- `TransE.py` (via **PyKEEN**) → `transe_embeddings.csv` / `transe_new_predictions.csv` (knowledge-graph embeddings; 64,755 candidate predictions generated).
- `predictivemodelling.py` trains a **RandomForest** link classifier on node2vec features against randomly generated negative samples; saved as `consensus_predictor_model.pkl`.
- `drugrepurposing.py` / `druganalysis.py` push DrugCentral drug→gene links into Neo4j and look for drug→disease "discovery deltas."

**Network-science analysis.**
- `louvain.py` — community detection + degree/betweenness/eigenvector centrality → `science_fair_data_report.txt` (606 nodes, 57 communities).
- `spectralanalysis.py` — Laplacian / Fiedler value (algebraic connectivity = 0.0137, "resilient") + spectral partition map.
- `temporalchart.py` — network growth 1995→2023 (`temporal_metrics.csv`, up to 440K nodes / 5.5M edges).
- `step3-1.py`…`step3-5.py`, `statssheet.py`, `networkwalker.py` — supporting stats and graph traversal.

**Interfaces & visualization.** `modelinterface.py` (Streamlit "Bio-Link Analysis Portal"), `consensus_graph_interactive.html`, `science_fair_graph_explorer.html`, and the `analysis_outputs_20260304_*` folder of publication-style figures (score distributions, relation-type shares, PageRank, PR-coverage, TransE-vs-strength plots).

---

## 5. Current results / deliverables

- A validated integrated KG: **440K nodes, 5.77M edges, 13 sources, ~20 relation types.**
- Two embedding sets (node2vec + TransE) and a trained RandomForest link predictor.
- **64,755 TransE candidate predictions** and a ranked top-20 hybrid candidate list (e.g. HGNC:10193 → MONDO:0971172 as rank 1).
- Community structure (57 communities on the analysis subset), spectral resilience metric, centrality hubs.
- A 29-year temporal growth curve of the network.
- Interactive explorers and a Streamlit demo.

These are real, presentable artifacts and make for a strong science-fair project. The sections below are about the distance between "strong science fair project" and "publishable in a peer-reviewed journal."

---

## 6. What is holding it back (critical issues)

These are ordered by how much they would matter to a reviewer.

### 6.1 The analysis runs on a tiny subset, not the real graph — **most important**
The integration produces 440K nodes / 5.77M edges, but the headline community, spectral, centrality, and "consensus" results are computed from `consensus similar table.csv`, which is **only 1,299 edges over 606 nodes** — 0.02% of the graph. The science-fair report's "57 communities" and "top hubs" therefore describe a sliver of the data, not the integrated resource the project claims to have built. A reviewer will immediately ask why the analysis abandons 99.98% of the data. Everything downstream inherits this problem.

### 6.2 No train/test split in the link-prediction evaluation — **disqualifying as written**
`TransE.py` passes `training=tf, testing=tf` — the model is **evaluated on the exact triples it trained on**. Any reported ranking metric is meaningless (guaranteed optimistic). The RandomForest in `predictivemodelling.py` does split, but its negatives are purely random node pairs, which makes the task artificially easy. There is no held-out test set, no time-sliced split, and no reported AUC/Hits@k/MRR anywhere in the outputs. Link prediction **cannot be published without a proper evaluation protocol.**

### 6.3 No quantitative evaluation or baselines
There are no metrics tables. Journals expect link-prediction work to report **MRR, Hits@1/3/10, AUROC, AUPRC** against **baselines** (random, degree/common-neighbors, and at least one standard KGE such as DistMult/ComplEx/RotatE). Right now the predictions are presented as ranked lists with no evidence they beat chance.

### 6.4 Predictions are not biologically validated
The top candidates (and the 64K TransE predictions) are never checked against held-out edges, newer database releases, or literature. Some predicted edges are even **self-loops** (e.g. `MONDO:0006877 → MONDO:0006877` appears in `transe_new_predictions.csv`), which signals the prediction step isn't filtering trivially invalid outputs. Without validation, the central claim ("we discover new gene–disease links") is unsupported.

### 6.5 The "consensus" score is under-justified and partly degenerate
The score in `step6-1.py` uses hand-picked weights (0.5/0.3/0.2) with no sensitivity analysis or calibration. In `consensus similar table.csv` many edges share **identical scores** (e.g. `3.2316589935…` repeated across many rows), suggesting the score is dominated by a single term rather than meaningfully ranking edges. The relation is collapsed to a single generic `CONSENSUS_SIMILAR` type, discarding the rich Biolink semantics the pipeline worked to preserve.

### 6.6 Data-leakage / overlap between databases not addressed
Many of these sources are **not independent** — Monarch aggregates several of the others, DisGeNET and Gene2Phenotype overlap, etc. Treating "supported by N datasets" as N independent votes overstates confidence. The paper needs to account for source dependence (provenance dedup, or weighting by genuinely independent primary sources).

### 6.7 Reproducibility and methods gaps
No README, no environment lockfile, no fixed random seeds across all scripts, hard-coded absolute Windows paths (`C:\Users\owenh\…`) and a hard-coded Neo4j password in `drugrepurposing.py` (remove before any public release). Source database **versions and download dates are not recorded**, which journals increasingly require. Embedding hyperparameters (TransE dim=32, 100 epochs, single seed) are untuned and undocumented as choices.

### 6.8 No clear novelty statement vs. existing resources
Integrated biomedical KGs already exist (Hetionet, PrimeKG, Monarch KG, OpenBioLink, BioKG). The project needs an explicit answer to "what is new here?" — e.g. a novel consensus-weighting method, a specific disease area where it outperforms, or a validated set of new predictions. Re-deriving an existing-style KG without a differentiator won't clear review.

### 6.9 Statistical rigor of the network analyses
"Resilient" from a Fiedler value, the centrality hubs, and the community count are reported as facts without **null models** (e.g. degree-preserving randomization), confidence intervals, or significance tests. Network claims need a comparison against randomized graphs to be meaningful.

---

## 7. How to improve it — roadmap to publication

### Tier 1 — Must fix before submitting anywhere
1. **Run all analysis on the full graph** (or a principled, documented filtered version), not the 606-node subset. If memory is the constraint, use the full edge list with sparse/streaming graph libraries (already partly in place with `networkx` + `scipy.sparse`); consider `igraph`/`graph-tool` for scale.
2. **Implement a real evaluation protocol** for link prediction: hold out a test set (ideally a **time-based split** using the temporal data you already have — train on pre-2020, test on 2020–2023 edges), and never evaluate on training triples.
3. **Report standard metrics with baselines**: MRR, Hits@{1,3,10}, AUROC, AUPRC vs. random, common-neighbors/Adamic-Adar, and ≥1 modern KGE (DistMult/ComplEx/RotatE) alongside TransE. PyKEEN gives you all of these for free.
4. **Filter invalid predictions** (self-loops, existing edges, type-incompatible pairs) before ranking.
5. **Strip secrets and absolute paths**; add config + relative paths.

### Tier 2 — Makes it a credible paper
6. **Biological validation of top predictions**: check predicted gene–disease links against a newer database release or literature (PubMed); report precision@k of recovered/confirmable links. Even a handful of literature-confirmed novel predictions is a compelling result.
7. **Justify and calibrate the consensus score**: sensitivity analysis over weights, calibration against held-out true edges, and account for **source non-independence** (collapse provenance from aggregators like Monarch to primary sources).
8. **Tune and document embeddings**: small hyperparameter sweep (dimension, epochs, negative sampling), fixed seeds, report variance over ≥3 runs.
9. **Add null-model comparisons** for the network metrics (configuration-model randomization) so "resilient," hub, and community claims carry statistical weight.
10. **Record data provenance**: exact source versions, download dates, license terms, and counts per source in a supplementary table.

### Tier 3 — Strengthens novelty & polish
11. **State the contribution crisply** and position against Hetionet/PrimeKG/Monarch KG/OpenBioLink in a related-work section.
12. **Pick a focused use case** (e.g. one disease family — the project already has rich rare-disease/Orphadata coverage) and demonstrate end-to-end discovery + validation there; depth beats breadth for a first paper.
13. **Release reproducibly**: public repo, environment lockfile, seeds, a `make`/single-command pipeline, and the figures-from-data scripts (the `analysis_outputs` figures are a good start).
14. **Choose a venue and format accordingly.** Realistic targets: *Bioinformatics* (Application Note if framed as a resource/tool), *BMC Bioinformatics*, *Database (Oxford)*, *GigaScience*, *PLOS Computational Biology*, or a workshop/conference (ISMB, PSB) as a first step. An Application Note is the most attainable framing for the current resource-heavy state.

---

## 8. Suggested target narrative for the paper

> *"We integrate 13 biomedical databases into a provenance-tracked knowledge graph of 440K entities and 5.8M relationships, introduce a consensus-confidence weighting that accounts for source dependence, and benchmark several knowledge-graph embedding methods for gene–disease link prediction under a time-based evaluation split. We validate the top predictions against newer database releases / literature and recover N% of held-out associations, plus M literature-supported novel candidates."*

That framing turns the current assets (which are real and substantial) into defensible claims, and every clause maps to one of the fixes above.

---

## 9. Bottom line

The **data engineering is publication-grade in ambition and largely sound in execution** — clean multi-format ingestion, real entity resolution, a validated 5.8M-edge graph, two embedding methods, and good visualization. What's missing is the **scientific rigor layer**: run the analysis on the actual graph, evaluate predictions honestly with held-out data and baselines, validate biologically, and document everything reproducibly. None of the blockers are fatal to the project — they're fixable with the data and code already in hand. The fastest credible path is to reframe as a resource/Application Note, fix the Tier-1 evaluation issues, and add one validated discovery use case.
