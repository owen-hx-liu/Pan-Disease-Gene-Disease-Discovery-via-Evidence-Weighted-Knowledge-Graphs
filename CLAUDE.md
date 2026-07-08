# CLAUDE.md

Guidance for AI assistants (and humans) working in this repository. Read this first.

---

## 1. Project overview

This is a **biomedical knowledge-graph (KG) and link-prediction project**. It integrates ~13 public biomedical databases into one large graph, resolves entities to canonical identifiers, embeds the graph, and predicts novel **gene–disease / disease–disease links** (framed as disease-gene discovery and drug repurposing). Network-science analyses (community detection, centrality, spectral, temporal) and interactive visualizations sit on top.

**The goal is to publish in a peer-reviewed journal.** Every decision should be evaluated against "would this survive peer review?" The data engineering is strong; the scientific-rigor layer (honest evaluation, baselines, biological validation, reproducibility) is what needs work. See `PROJECT_REPORT.md` for the full assessment and `FULL_GRAPH_ANALYSIS_REPORT.md` for the completed full-graph network analysis.

### Pipeline (phases, in `scripts/`)
1. **Ingestion/extraction** (`dataset_scanner.py`, `identifier.py`, `pipeline.py`, `auto_relationship_mapper.py`): stream raw CSV/TSV/JSON/XML/FASTA, detect ID columns, emit `(source, relation, target)` edges with provenance.
2. **Canonicalization** (`canonical_node_creator.py`, `canonical_edge_creator.py`, `entity_vocab*.py`, `id_priority.json`, `mondo_xrefs.tsv`): collapse 695,917 raw IDs → 458,475 canonical IDs.
3. **Cleaning/validation** (`remove_self_loops.py`, `step4-*.py`): dedup, drop self-loops, validate entity resolution → `edges_clean.csv`. (Validation went FATAL → PASS; see `data/processed/*validation_report*`.)
4. **Graph assembly** (`step5-1.py`, `step5-2.py`): category-per-node, write `relationships.csv` (Neo4j format).
5. **Consensus scoring** (`step6-1.py`): confidence weighting → `consensus similar table.csv`.
6. **Embedding/prediction**: `n2vextractor.py` (node2vec), `TransE.py`/`pykeen*.py` (PyKEEN TransE), `predictivemodelling.py` (RandomForest link classifier), `drugrepurposing.py`/`druganalysis.py` (Neo4j drug→disease).
7. **Analysis/UI**: `louvain.py`, `spectralanalysis.py`, `temporalchart.py`, `networkwalker.py`, `statssheet.py`, `modelinterface.py` (Streamlit), interactive HTML explorers.

### Key data facts (the real, full graph)
- **440,171 nodes / 5,767,668 directed edges** (5,494,276 unique undirected; 2,417 self-loops).
- Single giant connected component: **97.4 %** of nodes. Mean degree ~25, max 20,641. Mean clustering 0.28.
- Community structure: **Q ≈ 0.457, ~11,471 communities** (Louvain, full graph).
- Top hubs are **generic annotation/anatomy terms** (GO:0005515 "protein binding", GO roots, testis, liver), NOT disease IDs.
- Top edge types: `INTERACTS_WITH` (1.85M), `EXPRESSED_IN` (1.26M), `HAS_PHENOTYPE` (1.07M), `ORTHOLOGOUS_TO` (511K).
- Heavily **multi-organism**: mouse (MGI) and zebrafish (ZFIN) genes dominate over human (HGNC).

### Sources
Monarch, Orphadata, CiVIC, DGIdb, Gene2Phenotype, DIDA, ClinGen/ClinPGx, DisGeNET, DrugCentral, GWAS, LncRNADisease, LNCipedia. Note many are **not independent** (Monarch aggregates several others) — do not treat "supported by N datasets" as N independent votes.

---

## 2. Repository layout

```
data/raw/<Source><n>/data.raw.*      90 raw dumps across 13 databases (~7 GB)
data/processed/                       canonical_*, edges_*, relationships.csv, embeddings, models
scripts/                              the numbered pipeline + analysis scripts
analysis_outputs_20260304_*/          figure/CSV outputs from an earlier analysis run
full_graph_analysis/                  full-graph network results (degree/PR/eigenvector/communities/spectral)
fullgraph_*.png                       full-graph figures
*.html                                interactive graph explorers + Streamlit portal
*_report.txt, *_metrics.csv           assorted reports
PROJECT_REPORT.md                     overall state + publication roadmap (START HERE)
FULL_GRAPH_ANALYSIS_REPORT.md         results of re-running analyses on the full graph
venv/                                 a Windows virtualenv — do NOT use from the Linux sandbox
```

Canonical analysis input: `data/processed/relationships.csv` (`:START_ID,:END_ID,:TYPE,weight:int,dataset_sources`).

---

## 3. Environment, rules & constraints

**Compute sandbox (where code runs) is restrictive — plan around it:**
- **No internet / no package installs.** PyPI and apt are blocked. Only `numpy`, `pandas`, `matplotlib` are available. **`scipy`, `networkx`, `igraph`, `scikit-learn`, `pykeen` are NOT installed and cannot be installed here.** Graph algorithms had to be hand-written in numpy. Heavy scientific work should be run in the user's own full environment.
- **~45-second wall-clock limit per shell command.** Long jobs must be split and **checkpointed to disk** (filesystem persists between calls).
- **Background processes do NOT survive** between shell calls (fresh PID namespace each time). Don't rely on `nohup ... &` + poll. Run synchronously within the time budget, or checkpoint-and-resume across calls.
- **~3.8 GB RAM, 2 CPUs.** Prefer streaming, int32, CSR; avoid dense N×N anything (N = 440K).
- **Neo4j is part of the original workflow but is NOT available in the sandbox.** Some results (e.g. `drugrepurposing.py`) were produced in a local Neo4j instance and won't reproduce here. Do graph work directly on the CSVs in Python.
- The `venv/` in the repo is a **Windows** environment; it will not run under Linux. Ignore it.

**File-path rules:**
- File tools (Read/Write/Edit) use Windows paths: project root is `C:\Users\owenh\Downloads\ScienceFairYear2`.
- The shell sees the same folder at `/sessions/<id>/mnt/ScienceFairYear2/`. Translate accordingly; never write sandbox paths into deliverables.
- Save all final deliverables to the project folder so the user can see them; the scratch outputs dir is not user-visible.
- Scripts contain **hard-coded absolute Windows paths** and (in `drugrepurposing.py`) a **hard-coded Neo4j password** — must be removed before any public release.

**Correctness rules (learned the hard way):**
- **Build CSR adjacency by sorting**, not by vectorized scatter (`indices[fill[u]] = v` silently drops edges for repeated node IDs). Always verify `np.diff(indptr) == degree`.
- **Sanity-check every metric**: clustering ∈ [0,1]; λ_max of the Laplacian ≈ max degree; smallest Laplacian eigenvalue ≈ 0 for a connected component. A clustering value > 1 is what caught the CSR bug.
- Exact **betweenness** and the exact **Fiedler value** are not feasible in pure numpy at this scale — they need igraph/scipy. Current full-graph betweenness is a small sampled estimate; treat as indicative only.

---

## 4. Scientific constraints (publication bar)

- **Never evaluate a model on its training data.** `TransE.py` currently uses `training=tf, testing=tf` — this must be replaced with a proper held-out split (prefer a **time-based split** using the temporal data: train pre-2020, test 2020–2023).
- **Always report metrics with baselines.** Link prediction needs MRR, Hits@{1,3,10}, AUROC, AUPRC vs. random, common-neighbors/Adamic-Adar, and ≥1 modern KGE (DistMult/ComplEx/RotatE).
- **Filter invalid predictions** (self-loops, already-existing edges, type-incompatible pairs) before ranking.
- **Validate biologically**: check top predictions against a newer DB release or literature; report precision@k.
- **Account for source non-independence** when weighting consensus.
- **Down-weight or remove ultra-generic hubs** (GO:0005515, GO roots) for biologically meaningful prediction/centrality.
- **Reproducibility**: fixed seeds, recorded source versions/download dates, relative paths, no secrets, environment lockfile.
- Don't fabricate results. If a computation can't be completed in this environment, say so and specify what's needed (e.g. "run with scipy").

---

## 5. Goals & roadmap

**North star:** a publishable resource/method paper (realistic framing: an Application Note or resource paper for *Bioinformatics*, *Database (Oxford)*, *BMC Bioinformatics*, *GigaScience*, or *PLOS Computational Biology*).

Priority order (from `PROJECT_REPORT.md`):
- **Tier 1 (must-fix):** (1) run all analyses on the full graph — **DONE**, see `FULL_GRAPH_ANALYSIS_REPORT.md`; (2) real train/test (time-based) split; (3) metrics + baselines; (4) filter invalid predictions; (5) remove secrets/absolute paths.
- **Tier 2 (credibility):** biological validation of top predictions; justify/calibrate the consensus score and handle source overlap; tune+document embeddings; null-model tests for network metrics; record data provenance.
- **Tier 3 (novelty/polish):** crisp contribution statement vs. Hetionet/PrimeKG/Monarch KG/OpenBioLink; one focused validated use case (the rare-disease/Orphadata coverage is a good candidate); reproducible public release.

---

## 6. Working conventions for assistants

- Before substantive work, skim `PROJECT_REPORT.md` and `FULL_GRAPH_ANALYSIS_REPORT.md` so you don't re-derive context.
- When you produce results, write tables/figures to the project folder and a short markdown report; keep claims tied to actual computed numbers.
- Prefer editing/extending existing scripts over creating parallel ones; keep the numbered-phase structure.
- Flag, don't silently fix, anything that changes a scientific conclusion — explain what changed and why.
- Keep responses and docs concise and direct (the user prefers minimal verbosity).
