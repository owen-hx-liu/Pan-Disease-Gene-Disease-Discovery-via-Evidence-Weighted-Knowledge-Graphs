# Quantifying evaluation leakage in multi-species gene-disease link prediction

A graded leakage audit and a de-leaked benchmark for gene-disease link prediction on a
multi-species, Monarch-derived knowledge graph. This repository holds the code, the
de-leaked evaluation splits, the benchmark results, and the scripts that regenerate every
table and figure from the frozen result files.

> **Scope note.** This is a leakage-aware benchmark, not a database-integration resource.
> The knowledge graph is about 98.5% Monarch Initiative content (Monarch is itself an
> aggregator) augmented with seven curated non-Monarch sources. It is described honestly as
> a Monarch-derived graph, not an equal integration of many independent databases.

## What this is

Knowledge-graph embedding (KGE) models are widely used to predict gene-disease associations
and are usually reported with high ranking and AUROC scores. Those scores can be inflated by
evaluation leakage. This project measures how much, and from which sources, using a graded
set of controlled evaluation regimes run through one shared ranking harness.

Two KGE models (TransE, RotatE) and five topological baselines (Random, Common-Neighbors,
Adamic-Adar, Jaccard, Preferential-Attachment) are evaluated on 4,228 held-out human
gene-disease test edges across four regimes, with three seeds (42, 1, 7):

| Regime | Name | What it removes |
|---|---|---|
| R0 | standard | nothing (ordinary random split) |
| R1 | redundancy | direct held-pair edges (train copies of test facts) |
| R2 | degree-preserving null | genuine structure, while preserving the degree sequence exactly |
| R3 | orthology-blocked | cross-species ORTHOLOGOUS_TO bridges |

The graph has 452,012 nodes and 5,860,539 directed edges across 28 relation types and 8
node categories.

## Key result

Evaluation leakage is method-class-dependent:

- **Node-degree structure (R2) is the one universal leak.** Topological baselines lose about
  43 to 49% MRR and KGE about 21 to 24%. Preferential-Attachment (pure degree) is unchanged
  by construction, which confirms the null preserves the degree sequence exactly.
- **Train/test redundancy (R1) leaks only into multi-hop KGE** (about 6 to 9%) and is
  invisible to the overlap heuristics, because a direct gene-disease edge is embedded as
  proximity by KGE but is not a shared neighbor.
- **Cross-species orthology (R3) is a designed negative control that comes back null**
  (within about 0.5% for every method). This is the domain-specific hypothesis a biomedical
  reviewer would raise, and the data refutes it rather than confirming it.

A filtered full-ranking robustness check, in which each true disease is ranked against the
entire pool of about 14,300 diseases instead of 50 sampled negatives, reproduces and
strengthens the degree effect. Exact per-method numbers are in `tables/`: Table 2 for the
sampled-negative protocol and Table S1 for full ranking.

## Repository layout

```
scripts/                     pipeline and analysis code
  build_deleaked_splits.py   build the R0/R1/R3 splits and the split manifest
  build_degree_null.py       build the R2 degree-preserving null
  run_baselines.py           topological baselines across regimes
  run_kge.py                 TransE + RotatE (GPU); --full-rank for the robustness run
  run_robustness.py          baseline full-ranking robustness
  hetionet_audit.py          cross-graph replication on Hetionet
  lib_eval.py                the single shared ranking and scoring harness
  make_tables.py             Tables 1-4 and S1 from result JSONs (no hand-typed numbers)
  make_figures.py            Figures 1-3 from result JSONs
data/processed/
  splits/                    the de-leaked benchmark (train/valid/test, R1/R2/R3, manifest)
  results/                   per-regime baseline, KGE, robustness, and Hetionet results
  provenance/                per-source contribution table
tables/                      generated tables (Markdown and LaTeX)
figures/                     generated figures (PNG and PDF)
run_all.sh                   one-command end-to-end reproduction
requirements.txt             Python dependencies
```

## Reproducing the results

```bash
python -m venv venv
source venv/bin/activate            # Windows: venv\Scripts\activate
pip install -r requirements.txt
# GPU (RTX 50-series / Blackwell needs CUDA 12.8+):
pip install torch==2.11.0 --index-url https://download.pytorch.org/whl/cu128

./run_all.sh                 # full pipeline (the KGE step needs a CUDA GPU)
./run_all.sh --skip-kge      # everything except KGE, reusing the committed kge_summary.json
```

The tables and figures read only the frozen result JSONs, so you can regenerate them without
a GPU:

```bash
python scripts/make_tables.py    # tables/table1-4 and tableS1 (.md and .tex)
python scripts/make_figures.py   # figures/fig1-3 (.png and .pdf)
```

All runs use fixed seeds (42, 1, 7). Paths are repo-relative; there are no hard-coded machine
paths or secrets. Heavy steps checkpoint to disk. Reproducing the KGE arm from scratch needs
a CUDA GPU and is slow (about 40 minutes per model-by-regime-by-seed cell); the committed
result JSONs let you rebuild every table and figure without retraining.

## Data, provenance, and licensing

The knowledge graph is derived from the Monarch Initiative (about 98.5% of edges; CC BY 4.0)
augmented with seven curated sources: DGIdb, GWAS Catalog, DrugCentral, Orphadata,
Gene2Phenotype, CIViC, and DIDA. Per-source contributions are in
`data/processed/provenance/`. Because Monarch aggregates several primary sources, do not read
"supported by N datasets" as N independent votes.

Two licenses apply:

- **Code** (everything in `scripts/`): MIT. See [`LICENSE`](LICENSE).
- **Data** (the de-leaked splits and result files under `data/processed/`): CC BY 4.0, with
  attribution to the Monarch Initiative and the underlying primary sources. See
  [`LICENSE-DATA`](LICENSE-DATA).

## Citation

If you use this benchmark or code, please cite the archived release. Machine-readable metadata
is in [`CITATION.cff`](CITATION.cff) and [`.zenodo.json`](.zenodo.json); the Zenodo DOI will
be added on release.

## Caveats

- **No real temporal split.** The graph carries no genuine edge-discovery dates, so the
  benchmark controls optimism with the leakage regimes R0 to R3 on a random held-out split,
  not a time split.
- **Coverage is Monarch-dominated.** Even after integrating seven more sources, about 98.5%
  of edges are Monarch content.
- **Name the protocol.** Ranking metrics are reported for both a sampled-negative protocol
  (Table 2) and a filtered full-ranking protocol (Table S1); the two can rank methods
  differently, so always state which one a number comes from.
