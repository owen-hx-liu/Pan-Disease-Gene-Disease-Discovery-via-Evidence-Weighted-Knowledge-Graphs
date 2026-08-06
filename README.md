# Quantifying evaluation leakage in multi-species gene-disease link prediction

A graded leakage audit and a de-leaked benchmark for gene-disease link prediction on a
multi-species, Monarch-derived knowledge graph. This repository holds the code, the
de-leaked evaluation splits, the benchmark results, and the scripts that regenerate every
table and figure from the frozen result files.

> **Scope note.** This is a leakage-aware benchmark, not a database-integration resource.
> The knowledge graph is about 98.3% Monarch Initiative content (Monarch is itself an
> aggregator) augmented with seven curated non-Monarch sources. It is described honestly as
> a Monarch-derived graph, not an equal integration of many independent databases.

## What this is

Knowledge-graph embedding (KGE) models are widely used to predict gene-disease associations
and are usually reported with high ranking and AUROC scores. Those scores can be inflated by
evaluation leakage. This project measures how much, and from which sources, using a graded
set of controlled evaluation regimes run through one shared ranking harness.

Four KGE models spanning two families (translational: TransE, RotatE; bilinear: DistMult,
ComplEx) and five topological baselines (Random, Common-Neighbors, Adamic-Adar, Jaccard,
Preferential-Attachment) are evaluated on 4,228 held-out human gene-disease test edges
across four regimes, with three seeds (42, 1, 7):

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

- **Node-degree structure (R2) is the leak shared by every method tested.** Topological
  baselines lose about 43 to 49% MRR and the four KGE models about 21 to 37%.
  Preferential-Attachment (pure degree) is unchanged by construction, which confirms the null
  preserves the degree sequence exactly. Scope: this covers the five parameter-free
  topological heuristics and the four shallow embedding models; message-passing architectures
  could not be tested in the regime where the leak exists (see the boundary note below).
- **Train/test redundancy (R1) leaks only into the embedding models** (about 6 to 13%) and is
  invisible to the overlap heuristics, because a direct gene-disease edge is embedded as
  proximity by KGE but is not a shared neighbor.
- **Cross-species orthology (R3) is a designed negative control that comes back null** — it
  lowers no method's score, leaving the topological baselines within about 0.5% and raising
  the embedding models slightly (+0.2 to +3.0% under sampled negatives, +2 to +16% under full
  ranking). This is the domain-specific hypothesis a biomedical reviewer would raise, and the
  data refutes it rather than confirming it.

A filtered full-ranking robustness check, in which each true disease is ranked against the
entire pool of about 14,300 diseases instead of 50 sampled negatives, reproduces and
strengthens the degree effect: the KGE collapse deepens to 70 to 82% and the overlap
heuristics to 88 to 91%, while Preferential-Attachment still does not move. The check is
complete, covering all four KGE models and all five baselines across all four regimes and
three seeds. Exact per-method numbers are in `tables/`: Table 2 for the sampled-negative
protocol and Table S1 for full ranking.

Two controls bound the degree result, and both are reported. Table S2 shows the R2 collapse
is not an artifact of hyperparameters chosen on the standard split: the models fit the
rewired training graph almost as well as the real one while their held-out scores collapse,
and a matched learning-rate by negatives grid run on both regimes peaks at the same
configuration in each. Tables S3 and S4 report a message-passing (R-GCN) probe, which is a
boundary rather than a result: on this hardware the model trains only on subgraphs from which
the degree leak has already dissolved (the matched control retains 12% of its full-graph
collapse at the fraction we can train on, and the validity threshold is first cleared at a
fraction costing roughly a month of training), so the paper's degree claim is scoped to the
shallow-embedding and topological-heuristic families and message passing is stated as
untested.

A degree-stratified re-analysis (Table 5, Figure 4) decomposes each method's MRR by gene- and
disease-degree quartile, recomputed from the stored per-edge ranks with no retraining. It
localizes the leak: performance climbs steeply with disease degree (Preferential-Attachment
gains +0.90 MRR from the rarest to the most popular disease quartile), and the degree-null
penalty falls almost entirely on the rare-disease, low-degree strata. A worked case study
(Table 6) shows the practical cost: hub diseases such as insomnia are ranked first by degree
alone, while genuine rare-disease associations (for example RS1 to X-linked retinoschisis)
are pushed to the bottom of the list.

## Repository layout

```
scripts/                     pipeline and analysis code
  build_deleaked_splits.py   build the R0/R1/R3 splits and the split manifest
  build_degree_null.py       build the R2 degree-preserving null
  run_baselines.py           topological baselines across regimes
  run_kge.py                 TransE, RotatE, DistMult, ComplEx (GPU); --full-rank for the robustness run
  run_robustness.py          baseline full-ranking robustness
  hetionet_audit.py          cross-graph replication on Hetionet
  degree_stratified.py       degree-stratified re-analysis (gene vs disease degree; no retraining)
  case_study.py              worked case study (degree false-positives vs missed rare edges)
  lib_eval.py                the single shared ranking and scoring harness
  make_tables.py             Tables 1-6 and S1 from result JSONs (no hand-typed numbers)
  make_figures.py            Figures 1-4 from result JSONs
data/processed/
  splits/                    the de-leaked benchmark (train/valid/test, R1/R2/R3, manifest)
  results/                   per-regime baseline, KGE, robustness, Hetionet, degree-stratified, and case-study results
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
./run_all.sh --skip-kge      # everything except KGE, reusing the frozen result JSONs
```

A bare `git clone` only carries a partial subset of `data/processed/results/` (kept as a small,
browsable sample; GitHub has no clean way to track ~9 MB of JSON without also inviting the 160 MB
`splits.zip` past its 100 MB per-file limit). To regenerate every table and figure without a GPU,
first download and unzip `splits.zip` and `results.zip` from the
[Zenodo record](https://doi.org/10.5281/zenodo.21828739) into `data/processed/`, then:

```bash
python scripts/make_tables.py    # tables/table1-6 and tableS1 (.md and .tex)
python scripts/make_figures.py   # figures/fig1-4 (.png and .pdf)
```

All runs use fixed seeds (42, 1, 7). Paths are repo-relative; there are no hard-coded machine
paths or secrets. Heavy steps checkpoint to disk. Reproducing the KGE arm from scratch needs
a CUDA GPU and is slow (about 40 minutes per model-by-regime-by-seed cell); the frozen result
JSONs on Zenodo let you rebuild every table and figure without retraining.

## Data, provenance, and licensing

The knowledge graph is derived from the Monarch Initiative (about 98.3% of edges; CC BY 4.0)
augmented with seven curated sources: DGIdb, GWAS Catalog, DrugCentral, Orphadata,
Gene2Phenotype, CIViC, and DIDA. Per-source contributions are archived as
`provenance/source_provenance.csv` in `provenance.zip` on the
[Zenodo record](https://doi.org/10.5281/zenodo.21828739) (the `data/processed/provenance/`
path is git-ignored, not part of the GitHub tree). Because Monarch aggregates several primary
sources, do not read "supported by N datasets" as N independent votes.

Two licenses apply:

- **Code** (everything in `scripts/`): MIT. See [`LICENSE`](LICENSE).
- **Data** (the de-leaked splits and result files under `data/processed/`): CC BY 4.0, with
  attribution to the Monarch Initiative and the underlying primary sources. See
  [`LICENSE-DATA`](LICENSE-DATA).

## Citation

If you use this benchmark or code, please cite the archived release:
[10.5281/zenodo.21828739](https://doi.org/10.5281/zenodo.21828739). Machine-readable metadata
is in [`CITATION.cff`](CITATION.cff) and [`.zenodo.json`](.zenodo.json).

## Caveats

- **No real temporal split.** The graph carries no genuine edge-discovery dates, so the
  benchmark controls optimism with the leakage regimes R0 to R3 on a random held-out split,
  not a time split.
- **Coverage is Monarch-dominated.** Even after integrating seven more sources, about 98.3%
  of edges are Monarch content.
- **Name the protocol.** Ranking metrics are reported for both a sampled-negative protocol
  (Table 2) and a filtered full-ranking protocol (Table S1); the two can rank methods
  differently, so always state which one a number comes from.
