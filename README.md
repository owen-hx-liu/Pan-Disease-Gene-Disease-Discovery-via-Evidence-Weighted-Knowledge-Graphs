# Biomedical Knowledge Graph & Link Prediction

An integrated biomedical knowledge graph with honest, baseline-anchored gene–disease
link-prediction benchmarking and network-science analysis. Built for a peer-reviewed
resource/benchmark paper.

> **Read first:** [`CRITICAL_FINDINGS.md`](CRITICAL_FINDINGS.md) documents a
> publication-blocking data issue (the original graph was 100% Monarch-derived) and
> what has been done about it. [`PROJECT_REPORT.md`](PROJECT_REPORT.md) is the full
> assessment and publication roadmap.

## What's here

- **Integrated graph** — `data/processed/edges_clean_integrated.csv`
  (448,931 nodes / 5,814,437 edges). Monarch + three newly-ingested standalone
  databases (DGIdb, Orphadata, Gene2Phenotype); the original `edges_clean.csv` was
  Monarch-only.
- **Link-prediction benchmark** — TransE / DistMult / ComplEx / RotatE vs. Random /
  Common-Neighbors / Adamic-Adar under one protocol (held-out split, filtered
  all-entity ranking, shared AUROC/AUPRC with random **and type-matched hard**
  negatives, filtered novel-prediction generation).
- **Statistical network analysis** — structural metrics vs. a degree-preserving
  null model (z-scores).
- **Independence-aware consensus** — noisy-OR edge confidence that collapses
  aggregator sources and is sensitivity-tested.
- **Provenance** — per-source contribution table.

## Install

```bash
python -m venv venv && venv/Scripts/activate         # Windows; use source venv/bin/activate on *nix
pip install -r requirements.txt
# GPU (RTX 50-series / Blackwell needs CUDA 12.8+):
pip install torch==2.11.0 --index-url https://download.pytorch.org/whl/cu128
```

Paths are repo-relative (`scripts/config.py`); there are no hard-coded machine
paths or secrets. Neo4j scripts read credentials from `NEO4J_PASSWORD` (env var).

## Pipeline

```bash
# 1. Ingest the non-Monarch databases the original pipeline skipped, then merge
python scripts/extract_additional_sources.py        # -> data/processed/additional_edges.csv
python scripts/merge_integrated_edges.py            # -> data/processed/edges_clean_integrated.csv

# 2. Provenance / consensus / null-model analyses
python scripts/provenance_table.py    --edges data/processed/edges_clean_integrated.csv
python scripts/consensus_confidence.py
python scripts/null_model_tests.py    --edges data/processed/edges_clean_integrated.csv --realizations 5

# 3. Link-prediction benchmark (GPU strongly recommended)
python scripts/kge_benchmark.py --edges data/processed/edges_clean_integrated.csv \
    --out data/processed/benchmark_integrated --epochs 100 --dim 128
#   --skip-kge runs baselines only; outputs benchmark_table.md + novel_predictions_filtered.csv
```

All scripts use a fixed seed (42). Heavy steps checkpoint to disk.

## Reports

| File | Contents |
|---|---|
| `CRITICAL_FINDINGS.md` | Monarch-only data issue + what was fixed |
| `LINK_PREDICTION_BENCHMARK.md` | Benchmark protocol & results |
| `NETWORK_NULL_MODEL_REPORT.md` | Structure vs. null model (z-scores) |
| `CONSENSUS_CALIBRATION_REPORT.md` | Consensus score rebuild + sensitivity |
| `data/processed/provenance/source_provenance.md` | Per-source contribution table |

## Key caveats (state these in the paper)

- **No real temporal split.** The graph carries no genuine edge-discovery dates, so a
  time-based split is impossible; the benchmark controls optimism with the leakage
  regimes R0–R3 on a random held-out split, not a time split.
- **Coverage is still Monarch-dominated** (~99% of edges) even after integrating
  three more sources; the remaining standalone DBs (CiVIC, DisGeNET, DrugCentral,
  GWAS, ClinGen, …) need ingestion (they require cross-file joins or free-text →
  ontology mapping). See `extract_additional_sources.py` docstring.
- **Source non-independence**: Monarch aggregates several primary sources; the
  consensus score collapses aggregators rather than counting them as independent.
