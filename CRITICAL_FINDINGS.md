# Critical findings (2026-06-18)

This session ran the link-prediction work on the **real full graph** (PyKEEN + GPU,
not the sandbox) and audited data provenance. One finding is publication-blocking
and needs a scope decision before anything else matters.

---

## ✅ UPDATE 2026-06-21 — Finding 1 partially resolved (re-ingestion)

Resolution path #1 (re-ingest the other databases) has been substantially carried
out. `scripts/extract_additional_sources.py` now ingests **7 non-Monarch sources**
and the integrated graph (`data/processed/edges_clean_integrated.csv`) is
**5,860,539 edges / 452,012 nodes**:

| Source | Edges | Source | Edges |
|---|---|---|---|
| DGIdb | 42,389 | Gene2Phenotype | 3,319 |
| GWAS | 31,194 | CiVIC | 929 |
| DrugCentral | 13,944 | DIDA | 140 |
| Orphadata | 8,329 | **Monarch** | 5,773,610 |

A gene crosswalk (HGNC complete set + `entity_id_map`) maps symbol/Ensembl/Entrez/
UniProt → HGNC, and DOID/EFO/Orphanet→MONDO via `mondo_xrefs.tsv`, so the new
sources merge into the existing namespaces. **Caveat:** Monarch is still ~98.5% of
edges, so the honest framing is *"a Monarch-derived KG augmented with 7 curated
sources,"* not *"equal integration of 13 databases."*

**Still zero clean edges (genuinely not integrable from these dumps; see
`data/processed/provenance/source_provenance.md`):** DisGeNET (variant table, no
disease column), LNCipedia (FASTA only), ClinGen/PharmGKB (node/xref tables only,
no relationship table), LncRNADisease (ncRNA + disease both free text).

All downstream was regenerated on the rebuilt graph: provenance, KGE benchmark
(`benchmark_integrated/`), full-graph analysis (`full_graph_analysis_integrated/`),
null-model tests, Adamic-Adar predictions, consensus confidence.

---

## 🔴 FINDING 1 — The graph is 100% Monarch, not "13 integrated databases"

**Evidence (all reproducible):**

- `data/processed/edges_clean.csv` (5,767,668 edges) — the `dataset_sources`
  column contains **only** `Monarch*` tokens. Zero non-Monarch tokens:
  `cut -d, -f5 edges_clean.csv | grep -iv monarch` → empty.
- Same for `edges_canonicalized.csv`.
- The raw extraction `extracted_edges_20260205_213914.csv` (11,477,215 edges) has
  `dataset_id ∈ {Monarch, Monarch2 … Monarch19}` and **nothing else**.
- `data/registry/dataset_registry.csv` *did* scan all 13 databases (CiVIC, DGIdb,
  DIDA, ClinGen, Orphadata, …), so the folders are present and readable — but **no
  edges were ever extracted from them**. Only the 19 Monarch dumps became edges.

**Why the graph still looks diverse:** Monarch is itself a multi-species
*aggregator*. The MGI / ZFIN / RGD / HGNC / CHEBI / CLINVAR IDs in the graph all
come from *inside* the Monarch dumps — not from the standalone Orphadata, GWAS,
CiVIC, DGIdb, etc. folders, which contributed nothing.

**Why it's blocking:** the project's headline claim — *"we integrate ~13 public
biomedical databases"* (CLAUDE.md, PROJECT_REPORT.md §3, target narrative §8) —
is **not supported by the data as it currently stands**. A reviewer who runs the
one-line check above will reject the central framing. This is more important than
the train/test-split issue, because it concerns what the resource actually *is*.

**Two honest resolutions (mutually exclusive — this is the decision needed):**

1. **Re-ingest the other 12 databases.** Run the ingestion pipeline so Orphadata,
   GWAS, CiVIC, DGIdb, Gene2Phenotype, DIDA, ClinGen/ClinPGx, DisGeNET,
   DrugCentral, LncRNADisease, LNCipedia actually contribute edges, then
   re-canonicalize/clean/assemble. This makes the 13-source claim true but is a
   substantial pipeline effort (XML/JSON/FASTA parsing + per-schema relationship
   mapping for 12 heterogeneous sources) and **re-runs everything downstream**
   (the benchmark below would have to be recomputed on the new graph).

2. **Reframe honestly as a Monarch-derived resource.** Keep the current graph and
   describe it as *"a 440K-node / 5.8M-edge multi-species knowledge graph derived
   from the Monarch Initiative, benchmarked for gene–disease link prediction."*
   Still publishable (as a benchmark/analysis paper), honest, and immediately
   defensible — but it drops the "novel integration of 13 sources" novelty.

Everything else in this report (the KGE benchmark, baselines, null models,
consensus audit) is valid **for the current Monarch-derived graph** under either
choice; only the framing and the provenance table change.

---

## 🟠 FINDING 2 — "Consensus" / TransE artifacts are degenerate

- The old `transe_new_predictions.csv` contains **self-loops**
  (`MONDO:0006877 → MONDO:0006877`) and was produced by a model trained on the
  tiny `transetableuseforpykeen.csv` (the 606-node consensus subset), evaluated on
  its own training data. It is superseded by the filtered predictions from the new
  benchmark (`data/processed/benchmark/novel_predictions_filtered.csv`).

---

## ✅ What this session fixed / produced (valid regardless of the Finding-1 choice)

- **Removed all hard-coded secrets and absolute paths** → `scripts/config.py`
  (Neo4j creds from env; repo-relative paths). Verified zero remaining.
- **Real KGE benchmark** (`scripts/kge_benchmark.py`) now actually *runs* on this
  machine (RTX 5070, CUDA): TransE/DistMult/ComplEx/RotatE vs Random/CN/Adamic-Adar
  under one protocol — held-out split, filtered all-entity ranking, and **shared**
  AUROC/AUPRC with both random **and type-matched hard** negatives.
- **Filtered prediction generation** (no self-loops / known edges / type mismatches).
- **Provenance table** (`scripts/provenance_table.py`) — which surfaced Finding 1.
- See `LINK_PREDICTION_BENCHMARK.md` for the metrics once the run completes.
