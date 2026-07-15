# BUILD_GUIDE.md — Build the whole project, script by script, in order

*This is the construction manual. Each "build unit" below is one script (or one edit to an
existing script), in dependency order. For each: what it is, **who can run it**, inputs →
outputs, how it works, the acceptance test that proves it's correct, and a copy-paste prompt
you can hand to an AI to generate it. Do them top to bottom.*

**Status legend:** `NEW` = create it · `EXTEND` = modify an existing script · `REUSE` = already
works, just run it · `DONE` = already built in this session.

## Who runs what — the owner legend (read this)

The roles are about **the resource a task needs**, not fixed hardware. Every unit is tagged:

| Tag | Meaning | Which machine |
|---|---|---|
| 🟢 **A or C** | CPU-light: data, splits, figures, writing, coordination | **either** — run on A *during GPU downtime*, or on C. Pick whoever's free. |
| 🔴 **A-only (GPU)** | needs the RTX GPU (KGE training) | **A only** — no other machine can do it. |
| 🟠 **B (heavy CPU)** | RAM/CPU-heavy over 5.86M edges (baselines, leakage eval, LOSO) | **B** preferred. Can run on A *only when the GPU is idle* — never alongside KGE training (they fight for CPU/RAM). |
| 🟣 **split A+B** | parts go to different machines (see the unit) | the GPU part → A, the topology part → B |

**Why "A or C" exists:** A's GPU trains unattended, so while the KGE sweep runs, A's CPU and
you are free — that's exactly when to knock out the 🟢 light tasks. If you're short on people
or machines, **A can wear the C hat entirely**; you only truly need a separate C if a second
person will write/plot in parallel while A trains.

**Two hard rules (independent of which machine):**
1. **One owner per script.** Whoever builds/edits a script owns it — don't have two machines
   editing the same file (merge hell). Ownership can move between machines, but only one at a time.
2. **Only ONE machine ever regenerates the splits or the graph.** If A builds the splits, A
   owns split regeneration — C must not also rebuild them, or the hashes drift and results stop
   being comparable. (Everyone else *consumes* the hashed files.)

**Validated backbone numbers (what "correct" looks like):** 42,288 target edges → train
5,852,083 / valid 4,228 / test 4,228; R2 removes 520,198 bridge edges; 244 hub nodes.

---

## Dependency graph (build order)

```
STAGE 1  backbone                     STAGE 2  methods                STAGE 3  method + validation
  1. graph_stats.py ........ 🟢A/C     4. lib_eval.py ....... 🟢A/C    7. build_hybrid.py ..... 🟣A+B
  2. build_deleaked_splits [DONE] 🟢A/C 5. run_baselines.py .. 🟠B      8. loso_validation.py .. 🟠B(+🔴A)
  3. provenance fill ....... 🟢A/C     6. run_kge.py ........ 🔴A      9. newer_release_diff.py 🟢A/C
                                                                     10. literature_curate.py 🟢A/C
STAGE 4  outputs                      STAGE 5  release
 11. make_figures.py ....... 🟢A/C     13. run_all.ps1 + Zenodo 🟢A/C
 12. make_tables.py ........ 🟢A/C     14. manuscript ......... 🟢A/C
```

Nothing in Stage 2 can start until Unit 4 (`lib_eval.py`) and the split files (Unit 2) exist
and are hashed — that's why those two come first and have a single owner.

**The critical path is the 🔴/🟠 column (Units 6 → 7 → 8 on A and B).** Everything 🟢 can be
slotted into A's GPU-downtime or handed to C, so it never gates the timeline.

---

# STAGE 1 — Backbone (do first) 🟢 A or C

> Since **A already has the data**, the fastest path is for A to do this whole stage first
> (it unblocks A's own KGE runs and B's baselines), then start training. If a separate C exists,
> C can do it instead — just don't do it on both.

### Unit 1 — `scripts/graph_stats.py`  `NEW` · 🟢 A or C
- **Purpose:** produce Table 1 of the paper (graph size + composition) from the frozen graph,
  as JSON so no number is ever hand-typed.
- **In:** `data/processed/edges_clean_integrated.csv`
- **Out:** `data/processed/graph_stats.json`
- **Does:** counts nodes, directed edges, unique undirected edges, self-loops; per-namespace
  node counts (HGNC/MGI/ZFIN/MONDO/GO/…); per-relation edge counts; degree summary
  (mean/median/max), giant-component fraction (optional, needs a BFS).
- **Acceptance test:** prints `edges: 5860539`, `relations: 28`, and the top relation is
  `BIOLINK:INTERACTS_WITH` (~1.9M). Matches the numbers in `PAPER_SCOPE.md`.
- **Prompt:**
  ```
  Write scripts/graph_stats.py (pandas only). Input
  data/processed/edges_clean_integrated.csv with columns
  source_id,relation,target_id,weight,dataset_sources. Compute and save to
  data/processed/graph_stats.json: n_nodes (unique over source+target), n_edges_directed,
  n_edges_unique_undirected, n_self_loops, per-namespace node counts (split ID on ':' take
  prefix), per-relation edge counts, degree mean/median/max. Use scripts/kg_categories.py
  category_of to also report per-category node counts. Print a human-readable summary. Fixed,
  deterministic, no network libs.
  ```

### Unit 2 — `scripts/build_deleaked_splits.py`  `DONE` · 🟢 A or C (ONE owner — regenerates splits)
Already built and validated this session. **Run it for real** (drop `--dry-run`) to write the
split files, then:
```powershell
python scripts/build_deleaked_splits.py
Get-FileHash -Algorithm SHA256 data\processed\splits\*.csv
```
Add every printed sha256 to `ARTIFACT_HASHES.txt`, and **upload the whole `data/processed/splits/`
folder to the shared data channel** so A and B pull byte-identical files.
- ⚠️ **This is the one unit that must have a single owner forever** (rule 2). Whichever machine
  runs it owns all future re-runs. A is the natural choice (it has the data).
- **Acceptance test:** `split_manifest.json` shows train 5,852,083 / valid 4,228 / test 4,228,
  removed_for_R2 520,198, hub_nodes 244. No held-out edge appears in train (the script asserts this).
- **Later:** once Unit 9 gives you a MONDO label file, re-run with
  `--gene-symbols data/registry/hgnc_complete_set.txt --disease-labels data/registry/mondo_labels.tsv`
  to activate the R3 tautology filter, and re-hash.

### Unit 3 — Provenance finalize  `EXTEND`/data · 🟢 A or C (manual)
- **Purpose:** fill the two things reviewers require and you're missing.
- **Does:** run the existing `scripts/provenance_table.py`, then **manually** edit the output
  to add each source's exact upstream **release version + download date** and **license**
  (the current table flags `source_release_version` as a placeholder).
- **Acceptance test:** `data/processed/provenance/source_provenance.md` has no "placeholder"
  strings; every ingested source has a version, a date, and a license.

---

# STAGE 2 — Shared eval library + method runners

### Unit 4 — `scripts/lib_eval.py`  `NEW` · 🟢 A or C (shared library — the linchpin, ONE owner)
- **Purpose:** the single implementation of the ranking/metric protocol that **both**
  `run_baselines.py` (B) and `run_kge.py` (A) import, so every method is scored identically.
  One owner so it never diverges — build it wherever, but only one machine edits it.
- **In:** the split files (Unit 2).
- **Out:** importable functions (no CLI needed).
- **Must provide:**
  - `load_regime(regime, splits_dir)` → returns `(train_edges, test_edges, hub_set)` per the
    manifest map: R0/R1 use `train.csv`; R2/R3 use `train_R2_orthology_blocked.csv`; R3 uses
    `test_R3_deleaked.csv`; R1/R3 set `hub_filter=True`.
  - `rank_test_edges(score_fn, train_edges, test_edges, hub_set, hub_filter, n_neg=50, seed=42)`
    → for each test edge `(g, d_true)`: build a candidate set = `d_true` + `n_neg`
    **type-matched** disease negatives (same category, not an existing true disease of `g` →
    *filtered* setting); score all with `score_fn`; record the rank of `d_true`. This
    reproduces the existing `sampled-50neg` protocol in `benchmark_results.json`.
  - `ranking_metrics(ranks)` → MRR, Hits@1, Hits@3, Hits@10.
  - `classification_metrics(pos_scores, neg_scores)` → AUROC, AUPRC (support random AND
    type-matched negatives, like the current benchmark).
  - `bootstrap_ci(values, n=1000, seed=42)` → mean + 95% CI.
  - `paired_bootstrap_pvalue(rr_method_a, rr_method_b)` → significance of the MRR difference
    on per-edge reciprocal ranks.
- **Hub filter semantics (state in the code + paper):** for topological scorers, `hub_filter`
  removes hub nodes from the shared-neighbor evidence. For KGE the score is over embeddings, so
  R1 ≡ R0 for KGE (document this; R1 exists to test whether the *baselines'* win is hub-driven).
- **Acceptance test:** feed it a random scorer → MRR ≈ 1/((n_neg+1)/… ) ≈ chance (~0.09 for
  n_neg=50, matching the Random row already in `benchmark_results.json`); AUROC ≈ 0.50. If the
  Random control isn't at chance, the harness has a bug.
- **Prompt:**
  ```
  Write scripts/lib_eval.py: a shared evaluation library for a leakage-aware gene-disease link
  prediction benchmark. It reads fixed splits produced by build_deleaked_splits.py
  (data/processed/splits/: train.csv, valid.csv, test.csv, train_R2_orthology_blocked.csv,
  test_R3_deleaked.csv, hub_nodes.txt, split_manifest.json). Provide functions: load_regime(name)
  returning (train_edges, test_edges, hub_set) for R0/R1/R2/R3 per the manifest's regime map;
  rank_test_edges(score_fn, train_edges, test_edges, hub_set, hub_filter, n_neg=50, seed=42) that
  for each test (gene,true_disease) builds true+50 type-matched disease negatives (exclude the
  gene's other known diseases -> filtered), scores with score_fn, and returns the rank of the
  true disease; ranking_metrics(ranks)->MRR,Hits@1/3/10; classification_metrics(pos,neg)->
  AUROC,AUPRC; bootstrap_ci; paired_bootstrap_pvalue on per-edge reciprocal ranks. Use
  scripts/kg_categories.category_of for type matching. Deterministic (seed 42). Include a
  __main__ self-test that runs a RANDOM scorer and asserts MRR is near chance and AUROC near 0.5.
  ```

### Unit 5 — `scripts/run_baselines.py`  `NEW` (may fold in `bench_baselines.py`) · 🟠 B (heavy CPU)
- **Runs on B** (builds adjacency + ranks over 5.86M edges — RAM-heavy). Can run on A only if the
  GPU is idle; never during a KGE sweep.
- **Purpose:** Random, Common-Neighbors, Adamic-Adar over all four regimes R0–R3.
- **In:** split files + `lib_eval`.  **Out:** `data/processed/results/baselines_<regime>.json`
- **Does:** builds an undirected adjacency from the regime's train graph (networkx or dict of
  sets); defines `score_fn` for CN and AA (reuse the hub-capped AA logic from
  `predict_adamic_adar.py`); runs `lib_eval.rank_test_edges` for each regime × 3 seeds; saves
  metrics with CIs.
- **Acceptance test:** under R0, AdamicAdar MRR ≈ 0.69, Hits@10 ≈ 0.83 (matches
  `benchmark_results.json`). Under R2 (orthology-blocked) the MRR **drops** — that drop is the
  leakage signal; record it.
- **Prompt:**
  ```
  Write scripts/run_baselines.py. Import scripts/lib_eval.py. For each regime in [R0,R1,R2,R3]
  and each seed in [42,1,7]: load the regime, build an undirected graph from train_edges, define
  Random, CommonNeighbors, and Adamic-Adar scorers (skip shared neighbors with degree>2000 as in
  scripts/predict_adamic_adar.py; when hub_filter is on also exclude hub_set nodes from the
  neighbor evidence), evaluate with lib_eval.rank_test_edges + ranking_metrics +
  classification_metrics, and collect per-edge reciprocal ranks. Write
  data/processed/results/baselines_<regime>.json with mean±95%CI over seeds. Print a table of
  MRR per method per regime so I can see the R0->R3 drop.
  ```

### Unit 6 — `scripts/run_kge.py`  `EXTEND scripts/kge_benchmark.py` · 🔴 A-only (GPU)
- **Runs on A only** — this is the one thing no other machine can do. It's the critical path;
  start it as soon as Unit 4 lands and keep the GPU queue full.
- **Purpose:** TransE/DistMult/ComplEx/RotatE, tuned, ≥3 seeds, evaluated over R0/R2/R3 with the
  **same** `lib_eval` harness as the baselines.
- **In:** split files + `lib_eval` (+ PyKEEN/torch).  **Out:**
  `data/processed/results/kge_<model>_<regime>_seed<k>.json` and a merged `benchmark_results.json`.
- **Does (edits to the existing script):** (a) train on the regime's train graph instead of an
  internal random split; (b) after training, expose `score_fn(g,d)=model.score_hrt` to
  `lib_eval.rank_test_edges` so KGE and baselines share the exact ranking protocol; (c) loop
  seeds {42,1,7}; (d) a small sweep (dim∈{128,256}, epochs∈{100,300}); (e) **fix DistMult's
  collapse** (add regularization / matched loss+negative-sampler); (f) write every run's metrics
  to JSON (the current script only logs KGE metrics — this is the known gap).
- **Acceptance test:** `benchmark_results.json` contains KGE ranking+classification for every
  model × regime × seed with mean±sd; RotatE AUROC ~0.99 but MRR ≪ AdamicAdar under R0
  (the dissociation); DistMult no longer ~0.
- **Prompt:**
  ```
  Extend scripts/kge_benchmark.py into scripts/run_kge.py so it (1) trains each KGE model
  (TransE, DistMult, ComplEx, RotatE) on a regime's training graph from
  data/processed/splits/ (R0 uses train.csv; R2/R3 use train_R2_orthology_blocked.csv), (2)
  scores the held-out gene-disease test edges through scripts/lib_eval.rank_test_edges (the SAME
  harness the baselines use) rather than PyKEEN's internal evaluator, (3) loops seeds 42,1,7 and
  a small sweep dim in {128,256} epochs in {100,300}, (4) fixes DistMult collapsing to ~0 MRR via
  regularization or a matched loss/negative sampler and documents the change, (5) writes every
  run's ranking+classification metrics into data/processed/results/ and a merged
  benchmark_results.json with mean±sd. Keep GPU support, fixed seeds, config.py paths. Show me
  the diff before any long training run.
  ```

---

# STAGE 3 — The method (C4) + validation (C5)

### Unit 7 — `scripts/build_hybrid.py`  `NEW` · 🟣 split A+B
- **Split across machines:** the **KGE-scoring/embedding part runs on A (🔴 GPU)**; the
  **Adamic-Adar candidate-generation + calibrator training runs on B (🟠 CPU)**. They hand off
  via files. If A does both, run the AA part when the GPU is idle.
- **Purpose:** *your* method — a leakage-robust hybrid that stays strong under R3, especially for
  low-degree rare-disease nodes where orthology evidence is thin and KGE fails hardest.
- **Idea:** two stages. (1) **Candidate generation** with Adamic-Adar on the R2/R3 graph (cheap,
  high recall). (2) **Re-ranking** the top candidates with a KGE score and/or a small
  logistic/GBM calibrator over features [AA score, CN count, KGE score, gene degree, disease
  degree, rare-disease flag]. Train the calibrator on `valid.csv`, evaluate on `test.csv`.
- **Out:** `data/processed/results/hybrid_<regime>.json` + a trained calibrator.
- **Win condition (from PAPER_SCOPE §8):** best or tied-best MRR under **R3**, with a paired
  bootstrap p-value vs the best baseline. It does NOT need to win under R0.
- **Acceptance test:** hybrid R3 MRR ≥ best baseline R3 MRR, and the improvement is on the
  low-degree/rare-disease subset. If it doesn't beat baselines under R3, report that honestly
  and simplify the claim (the leakage finding stands regardless).
- **Prompt:**
  ```
  Write scripts/build_hybrid.py: a two-stage gene-disease predictor. Stage 1: Adamic-Adar on the
  regime training graph produces top-N candidate diseases per gene (reuse predict_adamic_adar.py
  logic). Stage 2: re-rank candidates with a scikit-learn calibrator (LogisticRegression or
  GradientBoosting) over features [AA score, common-neighbor count, KGE score from run_kge
  embeddings, gene degree, disease degree, is_rare_disease]. Train the calibrator on valid.csv,
  evaluate on test.csv through scripts/lib_eval. Run for R0 and R3. Report MRR/Hits@k with 95% CI
  and a paired bootstrap p-value vs the best baseline, overall and on the low-degree
  (rare-disease) subset. Seed 42.
  ```

### Unit 8 — `scripts/loso_validation.py`  `NEW` · 🟠 B (heavy CPU) + 🔴 A for the KGE retrain
- **Runs on B** for the Adamic-Adar retrain + scoring; the **"retrain best KGE" step needs A (GPU)**.
  Simplest: B does the AA-LOSO fully; A re-runs the one KGE model on the source-removed graph.
- **Purpose:** prospective-surrogate test #1 — predict edges a withheld source uniquely provided.
- **Does:** remove all edges whose sole source is `<SOURCE>` (start Gene2Phenotype: 2,139
  sole-source; then Orphadata: 4,654) using the `dataset_sources` column; retrain AA + best KGE;
  measure precision@{10,50,100} / recall on exactly those removed edges.
- **Out:** `data/processed/validation/loso_<source>.json`
- **Acceptance test:** precision@k > chance and reported for AA vs KGE; a real "predicted data it
  never saw" number.

### Unit 9 — `scripts/newer_release_diff.py`  `NEW` · 🟢 A or C
- **Purpose:** prospective-surrogate test #2 + the MONDO label file R3 needs.
- **Does:** download the *current* Monarch/Orphanet release; extract (a) gene–disease
  associations that appear ONLY in the newer release (= a prospective test set) and (b) a
  `mondo_id<TAB>label` table (feeds Unit 2's R3 tautology filter). Measure how many of your top-k
  predictions the newer release confirms.
- **Out:** `data/processed/validation/newer_release_recovered.json`, `data/registry/mondo_labels.tsv`
- **Acceptance test:** produces a non-empty label file and a recovery count; note the release
  versions/dates for the provenance table.

### Unit 10 — `scripts/literature_curate.py`  `NEW` + manual · 🟢 A or C
- **Purpose:** the curated top-k literature table (Table 4b), honestly de-tautologized.
- **Does:** take the top ~30 hybrid/AA novel predictions, resolve gene symbol (HGNC file) +
  disease label (Unit 9), auto-flag eponymous/obsolete, emit a worksheet for **manual** PubMed/OMIM
  checking. **Freeze this list before checking** (no cherry-picking — see WORK_SPLIT §5.4).
- **Out:** `data/processed/validation/literature_worksheet.csv` → you fill the evidence column.
- **Acceptance test:** ≥5 non-tautological, literature-confirmed novel candidates with citations.

---

# STAGE 4 — Figures & tables 🟢 A or C

### Unit 11 — `scripts/make_figures.py`  `NEW` · 🟢 A or C
- **Figure 1 (the money shot):** two panels sharing method order — AUROC (hard neg) vs MRR —
  showing KGE tall-left/tiny-right, AA the reverse.
- **Figure 2 (the thesis):** MRR per method across R0→R1→R2→R3 (lines dropping), with 95% CI
  bands — the leakage collapse.
- **Figure 3:** leakage decomposition (orthology-only vs hub-only vs tautology-only contribution).
- **Out:** `figures/fig1..3.{png,pdf}` at 300 dpi, colorblind-safe.
- **Acceptance test:** figures regenerate from the results JSONs with zero manual editing.

### Unit 12 — `scripts/make_tables.py`  `NEW` · 🟢 A or C
- Assembles Table 1 (graph_stats), Table 2 (benchmark mean±CI), Table 3 (R0–R3 drop), Table 4
  (LOSO/newer-release/literature) into markdown + LaTeX straight from JSON.
- **Acceptance test:** every table cell traces to a JSON file; no hand-typed numbers.

---

# STAGE 5 — Release & manuscript 🟢 A or C

### Unit 13 — `run_all.ps1` (+ `run_all.sh`) and Zenodo  `NEW` · 🟢 A or C
- **Purpose:** one-command reproduction on a clean clone; permanent DOIs.
- **Does:** ordered calls graph_stats → build_deleaked_splits → run_baselines → run_kge →
  build_hybrid → loso → figures → tables, with expected runtimes + hardware noted. Then upload
  data + tag a GitHub release → Zenodo DOIs for data and code.
- **Acceptance test:** on a machine that didn't build them, one command reproduces Table 2 +
  Figure 1.

### Unit 14 — Manuscript  · 🟢 A or C (or a dedicated writer)
Follow `PUBLICATION_ROADMAP.md` §6 (IMRaD) and the prompt library §7. Draft section-by-section,
verify every number against the JSONs, run the adversarial "Reviewer 2" pass (§7.7), get two
human readers + a mentor/co-author. Preprint to bioRxiv the day it's done; submit to JEI (parallel
guaranteed track) and BMC Bioinformatics / PLOS ONE (primary). **This is the ideal task for a
dedicated C (or second person) to run in parallel while A trains.**

---

## The single ordered checklist (print this)

Tags: 🔴 = A only (GPU) · 🟠 = B (heavy CPU) · 🟢 = A or C (light) · 🟣 = split A+B

- [ ] 1. `graph_stats.py` → graph_stats.json  🟢 A/C
- [ ] 2. `build_deleaked_splits.py` real run → splits/ + hashes + upload  🟢 A/C (one owner)  ✅ built
- [ ] 3. provenance versions/dates/licenses filled  🟢 A/C
- [ ] 4. `lib_eval.py` + Random-control passes  🟢 A/C  ← A & B blocked until this lands
- [ ] 5. `run_baselines.py` → baselines_R0..R3.json  🟠 B
- [ ] 6. `run_kge.py` → KGE in benchmark_results.json, DistMult fixed, ≥3 seeds  🔴 A
- [ ] 7. `build_hybrid.py` → hybrid best-or-tied under R3  🟣 A+B
- [ ] 8. `loso_validation.py` → precision@k  🟠 B (+🔴 A for KGE retrain)
- [ ] 9. `newer_release_diff.py` → recovery + mondo_labels.tsv → re-run Unit 2 with R3 filter  🟢 A/C
- [ ] 10. `literature_curate.py` → ≥5 confirmed novel candidates  🟢 A/C
- [ ] 11. `make_figures.py` → fig1–3  🟢 A/C
- [ ] 12. `make_tables.py` → tables 1–4  🟢 A/C
- [ ] 13. `run_all.ps1` + Zenodo DOIs  🟢 A/C
- [ ] 14. manuscript → preprint → submit  🟢 A/C

**Minimum viable team:** if it's just you, run **A + B** (A wears the C hat for all 🟢 work during
GPU downtime; B grinds the 🟠 heavy CPU jobs). A separate C only buys you parallel writing/plotting.

*Next action: whoever owns the backbone (A is easiest — it has the data) runs Unit 2 for real
(drop `--dry-run`), uploads `splits/`, then builds Unit 1 (`graph_stats.py`) and Unit 4
(`lib_eval.py`). Units 5 (B) and 6 (A) start the moment `lib_eval.py` passes its Random-control
self-test.*
