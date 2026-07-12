# BUILD_GUIDE.md (audit plan) — build the leakage audit, unit by unit, in order

*The construction manual for the pivoted paper. Each unit: what it is, **who runs it**, inputs →
outputs, the acceptance test that proves it's correct, and its origin (reuse / extend / new /
cut) relative to the old build guide. Do them top to bottom.*

**Status legend:** `REUSE` = already built, run it · `EXTEND` = modify an existing built script ·
`NEW` = create it · `CUT` = deleted from scope (do not build).

**Owner tags:** 🔴 A (KGE, GPU/Colab) · 🟠 B (heavy CPU) · 🟢 A/C (light: data, figures, writing).

**What you already built maps in like this:** old Units 1,3 → drop-in REUSE; old Unit 2 (splits) and
Unit 4 (`lib_eval`) → EXTEND (regimes redefined); old Unit 6 (`run_kge`) → EXTEND (trim to 2 models).
Old Units 7 (hybrid), 8 (LOSO), 9 (newer-release) → **CUT** (never built — pure savings).

---

## Dependency graph (build order + week)

```
STAGE 1 backbone (W1)              STAGE 2 regimes+methods (W1–4)        STAGE 3 robustness+stats (W5–6)
  1. graph_stats.py .... 🟢 REUSE    4. lib_eval regime map . 🟢 EXTEND    9.  hetionet_audit.py . 🟠 NEW
  2. build_deleaked_    🟢 EXTEND    5. run_baselines.py .... 🟠 NEW      10. make_figures.py ... 🟢 NEW
     splits.py (regimes)             6. R1 redundancy ....... 🟠 NEW      11. make_tables.py .... 🟢 NEW
  3. provenance ........ 🟢 REUSE    7. R2 degree-null ...... 🟠 NEW
                                     8. run_kge.py .......... 🔴 EXTEND   STAGE 4 release (W7–8)
                                        (2 models, filtered)             12. run_all + Zenodo .. 🟢 NEW
                                                                         13. manuscript ....... 🟢
```

Nothing in Stage 2's methods can start until **Unit 2 (regime files)** and **Unit 4 (lib_eval
regime map)** exist and are hashed. Those two come first and have a single owner (C).

**Backbone numbers to expect (from the built splits):** 42,288 target human gene→disease edges →
train 5,852,083 / valid 4,228 / test 4,228. These are unchanged by the pivot; only the R1/R2/R3
*variants* differ.

---

## Time estimates

Two clocks, because they're very different: **Build** = writing/adapting + debugging the code (your
active time); **Compute** = the run, mostly unattended (start it and go do the next unit). Estimates
assume **one person on home hardware** at a student's part-time pace; ranges cover "smooth" → "hit a
snag." REUSE units are near-zero build. Compute overlaps build, so wall-clock is far less than the sum
— especially split across 3 machines (B grinds compute while A trains KGE and C writes).

| Unit | Type | Build (active) | Compute (run) | Week | Notes |
|---|---|---|---|---|---|
| Setup (per machine) | — | 1–2 h each | +Hetionet download 15–30 min | W0 | one-time; +`pip install python-igraph` |
| 1 `graph_stats` | REUSE | ~1 h (Hetionet adapter) | 5–15 min | W1 | script already built |
| 2 `build_deleaked_splits` regimes | EXTEND | **1–2 days** | mins–1 h | W1 | the backbone; C owns it |
| 3 provenance | REUSE | ~0.5 day (manual lookup) | — | W1 | tedious, not hard |
| 4 `lib_eval` regime map | EXTEND | 0.5–1 day | seconds (self-test) | W1 | **blocks A & B until done** |
| 5 `run_baselines` (5 × R0–R3 × 3 seeds) | NEW | 1–2 days | **4–12 h on B** | W1–2 | biggest CPU scoring job |
| 6 R1 redundancy | NEW | 0.5–1 day | mins | W2 | can live inside Unit 2 |
| 7 R2 degree null | NEW | ~1 day | swaps mins–1 h/seed | W2–3 | igraph swaps are fast; re-scoring on the null is folded into Unit 5 |
| K `run_kge` (2 models, filtered) | EXTEND | 1–2 days | 6–15 h GPU, unattended | W4 | off critical path; Colab OK |
| 9 `hetionet_audit` | NEW | 1–2 days | 1–2 h | W5 | small graph, fast run |
| 10 `make_figures` | NEW | ~1 day | seconds | W6 | regenerate from JSON |
| 11 `make_tables` | NEW | 0.5–1 day | seconds | W6 | regenerate from JSON |
| 12 `run_all` + Zenodo | NEW | ~1 day | full-pipeline rerun overnight | W7 | proves reproducibility |
| 13 manuscript | — | **1–2 weeks** | — | W7–8 | the real time sink; start early |

**Totals:** roughly **3–4 weeks of active build/write effort** + a few days of mostly-unattended
compute, which is exactly why the calendar plan is **8 weeks** — the slack absorbs compute waits,
debugging, the manuscript, and buffer. **Critical path** is Unit 2 → 4 → 5 → 7 → figures → write;
KGE (Unit K) and Hetionet (Unit 9) run in parallel and never gate the timeline.

---

# STAGE 1 — Backbone (W1) 🟢 A or C

### Unit 1 — `scripts/graph_stats.py`  `REUSE` · 🟢
Already built. Run as-is on the Monarch graph for Table 1; **also run it on Hetionet** for the
second-graph column.
- **Out:** `data/processed/graph_stats.json`, `data/processed/graph_stats_hetionet.json`.
- **Acceptance:** prints `edges: 5860539`, top relation `BIOLINK:INTERACTS_WITH` (~1.9M) on Monarch.

### Unit 2 — `scripts/build_deleaked_splits.py`  `EXTEND` · 🟢 C (ONE owner — regenerates regimes)
Built for the old regimes; **rework the R1/R2/R3 definitions** to the audit set. Keep the reusable
machinery: target-edge extraction (42,288 human gene→disease), train/valid/test freeze, manifest,
hashing, and the "no held-out edge in train" assertion.
- **New regime outputs** (all frozen + hashed):
  - **R0 standard:** `train.csv`, `test.csv` — unchanged.
  - **R1 redundancy-controlled:** `train_R1_redundancy.csv` = train minus edges that leak a test
    edge via an inverse/duplicate/symmetric relation (e.g. the reverse of a held-out pair, or a
    parallel edge under a symmetric relation). Emit the count removed.
  - **R2 degree-controlled:** *this file is produced by Unit 7 (the null), not here* — Unit 2 just
    reserves the manifest slot and records the seed.
  - **R3 orthology-blocked (Monarch only):** `train_R3_orthology_blocked.csv` = train minus
    `ORTHOLOGOUS_TO` and `MODEL_OF` edges. (This is the old R2 file — reuse that logic directly.)
- **Acceptance:** `split_manifest.json` shows the four regimes with their file map + per-regime
  removed-edge counts; the leakage assertion passes for every regime; hashes written to
  `ARTIFACT_HASHES.txt`. Upload `data/processed/splits/` to the data channel.
- **Origin note:** the old R2-orthology-block code is exactly what new R3 needs — copy it over.

### Unit 3 — Provenance finalize  `REUSE`/data · 🟢 (manual)
Run the existing `scripts/provenance_table.py`; manually fill each source's upstream **release
version + download date + license**. Trim any "13-database integration" language — the audit paper
claims a Monarch-derived graph + a released benchmark, not a novel integration.
- **Acceptance:** `source_provenance.md` has no "placeholder" strings.

---

# STAGE 2 — Regime map + methods (W1–4)

### Unit 4 — `scripts/lib_eval.py` regime map  `EXTEND` · 🟢 C (shared library — ONE owner)
Already built and its ranking/metrics/bootstrap are correct — **keep all of that.** Two changes:
1. **Rewrite `load_regime`** to the new map: R0 → `train.csv`/`test.csv`; R1 →
   `train_R1_redundancy.csv`/`test.csv`; R2 → the degree-null train (Unit 7 output)/`test.csv`;
   R3 → `train_R3_orthology_blocked.csv`/`test.csv`. Drop the old hub_filter plumbing (or keep it
   dormant — R1 is now redundancy, not hubs).
2. **No new metric code needed for Jaccard/Preferential-Attachment** — they're just `score_fn`s
   passed by `run_baselines.py`; the harness is already generic. Confirm the docstring reflects the
   new regimes.
- **Acceptance (unchanged, must still pass):** a Random scorer gives MRR ≈ chance (~0.09 at
  n_neg=50) and AUROC ≈ 0.50. If not, the harness has a bug — fix before trusting any result.

### Unit 5 — `scripts/run_baselines.py`  `NEW` (fold in `bench_baselines.py` + `predict_adamic_adar.py`) · 🟠 B
The primary evidence. Five topological scorers over all regimes.
- **In:** regime files + `lib_eval`. **Out:** `results/baselines/baselines_<regime>.json`.
- **Does:** for each regime R0–R3 × seeds {42,1,7}: build an undirected adjacency from the regime's
  train graph; define **Random, Common-Neighbors, Adamic-Adar, Jaccard, Preferential-Attachment**
  score_fns (cap shared neighbors with degree > 2000 as in `predict_adamic_adar.py`); evaluate via
  `lib_eval.rank_test_edges` + `ranking_metrics` + `classification_metrics`; collect per-edge
  reciprocal ranks; save mean ± 95% CI.
- **Acceptance:** under R0, Adamic-Adar MRR ≈ 0.69, Hits@10 ≈ 0.83 (matches the existing
  `benchmark_results.json`). MRR **drops** R0 → R3 — that drop is the leakage signal; record it per method.

### Unit 6 — R1 redundancy control  `NEW` (may live inside Unit 2) · 🟠 B
Remove inverse/duplicate/symmetric edges that leak a held-out pair, then re-evaluate.
- **Does:** for each test edge `(g, d)`, drop from train any edge that is (a) the exact reverse
  `(d, g)` under a symmetric/invertible relation, (b) a parallel duplicate of the target relation, or
  (c) a trivially inverse relation type. Produces `train_R1_redundancy.csv` (feeds Unit 5's R1 run).
- **Acceptance:** count of removed edges is reported and non-trivial; R1 baseline MRR ≤ R0 MRR
  (removing leakage can only hurt or hold).

### Unit 7 — R2 degree-preserving permutation null  `NEW` (seed from `null_model_tests.py`) · 🟠 B
**The single heaviest job.** Attributes performance to degree alone.
- **Does:** build the undirected simple graph in **igraph** from the R0 train edges; run
  degree-preserving **double-edge swaps** (≈ 10× |E| swaps) to randomize topology while holding every
  node's degree fixed; write `train_R2_degree_null.csv`. Run ≥3 independent swap seeds. Assert the
  degree sequence is byte-identical before/after (this is the correctness guarantee).
- **In:** R0 train edges. **Out:** `train_R2_degree_null_seed<k>.csv`, `results/null/degree_null_meta.json`.
- **Compute:** checkpoint every N million swaps; keep int32 ids; this can take hours — start W2.
- **Acceptance:** post-swap degree sequence == pre-swap (assert); baselines scored on R2 collapse
  toward chance for any signal that was pure degree; the residual above chance is "real" structure.
  The **R2 → R3 gap** (add orthology block on top) is what tests whether orthology is distinct from degree.

### Unit K — `scripts/run_kge.py`  `EXTEND` (trim) · 🔴 A (GPU or Colab)
Already built for 4 models — **trim to two and retarget.**
- **Keep:** PyKEEN training loop, seed loop, scoring **through `lib_eval.rank_test_edges`** (same
  harness as baselines), per-run JSON output.
- **Change:** train only **TransE + ComplEx** (dim 64, standard cited hyperparameters + one small
  sweep on R0); train on the **filtered subgraph** (Gu 2024, ~11% of Monarch) for speed; drop the
  DistMult-collapse fix and RotatE. Add a Colab/Kaggle entry path (upload splits → run → download JSON).
- **In:** regime train files + `lib_eval`. **Out:** `results/kge/kge_<model>_<regime>_seed<k>.json`.
- **Acceptance:** JSON for TransE + ComplEx × {R0,R2,R3} × 3 seeds with mean±sd; each scored by the
  same harness as the baselines; RotatE/DistMult absent (out of scope now).
- **Runs in parallel, off the critical path** — start whenever the regime files land.

---

# STAGE 3 — Robustness + stats (W5–6)

### Unit 9 — `scripts/hetionet_audit.py`  `NEW` · 🟠 B
Zero-integration robustness: rerun the audit on a second, pre-built graph.
- **Does:** load Hetionet edges (`data/external/hetionet/`), map its Gene–Disease metaedge to the
  same `(gene, disease)` target format, build R0 + R1 (redundancy) + R2 (degree null); run the five
  baselines through `lib_eval`. **Skip R3 (orthology) — Hetionet lacks cross-species orthology;
  state this explicitly.** Do **not** merge Hetionet with Monarch.
- **Out:** `results/hetionet/baselines_<regime>.json`.
- **Acceptance:** the R0→R2 drop reproduces qualitatively on Hetionet (de-leaking drops methods on a
  second graph too) — this is the generality claim.

### Unit 10 — `scripts/make_figures.py`  `NEW` · 🟢 A/C
- **Figure 1:** AUROC vs MRR per method (the dissociation — KGE tall-left/tiny-right).
- **Figure 2 (the thesis):** MRR per method across R0→R1→R2→R3 with 95% CI bands — the collapse.
- **Figure 3:** leakage decomposition (redundancy-only vs degree-only vs orthology-only contribution).
- **Out:** `figures/fig1..3.{png,pdf}` at 300 dpi, colorblind-safe.
- **Acceptance:** regenerate from the results JSONs with zero manual editing.

### Unit 11 — `scripts/make_tables.py`  `NEW` · 🟢 A/C
Assemble Table 1 (both graphs), Table 2 (benchmark mean±CI, all methods × R0–R3), Table 3 (R0→R3
drop + permutation p-value + R2→R3 residual), Table 4 (Hetionet) → markdown + LaTeX from JSON.
- **Acceptance:** every cell traces to a JSON file; no hand-typed numbers.

**Optional stretch (only if W6 is on track):** temporal split using two dated Monarch releases —
train on older, test on newer-only edges. Adds a real prospective-flavored regime. Cut if behind.

---

# STAGE 4 — Release & manuscript (W7–8) 🟢 A/C

### Unit 12 — `run_all.ps1` (+ `run_all.sh`) and Zenodo  `NEW` · 🟢
Ordered calls graph_stats → build_deleaked_splits → run_baselines → R2 null → run_kge →
hetionet_audit → figures → tables, with expected runtimes + hardware noted. Then upload the
de-leaked splits + code, tag a GitHub release → **Zenodo DOIs** for data and code.
- **Acceptance:** on a machine that didn't build them, one command reproduces Table 2 + Figure 2.

### Unit 13 — Manuscript  · 🟢 (or a dedicated writer)
Resource-led IMRaD. Draft section-by-section; verify every number against the JSONs; **cite and
distinguish Gu 2024, Ranga 2025, and Alghamdi/Hoehndorf/Robinson 2022** in Related Work (novelty is
incremental methodological + resource — do not oversell). Run the adversarial "Reviewer 2" pass; get
a senior reader. bioRxiv the day it's done; submit to a resource-tier venue (GigaScience-DB /
Database / PeerJ / PLOS ONE / NAR-GB / Bioinformatics Advances).

---

## The single ordered checklist (print this)

Tags: 🔴 A (KGE) · 🟠 B (heavy CPU) · 🟢 A/C (light)

- [x] ~~0. Setup: env + Monarch graph + Hetionet downloaded & hashed~~ ✅ 2026-07-11
- [ ] 1. ~~`graph_stats.py` (Monarch → graph_stats.json)~~ ✅ — still need the **Hetionet** run  🟢  *(REUSE)*
- [x] ~~2. `build_deleaked_splits.py` reworked → R0/R1/R3 files + hashes~~ ✅ 2026-07-11 · R1 −1,661, R3 == old orthology, base split unchanged  🟢 C *(EXTEND)*
- [ ] 3. ~~provenance table generated~~ ✅ — still need versions/dates/licenses filled in  🟢  *(REUSE, manual — your task)*
- [x] ~~4. `lib_eval.py` regime map rewritten; Random-control passes~~ ✅ 2026-07-11 · R0 random MRR 0.088 = chance  🟢 C *(EXTEND)*
- [x] ~~5. `run_baselines.py` → 5 baselines × R0/R1/R3, 3 seeds~~ ✅ 2026-07-12 · **finding: R1 & R3 drops ≈ 0 for topological baselines** (orthology path is ≥3 hops, invisible to 2-hop methods); PA (degree) is best → R2 is decisive  🟠 B *(NEW)*
- [x] ~~6. R1 redundancy control → `train_R1_redundancy.csv`~~ ✅ folded into Unit 2  *(NEW)*
- [ ] 7. R2 degree-permutation null → `train_R2_degree_null_*.csv`  🟠 B  *(NEW, heaviest)*  ← **NEXT & DECISIVE** (only regime that hits PA/degree; make-or-break for the baseline leakage claim)
- [ ] K. `run_kge.py` trimmed → TransE+ComplEx, filtered subgraph, JSON  🔴 A  *(EXTEND)*
- [ ] 9. `hetionet_audit.py` → second-graph tables  🟠 B  *(NEW)*
- [ ] 10. `make_figures.py` → fig1–3  🟢  *(NEW)*
- [ ] 11. `make_tables.py` → tables 1–4  🟢  *(NEW)*
- [ ] 12. `run_all.ps1` + Zenodo DOIs  🟢  *(NEW)*
- [ ] 13. manuscript → preprint → submit  🟢

**CUT (do not build):** hybrid predictor, LOSO validation, newer-release recovery. Literature check
shrinks to a light spot-check of a few top predictions (optional, one worksheet).

**Minimum viable team:** if it's just you, run **A + B** — B grinds the heavy CPU jobs (baselines +
degree null + Hetionet), A trains the two KGE models (locally or on Colab) and wears the C hat for
all 🟢 work during downtime. A separate C only buys parallel writing/plotting.

*Next action: Machine C reworks Unit 2's regimes and Unit 4's `lib_eval` map, uploads the new
regime files, while B prototypes the Unit 7 degree swap on a 100k sample and A preps the filtered
subgraph. Units 5–7 (B) and K (A) start the moment the regime files + `lib_eval` land.*
