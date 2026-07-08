# WORK_SPLIT_3_MACHINES.md — Running this across 3 computers, and making it the best

*Companion to `PAPER_SCOPE.md`. This is the operations manual: how three machines divide the
work without colliding, and the practices that separate a paper that gets accepted from one
that gets desk-rejected.*

---

## 0. The one rule that makes multi-machine work possible

**Every machine must compute on the byte-identical frozen graph and the byte-identical
evaluation splits.** If machine A's KGE model and machine B's Adamic-Adar are evaluated on
different held-out edges, the numbers are not comparable and the whole benchmark is worthless.

So there are exactly **two shared "backbone" artifacts** that get built once, hashed, and
distributed to all three machines before anyone runs a method:

1. **The frozen graph:** `data/processed/edges_clean_integrated.csv` (record its `sha256`).
2. **The frozen splits:** train/valid/test edge lists **plus** the R0–R3 leakage-blocked
   variants, produced by `scripts/build_deleaked_splits.py` (Machine C builds this first).

Until both exist and are distributed, machines A and B are blocked on method-running (but can
do setup, environment, and their non-dependent analyses).

---

## 1. Shared infrastructure (set up on day 1, all three machines)

| Item | Choice | Notes |
|---|---|---|
| **Code sync** | One GitHub repo, three branches (`machine-a`, `machine-b`, `machine-c`) | Merge to `main` via PRs. Never two machines editing the same file on `main`. |
| **Big-data sync** | Zenodo (preferred, gives a DOI) or Google Drive / a shared NAS | The 360 MB graph + splits live here, NOT in git. Add to `.gitignore`. |
| **Integrity** | `sha256sum` on the graph + each split file, checked into git as `ARTIFACT_HASHES.txt` | Every machine verifies the hash before running. One mismatched byte = incomparable results. |
| **Environment** | Identical `requirements.txt` + pinned seeds (42, 1, 7) | Machine A also needs CUDA 12.8+ for the RTX GPU. B/C are CPU-only fine. |
| **Output namespacing** | Each machine writes ONLY to its own subfolder | e.g. `benchmark/kge/` (A), `benchmark/baselines/` + `leakage/` (B), `validation/` + `figures/` (C). No collisions → clean git merges. |
| **Daily sync** | 10-min standup or a shared checklist in `PROGRESS.md` | Who's running what, what's blocked, what landed. |

**Coordination discipline:** small result files (JSON/CSV tables, figures) go in git so everyone
sees them; large artifacts (embeddings, the graph) go in the data channel with a hash. If a
machine regenerates a backbone artifact, it bumps the hash file and pings the others to re-pull.

---

## 2. Machine roles (play to each machine's strength)

### Machine A — the GPU box (RTX 5070). *Role: heavy compute / the critical path.*
Owns everything that needs the GPU, because that's the scarce resource and the wall-clock
bottleneck.
- **Phase 2:** train all KGE models (TransE, DistMult, ComplEx, RotatE) × 3 seeds × the small
  hyperparameter sweep, on the frozen splits. Save per-run ranking + classification metrics to
  `benchmark_results.json`.
- Fix DistMult's collapse; report the best config found.
- **Phase 3/C3:** run each trained KGE model through the R0–R3 evaluation regimes.
- **C4:** train the embedding component of the hybrid predictor.
- Output → `benchmark/kge/`, `hybrid/embeddings/`. This is a serial GPU queue — start it early.

### Machine B — the CPU/RAM box. *Role: topology, leakage analysis, null models.*
Everything graph-algorithmic that does NOT need a GPU (networkx/scipy/numpy).
- **Phase 2 baselines:** Random, Common-Neighbors, Adamic-Adar on the frozen splits.
- **C1/C3 — the leakage analyses (the paper's core):** implement the orthology-path blocker,
  the hub blocker, tautology/eponym detection, and run every *topological* method through
  R0–R3. This is the biggest single work item and B owns it.
- **Phase 4 LOSO:** retrain baselines with a source removed; precision@k on sole-source edges.
- Null-model tests (already scripted) + graph characterization stats.
- Output → `benchmark/baselines/`, `leakage/`, `validation/loso/`, `analysis/null/`.

### Machine C — laptop / coordination. *Role: data steward, integration, writing.*
The lightest compute, the most coordination.
- **Backbone (do first, everyone waits on this):** build `scripts/build_deleaked_splits.py`
  → the frozen train/valid/test + R0–R3 variants; hash and publish them.
- `graph_stats.json` (Table 1); provenance-table finalization (versions/dates/licenses).
- **Phase 4:** newer-release download + diff (prospective-surrogate test set); curate the
  literature-validation table (de-tautologized).
- **Phase 5:** Zenodo upload + DOIs, one-command reproduction script, `DATA.md`, LICENSE.
- **Phase 6:** assemble figures, write the manuscript, run the adversarial pre-review, merge
  all machines' result files into the final tables.
- Output → `data/processed/*` (backbone), `figures/`, `manuscript/`.

> **If Machine B or C also has a usable GPU:** shard Machine A's `(model × seed × hp)` grid
> across them — the sweep is embarrassingly parallel. Split the config list into 2–3 disjoint
> chunks, one per GPU machine, each writing to `benchmark/kge/<machine>/`, then C concatenates.
> If only A has a GPU, A is the critical path — so B and C must keep A's queue *full* and do
> all non-GPU work in parallel, so total wall-clock ≈ A's queue length, not the sum.

---

## 3. The dependency-ordered job manifest

Run top-to-bottom; items on the same row run in parallel on different machines.

| Stage | Machine A (GPU) | Machine B (CPU) | Machine C (coord) | Blocks on |
|---|---|---|---|---|
| S0 setup | env + CUDA, pull graph | env, pull graph | build repo, hash graph, `PROGRESS.md` | — |
| S1 backbone | *(idle → start baseline env checks)* | *(prototype orthology blocker on a sample)* | **build & publish frozen splits R0–R3** | S0 |
| S2 methods | **KGE × seeds × sweep** | **baselines + leakage R0–R3** | graph_stats, provenance finalize | S1 |
| S3 analysis | KGE through R0–R3 regimes | LOSO + null models | newer-release diff, lit-validation | S2 |
| S4 method C4 | train hybrid embeddings | hybrid topology seeding + eval | figures from S2/S3 tables | S3 |
| S5 assemble | hand results to C | hand results to C | merge tables, write manuscript, Zenodo | S4 |
| S6 submit | — | — | preprint + submit + cover letter | S5 |

**Critical path = S1 (C) → S2 (A's KGE queue) → S3 → S4 → S5.** Keep the KGE GPU queue
saturated from the moment S1 lands; that queue, not the analysis, is what sets the timeline.

---

## 4. Concurrency hazards to avoid (learned the hard way)

- **Two machines regenerating a backbone artifact** → different hashes → incomparable results.
  *Rule:* only Machine C regenerates splits/graph; A and B consume, never regenerate.
- **Editing the same script on two branches** → merge conflicts on the code that matters.
  *Rule:* one owner per script (A owns `kge_benchmark.py`, B owns `build_deleaked_splits.py`'s
  consumers / leakage scripts, C owns figures + manuscript). Coordinate in `PROGRESS.md`.
- **Silent seed drift** → machine A uses seed 0, B uses 42 → the "same" split differs.
  *Rule:* seeds live in `config.py`, imported everywhere, never hard-typed per machine.
- **Committing 300 MB CSVs to git** → repo becomes unusable. *Rule:* big files → data channel
  only; enforce with `.gitignore` + a pre-commit size check.

---

## 5. "Make it the best" — the excellence practices that clear peer review

These are what turn a sound project into an *accepted* paper. Do all of them.

### 5.1 Pre-register the analysis (the highest-credibility move)
Before you run R3 and look at the final numbers, **write down the analysis plan and timestamp
it** — `PAPER_SCOPE.md` is most of it; add a short `PROTOCOL.md` stating the exact metrics,
regimes, seeds, and the hybrid's win-condition, then commit it (git timestamp) or post it to
OSF. This proves you didn't p-hack the leakage number after seeing it. Reviewers trust
pre-registered benchmarks far more, and it costs you nothing.

### 5.2 Statistics done properly, not eyeballed
- **≥3 seeds, report mean ± 95% CI** for every metric in every table.
- **Significance test the key comparisons:** paired bootstrap over per-edge reciprocal ranks
  for "AA vs KGE" and "hybrid vs best-baseline under R3." Report p or CI of the difference.
- **Negative controls that must pass:** Random method ≈ chance (AUROC 0.50, MRR ≈ H₅₁/51);
  λ_max(Laplacian) ≈ max degree; a shuffled-label run must destroy performance. If a control
  fails, a bug is hiding — fix before believing any result.

### 5.3 The leakage claim must be airtight (it's the whole paper)
- Show the R0→R3 drop **per method** with CIs, and whether the **method ranking changes**.
- Include a **path-level illustration**: for a handful of held-out edges, show the actual
  orthology→phenotype path the model exploited. One worked example makes the abstract concept
  concrete for reviewers.
- **Ablate the pieces of R3** (orthology-only vs hub-only vs tautology-only) so readers see
  which leakage source dominates — that decomposition is itself a result.

### 5.4 Don't cherry-pick the biology
- Freeze the top-k prediction list **before** literature checking; report precision@k on that
  frozen list, including the misses. Curating after looking is how face-validity becomes fraud.
- Exclude eponymous/tautological and obsolete-MONDO hits from the *headline* number, but
  *report* them separately and transparently.

### 5.5 Get a senior co-author / mentor
A named domain expert (a professor, a lab you can email, a teacher with a bio/CS PhD) who
reads the biology and co-signs **dramatically** lowers desk-reject odds at Tier-A venues and
catches errors you can't see. This is standard scientific practice, not a shortcut. Start
looking now; offer co-authorship for genuine intellectual contribution.

### 5.6 Adversarial self-review before submission
Run the "Reviewer 2" prompt (roadmap §7.7) on the full draft. Fix everything it raises. Then
have your two human readers (one bio, one non-expert) do the same. If a non-expert can't
restate your finding after reading the abstract, the abstract is not done.

### 5.7 Reproducibility as a feature, not a chore
Hashed frozen artifacts, one-command reproduction on a clean clone, Zenodo DOIs for data +
code, fixed seeds, pinned environment. At GigaScience/Database/PLOS this is *scored*; make it
your strength since the infrastructure is already 80% there.

### 5.8 Write limitations *proactively*
Put the entire honesty ledger (roadmap §1d) in the Limitations section, in your own words,
before a reviewer finds it. "We disclose that…" reads as rigor; the same fact discovered by a
reviewer reads as concealment. Honesty is not just ethical here — it is the winning strategy
for this specific paper, because honesty *is* the contribution.

---

## 6. Definition of done (the whole project)

- [ ] Backbone frozen + hashed; all three machines verified identical inputs.
- [ ] Table 1 (graph), Table 2 (benchmark, mean±CI), Figure 1 (R0→R3 leakage drop),
      Table 3 (leakage decomposition), Table 4 (LOSO/newer-release/literature).
- [ ] Hybrid (C4) best-or-tied under R3, with a significance test.
- [ ] Pre-registered protocol timestamped; negative controls pass; no fabricated numbers.
- [ ] ≥5 non-tautological literature-confirmed rare-disease candidates.
- [ ] Zenodo DOIs; one-command repro verified on a fresh clone by a machine that didn't build it.
- [ ] Senior co-author on board; adversarial review done; two human readers.
- [ ] bioRxiv preprint posted; primary venue submission formatted and sent; JEI version in parallel.

*First concrete step across the three machines: Machine C builds
`scripts/build_deleaked_splits.py` (the R0–R3 backbone) while A sets up the KGE sweep harness
and B prototypes the orthology-path blocker on a 100k-edge sample. Nothing else is comparable
until the backbone splits exist and are hashed.*
