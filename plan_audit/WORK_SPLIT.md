# WORK_SPLIT.md (audit plan) — running the leakage audit across 3 computers

*Companion to `SETUP_ENVIRONMENTS.md` and `BUILD_GUIDE.md` in this folder. Operations manual for
the pivoted paper: how three machines divide the audit without colliding. Key change from the old
plan: **KGE shrank (2 small models on a filtered subgraph, optionally on Colab), so the GPU box A
is no longer the critical path — the CPU/RAM box B is** (the degree-permutation null over 5.86M
edges is the single heaviest job).*

---

## 0. The one rule that makes multi-machine work possible

**Every machine computes on the byte-identical frozen graph and the byte-identical regime files.**
If A's KGE and B's Adamic-Adar are scored on different held-out edges, nothing is comparable.

Two shared "backbone" artifacts, built once, hashed, distributed **before** anyone runs a method:
1. **The frozen graphs:** `edges_clean_integrated.csv` (Monarch) + Hetionet — record each `sha256`.
2. **The frozen regimes:** train/valid/test + the R0–R3 variants from
   `scripts/build_deleaked_splits.py` (Machine C owns this — builds it first).

Until both exist and are distributed, A and B are blocked on method-running (but do setup + prototyping).

---

## 1. Shared infrastructure (day 1, all three)

| Item | Choice | Notes |
|---|---|---|
| **Code sync** | one GitHub repo, three branches, merge to `main` via PRs | never two machines editing one file. |
| **Big-data sync** | Google Drive / NAS (→ Zenodo at release) | the two graphs + regime files live here, not in git. |
| **Integrity** | `sha256` on every graph + regime file → `ARTIFACT_HASHES.txt` | verify before every run. |
| **Environment** | pinned `requirements.txt` + **igraph** + seeds `{42,1,7}` | A also needs CUDA *or* a Colab account. |
| **Output namespacing** | each machine writes only to its own subfolder | A→`results/kge/`, B→`results/baselines/`+`results/null/`+`results/hetionet/`, C→`splits/`+`figures/`+`tables/`. |
| **Daily sync** | short standup or a `PROGRESS.md` checklist | who's running what, what's blocked, what landed. |

Small result files (JSON/CSV/figures) go in git; large artifacts go in the data channel with a hash.

---

## 2. Machine roles (rebalanced for the audit)

### Machine B — CPU/RAM box. *Role: the heavy machine and the critical path.*
Owns the biggest compute now. Everything graph-algorithmic, no GPU.
- **All baselines** (Random, Common-Neighbors, Adamic-Adar, **Jaccard**, **Preferential-Attachment**)
  across R0–R3, ≥3 seeds.
- **R1 redundancy control** (remove inverse/duplicate/symmetric edges that leak the test edge).
- **R2 degree-permutation null** — degree-preserving double-edge swaps via **igraph** over 5.86M
  edges. *This is the single heaviest job in the project; start it early, checkpoint it.* Seed the
  implementation from the existing `scripts/null_model_tests.py`.
- **Hetionet audit** — rerun the whole baseline harness on the second graph (Monarch-only steps like
  orthology R3 are skipped there; state this).
- Output → `results/baselines/`, `results/null/`, `results/hetionet/`.

### Machine A — GPU box (or Colab). *Role: KGE only — now a side quest, not the bottleneck.*
- Train **TransE + ComplEx** (two families), dim 64, standard cited hyperparameters + one small
  sweep on R0, seeds `{42,1,7}`, on the **filtered subgraph** (Gu 2024) for speed.
- Score every trained model through the **same `lib_eval` harness** the baselines use (not PyKEEN's
  internal evaluator) so KGE and baselines are directly comparable.
- Runs on the local RTX **or** free Colab/Kaggle T4 — each run is tens of minutes.
- Output → `results/kge/`. Because it's light, **A also wears the C hat** (figures/writing) during any downtime.

### Machine C — laptop / coordination. *Role: data steward, regimes, figures, writing.*
- **Backbone (do first, everyone waits):** `scripts/build_deleaked_splits.py` → frozen
  train/valid/test + **the new R0–R3 regime files**; hash and publish them.
- `graph_stats.json` (Table 1) for **both** graphs; provenance finalize (versions/dates/licenses).
- Figures + tables from the results JSONs; the degree-vs-orthology decomposition figure.
- Release: Zenodo DOIs, one-command reproduction, `DATA.md`, LICENSE; manuscript + adversarial review.
- Output → `data/processed/splits/` (backbone), `figures/`, `tables/`, `manuscript/`.

> **If B or C also has a GPU:** shard A's tiny `(model × seed)` grid across them — it's
> embarrassingly parallel. But KGE is no longer the bottleneck, so the bigger win is a **second
> person on B's degree-null / Hetionet run** or on writing.

---

## 3. Dependency-ordered job manifest (maps to the 8-week timeline)

Run top-to-bottom; same-row items run in parallel on different machines.

| Stage / week | Machine A (KGE) | Machine B (heavy CPU) | Machine C (coord) | Blocks on |
|---|---|---|---|---|
| S0 setup (pre-W1) | env + CUDA/Colab, pull graphs | env + igraph, pull graphs | build repo, hash graphs, `PROGRESS.md` | — |
| S1 backbone (W1) | *(idle → filtered-subgraph prep)* | *(prototype R2 swap on a 100k sample)* | **build & publish R0–R3 regimes** + extend `lib_eval` regime map | S0 |
| S2 baselines (W1–2) | *(idle → C-hat: figures scaffold)* | **5 baselines on R0; then R1 redundancy + R2 null** | graph_stats (both graphs), provenance | S1 |
| S3 R3 + decomp (W3) | *(idle)* | **R3 orthology-block; R2→R3 residual + permutation p-value** | assemble R0–R3 baseline tables | S2 |
| S4 KGE (W4) | **TransE + ComplEx across regimes** (filtered) | *(help: Hetionet prep)* | merge KGE into benchmark tables | S1 (regimes) |
| S5 Hetionet (W5) | *(idle → C-hat)* | **full baseline audit on Hetionet** | second-graph tables | S4 harness stable |
| S6 stats/figures (W6) | hand KGE results to C | bootstrap CIs, paired tests, optional temporal split | **figures + CIs** | S3–S5 |
| S7 write/release (W7) | — | — | manuscript, Zenodo DOIs, one-command repro | S6 |
| S8 submit (W8) | — | — | preprint + submit + buffer | S7 |

**Critical path = S1 (C regimes) → S2/S3 (B's baselines + degree null) → S6 → S7.** KGE (S4) runs
in parallel off to the side and never gates the timeline. Keep B saturated.

---

## 4. Concurrency hazards (learned the hard way)

- **Two machines regenerating a backbone artifact** → different hashes → incomparable results.
  *Rule:* only C regenerates graphs/regimes; A and B consume, never regenerate.
- **Editing the same script on two branches** → merge conflicts. *Rule:* one owner per script.
- **Seed drift** (A uses 0, B uses 42) → the "same" split differs. *Rule:* seeds live in `config.py`,
  imported everywhere, never hard-typed.
- **Committing 300 MB CSVs to git** → repo unusable. *Rule:* big files → data channel only.
- **igraph vs networkx node-id mismatch** on B → the null swaps a different graph than the baselines
  score. *Rule:* build the null from the *same* edge list the regime files come from, and assert the
  degree sequence is preserved after swapping.

---

## 5. Excellence practices that clear peer review (adapted to the audit)

The old plan's hybrid-win-condition and heavy biological validation are gone. These remain and are
what get an audit/resource paper accepted:

### 5.1 Pre-register the analysis
Before running R3 and looking at the final numbers, timestamp a short `PROTOCOL.md` stating the exact
metrics, regimes, seeds, and the **decision rule for "orthology is distinct from degree"** (e.g. the
R2→R3 residual MRR drop and its permutation p-value threshold). Commit it (git timestamp) or post to
OSF. This proves you didn't pick the leakage story after seeing the number.

### 5.2 Statistics done properly, not eyeballed
- ≥3 seeds, **mean ± 95% CI (bootstrap)** for every metric in every table.
- **Paired bootstrap** over per-edge reciprocal ranks for method comparisons; **permutation p-value**
  for the degree control.
- **Negative controls that must pass:** Random ≈ chance (AUROC 0.50, MRR ≈ 1/((n_neg+1))); a
  shuffled-label run destroys performance. If a control fails, a bug is hiding — fix before believing anything.

### 5.3 The audit claim must be airtight (it's the whole paper)
- Show the R0→R3 drop **per method** with CIs, and whether the **method ranking changes**.
- **Decompose the drop** (redundancy-only vs degree-only vs orthology-only) — that decomposition *is*
  the result. The R2→R3 residual is the one analysis carrying the secondary novelty; report a null
  honestly in one sentence ("orthology leakage is largely a manifestation of degree bias").
- One **worked path example**: for a handful of held-out edges, show the actual
  orthology→phenotype path a method exploited. Makes the abstract concept concrete.

### 5.4 Reproducibility as the headline feature (this *is* the contribution)
Hashed frozen artifacts, one-command reproduction on a clean clone, Zenodo DOIs for data + code,
fixed seeds, pinned environment, SHA-256'd de-leaked splits. Database/GigaScience/PeerJ **score** this.

### 5.5 Honest positioning + proactive limitations
Cite and explicitly distinguish the close prior art (Gu 2024, Ranga 2025, **Alghamdi/Hoehndorf/Robinson
2022** ortholog-bias) in Related Work — novelty is *incremental methodological + resource*, not a
discovery. Put the whole honesty ledger in Limitations, in your own words, before a reviewer finds it.

### 5.6 Adversarial self-review + a senior reader
Run a "Reviewer 2" pass on the full draft; fix everything. A named domain expert who co-signs the
biology dramatically lowers desk-reject odds — start looking now.

---

## 6. Definition of done (whole project)

- [ ] Both graphs frozen + hashed; three machines verified identical inputs.
- [ ] R0–R3 regime files built, hashed, SHA-256'd, published.
- [ ] Table 1 (both graphs), Table 2 (benchmark, mean±CI, all methods × R0–R3), Figure 2 (R0→R3
      drop per method), Table 3 (leakage decomposition + permutation p-value).
- [ ] Hetionet robustness table (same audit, second graph).
- [ ] Pre-registered `PROTOCOL.md` timestamped; negative controls pass; no fabricated numbers.
- [ ] Zenodo DOIs; one-command repro verified on a fresh clone by a machine that didn't build it.
- [ ] Senior reader/co-author engaged; adversarial review done.
- [ ] bioRxiv preprint posted; resource-tier venue submission formatted and sent.

*First concrete step: Machine C builds the R0–R3 regime generator (BUILD_GUIDE Unit 2) and extends
`lib_eval`'s regime map, while B prototypes the degree-preserving swap on a 100k-edge sample and A
preps the filtered subgraph. Nothing is comparable until the regime files exist and are hashed.*
