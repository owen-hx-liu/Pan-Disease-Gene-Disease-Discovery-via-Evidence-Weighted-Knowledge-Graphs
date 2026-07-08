# From Science-Fair Project to Peer-Reviewed Paper — The Complete Roadmap

*Prepared 2026-07-07. This is the single document to work from. It contains: an honest
assessment of where the project stands, the exact reframing that makes it publishable,
target journals, and a phase-by-phase, copy-paste-ready execution plan. Read Sections
1–4 once, then live in Sections 6–8.*

---

## 0. TL;DR (read this even if you read nothing else)

**Where you are:** You have a genuinely strong *engineering* asset (a clean 5.86M-edge
biomedical knowledge graph + a reproducible link-prediction benchmark harness) and a
genuinely *interesting scientific result* (simple graph topology beats fashionable deep
KG-embedding models at gene–disease ranking). What you do **not** have yet is a paper, a
defensible novelty claim, or prospective validation. Much of the "rigor" work the old
reports list as TODO is actually **already done** (held-out split, baselines, null
models, secrets removed, consensus rebuilt).

**The one sentence that blocks you today:** the graph is **98.4% Monarch** (5,765,251 of
5,860,539 edges). So "we integrated 13 databases into a novel knowledge graph" is *false*
and would get you desk-rejected the moment a reviewer runs one `grep`. You cannot publish
the *resource* as novel. **But you can publish the *benchmark and the finding*.**

**Publishability, honest odds:**

| Framing | Odds as-is | Odds after this roadmap |
|---|---|---|
| "Novel 13-database KG resource" | ~2% (false claim) | still weak — don't |
| "Rigorous benchmark: topology beats KGE for gene–disease ranking" | ~15% | **60–75%** at the right venue |
| Student/high-school journal (JEI, JSR, NHSJS) | ~40% | **85–95%** |

**The move:** stop selling the graph, start selling **the honest benchmark + the
negative/cautionary result + a focused rare-disease validation**. That reframe is the
whole game. Everything below executes it.

---

## 1. Full audit — what's actually in the folder

I went through every report, the scripts, the data folders, and the outputs. Here is the
real state, separating *done* from *claimed* from *missing*.

### 1a. Assets that are real and strong ✅
- **Integrated graph** `data/processed/edges_clean_integrated.csv` — 448,931 nodes /
  5,860,539 edges. Entity resolution validated (`step_4_5_validation_report.txt` = PASS,
  after a documented FATAL→PASS iteration — that iteration is a *good* methods-section story).
- **A real, executed link-prediction benchmark** (`scripts/kge_benchmark.py`) with a
  proper held-out 80/10/10 split, filtered all-entity ranking, and shared AUROC/AUPRC vs
  **both** random and type-matched hard negatives. Baselines (Random, Common-Neighbors,
  Adamic-Adar) + KGE (TransE/DistMult/ComplEx/RotatE). *This is the crown jewel.*
- **The headline finding is genuine and interesting:** Adamic-Adar reaches MRR **0.69** /
  Hits@10 **0.83**; the best KGE (RotatE) reaches MRR **0.09** / Hits@10 **0.17** — an
  order of magnitude worse at *ranking*, yet KGE *wins* the easy AUROC classification task
  (0.88–0.99). That dichotomy is a real, teachable, publishable result.
- **Full-graph network analysis** with a **degree-preserving null model**
  (`null_model_tests.py`): clustering z ≈ +230, assortativity z ≈ −154. Statistically
  grounded, not hand-waved.
- **Consensus confidence rebuilt** as independence-aware noisy-OR, sensitivity-tested
  (Spearman ρ = 1.00 under ±0.10 prior perturbation).
- **Reproducibility largely fixed:** `scripts/config.py` (no absolute paths), no committed
  secrets, `requirements.txt` pinned, fixed seed 42, `README.md`.
- **Face-validity biological check:** top novel Adamic-Adar predictions form a clean MODY
  cluster (GCK/HNF1A/HNF4A = MODY2/3/1) plus other confirmable gene–disease pairs.

### 1b. Claims in old docs that are now outdated (don't act on these) ⚠️
- CLAUDE.md still says "no train/test split," "run on 606-node subset," "remove secrets."
  **All done.** The `PROJECT_REPORT.md` addendum and `CRITICAL_FINDINGS.md` supersede it.
- The "13 integrated databases" language throughout CLAUDE.md / README is the *false claim*
  — see 1d.

### 1c. Real, current blockers ❌ (this is the roadmap's target list)
1. **98.4% Monarch.** Provenance table confirms it: non-Monarch sources add only ~95K of
   5.86M edges (1.6%), and 4 of 12 (DisGeNET, ClinGen, LNCipedia, LncRNADisease) add
   **zero**. The resource is not a novel integration.
2. **No novelty of *method*.** The winning predictor is Adamic-Adar (2003). Nothing you
   *invented* is the best thing in the paper. A paper needs a "what's new" sentence.
3. **No genuine temporal / prospective validation.** `temporalchart.py` *fabricates*
   discovery years with `random.choices(...)` (lines 44–50). There are **no real edge
   dates**, so the "time-based split (train pre-2020)" that CLAUDE.md recommends is
   currently *impossible*. Validation so far is retrospective face-validity, some of it
   near-tautological (e.g. RPGR → "RPGR-related retinopathy").
4. **KGE numbers not fully re-saved on the integrated graph.** The integrated run trained
   TransE (`benchmark_integrated_run.txt`) but `benchmark_results.json` only stored the
   baselines; the KGE table in `LINK_PREDICTION_BENCHMARK.md` mixes runs. A reviewer asking
   "reproduce Table 1" would hit a gap.
5. **No manuscript exists.** There is not a single `.tex`/draft in the repo. The science is
   further along than the *writing*, which is at zero.
6. **Consensus score is near-inert:** 99.15% of edges are single-source, so it only
   discriminates ~4,851 edges. It's honest but it's a footnote, not a contribution.

### 1d. The honesty ledger (must appear verbatim-ish in the paper)
These are non-negotiable disclosures. Hiding any one of them is how good projects become
retracted projects. Keep this list; every item maps to a sentence in the Limitations section.
- Graph is ~98.4% Monarch-derived; it is *not* an equal integration of 13 sources.
- Monarch is itself an aggregator; "supported by N sources" ≠ N independent votes.
- Edge "years" are simulated; no time-based split was possible; the split is random.
- Best predictions come from Adamic-Adar, not from any KGE model.
- Some validated hits are eponymous/tautological (disease named after the gene).
- Betweenness is sampled/approximate; some MONDO nodes are obsolete.
- KGE models are single-run and lightly tuned; DistMult collapsed (a tuning failure).

---

## 2. The three hard truths (why a reviewer would currently reject)

**Truth 1 — "Integration" is your headline and it isn't true.** Reviewers of resource
papers *always* audit provenance. `cut -d, -f5 edges_clean.csv | grep -iv monarch` returns
empty for the original graph; the integrated one is still 98.4% Monarch. If the paper's
first contribution is "a novel multi-database KG," it dies here. **Fix by deletion:** stop
claiming integration as the contribution.

**Truth 2 — Your best model is old and off-the-shelf.** That is *fine* — but only if the
paper is *about that*. If the paper pretends KGE is the contribution, a reviewer sees MRR
0.09 and rejects. If the paper's thesis is "on hub-dominated biomedical graphs, KGE's
famous AUROC hides catastrophic ranking failure, and a 20-year-old heuristic wins," then
the old model is the *point*, not a weakness. **Fix by reframing:** make the comparison the
contribution.

**Truth 3 — Retrospective face-validity isn't validation.** "GCK causes MODY, and our model
ranked it high" proves nothing if GCK↔MODY was derivable from the training graph. Reviewers
want *prospective* evidence: predict, then confirm against data the model never saw. You
can approximate this without waiting a year — see Phase 4 (leave-one-source-out + a newer
release). **Fix by protocol:** a temporal-surrogate validation.

None of these are fatal. All three are converted from weaknesses to strengths by the
reframing in Section 3.

---

## 3. The vision — the paper you are actually going to write

> **Working title:** *"Simple topology beats knowledge-graph embeddings for gene–disease
> link prediction: a rigorous benchmark on a large multi-species biomedical graph."*
>
> **One-paragraph abstract (the North Star — everything serves this):**
> *Knowledge-graph embedding (KGE) models are widely used to predict gene–disease
> associations, and are usually reported with high AUROC. We assemble a 5.86M-edge
> multi-species biomedical knowledge graph (derived primarily from the Monarch Initiative
> and augmented with seven curated sources) and benchmark four popular KGE models
> (TransE, DistMult, ComplEx, RotatE) against classical topological heuristics
> (common-neighbors, Adamic–Adar) under a single, strict protocol: held-out edges,
> filtered all-entity ranking, and shared classification with both random and
> type-matched hard negatives. We find a sharp dissociation: KGE models achieve
> AUROC 0.88–0.99 yet rank the correct partner an order of magnitude worse than
> Adamic–Adar (MRR 0.09 vs 0.69; Hits@10 0.17 vs 0.83). Reporting only AUROC — as much
> of the literature does — would have hidden this. We further show the winning heuristic
> recovers held-out rare-disease gene associations under a leave-one-source-out protocol,
> and we release the graph, harness, and predictions for reproducibility.*

**Why this is publishable when the "resource paper" isn't:**
- The contribution is a **method-evaluation finding + a reproducible harness**, not a claim
  of novel data. Every hard truth above becomes supporting evidence.
- Negative/cautionary benchmark results are *wanted* by a specific set of venues (they fight
  publication bias).
- It's honest end-to-end, so it survives adversarial review.
- It's *finishable by one student* — no new large-scale data engineering required.

**Three framings, ranked (pick one in Phase 0):**
1. **★ Benchmark/cautionary-tale** (above) — recommended. Highest rigor-to-effort ratio.
2. **Rare-disease discovery** — "Adamic–Adar surfaces novel rare-disease gene candidates
   from a Monarch-derived KG; N validated against Orphadata release X." Narrower, needs more
   biology, but higher *impact* if the validation lands. Good as a *second* paper or as
   Section 4 of framing #1.
3. **Resource/Application Note** — "A reproducible, provenance-tracked, Monarch-derived KG +
   benchmark harness." Only viable at *Database (Oxford)* / *Bioinformatics* Application
   Note, and only if you're scrupulous that it's Monarch-derived, not novel. Weakest novelty.

**Recommendation:** Do **#1 as the paper**, and fold a scoped version of **#2** in as the
biological-validation section. That single paper is your target.

---

## 4. Where to send it (venue strategy, tiered by realism)

Match the venue to who you are (a pre-college / early-undergrad researcher) and to the
finding (an honest benchmark). Aim *slightly* high, keep a realistic fallback.

### Tier A — Real peer-reviewed venues that fit a benchmark/negative result
| Venue | Why it fits | Bar | Notes |
|---|---|---|---|
| **PLOS ONE** | Publishes methodologically-sound results regardless of "impact"; loves honest negatives | Soundness, not novelty | Realistic reach goal. ~$1.9k APC — check waivers. |
| **BMC Bioinformatics** | Benchmarks + methods; student-friendly reviewers | Solid eval + baselines | You already have the eval. Good primary target. |
| **Database (Oxford)** | Resource/benchmark home; provenance-focused | Reproducibility | Best if you lean resource+benchmark. |
| **GigaScience / GigaByte** | Rewards reproducibility + data release | Repro package | GigaByte is lighter-weight, faster. |
| **F1000Research** | Open post-publication review; publishes then reviews | Transparency | Fast, honest, good for a first paper; reviews are public. |
| **PeerJ** | Sound-science model like PLOS ONE | Soundness | Lower APC; strong for CS-bio. |

### Tier B — Workshops / preprints (do these *in parallel*, they're not exclusive)
- **bioRxiv / arXiv (cs.LG or q-bio.QM)** — post the preprint the day you finish the draft.
  Free, timestamps your priority, lets you cite it. **Do this regardless.**
- **ISMB / PSB / NeurIPS-GLB (Graph Learning on Biology) workshops** — a benchmark result is
  a perfect workshop paper; often a stepping stone to a journal version.

### Tier C — Student/high-school journals (high acceptance, real editorial process)
These are *legitimate* venues for your career stage and a near-certain win after this roadmap:
- **Journal of Emerging Investigators (JEI)** — Harvard-run, free, mentored review, made for
  exactly this. **Strongly recommended as a guaranteed publication + practice run.**
- **Journal of Student Research (JSR)**, **National High School Journal of Science (NHSJS)**,
  **STEM Fellowship Journal**.

**Avoid:** anything that emails you offering fast publication for a fee, isn't indexed, or
isn't on DOAJ/PubMed. Check every candidate against **DOAJ** and **Think-Check-Submit**.
Predatory journals will *destroy* the credibility you're building.

**Recommended plan:** Submit to **JEI** (near-certain, mentored, teaches you the process)
**while** preparing the **BMC Bioinformatics / PLOS ONE** version. Preprint on bioRxiv the
moment the draft is done.

---

## 5. How to use the rest of this document

Sections 6–8 are the execution plan:
- **Section 6** = the phased roadmap. Each phase has: *Goal · Why · Steps · Done-when ·
  Est. time*. Do them in order; phases 1–5 are the science, 6–7 are the paper.
- **Section 7** = the **copy-paste prompt library** — the "paste this into another chat and
  it gets done" pieces you asked for. Each is self-contained.
- **Section 8** = timeline, risk register, and the final submission checklist.

A realistic total: **~8–14 focused weeks** part-time for the Tier-A version; the JEI version
can be draft-ready in **~3–4 weeks**.

---

## 6. THE ROADMAP (phase by phase)

> Legend: 🎯 Goal · 🧠 Why it matters to a reviewer · 🔧 Steps · ✅ Done-when · ⏱ Estimate

### PHASE 0 — Lock the framing and freeze scope *(do first, 1 day)*

🎯 Commit to Framing #1 (benchmark + folded-in rare-disease validation) and write it down so
you never drift back to "we integrated 13 databases."

🧠 90% of unpublishable student papers fail because scope kept expanding. A one-page scope
lock is the cheapest insurance you will ever buy.

🔧 Steps:
1. Create `PAPER_SCOPE.md` with: the title, the one-paragraph abstract from Section 3, the
   *three* contributions (C1 harness, C2 the AA-beats-KGE finding, C3 leave-one-source-out
   rare-disease validation), and an explicit **"out of scope"** list (novel integration
   claim, Neo4j drug-repurposing, Streamlit UI, temporal analysis with fake years).
2. Pick the primary venue (recommend BMC Bioinformatics) + the guaranteed venue (JEI) +
   preprint server (bioRxiv). Note their length/format limits in `PAPER_SCOPE.md`.
3. Delete-or-quarantine the parts you're not using so they can't leak into the paper:
   move `temporalchart.py`, `drugrepurposing.py`, `modelinterface.py`, the `.html`
   explorers into an `archive/` folder (keep in git history, out of the paper's story).

✅ Done-when: `PAPER_SCOPE.md` exists, venues chosen, distractions archived.
⏱ 1 day.

---

### PHASE 1 — Finalize the graph & provenance (make the data claims bulletproof) *(2–4 days)*

🎯 Every number about the graph in the paper must be regenerable by one command, and the
provenance must say "Monarch-derived, augmented" everywhere.

🧠 Resource/benchmark reviewers audit provenance first. If your abstract says 5.86M edges
and a script says 5.81M, you look careless; if it says "integrated 13 databases," you look
dishonest. Kill both risks now.

🔧 Steps:
1. **Freeze one canonical graph file** and one node-count/edge-count. Use
   `data/processed/edges_clean_integrated.csv` (448,931 nodes / 5,860,539 edges). Write a
   tiny `scripts/graph_stats.py` that prints nodes, edges, unique undirected edges,
   self-loops, per-namespace node counts, per-relation edge counts, and saves them to
   `data/processed/graph_stats.json`. This single JSON becomes the source of truth for
   Table 1.
2. **Rewrite every "13 databases integrated" string** in README/CLAUDE/PROJECT_REPORT to
   "a Monarch-derived multi-species KG augmented with 7 curated sources (8 ingested, 4
   documented as non-integrable)." Grep for it: search `integrate`, `13`, `thirteen`.
3. **Finalize the provenance table** (`source_provenance.md` is 90% there). Add the two
   things journals now require and you're missing: (a) **exact upstream release version +
   download date per source** — go to each source's site, record the version string; put
   it in the table (the doc flags `source_release_version` as a placeholder — fill it).
   (b) **license per source** (Monarch = BSD-3/CC0 varies; check each) so your data-release
   is legally clean.
4. **Decide on the 4 zero-edge sources.** Either (a) get them contributing (DisGeNET needs a
   different dump with a disease column; ClinGen needs the gene-disease validity file, not
   the node table) — *only if easy*; or (b) **drop them from the count entirely** and say
   "8 sources ingested" not 12. Don't list sources that contributed nothing.

✅ Done-when: `graph_stats.json` regenerates Table 1; provenance table has versions + dates
+ licenses; no "13 databases" claim survives anywhere.
⏱ 2–4 days (mostly the manual version/license lookup).

---

### PHASE 2 — Make the benchmark reproducible and complete *(3–6 days, needs the GPU box)*

🎯 One command reproduces the **entire** Table 1 (baselines *and* all four KGE models) on
the frozen integrated graph, and the numbers are saved to disk with seeds.

🧠 The current `benchmark_results.json` only has baselines; the KGE numbers live in a
markdown table from a partly-different run. A reviewer (or the "Reproducibility" checkbox at
BMC) will ask you to regenerate Table 1. It must Just Work.

🔧 Steps:
1. **Re-run KGE on the frozen graph and save results to JSON**, not just to the log. Patch
   `kge_benchmark.py` so each KGE model's ranking + classification metrics are written into
   the same `benchmark_results.json` alongside the baselines. (Right now TransE trains but
   its metrics don't land in the integrated JSON.)
2. **Add ≥3 seeds and report mean ± sd** for every model. Single-run numbers are a
   guaranteed reviewer complaint. Loop seeds {42, 1, 7}; store all runs; report mean±sd in
   the table. This is the single highest-value rigor upgrade.
3. **Fix DistMult** (it collapsed to ~0 — a tuning failure, not a finding). Try: add
   regularization, use the same loss/negative-sampler as the others, or lower LR. If it
   still collapses after an honest attempt, report *that* as a documented tuning outcome.
   Either way it must not look like you fumbled a baseline.
4. **A small hyperparameter sweep** for the KGE models (dim ∈ {128, 256}, epochs ∈
   {100, 300}, negative sampling ∈ {basic, bernoulli}). You are not trying to make KGE win —
   you are pre-empting the reviewer who says "you just didn't tune it." Report the *best*
   KGE config you found; if AA still wins, the finding is bulletproof.
5. **Emit a `benchmark_table.md` and `benchmark_table.csv`** straight from the JSON so the
   paper's Table 1 is generated, never hand-typed.

✅ Done-when: `python scripts/kge_benchmark.py --edges <frozen> --out <dir> --seeds 42,1,7`
regenerates every cell of Table 1 (baselines + KGE, mean±sd) into JSON + a table file, and
the AA-beats-KGE-on-ranking conclusion holds under your best KGE config.
⏱ 3–6 days (GPU training time dominates; each KGE model ~15 min × models × seeds × sweep).

---

### PHASE 3 — Nail the headline analysis (the AUROC-vs-ranking dissociation) *(2–3 days)*

🎯 Turn the finding into a figure + a mechanism, not just a table row.

🧠 "AA > KGE" is a table. "AUROC hides a 7× ranking gap, *because* the graph is
hub-dominated and classification only needs to beat a random pair" is a *paper*. Reviewers
reward a mechanism.

🔧 Steps:
1. **Figure 1 — the dissociation plot.** Two panels: (left) AUROC (hard negatives) per
   method; (right) MRR per method. Same x-axis order. The visual of KGE tall on the left,
   tiny on the right, is your whole thesis in one image.
2. **Explain *why* mechanistically.** Compute and report: fraction of test edges touching a
   generic hub (GO:0005515 etc.); show that KGE ranking failure concentrates on high-degree
   tails; show classification is easy because random/type negatives are structurally far
   from true edges. Tie to the null-model result (clustering z≈+230 ⇒ neighborhood methods
   have signal to exploit; that's *why* AA wins).
3. **Stratify results** by relation type and by whether the disease is rare (MONDO under the
   Orphanet branch) vs common. A benchmark that reports where each method wins/loses is far
   stronger than one global number.
4. **Ablation: strip the ultra-generic hubs** (degree > 2,000) and re-rank. Show whether the
   gap narrows. This directly answers "is AA just exploiting hubs?" — and it's already
   half-done in `predict_adamic_adar.py` (hubs excluded as shared-neighbor evidence).

✅ Done-when: Figure 1 exists; you can state in one sentence *why* the dissociation happens,
backed by numbers; results are stratified (relation type, rare vs common, hub-stripped).
⏱ 2–3 days.

---

### PHASE 4 — Prospective-surrogate biological validation *(4–7 days — the credibility maker)*

🎯 Move from "these predictions look right" to "these predictions recover associations the
model never saw."

🧠 This is what separates a benchmark from a *useful* benchmark. It also converts Hard
Truth 3 into a strength. Two complementary protocols, neither needs you to wait a year:

🔧 Steps:
1. **Leave-one-source-out (LOSO) validation.** Retrain AA/KGE with one curated source fully
   *removed* (best candidate: Gene2Phenotype or Orphadata gene–disease edges), then measure
   precision@k / recall on exactly the edges that source uniquely contributed. This is a
   clean "predict data you didn't train on" test and is *feasible today* because you have
   sole-source counts in the provenance table (G2P sole-source = 2,139; Orphadata = 4,654).
   Report precision@{10,50,100} and compare AA vs best KGE — a second head-to-head that
   isn't circular.
2. **Newer-release temporal surrogate.** Download the *current* Monarch/Orphanet release,
   diff against your frozen version, and treat associations that appear *only* in the newer
   release as a prospective test set. Report how many your top-k predictions recovered.
   (This is the closest you can get to a real time-split given the fake-year problem — and
   you *say so* in the paper.)
3. **Literature confirmation of the top novel non-tautological predictions.** Take the top
   ~20 AA predictions, **remove eponymous/tautological ones** (disease named after the gene)
   and **obsolete MONDO terms**, then check each survivor against PubMed / OMIM / recent
   GWAS. Report a small precision@k with citations. Even 5–8 literature-confirmed novel
   candidates is a compelling result — and being upfront about excluding the tautological
   ones *builds* trust.
4. **One deep dive.** Pick the single most interesting validated novel prediction and write
   a paragraph: what the gene does, why the link is plausible, what evidence supports it.
   Reviewers remember the story, not the table.

✅ Done-when: LOSO precision@k table (AA vs KGE); newer-release recovery number; a curated,
de-tautologized top-k literature table with citations; one narrative deep-dive.
⏱ 4–7 days (the newer-release download + diff is the variable part).

---

### PHASE 5 — Reproducibility & public release (get your DOI-able artifact) *(2–3 days)*

🎯 A reviewer can clone your repo, run one command, and reproduce the tables/figures; the
data + code have permanent DOIs.

🧠 For BMC/PLOS/GigaScience/Database, the reproducibility package is *scored*. It's also the
easiest place to look professional, and it's mostly done.

🔧 Steps:
1. **Clean the public repo:** remove `venv/` (Windows env — already in `.gitignore`?),
   archived distractions out of the main path, ensure `scripts/config.py` relative paths
   only, confirm zero secrets (`git log -p | grep -i password` clean).
2. **One-command reproduction:** a `Makefile` or `run_all.sh` /
   `run_all.ps1` that runs provenance → null-model → benchmark → predictions → figures.
   Document expected runtime + hardware (RTX-class GPU for KGE; baselines run on CPU).
3. **Archive the data on Zenodo** (the 360 MB `edges_clean_integrated.csv` + provenance) →
   get a **DOI**. Zenodo integrates with GitHub to snapshot a release → get a **code DOI**.
   Cite both in the paper's Data/Code Availability section. *This is often the difference
   between "revise" and "accept" at reproducibility-focused venues.*
4. **Write the data card:** a short `DATA.md` — what the graph is, schema of the CSV
   (`:START_ID,:END_ID,:TYPE,weight,dataset_sources`), how it was built, license,
   Monarch-derived caveat.
5. **License the repo** (MIT or Apache-2.0 for code; match data license to the most
   restrictive source you include).

✅ Done-when: fresh clone + one command reproduces Table 1 and Figure 1; Zenodo DOIs minted
for data and code; `DATA.md` + LICENSE present.
⏱ 2–3 days.

---

### PHASE 6 — Write the manuscript *(2–3 weeks — the real work)*

🎯 A complete draft in the target journal's structure. Use the prompt library (Section 7)
to draft each section, then *you* edit for honesty and voice.

🧠 The science is done after Phase 5; now it's communication. Reviewers reject unclear
papers with good science more often than clear papers with modest science.

🔧 Structure (IMRaD, ~4,000–6,000 words for BMC Bioinformatics):
1. **Title + Abstract** — use Section 3's abstract as the seed; update with final numbers.
2. **Introduction** — the problem (gene–disease prediction matters), the gap (KGE reported
   with AUROC only; ranking under-examined; benchmarks often leak), your contributions
   (C1/C2/C3). End with an explicit "our contributions are…" bullet list.
3. **Related work** — position honestly vs Monarch KG, Hetionet, PrimeKG, OpenBioLink,
   BioKG. State plainly: *you use a Monarch-derived graph; your novelty is the evaluation
   and the finding, not the data.* This paragraph disarms the "not novel data" reviewer.
4. **Methods** — graph construction + provenance (Phase 1); the benchmark protocol (split,
   filtered ranking, hard negatives, seeds — Phase 2); baselines & KGE configs; the LOSO +
   newer-release validation (Phase 4); null-model tests. Include the FATAL→PASS entity-
   resolution story as a validation subsection — it shows rigor.
5. **Results** — Table 1 (benchmark, mean±sd), Figure 1 (dissociation), stratified results
   (Phase 3), LOSO/literature validation (Phase 4), network structure vs null (existing).
6. **Discussion** — the mechanism (why AUROC misleads, why AA wins), implications ("report
   ranking metrics; hard negatives; don't trust AUROC alone"), and **Limitations** = your
   entire honesty ledger (Section 1d), stated proactively.
7. **Conclusion**, **Data/Code Availability** (Zenodo DOIs), **Methods reproducibility**,
   **references** (use a manager; ~30–50 refs).

🔧 Process steps:
1. Draft section-by-section with the Section 7 prompts (paste real numbers in — never let
   the model invent numbers).
2. **Verify every number** against `graph_stats.json` / `benchmark_results.json`. One
   fabricated number found by a reviewer poisons the whole paper.
3. Make figures publication-quality (300 dpi, readable fonts, colorblind-safe palette).
4. Get **two humans** to read it: one bio-savvy, one not. If a non-expert can't state your
   finding after the abstract, rewrite the abstract.
5. Run it past a mentor/teacher/PI if you have one — a named senior co-author *dramatically*
   improves desk-reject odds at Tier-A venues (and is standard, not cheating).

✅ Done-when: a complete, numbers-verified draft in the journal's format, read by ≥2 people.
⏱ 2–3 weeks.

---

### PHASE 7 — Submit *(1 week + review cycles)*

🎯 Preprint up, manuscript submitted, cover letter written, reviewer responses handled.

🔧 Steps:
1. **Post to bioRxiv** the day the draft is final. Free, timestamps priority, citeable.
   Choose a CC-BY license.
2. **Format to the exact venue guidelines** (BMC has a template; PLOS has structural rules;
   JEI has its own). Length, section order, reference style, figure format — follow to the
   letter; editors reject on formatting.
3. **Write the cover letter** (Section 7 has a template): what the paper shows, why it fits
   *this* journal, why it matters, confirm it's original + not under review elsewhere,
   suggest 3–4 reviewers (and note any to exclude).
4. **Submit.** Then wait (weeks to months). 
5. **Revisions:** treat every reviewer comment as free improvement. Respond to *all* of
   them point-by-point in a response letter, even the ones you disagree with (politely
   explain). Most first decisions are "major revision," not accept — that's normal and good.
6. **If rejected:** it's data. Use the reviews, reformat, submit to the next venue down the
   list *the same week*. Papers get published by resubmission, not by being perfect first try.

✅ Done-when: submitted + preprint live.
⏱ 1 week to submit; 1–6 months of review cycles (venue-dependent).

---

## 7. Copy-paste prompt library ("enter this into another chat and get it done")

These are self-contained. Replace `<...>` placeholders. **Rule:** never let an AI *invent*
numbers — always paste your real JSON/CSV values in, and verify outputs against your files.

### 7.1 — Reproduce & complete the benchmark (Phase 2)
```
I have a link-prediction benchmark script at scripts/kge_benchmark.py that trains
TransE/DistMult/ComplEx/RotatE and compares them to Random/CommonNeighbors/AdamicAdar on a
biomedical knowledge graph (data/processed/edges_clean_integrated.csv, 448,931 nodes /
5,860,539 edges, columns :START_ID,:END_ID,:TYPE,weight,dataset_sources). Right now only the
baseline metrics get saved to benchmark_results.json; the KGE metrics only print to the log.

Please: (1) modify the script so every KGE model's ranking metrics (MRR, Hits@1/3/10) and
classification metrics (AUROC/AUPRC for both random and type-matched hard negatives) are
written into the same benchmark_results.json as the baselines; (2) add a --seeds argument
(default 42,1,7) that runs each KGE model once per seed and stores every run so I can report
mean±sd; (3) diagnose why DistMult collapses to ~0 MRR and try regularization / a matched
loss+negative-sampler to fix it, documenting what you changed; (4) emit benchmark_table.md
and benchmark_table.csv generated directly from the JSON. Keep the fixed-seed reproducibility
and the existing held-out filtered all-entity ranking protocol. Show me the diff before
running anything heavy.
```

### 7.2 — The dissociation figure + mechanism (Phase 3)
```
Using benchmark_results.json (fields: per-method MRR, Hits@1/3/10, AUROC_random,
AUROC_type, AUPRC_random, AUPRC_type), write a matplotlib script that produces a single
publication-quality figure (300 dpi, colorblind-safe, readable fonts) with two panels
sharing method order on the x-axis: left = AUROC (type-matched hard negatives), right = MRR.
The point is to visualize that KGE models are tall on AUROC but tiny on MRR while
Adamic-Adar is the reverse. Save as figures/fig1_dissociation.png and .pdf. Then write a
short analysis script that quantifies WHY: compute the fraction of held-out test edges that
touch a generic hub (degree > 2000), and report ranking metrics separately for
hub-touching vs non-hub-touching test edges, for AdamicAdar and the best KGE model. Output a
small markdown table I can drop into the Results section.
```

### 7.3 — Leave-one-source-out validation (Phase 4)
```
I need a leave-one-source-out (LOSO) validation. In data/processed/edges_clean_integrated.csv
the dataset_sources column tells which source(s) support each edge; data/processed/provenance/
source_provenance.md lists sole-source counts (e.g. Gene2Phenotype sole-source = 2,139,
Orphadata = 4,654). Write scripts/loso_validation.py that: (1) removes ALL edges whose sole
source is <SOURCE> (start with Gene2Phenotype) from the graph; (2) recomputes Adamic-Adar and
the best KGE model's gene->disease rankings on the reduced graph; (3) measures precision@k and
recall@k (k = 10, 50, 100) specifically on the held-out sole-source edges — i.e., did the
model predict associations it was never shown?; (4) writes a results table comparing AA vs
KGE. Filter to gene(HGNC)->disease(MONDO) pairs, exclude self-loops and generic hubs
(degree > 2000). Save outputs to data/processed/loso/. Explain the interpretation of the
numbers in plain language.
```

### 7.4 — Curated literature-validation table (Phase 4)
```
Here are my top 25 novel gene-disease predictions (Adamic-Adar) as gene HGNC IDs and disease
MONDO IDs: <paste top 25 rows of adamic_adar_gene_disease.csv, resolved to gene SYMBOL and
disease NAME if you can>. For each: (1) flag whether it is eponymous/tautological (the disease
is literally named after the gene) — these should be EXCLUDED from the headline precision
number; (2) flag obsolete MONDO terms; (3) for the survivors, tell me whether the gene-disease
association is supported in the literature (OMIM / recent GWAS / review articles), and give me
the specific evidence to check. Produce a markdown table with columns: rank, gene, disease,
tautological?, obsolete?, literature-supported?, evidence-to-verify. Be conservative — I will
manually verify every "supported" claim against PubMed before it goes in the paper; do not
assert support you can't point to a source for.
```

### 7.5 — Draft a manuscript section (Phase 6, repeat per section)
```
You are helping me write the <Methods> section of a bioinformatics benchmark paper for
<BMC Bioinformatics>. The paper's thesis: on a large hub-dominated multi-species biomedical
knowledge graph, KGE models achieve high AUROC but rank the correct gene-disease partner an
order of magnitude worse than the classical Adamic-Adar heuristic (MRR 0.09 vs 0.69), so
AUROC-only reporting is misleading. Here are the facts to write from — use ONLY these numbers,
do not invent any: <paste the relevant facts: graph_stats.json, the benchmark protocol from
README/LINK_PREDICTION_BENCHMARK.md, seeds, split ratios, KGE configs>. Write in clear,
precise scientific English, past tense, no hype. Include enough detail for exact
reproduction. Flag anywhere you need a number I haven't given you with <<NEED: ...>> rather
than guessing.
```

### 7.6 — Cover letter (Phase 7)
```
Write a submission cover letter to the Editor-in-Chief of <BMC Bioinformatics> for my
manuscript titled "<title>". Key points to convey: (1) what the paper shows — a strict
benchmark revealing that popular KGE models are outperformed an order of magnitude by
Adamic-Adar on gene-disease ranking despite higher AUROC, plus a reproducible harness and a
leave-one-source-out biological validation; (2) why it fits this journal's scope; (3) why it
matters (evaluation practice + a reusable resource); (4) a statement that the work is
original, not under consideration elsewhere, and that a preprint is on bioRxiv at <DOI>;
(5) 3–4 suggested reviewers with rationale. Keep it under one page, professional, no hype.
```

### 7.7 — Adversarial pre-review (do before submitting)
```
Act as a skeptical Reviewer 2 for a bioinformatics journal. Here is my manuscript draft:
<paste draft>. Find every weakness a reviewer would raise: unsupported claims, missing
baselines or ablations, statistical issues, reproducibility gaps, over-claims, and anywhere
the "novelty" is thin. For each, tell me the exact objection and the minimum change that
would satisfy it. Be harsh — I would rather hear it from you than in a rejection. Pay special
attention to: the Monarch-dominance of the graph, the fabricated edge years, single-run vs
multi-seed results, and whether the biological validation is truly prospective.
```

---

## 8. Timeline, risks, and the final checklist

### 8.1 Realistic timeline (part-time, one student)
| Weeks | Phases | Milestone |
|---|---|---|
| 1 | 0, start 1 | Scope locked, framing chosen, provenance versions gathered |
| 2 | 1, 2 | Graph frozen; benchmark re-runs with seeds into JSON |
| 3 | 2, 3 | KGE sweep done; DistMult fixed; Figure 1 + mechanism |
| 4 | 4 | LOSO + newer-release + literature validation |
| 5 | 5 | Repo cleaned; Zenodo DOIs; one-command repro |
| 6–8 | 6 | Full manuscript draft; internal reviews; adversarial pre-review |
| 9 | 7 | bioRxiv preprint + submit to JEI (guaranteed track) |
| 9–10 | 7 | Submit Tier-A version (BMC/PLOS); then review cycles |

*JEI-only fast path:* Phases 0–1 lite + existing results + write-up ≈ **3–4 weeks** to a
submittable JEI manuscript. Do this in parallel — it's a near-certain publication and teaches
you the whole pipeline while the Tier-A version matures.

### 8.2 Risk register
| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| Reviewer attacks Monarch-dominance | High | High | Pre-empt in Related Work + Limitations; novelty = evaluation not data |
| "You didn't tune KGE" | High | Med | Phase 2 sweep + best-config reporting; report mean±sd |
| Validation seen as retrospective-only | Med | High | Phase 4 LOSO + newer-release = genuine held-out signal |
| Single-seed results | Med | Med | ≥3 seeds, mean±sd everywhere |
| Fabricated years noticed | Med | High | Delete temporal analysis from paper; disclose in Limitations |
| Predatory-journal trap | Med | High | DOAJ + Think-Check-Submit; stick to the Section 4 list |
| Scope creep / never finishing | High | High | `PAPER_SCOPE.md`, archive distractions, JEI parallel track |
| No senior co-author → desk reject | Med | Med | Recruit a teacher/PI mentor as co-author early |

### 8.3 Pre-submission checklist (tick every box)
- [ ] Abstract states the finding + that the graph is Monarch-derived (no "novel integration").
- [ ] Every number in the paper traces to `graph_stats.json` / `benchmark_results.json`.
- [ ] Table 1 has all methods, mean±sd over ≥3 seeds, regenerated by one command.
- [ ] Figure 1 (dissociation) is 300 dpi, colorblind-safe, readable at print size.
- [ ] Baselines include Random + Common-Neighbors + Adamic-Adar; KGE includes ≥3 models.
- [ ] Ranking is filtered, all-entity/sampled-neg, documented; hard negatives used.
- [ ] Biological validation is LOSO/newer-release (not just eponymous face-validity).
- [ ] Tautological + obsolete predictions excluded from headline precision, disclosed.
- [ ] Null-model z-scores reported for network claims.
- [ ] Limitations section contains the full honesty ledger (Section 1d).
- [ ] Zenodo DOIs for data + code; one-command reproduction verified on a fresh clone.
- [ ] No secrets, no absolute paths (`git log -p | grep -i password` clean).
- [ ] License files present; per-source licenses compatible with release.
- [ ] Venue formatting followed exactly; cover letter written; preprint posted.
- [ ] ≥2 humans read it; adversarial pre-review done; a mentor/co-author is on board.

---

## 9. The bottom line

You do **not** have a "13-database integration" paper — and chasing that would waste months
and still fail review. You **do** have the raw material for an honest, publishable
**benchmark-and-cautionary-tale** paper whose central result (simple topology beats
fashionable KGE at gene–disease ranking, and AUROC hides it) is real, reproducible, and
genuinely useful to the field. The engineering that felt like the achievement is actually the
*supporting infrastructure*; the *finding* is the achievement.

Reframe around that finding (Section 3), execute Phases 0–7, publish the JEI version as a
guaranteed win while you push the BMC/PLOS version, and be scrupulously honest in the
Limitations. That is the realistic, concrete path from where this is today to a peer-reviewed
publication — and it is very achievable.

*Next action: create `PAPER_SCOPE.md` (Phase 0) and gather the per-source release
versions/dates (Phase 1, step 3). Those two things unblock everything else.*
