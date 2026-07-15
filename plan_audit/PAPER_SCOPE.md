# PAPER_SCOPE.md (audit plan) — locked scope for the leakage-audit + de-leaked benchmark

*Created 2026-07-10. This is the contract for the **pivoted** paper. If a task, figure, or
paragraph does not serve one of the contributions C1–C4 below, it is out of scope. Re-read before
adding anything. Supersedes `../PAPER_SCOPE.md` (the orthology-only framing) if the pivot is adopted.*

> **RECONCILED 2026-07-15 to the completed R0–R3 results.** The benchmark is finished
> (TransE + RotatE + 5 topological baselines, 3 seeds, 4,228 held-out human gene–disease edges).
> Two things this file said before the data came in are now **corrected in place** below, and the
> pre-registered framing is updated (this is a scientific-conclusion change — flagged, not silent):
>
> 1. **KGE models are TransE + RotatE, not TransE + ComplEx.** ComplEx is untrainable at this scale
>    (451k entities): every tractable negative-sampling recipe collapses to a trivial optimum, and
>    the only loss that avoids it is compute-intractable. RotatE was *not* cut. ComplEx's degeneration
>    is reported as an honest scale-limitation finding (pre-empts "why only one KGE family?").
> 2. **The headline is no longer "de-leaking drops *every* method."** Of the three candidate leakage
>    sources, **only degree (R2) is material.** Redundancy (R1) and cross-species orthology (R3) came
>    back **~null for every method** — so they are now *designed negative controls*, not load-bearing
>    claims (a strength: we tested the orthology hypothesis a biomedical reviewer would raise, and
>    refuted it). Measured drops, sampled-50neg MRR:
>
>    | R0→ (ΔMRR) | CommonNeigh | AdamicAdar | Jaccard | PrefAttach | TransE | RotatE | role |
>    |---|---|---|---|---|---|---|---|
>    | **R1 redundancy** | +0.3% | 0.0% | +0.3% | −0.8% | ~0% | ~0% | null control |
>    | **R2 degree null** | **−44%** | **−43%** | **−49%** | **+0.5%** | **−24%** | **−21%** | the effect |
>    | **R3 orthology** | +0.5% | +0.3% | +0.6% | −1.0% | +0.2% | +0.4% | null control |
>
>    PreferentialAttachment (pure degree) is unchanged under the degree null by construction — it
>    validates that R2 preserves the degree sequence exactly. See `manuscript/` and the memory notes
>    `benchmark-results-degree-driven` / `publication-strategy` for the locked framing.
>
> A **full-ranking robustness check** (rank the true disease against all ~14.3k diseases, not 50
> sampled negatives) on R0+R2 is in `data/processed/results/robustness/` + `ROBUSTNESS.md`; it
> hardens the centerpiece against the "sampled-50neg is optimistic" objection.

**Working title:** *Quantifying evaluation leakage in multi-species gene–disease link prediction:
a graded audit and a de-leaked benchmark.*

**Genre (the strategic decision):** a **benchmark + audit + resource** paper. Soundness-based venues
publish this genre on *correctness and usefulness*, not on a novel positive result — so acceptance
does not hinge on any single number coming out a particular way. This is deliberate: it removes the
single point of failure that the orthology-only draft carried.

---

## 0. Why this framing (and not the two alternatives)

**Not a plain benchmark.** "Adamic-Adar beats KGE and AUROC is misleading" is already known
(Crichton 2018; the AUROC-under-imbalance literature) and would be rejected as not novel.

**Not the orthology-only bet.** The earlier draft staked the entire paper on one make-or-break
number: does the orthology effect survive a degree control? If it doesn't, that draft collapses into
a restatement of known degree bias (Gu et al. 2024 already owns degree-bias-on-Monarch). Our own
strategy note rated a real-journal "novel finding" at ~25–40% on that bet.

**This plan** widens the scope into a full leakage *audit* and reframes the headline so it is robust
to (a) the orthology result and (b) reviewer complaints about KGE tuning:

> **Thesis (reconciled to results; still tuning-proof).** Under a graded leakage audit of a
> multi-species gene–disease benchmark, **degree structure is the one material source of evaluation
> leakage**: a degree-preserving null (R2) drops the overlap heuristics 43–49% and both KGE models
> 21–24% MRR, while redundancy (R1) and cross-species orthology (R3) leave every method unchanged
> (±1%). We quantify each candidate source, release the de-leaked splits, and show the degree effect
> is not a KGE-tuning artifact — it hits simple baselines and embeddings alike.

Because the *material* effect (degree) drops embeddings **and** simple baselines alike, the claim does
**not** depend on embeddings being well-tuned or on baselines "winning." The orthology and redundancy
regimes ride along as **negative controls** — both came back null, which is itself a reportable,
reviewer-anticipating result, not a load-bearing claim.

**Honesty note (state this in the paper).** Novelty here is *incremental and methodological*, not a
landmark discovery. The ortholog-annotation-bias phenomenon is known (Alghamdi/Hoehndorf/Robinson
2022); our delta is packaging it as a formal, *quantified* leakage source inside a KGE benchmark,
plus a released, reusable de-leaked resource. Do not oversell it.

---

## 1. The core content (what the paper actually establishes)

The prediction target is a real, human task, confirmed from
`data/processed/edges_clean_integrated.csv` (2026-07-07):

- Target relation `GENE_ASSOCIATED_WITH_CONDITION` = **42,288 edges, 100% human** (all source nodes
  are HGNC). The task is *human* gene → MONDO disease.
- The graph carries **510,703 `ORTHOLOGOUS_TO` edges (8.7% of all edges)** bridging human genes to
  mouse/fish/fly/worm orthologs, which hold the bulk of **1,066,114 `HAS_PHENOTYPE`** edges, plus
  **9,495 `MODEL_OF`** (disease ↔ model-organism) edges.

**A graded leakage taxonomy** — the audit measures how much of reported performance is each of:

| Source | What it is | Regime that removes it |
|---|---|---|
| **Redundancy** | inverse / duplicate / symmetric edges that restate a held-out pair | R1 |
| **Degree** | performance attributable to node degree alone (illegitimate feature) | R2 (permutation null) |
| **Orthology** | cross-species ortholog → shared-phenotype shortcut to the target | R3 (Monarch only) |

**The orthology leakage mechanism (the secondary result).** A held-out human edge
`HGNC:gene —GENE_ASSOCIATED_WITH_CONDITION→ MONDO:disease` can be recovered *without independent
human evidence* by traversing:

```
HGNC:gene --ORTHOLOGOUS_TO--> mouse/fish ortholog --HAS_PHENOTYPE--> phenotype P
MONDO:disease --HAS_PHENOTYPE--> phenotype P     (and/or  MONDO:disease --MODEL_OF--> ortholog)
```

Neighborhood methods and KGE can score the true partner highly purely through this shortcut — the
model looks like it is *discovering* human disease genes when it is *re-deriving* them from mouse.

**Why it matters (the "who cares").** A large fraction of computational gene–disease and
drug-repurposing work is benchmarked on aggregated, ortholog-rich graphs (Monarch KG, Hetionet-style,
PrimeKG-style). If headline performance is substantially leakage, the field is systematically
over-reporting its ability to find *novel human* associations. A rigorous audit **plus a corrected,
released protocol** is a genuine, useful, publishable contribution — the kind soundness venues want.

---

## 2. Working titles (pick one at draft time)

1. **"Quantifying evaluation leakage in multi-species gene–disease link prediction: a graded audit
   and a de-leaked benchmark."** *(recommended — names the deliverable)*
2. "De-leaking drops everyone: a graded leakage audit of gene–disease link-prediction benchmarks."
3. "Beyond AUROC: redundancy, degree, and orthology shortcuts in gene–disease link prediction, and a
   released de-leaked evaluation protocol."

---

## 3. Central hypotheses and the killer experiment

**H1 (primary, tuning-proof) — supported.** Of the graded leakage sources, the **degree-preserving
null (R2)** is the one that reduces measured ranking performance — for the overlap heuristics
(−43 to −49% MRR) and for both KGE models (−21 to −24%) alike; redundancy (R1) and orthology (R3) do
not (±1%). The per-source, per-method drop is the headline, together with whether the **method
ranking changes** under the degree null (PreferentialAttachment, pure degree, is unchanged by
construction and validates the null).

**H2 (secondary, Monarch-specific) — a clean null.** Does cross-species orthology contribute leakage
*beyond* degree? **No.** Removing all ~520k `ORTHOLOGOUS_TO` + `MODEL_OF` edges (R0→R3, and the
R2→R3 residual) leaves every method within ±1% (KGE +0.2–0.4%). Reported as the honest one-liner the
plan anticipated: orthology leakage is not a distinct source on this graph. This *refutes* the
domain-specific hypothesis a biomedical reviewer would raise — a strength, not a failed paper.

**Killer experiment (the figure that sells the paper):** every method's ranking (MRR, Hits@k) under
four regimes on the *same* held-out human gene–disease edges:

| Regime | What it removes | Tests for |
|---|---|---|
| **R0 Standard** | nothing (random split; current practice) | the inflated baseline everyone reports |
| **R1 Redundancy-controlled** | inverse / duplicate / symmetric edges leaking the test edge | trivial restatement leakage |
| **R2 Degree-controlled** | topology, via a degree-preserving permutation null | performance attributable to degree alone |
| **R3 Orthology-blocked** (Monarch only) | `ORTHOLOGOUS_TO` + `MODEL_OF` edges, measured *after* R2 | orthology leakage *distinct from* degree |

**Prediction (as it resolved):** performance falls sharply at **R2 only** (the degree null); R1 and
R3 are flat. The orthology residual (R2→R3) came back **small/null**, which is reported as an honest,
still-publishable finding (orthology leakage ≈ degree bias here). **Stretch regime (only if on
schedule):** a temporal split on two dated Monarch releases (train older, test newer-only edges).

*Honesty guardrail:* we report whatever the drops actually are. The paper's value is a *rigorous*
audit and a *reusable* resource, not a big number.

---

## 4. Contributions (the four things the paper delivers)

- **C1 — A graded leakage-audit protocol + open harness, and the released de-leaked splits.** A
  single scoring path (`lib_eval`) and reusable, SHA-256'd R0–R3 splits any multi-species
  gene–disease benchmark can adopt. **This released resource is the centerpiece.**
- **C2 — The measured, per-source leakage contribution across methods.** How much of each method's
  reported performance is redundancy vs degree vs orthology, under one strict protocol, with ≥3 seeds,
  bootstrap 95% CIs, paired-bootstrap comparisons, and a permutation p-value for the degree control.
  Includes the AUROC-vs-ranking dissociation as a *sub-result*, now explained mechanistically.
- **C3 — The degree-vs-orthology decomposition (secondary, Monarch-specific): a reported null.** The
  R2→R3 residual is ~0 — cross-species orthology is **not** a leakage source distinct from degree on
  this graph. Stated as one honest sentence; the audit + released resource stand regardless.
- **C4 — Cross-graph robustness (generality).** The same audit rerun, unmodified, on a second
  pre-built graph (**Hetionet**) — showing "de-leaking drops methods" is not an artifact of one graph.
  (Orthology R3 is Monarch-only; stated explicitly.)

**What the paper explicitly does NOT claim** (put these in Related Work / Limitations up front):
no novel dataset (the graph is Monarch-derived); no invented predictor (the hybrid is cut); no
discovery of ortholog bias (known since 2022) — only its formalization, quantification, and a
released de-leaked benchmark.

**One-sentence contribution statement for the intro:** *"We audit multi-species gene–disease
link-prediction benchmarks under a graded protocol that removes redundancy, degree, and cross-species
orthology leakage in turn, show that the **degree** component is the one that materially inflates every
method (redundancy and orthology are null), quantify each leakage source with significance tests,
verify the degree effect on a second graph, and release the de-leaked splits and scoring harness as a
reusable resource."*

---

## 5. How we defend the work against "this is already known / invalid"

| Anticipated reviewer objection | Our defense (must be in the paper) |
|---|---|
| "Your KGE was undertrained, so the comparison is invalid" | Headline is "the **degree null** drops *every* method — KGE (−21 to −24%) included," not "baselines beat KGE." We use cited standard hyperparameters + a small R0 sweep; the degree claim is invariant to KGE quality. |
| "Topology beats KGE is old news" (Crichton 2018) | Agreed and cited as a sub-result. Note Crichton found AA the *worst* baseline — our AA-strength is graph-specific, which is itself part of the audit. |
| "AUROC-is-misleading is known" | Cited; we go further — we show the *ranking* signal itself is partly leakage, and decompose it. |
| "Degree bias on Monarch is known" (Gu 2024) | We *reuse* their degree-permutation idea as our R2 control and their filtering as a compute-saver, and cite them explicitly. Our delta is the graded multi-source audit + released splits + the orthology regime they omit. |
| "Data leakage in KGE is already taxonomized" (Ranga 2025) | Cited. They cover redundancy/degree/cold-start on *drug* repurposing; we add a cross-species *orthology* regime for *gene–disease* and release cleaned splits rather than only measuring. |
| "Ortholog bias is already known" (Alghamdi/Hoehndorf/Robinson 2022) | Cited and explicitly distinguished — closest prior art on the *phenomenon*. Our contribution is the KGE-benchmark *leakage* framing, per-regime quantification, and a released de-leaked split, not the discovery. |
| "Only one KG" | Zero-integration robustness run on Hetionet (C4). |
| "Orthology is just degree bias" | Reported honestly via the R2→R3 residual; a null is one sentence, and the paper (audit + resource) still stands. |

---

## 6. In scope / out of scope (freeze this)

**In scope (build these):**
- The frozen Monarch-derived graph as substrate; **Hetionet** as a separate robustness graph
  (never merged).
- Human HGNC gene → MONDO disease as the single prediction task.
- Methods: **Random, Common-Neighbors, Adamic-Adar, Jaccard, Preferential-Attachment** (baselines,
  the primary evidence); **TransE + RotatE** (KGE, two families, dim 64, cited hyperparameters +
  one small R0 sweep, ≥3 seeds, on a filtered subgraph). *(ComplEx was attempted but is untrainable
  at 451k entities — reported as a scale-limitation finding, not a fourth model.)*
- Regimes R0–R3 (C1); measured per-source impact + significance (C2); degree-vs-orthology
  decomposition (C3); Hetionet robustness (C4).
- Released de-leaked splits + one-command reproduction + Zenodo DOIs.

**Out of scope (archive; do NOT let these into the paper):**
- The **hybrid predictor** (old C4) — cut; the paper claims no invented method.
- **Leave-one-source-out** and **newer-release recovery** validation — cut (near-circular; a later
  release may itself propagate via orthology).
- **Heavy biological validation** — reduced to a light literature spot-check of a few top predictions.
- ComplEx / DistMult (KGE is **TransE + RotatE**; ComplEx untrainable at this scale, DistMult not run).
- Any "we integrated 13 databases / novel resource" claim; temporal growth on fabricated years.
- Neo4j drug-repurposing, the Streamlit portal, the HTML explorers, the old 606-node consensus subset.

---

## 7. Frozen artifacts (single source of truth — hash these)

- **Graphs:** `data/processed/edges_clean_integrated.csv` (Monarch; sha256 in `ARTIFACT_HASHES.txt`)
  and `data/external/hetionet/` (Hetionet; hash it). 5,860,539 Monarch edges; node count per
  `graph_stats.json` (authoritative).
- **Regimes (the resource):** `data/processed/splits/` — `train.csv`, `valid.csv`, `test.csv`,
  `train_R1_redundancy.csv`, `train_R2_degree_null_seed*.csv`, `train_R3_orthology_blocked.csv`,
  `split_manifest.json`; every file SHA-256'd.
- **Benchmark:** `results/baselines/*.json`, `results/kge/*.json`, merged `benchmark_results.json`
  (per-seed) → Table 2.
- **Null / decomposition:** `results/null/` → Table 3 + Figure 3.
- **Robustness:** `results/hetionet/*.json` → Table 4.
- **Protocol:** a timestamped `PROTOCOL.md` (pre-registration) committed before R3 is examined.

---

## 8. Success criteria (what "done and publishable" means)

1. R0→R3 leakage drop measured, with bootstrap 95% CIs over ≥3 seeds, for every method (5 baselines
   + 2 KGE), on Monarch.
2. The degree control (R2) reported with a permutation p-value; the R2→R3 orthology residual reported
   honestly (large or null).
3. At least one pairwise method difference tested for significance (paired bootstrap on per-edge
   reciprocal ranks) — not eyeballed.
4. Negative controls pass: Random ≈ chance (AUROC 0.50, MRR ≈ 1/((n_neg+1))); a shuffled-label run
   destroys performance.
5. The Hetionet audit reproduces the qualitative R0→R2 drop on a second graph.
6. Every number regenerates from frozen, hashed inputs via one command on a clean clone; de-leaked
   splits released with Zenodo DOIs.
7. A domain expert (mentor/reader) has reviewed the framing and the biology.

---

## 9. Target venues (resource / soundness tier)

- **Resource-led (best odds if the released benchmark is the centerpiece):** GigaScience DB,
  *Database (Oxford)*.
- **Methods-led:** PLOS ONE, PeerJ, NAR Genomics & Bioinformatics, Bioinformatics Advances.
- **Guaranteed parallel track:** Journal of Emerging Investigators (mentored, honest, real).
- **Always:** bioRxiv preprint the day the draft is done. Avoid pay-to-publish/predatory venues.

---

## 10. Execution plan

Operations (roles, job manifest, coordination, excellence practices) live in
**`WORK_SPLIT.md`**; the unit-by-unit construction manual is **`BUILD_GUIDE.md`**; environment setup
is **`SETUP_ENVIRONMENTS.md`** — all in this `plan_audit/` folder. `README.md` holds the reuse/scrap
map from the old build.

*Next action after this file: Machine C reworks the R0–R3 regime generator
(`scripts/build_deleaked_splits.py`) and the `lib_eval` regime map — the backbone every machine
depends on — while B prototypes the degree-preserving permutation null and A preps the filtered
subgraph. Nothing is comparable until the regime files exist and are hashed.*
