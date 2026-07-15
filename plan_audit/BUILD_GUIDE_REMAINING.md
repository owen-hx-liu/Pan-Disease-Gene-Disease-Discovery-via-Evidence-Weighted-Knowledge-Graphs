# BUILD_GUIDE_REMAINING.md — remaining steps, in depth (details + prompts + timeline)

*Created 2026-07-12. Supersedes the forward half of `BUILD_GUIDE.md` (Units 7→13). Folds in the
interim results from Units 2/4/5, and for each remaining task gives: **why it matters**, the
**in-depth gotchas**, an **acceptance test**, and a **copy-pastable prompt** you can hand to an AI to
generate the script. The exact order (and what's parallel/anytime) is in §Timeline at the bottom.*

---

## Where we are (done + what it taught us)

**Done + verified:** Setup + Hetionet download/hash · Unit 2 (audit-regime splits: R1 −1,661 edges,
R3 = old orthology file) · Unit 4 (`lib_eval` regime map; Random control passes) · Unit 6 (R1,
folded into Unit 2) · **Unit 5 (5 baselines × R0/R1/R3 × 3 seeds).**

**The three findings that reshape the rest of the plan:**

1. **De-leaking barely moves the topological baselines** (mean MRR over 3 seeds):

   | Method | R0 | R1 | R3 | R0→R3 |
   |---|---|---|---|---|
   | Random | 0.091 | 0.091 | 0.091 | 0 |
   | CommonNeighbors | 0.441 | 0.442 | 0.444 | +0.002 |
   | AdamicAdar | 0.459 | 0.459 | 0.460 | +0.002 |
   | Jaccard | 0.397 | 0.398 | 0.399 | +0.002 |
   | PreferentialAttachment | **0.586** | 0.582 | 0.581 | −0.005 |

2. **Pure degree wins** (PreferentialAttachment) → degree is the dominant signal → **R2 is decisive.**
3. **Orthology (R3) is invisible to 2-hop baselines** (the leak is a ≥3-hop path) → the orthology
   claim lives in **KGE** (multi-hop), not the baselines.

**Refined thesis (what the data support):** a graded audit showing *which leakage source affects
which method class* — redundancy negligible; **degree affects everything**; **orthology affects only
multi-hop methods (KGE)**, with a mechanistic explanation — plus the released de-leaked benchmark.

**Operational reality:** everything so far ran on **Machine A alone** in minutes; igraph is
installed. B/C are optional (only to overlap KGE training or add a writer). Machine tags below are
for the multi-machine case; solo = do them all on A.

**The shared harness every scorer must use** (so numbers stay comparable):
```python
import lib_eval
reg = lib_eval.load_regime("R2", splits_dir)          # or R0/R1/R3
train, test, hubs = reg
pools = lib_eval._category_pools(train, hubs)          # build ONCE per graph, reuse
known = lib_eval._known_tails_by_head(train, test, test[:, 0])
ranks, pos, neg = lib_eval.rank_test_edges(            # score_fn(gene,disease)->float
    score_fn, train, test, hubs, reg.hub_filter,
    n_neg=50, seed=seed, return_scores=True, pools=pools, known=known)
rm = lib_eval.ranking_metrics(ranks)                   # MRR, Hits@1/3/10
cm = lib_eval.classification_metrics(pos, neg)         # AUROC, AUPRC
ci = lib_eval.bootstrap_ci(1.0 / ranks)                # 95% CI on MRR
p  = lib_eval.paired_bootstrap_pvalue(rr_a, rr_b)      # method-vs-method significance
```

---

## Leftover Stage-1 items

### Unit 1b — `graph_stats.py` on Hetionet  🟢 · *anytime*
**Why:** Table 1 needs the second-graph column.
**In depth:** the Hetionet edges file is `source⇥metaedge⇥target` (tab-separated), and every ID is
`Type::id` (e.g. `Gene::1234`, `Disease::DOID:0050156`). So (a) the column names differ from
Monarch's `source_id,relation,target_id`, and (b) `kg_categories.category_of` won't recognize the
`Type::` prefixes — the node "category" in Hetionet is simply the string before `::`.
**Accept:** ~47,031 nodes, 2,250,197 edges, ~24 metaedges; the Gene–Disease metaedge `DaG` ≈ 12,623.
**Prompt:**
```
Extend scripts/graph_stats.py (or write scripts/graph_stats_hetionet.py, pandas only) to also
handle Hetionet. Input data/external/hetionet/hetionet-v1.0-edges.sif, tab-separated with columns
source,metaedge,target where every ID is "Type::id". Add --sep '\t' and --columns
source,metaedge,target options, and derive each node's category as the substring before '::'
(NOT via kg_categories.category_of, which is Monarch-specific). Save
data/processed/graph_stats_hetionet.json with n_nodes, n_edges, per-metanode node counts,
per-metaedge edge counts, degree mean/median/max, and the DaG (Disease-associates-Gene) edge
count. Print a summary. Deterministic, no network libs.
```

### Unit 3 — provenance fill  🟢 · *anytime · manual (yours)*
**Why:** resource venues score provenance. **In depth:** for each ingested source open
`data/processed/provenance/source_provenance.md` and fill the real **release version + download
date + license** — these are facts to look up, not to generate. **Accept:** no "placeholder"
strings remain. *(No AI prompt — it's lookup. I can pre-format the table if you give me the values.)*

### Housekeeping  🟢 · *anytime*
- **Track result JSONs:** add `!data/processed/results/**/*.json` to `.gitignore` so the small
  result files version (they feed figures/tables). Then `git add -f` the existing baseline JSONs.
- **Reconcile AdamicAdar 0.459 vs old ~0.69:** read `benchmark_results.json`, find the protocol
  difference (likely easier/leakier negatives), record one sentence for the paper's Methods.

---

## Unit 7 — `build_degree_null.py` (R2)  🟠 B / runs on A  ← **NEXT & DECISIVE**

**Why:** the experiment the baseline story hinges on — decompose each method's R0 MRR into
**degree** vs **real structure**.

**In depth (the gotchas that make or break it):**
- **Type-preserving swap, NOT a global one.** A naive `igraph.Graph.rewire()` on the whole graph
  swaps edges blindly and will connect a gene to a gene where a gene–disease edge used to be —
  nonsense. Rewire **within each relation type separately** (Hetionet-XSwap / Gu-et-al. style) so
  node degrees *and* the type/bipartite structure are preserved and the null isolates *pure degree*.
- **Preserve exactly the right thing.** After swapping, assert the per-node degree sequence is
  **identical** to the original and the per-relation edge counts are unchanged. (CLAUDE.md's rule:
  a degree null with a changed degree sequence is a bug.)
- **Hold the evaluation fixed for a clean permutation test.** Score the baselines on each null
  replicate **using the R0 test set and the R0 negatives** (pass R0's `pools`/`known` into
  `rank_test_edges`), changing *only the scoring adjacency*. Otherwise you're changing the graph and
  the negatives at once and can't attribute the effect. Because `pools`/`known` are reused and each
  null adjacency rebuild is ~5s, ~10–20 replicates is only a few minutes.
- **Sanity anchor:** PreferentialAttachment reads only degree, which the null preserves → PA's null
  MRR must ≈ its R0 MRR (~0.58). If PA moves, the swap changed degrees (bug).
- **Output both:** one canonical `train_R2_degree_null_seed42.csv` (so `run_baselines --regimes R2`
  works as a point estimate) **and** N replicates for the null distribution + permutation p-value.

**Accept:** degree sequence preserved (assert); PA null ≈ PA R0; CN/AA/Jaccard get a real null MRR
with a permutation p-value = fraction of null replicates whose MRR ≥ the real R0 MRR (per method).
**Interpretation is reportable either way:** big real−null gap = genuine structure; small gap =
performance was mostly degree.

**Prompt:**
```
Write scripts/build_degree_null.py. Purpose: a degree-preserving, TYPE-preserving permutation null
of the R0 training graph for the leakage audit. Read data/processed/splits/train.csv
(source_id,relation,target_id). For each relation type separately, build an igraph graph and run
degree-preserving double-edge swaps (~10x the edge count) so every node keeps its exact degree and
edges never cross relation types; recombine into a full edge list. Generate N replicates (default
10) with seeds derived from --seed; assert for each that the per-node degree sequence is identical
to the original and per-relation edge counts are unchanged. Write the first replicate to
data/processed/splits/train_R2_degree_null_seed42.csv (so lib_eval's R2 regime resolves) and all
replicates to data/processed/splits/null/degree_null_rep<k>.csv. Then, importing scripts/lib_eval.py
and the scorers from scripts/run_baselines.py (Random/CommonNeighbors/AdamicAdar/Jaccard/
PreferentialAttachment), score each replicate on the R0 test edges while HOLDING THE EVALUATION
FIXED: build pools/known ONCE from R0 (lib_eval._category_pools / _known_tails_by_head on
train.csv) and pass them into rank_test_edges so only the scoring adjacency changes per replicate.
Report, per method: real R0 MRR, mean null MRR, and a permutation p-value = fraction of replicates
with null MRR >= real MRR. Save data/processed/results/null/degree_null.json. Deterministic; print a
table. Sanity-assert PreferentialAttachment's null MRR ~= its R0 MRR.
```

## Unit K — `run_kge.py` (TransE + ComplEx)  🔴 A (GPU) or Colab · *mostly anytime (parallel)*

**Why:** KGE is multi-hop, so **R3 (orthology) may drop KGE even though it didn't touch the
baselines** — the one place the orthology-leakage claim can be tested. Also the "embeddings inflate
too" evidence and the AUROC≫MRR dissociation.

**In depth:**
- **Trim** the existing 4-model `run_kge.py` to **TransE + ComplEx** only (two families), dim 64,
  standard cited hyperparameters + one small R0 sweep, seeds {42,1,7}.
- **Score through the shared harness.** After training, wrap the model as
  `score_fn(gene,disease) = model.score_hrt((gene_idx, TARGET_REL_idx, disease_idx))` and feed it to
  `lib_eval.rank_test_edges` — do NOT use PyKEEN's internal evaluator, or KGE and baselines won't be
  comparable. Map string IDs → PyKEEN entity/relation indices via the trained TriplesFactory; for
  a candidate not seen in training, return a low sentinel (score never skipped).
- **Filtered subgraph (Gu 2024) for speed:** train on ~11% of the graph (cite the filter) so each
  run is tens of minutes; full-graph dim-64 KGE is trainable but slow.
- **Regimes:** R0 (train.csv), R2 (degree null), R3 (orthology). R1≈R0 for KGE (optional).
- **Known gap in the old script:** it only logged KGE metrics; make it WRITE per-run JSON.

**Accept:** ranking+classification per model×regime×seed with mean±sd in
`data/processed/results/kge/`. Key comparison: **is KGE R3 MRR < KGE R0 MRR** (orthology bites
KGE)? Report either way. Expect AUROC high but filtered MRR low (the dissociation).

**Prompt:**
```
Refactor scripts/run_kge.py to (1) train ONLY TransE and ComplEx (PyKEEN, dim 64, fixed seeds
42,1,7, standard published hyperparameters + one small R0 sweep over dim in {64,128} epochs in
{100,300}); (2) train on a regime's training graph from data/processed/splits/ (R0=train.csv,
R2=train_R2_degree_null_seed42.csv, R3=train_R3_orthology_blocked.csv), optionally on an ~11%
filtered subgraph for speed (cite Gu et al. 2024); (3) score the held-out gene->disease test edges
through scripts/lib_eval.rank_test_edges by wrapping the model as score_fn(gene,disease)=
model.score_hrt for the GENE_ASSOCIATED_WITH_CONDITION relation, mapping IDs via the TriplesFactory
and returning a low sentinel for unseen entities; (4) loop seeds and regimes; (5) write every run's
ranking (MRR,Hits@1/3/10) and classification (AUROC,AUPRC) to
data/processed/results/kge/kge_<model>_<regime>_seed<k>.json plus a merged summary with mean+/-sd.
Keep GPU support and config.py paths. Print R0 vs R3 MRR per model so I can see if orthology removal
drops KGE. Show me the diff before any long training run.
```

## Unit 9 — `hetionet_audit.py` (robustness)  🟠 B / runs on A · *anytime (parallel)*

**Why:** rerun the audit on an independent graph → generality.

**In depth:**
- **ID remap is the trick.** Hetionet IDs are `Type::id`; `lib_eval` calls Monarch's `category_of`
  internally. Rather than fork the harness, **remap Hetionet IDs to prefixes `category_of`
  understands**: `Gene::1234 → NCBIGENE:1234` (gene), `Disease::DOID:x → DOID:x` (disease),
  `Compound::DB → DRUGBANK:DB`, etc. Then `lib_eval` runs unmodified.
- **Target = `DaG`** (Disease–associates–Gene); orient as (gene, disease). Hold out 10/10, build
  R0, R1 (redundancy), R2 (degree null via Unit 7's function). **No R3** — Hetionet has no
  cross-species orthology; state this explicitly. **Never merge** Hetionet with Monarch.
- The generality question: does **PA/degree dominate here too**, and does the **R2 degree null**
  behave the same? That's the claim, not the absolute numbers.

**Accept:** `data/processed/results/hetionet/baselines_<regime>.json`; the R0-vs-R2 pattern
reproduces (or is honestly reported not to).

**Prompt:**
```
Write scripts/hetionet_audit.py. Load data/external/hetionet/hetionet-v1.0-edges.sif
(tab-separated source,metaedge,target; IDs "Type::id"). Remap IDs to prefixes that
scripts/kg_categories.category_of understands: Gene::N->NCBIGENE:N, Disease::DOID:x->DOID:x,
Compound::DB->DRUGBANK:DB, and pass other types through as Type:id. Treat the DaG
(Disease-associates-Gene) metaedge as the gene->disease target task (gene=Gene node,
disease=Disease node). Reusing the SAME split logic as build_deleaked_splits.py and the SAME
harness lib_eval.py + the scorers from run_baselines.py, build regimes R0, R1 (remove direct
held-pair edges), and R2 (degree-preserving type-preserving null via build_degree_null); SKIP R3
(Hetionet has no orthology). Run the 5 baselines x 3 seeds and write
data/processed/results/hetionet/baselines_<regime>.json with mean+/-95%CI. Do not merge with
Monarch. Print an MRR table so I can compare the R0->R2 pattern to Monarch.
```

## Unit 10 — `make_figures.py`  🟢 A/C · *after 7, K, 9*

**Why:** the three figures that carry the paper. **In depth:** matplotlib only; read the results
JSONs (no hand-typed numbers); colorblind-safe; 300 dpi; deterministic regeneration.
- **Fig 1 — dissociation:** per method, AUROC vs MRR (KGE tall-AUROC/low-MRR; baselines reverse).
- **Fig 2 — audit:** MRR per method across R0→R1→R2→R3 with 95% CI bands. Honest shape: R1/R3 flat
  for baselines, the **R2 drop is the visible effect**; KGE panel shows any R3 drop.
- **Fig 3 — degree-vs-structure:** per method, R0 MRR split into degree component (null MRR) and
  structure component (R0−null), with the permutation p-value annotated.
**Accept:** `figures/fig1..3.{png,pdf}` regenerate from JSON with zero manual editing.
**Prompt:**
```
Write scripts/make_figures.py (matplotlib only, colorblind-safe, 300 dpi, deterministic). Read
data/processed/results/baselines/baselines_R*.json, .../kge/*.json, and .../null/degree_null.json.
Produce: Fig 1 = per-method AUROC vs MRR scatter (label methods; show the KGE-high-AUROC/low-MRR vs
baseline-reverse dissociation). Fig 2 = MRR per method across R0,R1,R2,R3 as lines with 95% CI
bands (separate panels for baselines and KGE). Fig 3 = per-method R0 MRR decomposed into degree
component (mean null MRR) and structure component (R0 minus null), stacked bars with the
permutation p-value annotated. Save figures/fig1..3.png and .pdf. No hand-typed numbers; every value
comes from the JSONs.
```

## Unit 11 — `make_tables.py`  🟢 A/C · *after 7, K, 9*
**In depth:** emit markdown + LaTeX from JSON; every cell traceable. **Accept:** no hand-typed
numbers. **Prompt:**
```
Write scripts/make_tables.py. From data/processed/graph_stats*.json and
data/processed/results/**/*.json, emit markdown + LaTeX for: Table 1 (graph stats, Monarch +
Hetionet), Table 2 (benchmark: every method x R0-R3, MRR/Hits@k/AUROC as mean+/-95%CI), Table 3
(degree-null decomposition: real MRR, null MRR, permutation p-value per method), Table 4 (Hetionet
R0/R1/R2). Write to tables/ as .md and .tex. Every number read from JSON; fail loudly if a source
file is missing rather than hard-coding.
```

## Unit 12 — `run_all.ps1` (+ `.sh`) + Zenodo  🟢 A/C · *after figures/tables*
**In depth:** ordered end-to-end script with runtimes/hardware noted; then a GitHub release wired to
Zenodo for DOIs. **Accept:** a clean clone reproduces Table 2 + Fig 2. **Prompt:**
```
Write run_all.ps1 and run_all.sh that reproduce the whole pipeline in order with echoed step
banners and expected runtimes: graph_stats (Monarch+Hetionet) -> build_deleaked_splits ->
build_degree_null -> run_baselines (R0,R1,R2,R3) -> run_kge -> hetionet_audit -> make_figures ->
make_tables. Verify data/processed/edges_clean_integrated.csv and the Hetionet files against
ARTIFACT_HASHES.txt before running and abort on mismatch. Print where each output lands. Assume a
fresh clone + the data channel; keep paths relative.
```

## Unit 13 — manuscript  🟢 · *start intro/methods ANYTIME; finalize after results*
**In depth:** resource-led IMRaD. **Lead with** the released de-leaked benchmark + the graded
decomposition. **Honest headline findings:** (i) degree dominates baseline performance (R2);
(ii) redundancy negligible; (iii) orthology removal doesn't affect 2-hop topological methods
(mechanistic, ≥3-hop) and its KGE effect is [Unit K]; (iv) replicates on Hetionet. **Cite +
distinguish** Gu 2024, Ranga 2025, Alghamdi/Hoehndorf/Robinson 2022 — novelty is incremental
methodological + resource, don't oversell. Adversarial self-review; a senior reader; bioRxiv on
completion; resource-tier venue.

---

## Timeline — exact order, and what's parallel/anytime

```
CRITICAL PATH (do in this order; each gates the next):
  Unit 7 (R2 degree null + baselines on R2)      <-- NEXT; the decisive result
     |
     v
  Unit 10 + 11 (figures + tables)  <-- need results from 7 AND K AND 9
     |
     v
  Unit 12 (run_all + Zenodo)
     |
     v
  Unit 13 finalize (numbers locked) -> preprint -> submit

PARALLEL LANE (start anytime; must LAND BEFORE Unit 10 figures):
  Unit K (KGE: TransE+ComplEx)   -- R0/R3 can run now; R2 after Unit 7. GPU/Colab, off critical path.
  Unit 9 (Hetionet audit)        -- fully independent of Monarch R2/KGE; needs Unit 7's null function.

ANYTIME / NO DEPENDENCY (whenever convenient; don't let them block the critical path):
  Unit 1b (Hetionet graph_stats)         -- Hetionet already downloaded
  Unit 3  (provenance fill)              -- manual lookup, yours
  Housekeeping (track result JSONs; reconcile AA 0.69)
  Unit 13 drafting (intro / related work / methods) -- write while compute runs
```

**Plain-English order:**
1. **Now:** Unit 7 (build the degree null, score baselines on R2). This is the one result the
   baseline story needs.
2. **In parallel with / right after 7:** Unit K (KGE — needs the GPU/Colab; R3 tests orthology) and
   Unit 9 (Hetionet — independent). These two + Unit 7 are the three that feed the figures.
3. **Whenever you have a spare moment** (they gate nothing): Unit 1b, Unit 3, housekeeping, and
   starting to *write* the manuscript's intro/methods/related-work.
4. **Once 7 + K + 9 are done:** Units 10–11 (figures + tables), then Unit 12 (release), then finalize
   Unit 13 and submit.

**Single-machine note:** solo on A, do the critical path in order and slot the parallel/anytime items
into gaps (e.g. run Unit 9 while KGE trains). Only split across B/C if you want KGE training and the
degree null/Hetionet running literally at the same time, or a second person writing.

---

## Updated checklist (remaining)

- [ ] 7.  `build_degree_null.py` → R2 files + decomposition  🟠/A  ← **NEXT & DECISIVE**
- [ ] K.  `run_kge.py` trimmed → TransE+ComplEx on R0/R2/R3 (orthology test)  🔴 A/Colab · *parallel*
- [ ] 9.  `hetionet_audit.py` → second-graph tables (R0/R1/R2)  🟠/A · *parallel*
- [ ] 1b. `graph_stats.py` on Hetionet  🟢 · *anytime*
- [ ] 3.  provenance versions/dates/licenses  🟢 · *anytime, manual*
- [ ] H.  housekeeping: track result JSONs + reconcile AA 0.69  🟢 · *anytime*
- [ ] 10. `make_figures.py` → fig1–3  🟢 · *after 7,K,9*
- [ ] 11. `make_tables.py` → tables 1–4  🟢 · *after 7,K,9*
- [ ] 12. `run_all.ps1` + Zenodo DOIs  🟢 · *after 10,11*
- [ ] 13. manuscript (draft anytime; finalize last) → preprint → submit  🟢
