# PAPER_SCOPE.md — Locked scope & the novelty that makes it matter

*Created 2026-07-07. This is the contract. If a task, figure, or paragraph does not serve
one of the contributions C1–C5 below, it is out of scope. Re-read before adding anything.*

---

## 0. The problem with the "safe" framing (why we are NOT doing a plain benchmark)

A paper that only says *"Adamic-Adar beats KGE for gene–disease link prediction and AUROC is
misleading"* will be **rejected as not novel**. Reviewer 2 will cite the existing literature
that already shows (a) simple topological baselines rival KG embeddings, and (b) AUROC is
optimistic under class imbalance. That version is fodder. We need a finding the field does
**not** already know and **should** act on.

We have one, and it is specific to this dataset.

---

## 1. The core insight (novel, important, and confirmed in our data)

**This graph is one of the most multi-species biomedical KGs in common use (Monarch-derived),
and that structure secretly leaks the answer to the gene–disease prediction task.**

Confirmed facts from `data/processed/edges_clean_integrated.csv` (2026-07-07):

- The prediction target — `GENE_ASSOCIATED_WITH_CONDITION` — is **42,288 edges, 100% human
  (all source nodes are HGNC)**. So the task is *human* gene → disease.
- The graph contains **510,703 `ORTHOLOGOUS_TO` edges (8.7% of all edges)** that bridge those
  human genes to model-organism orthologs: HGNC→FB 22,340, HGNC→WB 22,174, HGNC→ZFIN 18,019,
  Xenbase→HGNC 27,940, MGI/RGD/ZFIN cross-links, etc.
- Those orthologs carry the bulk of **1,066,114 `HAS_PHENOTYPE` edges** and there are **9,495
  `MODEL_OF`** (disease ↔ model-organism) edges.

**The leakage mechanism (the heart of the paper):** a held-out human edge
`HGNC:gene —GENE_ASSOCIATED_WITH_CONDITION→ MONDO:disease` can be recovered *without any
independent human evidence*, by traversing:

```
HGNC:gene --ORTHOLOGOUS_TO--> mouse/fish ortholog --HAS_PHENOTYPE--> phenotype P
MONDO:disease --HAS_PHENOTYPE--> phenotype P     (and/or  MONDO:disease --MODEL_OF--> ortholog)
```

Neighborhood methods (Adamic-Adar/common-neighbors) and message-passing/KGE models can score
the true partner highly purely through this orthology–phenotype shortcut. **The model looks
like it is "discovering" human disease genes when it is really re-deriving them from mouse.**
That is a form of **evaluation leakage that standard random splits do not control for**, and
— to our knowledge — it is **not quantified or corrected in the multi-species biomedical KG
link-prediction benchmarks the field relies on** (Monarch KG, Hetionet-style, PrimeKG-style
pipelines all inherit this structure).

**Why it matters (the "who cares"):** a large fraction of computational gene–disease and
drug-repurposing work is benchmarked on exactly this kind of aggregated, ortholog-rich graph.
If the headline performance is substantially orthology leakage, the field is **systematically
over-reporting** its ability to find *novel human* associations, and prioritizing candidates
for the wrong reasons. Showing this rigorously — and providing a corrected protocol — is a
genuine, useful, publishable contribution. It is a "the evaluation is leaky, here is proof and
a fix" paper, which good venues actively want.

---

## 2. Working titles (pick one at draft time)

1. **"Orthology leakage inflates gene–disease link prediction on multi-species knowledge
   graphs."** *(recommended — names the finding)*
2. "Are we discovering disease genes or re-deriving them from mice? A leakage-aware benchmark
   for biomedical knowledge-graph link prediction."
3. "Beyond AUROC: hub and orthology shortcuts in gene–disease link prediction, and a de-leaked
   evaluation protocol."

---

## 3. Central hypothesis and the killer experiment

**H1 (leakage).** A large share of held-out human gene–disease ranking performance is
attributable to cross-species orthology + shared-phenotype shortcuts, not independent
evidence.

**Killer experiment (the figure that sells the paper):** measure every method's ranking
(MRR, Hits@k) under four evaluation regimes on the *same* held-out human gene–disease edges:

| Regime | What it removes | Tests for |
|---|---|---|
| **R0 Standard** | nothing (random split; current practice) | the inflated baseline everyone reports |
| **R1 Hub-controlled** | generic hubs (deg > 2,000) as shared-neighbor evidence | hub shortcut |
| **R2 Orthology-blocked** | orthology + model-organism phenotype paths to the target | **the leakage claim** |
| **R3 De-leaked** | R1 + R2 + tautological/eponymous + source-non-independence | honest performance |

**Prediction:** performance falls sharply from R0 → R3 for *all* methods, and the **method
ranking can change** (a method that "won" under R0 may not under R3). The magnitude of the
R0→R3 drop is the paper's headline number. Even if the drop is moderate, quantifying it and
releasing the protocol is a contribution; if it is large, it is a strong result.

*Honesty guardrail:* we report whatever the drop actually is. If orthology leakage turns out
small, that is itself a publishable, reassuring finding ("multi-species KGs are safer than
feared") — we do not need a big number to have a paper, we need a *rigorous* one.

---

## 4. Contributions (the five things the paper delivers)

- **C1 — A leakage taxonomy for biomedical KG link prediction.** Formalize four leakage
  sources on aggregated multi-species KGs: (i) cross-species orthology shortcuts, (ii)
  generic-hub shortcuts, (iii) eponymous/tautological gene–disease edges, (iv) source
  non-independence (aggregator double-counting). Novel synthesis; (i) is the new one.
- **C2 — A de-leaked evaluation protocol + open harness.** Concrete, reusable splits and
  filters (R0–R3 above) that any multi-species KG benchmark can adopt. Released as code.
- **C3 — The measured impact.** How much each method's reported performance is leakage, under
  one strict protocol, with ≥3 seeds and significance tests. Includes the AUROC-vs-ranking
  dissociation as a *sub-result*, now explained mechanistically by the leakage, not just observed.
- **C4 — A leakage-robust method.** A topology-seeded, embedding-reranked hybrid built to
  degrade gracefully under de-leaking, especially for **low-degree rare-disease nodes** where
  orthology evidence is thin and current methods fail hardest. This is *our* method — it
  answers "what did you invent," and its win-condition is robustness under R3, not R0.
- **C5 — Prospective-surrogate rare-disease validation.** Leave-one-source-out + newer-release
  recovery + curated (de-tautologized) literature confirmation, focused on Orphanet rare
  diseases — showing the de-leaked pipeline still surfaces real, checkable novel candidates.

**One-sentence contribution statement for the intro:** *"We identify and quantify
cross-species orthology leakage as a previously under-examined confound in gene–disease link
prediction on multi-species knowledge graphs, release a de-leaked evaluation protocol under
which reported performance drops by X and the method ranking changes, and introduce a
leakage-robust hybrid predictor validated on rare-disease gene discovery."*

---

## 5. How we defend novelty against "this is already known"

| Anticipated reviewer objection | Our defense (must be in the paper) |
|---|---|
| "Topology beats KGE is old news" | We agree and cite it; it is a *sub-result*. Our novelty is *why* (leakage) and the *corrected protocol*, not the horse race. |
| "AUROC-is-misleading is known" | Cited; we go further — we show the *ranking* signal itself is largely leakage. |
| "Data leakage in ML is well studied" | Yes in general; **cross-species orthology leakage in KG link prediction is not** a controlled variable in the standard benchmarks. That specificity is the contribution. |
| "Just tune KGE harder" | Phase-2 sweep + best-config + mean±sd pre-empts it; and the point is invariance to method, not KGE-bashing. |
| "Not novel data" | Correct — we explicitly do **not** claim novel data; the graph is Monarch-derived. Novelty is evaluation + method. Stated plainly in Related Work. |

---

## 6. In scope / out of scope (freeze this)

**In scope (build these):**
- The frozen integrated graph as the substrate (Monarch-derived + 7 curated sources).
- Human HGNC gene → MONDO disease as the single prediction task.
- Methods: Random, Common-Neighbors, Adamic-Adar (baselines); TransE, DistMult, ComplEx,
  RotatE (KGE, tuned, ≥3 seeds); our hybrid (C4).
- Evaluation regimes R0–R3 (C2); leakage taxonomy (C1); measured impact + significance (C3).
- Validation: LOSO, newer-release recovery, curated literature (C5).
- Network structure vs null model — kept as a short "graph characterization" subsection only.

**Out of scope (archive; do NOT let these into the paper):**
- Any "we integrated 13 databases / novel resource" claim.
- Temporal growth analysis (fabricated years) and any time-based split.
- Neo4j drug-repurposing (`drugrepurposing.py`), the Streamlit portal, the HTML explorers.
- The old 606-node consensus subset and `transe_new_predictions.csv` (superseded/degenerate).
- Consensus confidence beyond one honest paragraph (it only discriminates ~4,851 edges).

---

## 7. Frozen artifacts (single source of truth — hash these)

- **Graph:** `data/processed/edges_clean_integrated.csv` — 448,931 nodes / 5,860,539 edges.
  Compute and record its `sha256`; every machine must use the byte-identical file.
- **Stats:** `data/processed/graph_stats.json` (to be generated in Phase 1) → Table 1.
- **Benchmark:** `benchmark_results.json` (baselines + KGE, per-seed) → Table 2.
- **Leakage:** `data/processed/leakage/` (R0–R3 results per method) → Figure 1 + Table 3.
- **Validation:** `data/processed/loso/`, `.../newer_release/`, `.../literature/` → Table 4.

---

## 8. Success criteria (what "done and publishable" means)

1. R0→R3 leakage drop is measured, with 95% CIs over ≥3 seeds, for every method.
2. At least one pairwise method difference is tested for significance (paired bootstrap on
   per-edge reciprocal ranks) — not eyeballed.
3. The hybrid (C4) is the best or tied-best method under **R3** (the honest regime), even if
   not under R0.
4. C5 yields ≥5 literature-confirmed, non-tautological novel rare-disease candidates *and* a
   quantitative newer-release/LOSO recovery number.
5. Every number regenerates from a frozen, hashed graph via one command on a clean clone.
6. A domain expert (mentor/co-author) has read and signed off on the biology.

---

## 9. Target venues (unchanged from roadmap, re-affirmed)

- **Reach / primary:** BMC Bioinformatics, PLOS ONE, PeerJ, *Bioinformatics* (if framed as a
  methods/benchmark advance), *Database (Oxford)*, GigaScience.
- **Guaranteed parallel track:** Journal of Emerging Investigators (mentored, honest, real).
- **Always:** bioRxiv preprint the day the draft is done; consider an ISMB/PSB or NeurIPS
  Graph-Learning-in-Bio workshop version.

---

## 10. The 3-machine execution plan

The operational plan for splitting this across three computers (roles, job manifest,
coordination, and the "make it the best" excellence practices) lives in
**`WORK_SPLIT_3_MACHINES.md`**. Read it before starting Phase 2.

*Next action after this file: generate `graph_stats.json` (Phase 1) and build the R0–R3 split
generator (`scripts/build_deleaked_splits.py`), which is the backbone every machine depends on.*
