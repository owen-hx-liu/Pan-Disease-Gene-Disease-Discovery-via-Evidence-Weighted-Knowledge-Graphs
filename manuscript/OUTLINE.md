# Manuscript outline — gene→disease link-prediction leakage audit

**Status:** working draft, untracked (do NOT commit to the public `release/v1.0.0` branch — this is an internal doc like the ones pruned in commit 67a52a2).
**Companion:** `MASTER_REFERENCE.md` (the big pull-from-here source file).
**DOI to embed:** `10.5281/zenodo.21398832`.

---

## Framing (hold every section to this)

**Thesis (one sentence):** Evaluation leakage in gene→disease link prediction is *method-class-dependent* — node-degree structure is the one universal leak, train/test redundancy leaks only into multi-hop KGE, and cross-species orthology (the intuitive culprit) is a null.

**Genre:** resource + audit paper. Soundness and honesty are the selling points, not a novel positive finding.

**Novelty posture (critical):** the core phenomenon (random-split evaluation is inflated by degree / study bias / leakage) is *already established* — PNAS (Yılmaz 2025), on Monarch itself (Bradshaw & Layer 2024 — VERIFY name, see MASTER §Citations), and in a near-identical leakage framework (Brière & Baudot 2025). **Lead Related Work with these three and frame the paper as careful consolidation + cross-graph replication + a degree-magnitude decomposition, not a discovery.** Over-claiming novelty is the most likely rejection cause.

**Realistic venue:** *Journal of Emerging Investigators (JEI)* + bioRxiv preprint as the guaranteed mentored track. A reproducibility/benchmark venue (*Database (Oxford)*, *NAR Genomics & Bioinformatics*, *Bioinformatics Advances*) is a stretch-but-possible IF the two flagged gaps (biological case study + low-degree/stratified eval) are added. Add a **senior/mentor co-author** either way (required for JEI).

---

## Estimated length

| Version | Main text | Display items | Refs | Typeset |
|---|---|---|---|---|
| Full (resource venue) | ~6,000–6,800 words | 4 tables + 3 figs + 1 supp table (all already generated) + 2 optional figs | ~30–35 | ~9–12 pp |
| JEI version | ~3,000–3,500 words | trim to 3 tables + 2 figs | ~20–25 | ~6–8 pp |

Per-section word budgets are in the headings below (full version).

---

## Front matter

- **Title:** "Quantifying evaluation leakage in multi-species gene–disease link prediction: a graded audit and a de-leaked benchmark." (matches `.zenodo.json`)
- **Authors:** Liu, Tang, Augustine + **mentor co-author (TBD)**.
- **Abstract (~250 words, structured):** problem (KGE GDA scores are inflated) → what we built (one shared harness; 2 KGE + 5 baselines; 4 graded regimes R0–R3; 3 seeds; two graphs) → method-class-dependent headline with the three numbers (R2 −43–49% baselines / −21–24% KGE; R1 KGE-only; R3 null) → de-leaked benchmark + DOI. **Put DOI `10.5281/zenodo.21398832` here.**
- **Keywords:** reuse the `.zenodo.json` list.
- **Author contributions + competing interests:** short, required by most venues (currently a ⚠️ gap).

## 1. Introduction (~700 words)
- KGE GDA scores drive real biological prioritization; several *distinct* leakage mechanisms inflate them, usually reported as one undifferentiated "optimism."
- Cite the general ML leakage crisis (Kapoor & Narayanan 2023) + the biological leakage taxonomy (Bernett 2024a).
- Gap: no *graded, controlled decomposition* separating the mechanisms on one graph under one harness, with a released de-leaked split and a formal degree-preserving null that measures *magnitude*.
- Contributions (enumerate 4): (1) graded regime framework R0–R3; (2) a degree-preserving null that isolates degree from genuine structure and quantifies each; (3) the method-class-dependent finding with orthology & redundancy as *negative controls*; (4) a reproducible, DOI-archived de-leaked benchmark, replicated on a second graph.
- **Honesty paragraph:** state up front this consolidates and sharpens known results.

## 2. Related work (~600 words) — LEAD with the competitors
- **2.1 The three direct competitors** (out-position, don't bury): Yılmaz PNAS 2025 (degree bias in PPI eval; PA-as-pure-bias), Bradshaw & Layer 2024 (KGE degree bias on Monarch — same graph/task), Brière & Baudot 2025 (3 leakage types on a KGE drug→disease pipeline). What each owns; what we add. → **Related-work differentiation table (MASTER Table D).**
- **2.2 Method/foundation lineage:** Akrami 2020 (origin of R1), Toutanova & Chen 2015 (FB15k-237 reversible-relation removal), Maslov & Sneppen 2002 + Bernett 2024b (degree-preserving permutation → R2), Aiyappa 2024 / Shomer 2023 (degree bias in link prediction / KGE).
- **2.3 The reconciliation to flag here and resolve in Discussion:** Brière's DL2 found *no* degree reliance; our R2 finds degree is *the* universal leak. Preview why (binary retrain vs magnitude decomposition; task/graph differences).

## 3. Materials: the knowledge graph (~700 words)
- Construction: **~98.5% Monarch-derived** + 7 curated sources (name them). Honest scope note — Monarch is itself an aggregator, so "N databases" ≠ N independent votes.
- Scale/composition → **Table 1** (452,012 nodes; 5,860,539 edges; 28 relations; 8 categories; degree mean 25.9 / median 2 / max 20,662; 42,288 gene–disease target edges; 510,703 orthology edges).
- Second graph: Hetionet v1.0 (Himmelstein 2017) for replication — stats in the same Table 1.
- Provenance & licensing (Monarch CC BY 4.0; per-source contribution table; orthology data basis for R3).
- **⚠️ Recommended figures (currently gaps):** schema/metagraph figure + degree-distribution figure (both flagged Med priority in your ref doc).

## 4. Methods (~1,400 words)
- **4.1 Graded leakage regimes.** R0 standard random split; R1 redundancy (1,661 direct held-pair edges removed, 0.03% of edges); R2 degree-preserving null; R3 orthology-blocked (~510–520K `ORTHOLOGOUS_TO` removed). Spell out exactly what each removes and why.
- **4.2 The degree-preserving null (R2).** Rewiring that holds the degree sequence exactly; verified by PreferentialAttachment being unchanged. Cite Maslov–Sneppen 2002; `build_degree_null.py`.
- **4.3 Models.** TransE + RotatE (d=64, 300 epochs, NSSA loss, seeds 42/1/7). **State ComplEx was attempted and is untrainable at this scale** (smoke tests only); DistMult likewise only smoke-tested — pre-empts "why not more KGE families?". Five topological baselines (Random, CN, AA, Jaccard, PrefAttach).
- **4.4 Shared evaluation harness.** One `lib_eval.py` ranks every method identically — the fairness backbone.
- **4.5 Ranking protocols.** Two, and NAME which every number comes from: sampled-50 type-matched negatives (main) and filtered full-ranking vs the ~14.3k-disease pool (robustness). Cold-start handling (30 genes → worst rank).
- **4.6 Statistics.** Paired bootstrap for MRR differences (1,000 resamples); 10-replicate permutation test for the degree decomposition; 3 seeds throughout; bootstrap CIs. **This rigor is the paper's biggest edge — foreground it.**

## 5. Results (~1,800 words) — 7 short subsections
- **5.1 The graded audit → Table 2 + Figure 2.** Headline MRR across R0–R3 for all 7 methods; lead with method-class-dependence.
- **5.2 Degree is the universal leak (R2) → Table 3 + Figure 3.** Overlap heuristics ≈ 55% degree / 45% structure (permutation p<0.001); PrefAttach 100% degree, unchanged by the null (built-in positive control).
- **5.3 Redundancy is a multi-hop-only leak (R1).** KGE −6 to −8.5%, 2-hop baselines flat (±0.3%). Mechanism: a direct gene→disease edge is embedded as proximity by KGE but is not a shared neighbor. Pre-empt "R1 is just less training data" (0.03% of edges).
- **5.4 Orthology is a null (R3).** All methods within ±0.5%. The refuted domain hypothesis / designed negative control.
- **5.5 Full-ranking robustness → Table S1.** TransE R0/R1/R2 full-rank: R2 collapse *strengthens* to −71%, R1 −10%. Name the protocol; note the MRR-vs-median dissociation nuance.
- **5.6 Independent-graph replication → Table 4.** Hetionet reproduces the degree collapse (−35 to −37% at R2; PrefAttach flat). Not a Monarch artifact; Monarch is *more* degree-dominated.
- **5.7 AUROC–MRR dissociation → Figure 1.** AUROC stays ~0.96–0.99 as MRR collapses. Why AUROC hides leakage and ranking metrics expose it (cite Saito & Rehmsmeier 2015).

## 6. Discussion (~900 words)
- Practice: report R2 as a mandatory sanity check; degree baselines belong in every GDA benchmark; AUROC alone is inadequate.
- **Reconcile Brière DL2** (binary retrain vs magnitude decomposition; their DL3 non-representative-test result actually supports us).
- Short recommendations checklist for the field.
- **Limitations (be generous):** no genuine temporal split; Monarch-dominated coverage; incremental novelty over PNAS/Brière/Bradshaw; KGE full-rank is TransE-only (justify as representative); sampled-vs-full protocol sensitivity; **two acknowledged gaps** — no biological case study, no low-degree/stratified eval — state whether added or flagged as future work.

## 7. Conclusion (~200 words)
- Restate the three-way method-class-dependent finding; point to the archived benchmark.

## 8. Data & code availability
- GitHub repo + **Zenodo DOI `10.5281/zenodo.21398832`**; MIT (code) / CC BY 4.0 (data); one-command `run_all.sh`; frozen result JSONs; split manifest SHA-256 `bc468864…`.

---

## Display-item inventory (all already generated)

| Item | File | Section |
|---|---|---|
| Table 1 — graph stats | `tables/table1_graph_stats.md` | §3 |
| Table 2 — leakage audit (sampled-50) | `tables/table2_audit.md` | §5.1 |
| Table 3 — degree-null decomposition | `tables/table3_degree_null.md` | §5.2 |
| Table 4 — Hetionet replication | `tables/table4_hetionet.md` | §5.6 |
| Table S1 — full-ranking robustness | `tables/tableS1_fullrank_robustness.md` | §5.5 |
| Figure 1 — AUROC–MRR dissociation | `figures/fig1_auroc_mrr_dissociation.*` | §5.7 |
| Figure 2 — leakage audit | `figures/fig2_leakage_audit.*` | §5.1 |
| Figure 3 — degree vs structure | `figures/fig3_degree_vs_structure.*` | §5.2 |
| Fig S? — schema/metagraph | ⚠️ not yet generated | §3 |
| Fig S? — degree distribution | ⚠️ not yet generated | §3 |

## What's left to DO (vs already done)
Everything numerical is frozen and reproducible. Net-new work is prose + positioning + (optionally) two experiments:
1. Write Intro, Related Work (lead with the 3 competitors), Discussion, Limitations.
2. **Recruit the mentor co-author** — biggest odds-mover; required for JEI.
3. Post **bioRxiv preprint**; cross-link its DOI into `.zenodo.json` (the one remaining metadata TODO).
4. **Optional, for a benchmark venue:** biological case study (rare-disease top-k + enrichment) + low-degree/stratified eval. Not needed for JEI.
5. **Optional:** add DistMult as a 3rd KGE family (feasible; smoke test exists). ComplEx is a dead end at this scale — do not promise it.
6. Verify every citation (see MASTER §Citations — esp. the Gu/Bradshaw name conflict and preprint→published-version updates).
