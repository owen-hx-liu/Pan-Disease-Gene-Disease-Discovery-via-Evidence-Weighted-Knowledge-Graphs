# MASTER REFERENCE — everything for the paper, sectioned by the outline

**Purpose:** one place to pull every fact, number, citation, talking point, and draft sentence from. Sectioned to match `OUTLINE.md`.
**Status:** internal working draft — untracked. Do NOT commit to `release/v1.0.0` (public branch).
**Provenance of numbers:** all result numbers below are taken from the frozen JSONs / generated tables in this repo (`tables/`, `data/processed/results/`). Where a fact comes from your literature review (`science fair year 2.md`), it is marked **[lit-doc]** and flagged VERIFY where the bibliographic detail is unconfirmed.
**DOI:** `10.5281/zenodo.21398832` · **Split manifest SHA-256:** `bc468864cdc3c395f7d5424c27c34edde7354d53284bfe876eb94a00010a08d1`

---

## 0. QUICK-REFERENCE: the numbers you will cite most

**Graph (Monarch, this work):** 452,012 nodes · 5,860,539 directed edges · 5,583,381 unique undirected · 28 relation types · 8 node categories · degree mean 25.9 / median 2 / max 20,662 · 336,460 gene nodes · 14,339 disease nodes · 42,288 gene–disease target edges · 510,703 orthology edges.

**Graph (Hetionet v1.0, robustness):** 47,031 nodes · 2,250,197 directed edges · 24 metaedge types · 11 metanodes · degree mean 99.7 / median 32 / max 23,879 · 19,145 gene nodes · 136 disease nodes · 12,623 gene–disease edges · 0 orthology (single-species).

**Task:** predict `BIOLINK:GENE_ASSOCIATED_WITH_CONDITION`; **4,228 held-out human gene–disease test edges**; seeds 42/1/7.

**Regimes:** R0 standard · R1 redundancy (**1,661 direct held-pair edges removed = 0.03%**) · R2 degree-preserving null · R3 orthology-blocked (**520,198 edges removed**, R0 train 5,852,083 → R3 train 5,331,885; Table 1 counts 510,703 ORTHOLOGOUS_TO — the ~10k gap is reverse-direction edges).

**Headline MRR (sampled-50 negatives, Table 2), R0 → R2 (%Δ), plus R1 / R3:**

| Method | R0 | R1 | R2 | R3 | R0→R2 |
|---|---|---|---|---|---|
| Random | 0.091 | 0.091 | 0.091 | 0.091 | +0% |
| CommonNeighbors | 0.441 | 0.442 | 0.247 | 0.443 | **−44%** |
| AdamicAdar | 0.459 | 0.459 | 0.262 | 0.460 | **−43%** |
| Jaccard | 0.397 | 0.398 | 0.201 | 0.399 | **−49%** |
| PreferentialAttachment | 0.587 | 0.582 | 0.589 | 0.581 | **+0%** (positive control) |
| TransE | 0.799 | 0.750 | 0.611 | 0.800 | **−24%** |
| RotatE | 0.827 | 0.757 | 0.655 | 0.830 | **−21%** |

**R1 KGE-only leak:** TransE −6.2%, RotatE −8.5%; all 2-hop baselines ±0.3%.
**R3 null:** every method within ±0.5% (KGE +0.2–0.4%).
**Degree decomposition (Table 3):** overlap heuristics ≈ 55% degree / 45% structure; structure residual CN +0.192, AA +0.194, Jaccard +0.196 (all permutation p<0.001); PrefAttach −0.002 (p=1.000), Random −0.000 (p=1.000).
**Full-ranking (Table S1, TransE, pool ~14.3k):** R0 0.079 → R1 0.071 (−10%) → R2 0.023 (**−71%**). 30 cold-start genes → worst rank.
**Hetionet (Table 4):** CN −35%, AA −37%, Jaccard −36% at R2; PrefAttach +0%; Random +0%.
**AUROC–MRR dissociation:** AUROC stays ~0.96–0.99 (type-matched ~0.88 at R2) while MRR collapses → AUROC hides the leak.

**KGE hyperparameters:** dim 64, 300 epochs, NSSA loss, sLCWA training loop, batch 16,384, 16 negatives, PyKEEN. Refs: Bordes 2013 (TransE), Sun 2019 (RotatE + NSSA).

---

## §Front matter

### Title options
- Primary (matches metadata): *Quantifying evaluation leakage in multi-species gene–disease link prediction: a graded audit and a de-leaked benchmark.*
- Alt (method-class angle): *Method-class-dependent evaluation leakage in gene–disease link prediction: degree, not orthology, drives inflated scores.*

### Abstract building blocks (assemble ~250 words)
1. **Motivation:** KGE models report high MRR/AUROC for gene–disease association (GDA) prediction; these numbers guide biological prioritization but are inflated by evaluation leakage of several distinct kinds.
2. **What we did:** a graded audit — five topological baselines (Random, Common-Neighbors, Adamic-Adar, Jaccard, Preferential-Attachment) and two KGE models (TransE, RotatE) through one shared ranking harness on 4,228 held-out human gene–disease edges, across four regimes (R0 standard, R1 redundancy, R2 degree-preserving null, R3 orthology-blocked), 3 seeds, on a 452k-node Monarch-derived graph, replicated on Hetionet.
3. **Findings:** leakage is method-class-dependent — degree structure (R2) is the one universal leak (baselines −43 to −49% MRR, KGE −21 to −24%; strengthening to −71% under full ranking); redundancy (R1) leaks only into multi-hop KGE (−6 to −8.5%); orthology (R3) is a null (±0.5%). PreferentialAttachment (pure degree) is unchanged by the null, validating it.
4. **Rigor + release:** confidence intervals, ≥3 seeds, significance tests, a formal degree-preserving null, and cross-graph replication — reported by almost no comparable paper. We release the de-leaked splits, result files, and one-command reproduction (DOI `10.5281/zenodo.21398832`).

### Keywords (from `.zenodo.json`)
knowledge graph · knowledge graph embedding · link prediction · gene–disease association · data leakage · benchmark evaluation · node degree bias · reproducibility · Monarch Initiative · TransE · RotatE

### Authorship / contributions (draft)
- Authors: Owen Liu (ORCID 0009-0003-0476-5270), Ryan Tang, Daniel Augustine, + **mentor co-author (TBD — recruit)**.
- Contributions template: O.L. conceived the audit, built the pipeline, ran experiments, wrote the draft; R.T. / D.A. contributed to [data integration / analysis / figures]; mentor supervised and revised.
- Competing interests: none declared. Funding: none / science-fair context.

---

## §1. Introduction — talking points & claim→citation map

- **Opening claim:** GDA link prediction on biomedical KGs routinely reports MRR > 0.8 and AUROC > 0.95; such scores feed disease-gene prioritization. → cite an application/resource example (TarKG 2024 [lit-doc]; Nunes 2023 [lit-doc]).
- **The leakage problem is general and under-appreciated in biology.** → Kapoor & Narayanan 2023 (ML leakage reproducibility crisis) [ADD]; Bernett 2024a (guiding questions to avoid data leakage in biological ML) [lit-doc].
- **Several distinct mechanisms, usually conflated:** train/test redundancy (Akrami 2020; Toutanova & Chen 2015 [lit-doc]); node-degree / popularity bias (Yılmaz 2025; Aiyappa 2024; Shomer 2023 [lit-doc]); non-representative test sets (Brière 2025 DL3 [lit-doc]).
- **The gap:** no graded, single-harness decomposition that (a) separates these mechanisms on one graph, (b) measures the *magnitude* each contributes via a formal degree-preserving null, and (c) releases the de-leaked splits.
- **Contributions (enumerate):**
  1. A graded regime framework (R0–R3) run through one shared harness for both topological and KGE methods.
  2. A degree-preserving null that isolates and *quantifies* degree vs. genuine structure (not a binary retrain).
  3. The method-class-dependent result, with orthology and redundancy as designed negative controls.
  4. A reproducible, DOI-archived de-leaked benchmark, replicated on a second independent graph.
- **Honesty sentence (write it early):** "The core phenomenon — that random-split evaluation is inflated by node degree and leakage — is established; our contribution is a controlled, quantitative consolidation across mechanisms and graphs, with a released resource."

---

## §2. Related work — the competitive landscape (pull heavily from lit-doc)

### The three DIRECT competitors (lead with these; out-position, don't bury)

**A. Yılmaz, Yorgancioglu & Koyutürk 2025 — PNAS 122(24) e2416646122 [lit-doc, read-in-full].**
- Showed random-edge-sampling evaluation is biased toward high-degree ("rich") nodes (Matthew effect), persisting across temporal snapshots. PPI networks (Biogrid + STRING). Preferential Attachment as the pure-bias reference. Proposed AWARE (weighted/stratified low-degree eval).
- Headline: "a simple predictive model that uses bias as the only feature outperforms sophisticated models."
- Stakes quote: "95% of all life science publications focus on ~5,000 well-studied proteins."
- **How we differ:** heterogeneous multi-relational KG (not PPI); KGE models included; permutation-null *magnitude* decomposition; the orthology regime. **Our sibling in a top venue — cite respectfully, claim the KG/KGE + decomposition + orthology daylight.**

**B. Bradshaw & Layer 2024 — bioRxiv 10.1101/2024.06.10.598277, CU Boulder [lit-doc, read-in-full]. ⚠️ NAME CONFLICT — see §Citations.**
- KGE link prediction for rare-disease gene→disease **on the Monarch KG (our graph)**. Filtering Monarch to ~11% of size *improves* performance ("quality > quantity"); names KGE degree bias; coins the "bias cycle."
- **How we differ:** same graph/task — MUST cite. We add the formal 4-regime audit + degree-null + second graph + topological baselines/anchors. (Tension with Nunes 2023, who found richer graphs help.)

**C. Brière & Baudot 2025 — bioRxiv 10.1101/2025.01.23.634511, Aix-Marseille/INSERM/CNRS [lit-doc, read-in-full].**
- Systematic framework for **three leakage types** in KGE link prediction on **ShepherdKG** (PrimeKG + Orphanet; ~5.5M triples), drug→disease (repurposing), 7 KGE models, filtered MRR only.
- DL1 = train/test separation (remove near-dup/near-reverse/Cartesian relations, per Akrami 2020) ≈ **our R1**. DL2 = degree-permutation retrain ≈ **our null (but binary)**. DL3 = non-representative test set.
- Key numbers: DL1 inflates MRR most for ComplEx (+24.9%), ANALOGY (+19.5%), HolE (+11.4%); TransE/TransD/DistMult/RESCAL insensitive. **DL2: found NO evidence of degree reliance.** DL3: test MRR ~0.2 vs real-inference MRR ~0.01.
- **How we differ + THE reconciliation:** topological baselines + Random/PA anchors; a degree-null that *decomposes magnitude* (not binary retrain); CIs/seeds/significance; two graphs; gene→disease. **Their DL2 "no degree" must be reconciled with our R2 "degree is THE leak" — see §6.**

### Method / foundation lineage (Cluster 3 in lit-doc)
- **Akrami 2020** — Realistic Re-evaluation of KG Completion (SIGMOD). Removing near-duplicate/near-reverse/Cartesian relations sharply drops KGE scores on FB15k/WN18/YAGO3-10. **Origin of our R1 and Brière's DL1.** [lit-doc]
- **Toutanova & Chen 2015** — created FB15k-237 by removing reversible relations. ML precedent for R1. [lit-doc]
- **Bernett 2024a** — Nature Methods 21:1444, three-source leakage taxonomy. [lit-doc]
- **Bernett 2024b** — Brief. Bioinform. 25(2) bbae076; sequence-based PPI models exploit degree, exposed by degree-preserving permutation. **Ancestor of our null.** [lit-doc]
- **Maslov & Sneppen 2002** — Science 296:910, degree-preserving edge randomization (configuration/edge-swap). **The primitive behind R2 — cite this for the null.** [ADD — currently missing]

### Degree / popularity-bias background
- **Shomer 2023** — Toward Degree Bias in Embedding-Based KG Completion (WWW). KGE learns poor low-degree reps; KG-Mixup. [lit-doc]
- **Aiyappa 2024** — Implicit degree bias in link prediction (arXiv 2405.14985); a degree-only baseline performs near-optimally; proposes a degree-corrected task. [lit-doc]

### Application / convention examples (show you fit the genre)
- **Gualdi 2024** (DLemb/BioKG2Vec) — new embeddings; validated by predicting IDD genes + enrichment; compared vs a null. **Template for the case study we lack.** [lit-doc]
- **Nunes 2023** — multi-ontology KGE; richer interconnected KGs improve GDA (counterpoint to Bradshaw). [lit-doc]
- **TarKG 2024** — resource KG + 6 KGE + Alzheimer's top-10 case study; no stats. [lit-doc, summary]
- **Diabetes gene discovery 2021** (Front. Genet. 12:779186) — embedding+classifier + enrichment; no stats. [lit-doc, summary]
- **HGNN benchmark 2025** (Bioinformatics Advances vbaf187) — 9 GNNs vs 33 baselines, 8 networks; no stats. Closest to our benchmark structure. [lit-doc, summary]
- **Systematic KGE for GDA 2025** (arXiv 2504.08445) — link-pred vs node-pair classification on DisGeNET + ontologies; Hits@k only; no stats. [lit-doc, summary]

### Related-work differentiation table (MASTER Table D — drop into §2 or §6)
| Competitor | What they did | What we add |
|---|---|---|
| Brière 2025 | 3 leakage types, KGE-only, MRR-only, single runs, 1 KG (drug repurposing) | topo+KGE one harness; degree-null *magnitude* decomposition; CIs/seeds/sig; two graphs; gene→disease |
| Yılmaz PNAS 2025 | Degree bias in PPI eval; PA-as-bias; stratified metrics; no KG/KGE | heterogeneous KG + KGE; permutation null; orthology regime |
| Bradshaw 2024 | Degree bias + "bias cycle" on Monarch; filtering helps | formal 4-regime audit + null + second graph + topo baselines/anchors |
| Gualdi 2024 | New algorithms + single-disease validation | we audit *existing* methods across regimes/graphs (different goal) |

---

## §3. Materials — the knowledge graph

### Composition & honest scope
- **~98.5% Monarch Initiative content** (Putman 2024, NAR 52(D1):D938) [lit-doc] + 7 curated non-Monarch sources: DGIdb, GWAS Catalog, DrugCentral, Orphadata, Gene2Phenotype, CIViC, DIDA.
- **Scope discipline (from README, keep it):** Monarch is itself an aggregator; describe as a *Monarch-derived* graph, NOT an equal integration of many independent databases. Do not read "supported by N datasets" as N independent votes. Per-source contributions: `data/processed/provenance/`.
- Licensing: code MIT; data CC BY 4.0 with attribution to Monarch + primary sources (`LICENSE`, `LICENSE-DATA`).

### Table 1 (verbatim source: `tables/table1_graph_stats.md`)
| Property | Monarch (this work) | Hetionet v1.0 (robustness) |
|---|---|---|
| Nodes | 452,012 | 47,031 |
| Edges (directed) | 5,860,539 | 2,250,197 |
| Edges (unique undirected) | 5,583,381 | 2,107,709 |
| Node categories / metanodes | 8 | 11 |
| Relation / metaedge types | 28 | 24 |
| Degree mean | 25.9 | 99.7 |
| Degree median | 2 | 32 |
| Degree max | 20,662 | 23,879 |
| Gene nodes | 336,460 | 19,145 |
| Disease nodes | 14,339 | 136 |
| Gene–disease target edges | 42,288 | 12,623 |
| Orthology edges | 510,703 | 0 (single-species) |

- **Key contrast to narrate:** Monarch has degree median 2 (very heavy-tailed, many singletons) vs Hetionet median 32 — Monarch is *more* degree-dominated, which is why the degree collapse is larger there.
- Second graph provenance: Hetionet (Himmelstein 2017, eLife 6:e26726) [lit-doc]. PrimeKG (Chandak 2023) is Brière's source, for comparison [lit-doc].
- **Orthology basis for R3** (cite the orthology data source — Monarch draws cross-species orthology; verify whether via PANTHER / Alliance of Genome Resources and cite it).

### ⚠️ Figures to generate (currently gaps, Med priority)
- Schema/metagraph figure (node categories + relation types).
- Degree-distribution figure (log-log; shows the heavy tail that drives everything).

---

## §4. Methods — full detail

### 4.1 The graded leakage regimes
- **R0 standard:** ordinary random train/valid/test split. The baseline everyone reports.
- **R1 redundancy:** remove **1,661 direct held-pair edges** (train copies of test facts) — 0.03% of edges. Built at `scripts/build_deleaked_splits.py:215`. NOT a training-size manipulation.
- **R2 degree-preserving null:** rewire the graph to destroy genuine structure while holding the degree sequence *exactly*. Built by `scripts/build_degree_null.py`. Validated because PreferentialAttachment (pure degree) is unchanged.
- **R3 orthology-blocked:** remove cross-species `ORTHOLOGOUS_TO` bridges (**520,198 edges**: R0 train 5,852,083 → R3 train 5,331,885). The domain-specific hypothesis a biomedical reviewer would raise.

### 4.2 The degree-preserving null (centerpiece method)
- Configuration-model / edge-swap randomization preserving each node's degree (Maslov–Sneppen 2002 [ADD]; kin to Bernett 2024b's PPI permutation).
- Decomposition logic: score on real R0 vs the degree-null gives *structure = real − null*; a positive, significant residual = genuine signal beyond degree. PrefAttach residual ≈ 0 by construction = the built-in positive control.
- 10 permutation replicates; R0 evaluation held fixed; permutation p-values.

### 4.3 Models & baselines
- **KGE:** TransE (Bordes 2013) + RotatE (Sun 2019, with NSSA loss). dim 64, 300 epochs, sLCWA, batch 16,384, 16 negatives, seeds 42/1/7. Toolkit: PyKEEN (Ali 2021) [lit-doc].
- **ComplEx / DistMult status (IMPORTANT, be honest):** ComplEx (Trouillon 2016 [ADD]) was attempted repeatedly (smoke + 3 "rescue" runs) and is **untrainable/collapses at 452k entities** — only smoke-test JSONs exist (`kge_smoke*`, `kge_ComplEx_R0_standard_seed42_d128_e5.json`). DistMult (Yang 2015 [ADD]) has a single d128/e5 smoke test only. **Do not claim full ComplEx/DistMult runs.** If a reviewer asks for more KGE families, DistMult is the feasible add; ComplEx is a documented dead end at this scale.
- **Topological baselines:** Random, Common-Neighbors, Adamic-Adar, Jaccard, Preferential-Attachment (Liben-Nowell & Kleinberg 2007; PA from Barabási & Albert 1999) [lit-doc]. PrefAttach doubles as the degree positive-control.

### 4.4 Shared evaluation harness
- `scripts/lib_eval.py` ranks *every* method identically (target relation `BIOLINK:GENE_ASSOCIATED_WITH_CONDITION`). Single-harness fairness is a core design claim — topo and KGE are never scored by different code.

### 4.5 Ranking protocols (name which every number is from)
- **Sampled-50 negatives** (main, Table 2): rank the true disease against 50 type-matched sampled negatives; `rank_seed=42`.
- **Filtered full ranking** (robustness, Table S1): rank the true disease against the entire ~14,307-disease pool (14,252 for R1). Fast path verified exactly against brute force.
- **Cold-start:** 30 test genes absent from the training vocabulary receive the worst possible rank (stated explicitly).
- **Protocol sensitivity to disclose:** "PrefAttach = strongest single baseline" is true by sampled-50 MRR and by full-ranking median/Hits@100, but FALSE by full-ranking MRR (there PrefAttach is weakest, AdamicAdar strongest). Always name the metric — the MRR-vs-median dissociation is itself a reportable sub-result.

### 4.6 Statistics (the paper's biggest edge — foreground)
- Bootstrap **confidence intervals** on every metric.
- **≥3 seeds** (42, 1, 7), mean ± SD reported.
- **Significance:** paired bootstrap MRR differences, 1,000 resamples, n=4,228 (Table 2 note).
- **Null model:** the formal degree-preserving permutation with p-values.
- Almost no comparable paper reports all four (see §Field patterns) — say so.

---

## §5. Results — every number, subsection by subsection

### 5.1 The graded audit → Table 2 + Figure 2
Source: `tables/table2_audit.md` (mean ± SD over seeds 42/1/7; 50 type-matched negatives).

| Method | MRR R0 | MRR R1 | MRR R2 | MRR R3 | ΔR0→R2 | AUROC R0 | Hits@10 R0 |
|---|---|---|---|---|---|---|---|
| Random | 0.091±.003 | 0.091±.003 | 0.091±.003 | 0.091±.003 | +0% | 0.503±.001 | 0.199±.002 |
| CommonNeighbors | 0.441±.001 | 0.442±.001 | 0.247±.001 | 0.443±.002 | −44% | 0.736 | 0.522 |
| AdamicAdar | 0.459 | 0.459±.002 | 0.262 | 0.460±.003 | −43% | 0.738 | 0.528 |
| Jaccard | 0.397±.001 | 0.398±.002 | 0.201±.001 | 0.399±.001 | −49% | 0.690±.001 | 0.474±.002 |
| PreferentialAttachment | 0.587±.002 | 0.582±.001 | 0.589±.002 | 0.581±.002 | +0% | 0.794 | 0.774±.001 |
| TransE | 0.799 | 0.750±.004 | 0.611±.004 | 0.800±.002 | −24% | 0.962±.001 | 0.975±.001 |
| RotatE | 0.827±.007 | 0.757±.004 | 0.655±.004 | 0.830±.004 | −21% | 0.967±.001 | 0.974 |

Paired-bootstrap (AdamicAdar − method; R0, seed 42, 1,000 resamples, n=4,228): CommonNeighbors Δ=+0.017 (p<0.001); Jaccard Δ=+0.060 (p<0.001); PreferentialAttachment Δ=−0.130 (p<0.001); Random Δ=+0.365 (p<0.001). All significant at p<0.001 (two-sided).

**Sentence:** "Removing genuine structure while preserving degree (R2) collapses every overlap heuristic by 43–49% and both KGE models by 21–24%, but leaves Preferential-Attachment — which uses degree alone — unchanged, isolating node degree as the single leakage source common to all method classes."

### 5.2 Degree is the universal leak → Table 3 + Figure 3
Source: `tables/table3_degree_null.md` (10 permutation replicates; R0 held fixed).
PrefAttach sanity: real 0.589 vs null 0.591 (unchanged → the null preserves the degree sequence).

| Method | Real R0 MRR | Degree-null MRR | Structure (real−null) | Perm. p | Structure>degree? |
|---|---|---|---|---|---|
| CommonNeighbors | 0.441 | 0.249±.003 | +0.192 | <0.001 | yes |
| AdamicAdar | 0.458 | 0.264±.004 | +0.194 | <0.001 | yes |
| Jaccard | 0.398 | 0.202±.003 | +0.196 | <0.001 | yes |
| PreferentialAttachment | 0.589 | 0.591±.000 | −0.002 | 1.000 | no |
| Random | 0.094 | 0.094±.000 | −0.000 | 1.000 | no |

**Sentence:** "For the overlap heuristics roughly 55% of R0 performance is attributable to degree and ~45% to genuine structure (permutation p<0.001); for Preferential-Attachment the structure residual is statistically zero, confirming it is a pure-degree predictor."

### 5.3 Redundancy is a multi-hop-only leak (R1)
- KGE: TransE −6.2%, RotatE −8.5% (sampled-50). Full-rank TransE R1 −10%.
- 2-hop baselines: CN/AA/Jaccard flat (±0.3%).
- **Mechanism:** a direct gene→disease edge is embedded as *proximity* by KGE (it shortens the translation/rotation distance) but is NOT a shared neighbor, so it never enters CN/AA/Jaccard. Confirmed at `build_deleaked_splits.py:215`.
- **Pre-empt "R1 is just less training data":** only 1,661 edges = 0.03% of the graph; the drop is structural, not a data-volume artifact.
- **Reverse-of-pre-registration note (honesty point):** the pre-registered guess was that *orthology* would be the multi-hop leak; it turned out to be *redundancy*. Report this honestly — it strengthens credibility.

### 5.4 Orthology is a null (R3)
- Every method within ±0.5% of R0 (KGE +0.2–0.4%). Removing ~520K orthology edges does NOT reduce performance.
- **Framing:** the designed negative control / refuted domain hypothesis. This is the biomedical reviewer's intuitive culprit, and the data refutes it. A *clean null is a result*, not a failure.

### 5.5 Full-ranking robustness → Table S1
Source: `tables/tableS1_fullrank_robustness.md` (TransE, R0/R1/R2, seeds 1/7/42, filtered full ranking vs whole disease pool).

| Model | Regime | Pool D | MRR | Hits@1 | Hits@3 | Hits@10 | AUROC(type) | ΔMRR vs R0 |
|---|---|---|---|---|---|---|---|---|
| TransE | R0 | 14,307 | 0.079±.001 | 0.036 | 0.072 | 0.150 | 0.962 | ref |
| TransE | R1 | 14,252 | 0.071±.003 | 0.031 | 0.066 | 0.138 | 0.954 | −10% |
| TransE | R2 | 14,307 | 0.023±.002 | 0.005 | 0.014 | 0.042 | 0.879 | **−71%** |

- Baseline full-ranking (R0+R2) also exists (`data/processed/results/robustness/robustness_R{0,2}.json`): overlap heuristics drop far more under full ranking (median rank → pool midpoint ≈ chance); PrefAttach unchanged.
- **Sentence:** "Under the stricter filtered full-ranking protocol the degree collapse does not merely persist — it strengthens, from −24% to −71% MRR for TransE — showing the sampled-negative protocol *understates* degree dominance."
- **Justify TransE-only:** TransE is representative of the KGE class (tracks RotatE within a few points throughout Table 2); full-ranking RotatE and R3 were deliberately not run because R3 is already a null and the point is the R2 collapse. Flag RotatE full-rank as trivial future work if a reviewer insists.

### 5.6 Independent-graph replication → Table 4
Source: `tables/table4_hetionet.md` (R3 omitted — Hetionet is single-species, no orthology).

| Method | MRR R0 | MRR R1 | MRR R2 | ΔR0→R2 |
|---|---|---|---|---|
| Random | 0.092±.001 | 0.092±.001 | 0.092±.001 | +0% |
| CommonNeighbors | 0.269±.003 | 0.268±.002 | 0.173±.001 | −35% |
| AdamicAdar | 0.285±.005 | 0.284±.004 | 0.178±.002 | −37% |
| Jaccard | 0.162±.002 | 0.162±.001 | 0.104±.001 | −36% |
| PreferentialAttachment | 0.231±.002 | 0.231±.002 | 0.231±.002 | −0% |

- **Sentence:** "The same degree collapse (−35 to −37% at R2, Preferential-Attachment unchanged) reproduces on Hetionet, an independently constructed graph with no shared pipeline, establishing that the effect is a property of biomedical KG topology rather than of Monarch or of our build."
- Monarch collapse (−43 to −49%) > Hetionet (−35 to −37%): Monarch is more degree-dominated (median degree 2 vs 32).

### 5.7 AUROC–MRR dissociation → Figure 1
- AUROC (random negatives) stays ~0.99 across all regimes; AUROC (type-matched) ~0.96 at R0, ~0.88 at R2 — i.e. even where MRR collapses by tens of percent, AUROC barely moves.
- **Sentence:** "Because AUROC against easy negatives saturates near 1.0 regardless of regime, it is blind to the leakage that ranking metrics expose; benchmarks that report only AUROC will not detect degree inflation." Cite Saito & Rehmsmeier 2015 (PR vs ROC under imbalance) [ADD].

---

## §6. Discussion — arguments to make

### Implications for practice
- Report an **R2 degree-null sanity check** as standard in GDA benchmarks.
- Include **degree baselines (PrefAttach)** — if your KGE barely beats pure degree, say so (RotatE 0.655 vs PrefAttach 0.589 at R2: KGE keeps only a small increment above degree-alone).
- **Do not report AUROC alone** — pair with ranking metrics and hard/type-matched negatives.
- Name the ranking protocol on every number (sampled vs full).

### THE reconciliation: our R2 vs Brière DL2
- Brière DL2 (degree-permute + retrain) found *no* degree reliance; our R2 finds degree is *the* universal leak. Reconcile via:
  1. **Binary vs magnitude:** DL2 asks "does removing degree destroy performance?" (a yes/no retrain); R2 *decomposes* how much of the score is degree vs structure. A model can retain rank-order while still owing ~half its score to degree.
  2. **Task/graph differences:** DL2 is drug→disease on ShepherdKG/PrimeKG; ours is gene→disease on Monarch/Hetionet, which are more degree-dominated (median degree 2).
  3. **Their DL3 supports us:** Brière's own non-representative-test result (test MRR ~0.2 vs real-inference ~0.01; well-studied drugs predicted far better) is exactly a study/degree-bias signature.
- **Net message:** the two papers are consistent once you separate "is degree used at all" from "how much of the number is degree."

### Field-pattern framing (from lit-doc Table C — 9 comparison papers)
- ~1/9 report CIs; ~2/9 multiple seeds; ~2/9 significance; ~3/9 a null model → **our rigor is the edge.**
- ~6/9 include a biological case study → **the main thing they have that we don't.**
- 3/9 address degree/leakage (PNAS, Brière, Bradshaw) → **our topic and our competition.**

### Recommendations checklist (put as a boxed list)
1. Publish a degree-preserving null result alongside headline numbers.
2. Report PrefAttach (and Random) anchors.
3. Report ranking metrics + hard negatives, not AUROC alone.
4. State the ranking protocol explicitly.
5. Release de-leaked splits + seeds + hashes for reproduction.

---

## §7. Limitations (be generous — this is a soundness paper)
- **No genuine temporal split** — the graph carries no edge-discovery dates; regimes control optimism on a random split, not time.
- **Monarch-dominated coverage** — ~98.5% Monarch even after 7 more sources; not an independent multi-database integration.
- **Incremental novelty** — the phenomenon is established (PNAS 2025 / Brière 2025 / Bradshaw 2024); contribution is consolidation + decomposition + release.
- **KGE full-ranking is TransE-only** — representative of the class; RotatE full-rank is straightforward future work.
- **Two KGE families only** — ComplEx is untrainable at this scale (documented); DistMult only smoke-tested.
- **Protocol sensitivity** — sampled-50 vs full ranking can reorder methods (the PrefAttach MRR reversal); we report both.
- **Acknowledged gaps (state as future work or add before a benchmark venue):** no biological/literature case study; no low-degree/stratified (held-out-disease, cold-start) evaluation.

---

## §8. Data & code availability (draft text)
"All code (MIT) and de-leaked data (CC BY 4.0) are archived at Zenodo (DOI 10.5281/zenodo.21398832) and on GitHub (github.com/owen-hx-liu/gene-disease-leakage-benchmark). The release includes the de-leaked splits with manifest (SHA-256 bc468864…), per-regime baseline/KGE/robustness/Hetionet result JSONs, table and figure generators that read only the frozen result files, and a one-command reproduction script (`run_all.sh`). All runs use fixed seeds (42, 1, 7)."

---

## §Citations — master list + how to cite + what to ADD

### How to build the bibliography (workflow)
1. **Use a reference manager** (Zotero, free). For each DOI/arXiv id, "Add by identifier" auto-pulls metadata; export BibTeX or CSL-JSON.
2. **Cite at the point of claim** — every number/claim in the master file above has its citation named inline; carry those into the draft.
3. **Match the venue's style** — JEI has its own reference format (check its author guide); *Database (Oxford)* / *NAR GB* use numbered Vancouver. Set the CSL style in Zotero once and it reformats all.
4. **Update preprints to published versions** — several 2024/2025 bioRxiv/arXiv items may now be peer-reviewed; cite the journal version if it exists (check the DOI landing page). Especially Brière, Bradshaw, the arXiv items.
5. **Verify before submission** — confirm author list, year, venue, and DOI for every reference you did not read in full (the 🟨/🔵 items in your lit-doc).

### ⚠️ Conflicts / things to resolve
- **Gu et al. vs Bradshaw & Layer (Monarch degree-bias paper).** My prior project note referred to the closest competitor as "Gu et al. 2024"; your lit-doc calls it "Bradshaw & Layer 2024" (CU Boulder, bioRxiv 2024.06.10.598277). These may be the same paper (mis-remembered on my side) or two different ones. **Resolve before writing Related Work — getting the lead author wrong on your nearest competitor is a red flag to reviewers.** Action: open the bioRxiv DOI, confirm the author list, and standardize.

### Already in your lit-doc (23 — keep, verify details)
Direct competitors: Brière & Baudot 2025; Yılmaz et al. 2025 (PNAS); Bradshaw & Layer 2024. Applications/conventions: Gualdi 2024; Nunes 2023; TarKG 2024; Diabetes-genes 2021; HGNN benchmark 2025; Systematic-KGE 2025. Foundations: Akrami 2020; Bernett 2024a; Bernett 2024b; Toutanova & Chen 2015. Degree bias: Shomer 2023; Aiyappa 2024. Data/models: Putman 2024 (Monarch); Himmelstein 2017 (Hetionet); Chandak 2023 (PrimeKG); Bordes 2013 (TransE); Sun 2019 (RotatE); Ali 2021 (PyKEEN); Liben-Nowell & Kleinberg 2007 (topo baselines); Barabási & Albert 1999 (PA).

### RECOMMENDED ADDITIONS (not in your 23) — with why & priority

| # | Citation (verify exact details) | Why you need it | Priority |
|---|---|---|---|
| A1 | **Maslov & Sneppen 2002**, Science 296:910 — degree-preserving edge randomization | The methodological *primitive* behind your R2 null; right now R2's ancestry is only cited via a PPI paper. Cite the source. | **High** |
| A2 | **Kapoor & Narayanan 2023**, Patterns 4(9):100804 — "Leakage and the reproducibility crisis in ML-based science" | High-profile, general framing for the whole leakage motivation in §1. | **High** |
| A3 | **Trouillon et al. 2016**, ICML — ComplEx | You discuss ComplEx (untrainable / vs Brière); cite the model paper. | **High** |
| A4 | **Yang et al. 2015**, ICLR — DistMult | Same; needed if you mention/add DistMult. | Med |
| A5 | **Saito & Rehmsmeier 2015**, PLoS ONE 10(3):e0118432 — PR vs ROC under class imbalance | Backs the §5.7 AUROC-hides-leakage argument. | Med |
| A6 | **Sun et al. 2020**, ICLR — "A Re-evaluation of KG Completion Methods" (scoring-tie test leakage) | Reinforces the evaluation-protocol critique alongside Akrami. | Med |
| A7 | **Hu et al. 2020**, NeurIPS — Open Graph Benchmark (OGB) | Standard link-prediction evaluation conventions (filtered ranking, negatives). | Low |
| A8 | Orthology source (e.g. **Alliance of Genome Resources 2024** and/or **PANTHER, Mi et al.**) — VERIFY what Monarch uses | Grounds the R3 orthology regime in its data provenance. | Med |
| A9 | **Rossi et al. 2021**, ACM TKDD — KGE comparative analysis | Optional benchmark-conventions reference. | Low |
| A10 | A GDA/DisGeNET reference (**Piñero et al. 2020**, NAR) if you reference DisGeNET | Convention; only if used. | Low |

---

## §Open issues / TODO before the paper is submission-ready
1. **Resolve the Gu/Bradshaw citation conflict** (above).
2. **Recruit a mentor/senior co-author** (required for JEI; big odds-mover generally).
3. **Post bioRxiv preprint**, then cross-link its DOI into `.zenodo.json` `related_identifiers`.
4. **Decide the venue path:** JEI now (paper is complete for it) vs. benchmark venue (add the two gaps first).
5. **Optional experiments for a benchmark venue:** (a) rare-disease top-k + enrichment case study (template: Gualdi 2024); (b) low-degree / held-out-disease / cold-start stratified eval (differentiates from PNAS & Brière); (c) DistMult as a 3rd KGE family (feasible).
6. **Generate the two missing figures:** schema/metagraph + degree distribution.
7. **Write the Brière-DL2 reconciliation** paragraph (§6) — a likely reviewer question.
8. **Verify all 🟨/🔵 citations** and update preprints to published versions.
