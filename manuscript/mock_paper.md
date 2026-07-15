# ⚠️ MOCK PAPER — ALL NUMBERS ARE FAKE PLACEHOLDERS ⚠️

*This is a dummy manuscript to show what the final paper would look like. Every number in
`[[double brackets]]` is INVENTED and must be replaced with real computed results. Do not
submit, cite, or share this — it is a layout mockup only. The Related Work section cites REAL
papers (found in the 2026-07 literature search) so the honest positioning is baked in.*

---

# Cross-species orthology inflates gene–disease link-prediction benchmarks on multi-species knowledge graphs

**Authors:** [[Your Name]]¹, [[Mentor/Co-author Name]]²
¹[[School / affiliation]] · ²[[Institution]]
**Corresponding author:** [[email]]

---

## Abstract

Knowledge-graph embedding (KGE) models are widely used to predict gene–disease associations
and are typically reported with high area-under-the-curve scores. Recent work has shown that
such benchmarks are vulnerable to data leakage from train/test redundancy and node-degree bias
(Ranga et al., 2025; Gu et al., 2024). We identify and quantify an additional, biology-specific
leakage source that existing leakage taxonomies do not address: **cross-species orthology**. On
a multi-species knowledge graph derived from the Monarch Initiative (448,931 nodes; 5,860,539
edges), human gene→disease facts can be recovered through a held-out gene's model-organism
ortholog and its phenotype profile, without independent human evidence. We introduce a graded
evaluation protocol (R0–R3) that removes hub, orthology, and tautological shortcuts in turn, and
we benchmark four KGE models against topological baselines under each regime. Under the standard
random split, Adamic–Adar reaches MRR [[0.69]] versus [[0.09]] for the best KGE model; after
removing orthology bridges, all methods fall to MRR [[0.xx]]–[[0.xx]], a relative drop of
[[XX%]], and the method ranking [[changes/holds]]. Crucially, [[XX%]] of the orthology effect
[[remains / disappears]] after controlling for node degree, indicating orthology leakage is
[[a distinct confound / largely explained by degree]]. We further validate a leakage-robust
predictor on held-out rare-disease associations, recovering [[N]] links absent from the training
release. Code, data, and the de-leaked splits are released for reproducibility.

**Keywords:** knowledge graph, link prediction, gene–disease association, data leakage,
orthology, benchmark evaluation

---

## 1. Introduction

Identifying which genes underlie which diseases is a central problem in genomics, and exhaustive
experimental testing of candidate gene–disease pairs is infeasible. Computational prioritization
from integrated biomedical knowledge graphs (KGs) is therefore widely used, most prominently via
knowledge-graph embedding (KGE) models such as TransE, DistMult, ComplEx, and RotatE.

A recurring concern is that reported performance may be inflated by **data leakage** — situations
where the test set is recoverable from the training data through trivial or non-biological
shortcuts rather than genuine inference. In the general KG literature, inverse-relation leakage
motivated the FB15k-237 and WN18RR benchmarks (Toutanova & Chen, 2015; Dettmers et al., 2018).
In biomedicine, OpenBioLink (Breit et al., 2020) removed reverse/symmetric-relation leakage;
more recently, Ranga et al. (2025) benchmarked leakage from train/test redundancy, node degree,
and cold-start conditions in drug-repurposing prediction, and Gu et al. (2024) showed that KGE
models on the Monarch KG are strongly degree-biased and that aggressive graph filtering improves
performance.

We observe that a leakage source specific to **multi-species** biomedical KGs has not been
isolated or quantified: **cross-species orthology**. Aggregated resources such as Monarch
deliberately propagate phenotype and disease information across species via orthology to increase
coverage (from ~51% to ~89% of human genes; Mungall et al., 2017). This is a genuine biological
asset for discovery. However, for *evaluation*, it creates a shortcut: a held-out human
`gene —associated_with→ disease` edge can be scored highly through the gene's mouse or zebrafish
ortholog and that ortholog's phenotype annotations, independent of any human evidence. Under a
naive random split, this inflates apparent performance and may misrepresent a method's ability to
make genuinely novel *human* predictions.

**Contributions.** (C1) We characterize cross-species orthology as a leakage source distinct from
the degree and redundancy leakage in existing taxonomies. (C2) We release a graded, de-leaked
evaluation protocol (R0–R3) and an open harness. (C3) We quantify how much reported gene–disease
ranking performance is attributable to orthology, whether it survives degree control, and whether
it changes the KGE-versus-baseline verdict. (C4) We introduce a leakage-robust hybrid predictor.
(C5) We validate it on held-out rare-disease associations. We make no novelty claim about the
graph itself, which is Monarch-derived.

---

## 2. Related work

**KG-embedding link prediction and leakage.** Inverse-relation leakage in FB15k/WN18 (Toutanova
& Chen, 2015) motivated the widely used FB15k-237/WN18RR benchmarks. In biomedicine, OpenBioLink
(Breit et al., 2020) explicitly filtered reverse and symmetric-relation leakage and restricted
test entities to those seen in training. Ranga et al. (2025) provide the most directly comparable
recent work, a framework for detecting and controlling three leakage types — train/test
redundancy, illegitimate features (node degree), and cold-start — for **drug–disease** repurposing
prediction. Our work is complementary: we target **gene–disease** prediction and a **fourth**
leakage source (orthology) that their taxonomy does not include.

**Topology and degree bias.** Gu et al. (2024) analyze the Monarch KG and show KGE models are
degree-biased, introducing a permutation framework that attributes much of measured performance
to degree alone, and demonstrating that filtering the KG improves prediction. Because orthology
hubs are high-degree, we explicitly separate the orthology effect from degree (Section 5.3) so as
not to re-report a known result. Topological imbalance in biomedical KGs has also been studied by
[[Bonner et al., 2021]].

**Topological baselines.** Simple neighborhood heuristics (common-neighbors, Adamic–Adar, Jaccard)
are known to be competitive with embedding methods on some biomedical graphs (Crichton et al.,
2018), though the strongest baseline is graph-dependent; we therefore report the full baseline set
rather than a single heuristic.

**Multi-species inference.** Cross-species phenotype transfer (phenologs; McGary et al., 2010) and
ortholog-based disease-gene inference are established, and Monarch is explicitly built to exploit
them. We treat this signal not as invalid but as a source of *evaluation* leakage that must be
controlled to obtain honest estimates of novel-prediction performance.

---

## 3. Data

We use the integrated Monarch-derived knowledge graph (`edges_clean_integrated.csv`): 448,931
nodes and 5,860,539 directed edges (5,494,276 unique undirected; 0 self-loops after cleaning),
across 28 Biolink relation types. Approximately 98.4% of edges derive from Monarch, augmented with
seven curated sources contributing ~1.6% (Table 1). The prediction target is the human
`gene —GENE_ASSOCIATED_WITH_CONDITION→ disease` relation (42,288 edges; all gene endpoints are
human HGNC identifiers). The graph contains 510,703 `ORTHOLOGOUS_TO` edges bridging human genes to
model-organism orthologs and 9,495 `MODEL_OF` edges. Provenance, versions, and licenses for all
sources are in Supplementary Table S1. Data QC (no malformed identifiers, duplicates, or field
errors) is reported in Supplementary Section S2.

**Table 1. Graph composition.** *(fill from `graph_stats.json`)*

| Property | Value |
|---|---|
| Nodes | 448,931 |
| Directed edges | 5,860,539 |
| Unique undirected edges | 5,494,276 |
| Relation types | 28 |
| Target edges (human gene→disease) | 42,288 |
| ORTHOLOGOUS_TO edges | 510,703 |
| MODEL_OF edges | 9,495 |
| Monarch-derived fraction | 98.4% |

---

## 4. Methods

**Prediction task and splits.** We hold out 10% of target gene→disease edges as a test set and
10% as validation (seed 42); all remaining edges form the training graph. We evaluate four regimes:

- **R0 (standard):** full training graph; standard random-split test.
- **R1 (hub-controlled):** generic hubs (degree > 2,000) excluded as shared-neighbor evidence.
- **R2 (orthology-blocked):** `ORTHOLOGOUS_TO` and `MODEL_OF` edges removed from training (520,198
  edges; 8.9%), severing cross-species shortcuts while preserving graph scale.
- **R3 (de-leaked):** R1 + R2 + removal of eponymous/tautological and obsolete test targets.

**Methods compared.** Baselines: Random, Common-Neighbors, Adamic–Adar. KGE: TransE, DistMult,
ComplEx, RotatE (PyKEEN; bilinear models trained with 1-vs-all cross-entropy and Lp
regularization to avoid collapse). All methods are scored through one shared harness: for each
test edge, the true tail plus 50 type-matched disease negatives (excluding the gene's other known
diseases) are ranked; we report MRR and Hits@{1,3,10}, plus AUROC/AUPRC against random and
type-matched negatives. All results are mean ± 95% CI over seeds {42, 1, 7}.

**Degree decomposition.** To separate orthology leakage from degree bias, we compare the R0→R2
drop against a degree-preserving permutation null (following Gu et al., 2024) and report the
orthology effect residual to degree.

**Hybrid predictor.** A two-stage model: Adamic–Adar candidate generation followed by a
gradient-boosted re-ranker over [AA, common-neighbor count, KGE score, gene degree, disease
degree, rare-disease flag], trained on validation and evaluated on test.

**Validation.** Leave-one-source-out (retrain with one curated source removed; measure precision@k
on its sole-source edges) and newer-release recovery (predictions confirmed by a subsequent Monarch
release).

---

## 5. Results

### 5.1 The AUROC–ranking dissociation (R0)
Under the standard regime, KGE models achieve high classification scores (RotatE AUROC
[[0.99]]) yet rank the true disease far worse than topological baselines (Adamic–Adar MRR
[[0.69]], Hits@10 [[0.83]]; best KGE MRR [[0.09]], Hits@10 [[0.17]]). Reporting AUROC alone would
obscure this gap (Figure 1).

**Table 2. Benchmark under R0 (mean ± 95% CI, fake numbers).**

| Method | MRR | Hits@1 | Hits@10 | AUROC (hard) |
|---|---|---|---|---|
| Random | [[0.09±0.00]] | [[0.02]] | [[0.20]] | [[0.50]] |
| Common-Neighbors | [[0.65±0.01]] | [[0.52]] | [[0.83]] | [[0.87]] |
| Adamic–Adar | [[0.69±0.01]] | [[0.60]] | [[0.83]] | [[0.87]] |
| TransE | [[0.04]] | [[0.01]] | [[0.10]] | [[0.89]] |
| DistMult | [[0.11]] | [[0.05]] | [[0.19]] | [[0.86]] |
| ComplEx | [[0.10]] | [[0.04]] | [[0.18]] | [[0.88]] |
| RotatE | [[0.12]] | [[0.06]] | [[0.22]] | [[0.99]] |

### 5.2 Leakage collapse across regimes
Performance falls monotonically from R0 to R3 (Figure 2). Adamic–Adar drops from MRR [[0.69]]
(R0) to [[0.41]] (R2) to [[0.33]] (R3); the best KGE model from [[0.12]] to [[0.08]] to [[0.07]].
The orthology block (R2) accounts for the largest single drop ([[XX%]] relative), larger than the
hub block ([[YY%]]).

**Table 3. MRR by regime (fake numbers).**

| Method | R0 | R1 hub | R2 ortho | R3 de-leaked |
|---|---|---|---|---|
| Adamic–Adar | [[0.69]] | [[0.61]] | [[0.41]] | [[0.33]] |
| Best KGE (RotatE) | [[0.12]] | [[0.11]] | [[0.08]] | [[0.07]] |
| Hybrid (ours) | [[0.70]] | [[0.63]] | [[0.46]] | [[0.39]] |

### 5.3 Is orthology leakage separable from degree bias?
After the degree-preserving permutation control, [[XX%]] of the R0→R2 drop remains attributable to
orthology specifically (permutation p [[< 0.001 / = 0.12]]). *(This is the make-or-break result:
if the residual is large and significant, orthology is a distinct confound and the paper's core
claim holds; if it vanishes, we honestly report that orthology leakage is largely a manifestation
of the degree bias described by Gu et al. (2024).)*

### 5.4 A leakage-robust hybrid
Under R3, the hybrid attains MRR [[0.39±0.02]], exceeding the best baseline ([[0.33]]; paired
bootstrap p [[< 0.01]]), with the largest gains on low-degree rare-disease targets.

### 5.5 Biological validation
Leave-one-source-out: precision@50 = [[0.xx]] for held-out Gene2Phenotype edges. Newer-release:
[[N]] of the top [[100]] predictions were confirmed by a subsequent Monarch release. After
excluding eponymous and obsolete targets, [[M]] literature-supported novel candidates remain
(Supplementary Table S3), including [[example gene → disease]].

---

## 6. Discussion

Our results indicate that a substantial fraction of apparent gene–disease ranking performance on
multi-species KGs is [[attributable to / partly explained by]] cross-species orthology shortcuts.
[[If the residual-to-degree effect is significant:]] because this effect survives degree control,
it is a leakage source distinct from those in existing taxonomies (Ranga et al., 2025; Gu et al.,
2024) and specific to multi-species graphs. Practically, benchmarks on Monarch-scale graphs should
report an orthology-blocked regime, and novel-human-prediction claims should be evaluated with
cross-species evidence removed.

We emphasize that orthology is a legitimate biological signal, not an artifact; the issue is purely
one of *evaluation*. R0 (standard) is an optimistic upper bound and R3 a conservative lower bound
on honest performance; genuine prospective discovery lies between them.

## 7. Limitations
*(See `manuscript/limitations.md` — Monarch dominance and non-independence; four databases
evaluated but not integrated; source vintage incl. DrugCentral 2021; random rather than temporal
split; residual non-orthology leakage paths; predictions from the heuristic/hybrid rather than
KGE; sampled-negative ranking; single-run-adjacent KGE tuning.)*

## 8. Conclusion
Cross-species orthology is a quantifiable, biology-specific leakage source in gene–disease link
prediction on multi-species knowledge graphs. Controlling for it [[changes / does not change]]
the ranking of methods and reduces apparent performance by [[XX%]]. We release a de-leaked
protocol and harness to support honest evaluation.

## Data & code availability
Graph, splits (SHA-256 verified), and code: [[Zenodo DOI]] / [[GitHub URL]]. Fixed seed 42;
one-command reproduction in `run_all`.

## References *(real papers — verify and format to venue style)*
- Toutanova & Chen (2015), *Observed vs latent features for KG and text inference.*
- Dettmers et al. (2018), ConvE / WN18RR.
- Breit et al. (2020), *OpenBioLink*, Bioinformatics.
- Ranga et al. (2025), *Benchmarking the Impact of Data Leakage … Biomedical Link Prediction*, bioRxiv 2025.01.23.634511.
- Gu et al. (2024), *The effects of biological knowledge graph topology on embedding-based link prediction*, bioRxiv 2024.06.10.598277.
- Crichton et al. (2018), *Neural networks for link prediction in realistic biomedical graphs*, BMC Bioinformatics.
- Mungall et al. (2017), *The Monarch Initiative*, Nucleic Acids Research.
- McGary et al. (2010), *phenologs*, PNAS.
- [[verify all author names/years — placeholders above are approximate]]
