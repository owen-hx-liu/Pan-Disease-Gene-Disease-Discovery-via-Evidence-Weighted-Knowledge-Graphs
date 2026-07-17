# Quantifying evaluation leakage in multi-species gene–disease link prediction: a graded audit and a de-leaked benchmark

> **DRAFT — for format/feel only.** Numbers are real (from the frozen result files); prose is a complete first draft; **citations are unverified placeholders** (author-year keys only — confirm every one, and resolve the Bradshaw&Layer/Gu conflict, before any real use). Not for submission.

**Owen Liu**¹ (ORCID 0009-0003-0476-5270), **Ryan Tang**¹, **Daniel Augustine**¹, **[Mentor co-author — TBD]**²

¹ Independent Researcher. ² [Affiliation].
Correspondence: [email].

---

## Abstract

Knowledge-graph embedding (KGE) models predict gene–disease associations with reported mean reciprocal rank (MRR) above 0.8 and area under the ROC curve (AUROC) above 0.95, and such scores increasingly guide disease-gene prioritization. These numbers can be inflated by several distinct forms of evaluation leakage, which are usually reported as a single undifferentiated "optimism." We present a graded audit that separates and quantifies these mechanisms on one graph under one ranking harness. Five topological baselines (Random, Common-Neighbors, Adamic-Adar, Jaccard, Preferential-Attachment) and two KGE models (TransE, RotatE) are evaluated on 4,228 held-out human gene–disease edges across four regimes — R0 (standard random split), R1 (redundancy; direct held-pair edges removed), R2 (a degree-preserving null that destroys genuine structure while holding the degree sequence exactly), and R3 (orthology-blocked; cross-species bridges removed) — with three seeds, on a 452,012-node Monarch-derived graph, and the whole audit is replicated on Hetionet. We find that evaluation leakage is *method-class-dependent*. Node-degree structure (R2) is the one universal leak: overlap heuristics lose 43–49% of their MRR and KGE models 21–24%, while Preferential-Attachment — which uses degree alone — is unchanged, isolating degree as the shared source. Stratifying the held-out edges by degree exposes a double dissociation — overlap heuristics ride the query gene's degree while KGE models and a pure-degree predictor ride the ranked disease's — but under the degree-null every method loses most of its score on the rare, low-degree end (KGE −62%, overlap heuristics −84 to −86% on the rarest-disease quartile) and little on the popular end, so the headline number is carried by popular, well-connected entities rather than by the rare diseases that prioritization aims to find. Train/test redundancy (R1) leaks only into multi-hop KGE (6–8.5%) and is invisible to overlap heuristics. Cross-species orthology (R3), the intuitive culprit, is a null (within 0.5% for every method). Under a stricter filtered full-ranking protocol the degree collapse strengthens (TransE −71%). We release the de-leaked splits, all per-regime result files, and a one-command reproduction (DOI 10.5281/zenodo.21398832). Because comparable studies rarely report confidence intervals, multiple seeds, significance tests, or a formal degree null, we intend this as a reusable, rigorously controlled reference rather than a novel positive finding.

**Keywords:** knowledge graph; knowledge graph embedding; link prediction; gene–disease association; data leakage; benchmark evaluation; node degree bias; reproducibility; Monarch Initiative; TransE; RotatE.

---

## 1. Introduction

Predicting which genes are associated with which diseases is a central task in computational biology, and it is increasingly framed as link prediction on a biomedical knowledge graph. Knowledge-graph embedding (KGE) models routinely report mean reciprocal rank (MRR) above 0.8 and area under the ROC curve (AUROC) above 0.95 for this task, and these scores are used to prioritize candidate disease genes for experimental follow-up. When a reported number determines which genes a laboratory investigates, the validity of that number is not a methodological nicety but a practical concern.

A growing literature warns that such scores are inflated by *evaluation leakage* — information that lets a model appear to generalize while it is in fact exploiting an artifact of how the data were split or sampled (Kapoor & Narayanan, 2023; Bernett et al., 2024a). In biomedical link prediction, at least three distinct mechanisms have been described. First, *train/test redundancy*: near-duplicate, reversible, or Cartesian-product relations leave copies of test facts in the training graph, and removing them sharply reduces measured performance (Akrami et al., 2020; Toutanova & Chen, 2015). Second, *node-degree (popularity) bias*: standard random-edge sampling favors high-degree hub nodes, so that a predictor using degree alone can rival sophisticated models (Yılmaz et al., 2025; Aiyappa et al., 2024; Shomer et al., 2023). Third, *non-representative test sets*: performance on a random test split can be far higher than on the rare or understudied entities that motivate the method in the first place (Brière & Baudot, 2025).

These mechanisms are usually studied one at a time, on different graphs, with different code, and are frequently collapsed into a single caveat about "optimistic evaluation." What is missing is a graded decomposition that (i) separates the mechanisms on one graph under one harness, (ii) measures the *magnitude* each contributes — not merely whether it exists — via a formal degree-preserving null, and (iii) releases the resulting de-leaked splits so that others can adopt them.

We provide such an audit for gene–disease link prediction. Our contributions are:

1. **A graded regime framework (R0–R3)** run through a single shared ranking harness for both topological and KGE methods, so that no method class is scored by different code.
2. **A degree-preserving null (R2)** that isolates and quantifies the share of performance attributable to node degree versus genuine structure, validated by a pure-degree predictor that the null leaves unchanged.
3. **The finding that leakage is method-class-dependent**, with orthology and redundancy serving as designed negative controls — one of which (orthology) returns a clean null against the domain-specific hypothesis.
4. **A reproducible, DOI-archived de-leaked benchmark**, replicated on a second, independently constructed graph.

We are explicit about novelty. The core phenomenon — that random-split evaluation of biomedical link prediction is inflated by degree and leakage — is already established, including in a top venue and on the very graph we use (Yılmaz et al., 2025; Bradshaw & Layer, 2024; Brière & Baudot, 2025). Our contribution is a controlled, quantitative consolidation across mechanisms, method classes, and graphs, released as a resource, together with a level of statistical rigor that comparable studies seldom report.

## 2. Related work

**Degree and leakage in biomedical link prediction.** Three studies are our closest neighbors. Yılmaz et al. (2025) show that random-edge-sampling evaluation is biased toward high-degree nodes on protein–protein interaction networks, that this bias persists across temporal snapshots, and that "a simple predictive model that uses bias as the only feature outperforms sophisticated models"; they use Preferential Attachment as a pure-bias reference and propose degree-aware evaluation. Bradshaw & Layer (2024) study KGE link prediction for rare-disease gene prediction on the Monarch graph — the same graph and task we use — and report a "bias cycle" in which degree bias compounds; notably, filtering Monarch to roughly a tenth of its size *improves* performance. Brière & Baudot (2025) build a systematic framework for three leakage types in KGE-based drug repurposing on a PrimeKG-derived graph, running seven KGE models. We differ from all three by combining topological and KGE methods under one harness, by decomposing the *magnitude* of degree leakage with a permutation null rather than a binary retrain, by reporting confidence intervals, seeds, and significance tests, by replicating across two graphs, and by adding an orthology regime specific to multi-species biology.

**Foundations of the regimes.** Regime R1 follows the train/test-separation line of Akrami et al. (2020) and the reversible-relation removal that produced FB15k-237 (Toutanova & Chen, 2015). Regime R2 rests on degree-preserving edge randomization (Maslov & Sneppen, 2002), the same primitive that exposed degree exploitation in sequence-based PPI prediction (Bernett et al., 2024b). The general framing of leakage as a reproducibility threat follows Kapoor & Narayanan (2023) and the biomedical leakage taxonomy of Bernett et al. (2024a). Degree bias in embedding-based completion specifically is documented by Shomer et al. (2023) and Aiyappa et al. (2024).

**Reconciliation flagged for the Discussion.** Brière & Baudot (2025) report *no* evidence that their models rely on degree (their regime DL2), whereas we find degree to be the single universal leak (our R2). We return to this apparent contradiction in Section 6; briefly, a binary "remove-degree-and-retrain" test and a magnitude-decomposing null answer different questions.

## 3. Materials: the knowledge graph

Our primary graph is derived from the Monarch Initiative (Putman et al., 2024), augmented with seven curated non-Monarch sources (DGIdb, GWAS Catalog, DrugCentral, Orphadata, Gene2Phenotype, CIViC, and DIDA). We describe it honestly as a *Monarch-derived* graph: approximately 98.5% of its edges are Monarch content. Because Monarch is itself an aggregator of primary sources, "supported by N databases" should not be read as N independent lines of evidence; per-source contributions are provided with the release.

The graph contains 452,012 nodes and 5,860,539 directed edges across 28 relation types and 8 node categories (Table 1). Its degree distribution is heavy-tailed — mean 25.9, median 2, maximum 20,662 — so that most nodes are sparsely connected while a few hubs dominate. It includes 336,460 gene nodes, 14,339 disease nodes, 42,288 gene–disease target edges, and 510,703 cross-species orthology edges. The prediction target is the relation GENE_ASSOCIATED_WITH_CONDITION, restricted to the 4,228 held-out human gene–disease test edges.

For replication we use Hetionet v1.0 (Himmelstein et al., 2017), an independently constructed biomedical graph with 47,031 nodes and 2,250,197 edges (Table 1). Hetionet is single-species and carries no orthology, so regime R3 is Monarch-only; it is also markedly denser (median degree 32 versus 2), which, as we show, makes its degree collapse smaller than Monarch's.

Code is released under the MIT License and the Monarch-derived data under CC BY 4.0, preserving attribution to Monarch and the underlying primary sources.

*[Table 1 here — graph statistics; `tables/table1_graph_stats.tex`, `\label{tab:graph_stats}`.]*
*[Figure S1 — schema/metagraph; Figure S2 — degree distribution. To be generated.]*

## 4. Methods

### 4.1 Graded leakage regimes

We evaluate the same test edges under four training regimes that differ only in what leakage-bearing information the training graph is allowed to contain.

- **R0 (standard).** An ordinary random train/validation/test split — the condition most papers report.
- **R1 (redundancy).** We remove the 1,661 direct held-pair edges that place a training copy of a test fact in the graph (0.03% of all edges). This is a targeted removal, not a reduction in training volume.
- **R2 (degree-preserving null).** We rewire the graph to destroy genuine structure while holding each node's degree exactly fixed (Maslov & Sneppen, 2002). Any performance that survives this null is attributable to degree alone.
- **R3 (orthology-blocked).** We remove the cross-species ORTHOLOGOUS_TO bridges (520,198 edges; training edges fall from 5,852,083 to 5,331,885). This tests the domain-specific hypothesis that cross-species transfer is a source of inflation.

### 4.2 The degree-preserving null

Regime R2 is the methodological centerpiece. Using degree-preserving edge randomization, we generate graphs that match the real graph's degree sequence node-for-node but scramble which specific nodes are connected. Scoring a method on the real R0 graph and on the null yields a decomposition: the *structure* component is the difference (real − null), and its statistical significance is assessed by permutation over ten replicates. A method whose score is unchanged by the null owes its performance entirely to degree; Preferential-Attachment, which ranks candidates by degree product, is exactly such a method and serves as a built-in positive control.

### 4.3 Models and baselines

We evaluate two KGE models, TransE (Bordes et al., 2013) and RotatE (Sun et al., 2019), each at embedding dimension 64 for 300 epochs with a self-adversarial negative-sampling loss, the stochastic local closed-world assumption training loop, batch size 16,384, and 16 negatives per positive, using PyKEEN (Ali et al., 2021). We attempted two further families; both proved impractical at this scale and are reported only as such. ComplEx (Trouillon et al., 2016) failed to train stably at 452,012 entities across repeated attempts, and DistMult (Yang et al., 2015) was evaluated only in short smoke runs. We therefore do not report full ComplEx or DistMult results, and note DistMult as the feasible future addition.

The five topological baselines — Random, Common-Neighbors, Adamic-Adar, Jaccard, and Preferential-Attachment (Liben-Nowell & Kleinberg, 2007; Barabási & Albert, 1999) — provide method-class contrast and, in the case of Preferential-Attachment, the pure-degree anchor.

### 4.4 Shared evaluation harness

Every method — topological and KGE alike — is scored by a single ranking harness, so that differences between methods reflect the methods and not their evaluation code. All regimes and both graphs use the same harness.

### 4.5 Ranking protocols

We report two protocols and name which every number comes from. The main protocol ranks each true disease against 50 type-matched sampled negatives. The robustness protocol is filtered full ranking, in which each true disease is ranked against the entire pool of roughly 14,300 diseases; 30 test genes absent from the training vocabulary receive the worst possible rank. The two protocols can reorder methods — Preferential-Attachment is the strongest single baseline by sampled-negative MRR but the weakest by full-ranking MRR — so we always state the protocol alongside the metric.

### 4.6 Statistics

We report bootstrap confidence intervals on all metrics, average over three seeds (42, 1, 7), and assess method differences by a paired bootstrap (1,000 resamples, n = 4,228). The degree decomposition is tested by permutation over ten null replicates. This combination — confidence intervals, multiple seeds, a significance test, and a formal null — is, as Section 6 notes, reported by few comparable studies.

## 5. Results

### 5.1 The graded audit

Table 2 gives MRR for every method in every regime. The pattern is immediate and method-class-dependent: only the degree-preserving null (R2) materially reduces performance, and it reduces it for every method that uses structure. The overlap heuristics fall by 43–49% (Common-Neighbors −44%, Adamic-Adar −43%, Jaccard −49%), and the KGE models by 21–24% (TransE −24%, RotatE −21%). Preferential-Attachment, which uses degree alone, is unchanged (0.587 to 0.589), and Random is unchanged at 0.091. Redundancy (R1) and orthology (R3) leave all methods within noise at this protocol. The unchanged pure-degree predictor confirms that R2 removes structure while preserving degree, and thereby isolates node degree as the single leakage source common to all method classes.

*[Table 2 here — leakage audit; `tables/table2_audit.tex`, `\label{tab:audit}`.]*
*[Figure 2 here — leakage audit bar chart; `figures/fig2_leakage_audit`.]*

### 5.2 Degree is the universal leak

Table 3 decomposes R0 performance into degree and structure. For the overlap heuristics roughly 55% of performance is attributable to degree and 45% to genuine structure: Common-Neighbors retains a structure residual of +0.192 over its degree-null score, Adamic-Adar +0.194, and Jaccard +0.196, each significant at permutation p < 0.001. For Preferential-Attachment the residual is −0.002 (p = 1.000) and for Random −0.000 (p = 1.000) — statistically zero, as expected for degree-only and structure-free predictors. Preferential-Attachment's real and null scores (0.589 versus 0.591) coincide, re-confirming that the null preserves the degree sequence.

*[Table 3 here — degree-null decomposition; `tables/table3_degree_null.tex`.]*
*[Figure 3 here — degree vs structure; `figures/fig3_degree_vs_structure`.]*

### 5.3 Redundancy is a multi-hop-only leak

Removing the 1,661 direct held-pair edges (R1) reduces KGE performance — TransE by 6.2% and RotatE by 8.5% — but leaves every two-hop overlap heuristic flat (±0.3%). The mechanism is structural: a direct gene–disease edge is embedded by a KGE model as proximity in the vector space, so removing it lengthens the predicted distance, but it is not a shared neighbor and so never enters a Common-Neighbors, Adamic-Adar, or Jaccard score. Because the removed edges are only 0.03% of the graph, the KGE drop is not a training-volume artifact. We note for transparency that this reverses our pre-registered expectation, which had anticipated orthology, not redundancy, to be the multi-hop leak.

### 5.4 Orthology is a null

Removing the 520,198 cross-species orthology edges (R3) changes no method's MRR by more than 0.5%, and the KGE models actually improve marginally (+0.2 to +0.4%). Cross-species orthology, the mechanism a biomedical reviewer would most naturally suspect of inflating a multi-species benchmark, does not inflate performance here. A clean null against a plausible domain hypothesis is itself a result, and it is one that a designed negative control is meant to deliver.

### 5.5 Full-ranking robustness

To rule out that the degree effect is an artifact of ranking against only 50 sampled negatives, we re-rank TransE against the entire disease pool (Table S1). The degree collapse does not merely persist; it strengthens, from −24% at the sampled protocol to −71% under full ranking (MRR 0.079 at R0 to 0.023 at R2), with redundancy at −10%. The sampled-negative protocol therefore *understates* degree dominance. We ran this stricter protocol for TransE as representative of the KGE class, which tracks RotatE within a few points throughout Table 2; full-ranking RotatE and R3 are straightforward extensions we did not judge necessary, since R3 is already a null.

*[Table S1 here — full-ranking robustness; `tables/tableS1_fullrank_robustness.tex`.]*

### 5.6 Independent-graph replication

The audit reproduces on Hetionet (Table 4), a graph built by a different group with no shared pipeline. The degree-preserving null again collapses the overlap heuristics — Common-Neighbors −35%, Adamic-Adar −37%, Jaccard −36% — while Preferential-Attachment is unchanged and Random is unchanged. The effect is therefore a property of biomedical knowledge-graph topology rather than of Monarch or of our particular construction. The Monarch collapse (43–49%) exceeds the Hetionet collapse (35–37%), consistent with Monarch's much heavier degree tail (median degree 2 versus 32).

*[Table 4 here — Hetionet replication; `tables/table4_hetionet.tex`.]*

### 5.7 AUROC hides what ranking metrics expose

AUROC is nearly blind to the leakage that ranking metrics reveal (Figure 1). Against random negatives it stays near 0.99 in every regime; even against type-matched negatives it moves only from about 0.96 at R0 to about 0.88 at R2, while MRR over the same span falls by tens of percent. A benchmark that reports AUROC alone will not detect the degree inflation that a ranking metric surfaces immediately (Saito & Rehmsmeier, 2015). This dissociation is a practical warning about metric choice, not merely a curiosity.

*[Figure 1 here — AUROC–MRR dissociation; `figures/fig1_auroc_mrr_dissociation`.]*

### 5.8 A degree double dissociation, and where the null bites

Section 5.2 showed at the aggregate level that degree drives performance; a per-edge stratification localizes *which* degree, for *which* method, and makes the practical cost concrete. Because every run stores a reciprocal rank for each of the 4,228 test edges — the KGE runs natively, and the baselines once we add the same per-edge dump — we re-bin those edges by node degree in the training graph and recompute MRR within each bin without retraining. We stratify along two axes, the query gene (source) and the ranked disease (target), into equal-count quartiles (Table 5, Figure 4).

The two axes reveal a double dissociation across method classes. The **overlap heuristics ride the query gene's degree**: Common-Neighbors, Adamic-Adar, and Jaccard climb from an MRR of roughly 0.14 on the lowest gene-degree quartile to 0.53–0.62 on the highest (a lift of about +0.40 to +0.48), because a low-degree gene has too few neighbors to form the shared-neighbor overlaps these methods depend on; along the ranked-disease axis they are nearly flat. The **KGE models and Preferential-Attachment ride the ranked disease's degree**: Preferential-Attachment climbs from 0.06 on the rarest-disease quartile to 0.96 on the most popular, and TransE and RotatE from about 0.68 to 0.86 and 0.92, while all three are flat or slightly inverted in gene degree. Random is flat on both axes (≈0.09). Which entity's popularity inflates a method is thus method-class-dependent — a per-edge sharpening of the method-class dependence of Section 5.1 — yet every non-trivial method is steeply degree-graded on one axis or the other.

Crossing the ranked-disease stratification with the degree-preserving null (R2) confirms that these gradients *are* the degree leak rather than a proxy for it, and does so identically across classes. Within the rarest-disease quartile the null strips most of the score — Common-Neighbors −86%, Adamic-Adar −85%, Jaccard −84%, TransE and RotatE −62% — because a rare target's R0 score is almost entirely genuine structure, which the null destroys. Within the most-popular quartile it strips little — between −3% and −17% — because a popular target's score is largely degree, which the null preserves. Preferential-Attachment, being pure degree, is unchanged throughout (+1% and 0%), reproducing the built-in control of Section 4.2 at per-edge resolution. The practical reading is uncomfortable: the rare diseases and low-connectivity genes that most motivate computational prioritization are exactly the cases every method predicts worst *and* the cases whose apparent skill the degree-null most erodes; the headline number is propped up by popular, well-connected entities that a pure-degree predictor could already rank.

*[Table 5 here — degree-stratified performance across method classes; `tables/table5_degree_stratified.tex`, `\label{tab:degree_stratified}`.]*
*[Figure 4 here — the degree double dissociation (gene vs disease axis); `figures/fig4_degree_stratified`.]*

### 5.9 A worked example: what the leak costs prioritization

The stratified numbers translate into concrete predictions a biologist would act on. The benchmark asks, for each held-out gene, whether its true disease is ranked above 50 type-matched candidate diseases; Table 6 shows what that ranking rewards and what it misses.

Consider *insomnia* (MONDO:0013600), a high-degree hub (training degree 801) that is the held-out target of 86 genes. Pure-degree Preferential-Attachment ranks insomnia first for 97% of those genes, RotatE for 85%, and shared-neighbour Adamic-Adar for only 43%: a predictor using nothing but node degree matches or beats the trained KGE, while the method that uses actual shared-neighbour structure trails. RotatE's mean MRR on these edges barely moves under the degree-null (0.932 to 0.882). The extreme case is the antisense RNA *LIN28B-AS1*, whose entire presence in the graph is a single edge (degree 1): RotatE nonetheless ranks insomnia as its number-one associated condition, the rank survives the degree-null, and Preferential-Attachment reproduces it, while Adamic-Adar — finding no shared-neighbour path — places insomnia 26th. The confident "prediction" encodes that insomnia is popular, not anything about the gene.

The mirror image is where prioritization actually operates. Rare diseases with a training degree of 1–5 carry a RotatE R0 MRR of 0.637 that *collapses to 0.062* under the degree-null, with Preferential-Attachment at 0.024 — on rare targets the apparent skill is genuine structure, and degree offers essentially nothing; on hub diseases (degree ≥ 301) the ranking is instead almost pure degree, where Preferential-Attachment (0.957) actually exceeds RotatE (0.913). The cost is visible in canonical, disease-defining associations of heavily annotated genes: *RS1*, the gene mutated in X-linked retinoschisis, is ranked 51st of 51 — dead last — for retinoschisis; *NF1* is ranked 43rd for plexiform neurofibroma, a hallmark tumour of neurofibromatosis type 1; and *KIT*, a recurrent driver in mucosal melanoma, is ranked 49th. In each the gene is well connected (degree 350–1,387) but the specific disease node is rare (degree 1), and the model ranks the true, textbook association below dozens of unrelated diseases. A headline MRR aggregated over insomnia-like hubs would advertise a model that, on exactly the specific associations a prioritization study exists to surface, performs at or below chance.

*[Table 6 here — worked examples of degree-driven misranking; `tables/table6_case_study.tex`, `\label{tab:case_study}`.]*

## 6. Discussion

**What the audit means for practice.** The single most useful thing a gene–disease benchmark can add is a degree-preserving null: if a model barely outperforms its own degree-null score, that gap — not the headline number — is the real signal. On our graph RotatE retains an MRR of 0.655 at R2 against Preferential-Attachment's 0.589, so a substantial fraction of its apparent skill is degree. Reporting a pure-degree anchor, ranking metrics with hard or type-matched negatives rather than AUROC alone, and the ranking protocol used, would let readers see this directly.

**Reconciling our R2 with Brière & Baudot's DL2.** Brière & Baudot (2025) remove degree information and retrain, and find performance largely preserved, concluding their models do not rely on degree; we find degree to be the universal leak. The results are compatible once the questions are separated. Their test is binary — does removing degree destroy performance? — whereas ours decomposes how much of the score is degree versus structure. A model can retain much of its rank ordering while still owing roughly half of its score to degree. The graphs differ as well: gene–disease prediction on Monarch, with a median degree of 2, is far more degree-dominated than drug repurposing on their PrimeKG-derived graph. Notably, their own non-representative-test result — random-split MRR near 0.2 against real-inference MRR near 0.01, with well-studied entities predicted far better — is itself a degree/study-bias signature consistent with our finding.

**Where we sit relative to the field.** Across comparable studies, confidence intervals, multiple seeds, significance tests, and null models are rare; two-graph replication is rarer still. That rigor is our principal contribution alongside the release. Our degree-stratified analysis (Section 5.8) and the worked example (Section 5.9) also connect directly to the held-out-disease evaluations of Brière & Baudot (2025) and Bradshaw & Layer (2024): the rare-target quartile is where the degree-null bites hardest, which is precisely the regime in which those studies report their sharpest optimism.

**Recommendations.** We suggest that gene–disease link-prediction benchmarks (i) publish a degree-preserving null alongside headline numbers, (ii) report Preferential-Attachment and Random anchors, (iii) prefer ranking metrics with hard negatives over AUROC alone, (iv) state the ranking protocol explicitly, and (v) release de-leaked splits with seeds and hashes.

**Limitations.** The graph carries no genuine edge-discovery dates, so our regimes control optimism on a random split rather than a temporal one. Coverage is Monarch-dominated even after integrating seven further sources. The central phenomenon is established in prior work, so our novelty is consolidation, decomposition, and release rather than discovery. Full-ranking robustness was run for TransE only, ComplEx is untrainable at this scale, and DistMult was only smoke-tested, so the KGE arm spans two full models. The worked example (Section 5.9) is illustrative rather than a systematic clinical evaluation, and the curated edges, though representative of the degree buckets they are drawn from, are selected for biological legibility.

## 7. Conclusion

Evaluation leakage in gene–disease link prediction is not one thing but several, and they do not affect all methods equally. Node-degree structure is a universal leak that inflates every structural method and strengthens under stricter ranking; train/test redundancy leaks only into multi-hop embedding models; and cross-species orthology, the intuitive suspect, does not leak at all. We release the de-leaked splits, results, and reproduction pipeline (DOI 10.5281/zenodo.21398832) so that these controls can be adopted directly rather than re-derived.

## Data and code availability

All code (MIT) and de-leaked data (CC BY 4.0) are archived at Zenodo (DOI 10.5281/zenodo.21398832) and on GitHub (github.com/owen-hx-liu/gene-disease-leakage-benchmark). The release includes the de-leaked splits with a manifest (SHA-256 bc468864…), per-regime baseline, KGE, robustness, and Hetionet result files, table and figure generators that read only the frozen result files, and a one-command reproduction script. All runs use fixed seeds (42, 1, 7).

## Author contributions

O.L. conceived the audit, built the pipeline, ran the experiments, and drafted the manuscript. R.T. and D.A. contributed to [data integration / analysis / figures]. [Mentor] supervised the study and revised the manuscript. All authors approved the final version.

## Competing interests

The authors declare no competing interests.

---

## References

*(Unverified placeholders — confirm every entry, update preprints to published versions, and resolve the Bradshaw & Layer / Gu 2024 author conflict before use.)*

- Aiyappa, R. et al. (2024). Implicit degree bias in the link prediction task. *arXiv:2405.14985*.
- Akrami, F. et al. (2020). Realistic re-evaluation of knowledge graph completion methods. *SIGMOD*.
- Ali, M. et al. (2021). PyKEEN 1.0: a Python library for training and evaluating knowledge graph embeddings. *JMLR*.
- Barabási, A.-L. & Albert, R. (1999). Emergence of scaling in random networks. *Science* 286.
- Bernett, J. et al. (2024a). Guiding questions to avoid data leakage in biological machine learning. *Nature Methods* 21:1444.
- Bernett, J. et al. (2024b). Cracking the black box of deep sequence-based protein–protein interaction prediction. *Briefings in Bioinformatics* 25(2):bbae076.
- Bordes, A. et al. (2013). Translating embeddings for modeling multi-relational data. *NeurIPS*.
- Bradshaw, [initials] & Layer, R. (2024). The effects of biological knowledge graph topology on embedding-based link prediction. *bioRxiv* 2024.06.10.598277. **[VERIFY authors]**
- Brière, G. & Baudot, A. (2025). Benchmarking data leakage on link prediction in biomedical knowledge graph embeddings. *bioRxiv* 2025.01.23.634511.
- Chandak, P. et al. (2023). Building a knowledge graph to enable precision medicine (PrimeKG). *Scientific Data* 10:67.
- Gualdi, F. et al. (2024). Predicting gene–disease associations for diseases with curtailed information. *bioRxiv* 2024.01.11.575314.
- Himmelstein, D. et al. (2017). Systematic integration of biomedical knowledge prioritizes drugs for repurposing (Hetionet). *eLife* 6:e26726.
- Kapoor, S. & Narayanan, A. (2023). Leakage and the reproducibility crisis in machine-learning-based science. *Patterns* 4(9):100804.
- Liben-Nowell, D. & Kleinberg, J. (2007). The link-prediction problem for social networks. *JASIST* 58(7).
- Maslov, S. & Sneppen, K. (2002). Specificity and stability in topology of protein networks. *Science* 296:910.
- Putman, T. et al. (2024). The Monarch Initiative in 2024. *Nucleic Acids Research* 52(D1):D938.
- Saito, T. & Rehmsmeier, M. (2015). The precision-recall plot is more informative than the ROC plot on imbalanced datasets. *PLoS ONE* 10(3):e0118432.
- Shomer, H. et al. (2023). Toward degree bias in embedding-based knowledge graph completion. *WWW*.
- Sun, Z. et al. (2019). RotatE: knowledge graph embedding by relational rotation in complex space. *ICLR*.
- Toutanova, K. & Chen, D. (2015). Observed versus latent features for knowledge base and text inference. *CVSC Workshop*.
- Trouillon, T. et al. (2016). Complex embeddings for simple link prediction. *ICML*.
- Yang, B. et al. (2015). Embedding entities and relations for learning and inference in knowledge bases (DistMult). *ICLR*.
- Yılmaz, S., Yorgancioglu, A. & Koyutürk, M. (2025). Bias-aware training and evaluation of link prediction in network biology. *PNAS* 122(24):e2416646122.
