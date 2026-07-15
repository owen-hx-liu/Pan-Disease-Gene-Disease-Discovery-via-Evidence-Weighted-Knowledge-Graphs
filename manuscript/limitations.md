# Limitations (manuscript draft — Discussion section)

*Drop-in prose for the paper's Limitations subsection. Written proactively (disclosing a
weakness reads as rigor; the same fact found by a reviewer reads as concealment). Replace
`<<...>>` with final numbers once the runs complete. Keep every claim scoped to what we
actually did — never overstate what a source or method can do.*

---

## Limitations

**Knowledge-graph provenance.** Although our graph incorporates eight biomedical databases,
it is not an equal-weight integration: approximately 98.4% of edges derive from the Monarch
Initiative, with seven curated sources (DGIdb, GWAS, DrugCentral, Orphadata, Gene2Phenotype,
CiVIC, DIDA) contributing the remaining ~1.6%. We therefore describe the resource as a
Monarch-derived multi-species knowledge graph augmented with curated sources, rather than a
novel multi-database integration, and we make no novelty claim about the graph itself; our
contribution is the leakage-aware evaluation and the associated method. Because Monarch is
itself an aggregator of primary resources (HPO, MGI, ZFIN, ClinGen, OMIM, Orphanet, and
others), edge support across sources is not independent, and "supported by N databases" does
not represent N independent observations. Our consensus-confidence model accounts for this by
collapsing aggregator votes (Methods), but the practical consequence is that multi-source
corroboration is available for only a small subset of edges.

**Sources considered but not integrated.** Four additional databases were evaluated but not
incorporated in this release, owing to limitations of the specific dumps we obtained rather
than of the databases themselves: the DisGeNET dump was a variant-level annotation table
without a gene–disease-association column; the ClinGen/PharmGKB dump contained node and
cross-reference tables but no relationship table; LNCipedia provided sequence (FASTA) data
only; and LncRNADisease encoded both partners as free text with no offline ontology index for
unambiguous mapping. Their exclusion does not affect the human gene–disease content that our
analyses target. Provenance, versions, download dates, and licenses for all sources — included
and excluded — are reported in the supplementary provenance table.

**Source vintage.** The constituent sources span 2019–2026; notably, the DrugCentral export is
a 2021 release. Because the graph functions here as a fixed benchmark testbed rather than a
currency-focused resource, and because the leakage phenomena we study are structural properties
of the graph's cross-species topology (independent of any single source's release date), source
vintage does not affect our central findings. Drug–target edges in particular are peripheral to
the gene→disease prediction task. We record exact versions and dates so that the analysis is
fully reproducible against the same snapshot.

**Evaluation split.** The dataset contains no genuine edge-discovery dates (an earlier temporal
analysis simulated years and is not used here), so a time-based train/test split — the ideal
guard against retrospective leakage — is not possible. We therefore use a fixed random held-out
split (seed 42) and mitigate leakage through the R0–R3 regimes and the leave-one-source-out and
newer-release validations instead. A genuine prospective evaluation must await datasets with
reliable association dates.

**Residual leakage.** Our de-leaked regime (R3) blocks the dominant cross-species pathway by
removing orthology and model-organism bridges and controlling for generic hubs, tautological
targets, and source non-independence. It cannot guarantee the removal of every indirect
information path (e.g., leakage mediated by shared pathway or protein-complex membership), so
the R3 numbers should be read as a conservative-but-not-exhaustive lower bound on honest
performance. We report the R0→R3 drop of <<X>> as evidence of the magnitude, not as a claim of
zero residual leakage.

**Prediction method and validation.** The strongest predictor under our protocol is a
topological heuristic (Adamic–Adar), not a knowledge-graph embedding model; our reported novel
hypotheses therefore derive from the heuristic (and the hybrid), and the KGE models are
benchmarked rather than used for discovery. The KGE models were tuned over a limited grid
(dimensions <<{128,256}>>, epochs <<{100,300}>>) with <<three>> seeds; a more exhaustive sweep
could narrow — though, given the order-of-magnitude ranking gap, is unlikely to close — the
difference. <<DistMult required [regularization/loss change] to train stably; we report the
corrected configuration.>> Biological face-validity is retrospective: we exclude eponymous or
tautological gene–disease pairs (where the disease is defined by the gene) and obsolete disease
terms from the headline precision, and report them separately.

**Computational approximations.** For tractability at 5.86M edges, ranking is evaluated against
sampled type-matched negatives rather than the full entity set for every method, and selected
network metrics (e.g., betweenness) are estimated from samples; the reported effect sizes are
large relative to the resulting sampling error, but exact values would require a compiled graph
library.
