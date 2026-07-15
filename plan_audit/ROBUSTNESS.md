# Full-ranking robustness check (R0 vs R2)

*Generated 2026-07-15. Data: `data/processed/results/robustness/robustness_R{0,2}.json`.
Harness: `scripts/run_robustness.py`. Reconciles with `PAPER_SCOPE.md` and the memory notes
`benchmark-results-degree-driven` / `publication-strategy`.*

## Why this exists

Table 2 ranks each held-out true disease against only **50 type-matched sampled negatives**
(`lib_eval.rank_test_edges`, `n_neg=50`). That is standard and keeps KGE tractable, but the first
reviewer objection is: *"sampled-50neg MRR is optimistic; the real test is ranking against all
candidate entities."* This check answers that objection for the topological baselines — the methods
that carry the degree story — by re-ranking the true disease against the **entire filtered disease
pool** (D = 14,307 disease-category nodes; mean 14,299 negatives per edge after filtering the
gene's known tails). Same graph, same scorers (`run_baselines.make_scorers`, hub-cap 2000), same
tie-averaged rank rule as `lib_eval`; strictly harder ranking.

Correctness: the fast per-method ranking is **verified against a brute-force reference** (loop the
exact scorer over the full pool) on 150 test edges × 4 deterministic methods per regime — all ranks
agree exactly. Random is the analytic chance floor (true edge placed uniformly among D+1 slots).

**Scope:** topological baselines on R0 (standard) and R2 (degree-preserving null) — the pair the
robustness claim is about. KGE full-ranking needs the trained TransE/RotatE embeddings, which are
not on disk in this environment; it is a noted GPU-env follow-up (see *Limitations* below).

## Results

Full filtered ranking, 4,228 held-out human gene→disease edges, ~14,299-way per edge.

### MRR (headline)

| Method | R0 | R2 | Δ vs R0 | (sampled-50neg Δ, for reference) |
|---|---|---|---|---|
| Random (chance floor) | 0.0007 | 0.0007 | 0% | 0% |
| CommonNeighbors | 0.1177 | 0.0145 | **−88%** | −44% |
| AdamicAdar | 0.1326 | 0.0144 | **−89%** | −43% |
| Jaccard | 0.1061 | 0.0091 | **−91%** | −49% |
| **PreferentialAttachment** | 0.0251 | 0.0255 | **+1%** | +0.5% |

### Median rank (central tendency) and Hits@100

| Method | R0 MedRank | R2 MedRank | R0 Hits@100 | R2 Hits@100 |
|---|---|---|---|---|
| Random | 7152 | 7152 | 0.007 | 0.007 |
| CommonNeighbors | 1548 | **7154** | 0.361 | 0.150 |
| AdamicAdar | 1309 | **7154** | 0.378 | 0.160 |
| Jaccard | 4398 | **7152** | 0.286 | 0.086 |
| **PreferentialAttachment** | **156** | **147** | **0.441** | **0.440** |

(Pool midpoint ≈ 14,299 / 2 ≈ 7,150. A median rank there = chance.)

### Hits@1 / Hits@10 at R0

| Method | Hits@1 | Hits@10 | Mean rank |
|---|---|---|---|
| CommonNeighbors | 0.076 | 0.192 | 3605 |
| AdamicAdar | 0.091 | 0.214 | 3576 |
| Jaccard | 0.079 | 0.150 | 4086 |
| PreferentialAttachment | 0.000 | 0.081 | 1791 |

## What it means

1. **The degree pattern holds — and is stronger under full ranking.** Stripping degree-correlated
   structure (R2) drops the overlap heuristics **88–91%** in MRR (vs 43–49% at sampled-50neg) and
   pushes their **median rank to the pool midpoint (~7,153) = chance**. The residual "genuine
   structure" that looked substantial at sampled-50neg (overlap MRR ~0.25 at R2) is largely an
   artifact of the easy 50-negative setting: against 14,299 real candidates the overlap methods keep
   only a small top-hit residual (R2 MRR ≈ 0.014, ≈ 20× chance; R2 Hits@100 ≈ 0.15, ≈ 21× chance),
   while their bulk ranking is at chance. Their apparent skill is **overwhelmingly degree-driven**.

2. **PreferentialAttachment is the built-in positive control and is invariant.** Pure degree:
   MRR 0.025→0.026, median 156→147, Hits@100 0.441→0.440 — unchanged across the degree-preserving
   null, confirming R2 preserves the degree sequence exactly. It carries the true disease into a good
   neighborhood (median rank ~150 / 14,299; 44% within top-100) but essentially never ranks it #1
   (Hits@1 = 0.000), because ordering diseases by raw degree rarely puts the *specific* partner
   first.

3. **A metric dissociation worth reporting (sub-result).** Under full ranking at R0, MRR and
   central-tendency metrics rank the baselines **oppositely**: AdamicAdar wins by MRR (0.133) but is
   near-worst by median rank (1,309), while PreferentialAttachment is worst by MRR (0.025) but far
   best by median rank (156) and Hits@100 (0.44). The overlap heuristics occasionally nail the exact
   partner (higher Hits@1/MRR) but bury most; degree ranks everything moderately well but nothing
   first. This reinforces the paper's broader point that a single top-heavy headline metric — and
   especially few-negative sampling — misleads. **Caveat for the write-up:** the "PreferentialAttachment
   is the strongest single baseline" statement is *protocol-dependent* (true by sampled-50neg MRR and
   by full-ranking median/Hits@100; false by full-ranking MRR) and must be stated with the metric named.

**Bottom line:** the centerpiece — overlap-heuristic performance is largely degree, and removing
degree collapses it while pure degree is untouched — is **confirmed and reinforced** by full ranking.
The sampled-50neg protocol, if anything, *understated* the degree dominance.

## Limitations / follow-ups

- **KGE not yet full-ranked.** TransE/RotatE full-ranking requires the trained d64 embeddings (not on
  disk here). Sampled-50neg already shows the same-direction degree drop (TransE −24%, RotatE −21%
  R0→R2); a full-ranking KGE pass in the GPU environment would close the loop and is the one
  remaining piece. `run_robustness.py` scores from a `score_fn(gene, disease)`, so plugging in a KGE
  scorer is a small addition.
- **R1/R3 not full-ranked** (out of scope for this check — they are null at sampled-50neg, so a
  robustness pass is low value; `--regimes R0,R1,R2,R3` runs them if wanted).

## Reproduce

```
python scripts/run_robustness.py --regimes R0,R2   # ~3 min, writes robustness_R{0,2}.json
```
