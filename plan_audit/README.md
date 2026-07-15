# plan_audit/ — the pivoted plan (Leakage Audit + De-Leaked Benchmark)

*Created 2026-07-10. This folder holds the operational docs for the **pivoted** paper:
a leakage **audit + de-leaked benchmark + resource**, replacing the old orthology-only
framing. The three originals (`../PAPER_SCOPE.md`, `../WORK_SPLIT_3_MACHINES.md`,
`../BUILD_GUIDE.md`, `../SETUP_ENVIRONMENTS.md`) are left intact. If you commit to this
pivot, archive them (move to `archive/`) so there is one source of truth.*

**Paper thesis (tuning-proof):** *De-leaking a multi-species gene–disease benchmark reduces
the measured performance of **every** method — topological heuristics and KGE alike. We
quantify each leakage source and release the de-leaked splits.* The orthology finding rides
along as a secondary, Monarch-specific result, no longer load-bearing.

Read order: `SETUP_ENVIRONMENTS.md` → `WORK_SPLIT.md` → `BUILD_GUIDE.md`.

---

## Why this pivot (one screen)

The old plan staked the whole paper on one number: does orthology survive a degree control?
`memory/publication-strategy.md` put that at ~25–40% for a real journal. This plan makes the
headline "de-leaking drops everyone" (result-independent), turns the make-or-break test into a
**control** (R2 degree-permutation null), and pre-empts the two killer reviews (under-tuned KGE;
"orthology is just degree"). Trade-off: lower ceiling, and it **drops the hybrid predictor**
(your one invented method). Same venue tier, but now reliably reachable.

---

## What carries over from work already done (units 1–4, 6 of the old build guide)

Nothing built is scrapped. Cuts fall on *future* work that was never started.

| Old unit (built) | Fate here | Notes |
|---|---|---|
| 1 `graph_stats.py` | **Reuse as-is** | Table 1; also run on Hetionet. |
| 2 `build_deleaked_splits.py` | **Reuse ~60%, rewrite regimes** | Keep freeze/hash/manifest/assert; R1/R2/R3 redefined. |
| 3 provenance | **Reuse, trim** | Still scored by resource venues. |
| 4 `lib_eval.py` | **Reuse ~80%, extend** | Keep protocol/metrics/bootstrap; new regime map + Jaccard/PA scorers. |
| 6 `run_kge.py` | **Reuse ~70%, trim** | 4 models → 2 (TransE+ComplEx); filtered subgraph + Colab path. |

**Cut (never built — pure savings):** old Unit 7 hybrid, Unit 8 LOSO, Unit 9 newer-release.
Old Unit 10 literature → shrinks to a light spot-check.

**New work added:** R1 redundancy control, R2 degree-permutation null (seed from
`../scripts/null_model_tests.py`), Jaccard + Preferential-Attachment baselines, Hetionet run.

---

## Regime redefinition (the one thing that changes in code)

| | Old (orthology-load-bearing) | **New (audit)** |
|---|---|---|
| R0 | standard random split | standard random split *(unchanged)* |
| R1 | hub-controlled | **redundancy-controlled** (inverse/duplicate/symmetric edges removed) |
| R2 | orthology-blocked | **degree-controlled** (degree-preserving permutation null) |
| R3 | fully de-leaked | **orthology-blocked** (Monarch only; measured *after* R2) |
| stretch | — | temporal split (two dated Monarch releases) |

The R2→R3 residual (orthology drop *after* degree control) carries the secondary novelty. A
null there is now a sentence, not a failure.
