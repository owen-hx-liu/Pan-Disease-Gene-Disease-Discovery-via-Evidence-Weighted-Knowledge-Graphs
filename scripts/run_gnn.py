#!/usr/bin/env python3
"""
run_gnn.py -- the MESSAGE-PASSING arm of the leakage-aware benchmark (R-GCN).

WHY THIS ARM EXISTS
===================
Not because "KGEs are not state of the art" -- this paper is an evaluation audit, not a
leaderboard, so that objection alone does not land. The reason is internal to our own
Related Work (OUTLINE.md 2.3 / DRAFT.md "Reconciliation flagged for the Discussion"):

    Briere et al. (2025) remove degree information and retrain (their DL2 regime) and find
    NO evidence their models rely on degree. Our R2 degree-preserving null finds degree to
    be THE universal leak. We currently resolve that conflict by ARGUMENT (binary retrain
    vs magnitude decomposition; different graph; different task).

Our thesis is that evaluation leakage is METHOD-CLASS-DEPENDENT. We currently test two
classes: parameter-free topological heuristics and shallow KGE. Message passing is the
untested third class, and the one most likely to behave differently -- a GNN aggregates
over a node's neighbourhood, so its dependence on degree is architectural rather than
incidental. Measuring the R0 -> R2 drop for an R-GCN either extends the "degree is
universal" claim to a third class (strengthening it against DL2) or bounds it (which is
also a real result, honestly reported).

SCOPE (deliberately narrow)
===========================
  * ONE model: R-GCN (Schlichtkrull et al., 2018). No NBFNet -- it costs far more and adds
    nothing to this specific argument.
  * Regimes R0 and R2 only; seeds 42/1/7; sampled-50 protocol only (no full-rank).
    => 6 production cells.

WHY THIS IS A SEPARATE SCRIPT THAT IMPORTS run_kge
==================================================
Every piece of the evaluation contract -- split loading, the sampled-50 type-matched
negative protocol, the fixed rank_seed, cold-start -> SENTINEL -> worst rank, per-edge
reciprocal-rank storage, the run-JSON schema, seed merging -- is reused by IMPORTING
run_kge and calling its functions (``load_regimes``, ``build_subgraph``, ``KGEScorer``,
``evaluate_regime``, ``train_eval_write``, ``merge_results``). Nothing is reimplemented,
so protocol parity holds by construction rather than by inspection.

run_kge.py itself is NOT modified. The frozen KGE arm (results/kge/*.json, the manuscript's
Table 2/3/5 source) is on that code path, and this arm has no business perturbing it. The
one thing an R-GCN needs that ``run_kge._build_model`` does not pass through --
architecture kwargs (num_layers, edge_dropout, ...) -- is supplied by making ``cfg["cls"]``
a ``functools.partial`` of the model class, which ``_build_model`` calls exactly as it
calls a bare class. Zero edits, zero monkey-patching.

Per-run files are written as ``kge_RGCN_<regime>_seed<k>.json`` under a SEPARATE results
directory (default ``results/gnn/``). The ``kge_`` prefix is the harness's per-run record
convention -- it is what ``merge_results`` and ``degree_stratified.py`` glob for -- not a
claim that an R-GCN is a KGE. The separate directory is what guarantees this arm can never
land in ``results/kge/kge_summary.json``.

FEASIBILITY: WHY A SUBGRAPH IS MANDATORY
========================================
PyKEEN's R-GCN re-runs message passing over the WHOLE training graph on every optimizer
step (``RGCNRepresentation`` caches enriched embeddings and ``post_parameter_update``
invalidates the cache once per batch). Cost per epoch is therefore
``n_batches x one_full_propagation``, and peak memory is set by the propagation alone --
measured flat at 5.33-5.36 GB across batch 4,096 / 16,384 / 65,536 at an 11.6%-target
292k-edge subgraph. On this 8.5 GB RTX 5070 Laptop GPU (see FINDINGS_gnn_feasibility.md
for the measured table) the full 5.85M-edge graph is not tractable, and 0.11 already sits
at the memory cliff and goes superlinear. ~5% (--subgraph-frac 0.05) is the operating
point.

One inherited behaviour to know about: ``run_kge.train_model`` reacts to CUDA OOM by
halving the batch and retrying. For a KGE that works; for this R-GCN it does NOT, because
peak memory is set by the propagation and is flat in batch size (measured), while halving
the batch DOUBLES the number of propagations per epoch. If an OOM happens here, lower
--subgraph-frac rather than relying on the retry.

Subsampling is scientifically valid HERE because the claim is WITHIN-model ACROSS regimes
(R0 vs R2), and R0 and R2 are given an identical subgraph procedure and seed -- asserted at
runtime by ``assert_subgraph_parity``, not assumed. What subsampling does NOT license is
comparing these absolute numbers against the full-graph KGE numbers: the negative candidate
pool is drawn from the retained subgraph, so it is a DIFFERENT pool than the frozen KGE runs
used, and the two are not paired edge-for-edge. That is why this script runs a matched
reference model (default DistMult) on the byte-identical subgraph -- see below.

WHY THE REFERENCE MODEL IS DistMult
===================================
PyKEEN's R-GCN is, by construction, a DistMult decoder on top of a message-passing encoder
(``interaction="DistMult"``, the published Schlichtkrull setup). DistMult trained on the
IDENTICAL subgraph is therefore the exact control: same scoring function, same graph, same
negatives, same harness -- differing ONLY by the presence of the message-passing encoder.
It makes the R-GCN's R0 number interpretable (which the convergence gate needs) and makes
"does the message-passing class leak differently?" a controlled comparison instead of a
cross-graph guess. It is cheap (DistMult on 292k edges is seconds per epoch).

CONVERGENCE GATE
================
An undertrained GNN scoring near chance at BOTH R0 and R2 would read as "no degree
reliance" -- a false finding that would damage the paper by appearing to support DL2. So
the R2 delta is GATED (``convergence_gate``): unless every R0 seed clears an absolute MRR
floor, an Hits@10 floor, a real training-loss decrease, and a fraction of the matched
reference model's R0 MRR, this script REFUSES to report an R0->R2 delta and says the
baseline did not converge. Thresholds are flags, and the verdict (with reasons) is written
to ``gnn_gate.json`` next to the results.

DEGREE-STRATIFIED REANALYSIS (Table 5 / Figure 4) COMES FREE
============================================================
``evaluate_regime`` already stores one tie-averaged reciprocal rank per test edge in
test.csv row order, so the degree-stratified reanalysis needs NO retraining. With
``--stratify`` this script invokes ``degree_stratified.py`` against its own results
directory, writing to ``<results-dir>/degree_stratified.json`` and leaving the frozen
production artifact untouched. It is skipped automatically when --limit-test truncated the
test set (the per-edge arrays would not align to test.csv).

Usage
-----
    # wiring check, no GPU, no training (random scorer through the real harness):
    python scripts/run_gnn.py --dry-run --limit-test 200

    # feasibility smoke: one seed, few epochs, small subgraph, ~minutes
    python scripts/run_gnn.py --smoke

    # production (6 cells: R0/R2 x seeds 42/1/7) + matched DistMult reference + stratify
    python scripts/run_gnn.py --epochs 300 --subgraph-frac 0.05 \
        --reference-models DistMult --stratify
"""
from __future__ import annotations

import argparse
import functools
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

try:
    import lib_eval as le
    import run_kge as rk
except ImportError:  # running from repo root
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import lib_eval as le
    import run_kge as rk

log = rk.log
TARGET_REL = rk.TARGET_REL
RESULTS_DIR = rk.RESULTS_DIR.parent / "gnn"

# Chance MRR under the shared sampled-50 protocol (51-way ranking, tie-averaged).
# run_kge --dry-run prints this empirically as ~0.089.
CHANCE_MRR = 0.089
CHANCE_HITS10 = 10.0 / 51.0  # ~0.196

# Full-graph, 300-epoch, sampled-50 R0->R2 MRR drop from the frozen KGE arm
# (data/processed/results/kge/kge_summary.json). Used ONLY as the yardstick for the
# subgraph-validity check below -- never as a target for the GNN to match.
FULLGRAPH_R2_DROP = {"DistMult": -0.3632, "ComplEx": -0.3712,
                     "TransE": -0.2358, "RotatE": -0.2084}


# --------------------------------------------------------------------------- #
# Model config: R-GCN, published recipe, shaped for run_kge._build_model
# --------------------------------------------------------------------------- #
def rgcn_config(num_layers=2, edge_dropout=0.4, self_loop_dropout=0.2,
                decomposition=None, num_bases=None, edge_weighting="symmetric",
                num_negs=16, loss="softplus", lr=None):
    """Return {"RGCN": cfg} in the shape ``run_kge.train_model`` expects.

    Everything architectural is bound into ``cls`` as a ``functools.partial`` so
    ``run_kge._build_model`` -- which only forwards triples_factory / embedding_dim /
    random_seed / loss / regularizer -- can stay untouched.

    Defaults are PyKEEN's, which ARE the published R-GCN link-prediction recipe
    (Schlichtkrull et al., 2018): 2 layers, ReLU (PyKEEN resolves ``activation=None`` to
    ReLU, and suppresses it on the last layer as in the reference implementation), basis
    decomposition of the relation-specific weight matrices, symmetric ``1/c_{i,r}`` edge
    normalisation, and a DistMult decoder.

    Two deliberate choices, both for parity with the existing KGE arm rather than tuning:

    * ``loss="softplus"`` instead of PyKEEN's RGCN default (MarginRankingLoss). The decoder
      IS DistMult, and the bilinear models in ``run_kge.model_config()`` use the
      logistic/softplus objective (Yang et al., 2015; Trouillon et al., 2016); the
      published R-GCN also trains with a sigmoid cross-entropy against sampled negatives.
      Using softplus keeps R-GCN and its DistMult control on the same objective, so the
      comparison isolates the encoder.
    * 16 negatives per positive under SLCWA -- identical to every model in the KGE arm.

    ``regularizer`` is deliberately NOT set: unlike DistMult and ComplEx, PyKEEN attaches
    no default regularizer to RGCN (``RGCN.regularizer_default is None``, verified), so the
    Lp-collapse trap documented in ``run_kge.model_config()`` does not apply here and there
    is nothing to correct.
    """
    from pykeen.models import RGCN

    arch = dict(num_layers=num_layers, edge_dropout=edge_dropout,
                self_loop_dropout=self_loop_dropout, edge_weighting=edge_weighting)
    if decomposition:
        arch["decomposition"] = decomposition
        if num_bases:
            arch["decomposition_kwargs"] = dict(num_bases=num_bases)
    cfg = dict(
        cls=functools.partial(RGCN, **arch),
        num_negs=num_negs,
        loss=loss,
        ref="Schlichtkrull+2018 (R-GCN); DistMult decoder; softplus loss for parity "
            "with the bilinear KGE arm",
    )
    if lr is not None:
        cfg["lr"] = lr
    # Recorded verbatim into every run JSON (a partial's repr is useless for provenance).
    cfg["_gnn_meta"] = dict(
        method_class="message_passing",
        encoder="RGCN", decoder="DistMult (PyKEEN RGCN interaction default)",
        num_layers=num_layers, activation="relu (PyKEEN default; none on last layer)",
        edge_dropout=edge_dropout, self_loop_dropout=self_loop_dropout,
        edge_weighting=edge_weighting,
        decomposition=(decomposition or "bases (PyKEEN default)"),
        num_bases=num_bases, loss=loss, num_negs=num_negs,
        propagation="full-graph over the retained subgraph, recomputed every optimizer step",
        reference="Schlichtkrull et al., ESWC 2018",
    )
    return {"RGCN": cfg}


# --------------------------------------------------------------------------- #
# PARITY: assert, do not assume. Any deviation is logged LOUDLY.
# --------------------------------------------------------------------------- #
def _degree_map(edges):
    import collections
    c = collections.Counter(edges[:, 0].tolist())
    c.update(edges[:, 2].tolist())
    return c


def _degree_preservation(regimes, dev) -> dict:
    """Measure how well the retained subgraph preserves the R2 null's defining property.

    Reported on three levels, because they behave differently and only the third is
    degraded by subsampling:
      * task relation, per node -- must be EXACT (every target edge is force-kept);
      * whole subgraph -- approximate, quantified by Pearson r / mean |delta| over the
        entities present in both regimes;
      * the test genes and test diseases -- the entities the ranking protocol actually
        exposes; the ranked-disease axis is the one Table 5 shows carries the leak.
    """
    out: dict = {}
    base = regimes[0]
    dm = {rd["short"]: _degree_map(rd["train_edges"]) for rd in regimes}
    tm = {rd["short"]: _degree_map(rd["train_edges"][rd["train_edges"][:, 1] == TARGET_REL])
          for rd in regimes}
    b = base["short"]
    for rd in regimes[1:]:
        n = rd["short"]
        # (2) exact on the task relation
        exact = set(tm[b]) == set(tm[n]) and all(tm[b][k] == tm[n][k] for k in tm[b])
        maxd = (0 if exact else
                max(abs(tm[b].get(k, 0) - tm[n].get(k, 0)) for k in set(tm[b]) | set(tm[n])))
        out[f"target_relation_{b}_vs_{n}"] = {
            "n_nodes": len(tm[b]), "per_node_degree_exact": bool(exact),
            "max_abs_delta": int(maxd)}
        if not exact:
            dev.append(f"CRITICAL: target-relation per-node degree is NOT preserved between "
                       f"{b} and {n} (max |delta| {maxd}). The R2 null is supposed to hold "
                       f"degree fixed exactly on the task relation; this invalidates the "
                       f"degree-vs-structure attribution.")
        # (3) approximate on the whole subgraph + on the evaluated entities
        block = {}
        shared = sorted(set(dm[b]) & set(dm[n]))
        x = np.array([dm[b][k] for k in shared], float)
        y = np.array([dm[n][k] for k in shared], float)
        block["whole_subgraph"] = {
            "n_entities": {b: len(dm[b]), n: len(dm[n])}, "n_shared": len(shared),
            "pearson_r": round(float(np.corrcoef(x, y)[0, 1]), 6) if len(shared) > 1 else None,
            "mean_abs_delta": round(float(np.abs(x - y).mean()), 4),
            "max_abs_delta": int(np.abs(x - y).max()) if len(shared) else 0,
            "mean_degree": {b: round(float(x.mean()), 4), n: round(float(y.mean()), 4)},
        }
        te = base["test_edges"]
        for lbl, col in (("test_gene_query", 0), ("test_disease_ranked_target", 2)):
            ids = np.unique(te[:, col])
            u = np.array([dm[b].get(i, 0) for i in ids], float)
            v = np.array([dm[n].get(i, 0) for i in ids], float)
            block[lbl] = {
                "n": int(len(ids)),
                "pearson_r": round(float(np.corrcoef(u, v)[0, 1]), 6) if len(ids) > 1 else None,
                "mean_abs_delta": round(float(np.abs(u - v).mean()), 4),
                "mean_degree": {b: round(float(u.mean()), 3), n: round(float(v.mean()), 3)},
            }
        out[f"context_graph_{b}_vs_{n}"] = block
        r_dis = block["test_disease_ranked_target"]["pearson_r"]
        if r_dis is not None and r_dis < 0.99:
            dev.append(f"subgraph degree preservation on the RANKED DISEASE axis has "
                       f"degraded ({b} vs {n} Pearson r={r_dis:.4f} < 0.99). This is the "
                       f"axis the degree leak runs on, so part of any R0->R2 drop could be "
                       f"a degree-distribution change rather than structure loss. Raise "
                       f"--subgraph-frac.")
        else:
            dev.append(f"subgraph degree preservation ({b} vs {n}): EXACT per node on the "
                       f"task relation; approximate on the subsampled context graph "
                       f"(whole-graph Pearson r="
                       f"{out[f'context_graph_{b}_vs_{n}']['whole_subgraph']['pearson_r']}, "
                       f"ranked-disease axis r={r_dis}). Expected consequence of edge "
                       f"subsampling, not of the null; quantified in every run JSON.")
    return out


def assert_subgraph_parity(regimes, args) -> dict:
    """Verify that every regime got the SAME subgraph procedure, and quantify what the
    subgraph did to the evaluation. Returns a parity block for the run JSONs.

    The within-model / across-regime comparison is only valid if R0 and R2 were handed
    graphs that differ ONLY by the degree-preserving rewiring -- not by which edges the
    subsampler happened to keep. ``run_kge.build_subgraph`` selects by ROW INDEX with a
    fixed seed after force-keeping every target edge, so identical (frac, seed) on two
    files with the same row count selects the same index set. That is checked here rather
    than trusted, by three invariants:

    1. RELATION COMPOSITION. ``build_degree_null.py`` rewires WITHIN each relation type and
       writes its output in the input's row order (verified: the relation column is
       identical row-for-row between train.csv and train_R2_degree_null_seed42.csv), so
       index-based sampling must yield identical per-relation edge counts. If it does not,
       the regimes differ in graph COMPOSITION and the comparison is confounded.
       NOTE: the retained edges' ENDPOINTS differ almost completely between R0 and R2
       (Jaccard ~0.002). That is the null working as designed -- it rewires every relation,
       not just the target one -- and is not a parity failure.
    2. DEGREE PRESERVATION ON THE TASK RELATION. Every target (gene->disease) edge is
       force-kept in both regimes and the rewiring is degree-preserving within a relation,
       so per-node target-relation degree must be preserved EXACTLY (measured: max |delta|
       = 0 over 17,139 nodes at frac 0.05).
    3. DEGREE PRESERVATION ON THE CONTEXT GRAPH. This one is only APPROXIMATE and it is a
       genuine cost of subsampling: an edge sample of a rewired graph does not reproduce
       the degree sequence of the same sample of the original, because the rewiring moved
       which edges attach where. Quantified rather than hidden (measured at frac 0.05:
       Pearson r = 0.990, mean |delta| = 0.9, identical mean degree; for the ranked test
       diseases -- the axis the sampled-negative protocol actually exposes -- r = 0.9999,
       mean |delta| = 0.44). Reported so a reader can judge it; flagged as a deviation if
       it degrades.

    Deviations from the frozen full-graph KGE protocol are collected and returned so they
    end up in every run JSON and in the console log. The important one: the negative
    candidate pool is built from the RETAINED subgraph, so it is a different pool than the
    frozen KGE runs used. Within-subgraph comparisons are paired and exact; comparisons
    against results/kge/ absolute numbers are trend-level only.
    """
    dev: list[str] = []
    info: dict = {}

    fracs = {rd["short"]: rd["subgraph_frac"] for rd in regimes}
    if len(set(fracs.values())) > 1:
        dev.append(f"CRITICAL: regimes got DIFFERENT subgraph fractions {fracs} -- the "
                   f"across-regime comparison is invalid")
    tests = {rd["short"]: rd["test_file"] for rd in regimes}
    if len(set(tests.values())) > 1:
        dev.append(f"regimes use different test files {tests}")
    info["subgraph_frac"] = fracs
    info["subgraph_seed"] = args.subgraph_seed
    info["test_files"] = tests
    info["n_train_edges"] = {rd["short"]: int(len(rd["train_edges"])) for rd in regimes}
    info["n_train_full"] = {rd["short"]: int(rd["n_train_full"]) for rd in regimes}

    # Row-count parity of the SOURCE graphs: index-based subsampling only lines up if the
    # regime files have the same number of rows.
    if len({rd["n_train_full"] for rd in regimes}) > 1:
        dev.append(f"regime source graphs have different edge counts "
                   f"{info['n_train_full']} -- index-based subsampling does NOT select the "
                   f"same non-target edges across regimes")

    # Target-edge counts per retained subgraph: build_subgraph force-keeps all of them, so
    # these must agree across regimes (R2 rewires target edges, it does not add/remove).
    n_tgt = {rd["short"]: int((rd["train_edges"][:, 1] == TARGET_REL).sum()) for rd in regimes}
    info["n_target_edges_kept"] = n_tgt
    if len(set(n_tgt.values())) > 1:
        dev.append(f"retained target-edge counts differ across regimes {n_tgt}")

    # Invariant 1: identical per-relation edge counts => the same ROWS were sampled.
    # (Endpoint differences are the null doing its job; see the docstring.)
    if len(regimes) > 1:
        import collections
        comps = {rd["short"]: collections.Counter(rd["train_edges"][:, 1].tolist())
                 for rd in regimes}
        base_name, base = next(iter(comps.items()))
        for name, c in comps.items():
            if name == base_name or c == base:
                continue
            diff = {k: (base.get(k, 0), c.get(k, 0)) for k in set(base) | set(c)
                    if base.get(k, 0) != c.get(k, 0)}
            dev.append(f"CRITICAL: retained per-relation edge counts differ between "
                       f"{base_name} and {name} ({len(diff)} relation(s), e.g. "
                       f"{dict(list(diff.items())[:3])}). The subsampler did NOT select the "
                       f"same rows, so the regimes differ in graph composition and the "
                       f"across-regime comparison is confounded.")
        info["relation_composition_identical_across_regimes"] = all(
            c == base for c in comps.values())
        info["n_relations_retained"] = {k: len(v) for k, v in comps.items()}

        # Invariants 2 & 3: how well the SUBGRAPH preserves degree. Exact on the task
        # relation (all its edges are kept); approximate on the subsampled context graph.
        info["degree_preservation"] = _degree_preservation(regimes, dev)

    # Protocol constants that must match the frozen KGE runs exactly.
    if args.n_neg != 50:
        dev.append(f"n_neg={args.n_neg} != 50 (frozen protocol is sampled-50)")
    if args.rank_seed != 42:
        dev.append(f"rank_seed={args.rank_seed} != 42 (frozen runs used 42; per-edge ranks "
                   f"will not be paired-comparable)")
    if args.limit_test:
        dev.append(f"--limit-test {args.limit_test}: test set TRUNCATED; per-edge ranks do "
                   f"not align to test.csv and the degree-stratified reanalysis is skipped")
    if args.full_rank:
        dev.append("--full-rank set; the GNN arm is specified as sampled-50 only")

    # Cold start: the subgraph makes far more test entities out-of-vocabulary than the
    # 30 genes of the full-graph runs. Quantified per regime, never absorbed silently.
    for rd in regimes:
        vocab = set(np.unique(np.concatenate([rd["train_edges"][:, 0],
                                              rd["train_edges"][:, 2]])).tolist())
        te = rd["test_edges"]
        n_g = int(sum(1 for g in np.unique(te[:, 0]) if g not in vocab))
        n_d = int(sum(1 for d in np.unique(te[:, 2]) if d not in vocab))
        n_edge_cold = int(sum(1 for g, _r, d in te if g not in vocab or d not in vocab))
        info[f"cold_start_{rd['short']}"] = {
            "n_test_genes_oov": n_g, "n_test_genes": int(len(np.unique(te[:, 0]))),
            "n_test_diseases_oov": n_d, "n_test_diseases": int(len(np.unique(te[:, 2]))),
            "n_test_edges_touching_oov": n_edge_cold, "n_test_edges": int(len(te)),
            "frac_test_edges_touching_oov": round(n_edge_cold / max(len(te), 1), 4),
        }
    cold = [info[f"cold_start_{rd['short']}"]["frac_test_edges_touching_oov"] for rd in regimes]
    if max(cold, default=0) > 0.02:
        dev.append(f"subgraph cold start: up to {100*max(cold):.1f}% of test edges touch an "
                   f"out-of-vocabulary entity (full-graph KGE runs: 30 genes). Handled by "
                   f"the SAME SENTINEL -> worst-rank rule, but it depresses absolute MRR, "
                   f"so absolute numbers are NOT comparable to results/kge/")

    dev.append("negative candidate pool is drawn from the RETAINED subgraph, so it differs "
               "from the pool the frozen full-graph KGE runs used: within-subgraph "
               "comparisons (RGCN R0 vs R2, RGCN vs the matched reference model) are "
               "paired and protocol-identical; comparisons against results/kge/ absolute "
               "numbers are trend-level only")

    info["deviations"] = dev
    info["harness"] = "scripts/lib_eval.py via scripts/run_kge.py (imported, not reimplemented)"
    info["protocol"] = f"sampled-{args.n_neg}neg type-matched, rank_seed={args.rank_seed}"

    log("=== PROTOCOL PARITY ===")
    log(f"  harness      : {info['harness']}")
    log(f"  protocol     : {info['protocol']}")
    log(f"  subgraph     : frac={fracs} seed={args.subgraph_seed}")
    log(f"  train edges  : {info['n_train_edges']} (of {info['n_train_full']})")
    log(f"  target edges : {n_tgt} (all force-kept)")
    for rd in regimes:
        c = info[f"cold_start_{rd['short']}"]
        log(f"  cold start {rd['short']}: {c['n_test_genes_oov']}/{c['n_test_genes']} test genes, "
            f"{c['n_test_diseases_oov']}/{c['n_test_diseases']} test diseases OOV -> "
            f"{c['n_test_edges_touching_oov']}/{c['n_test_edges']} test edges "
            f"({100*c['frac_test_edges_touching_oov']:.1f}%) touch an OOV entity")
    if "relation_composition_identical_across_regimes" in info:
        log(f"  relation composition identical across regimes: "
            f"{info['relation_composition_identical_across_regimes']} "
            f"(=> the subsampler selected the same rows)")
    for k, v in (info.get("degree_preservation") or {}).items():
        if k.startswith("target_relation_"):
            log(f"  degree, task relation ({k[len('target_relation_'):]}): "
                f"per-node EXACT={v['per_node_degree_exact']} over {v['n_nodes']:,} nodes "
                f"(max |delta| {v['max_abs_delta']})")
        else:
            w = v["whole_subgraph"]
            d = v["test_disease_ranked_target"]
            log(f"  degree, context graph ({k[len('context_graph_'):]}): whole-subgraph "
                f"r={w['pearson_r']} mean|delta|={w['mean_abs_delta']} "
                f"mean deg {w['mean_degree']}; ranked-disease axis r={d['pearson_r']} "
                f"mean|delta|={d['mean_abs_delta']}")
    if dev:
        log(f"  !! {len(dev)} PARITY NOTE(S)/DEVIATION(S) -- recorded in every run JSON:")
        for d in dev:
            log(f"     - {d}")
    else:
        log("  no deviations")
    return info


# --------------------------------------------------------------------------- #
# CONVERGENCE GATE
# --------------------------------------------------------------------------- #
def _load_runs(results_dir, model):
    out = {}
    for p in sorted(Path(results_dir).glob(f"kge_{model}_*_seed*.json")):
        try:
            r = json.loads(p.read_text())
        except Exception:
            continue
        if r.get("model") == model and r.get("regime"):
            out.setdefault(r["regime"], []).append(r)
    return out


def convergence_gate(results_dir, args, ref_r0_mrr=None) -> dict:
    """Decide whether the R-GCN baseline converged well enough for an R2 delta to MEAN
    anything. Returns a verdict dict; also written to ``gnn_gate.json``.

    The failure mode this exists to catch: a model that never learned scores ~chance at R0
    AND at R2, which reads as "the GNN shows no degree reliance" and would look like
    support for Briere's DL2. That would be a false finding, so the delta is withheld
    unless the R0 baseline is demonstrably a working model.

    Four checks, all on the R0 cells, all per seed:
      1. training loss actually decreased (catches the flat-loss collapse mode that the
         bilinear models exhibited: loss pinned at its zero-gap value, MRR ~ chance);
      2. MRR clears an absolute floor well above the ~0.089 chance level;
      3. Hits@10 clears a floor well above its ~0.196 chance level;
      4. MRR reaches a fraction of the MATCHED-SUBGRAPH reference model's R0 MRR -- the
         only fair "near our KGE numbers" test, since the frozen KGE numbers were measured
         on the full graph with a different negative pool.
    Check 4 is skipped (and said to be skipped) when no reference model was run.
    """
    runs = _load_runs(results_dir, "RGCN").get("R0", [])
    checks, per_seed = [], []
    if not runs:
        verdict = {"verdict": "fail", "reason": "no RGCN R0 runs found", "checks": [],
                   "n_r0_runs": 0}
    else:
        for r in sorted(runs, key=lambda x: x.get("seed", 0)):
            seed = r.get("seed")
            losses = r.get("epoch_losses") or []
            c: list[dict] = []
            if len(losses) >= 2:
                drop_ok = losses[-1] <= args.gate_loss_drop * losses[0]
                c.append({"check": "training loss decreased",
                          "detail": f"loss {losses[0]:.4f} -> {losses[-1]:.4f} "
                                    f"(need <= {args.gate_loss_drop:g} x first)",
                          "pass": bool(drop_ok)})
            else:
                c.append({"check": "training loss decreased",
                          "detail": f"only {len(losses)} epoch loss(es) recorded",
                          "pass": len(losses) > 0})
            c.append({"check": "MRR above absolute floor",
                      "detail": f"MRR {r['MRR']:.4f} vs floor {args.gate_min_mrr:g} "
                                f"(chance ~{CHANCE_MRR})",
                      "pass": bool(r["MRR"] >= args.gate_min_mrr)})
            c.append({"check": "Hits@10 above absolute floor",
                      "detail": f"Hits@10 {r['Hits@10']:.4f} vs floor "
                                f"{args.gate_min_hits10:g} (chance ~{CHANCE_HITS10:.3f})",
                      "pass": bool(r["Hits@10"] >= args.gate_min_hits10)})
            if ref_r0_mrr is not None and ref_r0_mrr > 0:
                need = args.gate_ref_ratio * ref_r0_mrr
                c.append({"check": "MRR vs matched-subgraph reference",
                          "detail": f"MRR {r['MRR']:.4f} vs {args.gate_ref_ratio:g} x "
                                    f"reference {ref_r0_mrr:.4f} = {need:.4f}",
                          "pass": bool(r["MRR"] >= need)})
            else:
                c.append({"check": "MRR vs matched-subgraph reference",
                          "detail": "SKIPPED: no matched reference model run "
                                    "(--reference-models). The 'near our KGE numbers' "
                                    "criterion cannot be evaluated against full-graph "
                                    "KGE numbers, which used a different negative pool.",
                          "pass": None})
            per_seed.append({"seed": seed, "MRR": r["MRR"], "Hits@10": r["Hits@10"],
                             "epochs": r.get("epochs"), "checks": c,
                             "pass": all(x["pass"] for x in c if x["pass"] is not None)})
            checks.extend(c)
        failed = [f"seed {s['seed']}: " + "; ".join(
            f"{x['check']} ({x['detail']})" for x in s["checks"] if x["pass"] is False)
            for s in per_seed if not s["pass"]]
        verdict = {
            "verdict": "pass" if not failed else "fail",
            "n_r0_runs": len(per_seed),
            "reference_r0_mrr": ref_r0_mrr,
            "thresholds": {"min_MRR": args.gate_min_mrr, "min_Hits@10": args.gate_min_hits10,
                           "loss_drop_factor": args.gate_loss_drop,
                           "ref_ratio": args.gate_ref_ratio},
            "chance": {"MRR": CHANCE_MRR, "Hits@10": round(CHANCE_HITS10, 4)},
            "per_seed": per_seed,
            "failures": failed,
            "reason": ("all R0 seeds cleared every applicable check" if not failed
                       else "; ".join(failed)),
        }
    verdict["generated"] = time.strftime("%Y-%m-%d %H:%M:%S")
    dest = Path(results_dir) / "gnn_gate.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(verdict, indent=2))

    log("=== CONVERGENCE GATE ===")
    for s in per_seed:
        log(f"  seed {s['seed']}: MRR={s['MRR']:.4f} Hits@10={s['Hits@10']:.4f} "
            f"-> {'PASS' if s['pass'] else 'FAIL'}")
        for x in s["checks"]:
            mark = {True: "ok  ", False: "FAIL", None: "skip"}[x["pass"]]
            log(f"      [{mark}] {x['check']}: {x['detail']}")
    log(f"  VERDICT: {verdict['verdict'].upper()} -> {dest.name}")
    return verdict


def report_regime_delta(results_dir, gate, models) -> dict:
    """R0 -> R2 MRR delta per model -- the number this arm exists to produce.

    Written ONLY if the convergence gate passed. On failure the payload records
    ``reported: false`` plus the gate's reason, so a failed run leaves an explicit
    "did not converge" artifact rather than a delta that could be mistaken for evidence
    that message passing does not rely on degree.
    """
    out = {"gate": gate["verdict"], "gate_reason": gate["reason"], "models": {},
           "generated": time.strftime("%Y-%m-%d %H:%M:%S")}
    if gate["verdict"] != "pass":
        out["reported"] = False
        out["note"] = ("R0 -> R2 delta WITHHELD: the R-GCN baseline did not converge, so a "
                       "small delta cannot be distinguished from a model that never "
                       "learned. Do NOT read this run as evidence about degree reliance.")
        Path(results_dir, "gnn_regime_delta.json").write_text(json.dumps(out, indent=2))
        log("=" * 78)
        log("!! BASELINE DID NOT CONVERGE -- REFUSING TO REPORT AN R0->R2 DELTA")
        log(f"!! {gate['reason']}")
        log("!! A near-chance model scores near chance under BOTH regimes; that is not a")
        log("!! finding about degree reliance. Train longer / raise --subgraph-frac / see")
        log("!! FINDINGS_gnn_feasibility.md before drawing any conclusion.")
        log("=" * 78)
        return out

    out["reported"] = True
    out["subgraph_validity"] = None
    for model in models:
        by = _load_runs(results_dir, model)
        cell = {}
        for reg in ("R0", "R2"):
            rs = by.get(reg, [])
            if not rs:
                continue
            m = np.array([r["MRR"] for r in rs], float)
            cell[reg] = {"MRR_mean": float(m.mean()),
                         "MRR_sd": float(m.std(ddof=1)) if len(m) > 1 else 0.0,
                         "n_seeds": len(m), "seeds": sorted(r["seed"] for r in rs)}
        if "R0" in cell and "R2" in cell:
            a, b = cell["R0"]["MRR_mean"], cell["R2"]["MRR_mean"]
            cell["delta_MRR"] = float(b - a)
            cell["pct_change"] = float(100.0 * (b - a) / a) if a else float("nan")
        out["models"][model] = cell

    out["subgraph_validity"] = _subgraph_validity(out, args_min_ratio=0.5)

    Path(results_dir, "gnn_regime_delta.json").write_text(json.dumps(out, indent=2))
    log("=== R0 -> R2 (degree-preserving null) on the MATCHED subgraph ===")
    for model, cell in out["models"].items():
        if "delta_MRR" in cell:
            log(f"  {model:10s} R0 MRR={cell['R0']['MRR_mean']:.4f} "
                f"-> R2 MRR={cell['R2']['MRR_mean']:.4f}  "
                f"({cell['pct_change']:+.1f}%)")
        else:
            log(f"  {model:10s} incomplete (need both R0 and R2): {sorted(cell)}")
    sv = out["subgraph_validity"]
    if sv and sv.get("checked"):
        log("=== SUBGRAPH VALIDITY (does the control still show the paper's leak?) ===")
        for m, r in sv["reference_models"].items():
            log(f"  {m}: subgraph R2 drop {100*r['subgraph_drop']:+.1f}% vs full-graph "
                f"{100*r['fullgraph_drop']:+.1f}% -> retained "
                f"{100*r['retained_fraction']:.0f}% of the known collapse "
                f"[{'OK' if r['pass'] else 'ATTENUATED'}]")
        if not sv["pass"]:
            log("=" * 78)
            log("!! SUBGRAPH DOES NOT REPRODUCE THE PAPER'S DEGREE COLLAPSE")
            log("!! The matched reference model's R0->R2 drop is far smaller here than on")
            log("!! the full graph. The subgraph regime is therefore NOT the regime the R2")
            log("!! result was established in, so the RGCN delta above CANNOT be read as a")
            log("!! statement about degree leakage in this paper's sense -- with a weak leak")
            log("!! present there is little for a method class to depend on either way.")
            log("!! Likely causes, in order: (a) undertraining (the full-graph reference ran")
            log("!! 300 epochs); (b) build_subgraph force-keeps ALL target edges, raising the")
            log("!! task relation from 0.6% to ~11.6% of edges (~20x enrichment) and cutting")
            log("!! mean degree ~5x, which removes much of the popularity signal the leak")
            log("!! rides on. Fix before production: raise --epochs and --subgraph-frac until")
            log("!! the reference recovers its collapse. See FINDINGS_gnn_feasibility.md.")
            log("=" * 78)
    return out


def _subgraph_validity(delta_payload, args_min_ratio=0.5) -> dict:
    """Is the retained subgraph still a graph in which the paper's R2 leak EXISTS?

    A second failure mode, distinct from the convergence gate and just as capable of
    manufacturing a false 'message passing does not rely on degree' result: if the subgraph
    (or an undertrained model) has already destroyed most of the degree leak, then EVERY
    method will show a small R0->R2 drop, and the R-GCN's small drop says nothing about
    method classes. The matched reference model is the instrument: its full-graph R2 drop is
    known and frozen, so comparing its SUBGRAPH drop against that value measures how much
    of the phenomenon survived the subsample.
    """
    refs = {}
    for model, cell in delta_payload.get("models", {}).items():
        known = FULLGRAPH_R2_DROP.get(model)
        if known is None or "pct_change" not in cell:
            continue  # RGCN has no full-graph counterpart: it is the thing being measured
        sub = cell["pct_change"] / 100.0
        retained = sub / known if known else float("nan")
        refs[model] = {"subgraph_drop": round(sub, 4), "fullgraph_drop": known,
                       "retained_fraction": round(retained, 4),
                       "pass": bool(retained >= args_min_ratio)}
    if not refs:
        return {"checked": False, "pass": None,
                "reason": "no reference model with a known full-graph R2 drop was run; "
                          "pass --reference-models DistMult to enable this check"}
    return {"checked": True, "min_retained_fraction": args_min_ratio,
            "pass": all(r["pass"] for r in refs.values()),
            "reference_models": refs,
            "interpretation": ("The reference model's R0->R2 drop on this subgraph, as a "
                               "fraction of its frozen full-graph drop. Well below 1.0 means "
                               "the subgraph regime has attenuated the very leak under "
                               "study, so no method's delta measured here transfers to the "
                               "paper's full-graph claim.")}


# --------------------------------------------------------------------------- #
# Degree-stratified reanalysis (free: per-edge ranks are already stored)
# --------------------------------------------------------------------------- #
def run_stratify(results_dir, models, regimes, seeds, parity) -> bool:
    """Invoke degree_stratified.py against THIS arm's results directory.

    Requires the per-edge arrays to align to test.csv (i.e. no --limit-test), which the
    caller checks. Output goes to <results_dir>/degree_stratified.json so the frozen
    production artifact under results/degree_stratified/ is never touched.
    """
    script = Path(__file__).resolve().parent / "degree_stratified.py"
    out = Path(results_dir) / "degree_stratified.json"
    cmd = [sys.executable, str(script),
           "--kge-dir", str(results_dir),
           "--models", *models,
           "--regimes", *regimes,
           "--seeds", *[str(s) for s in seeds],
           "--out", str(out),
           "--no-baselines"]
    log(f"=== degree-stratified reanalysis (no retraining) ===")
    log(f"  $ {' '.join(cmd[1:])}")
    r = subprocess.run(cmd, capture_output=True, text=True)
    for line in (r.stdout or "").rstrip().splitlines():
        log(f"  | {line}")
    if r.returncode != 0:
        log(f"  WARNING degree_stratified.py exited {r.returncode}")
        for line in (r.stderr or "").rstrip().splitlines()[-15:]:
            log(f"  ! {line}")
        return False
    # Stamp the subgraph caveat into the artifact: the stratification axis is FULL-graph R0
    # degree (identical bins to the frozen Table 5, which is the point), while this model
    # only ever saw the subgraph.
    try:
        payload = json.loads(out.read_text())
        payload.setdefault("meta", {})["gnn_subgraph_caveat"] = (
            "Stratification axis is node degree in the FULL R0 training graph "
            "(splits/node_degree.csv), deliberately identical to the frozen Table 5 bins so "
            "the message-passing column lands in the same strata as the KGE columns. This "
            "model was trained on a subgraph, so the degree it actually SAW is smaller; the "
            "bins are a measure of true entity popularity, not of the model's input.")
        payload["meta"]["gnn_parity"] = parity
        out.write_text(json.dumps(payload, indent=2))
    except Exception as e:
        log(f"  WARNING could not stamp caveat into {out.name}: {type(e).__name__}: {e}")
    return True


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--splits-dir", default=str(rk.SPLITS_DIR))
    ap.add_argument("--results-dir", default=str(RESULTS_DIR))
    ap.add_argument("--regimes", nargs="+", default=["R0", "R2"],
                    help="R0 and R2 only: this arm tests the degree null, not the full grid")
    ap.add_argument("--seeds", nargs="+", type=int, default=[42, 1, 7])
    ap.add_argument("--dim", type=int, default=64, help="embedding dim (KGE arm parity)")
    ap.add_argument("--epochs", type=int, default=300)
    ap.add_argument("--reference-models", nargs="*", default=[],
                    help="KGE models to train on the IDENTICAL subgraph as controls. "
                         "DistMult is the exact control: it is the R-GCN's own decoder.")

    # Architecture
    ap.add_argument("--num-layers", type=int, default=2)
    ap.add_argument("--edge-dropout", type=float, default=0.4)
    ap.add_argument("--self-loop-dropout", type=float, default=0.2)
    ap.add_argument("--decomposition", default=None, choices=["bases", "block"],
                    help="default None = PyKEEN's basis decomposition (the published recipe)")
    ap.add_argument("--num-bases", type=int, default=None)
    ap.add_argument("--edge-weighting", default="symmetric")
    ap.add_argument("--num-negs", type=int, default=16)
    ap.add_argument("--loss", default="softplus")

    # Feasibility knobs
    ap.add_argument("--subgraph-frac", type=float, default=0.05,
                    help="REQUIRED in practice: PyKEEN's RGCN propagates over the whole "
                         "training graph every optimizer step, so the full 5.85M-edge "
                         "graph does not fit this GPU. Identical for every regime + seed.")
    ap.add_argument("--subgraph-seed", type=int, default=42)
    ap.add_argument("--slcwa-batch", type=int, default=16384,
                    help="per-epoch cost is (n_batches x one full-graph propagation), so a "
                         "LARGER batch means fewer propagations per epoch; peak memory is "
                         "set by the propagation and is flat in this value")
    ap.add_argument("--lr", type=float, default=None)

    # Protocol (defaults MUST match the frozen KGE runs)
    ap.add_argument("--n-neg", type=int, default=50)
    ap.add_argument("--rank-seed", type=int, default=42)
    ap.add_argument("--class-cap", type=int, default=None)
    ap.add_argument("--limit-test", type=int, default=None)
    ap.add_argument("--train-fit-n", type=int, default=0,
                    help="also score N sampled TRAINING edges (optimization-vs-"
                         "generalization diagnostic; see run_kge.sample_train_fit_edges)")
    ap.add_argument("--train-fit-seed", type=int, default=42)

    # Convergence gate
    ap.add_argument("--gate-min-mrr", type=float, default=0.35,
                    help="R0 MRR floor. Chance is ~0.089; the full-graph KGE R0 band is "
                         "0.669-0.827. 0.35 is ~4x chance and roughly the midpoint to the "
                         "weakest KGE, i.e. unambiguously a working model while not "
                         "demanding full-graph performance from a 5%% subgraph.")
    ap.add_argument("--gate-min-hits10", type=float, default=0.50,
                    help="R0 Hits@10 floor (chance ~0.196)")
    ap.add_argument("--gate-loss-drop", type=float, default=0.9,
                    help="require last epoch loss <= this x first epoch loss")
    ap.add_argument("--gate-ref-ratio", type=float, default=0.60,
                    help="require R0 MRR >= this x the matched-subgraph reference model's "
                         "R0 MRR (skipped if no reference model was run)")

    ap.add_argument("--stratify", action="store_true",
                    help="run the degree-stratified reanalysis afterwards (free: reads the "
                         "per-edge ranks this run already stored)")
    ap.add_argument("--smoke", action="store_true",
                    help="feasibility smoke preset: 1 seed, 3 epochs, frac 0.02, "
                         "results/gnn_smoke/, DistMult reference, relaxed gate")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--dry-run", action="store_true",
                    help="validate harness wiring with a random scorer; no GPU, no training")
    ap.add_argument("--progress-secs", type=float, default=2.0)

    # Accepted so this namespace can be handed straight to run_kge.train_eval_write.
    ap.add_argument("--full-rank", action="store_true",
                    help="NOT part of this arm's scope (sampled-50 only); flagged as a "
                         "parity deviation if set")
    ap.add_argument("--lcwa-slice-size", type=int, default=20000)
    ap.add_argument("--validate-fullrank", type=int, default=0)
    args = ap.parse_args()

    if args.smoke:
        args.seeds = [42]
        args.epochs = 3 if args.epochs == 300 else args.epochs
        args.subgraph_frac = 0.02 if args.subgraph_frac == 0.05 else args.subgraph_frac
        if Path(args.results_dir) == RESULTS_DIR:
            args.results_dir = str(RESULTS_DIR.parent / "gnn_smoke")
        if not args.reference_models:
            args.reference_models = ["DistMult"]
        args.gate_min_mrr = min(args.gate_min_mrr, 0.20)
        log("--smoke: 1 seed, few epochs, small subgraph, relaxed gate, "
            f"-> {args.results_dir}. NOT a production run.")

    os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")
    results_dir = Path(args.results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    regimes = rk.load_regimes(args.regimes, args.splits_dir,
                              args.subgraph_frac, args.subgraph_seed)
    if not regimes:
        log("no regimes available (training files missing). nothing to do.")
        return
    if args.limit_test:
        for rd in regimes:
            rd["test_edges"] = rd["test_edges"][:args.limit_test]

    parity = assert_subgraph_parity(regimes, args)

    # ----- dry run: exercise the real harness with a random scorer ----- #
    if args.dry_run:
        rng = np.random.default_rng(0)
        for rd in regimes:
            ranks = le.rank_test_edges(lambda g, d: float(rng.random()),
                                       rd["train_edges"], rd["test_edges"],
                                       rd["hub_set"], rd["hub_filter"],
                                       n_neg=args.n_neg, seed=args.rank_seed)
            rm = le.ranking_metrics(ranks)
            log(f"[dry-run] {rd['short']}: n_test={rm['n']} MRR={rm['MRR']:.4f} "
                f"(chance ~{CHANCE_MRR}) Hits@10={rm['Hits@10']:.4f}")
        log("dry-run OK: the shared harness ranks this arm's regimes. No model trained.")
        return

    # ----- real runs ----- #
    import torch
    from pykeen.triples import TriplesFactory
    try:
        from keep_awake import keep_awake; keep_awake()
    except Exception:
        pass
    device = "cuda" if torch.cuda.is_available() else "cpu"
    log(f"device = {device}")
    if device != "cuda":
        log("WARNING no CUDA: full-graph message passing on CPU is not tractable here.")

    cfgs = rgcn_config(num_layers=args.num_layers, edge_dropout=args.edge_dropout,
                       self_loop_dropout=args.self_loop_dropout,
                       decomposition=args.decomposition, num_bases=args.num_bases,
                       edge_weighting=args.edge_weighting, num_negs=args.num_negs,
                       loss=args.loss, lr=args.lr)
    kge_cfgs = rk.model_config()
    for m in args.reference_models:
        if m not in kge_cfgs:
            log(f"WARNING unknown --reference-models entry {m!r}; known: {sorted(kge_cfgs)}")
            continue
        cfgs[m] = dict(kge_cfgs[m])
        if args.lr is not None:
            cfgs[m]["lr"] = args.lr
    models = ["RGCN"] + [m for m in args.reference_models if m in kge_cfgs]
    log(f"models: {models} (RGCN + matched-subgraph reference control(s))")

    poolcache: dict = {}
    plan = [(rd, seed, model) for rd in regimes for seed in args.seeds for model in models]
    durations: list = []
    t_grid = time.time()
    for idx, (rd, seed, model_name) in enumerate(plan, 1):
        if rd["short"] not in poolcache:
            poolcache[rd["short"]] = (
                le._category_pools(rd["train_edges"], rd["hub_set"]),
                le._known_tails_by_head(rd["train_edges"], rd["test_edges"],
                                        rd["test_edges"][:, 0]),
            )
        pools, known = poolcache[rd["short"]]
        cfg = dict(cfgs[model_name])
        gnn_meta = cfg.pop("_gnn_meta", None)
        outp = rk.main_run_path(results_dir, model_name, rd["short"], seed)
        eta = (float(np.mean(durations)) * (len(plan) - idx + 1) / 60.0) if durations else float("nan")
        log(f"[grid {idx}/{len(plan)}] {model_name} {rd['short']} seed{seed} "
            f"| elapsed {(time.time()-t_grid)/60:.1f}min | "
            + (f"eta ~{eta:.0f}min" if durations else "eta (measuring)"))

        t0 = time.time()
        trained = rk.train_eval_write(cfg, model_name, rd, args.dim, args.epochs, seed,
                                      device, args, outp, pools, known, TriplesFactory, torch)
        if trained:
            durations.append(time.time() - t0)
        # Stamp provenance the shared writer does not know about. Done by patching the
        # JSON rather than by editing run_kge's writer, so the frozen KGE path is untouched.
        if outp.exists():
            try:
                r = json.loads(outp.read_text())
                r["method_class"] = (gnn_meta or {}).get("method_class", "kge")
                r["parity"] = parity
                r["subgraph_seed"] = args.subgraph_seed
                if gnn_meta:
                    r["gnn"] = gnn_meta
                elif model_name in kge_cfgs:
                    r["role"] = ("matched-subgraph reference control for RGCN "
                                 "(same graph, same negatives, same harness)")
                outp.write_text(json.dumps(r, indent=2))
            except Exception as e:
                log(f"    WARNING could not stamp provenance into {outp.name}: "
                    f"{type(e).__name__}: {e}")

    payload = rk.merge_results(results_dir, results_dir / "gnn_summary.json")

    # Reference R0 MRR (mean over seeds) for gate check 4.
    ref_r0 = None
    for m in args.reference_models:
        rs = _load_runs(results_dir, m).get("R0", [])
        if rs:
            ref_r0 = float(np.mean([r["MRR"] for r in rs]))
            log(f"matched-subgraph reference {m}: R0 MRR = {ref_r0:.4f} "
                f"(n_seeds={len(rs)}) -- the fair yardstick for the gate")
            break

    gate = convergence_gate(results_dir, args, ref_r0_mrr=ref_r0)
    report_regime_delta(results_dir, gate, models)

    if args.stratify:
        n_test = int(np.median([len(rd["test_edges"]) for rd in regimes]))
        if args.limit_test:
            log("skip degree-stratified reanalysis: --limit-test truncated the test set, so "
                "the per-edge ranks do not align to test.csv row order.")
        elif gate["verdict"] != "pass":
            log("skip degree-stratified reanalysis: the convergence gate FAILED, so "
                "stratified numbers from this model would be uninterpretable.")
        else:
            regs = [rd["short"] for rd in regimes]
            run_stratify(results_dir, models, regs, args.seeds, parity)
            log(f"  (stratified over {n_test} test edges)")

    log(f"done. {len(payload.get('runs', []))} run(s) under {results_dir}")


if __name__ == "__main__":
    main()
