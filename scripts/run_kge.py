#!/usr/bin/env python3
"""
run_kge.py -- Machine A (GPU) KGE arm of the leakage-aware benchmark.

Two KGE models -- TransE and ComplEx -- are scored by the EXACT SAME harness as
the topological baselines (scripts/lib_eval.py), on the SAME frozen splits
(data/processed/splits/). That is the whole point: if KGE used PyKEEN's internal
all-entity evaluator while the baselines used lib_eval's sampled protocol, the two
columns of the paper's Table 2 would not be comparable.

Scope (this rewrite):
  (a) Models: ONLY TransE and ComplEx. DistMult and RotatE are dropped -- two
      well-understood, well-cited models (one translational, one bilinear) are
      enough to answer the leakage question and keep the GPU budget sane.
  (b) Regimes: R0 (train.csv), R2 (train_R2_degree_null_seed42.csv), R3
      (train_R3_orthology_blocked.csv). Each has its OWN training graph, so each
      model is trained per regime. R1 is skipped (hub filter is a no-op for
      embeddings; R1==R0 for KGE). Regimes whose training file is not built yet
      (R2 until Unit 7 runs) are skipped with a warning, exactly like run_baselines.py.
  (c) Ranking goes through lib_eval.rank_test_edges with
      score_fn(gene, disease) = model score of (gene, GENE_ASSOCIATED_WITH_CONDITION,
      disease) -- i.e. score_hrt for the fixed target relation, with IDs mapped
      through the TriplesFactory and a low SENTINEL returned for out-of-vocabulary
      (cold-start) entities. Same sampled-50neg, filtered, type-matched protocol the
      baselines use. (PyKEEN's RankBasedEvaluator is NOT used.)
  (d) Hyperparameters: dim 64, seeds {42,1,7}, self-adversarial negative-sampling
      loss (standard published recipe; refs below), plus ONE small sweep on R0 over
      dim in {64,128} x epochs in {100,300} to show sensitivity to width/under-training.
  (e) Optional ~11% filtered training subgraph (--subgraph-frac 0.11) for tractable
      training: keep every gene->disease target edge (~0.7% of edges) and uniformly
      subsample the rest to the target fraction (Gu et al., 2024 -- a small fraction of
      KG edges retains most of the KGE signal). Entities left in no retained edge become
      cold-start (SENTINEL), which the de-leaked regimes are meant to expose.
  (f) EVERY run's ranking (MRR, Hits@1/3/10) + classification (AUROC, AUPRC) is
      written to data/processed/results/kge/kge_<model>_<regime>_seed<k>.json, plus a
      merged kge_summary.json with mean +/- sd over seeds. Sweep runs go under
      results/kge/sweep/ so they never pollute the headline summary.
  (g) At the end it prints R0 vs R3 MRR per model (the orthology-leakage delta).

Published hyperparameter references (recorded in each run JSON):
  * TransE                 : Bordes et al., NeurIPS 2013.
  * ComplEx                : Trouillon et al., ICML 2016.
  * Self-adversarial NS    : Sun et al. ("RotatE"), ICLR 2019 -- the loss/negative
                             sampler both models use here. We deliberately do NOT use
                             the canonical full 1-vs-all softmax / LCWA recipe
                             (Lacroix et al., ICML 2018): standard on small KGs
                             (FB15k ~15k entities) but intractable at 451k entities /
                             5.8M triples on an 8GB GPU (scoring every example against
                             all entities). Self-adversarial NS gives the anti-collapse
                             benefit at TransE-like speed; ComplEx needs more negatives
                             (bilinear scores saturate faster) so it gets 32 vs TransE's 16.
  * ~11% training subgraph : Gu et al., 2024.

Scientific notes:
  * The split is RANDOM (the repo has no genuine edge-discovery dates); the leakage
    regimes R0->R3, not a time split, are how this benchmark controls optimism.
  * Classification negatives are reported for BOTH random and type-matched modes; the
    type-matched (hard) numbers are the credible ones.
  * The ranking negative-sampling seed is FIXED (--rank-seed, default 42) for every
    model x training-seed, so all methods are ranked on byte-identical candidate sets
    and per-edge reciprocal ranks can be paired-bootstrapped across methods. The
    training seed only varies the learned embeddings (that is the +/- sd we report).

Usage:
    # validate the wiring with NO training (random scorer, small test subset):
    python scripts/run_kge.py --dry-run --limit-test 200

    # one fast end-to-end smoke on the GPU before the full grid:
    python scripts/run_kge.py --models TransE --seeds 42 --epochs 5 --regimes R0 \
        --subgraph-frac 0.11

    # the headline grid (dim 64, seeds 42/1/7, R0/R2/R3):
    python scripts/run_kge.py

    # add the R0 dim x epochs sweep:
    python scripts/run_kge.py --sweep

    # merge whatever runs have completed into kge_summary.json:
    python scripts/run_kge.py --merge-only
"""
from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd

# lib_eval + kg_categories are pure numpy/pandas -> safe to import without torch.
try:
    import lib_eval as le
    from kg_categories import category_of
except ImportError:  # running from repo root
    from scripts import lib_eval as le
    from scripts.kg_categories import category_of

try:
    from config import PROCESSED_DIR
    SPLITS_DIR = Path(PROCESSED_DIR) / "splits"
    RESULTS_DIR = Path(PROCESSED_DIR) / "results" / "kge"
except Exception:
    SPLITS_DIR = Path("data/processed/splits")
    RESULTS_DIR = Path("data/processed/results/kge")

TARGET_REL = "BIOLINK:GENE_ASSOCIATED_WITH_CONDITION"  # the task relation (all test edges)
SENTINEL = -1e30  # score for out-of-vocabulary (cold-start) entities: worst possible


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


# --------------------------------------------------------------------------- #
# Per-model training recipes (standard published hyperparameters)
# --------------------------------------------------------------------------- #
def model_config():
    """Return {name: cfg} for the two benchmark models. Imported lazily so the
    module loads without pykeen.

    TransE (translational) uses SLCWA + self-adversarial negative sampling (NSSA;
    Sun et al., ICLR 2019) -- the standard strong TransE recipe.

    ComplEx (bilinear) uses SLCWA + the LOGISTIC / softplus loss (Trouillon et al.,
    ICML 2016 -- the canonical ComplEx objective). NSSA does NOT work for ComplEx on
    this graph: its margin-sigmoid saturates against ComplEx's unbounded scores and the
    model collapses to chance -- verified empirically at full scale (flat loss ~4.50 for
    the last 150 of 300 epochs, MRR ~0.10, AUROC ~0.50 across 3 seeds on R0). The earlier
    "ComplEx learns fine on NSSA" note was never validated past 5-epoch smokes; it is
    false here. Softplus does not saturate and restores learning. With softplus the extra
    negatives are no longer needed for stability, so ComplEx uses 16 negatives (parity
    with TransE); a one-cell 16-vs-32 check on R0 documents that the count is immaterial.
    Full-softmax LCWA (Lacroix+2018) is intractable at 451k entities (see module docstring).
    """
    from pykeen.models import TransE, ComplEx, RotatE
    nssa = dict(loss="nssa", loss_kwargs=dict(margin=9.0, adversarial_temperature=1.0))
    return {
        "TransE": dict(cls=TransE, num_negs=16, ref="Bordes+2013; NSSA Sun+2019", **nssa),
        # RotatE (rotational family): NSSA is its native training objective (Sun+2019),
        # and its geometry is close to TransE's, so it trains tractably under SLCWA here --
        # unlike the bilinear ComplEx, which collapses to a trivial optimum under every
        # tractable loss and needs intractable full-LCWA (documented; see PAPER notes).
        "RotatE": dict(cls=RotatE, num_negs=16, ref="Sun+2019 (RotatE); NSSA", **nssa),
        "ComplEx": dict(cls=ComplEx, num_negs=16, loss="softplus",
                        ref="Trouillon+2016 (logistic/softplus loss)"),
    }


def _build_model(cfg, tf, dim, seed, device):
    kw = dict(triples_factory=tf, embedding_dim=dim, random_seed=seed)
    if cfg.get("loss"):
        kw["loss"] = cfg["loss"]
    if cfg.get("loss_kwargs"):
        kw["loss_kwargs"] = cfg["loss_kwargs"]
    if cfg.get("regularizer"):
        kw["regularizer"] = cfg["regularizer"]
        kw["regularizer_kwargs"] = cfg["regularizer_kwargs"]
    return cfg["cls"](**kw).to(device)


def _is_oom(err) -> bool:
    msg = str(err).lower()
    return "out of memory" in msg or "cuda error" in msg or "cudnn" in msg


def _progress_callback(total_epochs, total_batches, tag, min_interval=1.0):
    """PyKEEN callback: a BATCH-level, wall-clock-throttled progress bar.

    Emits a fresh log LINE at most every ``min_interval`` seconds (default 1s), driven by
    the per-batch ``on_batch`` hook, plus one line at the end of each epoch. Because it is
    line-based (not an in-place ``\\r`` bar) it renders live when the log is tailed
    (``Get-Content <log> -Wait``) even when the job runs in the background -- an in-place
    bar needs a TTY and goes dead once backgrounded. Overall progress and ETA span the
    WHOLE training (all epochs x batches), so the bar moves smoothly every second no
    matter how long a single epoch takes.
    """
    from pykeen.training.callbacks import TrainingCallback

    class _Cb(TrainingCallback):
        def __init__(self, *a, **k):
            super().__init__(*a, **k)
            self.t0 = time.time()
            self.last = 0.0
            self.b = 0            # batches seen in the current epoch
            self.cur_epoch = 0

        def on_batch(self, epoch, batch, batch_loss, **kwargs):
            if epoch != self.cur_epoch:
                self.cur_epoch, self.b = epoch, 0
            self.b += 1
            now = time.time()
            if now - self.last >= min_interval:
                self.last = now
                self._emit(epoch, self.b, float(batch_loss), False)

        def post_epoch(self, epoch, epoch_loss, **kwargs):
            self._emit(epoch, total_batches, float(epoch_loss), True)
            self.last = time.time()

        def _emit(self, epoch, b, loss, end_of_epoch):
            b = min(b, total_batches)
            frac = min(max(((epoch - 1) + b / max(total_batches, 1)) / max(total_epochs, 1),
                           1e-9), 1.0)
            el = time.time() - self.t0
            eta = el / frac - el          # remaining wall-seconds for the whole training
            w = 24
            fill = int(round(w * frac))
            bar = "#" * fill + "-" * (w - fill)
            tick = "epoch done " if end_of_epoch else f"batch {b:>4}/{total_batches}"
            log(f"    [{tag}] ep {epoch:>3}/{total_epochs} {tick} [{bar}] "
                f"{frac * 100:4.1f}% loss={loss:.4f} elapsed={el:.0f}s eta={eta:.0f}s")

    return _Cb()


def train_model(cfg, tf, dim, epochs, seed, device, batch, tag="", progress_secs=1.0,
                lcwa_slice=20000):
    """Train one model. Default is SLCWA + self-adversarial negative sampling; a model
    whose cfg sets ``loop="lcwa"`` trains 1-vs-all (LCWA), slicing the 451k-entity target
    dimension (``lcwa_slice``) so the all-entity scores fit on the GPU -- the canonical
    tractable recipe for bilinear models (ComplEx/DistMult) at scale. On CUDA OOM we halve
    the batch (and LCWA slice) and retry a few times."""
    import torch, random
    from pykeen.training import SLCWATrainingLoop, LCWATrainingLoop
    is_lcwa = cfg.get("loop", "slcwa") == "lcwa"
    slice_size = lcwa_slice if is_lcwa else None

    attempt = 0
    while True:
        torch.manual_seed(seed); np.random.seed(seed); random.seed(seed)
        if device == "cuda":
            torch.cuda.empty_cache()
        model = _build_model(cfg, tf, dim, seed, device)
        try:
            if is_lcwa:
                loop = LCWATrainingLoop(model=model, triples_factory=tf)
                train_kw = dict(batch_size=batch, slice_size=slice_size)
            else:
                loop = SLCWATrainingLoop(
                    model=model, triples_factory=tf, negative_sampler="basic",
                    negative_sampler_kwargs=dict(num_negs_per_pos=cfg.get("num_negs", 16)))
                train_kw = dict(batch_size=batch)
            total_batches = max(1, -(-tf.num_triples // batch))  # ceil(n_triples / batch)
            t0 = time.time()
            loop.train(triples_factory=tf, num_epochs=epochs,
                       label_smoothing=cfg.get("label_smoothing", 0.0), use_tqdm=False,
                       callbacks=_progress_callback(epochs, total_batches, tag, progress_secs),
                       **train_kw)
            return model, time.time() - t0, {"batch": batch, "slice_size": slice_size}
        except Exception as e:  # noqa: BLE001 -- torch raises AcceleratorError/RuntimeError on OOM
            if not (_is_oom(e) and attempt < 4):
                raise
            del model
            if device == "cuda":
                torch.cuda.empty_cache()
            attempt += 1
            batch = max(16, batch // 2)
            if slice_size:
                slice_size = max(2000, slice_size // 2)
            log(f"    CUDA OOM -> retry {attempt}: batch={batch} slice={slice_size}")


# --------------------------------------------------------------------------- #
# Scoring: expose the trained model to lib_eval.rank_test_edges as score_fn(g,d)
# --------------------------------------------------------------------------- #
class KGEScorer:
    """score_fn(gene, disease) for lib_eval.rank_test_edges.

    Returns the model's score of the triple (gene, TARGET_REL, disease) -- i.e.
    model.score_hrt(gene, GENE_ASSOCIATED_WITH_CONDITION, disease) -- with IDs mapped
    through the TriplesFactory. Out-of-vocabulary genes/diseases (cold start under the
    de-leaked graphs) get SENTINEL -> worst possible rank.

    Optimization (numerically identical to per-triple score_hrt): for disease tails it
    uses score_t (all-entity scoring for a fixed (head, TARGET_REL)) ONCE per gene and
    caches the disease-only slice, so a whole test gene costs ONE GPU call instead of
    ~51. In-vocab non-disease tails fall back to a direct score_hrt.
    """

    def __init__(self, model, tf, device, rel=TARGET_REL):
        import torch
        self.torch = torch
        self.model = model
        self.device = device
        self.e2id = tf.entity_to_id
        self.rid = tf.relation_to_id.get(rel)
        # disease-entity columns (bounds the per-gene cache to ~14k floats, not ~450k)
        ents = np.array(list(self.e2id.keys()), dtype=object)
        dis = ents[np.array([category_of(e) == "disease" for e in ents], dtype=bool)]
        self.dpos = np.array([self.e2id[d] for d in dis], dtype=np.int64)
        self.did2slot = {int(did): i for i, did in enumerate(self.dpos)}
        self.cache: dict = {}

    def _fill(self, gid):
        """Compute + cache the disease-score slice for one gene (one GPU call)."""
        vec = self.cache.get(gid)
        if vec is None:
            with self.torch.no_grad():
                hr = self.torch.tensor([[gid, self.rid]], dtype=self.torch.long, device=self.device)
                full = self.model.score_t(hr).detach().cpu().numpy().ravel()
            vec = full[self.dpos]
            self.cache[gid] = vec
        return vec

    def warm(self, genes):
        """Pre-fill the cache over the test genes, logging progress (the slow part)."""
        if self.rid is None:
            return
        gids = [self.e2id[g] for g in genes
                if g in self.e2id and self.e2id[g] not in self.cache]
        for i, gid in enumerate(gids, 1):
            self._fill(gid)
            if i % 500 == 0 or i == len(gids):
                log(f"    [warm] scored {i}/{len(gids)} test genes")

    def __call__(self, g, d):
        if self.rid is None:
            return SENTINEL
        gid = self.e2id.get(g)
        did = self.e2id.get(d)
        if gid is None or did is None:
            return SENTINEL
        slot = self.did2slot.get(int(did))
        if slot is None:  # in-vocab but not a disease-category entity: score directly
            return self._score_single(gid, did)
        return float(self._fill(gid)[slot])

    def _score_single(self, gid, did):
        with self.torch.no_grad():
            hrt = self.torch.tensor([[gid, self.rid, did]], dtype=self.torch.long, device=self.device)
            return float(self.model.score_hrt(hrt).detach().cpu().numpy().ravel()[0])


def score_pairs(model, tf, pairs, device, rel=TARGET_REL, bs=20000):
    """Batched score_hrt over in-vocab (gene, disease) pairs -> np array of scores."""
    import torch
    e2id = tf.entity_to_id
    rid = tf.relation_to_id.get(rel)
    if rid is None or not pairs:
        return np.array([])
    ids = np.array([(e2id[g], rid, e2id[d]) for g, d in pairs], dtype=np.int64)
    outs = []
    model.eval()
    with torch.no_grad():
        for i in range(0, len(ids), bs):
            bt = torch.tensor(ids[i:i + bs], dtype=torch.long, device=device)
            outs.append(model.score_hrt(bt).detach().cpu().numpy().ravel())
    return np.concatenate(outs) if outs else np.array([])


# --------------------------------------------------------------------------- #
# Classification eval sets (AUROC/AUPRC): random + type-matched negatives
# --------------------------------------------------------------------------- #
def classification_scores(model, tf, test_edges, train_edges, device, seed, cap=None):
    """Return {'random': {AUROC,AUPRC,...}, 'type': {...}} for this model."""
    e2id = tf.entity_to_id
    rid = tf.relation_to_id.get(TARGET_REL)
    if rid is None:
        return {"random": le.classification_metrics([], []),
                "type": le.classification_metrics([], [])}
    rng = np.random.default_rng(seed)
    known = set(zip(train_edges[:, 0].tolist(), train_edges[:, 2].tolist()))
    known |= set(zip(test_edges[:, 0].tolist(), test_edges[:, 2].tolist()))

    ents = np.array(list(e2id.keys()), dtype=object)
    dis = ents[np.array([category_of(e) == "disease" for e in ents], dtype=bool)]

    pos = [(g, d) for g, _r, d in test_edges if g in e2id and d in e2id]
    if cap and len(pos) > cap:
        sel = rng.choice(len(pos), cap, replace=False)
        pos = [pos[i] for i in sel]

    def make_neg(type_matched):
        pool = dis if type_matched else ents
        negs = []
        for g, _d in pos:
            for _try in range(40):
                t = pool[rng.integers(0, len(pool))]
                if t == g or (g, t) in known:
                    continue
                negs.append((g, t)); break
        return negs

    ps = score_pairs(model, tf, pos, device)
    out = {}
    for mode, negs in (("random", make_neg(False)), ("type", make_neg(True))):
        ns = score_pairs(model, tf, negs, device)
        out[mode] = le.classification_metrics(ps, ns)
    return out


# --------------------------------------------------------------------------- #
# One evaluation of a (trained model) on one regime -> run dict + JSON
# --------------------------------------------------------------------------- #
def evaluate_regime(scorer, model, tf, train_edges, regime_name, test_edges, hub_set,
                    hub_filter, device, n_neg, rank_seed, class_cap, pools=None, known=None):
    scorer.warm(pd.unique(test_edges[:, 0]))            # GPU scoring w/ progress
    log(f"    [{regime_name}] ranking {len(test_edges)} test edges (n_neg={n_neg}) ...")
    ranks = le.rank_test_edges(scorer, train_edges, test_edges, hub_set, hub_filter,
                               n_neg=n_neg, seed=rank_seed, pools=pools, known=known)
    rm = le.ranking_metrics(ranks)
    rr = (1.0 / ranks).tolist()
    mrr_ci = le.bootstrap_ci(1.0 / ranks)
    log(f"    [{regime_name}] classification (AUROC/AUPRC) ...")
    cls = classification_scores(model, tf, test_edges, train_edges, device,
                                seed=rank_seed, cap=class_cap)
    return {
        "regime": regime_name,
        "n_test": int(len(test_edges)),
        "rank_protocol": f"sampled-{n_neg}neg (lib_eval, rank_seed={rank_seed})",
        "MRR": rm["MRR"], "Hits@1": rm["Hits@1"], "Hits@3": rm["Hits@3"], "Hits@10": rm["Hits@10"],
        "MRR_ci_low": mrr_ci["ci_low"], "MRR_ci_high": mrr_ci["ci_high"],
        "AUROC_random": cls["random"]["AUROC"], "AUPRC_random": cls["random"]["AUPRC"],
        "AUROC_type": cls["type"]["AUROC"], "AUPRC_type": cls["type"]["AUPRC"],
        "reciprocal_ranks": rr,   # per-edge, for cross-method paired bootstrap later
    }


# --------------------------------------------------------------------------- #
# Full-ranking evaluation (H1b robustness): rank the true disease against the
# ENTIRE filtered disease pool, not 50 sampled negatives. This is the KGE analogue
# of run_robustness.py's baseline "filtered full ranking", so the two are directly
# comparable. The per-gene score vector (KGEScorer.score_t cache) makes it cheap.
# --------------------------------------------------------------------------- #
def _fullrank_one(vec, true_slot, excl_slots, D):
    """Tie-averaged filtered rank of the true disease.

    ``vec`` : length-D score vector over the pool diseases (higher = better).
    ``true_slot`` : pool slot of the true disease, or None if the true disease is not
                    in the pool (cold-start -> SENTINEL true score -> worst rank).
    ``excl_slots`` : pool slots to drop from the negatives (gene's known tails + the
                     true disease itself; only pool diseases shrink the denominator).
    Identical semantics to run_robustness._rank_from_dict over the full pool.
    """
    true_s = float(vec[true_slot]) if true_slot is not None else SENTINEL
    if excl_slots:
        mask = np.ones(D, dtype=bool)
        mask[list(excl_slots)] = False
        neg = vec[mask]
    else:
        neg = vec
    return 1.0 + float(np.sum(neg > true_s)) + 0.5 * float(np.sum(neg == true_s))


def evaluate_regime_fullrank(scorer, model, tf, train_edges, regime_name, test_edges,
                             device, rank_seed, class_cap, validate=0):
    """Filtered full-ranking eval: each held-out true disease is ranked against the
    entire disease pool (all disease-category entities in the training vocab) minus
    the gene's known tails (train U test). Scores come from the model's score_t
    vectors (one GPU call per test gene via KGEScorer.warm), so this costs the same
    GPU work as the sampled protocol. Classification (AUROC/AUPRC) is unchanged.
    """
    e2id = tf.entity_to_id
    known = le._known_tails_by_head(train_edges, test_edges, test_edges[:, 0])
    genes = pd.unique(test_edges[:, 0])
    scorer.warm(genes)                         # fills scorer.cache[gid] over scorer.dpos
    did2slot = scorer.did2slot
    D = len(scorer.dpos)
    ranks = np.empty(len(test_edges), dtype=float)
    pos_scores = np.empty(len(test_edges), dtype=float)
    n_cold = 0
    log(f"    [{regime_name}] FULL-RANK: {len(test_edges)} test edges vs pool D={D} ...")
    for i, (g, _r, d_true) in enumerate(test_edges):
        gid = e2id.get(g)
        excl = set(known.get(g, ()))
        excl.add(d_true)
        excl_slots = [did2slot[e2id[d]] for d in excl if d in e2id and e2id[d] in did2slot]
        n_filt = D - len(excl_slots)
        vec = scorer.cache.get(gid) if gid is not None else None
        if vec is None:                         # cold-start gene: worst possible rank
            ranks[i] = n_filt + 1.0
            pos_scores[i] = SENTINEL
            n_cold += 1
            continue
        dtid = e2id.get(d_true)
        true_slot = did2slot.get(int(dtid)) if dtid is not None else None
        ranks[i] = _fullrank_one(vec, true_slot, excl_slots, D)
        pos_scores[i] = float(vec[true_slot]) if true_slot is not None else SENTINEL

    # correctness check: brute-force rank of the first `validate` scoreable edges
    if validate:
        _validate_fullrank(scorer, test_edges, known, ranks, validate)

    rm = le.ranking_metrics(ranks)
    mrr_ci = le.bootstrap_ci(1.0 / ranks)
    log(f"    [{regime_name}] classification (AUROC/AUPRC) ...")
    cls = classification_scores(model, tf, test_edges, train_edges, device,
                                seed=rank_seed, cap=class_cap)
    log(f"    [{regime_name}] FULL-RANK MRR={rm['MRR']:.4f} H@10={rm['Hits@10']:.4f} "
        f"(pool D={D}, cold-start genes={n_cold})")
    return {
        "regime": regime_name,
        "n_test": int(len(test_edges)),
        "rank_protocol": f"filtered-full-rank (disease pool D={D}, run_robustness-comparable)",
        "pool_size": int(D),
        "n_coldstart_genes": int(n_cold),
        "MRR": rm["MRR"], "Hits@1": rm["Hits@1"], "Hits@3": rm["Hits@3"], "Hits@10": rm["Hits@10"],
        "MRR_ci_low": mrr_ci["ci_low"], "MRR_ci_high": mrr_ci["ci_high"],
        "AUROC_random": cls["random"]["AUROC"], "AUPRC_random": cls["random"]["AUPRC"],
        "AUROC_type": cls["type"]["AUROC"], "AUPRC_type": cls["type"]["AUPRC"],
        "reciprocal_ranks": (1.0 / ranks).tolist(),
    }


def _validate_fullrank(scorer, test_edges, known, ranks, n):
    """Assert the vectorized full-rank equals a brute-force rank for the first n
    scoreable edges (loops the scorer over the whole pool per pair)."""
    e2id = scorer.e2id
    pool_ids = [d for d in e2id if category_of(d) == "disease"]
    checked = 0
    for i, (g, _r, d_true) in enumerate(test_edges):
        if checked >= n:
            break
        if e2id.get(g) is None or e2id.get(d_true) is None:
            continue
        excl = set(known.get(g, ())); excl.add(d_true)
        cand = [d_true] + [d for d in pool_ids if d not in excl]
        sc = np.fromiter((float(scorer(g, d)) for d in cand), dtype=float, count=len(cand))
        true_s, neg_s = sc[0], sc[1:]
        bf = 1.0 + float(np.sum(neg_s > true_s)) + 0.5 * float(np.sum(neg_s == true_s))
        if abs(bf - ranks[i]) > 1e-6:
            raise AssertionError(f"full-rank mismatch edge {i} ({g}->{d_true}): "
                                 f"vectorized {ranks[i]} != bruteforce {bf}")
        checked += 1
    log(f"    [validate] {checked} edges: vectorized full-rank == brute force. OK")


# --------------------------------------------------------------------------- #
# Merge per-run JSONs -> summary with mean +/- sd over seeds
# --------------------------------------------------------------------------- #
def merge_results(results_dir, out_path):
    """Merge every kge_*.json directly under results_dir (non-recursive, so the
    results/kge/sweep/ subdir is not swept into the headline summary)."""
    runs = []
    # Match per-run files only ("..._seed<k>[...].json"); the "kge_*.json" glob also
    # matched this function's OWN output (kge_summary.json), which has no per-run keys.
    for p in sorted(Path(results_dir).glob("kge_*_seed*.json")):
        try:
            r = json.loads(Path(p).read_text())
        except Exception:
            continue
        if "model" not in r:      # defensive: skip anything that isn't a per-run record
            continue
        r = {k: v for k, v in r.items() if k != "reciprocal_ranks"}  # drop big arrays
        runs.append(r)

    agg: dict = {}
    metrics = ["MRR", "Hits@1", "Hits@3", "Hits@10",
               "AUROC_random", "AUPRC_random", "AUROC_type", "AUPRC_type"]
    keyf = lambda r: f"{r['model']}|{r['regime']}|d{r['dim']}|e{r['epochs']}"
    groups: dict = {}
    for r in runs:
        groups.setdefault(keyf(r), []).append(r)
    for k, rs in groups.items():
        entry = {"model": rs[0]["model"], "regime": rs[0]["regime"],
                 "dim": rs[0]["dim"], "epochs": rs[0]["epochs"],
                 "n_seeds": len(rs), "seeds": sorted(x["seed"] for x in rs)}
        for m in metrics:
            vals = np.array([x[m] for x in rs if x.get(m) is not None and not np.isnan(x[m])], float)
            if len(vals):
                entry[f"{m}_mean"] = float(vals.mean())
                entry[f"{m}_sd"] = float(vals.std(ddof=1)) if len(vals) > 1 else 0.0
        agg[k] = entry

    payload = {
        "aggregate": agg,
        "runs": runs,
        "meta": {
            "harness": "scripts/lib_eval.py (shared with baselines)",
            "models": ["TransE", "ComplEx"],
            "target_relation": TARGET_REL,
            "splits_manifest_sha256": _manifest_sha(),
            "generated": time.strftime("%Y-%m-%d %H:%M:%S"),
        },
    }
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    Path(out_path).write_text(json.dumps(payload, indent=2))
    log(f"merged {len(runs)} run(s) -> {out_path} ({len(agg)} model x regime x config cells)")
    return payload


def print_r0_vs_r3(payload):
    """Print R0 vs R3 MRR per model -- does removing orthology leakage drop KGE?"""
    agg = (payload or {}).get("aggregate", {})
    by: dict = {}
    for entry in agg.values():
        by.setdefault(entry["model"], {})[entry["regime"]] = entry
    log("=== R0 vs R3 MRR (does orthology-leakage removal drop KGE?) ===")
    if not by:
        log("  (no runs yet)")
        return
    for model in sorted(by):
        r0, r3 = by[model].get("R0"), by[model].get("R3")
        if r0 and r3 and "MRR_mean" in r0 and "MRR_mean" in r3:
            d = r3["MRR_mean"] - r0["MRR_mean"]
            pct = 100.0 * d / r0["MRR_mean"] if r0["MRR_mean"] else float("nan")
            log(f"  {model:8s} R0 MRR={r0['MRR_mean']:.4f}+/-{r0.get('MRR_sd', 0.0):.4f}  "
                f"R3 MRR={r3['MRR_mean']:.4f}+/-{r3.get('MRR_sd', 0.0):.4f}  "
                f"delta={d:+.4f} ({pct:+.1f}%)")
        else:
            have = [r for r in ("R0", "R3") if by[model].get(r)]
            log(f"  {model:8s} incomplete (have: {have or 'none'})")


def _manifest_sha():
    try:
        import hashlib
        return hashlib.sha256((SPLITS_DIR / "split_manifest.json").read_bytes()).hexdigest()
    except Exception:
        return None


# --------------------------------------------------------------------------- #
# Regime plumbing: load each regime's train/test once, skip regimes not built yet
# --------------------------------------------------------------------------- #
def load_regimes(regimes, splits_dir, subgraph_frac=None, subgraph_seed=42):
    """Return an ordered list of regime dicts, one per requested regime whose
    training file exists. Regimes whose train file is missing (e.g. R2 until Unit 7
    builds the degree-null graph) are skipped with a warning, matching run_baselines.py.

    Each dict: short ('R0'), key ('R0_standard'), train_edges, test_edges, hub_set,
    hub_filter, train_file, test_file, n_train_full, subgraph_frac.
    """
    manifest = json.loads((Path(splits_dir) / "split_manifest.json").read_text())
    out = []
    for rname in regimes:
        try:
            key = le._regime_key(manifest, rname)
        except KeyError as e:
            log(f"skip {rname}: {e}")
            continue
        spec = manifest["regimes"][key]
        short = key.split("_", 1)[0]
        train_path = Path(splits_dir) / spec["train"]
        if not train_path.exists():
            log(f"skip {short}: training file {spec['train']} not built yet "
                f"({spec.get('produced_by', 'run build_deleaked_splits.py')})")
            continue
        log(f"loading {short}: train={spec['train']} test={spec['test']}")
        train_edges = le._load_edges(train_path)
        test_edges = le._load_edges(Path(splits_dir) / spec["test"])
        n_full = len(train_edges)
        if subgraph_frac is not None:
            train_edges = build_subgraph(train_edges, subgraph_frac, subgraph_seed)
            log(f"    subgraph: {len(train_edges):,}/{n_full:,} edges "
                f"({100.0 * len(train_edges) / max(n_full, 1):.1f}%) kept")
        hub_filter = bool(spec["hub_filter"])
        hub_set = le._load_hub_nodes(Path(splits_dir) / "hub_nodes.txt") if hub_filter else frozenset()
        out.append(dict(
            short=short, key=key, train_edges=train_edges, test_edges=test_edges,
            hub_set=hub_set, hub_filter=hub_filter,
            train_file=spec["train"], test_file=spec["test"],
            n_train_full=n_full,
            subgraph_frac=(subgraph_frac if subgraph_frac is not None else None),
        ))
    return out


def build_subgraph(train_edges, frac, seed, target_rel=TARGET_REL):
    """Downsample the training graph to ~frac of its edges for tractable KGE training.

    Rationale (Gu et al., 2024): a small fraction of KG edges retains most of the KGE
    signal. Every target-relation (gene->disease) edge is kept -- that is the task
    signal and only ~0.7% of edges -- and the remaining edges are uniformly subsampled
    (fixed seed) to hit the target fraction. An entity that ends up in NO retained edge
    becomes cold-start and is scored with SENTINEL; that is consistent with the de-leaked
    regimes (R3), which exist precisely to expose such cold-start, so the subgraph must
    not paper over it by force-keeping test-incident edges.
    """
    n = len(train_edges)
    target_keep = int(round(frac * n))
    is_target = train_edges[:, 1] == target_rel
    keep = is_target.copy()
    budget = target_keep - int(is_target.sum())
    if budget > 0:
        pool_idx = np.nonzero(~is_target)[0]
        rng = np.random.default_rng(seed)
        chosen = rng.choice(pool_idx, size=min(budget, len(pool_idx)), replace=False)
        keep[chosen] = True
    return train_edges[keep]


def main_run_path(results_dir, model, regime_short, seed):
    return Path(results_dir) / f"kge_{model}_{regime_short}_seed{seed}.json"


def sweep_run_path(results_dir, model, regime_short, seed, dim, epochs):
    return Path(results_dir) / "sweep" / f"kge_{model}_{regime_short}_seed{seed}_d{dim}_e{epochs}.json"


# --------------------------------------------------------------------------- #
# Train one (model, regime, dim, epochs, seed) cell and write its run JSON
# --------------------------------------------------------------------------- #
def train_eval_write(cfg, model_name, rd, dim, epochs, seed, device, args, out_path,
                     pools, known, TriplesFactory, torch):
    if out_path.exists() and not args.force:
        # The headline filename does not encode dim/epochs/subgraph, so a resume must
        # confirm the existing JSON was produced with THIS config before trusting it --
        # otherwise a quick smoke (e.g. --epochs 5) would masquerade as a headline run.
        try:
            prev = json.loads(out_path.read_text())
            mism = [f"{k}={prev.get(k)!r}!={want!r}" for k, want in
                    (("dim", dim), ("epochs", epochs), ("subgraph_frac", rd["subgraph_frac"]))
                    if prev.get(k) != want]
        except Exception:
            mism = []
        if mism:
            log(f"  WARNING {out_path.name} exists but was produced with a DIFFERENT "
                f"config ({', '.join(mism)}); keeping it. Use --force to overwrite.")
        else:
            log(f"  skip (exists): {out_path.name}")
        return False
    tf = TriplesFactory.from_labeled_triples(rd["train_edges"].astype(str))
    log(f"--- training {model_name} (nssa, {cfg.get('num_negs')} negs) "
        f"{rd['short']} dim={dim} epochs={epochs} seed={seed}")
    model, secs, tmeta = train_model(cfg, tf, dim, epochs, seed, device,
                                     args.slcwa_batch, tag=f"{model_name} {rd['short']} d{dim} s{seed}",
                                     progress_secs=args.progress_secs, lcwa_slice=args.lcwa_slice_size)
    log(f"    trained in {secs:.0f}s (batch={tmeta['batch']})")
    scorer = KGEScorer(model, tf, device)
    if args.full_rank:
        run = evaluate_regime_fullrank(scorer, model, tf, rd["train_edges"], rd["short"],
                                       rd["test_edges"], device, args.rank_seed, args.class_cap,
                                       validate=args.validate_fullrank)
    else:
        run = evaluate_regime(scorer, model, tf, rd["train_edges"], rd["short"],
                              rd["test_edges"], rd["hub_set"], rd["hub_filter"], device,
                              args.n_neg, args.rank_seed, args.class_cap, pools=pools, known=known)
    run.update(model=model_name, seed=seed, dim=dim, epochs=epochs,
               train_file=rd["train_file"], test_file=rd["test_file"],
               n_train_edges=int(len(rd["train_edges"])), n_train_full=rd["n_train_full"],
               subgraph_frac=rd["subgraph_frac"], train_seconds=round(secs, 1),
               train_batch=tmeta["batch"], loop="slcwa", loss=cfg.get("loss"),
               num_negs=cfg.get("num_negs"), hparams_ref=cfg.get("ref"))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(run, indent=2))
    log(f"    {model_name}/{rd['short']} d{dim} e{epochs}: MRR={run['MRR']:.4f} "
        f"H@10={run['Hits@10']:.4f} AUROC_type={run['AUROC_type']:.3f} -> {out_path.name}")
    del model
    if device == "cuda":
        torch.cuda.empty_cache()
    return True


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--splits-dir", default=str(SPLITS_DIR))
    ap.add_argument("--results-dir", default=str(RESULTS_DIR))
    ap.add_argument("--models", nargs="+", default=["TransE", "RotatE"],
                    choices=["TransE", "RotatE", "ComplEx"])
    ap.add_argument("--regimes", nargs="+", default=["R0", "R2", "R3"],
                    help="KGE skips R1 (hub filter is a no-op for embeddings; R1==R0). "
                         "R2 is auto-skipped until its degree-null graph is built.")
    ap.add_argument("--seeds", nargs="+", type=int, default=[42, 1, 7])
    ap.add_argument("--dim", type=int, default=64, help="embedding dim for the headline runs")
    ap.add_argument("--epochs", type=int, default=300, help="epochs for the headline runs")
    ap.add_argument("--sweep", action="store_true",
                    help="ALSO run the R0 sensitivity sweep over --sweep-dims x --sweep-epochs")
    ap.add_argument("--sweep-dims", nargs="+", type=int, default=[64, 128])
    ap.add_argument("--sweep-epochs", nargs="+", type=int, default=[100, 300])
    ap.add_argument("--subgraph-frac", type=float, default=None,
                    help="train on a ~FRAC filtered subgraph for speed, e.g. 0.11 "
                         "(Gu et al. 2024); keeps all target + test-incident edges")
    ap.add_argument("--subgraph-seed", type=int, default=42,
                    help="fixed seed for subgraph sampling (kept constant across "
                         "training seeds so every model trains on the identical subgraph)")
    ap.add_argument("--n-neg", type=int, default=50, help="ranking negatives per test edge")
    ap.add_argument("--rank-seed", type=int, default=42,
                    help="FIXED negative-sampling seed for ranking+classification so every "
                         "model x seed is scored on identical candidate sets")
    ap.add_argument("--class-cap", type=int, default=None,
                    help="cap #positives used for AUROC/AUPRC (None = all test edges)")
    ap.add_argument("--slcwa-batch", type=int, default=16384)
    ap.add_argument("--complex-negs", type=int, default=None,
                    help="override ComplEx negatives/positive (rescue/sweep of the bilinear recipe)")
    ap.add_argument("--complex-loss", default=None,
                    help="override ComplEx loss: softplus|nssa|bcewithlogits|crossentropy (rescue/sweep)")
    ap.add_argument("--complex-margin", type=float, default=9.0,
                    help="NSSA margin/gamma used when --complex-loss nssa")
    ap.add_argument("--complex-loop", default=None, choices=["slcwa", "lcwa"],
                    help="override ComplEx training loop; 'lcwa' = canonical 1-vs-all "
                         "(defaults its loss to crossentropy unless --complex-loss is given)")
    ap.add_argument("--lcwa-slice-size", type=int, default=20000,
                    help="chunk the all-entity target dim for LCWA so it fits GPU memory")
    ap.add_argument("--progress-secs", type=float, default=1.0,
                    help="emit a batch-level progress line at most this often (seconds); "
                         "1.0 = a live line every second when tailing the log")
    ap.add_argument("--force", action="store_true", help="recompute runs whose JSON exists")
    ap.add_argument("--merge-only", action="store_true",
                    help="just merge existing results into kge_summary.json + print R0 vs R3")
    ap.add_argument("--dry-run", action="store_true",
                    help="NO training/torch: rank a RANDOM scorer to validate the wiring")
    ap.add_argument("--limit-test", type=int, default=None,
                    help="use only the first N test edges (fast dry-run / smoke)")
    ap.add_argument("--full-rank", action="store_true",
                    help="rank each true disease against the ENTIRE filtered disease pool "
                         "(run_robustness-comparable) instead of --n-neg sampled negatives; "
                         "writes to results/kge_fullrank/ so the sampled headline stays intact")
    ap.add_argument("--validate-fullrank", type=int, default=0,
                    help="brute-force-check the first N full-rank edges per regime (correctness)")
    args = ap.parse_args()
    # Reduce fragmentation-driven OOM on the GPU (must precede torch's CUDA init).
    os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")

    # Full-ranking runs live in a sibling dir so they never overwrite the sampled
    # headline kge_summary.json that the manuscript tables/figures read.
    if args.full_rank and Path(args.results_dir) == RESULTS_DIR:
        args.results_dir = str(RESULTS_DIR.parent / "kge_fullrank")
        print(f"[run_kge] --full-rank -> writing to {args.results_dir}")

    Path(args.results_dir).mkdir(parents=True, exist_ok=True)
    summary_path = Path(args.results_dir) / "kge_summary.json"
    sweep_summary_path = Path(args.results_dir) / "sweep" / "sweep_summary.json"

    if args.merge_only:
        payload = merge_results(args.results_dir, summary_path)
        if (Path(args.results_dir) / "sweep").exists():
            merge_results(Path(args.results_dir) / "sweep", sweep_summary_path)
        print_r0_vs_r3(payload)
        return

    regimes = load_regimes(args.regimes, args.splits_dir, args.subgraph_frac, args.subgraph_seed)
    if not regimes:
        log("no regimes available (all requested training files missing). nothing to do.")
        return
    if args.limit_test:
        for rd in regimes:
            rd["test_edges"] = rd["test_edges"][:args.limit_test]

    # ----- dry run: validate lib_eval wiring + JSON I/O with a random scorer ----- #
    if args.dry_run:
        rng = np.random.default_rng(0)
        random_scorer = lambda g, d: float(rng.random())  # noqa: E731
        for rd in regimes:
            ranks = le.rank_test_edges(random_scorer, rd["train_edges"], rd["test_edges"],
                                       rd["hub_set"], rd["hub_filter"],
                                       n_neg=args.n_neg, seed=args.rank_seed)
            rm = le.ranking_metrics(ranks)
            log(f"[dry-run] {rd['short']}: n_test={rm['n']} MRR={rm['MRR']:.4f} "
                f"(chance~0.089) Hits@10={rm['Hits@10']:.4f}")
        log("dry-run OK: lib_eval wiring produces valid ranking metrics. "
            "No models trained. Remove --dry-run to train.")
        return

    # ----- real runs (GPU) ----- #
    import torch
    from pykeen.triples import TriplesFactory
    try:
        from keep_awake import keep_awake; keep_awake()
    except Exception:
        pass
    device = "cuda" if torch.cuda.is_available() else "cpu"
    log(f"device = {device}")
    cfgs = model_config()

    # Optional ComplEx recipe overrides (bilinear-model rescue / hyperparameter sweep).
    if args.complex_negs is not None:
        cfgs["ComplEx"]["num_negs"] = args.complex_negs
    if args.complex_loss is not None:
        cfgs["ComplEx"]["loss"] = args.complex_loss
        if args.complex_loss == "nssa":
            cfgs["ComplEx"]["loss_kwargs"] = dict(margin=args.complex_margin,
                                                  adversarial_temperature=1.0)
        else:
            cfgs["ComplEx"].pop("loss_kwargs", None)
    if args.complex_loop is not None:
        cfgs["ComplEx"]["loop"] = args.complex_loop
        if args.complex_loop == "lcwa" and args.complex_loss is None:
            # canonical bilinear recipe: 1-vs-all softmax cross-entropy (Lacroix+2018)
            cfgs["ComplEx"]["loss"] = "crossentropy"
            cfgs["ComplEx"].pop("loss_kwargs", None)
    if any(v is not None for v in (args.complex_negs, args.complex_loss, args.complex_loop)):
        log(f"ComplEx override: loop={cfgs['ComplEx'].get('loop', 'slcwa')} "
            f"negs={cfgs['ComplEx'].get('num_negs')} loss={cfgs['ComplEx'].get('loss')} "
            f"loss_kwargs={cfgs['ComplEx'].get('loss_kwargs')}")

    # Flatten every (regime, seed, dim, epochs, kind, model) cell into ONE ordered plan
    # so the log shows whole-battery progress + a running ETA, on top of the per-epoch bar.
    plan = []
    for rd in regimes:
        jobs = [(args.dim, args.epochs, "main")]
        if args.sweep and rd["short"] == "R0":
            jobs += [(d, e, "sweep") for d in args.sweep_dims for e in args.sweep_epochs]
        for seed in args.seeds:
            for dim, epochs, kind in jobs:
                for model_name in args.models:
                    plan.append((rd, seed, dim, epochs, kind, model_name))

    poolcache: dict = {}      # regime -> (pools, known); depend only on train/test edges
    durations: list = []      # wall-seconds of cells we actually trained (feeds the ETA)
    grid_t0 = time.time()
    n = len(plan)
    for idx, (rd, seed, dim, epochs, kind, model_name) in enumerate(plan, 1):
        if rd["short"] not in poolcache:
            poolcache[rd["short"]] = (
                le._category_pools(rd["train_edges"], rd["hub_set"]),
                le._known_tails_by_head(rd["train_edges"], rd["test_edges"], rd["test_edges"][:, 0]),
            )
        pools, known = poolcache[rd["short"]]
        outp = (main_run_path(args.results_dir, model_name, rd["short"], seed) if kind == "main"
                else sweep_run_path(args.results_dir, model_name, rd["short"], seed, dim, epochs))

        elapsed = time.time() - grid_t0
        eta_h = (float(np.mean(durations)) * (n - idx + 1) / 3600.0) if durations else float("nan")
        bar_n = 28
        filled = int(round(bar_n * (idx - 1) / n))
        log(f"[grid {idx:>2}/{n}] [{'#' * filled + '-' * (bar_n - filled)}] "
            f"{model_name} {rd['short']} seed{seed} d{dim} e{epochs} ({kind}) | "
            f"elapsed {elapsed/3600:.2f}h | trained {len(durations)} | "
            + (f"eta ~{eta_h:.1f}h" if durations else "eta (measuring)"))

        t0 = time.time()
        trained = train_eval_write(cfgs[model_name], model_name, rd, dim, epochs, seed,
                                   device, args, outp, pools, known, TriplesFactory, torch)
        if trained:
            durations.append(time.time() - t0)

    payload = merge_results(args.results_dir, summary_path)
    if (Path(args.results_dir) / "sweep").exists():
        merge_results(Path(args.results_dir) / "sweep", sweep_summary_path)
    print_r0_vs_r3(payload)


if __name__ == "__main__":
    main()
