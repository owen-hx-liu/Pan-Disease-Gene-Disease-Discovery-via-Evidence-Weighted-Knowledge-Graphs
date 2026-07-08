#!/usr/bin/env python3
"""
kge_benchmark.py -- Honest, publication-grade link-prediction benchmark for the
biomedical KG.

Replaces the train==test evaluation in TransE.py with a proper held-out split and
reports STANDARD metrics WITH BASELINES under one comparable protocol:

  Ranking (filtered, both-sides, all-entity -- PyKEEN gold standard for KGE):
      MRR, Hits@1, Hits@3, Hits@10
  Classification (test positives vs. negatives, SHARED across all methods so the
  baselines and KGE models are directly comparable):
      AUROC, AUPRC

  Models:    TransE, DistMult, ComplEx, RotatE   (PyKEEN, GPU if available)
  Baselines: Random, Common-Neighbors, Adamic-Adar   (graph topology)

After training, the best model (by MRR) is used to generate FILTERED novel
predictions for a target gene->disease relation (default BIOLINK:CAUSES):
self-loops, already-known edges, and type-incompatible pairs are removed before
ranking (Tier-1 requirement; the old transe_new_predictions.csv contained
self-loops such as MONDO:0006877 -> MONDO:0006877).

Scientific notes (see CLAUDE.md):
  * Split is RANDOM. The repo has NO real edge-discovery dates (temporalchart.py
    simulates years with random.choices), so a time-based split would be
    meaningless. --year-col enables a real time split if a genuine date column is
    ever added.
  * Classification negatives: --neg-mode {random,type} . Random negatives make
    AUPRC optimistic (most node pairs are trivially non-edges); type-matched
    negatives corrupt the tail only with entities of the SAME category, a much
    harder and more credible test. The run reports BOTH by default.
  * Fixed seed, relative paths (scripts/config.py), no secrets.

Usage (full environment with PyTorch/PyKEEN):
    python scripts/kge_benchmark.py --epochs 100 --dim 128
    python scripts/kge_benchmark.py --skip-kge            # baselines only
"""
import argparse, json, math, os, random, time
import numpy as np
import pandas as pd

try:
    from kg_categories import category_of
except ImportError:  # allow running from repo root
    from scripts.kg_categories import category_of


def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


# --------------------------------------------------------------------------- #
# Data loading / splitting
# --------------------------------------------------------------------------- #
def load_triples(path, year_col=None):
    cols = ["source_id", "relation", "target_id"]
    usecols = cols + ([year_col] if year_col else [])
    df = pd.read_csv(path, usecols=usecols, dtype=str)
    df = df[df["source_id"] != df["target_id"]]            # drop self-loops
    df = df.dropna(subset=cols).drop_duplicates(subset=cols)
    return df.rename(columns={"source_id": "h", "relation": "r", "target_id": "t"})


def split_triples(df, seed, test_frac=0.1, valid_frac=0.1, year_col=None, split_year=None):
    if year_col and split_year is not None:
        yr = pd.to_numeric(df[year_col], errors="coerce")
        train = df[yr < split_year]; held = df[yr >= split_year]
        valid = held.sample(frac=0.5, random_state=seed); test = held.drop(valid.index)
        log(f"TIME split @ {split_year}: train={len(train)} valid={len(valid)} test={len(test)}")
    else:
        rng = np.random.default_rng(seed)
        idx = rng.permutation(len(df))
        n_test = int(len(df) * test_frac); n_val = int(len(df) * valid_frac)
        test = df.iloc[idx[:n_test]]; valid = df.iloc[idx[n_test:n_test + n_val]]
        train = df.iloc[idx[n_test + n_val:]]
        log(f"RANDOM split: train={len(train)} valid={len(valid)} test={len(test)}")
    cols = ["h", "r", "t"]
    return (train[cols].to_numpy(), valid[cols].to_numpy(), test[cols].to_numpy())


# --------------------------------------------------------------------------- #
# Metric helpers (no sklearn dependency; unit-checked)
# --------------------------------------------------------------------------- #
def auroc(pos, neg):
    s = np.concatenate([pos, neg])
    y = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])
    order = np.argsort(s, kind="stable"); ranks = np.empty(len(s)); sx = s[order]
    i = 0
    while i < len(s):
        j = i
        while j + 1 < len(s) and sx[j + 1] == sx[i]:
            j += 1
        ranks[order[i:j + 1]] = (i + j) / 2.0 + 1.0; i = j + 1
    n_p, n_n = len(pos), len(neg)
    return float((ranks[y == 1].sum() - n_p * (n_p + 1) / 2.0) / (n_p * n_n))


def auprc(pos, neg):
    s = np.concatenate([pos, neg])
    y = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])
    order = np.argsort(-s, kind="stable"); y = y[order]
    tp = np.cumsum(y); fp = np.cumsum(1 - y)
    prec = tp / (tp + fp); rec = tp / y.sum()
    rec = np.concatenate([[0], rec]); prec = np.concatenate([[1], prec])
    return float(np.sum((rec[1:] - rec[:-1]) * prec[1:]))


# --------------------------------------------------------------------------- #
# Shared classification eval sets (positives + negatives) -- used by ALL methods
# --------------------------------------------------------------------------- #
def build_classification_sets(train, valid, test, seed, n_pos, neg_per_pos):
    """Sample test positives and build random + type-matched negatives, filtered
    against ALL known triples. Returns dict with positives and both negative sets.
    Each entry is an array of (h, r, t) string triples."""
    rng = np.random.default_rng(seed)
    known = set(map(tuple, np.concatenate([train, valid, test])[:, :3]))
    entities = np.unique(np.concatenate([train[:, 0], train[:, 2]]))
    # bucket entities by category for type-matched negatives
    cats = np.array([category_of(e) for e in entities])
    by_cat = {c: entities[cats == c] for c in np.unique(cats)}

    te = test if len(test) <= n_pos else test[rng.choice(len(test), n_pos, replace=False)]
    pos = te

    def sample_negs(type_matched):
        negs = []
        for h, r, t in te:
            pool = by_cat.get(category_of(t)) if type_matched else entities
            if pool is None or len(pool) == 0:
                pool = entities
            got = 0; tries = 0
            while got < neg_per_pos and tries < neg_per_pos * 40:
                tries += 1
                tt = pool[rng.integers(0, len(pool))]
                if tt == h or (h, r, tt) in known:
                    continue
                negs.append((h, r, tt)); got += 1
        return np.array(negs, dtype=object)

    return {"pos": pos, "neg_random": sample_negs(False), "neg_type": sample_negs(True),
            "known": known}


# --------------------------------------------------------------------------- #
# Topological baselines on the SHARED classification sets + sampled ranking
# --------------------------------------------------------------------------- #
def baseline_metrics(train, evalsets, seed, neg_modes, rank_neg=50, rank_eval=4000):
    import networkx as nx
    rng = np.random.default_rng(seed)
    log("baselines: building train graph (networkx)...")
    G = nx.Graph(); G.add_edges_from(train[:, [0, 2]])
    nodes = np.array(list(G.nodes())); known = evalsets["known"]

    def cn(u, v):
        if u not in G or v not in G: return 0.0
        return float(len(set(G[u]) & set(G[v])))

    def aa(u, v):
        if u not in G or v not in G: return 0.0
        return float(sum(1.0 / math.log(G.degree(w))
                         for w in (set(G[u]) & set(G[v])) if G.degree(w) > 1))

    scorers = {"CommonNeighbors": cn, "AdamicAdar": aa, "Random": None}
    out = {}
    for name, fn in scorers.items():
        row = {}
        # classification AUROC/AUPRC on shared sets (per negative mode)
        for mode in neg_modes:
            negs = evalsets["neg_random" if mode == "random" else "neg_type"]
            if fn is None:
                ps = rng.random(len(evalsets["pos"])); ns = rng.random(len(negs))
            else:
                ps = np.array([fn(h, t) for h, _, t in evalsets["pos"]])
                ns = np.array([fn(h, t) for h, _, t in negs])
            row[f"AUROC_{mode}"] = auroc(ps, ns); row[f"AUPRC_{mode}"] = auprc(ps, ns)
        # sampled filtered ranking (50 negatives) -> MRR / Hits
        te = evalsets["pos"]
        mask = np.array([G.has_node(h) and G.has_node(t) for h, _, t in te])
        ev = te[mask]
        if len(ev) > rank_eval:
            ev = ev[rng.choice(len(ev), rank_eval, replace=False)]
        rr = []; hits = {1: 0, 3: 0, 10: 0}
        for h, r, t in ev:
            ts = rng.random() if fn is None else fn(h, t)
            cand = nodes[rng.integers(0, len(nodes), rank_neg * 3)]
            cand = [c for c in cand if c != h and (h, r, c) not in known][:rank_neg]
            ns = rng.random(len(cand)) if fn is None else np.array([fn(h, c) for c in cand])
            rank = 1 + int(np.sum(ns > ts)) + 0.5 * int(np.sum(ns == ts))
            rr.append(1.0 / rank)
            for k in hits: hits[k] += rank <= k
        n = len(ev)
        row.update({"MRR": float(np.mean(rr)), "Hits@1": hits[1] / n,
                    "Hits@3": hits[3] / n, "Hits@10": hits[10] / n,
                    "rank_protocol": f"sampled-{rank_neg}neg"})
        out[name] = row
        log(f"  {name}: { {k: round(v,3) for k,v in row.items() if isinstance(v,float)} }")
    return out


# --------------------------------------------------------------------------- #
# KGE models via PyKEEN
# --------------------------------------------------------------------------- #
def kge_metrics(train, valid, test, evalsets, models, dim, epochs, batch, seed,
                neg_modes, rank_eval, out_dir):
    import torch
    from pykeen.triples import TriplesFactory
    from pykeen.models import TransE, DistMult, ComplEx, RotatE
    from pykeen.training import SLCWATrainingLoop
    from pykeen.evaluation import RankBasedEvaluator

    device = "cuda" if torch.cuda.is_available() else "cpu"
    log(f"KGE device = {device}")
    model_cls = {"TransE": TransE, "DistMult": DistMult, "ComplEx": ComplEx, "RotatE": RotatE}

    tf_train = TriplesFactory.from_labeled_triples(train)
    e2id, r2id = tf_train.entity_to_id, tf_train.relation_to_id

    def map_triples(arr):
        keep = [(e2id[h], r2id[r], e2id[t]) for h, r, t in arr
                if h in e2id and r in r2id and t in e2id]
        return torch.tensor(keep, dtype=torch.long)

    valid_mt = map_triples(valid); test_mt = map_triples(test)
    # ranking eval on a sample of test triples (filtered against ALL known)
    rng = np.random.default_rng(seed)
    if len(test_mt) > rank_eval:
        sel = rng.choice(len(test_mt), rank_eval, replace=False)
        test_rank = test_mt[sel]
    else:
        test_rank = test_mt
    all_known_mt = torch.cat([map_triples(train), valid_mt, test_mt], dim=0)

    # pre-map shared classification sets to ids (skip any with unknown ids)
    def map_class(arr):
        return np.array([(e2id[h], r2id[r], e2id[t]) for h, r, t in arr
                         if h in e2id and r in r2id and t in e2id], dtype=np.int64)
    cls = {"pos": map_class(evalsets["pos"]),
           "neg_random": map_class(evalsets["neg_random"]),
           "neg_type": map_class(evalsets["neg_type"])}

    def score_batched(model, triples_np, bs=20000):
        model.eval(); outs = []
        with torch.no_grad():
            for i in range(0, len(triples_np), bs):
                batch_t = torch.tensor(triples_np[i:i + bs], dtype=torch.long, device=device)
                outs.append(model.score_hrt(batch_t).cpu().numpy().ravel())
        return np.concatenate(outs) if outs else np.array([])

    # resume: load already-completed models so interruptions (e.g. the machine
    # sleeping mid-run) don't force re-training finished models.
    results = {}
    partial_path = os.path.join(out_dir, "kge_partial.json")
    if os.path.exists(partial_path):
        try:
            results = json.load(open(partial_path))
            if results:
                log(f"resuming: {list(results)} already done, skipping their training")
        except Exception:
            results = {}
    best = (None, -1.0)
    for name in models:
        if name in results:
            continue
        log(f"=== training {name} (dim={dim}, epochs={epochs}, bs={batch}) ===")
        torch.manual_seed(seed); np.random.seed(seed); random.seed(seed)
        model = model_cls[name](triples_factory=tf_train, embedding_dim=dim,
                                random_seed=seed).to(device)
        loop = SLCWATrainingLoop(model=model, triples_factory=tf_train)
        t0 = time.time()
        loop.train(triples_factory=tf_train, num_epochs=epochs, batch_size=batch,
                   use_tqdm=False)
        log(f"  trained in {time.time()-t0:.0f}s")

        # batch_size=512 is proven safe on the 8GB GPU for filtered all-entity
        # ranking (1024 risks OOM). AMO in this PyKEEN version is a kwarg of
        # evaluate(), not the constructor.
        evaluator = RankBasedEvaluator()
        rr = evaluator.evaluate(model=model, mapped_triples=test_rank,
                                additional_filter_triples=[all_known_mt],
                                batch_size=512, use_tqdm=False)
        def g(metric): return float(rr.get_metric(metric))
        row = {"MRR": g("both.realistic.inverse_harmonic_mean_rank"),
               "Hits@1": g("both.realistic.hits_at_1"),
               "Hits@3": g("both.realistic.hits_at_3"),
               "Hits@10": g("both.realistic.hits_at_10"),
               "rank_protocol": f"filtered-all-entity (n={len(test_rank)})"}
        ps = score_batched(model, cls["pos"])
        for mode in neg_modes:
            ns = score_batched(model, cls["neg_random" if mode == "random" else "neg_type"])
            row[f"AUROC_{mode}"] = auroc(ps, ns); row[f"AUPRC_{mode}"] = auprc(ps, ns)
        results[name] = row
        log(f"  {name}: { {k: round(v,3) for k,v in row.items() if isinstance(v,float)} }")
        if row["MRR"] > best[1]:
            best = (name, row["MRR"], model, tf_train)
        # checkpoint after each model so a later crash doesn't lose finished ones
        json.dump(results, open(os.path.join(out_dir, "kge_partial.json"), "w"), indent=2)
        del model
        if device == "cuda":
            torch.cuda.empty_cache()
    return results, best


# --------------------------------------------------------------------------- #
# Filtered novel-prediction generation from the best model
# --------------------------------------------------------------------------- #
def generate_predictions(best, train, valid, test, relation, head_cat, tail_cat,
                         top_k, max_heads, seed, out_csv):
    import torch
    name, _mrr, model, tf = best
    e2id = tf.entity_to_id; r2id = tf.relation_to_id
    id2e = {v: k for k, v in e2id.items()}
    if relation not in r2id:
        log(f"prediction: relation {relation} not in graph; skipping")
        return None
    device = next(model.parameters()).device
    rid = r2id[relation]
    ents = np.array(list(e2id.keys()))
    heads = np.array([e for e in ents if category_of(e) == head_cat])
    tails = np.array([e for e in ents if category_of(e) == tail_cat])
    if len(heads) == 0 or len(tails) == 0:
        log(f"prediction: no {head_cat} heads or {tail_cat} tails; skipping")
        return None
    rng = np.random.default_rng(seed)
    if len(heads) > max_heads:
        heads = heads[rng.choice(len(heads), max_heads, replace=False)]
    known = set(map(tuple, np.concatenate([train, valid, test])[:, :3]))
    tail_ids = torch.tensor([e2id[t] for t in tails], dtype=torch.long, device=device)
    tail_labels = tails
    log(f"prediction: scoring {len(heads)} {head_cat} heads x {len(tails)} {tail_cat} "
        f"tails for relation {relation} (model={name})")
    rows = []
    model.eval()
    with torch.no_grad():
        for h in heads:
            hid = e2id[h]
            hr = torch.tensor([[hid, rid]], dtype=torch.long, device=device)
            scores = model.score_t(hr).cpu().numpy().ravel()  # over all entities
            s_tail = scores[tail_ids.cpu().numpy()]
            order = np.argsort(-s_tail)
            taken = 0
            for j in order:
                t = tail_labels[j]
                if t == h or (h, relation, t) in known:
                    continue
                rows.append((h, relation, t, float(s_tail[j])))
                taken += 1
                if taken >= 50:   # keep top-50 per head, global top-k filtered after
                    break
    pred = pd.DataFrame(rows, columns=["head", "relation", "tail", "score"])
    pred = pred.sort_values("score", ascending=False).head(top_k).reset_index(drop=True)
    pred.insert(0, "rank", pred.index + 1)
    pred["head_category"] = head_cat; pred["tail_category"] = tail_cat
    pred.to_csv(out_csv, index=False)
    log(f"prediction: wrote {len(pred)} filtered novel predictions -> {out_csv}")
    return out_csv


# --------------------------------------------------------------------------- #
def write_table(results, path, neg_modes):
    cols = ["MRR", "Hits@1", "Hits@3", "Hits@10"]
    for m in neg_modes:
        cols += [f"AUROC_{m}", f"AUPRC_{m}"]
    hdr = "| Method | " + " | ".join(cols) + " |"
    sep = "|" + "---|" * (len(cols) + 1)
    lines = [hdr, sep]
    order = ["Random", "CommonNeighbors", "AdamicAdar", "TransE", "DistMult", "ComplEx", "RotatE"]
    for name in order:
        if name not in results: continue
        r = results[name]
        cells = [f"{r.get(c, float('nan')):.3f}" for c in cols]
        lines.append(f"| {name} | " + " | ".join(cells) + " |")
    open(path, "w").write("\n".join(lines) + "\n")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--edges", default="data/processed/edges_clean.csv")
    ap.add_argument("--out", default="data/processed/benchmark")
    ap.add_argument("--models", nargs="+", default=["TransE", "DistMult", "ComplEx", "RotatE"])
    ap.add_argument("--dim", type=int, default=128)
    ap.add_argument("--epochs", type=int, default=100)
    ap.add_argument("--batch", type=int, default=16384)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--neg-modes", nargs="+", default=["random", "type"])
    ap.add_argument("--n-class", type=int, default=10000, help="# test positives for AUROC/AUPRC")
    ap.add_argument("--neg-per-pos", type=int, default=1)
    ap.add_argument("--rank-eval", type=int, default=10000, help="# test triples for ranking")
    ap.add_argument("--year-col", default=None)
    ap.add_argument("--split-year", type=int, default=2020)
    ap.add_argument("--skip-kge", action="store_true")
    ap.add_argument("--skip-baselines", action="store_true")
    # prediction generation
    ap.add_argument("--pred-relation", default="BIOLINK:CAUSES")
    ap.add_argument("--pred-head-cat", default="gene")
    ap.add_argument("--pred-tail-cat", default="disease")
    ap.add_argument("--pred-top-k", type=int, default=500)
    ap.add_argument("--pred-max-heads", type=int, default=3000)
    ap.add_argument("--no-predict", action="store_true")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    random.seed(args.seed); np.random.seed(args.seed)
    try:
        from keep_awake import keep_awake; keep_awake()   # prevent sleep killing the run
    except Exception:
        pass

    df = load_triples(args.edges, args.year_col)
    log(f"loaded {len(df)} unique non-self-loop triples")
    train, valid, test = split_triples(df, args.seed, year_col=args.year_col,
                                       split_year=args.split_year)

    log(f"building shared classification eval sets (n_pos={args.n_class}, "
        f"neg/pos={args.neg_per_pos}, modes={args.neg_modes})")
    evalsets = build_classification_sets(train, valid, test, args.seed,
                                         args.n_class, args.neg_per_pos)

    results = {}
    meta = {"seed": args.seed, "dim": args.dim, "epochs": args.epochs,
            "batch": args.batch, "n_triples": len(df),
            "n_train": len(train), "n_valid": len(valid), "n_test": len(test),
            "neg_modes": args.neg_modes, "n_class": args.n_class,
            "neg_per_pos": args.neg_per_pos, "rank_eval": args.rank_eval}

    if not args.skip_baselines:
        results.update(baseline_metrics(train, evalsets, args.seed, args.neg_modes,
                                        rank_eval=args.rank_eval))
        json.dump({"baselines": results, "meta": meta},
                  open(os.path.join(args.out, "benchmark_results.json"), "w"), indent=2)

    best = None
    if not args.skip_kge:
        kge, best = kge_metrics(train, valid, test, evalsets, args.models, args.dim,
                                args.epochs, args.batch, args.seed, args.neg_modes,
                                args.rank_eval, args.out)
        results.update(kge)

    json.dump({"results": results, "meta": meta},
              open(os.path.join(args.out, "benchmark_results.json"), "w"), indent=2)
    table = write_table(results, os.path.join(args.out, "benchmark_table.md"), args.neg_modes)
    log("wrote benchmark_results.json and benchmark_table.md")
    print("\n" + table + "\n")

    if best and best[0] is not None and len(best) >= 4 and not args.no_predict:
        generate_predictions(best, train, valid, test, args.pred_relation,
                             args.pred_head_cat, args.pred_tail_cat, args.pred_top_k,
                             args.pred_max_heads, args.seed,
                             os.path.join(args.out, "novel_predictions_filtered.csv"))
    elif not args.no_predict:
        log("predictions skipped: best model was resumed (no live model object). "
            "Re-run a single model fresh, or use --models <best> to regenerate.")


if __name__ == "__main__":
    main()
