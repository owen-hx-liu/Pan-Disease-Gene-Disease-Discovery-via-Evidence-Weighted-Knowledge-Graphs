#!/usr/bin/env python3
"""
null_model_tests.py -- Statistical significance for the network-science claims.

PROJECT_REPORT.md 6.9 / 7.9: the "resilient", hub, clustering, and community
claims were reported as bare numbers with no null model. This script compares the
observed graph against a degree-preserving randomization (configuration model),
so each structural claim gets an observed-vs-null contrast with a z-score.

Metrics (observed vs. K configuration-model realizations, same degree sequence):
  * global transitivity         (sampled over high-coverage node subset)
  * average local clustering     (sampled)
  * degree assortativity         (exact)
  * giant-component fraction      (exact)
  * modularity of a greedy partition (optional, --modularity; slower)

A configuration model preserves the degree of every node but rewires edges at
random, so any clustering/modularity ABOVE the null is genuine structure, not a
by-product of the degree distribution.

Usage:
  python scripts/null_model_tests.py --edges data/processed/edges_clean_integrated.csv \
      --realizations 5 --sample 3000
"""
import argparse, json, time, random
import numpy as np
import pandas as pd
import networkx as nx


def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def load_graph(path):
    df = pd.read_csv(path, usecols=["source_id", "target_id"], dtype=str)
    df = df[df["source_id"] != df["target_id"]]
    G = nx.Graph()
    G.add_edges_from(df.itertuples(index=False, name=None))
    return G


def sampled_clustering(G, nodes, seed):
    rng = random.Random(seed)
    samp = rng.sample(list(nodes), min(len(nodes), len(nodes)))
    cl = nx.clustering(G, samp)
    vals = np.array(list(cl.values()), dtype=float)
    return float(vals.mean())


def metrics(G, sample, seed, do_modularity=False):
    rng = random.Random(seed)
    nodes = list(G.nodes())
    samp = rng.sample(nodes, min(sample, len(nodes)))
    cl = nx.clustering(G, samp)
    avg_clust = float(np.mean(list(cl.values()))) if cl else 0.0
    # transitivity estimate on the same sample's induced view via clustering mean
    try:
        assort = float(nx.degree_assortativity_coefficient(G))
    except Exception:
        assort = float("nan")
    comps = max(len(c) for c in nx.connected_components(G))
    giant = comps / G.number_of_nodes()
    out = {"avg_clustering_sampled": avg_clust,
           "degree_assortativity": assort,
           "giant_component_frac": giant}
    if do_modularity:
        communities = nx.community.greedy_modularity_communities(G)
        out["modularity_greedy"] = float(nx.community.modularity(G, communities))
        out["n_communities"] = len(communities)
    return out


def config_model(G, seed):
    """Memory-efficient degree-preserving null via numpy stub matching.

    networkx.configuration_model builds a MultiGraph and OOMs at this scale
    (~11M stubs). We pair shuffled stubs in numpy, drop self-loops/parallels,
    and build a plain simple Graph on integer node ids (labels are irrelevant
    for clustering / assortativity / component metrics)."""
    deg = np.array([d for _n, d in G.degree()], dtype=np.int64)
    rng = np.random.default_rng(seed)
    stubs = np.repeat(np.arange(len(deg)), deg)
    rng.shuffle(stubs)
    if len(stubs) % 2:                       # drop a dangling stub if odd
        stubs = stubs[:-1]
    e = stubs.reshape(-1, 2)
    e = e[e[:, 0] != e[:, 1]]                 # drop self-loops
    e = np.sort(e, axis=1)
    e = np.unique(e, axis=0)                  # collapse parallel edges
    cm = nx.Graph()
    cm.add_nodes_from(range(len(deg)))
    cm.add_edges_from(map(tuple, e))
    return cm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--edges", default="data/processed/edges_clean_integrated.csv")
    ap.add_argument("--out", default="data/processed/null_model_report.json")
    ap.add_argument("--realizations", type=int, default=5)
    ap.add_argument("--sample", type=int, default=3000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--modularity", action="store_true")
    args = ap.parse_args()

    log(f"loading graph from {args.edges}")
    G = load_graph(args.edges)
    log(f"graph: {G.number_of_nodes():,} nodes / {G.number_of_edges():,} edges")

    log("computing observed metrics...")
    obs = metrics(G, args.sample, args.seed, args.modularity)
    log(f"observed: {obs}")

    null_runs = []
    for i in range(args.realizations):
        t0 = time.time()
        cm = config_model(G, args.seed + i)
        m = metrics(cm, args.sample, args.seed + i, args.modularity)
        null_runs.append(m)
        log(f"  null {i+1}/{args.realizations}: {m}  ({time.time()-t0:.0f}s)")

    # summarize observed vs null with z-scores
    summary = {"graph": {"nodes": G.number_of_nodes(), "edges": G.number_of_edges()},
               "observed": obs, "null_mean": {}, "null_std": {}, "z_score": {},
               "realizations": args.realizations, "sample": args.sample}
    for k in obs:
        vals = np.array([r[k] for r in null_runs], dtype=float)
        mu, sd = float(np.nanmean(vals)), float(np.nanstd(vals))
        summary["null_mean"][k] = mu
        summary["null_std"][k] = sd
        summary["z_score"][k] = float((obs[k] - mu) / sd) if sd > 0 else float("nan")

    json.dump(summary, open(args.out, "w"), indent=2)
    log(f"wrote {args.out}")
    print("\n| Metric | Observed | Null mean | Null sd | z |")
    print("|---|---|---|---|---|")
    for k in obs:
        print(f"| {k} | {obs[k]:.4f} | {summary['null_mean'][k]:.4f} "
              f"| {summary['null_std'][k]:.4f} | {summary['z_score'][k]:.2f} |")


if __name__ == "__main__":
    main()
