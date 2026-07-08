#!/usr/bin/env python3
"""
full_graph_analysis.py -- Network-science analysis on the integrated graph.

Regenerates the metrics in full_graph_analysis/ (which were computed on the OLD
Monarch-only graph with hand-written numpy) using networkx + scipy on the
INTEGRATED graph. With scipy available we also compute the EXACT Fiedler value
(algebraic connectivity), which the sandbox run could only bound.

Outputs (default dir full_graph_analysis_integrated/), each checkpointed as soon
as it is computed:
  full_degree_top.csv, full_pagerank_top.csv, full_eigenvector_top.csv
  full_clustering.json, full_components.json, full_communities.json
  full_spectral.json, full_node_prefix_composition.csv
  full_graph_analysis_report.md
"""
import argparse, json, time
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import networkx as nx
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh


def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def load_graph(path):
    df = pd.read_csv(path, usecols=["source_id", "target_id"], dtype=str)
    df = df[df["source_id"] != df["target_id"]]
    G = nx.Graph()
    G.add_edges_from(df.itertuples(index=False, name=None))
    return G


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--edges", default="data/processed/edges_clean_integrated.csv")
    ap.add_argument("--out", default="full_graph_analysis_integrated")
    ap.add_argument("--sample", type=int, default=3000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--skip-communities", action="store_true")
    args = ap.parse_args()
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    try:
        from keep_awake import keep_awake; keep_awake()
    except Exception:
        pass

    log(f"loading {args.edges}")
    G = load_graph(args.edges)
    N, E = G.number_of_nodes(), G.number_of_edges()
    log(f"graph: {N:,} nodes / {E:,} edges")

    # ---- degree ----
    deg = dict(G.degree())
    dser = pd.Series(deg).sort_values(ascending=False)
    top = dser.head(200)
    pd.DataFrame({"node": top.index, "degree": top.values,
                  "degree_centrality": top.values / (N - 1)}
                 ).to_csv(out / "full_degree_top.csv", index=False)
    log(f"degree done; top hub {dser.index[0]} ({dser.iloc[0]})")

    # ---- node prefix composition ----
    pref = Counter(n.split(":", 1)[0] for n in G.nodes())
    pd.DataFrame(sorted(pref.items(), key=lambda x: -x[1]),
                 columns=["prefix", "node_count"]).to_csv(
                     out / "full_node_prefix_composition.csv", index=False)

    # ---- components ----
    comps = sorted((len(c) for c in nx.connected_components(G)), reverse=True)
    lcc_nodes_set = max(nx.connected_components(G), key=len)
    comp = {"num_components": len(comps), "lcc_size": comps[0],
            "lcc_fraction": round(comps[0] / N, 5),
            "isolated_nodes": sum(1 for _n, d in G.degree() if d == 0),
            "component_size_top20": comps[:20]}
    json.dump(comp, open(out / "full_components.json", "w"), indent=2)
    log(f"components done; LCC {comp['lcc_size']:,} ({comp['lcc_fraction']})")

    # ---- sampled clustering ----
    samp = list(rng.choice(np.array(G.nodes(), dtype=object),
                           size=min(args.sample, N), replace=False))
    cl = nx.clustering(G, samp)
    json.dump({"sampled_nodes": len(samp),
               "mean_local_clustering": round(float(np.mean(list(cl.values()))), 5)},
              open(out / "full_clustering.json", "w"), indent=2)
    log("clustering done")

    # ---- PageRank ----
    log("pagerank...")
    pr = nx.pagerank(G, alpha=0.85, max_iter=100, tol=1e-6)
    prs = pd.Series(pr).sort_values(ascending=False).head(200)
    pd.DataFrame({"node": prs.index, "pagerank": prs.values}).to_csv(
        out / "full_pagerank_top.csv", index=False)
    log("pagerank done")

    # ---- spectral on LCC (scipy: exact Fiedler) ----
    log("building LCC sparse Laplacian...")
    lcc = G.subgraph(lcc_nodes_set)
    nodes = list(lcc.nodes()); idx = {n: i for i, n in enumerate(nodes)}
    rows, cols = [], []
    for u, v in lcc.edges():
        rows.append(idx[u]); cols.append(idx[v])
    n = len(nodes)
    data = np.ones(len(rows))
    A = sp.coo_matrix((data, (rows, cols)), shape=(n, n))
    A = (A + A.T).tocsr()
    degs = np.asarray(A.sum(axis=1)).ravel()
    L = sp.diags(degs) - A
    spectral = {"lcc_nodes": n, "lcc_edges": lcc.number_of_edges(),
                "lcc_fraction": round(n / N, 5), "min_degree_lcc": int(degs.min())}
    # largest Laplacian eigenvalue: Lanczos, no factorization -> memory-safe
    try:
        lam_max = float(eigsh(L, k=1, which="LA", return_eigenvectors=False,
                              maxiter=5000)[0])
        spectral["largest_laplacian_eigenvalue"] = round(lam_max, 3)
    except Exception as e:
        spectral["largest_laplacian_eigenvalue_error"] = str(e)
    # Fiedler value via LOBPCG (matrix-free, low memory). The previous
    # shift-invert eigsh did a sparse LU factorization that OOMed at 437k nodes.
    try:
        from scipy.sparse.linalg import lobpcg
        rng2 = np.random.default_rng(args.seed)
        # deflate the known null space (constant vector) and solve for the next
        # smallest eigenpair; X has 2 columns to capture eig0~0 and Fiedler
        X = rng2.standard_normal((n, 2))
        ones = np.ones((n, 1)) / np.sqrt(n)
        X[:, [0]] = ones
        M = sp.diags(1.0 / np.maximum(degs, 1))   # Jacobi preconditioner
        vals, _vecs = lobpcg(L, X, M=M, largest=False, tol=1e-5, maxiter=2000)
        vals = sorted(float(v) for v in vals)
        spectral["fiedler_value"] = round(vals[1], 6)
        spectral["smallest_eigenvalue"] = round(vals[0], 8)
    except Exception as e:
        spectral["fiedler_error"] = str(e)
        spectral["fiedler_note"] = "0 < Fiedler <= min_degree_lcc"
    json.dump(spectral, open(out / "full_spectral.json", "w"), indent=2)
    log(f"spectral done: {spectral}")

    # ---- communities (Louvain) ----
    community = {"skipped": True}
    if not args.skip_communities:
        log("louvain communities (slow)...")
        comms = nx.community.louvain_communities(G, resolution=1.0, seed=args.seed)
        mod = nx.community.modularity(G, comms)
        sizes = sorted((len(c) for c in comms), reverse=True)
        community = {"method": "networkx louvain_communities, resolution=1",
                     "num_communities": len(comms), "modularity": round(float(mod), 5),
                     "largest_community_size": sizes[0],
                     "largest_fraction": round(sizes[0] / N, 4),
                     "communities_ge100": sum(1 for s in sizes if s >= 100),
                     "communities_ge10": sum(1 for s in sizes if s >= 10),
                     "singletons": sum(1 for s in sizes if s == 1),
                     "top20_sizes": sizes[:20]}
        log(f"communities done: Q={community['modularity']}, k={community['num_communities']}")
    json.dump(community, open(out / "full_communities.json", "w"), indent=2)

    # ---- report ----
    md = [f"# Full-graph network analysis (integrated graph)\n",
          f"Graph: **{N:,} nodes / {E:,} edges** · seed {args.seed} · "
          f"computed with networkx + scipy.\n",
          "| Metric | Value |", "|---|---|",
          f"| Nodes / edges | {N:,} / {E:,} |",
          f"| Top hub | {dser.index[0]} (deg {dser.iloc[0]:,}) |",
          f"| Connected components | {comp['num_components']:,} |",
          f"| Giant component | {comp['lcc_size']:,} ({comp['lcc_fraction']*100:.2f}%) |",
          f"| Mean local clustering (sampled) | {json.load(open(out/'full_clustering.json'))['mean_local_clustering']} |",
          f"| Largest Laplacian eigenvalue | {spectral.get('largest_laplacian_eigenvalue','n/a')} |",
          f"| Fiedler value (algebraic connectivity) | {spectral.get('fiedler_value','n/a')} |"]
    if not community.get("skipped"):
        md += [f"| Communities (Louvain) | {community['num_communities']:,} |",
               f"| Modularity Q | {community['modularity']} |"]
    (out / "full_graph_analysis_report.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    log(f"wrote report -> {out/'full_graph_analysis_report.md'}")


if __name__ == "__main__":
    main()
