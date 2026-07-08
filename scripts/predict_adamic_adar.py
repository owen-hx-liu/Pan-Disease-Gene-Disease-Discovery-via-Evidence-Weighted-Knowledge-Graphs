#!/usr/bin/env python3
"""
predict_adamic_adar.py -- Novel gene->disease hypotheses from Adamic-Adar.

The benchmark shows Adamic-Adar is the strongest link-prediction method on this
graph (MRR 0.69 vs <=0.09 for every KGE model). So the actual novel-link
hypotheses should come from Adamic-Adar, not the weak KGE models.

For each human gene (HGNC) we score candidate gene->disease (MONDO) pairs that are
NOT already edges by the Adamic-Adar index
    AA(g, d) = sum_{w in N(g) ∩ N(d)} 1 / log(deg(w))
accumulated efficiently over shared neighbors. Generic hubs (degree > --hub-cap)
are skipped as intermediates (they carry little specific signal and dominate cost;
CLAUDE.md recommends down-weighting ultra-generic hubs). Self-loops, already-known
edges, and non gene->disease pairs are filtered out.

Output: data/processed/predictions/adamic_adar_gene_disease.csv  (rank, gene,
disease, adamic_adar, n_shared_neighbors)
"""
import argparse, math, heapq, time
from collections import defaultdict
from pathlib import Path

import pandas as pd
import networkx as nx

from kg_categories import category_of
from config import PROCESSED_DIR


def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--edges", default="data/processed/edges_clean_integrated.csv")
    ap.add_argument("--out", default=str(PROCESSED_DIR / "predictions" / "adamic_adar_gene_disease.csv"))
    ap.add_argument("--head-prefix", default="HGNC", help="restrict gene heads to this prefix")
    ap.add_argument("--top-k", type=int, default=1000)
    ap.add_argument("--hub-cap", type=int, default=2000, help="skip shared neighbors with degree above this")
    ap.add_argument("--max-heads", type=int, default=0, help="0 = all matching heads")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    out = Path(args.out); out.parent.mkdir(parents=True, exist_ok=True)
    try:
        from keep_awake import keep_awake; keep_awake()
    except Exception:
        pass

    log(f"loading {args.edges}")
    df = pd.read_csv(args.edges, usecols=["source_id", "target_id"], dtype=str)
    df = df[df["source_id"] != df["target_id"]]
    G = nx.Graph(); G.add_edges_from(df.itertuples(index=False, name=None))
    log(f"graph: {G.number_of_nodes():,} nodes / {G.number_of_edges():,} edges")
    deg = dict(G.degree())

    genes = [n for n in G if n.split(":", 1)[0].upper() == args.head_prefix.upper()
             and category_of(n) == "gene"]
    if args.max_heads and len(genes) > args.max_heads:
        import random; random.Random(args.seed).shuffle(genes); genes = genes[:args.max_heads]
    log(f"scoring Adamic-Adar for {len(genes):,} {args.head_prefix} gene heads "
        f"(hub-cap deg<={args.hub_cap})")

    heap = []  # min-heap of (score, gene, disease, n_shared)
    for i, g in enumerate(genes):
        if i % 2000 == 0 and i:
            log(f"  {i}/{len(genes)} genes; heap={len(heap)}")
        known = set(G[g])
        aa = defaultdict(float); shared = defaultdict(int)
        for w in G[g]:
            dw = deg[w]
            if dw <= 1 or dw > args.hub_cap:
                continue
            contrib = 1.0 / math.log(dw)
            for d in G[w]:
                if d in known or d == g:
                    continue
                if category_of(d) != "disease":
                    continue
                aa[d] += contrib; shared[d] += 1
        for d, s in aa.items():
            item = (s, g, d, shared[d])
            if len(heap) < args.top_k:
                heapq.heappush(heap, item)
            elif s > heap[0][0]:
                heapq.heapreplace(heap, item)

    rows = sorted(heap, key=lambda x: -x[0])
    pred = pd.DataFrame([(r+1, g, d, round(s, 4), n) for r, (s, g, d, n) in enumerate(rows)],
                        columns=["rank", "gene", "disease", "adamic_adar", "n_shared_neighbors"])
    pred.to_csv(out, index=False)
    log(f"wrote {len(pred)} filtered novel gene->disease predictions -> {out}")
    if len(pred):
        print(pred.head(15).to_string(index=False))


if __name__ == "__main__":
    main()
