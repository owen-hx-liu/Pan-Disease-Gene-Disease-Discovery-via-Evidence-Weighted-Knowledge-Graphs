#!/usr/bin/env python3
"""
build_degree_null.py -- Unit 7: the R2 degree-preserving null training graph.

Writes data/processed/splits/train_R2_degree_null_seed<seed>.csv: a randomized copy of
the R0 training graph (train.csv) in which every node's in- and out-degree is preserved
EXACTLY, per relation, while the specific head->tail wiring is destroyed. Training a model
on this null and evaluating it on the REAL held-out test edges isolates how much of a
method's score is explained by degree structure alone -- the R0 - R2 drop is the signal
that is NOT attributable to degree. (On this task pure PreferentialAttachment already
scores MRR ~0.59, so degree carries a lot; R2 is the control that quantifies it.)

Method -- per-relation directed double-edge swap (Maslov & Sneppen 2002):
    For two edges (h1, r, t1) and (h2, r, t2) of the SAME relation r, propose
    (h1, r, t2) and (h2, r, t1); accept iff it makes no self-loop and no parallel edge.
    Only tail endpoints are permuted among a relation's edges, so the head multiset is
    untouched (head out-degree preserved) and the tail multiset is a permutation (tail
    in-degree preserved) -- degrees are preserved by construction, exactly, and the
    relation-type histogram is identical to R0. ~--swaps-per-edge x |E_r| swaps per
    relation give good mixing.

Only the TRAINING graph is randomized; valid.csv / test.csv are untouched, so the held-out
target edges are unchanged. The seed is fixed and split_manifest.json is updated with the
new file's edge count + sha256 so the split set stays reproducible.

Live progress: a line-based bar (a fresh line at most every --progress-secs seconds) over
total swap attempts, with accept-rate + ETA. It is line-based on purpose so it renders when
the log is tailed (Get-Content <log> -Wait) even under a backgrounded run.

Usage:
    python scripts/build_degree_null.py                          # full build, seed 42
    python scripts/build_degree_null.py --swaps-per-edge 10
    python scripts/build_degree_null.py --limit-relations 3 --dry-run   # fast calibration
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from config import PROCESSED_DIR
    SPLITS_DIR = Path(PROCESSED_DIR) / "splits"
except Exception:
    SPLITS_DIR = Path("data/processed/splits")


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


class Progress:
    """Line-based, wall-clock-throttled progress bar over swap attempts."""

    def __init__(self, total, tag="swap", secs=1.0):
        self.total = max(int(total), 1)
        self.tag, self.secs = tag, secs
        self.t0 = time.time()
        self.last = 0.0
        self.done = 0
        self.accepted = 0

    def add(self, k, accepted):
        self.done += k
        self.accepted += accepted
        now = time.time()
        if now - self.last >= self.secs or self.done >= self.total:
            self.last = now
            self._emit()

    def _emit(self):
        frac = min(max(self.done / self.total, 1e-9), 1.0)
        el = time.time() - self.t0
        eta = el / frac - el
        w = 28
        fill = int(round(w * frac))
        bar = "#" * fill + "-" * (w - fill)
        acc = 100.0 * self.accepted / max(self.done, 1)
        log(f"[{self.tag}] [{bar}] {frac * 100:5.1f}% "
            f"swaps {self.done:,}/{self.total:,} accept={acc:4.1f}% "
            f"elapsed={el:.0f}s eta={eta:.0f}s")


def swap_relation(head_ids, tail_ids, n_nodes, n_swaps, rng, prog):
    """Degree-preserving double-edge swaps on one relation.

    head_ids/tail_ids are int32 arrays (global node ids) for the edges of one relation.
    Returns a new tails array (heads are never moved). Packs (h, t) -> h*n_nodes + t
    (fits int64: n_nodes ~ 4.5e5, so h*n_nodes+t < 2e11 << 9.2e18).
    """
    m = len(head_ids)
    heads = head_ids
    tails = tail_ids.copy()
    if m < 2 or n_swaps <= 0:
        prog.add(int(n_swaps), 0)
        return tails, 0

    N = n_nodes
    edgeset = set((int(h) * N + int(t)) for h, t in zip(heads.tolist(), tails.tolist()))
    accepted = 0
    CH = 100_000
    done = 0
    while done < n_swaps:
        k = int(min(CH, n_swaps - done))
        I = rng.integers(0, m, size=k)
        J = rng.integers(0, m, size=k)
        acc_chunk = 0
        for a in range(k):
            i = int(I[a]); j = int(J[a])
            if i == j:
                continue
            h1 = int(heads[i]); t1 = int(tails[i])
            h2 = int(heads[j]); t2 = int(tails[j])
            if h1 == h2 or t1 == t2:      # degenerate (no-op / guaranteed parallel)
                continue
            if h1 == t2 or h2 == t1:      # would create a self-loop
                continue
            e1 = h1 * N + t2
            e2 = h2 * N + t1
            if e1 in edgeset or e2 in edgeset:   # would create a parallel edge
                continue
            edgeset.discard(h1 * N + t1)
            edgeset.discard(h2 * N + t2)
            edgeset.add(e1)
            edgeset.add(e2)
            tails[i] = t2
            tails[j] = t1
            accepted += 1
            acc_chunk += 1
        done += k
        prog.add(k, acc_chunk)
    return tails, accepted


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--splits-dir", default=str(SPLITS_DIR))
    ap.add_argument("--train", default="train.csv", help="R0 training graph to randomize")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--swaps-per-edge", type=int, default=10,
                    help="double-edge swap attempts per edge, per relation (mixing)")
    ap.add_argument("--progress-secs", type=float, default=1.0,
                    help="emit a progress line at most this often (seconds)")
    ap.add_argument("--limit-relations", type=int, default=None,
                    help="only process the N largest relations (fast calibration)")
    ap.add_argument("--dry-run", action="store_true",
                    help="swap + verify but write NO file and do NOT touch the manifest")
    args = ap.parse_args()
    splits = Path(args.splits_dir)
    rng = np.random.default_rng(args.seed)

    log(f"loading {args.train} ...")
    df = pd.read_csv(splits / args.train, usecols=["source_id", "relation", "target_id"], dtype=str)
    df = df.dropna(subset=["source_id", "relation", "target_id"])
    n0 = len(df)
    # Exact-duplicate (h,r,t) rows would corrupt the per-relation edge set; drop them
    # (a clean split has none). Report so any count change is explicit in the manifest.
    df = df.drop_duplicates(subset=["source_id", "relation", "target_id"]).reset_index(drop=True)
    if len(df) != n0:
        log(f"  dropped {n0 - len(df):,} exact-duplicate triples (kept {len(df):,})")

    # Integer-encode node labels once (shared head/tail id space).
    nodes = pd.unique(pd.concat([df["source_id"], df["target_id"]], ignore_index=True))
    id_of = {lbl: i for i, lbl in enumerate(nodes)}
    n_nodes = len(nodes)
    label_of = np.asarray(nodes, dtype=object)
    head_ids = df["source_id"].map(id_of).to_numpy(np.int64)
    tail_ids = df["target_id"].map(id_of).to_numpy(np.int64)
    rels = df["relation"].to_numpy(object)
    log(f"graph: {len(df):,} edges, {n_nodes:,} nodes, {len(np.unique(rels)):,} relations")

    # Order relations by size (largest first); optionally limit for calibration.
    uniq, counts = np.unique(rels, return_counts=True)
    order = np.argsort(-counts)
    uniq, counts = uniq[order], counts[order]
    if args.limit_relations:
        uniq, counts = uniq[:args.limit_relations], counts[:args.limit_relations]
        log(f"  calibration: only the {len(uniq)} largest relation(s)")

    total_swaps = int(sum(args.swaps_per_edge * c for c in counts if c >= 2))
    prog = Progress(total_swaps, tag="swap", secs=args.progress_secs)
    log(f"target swaps: {total_swaps:,} ({args.swaps_per_edge}x per edge, per relation)")

    new_tails = tail_ids.copy()
    total_accepted = 0
    for rel, c in zip(uniq.tolist(), counts.tolist()):
        idx = np.nonzero(rels == rel)[0]
        h = head_ids[idx]
        t = tail_ids[idx]
        # guard against pre-existing intra-relation parallels (would break the edge set)
        if len(idx) != len({(int(a), int(b)) for a, b in zip(h.tolist(), t.tolist())}):
            log(f"  NOTE relation {rel} has intra-relation parallels; edge set dedups them")
        swapped, acc = swap_relation(h, t, n_nodes, args.swaps_per_edge * c, rng, prog)
        new_tails[idx] = swapped
        total_accepted += acc
    log(f"done: {total_accepted:,} swaps accepted "
        f"({100.0 * total_accepted / max(total_swaps, 1):.1f}% of attempts)")

    # Degree preservation is guaranteed by construction (heads are never moved; tails are
    # only permuted within a relation) -- assert the tail multiset really is a permutation
    # as a hard sanity check that in-degree is preserved exactly.
    assert np.array_equal(np.bincount(tail_ids, minlength=n_nodes),
                          np.bincount(new_tails, minlength=n_nodes)), "tail in-degree changed"
    # And confirm we actually rewired something (unless a tiny calibration).
    frac_moved = float(np.mean(tail_ids != new_tails))
    log(f"verified: in/out degree preserved exactly; {frac_moved * 100:.1f}% of edges rewired")

    if args.dry_run:
        log("dry-run: no file written, manifest untouched. "
            f"(extrapolated full build ~ scale by total_edges/processed_edges.)")
        return

    out_name = f"train_R2_degree_null_seed{args.seed}.csv"
    out_path = splits / out_name
    out = pd.DataFrame({
        "source_id": label_of[head_ids],
        "relation": rels,
        "target_id": label_of[new_tails],
    })
    out.to_csv(out_path, index=False)
    sha = hashlib.sha256(out_path.read_bytes()).hexdigest()
    log(f"wrote {out_path} ({len(out):,} edges)  sha256={sha[:16]}...")

    # Update the manifest: edge count + sha256 (leave everything else intact).
    man_path = splits / "split_manifest.json"
    man = json.loads(man_path.read_text())
    man.setdefault("counts", {})["train_R2_degree_null"] = int(len(out))
    man.setdefault("sha256", {})[out_name] = sha
    man.setdefault("params", {})["R2_swaps_per_edge"] = args.swaps_per_edge
    man.setdefault("params", {})["R2_seed"] = args.seed
    man_path.write_text(json.dumps(man, indent=2))
    log(f"updated manifest: counts.train_R2_degree_null={len(out):,}, sha256[{out_name}]")


if __name__ == "__main__":
    main()
