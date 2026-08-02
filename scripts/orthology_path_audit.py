#!/usr/bin/env python
"""
Why is R3 (orthology-blocked) a null?

Removing 520,198 cross-species bridges leaves every method's score untouched,
which surprises biologists who rely on model-organism evidence. The obvious
mechanistic question is whether those bridges ever carried information about a
held-out human gene-disease pair in the first place. A cross-species path runs

    human gene --ORTHOLOGOUS_TO--> non-human gene --(any relation)--> disease

so this script counts, for the 4,228 held-out test pairs, how many are reachable
by such a path at all, and of those, how many are *also* reachable by a two-hop
path that uses no cross-species edge. A pair in the second group gains nothing
from the bridge: the bridge is structurally redundant with evidence the training
graph already carries.

Reads the frozen R0 training graph and the frozen test split; writes one JSON.

    python scripts/orthology_path_audit.py
"""
from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPLITS = ROOT / "data" / "processed" / "splits"
OUT = ROOT / "data" / "processed" / "results" / "null" / "orthology_path_audit.json"

TRAIN = SPLITS / "train.csv"          # R0 training graph
TEST = SPLITS / "test.csv"

ORTHO = "BIOLINK:ORTHOLOGOUS_TO"
HUMAN = "HGNC:"
DISEASE = "MONDO:"


def main() -> None:
    # ortholog[g] = set of genes g is declared orthologous to, in either direction
    ortholog: dict[str, set[str]] = defaultdict(set)
    # out_disease[n] = diseases node n points at (or is pointed at by) directly
    disease_nbr: dict[str, set[str]] = defaultdict(set)
    # non-orthology adjacency, for the redundancy check
    plain_nbr: dict[str, set[str]] = defaultdict(set)

    n_edges = 0
    with TRAIN.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            s, rel, t = row["source_id"], row["relation"], row["target_id"]
            n_edges += 1
            if rel == ORTHO:
                ortholog[s].add(t)
                ortholog[t].add(s)
                continue
            plain_nbr[s].add(t)
            plain_nbr[t].add(s)
            if t.startswith(DISEASE):
                disease_nbr[s].add(t)
            if s.startswith(DISEASE):
                disease_nbr[t].add(s)

    test = []
    with TEST.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            test.append((row["source_id"], row["target_id"]))

    n_reachable = 0          # pair reachable through an orthology bridge
    n_reachable_nonhuman = 0 # ... where the bridging ortholog is non-human
    n_also_plain = 0         # ... and also reachable without any orthology edge
    n_plain_only = 0         # reachable only without orthology
    for gene, dis in test:
        orths = ortholog.get(gene, ())
        via_ortho = any(dis in disease_nbr.get(o, ()) for o in orths)
        via_ortho_nonhuman = any(
            not o.startswith(HUMAN) and dis in disease_nbr.get(o, ()) for o in orths)
        # two-hop path using no orthology edge at all
        via_plain = bool(plain_nbr.get(gene, set()) & plain_nbr.get(dis, set()))
        if via_ortho:
            n_reachable += 1
            if via_ortho_nonhuman:
                n_reachable_nonhuman += 1
            if via_plain:
                n_also_plain += 1
        elif via_plain:
            n_plain_only += 1

    # How narrow is the cross-species landing zone? Count the distinct diseases
    # any non-human gene points at directly.
    nonhuman_targets = set()
    for node, ds in disease_nbr.items():
        if not node.startswith(HUMAN) and not node.startswith(DISEASE):
            nonhuman_targets |= ds
    all_diseases = set()
    for ds in disease_nbr.values():
        all_diseases |= ds

    out = {
        "description": "Cross-species reachability of the held-out test pairs on the R0 graph.",
        "train_file": TRAIN.name,
        "n_train_edges": n_edges,
        "n_test_pairs": len(test),
        "reachable_via_orthology_bridge": n_reachable,
        "reachable_via_orthology_bridge_frac": n_reachable / len(test),
        "reachable_via_nonhuman_ortholog": n_reachable_nonhuman,
        "of_those_also_reachable_without_orthology": n_also_plain,
        "reachable_only_without_orthology": n_plain_only,
        "distinct_diseases_reachable_from_nonhuman_genes": len(nonhuman_targets),
        "distinct_diseases_with_any_direct_neighbour": len(all_diseases),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
