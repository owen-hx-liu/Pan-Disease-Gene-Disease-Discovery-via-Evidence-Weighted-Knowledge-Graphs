#!/usr/bin/env python3
"""
qc_graph.py -- Data-quality audit of the integrated edge list.

Checks for structural problems that would silently harm the benchmark or get caught
in review:
  1. CSV field-count integrity (embedded/unquoted delimiters -> row corruption)
  2. Null / empty source, relation, target
  3. Self-loops (source == target)
  4. Exact duplicate edges (source, relation, target)
  5. Reverse-duplicate edges for symmetric relations (train/test LEAKAGE risk)
  6. Malformed IDs: duplicated prefixes (HGNC:HGNC:123), missing prefix (no ':')
  7. Relation prefix consistency (all BIOLINK:?)
  8. Type sanity on the target relation (gene -> disease)
  9. weight / dataset_sources population

Writes data/processed/qc_graph_report.json and prints a summary.
Usage: python scripts/qc_graph.py [--edges PATH]
"""
import argparse, json, re, time
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from kg_categories import category_of
except ImportError:
    from scripts.kg_categories import category_of

SYMMETRIC = {"INTERACTS_WITH", "ORTHOLOGOUS_TO", "GENETICALLY_INTERACTS_WITH",
             "GENETICALLY_ASSOCIATED_WITH", "RELATED_TO"}


def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def suffix(r): return r.split(":", 1)[1] if ":" in r else r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--edges", default="data/processed/edges_clean_integrated.csv")
    args = ap.parse_args()
    report = {}

    # -- 1. raw field-count integrity (before pandas can hide misalignment) ---- #
    log("scanning raw lines for field-count integrity...")
    total_lines = 0
    bad_fieldcount = 0
    bad_examples = []
    with open(args.edges, "r", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n")
        n_cols = header.count(",") + 1
        for line in f:
            total_lines += 1
            c = line.count(",")
            if c != n_cols - 1:
                bad_fieldcount += 1
                if len(bad_examples) < 5:
                    bad_examples.append(line.rstrip("\n")[:160])
    report["header"] = header
    report["expected_cols"] = n_cols
    report["data_rows_raw"] = total_lines
    report["rows_wrong_fieldcount"] = bad_fieldcount
    report["wrong_fieldcount_examples"] = bad_examples
    log(f"  {total_lines:,} data rows; {bad_fieldcount:,} rows with != {n_cols} fields")

    # -- load with pandas -------------------------------------------------- #
    log("loading with pandas...")
    df = pd.read_csv(args.edges, dtype=str, keep_default_na=False, na_values=[])
    report["data_rows_pandas"] = len(df)
    cols = list(df.columns)
    report["columns"] = cols
    s, r, t = "source_id", "relation", "target_id"

    # -- 2. null / empty --------------------------------------------------- #
    empt = {c: int((df[c].str.strip() == "").sum()) for c in [s, r, t]}
    report["empty_fields"] = empt
    log(f"  empty core fields: {empt}")

    # -- 3. self-loops ----------------------------------------------------- #
    selfloops = int((df[s] == df[t]).sum())
    report["self_loops"] = selfloops
    log(f"  self-loops: {selfloops:,}")

    # -- 4. exact duplicate edges ------------------------------------------ #
    dup = int(df.duplicated(subset=[s, r, t]).sum())
    report["exact_duplicate_edges"] = dup
    log(f"  exact duplicate (s,r,t) rows: {dup:,}")

    # -- 5. reverse-duplicate edges on symmetric relations ----------------- #
    rel_suf = df[r].map(suffix)
    sym_mask = rel_suf.isin(SYMMETRIC)
    sym = df[sym_mask]
    # undirected key = sorted(src,tgt)+rel ; count directed edges whose reverse also exists
    a = np.minimum(sym[s].values, sym[t].values)
    b = np.maximum(sym[s].values, sym[t].values)
    keys = pd.Series([f"{x}\t{y}\t{z}" for x, y, z in zip(a, b, rel_suf[sym_mask].values)])
    vc = keys.value_counts()
    both_dir = int((vc >= 2).sum())
    report["symmetric_edges_total"] = int(sym_mask.sum())
    report["symmetric_pairs_with_both_directions"] = both_dir
    report["symmetric_reverse_dupe_directed_edges"] = int((vc[vc >= 2]).sum())
    log(f"  symmetric relations: {int(sym_mask.sum()):,} edges; "
        f"{both_dir:,} undirected pairs present in BOTH directions "
        f"(leakage-relevant)")

    # -- 6. malformed IDs -------------------------------------------------- #
    def dup_prefix(series):
        # 'HGNC:HGNC:123' -> prefix repeated
        parts = series.str.split(":", n=2)
        return parts.map(lambda p: len(p) >= 3 and p[0] == p[1])
    dp_s = int(dup_prefix(df[s]).sum()); dp_t = int(dup_prefix(df[t]).sum())
    noprefix_s = int((~df[s].str.contains(":", regex=False)).sum())
    noprefix_t = int((~df[t].str.contains(":", regex=False)).sum())
    report["duplicated_prefix_ids"] = {"source": dp_s, "target": dp_t}
    report["ids_without_prefix"] = {"source": noprefix_s, "target": noprefix_t}
    log(f"  duplicated-prefix IDs: source={dp_s:,} target={dp_t:,}")
    log(f"  IDs without ':' prefix: source={noprefix_s:,} target={noprefix_t:,}")

    # -- 7. relation prefix consistency ------------------------------------ #
    non_biolink = int((~df[r].str.startswith("BIOLINK:")).sum())
    report["relations_not_biolink_prefixed"] = non_biolink
    report["n_relations"] = int(df[r].nunique())
    log(f"  relations: {df[r].nunique()} distinct; {non_biolink:,} not BIOLINK:-prefixed")

    # -- 8. type sanity on target relation --------------------------------- #
    gac = df[rel_suf == "GENE_ASSOCIATED_WITH_CONDITION"]
    src_cat = gac[s].map(category_of).value_counts().to_dict()
    tgt_cat = gac[t].map(category_of).value_counts().to_dict()
    report["GENE_ASSOCIATED_WITH_CONDITION"] = {
        "n": int(len(gac)),
        "source_category": src_cat,
        "target_category": tgt_cat,
    }
    log(f"  GENE_ASSOCIATED_WITH_CONDITION: {len(gac):,} edges; "
        f"source cats={src_cat}; target cats={tgt_cat}")

    # -- 9. weight / dataset_sources --------------------------------------- #
    if "weight" in cols:
        w = pd.to_numeric(df["weight"], errors="coerce")
        report["weight"] = {"null_or_nonnumeric": int(w.isna().sum()),
                            "min": None if w.isna().all() else float(w.min()),
                            "max": None if w.isna().all() else float(w.max()),
                            "n_distinct": int(w.nunique())}
    if "dataset_sources" in cols:
        report["dataset_sources_empty"] = int((df["dataset_sources"].str.strip() == "").sum())

    # -- verdict ----------------------------------------------------------- #
    blockers = []
    if bad_fieldcount: blockers.append(f"{bad_fieldcount} rows with wrong field count")
    if empt[s] or empt[r] or empt[t]: blockers.append("empty core fields present")
    if dup: blockers.append(f"{dup} exact duplicate edges")
    if dp_s or dp_t: blockers.append("duplicated-prefix IDs present")
    if noprefix_s or noprefix_t: blockers.append("IDs without prefix present")
    report["blockers"] = blockers

    out = Path("data/processed/qc_graph_report.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2))
    log(f"wrote {out}")
    print("\n===== QC VERDICT =====")
    if not blockers:
        print("No hard blockers found. Review the warnings below before finalizing.")
    else:
        print("ISSUES TO ADDRESS:")
        for b in blockers:
            print(f"  - {b}")
    print(f"self-loops={selfloops:,}  exact-dupes={dup:,}  "
          f"symmetric-both-dir-pairs={both_dir:,}")


if __name__ == "__main__":
    main()
