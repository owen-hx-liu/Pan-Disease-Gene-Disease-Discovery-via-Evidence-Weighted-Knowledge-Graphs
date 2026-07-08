#!/usr/bin/env python3
"""
provenance_table.py -- Build a data-provenance / source-contribution table.

Addresses the Tier-2 reproducibility requirement (CLAUDE.md, PROJECT_REPORT.md
6.7/7.10): record per-source counts, formats, raw sizes, and -- importantly --
the SOURCE NON-INDEPENDENCE that makes "supported by N datasets" not equal to N
independent votes (Monarch aggregates several of the other sources).

Outputs:
  data/processed/provenance/source_provenance.csv
  data/processed/provenance/source_provenance.md

Note on dates: the raw dumps carry no machine-readable source-release version, so
we record the local file modification time as the INGESTION date and flag that the
upstream release version/date must be filled in manually before submission.
"""
import os, re, json, csv, sys, argparse
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd

from config import EDGES_CLEAN_CSV, RAW_DIR, PROCESSED_DIR

OUT_DIR = PROCESSED_DIR / "provenance"

# Which sources are aggregators (re-distribute other primary sources) -> their
# edges should NOT be counted as independent of the sources they ingest.
AGGREGATORS = {
    "Monarch": "Aggregates HPO, MGI, ZFIN, ClinGen, OMIM, Orphanet, etc. "
               "Counts overlap heavily with primary sources.",
    "DisGeNET": "Aggregates curated + text-mined gene-disease sources "
                "(overlaps ClinGen, Orphanet, CTD).",
    "DGIDB": "Aggregates ~30 drug-gene sources.",
}

FORMAT_HINTS = {  # default format per DB (from data/raw inspection)
    "Monarch": "TSV", "Orphadata": "JSON/XML", "CiVIC": "TSV", "DGIdb": "TSV",
    "Gene2Phenotype": "CSV", "DIDA": "CSV", "ClinGen": "TSV", "ClinPGX": "TSV",
    "DisGeNET": "TSV", "DrugCentral": "TSV", "GWAS": "TSV",
    "LncRNADisease": "TSV", "LNCipedia": "FASTA",
}


# canonicalize source-name variants so edge tokens and raw-folder names merge
NAME_CANON = {"DGIDB": "DGIdb", "DGIDB2": "DGIdb", "CLINGPX": "ClinGen",
              "CLINPGX": "ClinGen", "GENE2PHENOTYPE": "Gene2Phenotype"}


def db_of(token: str) -> str:
    """Strip a trailing dump index and canonicalize:
    'Monarch16' -> 'Monarch', 'DGIDB2' -> 'DGIdb'."""
    base = re.sub(r"\d+$", "", token.strip())
    return NAME_CANON.get(base.upper(), base)


def dir_size_mb(path: Path) -> float:
    total = 0
    for root, _dirs, files in os.walk(path):
        for f in files:
            try:
                total += os.path.getsize(os.path.join(root, f))
            except OSError:
                pass
    return round(total / 1e6, 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--edges", default=str(EDGES_CLEAN_CSV),
                    help="edge CSV with a dataset_sources column")
    args = ap.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # 1) edges contributed per DB (an edge counts toward each DB that supports it)
    edge_counts = defaultdict(int)
    sole_source_counts = defaultdict(int)  # edges where this DB is the ONLY source
    total_edges = 0
    for chunk in pd.read_csv(args.edges, usecols=["dataset_sources"],
                             dtype=str, chunksize=500_000):
        for src in chunk["dataset_sources"].dropna():
            total_edges += 1
            dbs = {db_of(t) for t in str(src).split(";") if t.strip()}
            for db in dbs:
                edge_counts[db] += 1
            if len(dbs) == 1:
                sole_source_counts[next(iter(dbs))] += 1

    # 2) raw-folder metadata per DB
    raw_meta = defaultdict(lambda: {"dumps": 0, "size_mb": 0.0, "mtime": None})
    if RAW_DIR.exists():
        for d in sorted(RAW_DIR.iterdir()):
            if not d.is_dir():
                continue
            db = db_of(d.name)
            raw_meta[db]["dumps"] += 1
            raw_meta[db]["size_mb"] += dir_size_mb(d)
            mt = datetime.fromtimestamp(d.stat().st_mtime, tz=timezone.utc).date().isoformat()
            cur = raw_meta[db]["mtime"]
            raw_meta[db]["mtime"] = mt if cur is None else max(cur, mt)

    dbs = sorted(set(edge_counts) | set(raw_meta), key=lambda x: -edge_counts.get(x, 0))
    rows = []
    for db in dbs:
        m = raw_meta.get(db, {})
        rows.append({
            "source": db,
            "dumps": m.get("dumps", 0),
            "format": FORMAT_HINTS.get(db, "?"),
            "raw_size_mb": round(m.get("size_mb", 0.0), 1),
            "ingestion_date_utc": m.get("mtime") or "",
            "edges_supported": edge_counts.get(db, 0),
            "edges_sole_source": sole_source_counts.get(db, 0),
            "is_aggregator": db in AGGREGATORS,
            "independence_note": AGGREGATORS.get(db, "primary source"),
            "source_release_version": "TODO-fill-before-submission",
        })

    # write CSV
    csv_path = OUT_DIR / "source_provenance.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    # write markdown
    md = ["# Source provenance & contribution\n",
          f"Total cleaned edges: **{total_edges:,}**. "
          "`edges_supported` = edges citing this source (an edge may cite several); "
          "`edges_sole_source` = edges where it is the ONLY source.\n",
          "| Source | Dumps | Format | Raw MB | Ingested | Edges supported | "
          "Edges sole-source | Aggregator? |",
          "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        md.append(f"| {r['source']} | {r['dumps']} | {r['format']} | {r['raw_size_mb']} "
                  f"| {r['ingestion_date_utc']} | {r['edges_supported']:,} | "
                  f"{r['edges_sole_source']:,} | {'yes' if r['is_aggregator'] else 'no'} |")
    md += ["\n## Non-independence (read before weighting consensus)\n"]
    for db, note in AGGREGATORS.items():
        md.append(f"- **{db}**: {note}")
    md.append("\n> `source_release_version` is a placeholder: upstream release "
              "versions/dates are not embedded in the dumps and must be recorded "
              "manually before submission (journals increasingly require this).")
    (OUT_DIR / "source_provenance.md").write_text("\n".join(md) + "\n", encoding="utf-8")

    print(f"wrote {csv_path}")
    print(f"wrote {OUT_DIR / 'source_provenance.md'}")
    print("\n".join(md[2:]))


if __name__ == "__main__":
    main()
