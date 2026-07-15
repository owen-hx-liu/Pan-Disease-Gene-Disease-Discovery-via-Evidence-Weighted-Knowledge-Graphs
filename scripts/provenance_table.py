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

# Upstream release version + license per source. Versions were recovered from the
# dumps themselves (embedded dates / preserved file timestamps) and cross-checked
# against the upstream download portals; licenses were confirmed on each project's
# terms page (see the reference URLs at the bottom of source_provenance.md).
# The download date is the local directory mtime (reported as `ingestion_date_utc`).
SOURCE_META = {
    "Monarch": {
        "release_version": "Monarch KG monthly release, downloaded within the "
                           "2025-12 / 2026-01-11 release window",
        "license": "BSD-3-Clause (Monarch KG & app; constituent primary sources "
                   "retain their own licenses)",
    },
    "DGIdb": {
        "release_version": "5.0.11 (upstream release dated 2026-01-13)",
        "license": "MIT (code/app); aggregated interaction data carries per-source "
                   "licenses, some redistribution-restricted (see dgidb.org sources)",
    },
    "GWAS": {
        "release_version": "NHGRI-EBI GWAS Catalog associations, v1.0.2 schema "
                           "(latest association added 2025-12-17)",
        "license": "EMBL-EBI Terms of Use (curated associations, freely available); "
                   "post-2021 summary statistics are CC0",
    },
    "DrugCentral": {
        "release_version": "2021 release, drug-target interaction export "
                           "(upstream file dated 2021-10-29)",
        "license": "CC BY-SA 4.0",
    },
    "Orphadata": {
        "release_version": "Orphadata product / Orphanet knowledge-base extract "
                           "(upstream file dated 2025-12-09)",
        "license": "CC BY 4.0",
    },
    "Gene2Phenotype": {
        "release_version": "EBI G2P panels (latest panel review 2025-07-23)",
        "license": "EMBL-EBI Terms of Use (freely available; EBI prefers CC0)",
    },
    "CiVIC": {
        "release_version": "CiVIC monthly data release, Jan-2026 build "
                           "(content through 2025-12-26)",
        "license": "CC0 1.0 (Public Domain Dedication)",
    },
    "DIDA": {
        "release_version": "v1 (2015; database static, no updates since ~2016)",
        "license": "CC BY 4.0 (per DIDA, Nucleic Acids Research 2016)",
    },
    # --- scanned but yield 0 clean edges from these dumps (not integrated) ---
    "DisGeNET": {
        "release_version": "variant-annotation table dump (NOT the standard "
                           "gene-disease-association release; ~v24 era)",
        "license": "CC BY-NC-SA 4.0 (academic / non-commercial)",
    },
    "LncRNADisease": {
        "release_version": "v3.0 (2023)",
        "license": "CC BY-NC (Attribution-NonCommercial)",
    },
    "ClinGen": {
        "release_version": "ClinGen + ClinPGx/PharmGKB node & xref tables "
                           "(no relationship table present in the dump)",
        "license": "ClinGen: CC0 1.0; PharmGKB/ClinPGx: CC BY-SA 4.0",
    },
    "LNCipedia": {
        "release_version": "5.2 (2019; current release)",
        "license": "Academic / non-commercial only (commercial use requires "
                   "written approval from Ghent University)",
    },
}

# Reference URLs backing the version/license values above (recorded for the paper's
# data-availability statement).
LICENSE_REFS = [
    ("Monarch Initiative", "https://monarchinitiative.org/  (KG BSD-3-Clause)"),
    ("DGIdb", "https://www.dgidb.org/downloads  (v5.0.11; per-source licenses)"),
    ("GWAS Catalog", "https://www.ebi.ac.uk/gwas/docs/about  (EMBL-EBI Terms of Use)"),
    ("DrugCentral", "https://drugcentral.org/  (CC BY-SA 4.0)"),
    ("Orphadata", "https://www.orphadata.com/legal-notice/  (CC BY 4.0)"),
    ("Gene2Phenotype", "https://www.ebi.ac.uk/gene2phenotype/  (EMBL-EBI Terms of Use)"),
    ("CiVIC", "https://docs.civicdb.org/en/latest/about.html  (CC0 1.0)"),
    ("DIDA", "https://doi.org/10.1093/nar/gkv1068  (CC BY 4.0)"),
    ("DisGeNET", "https://www.disgenet.org/  (CC BY-NC-SA 4.0)"),
    ("LncRNADisease", "http://www.rnanut.net/lncrnadisease/  (CC BY-NC)"),
    ("ClinGen / PharmGKB", "https://www.pharmgkb.org/  (CC0 / CC BY-SA 4.0)"),
    ("LNCipedia", "https://lncipedia.org/download  (academic / non-commercial)"),
]


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
        meta = SOURCE_META.get(db, {})
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
            "source_release_version": meta.get("release_version",
                                               "unknown (not embedded in dump)"),
            "license": meta.get("license", "unknown (verify upstream)"),
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
    # upstream release version / download date / license (data-availability table)
    md += ["\n## Upstream release version, download date & license\n",
           "Download date = local ingestion (directory mtime). Versions were "
           "recovered from the dumps (embedded dates / preserved file timestamps) "
           "and cross-checked against the upstream portals; licenses were confirmed "
           "on each project's terms page (URLs below). Sources with 0 edges were "
           "scanned but not integrated from these dumps.\n",
           "| Source | Upstream release version | Download date (UTC) | License |",
           "|---|---|---|---|"]
    for r in rows:
        md.append(f"| {r['source']} | {r['source_release_version']} | "
                  f"{r['ingestion_date_utc']} | {r['license']} |")

    md += ["\n### Reference URLs\n"]
    for name, url in LICENSE_REFS:
        md.append(f"- **{name}**: {url}")

    md += ["\n## Non-independence (read before weighting consensus)\n"]
    for db, note in AGGREGATORS.items():
        md.append(f"- **{db}**: {note}")
    (OUT_DIR / "source_provenance.md").write_text("\n".join(md) + "\n", encoding="utf-8")

    print(f"wrote {csv_path}")
    print(f"wrote {OUT_DIR / 'source_provenance.md'}")
    print("\n".join(md[2:]))


if __name__ == "__main__":
    main()
