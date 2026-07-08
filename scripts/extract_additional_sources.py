#!/usr/bin/env python3
"""
extract_additional_sources.py -- Ingest the non-Monarch databases into edges.

Background: the canonical graph (edges_clean.csv) was built ONLY from Monarch,
because the original extractor (step3-4.py) only recognized Monarch's KGX
subject/predicate/object schema. The other databases use bespoke schemas and were
silently skipped (see CRITICAL_FINDINGS.md).

This script adds the databases whose endpoints resolve to CLEAN CURIEs that
integrate with the existing HGNC/MONDO/ORPHA/CHEMBL namespaces. Quality-first:
we only emit an edge when BOTH endpoints are resolvable identifiers (free-text
disease/drug names are skipped rather than guessed, to avoid polluting the graph).

To make the bespoke gene identifiers (symbols / Ensembl / Entrez / UniProt) used
by these databases resolve to the graph's dominant HGNC namespace, we build a gene
crosswalk from the HGNC complete set (data/registry/hgnc_complete_set.txt) plus the
project's own entity_id_map. Diseases coded as DOID / EFO / Orphanet / OMIM are
mapped to MONDO via mondo_xrefs.tsv so they merge with existing MONDO nodes.

Integrated here:
  * Gene2Phenotype  HGNC gene  --GENE_ASSOCIATED_WITH_CONDITION-->  MONDO disease
  * DGIdb           HGNC/NCBIGene gene --INTERACTS_WITH--> CHEMBL/DRUGBANK/RXCUI drug
  * Orphadata(36)   HGNC gene  --CAUSES / GENE_ASSOCIATED_WITH_CONDITION--> MONDO/ORPHA
  * CiVIC           NCBIGene/HGNC gene --GENE_ASSOCIATED_WITH_CONDITION--> MONDO
                    (molecular_profile_id join: variants->entrez x evidence->DOID)
  * GWAS catalog    HGNC gene  --GENE_ASSOCIATED_WITH_CONDITION--> MONDO
                    (SNP_GENE_IDS Ensembl x MAPPED_TRAIT_URI EFO->MONDO)
  * DrugCentral     HGNC gene  --INTERACTS_WITH--> DrugCentral drug (STRUCT_ID)
  * DIDA            HGNC gene  --GENETICALLY_INTERACTS_WITH--> HGNC gene (digenic)

Genuinely NOT integrable here (documented, no clean edges available):
  DisGeNET   -- this dump is a variant annotation table (no disease column at all)
  LNCipedia  -- FASTA nucleotide sequences only, no relationships
  ClinGen / PharmGKB (ClinGPX*, ClinPGX) -- node/xref tables only (genes, drugs,
             phenotypes, variants); the relationship tables are not in this dump
  LncRNADisease -- ncRNA symbol + disease are BOTH free text; no offline MONDO
             label index is available to resolve disease names without guessing

Output: data/processed/additional_edges.csv in edges_clean format
        (source_id,relation,target_id,weight,dataset_sources)
"""
import csv
import re
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import pandas as pd

from config import RAW_DIR, PROCESSED_DIR, DATA_DIR

# csv field size: GWAS/CiVIC have very long free-text cells
csv.field_size_limit(10_000_000)

OUT_CSV = PROCESSED_DIR / "additional_edges.csv"
MONDO_XREFS = Path(__file__).resolve().parent / "mondo_xrefs.tsv"
HGNC_TSV = DATA_DIR / "registry" / "hgnc_complete_set.txt"
ENTITY_ID_MAP = PROCESSED_DIR / "entity_id_map.csv"


# --------------------------------------------------------------------------- #
# Gene crosswalk:  symbol / Ensembl / Entrez / UniProt  ->  HGNC:n
# --------------------------------------------------------------------------- #
class GeneXwalk:
    def __init__(self):
        self.sym = {}      # SYMBOL (upper) -> HGNC:n   (primary symbols only)
        self.alias = {}    # alias/prev SYMBOL (upper) -> HGNC:n (only if unambiguous)
        self.ens = {}      # ENSG... -> HGNC:n
        self.entrez = {}   # entrez str -> HGNC:n
        self.uniprot = {}  # UniProt acc -> HGNC:n

    def _add_alias(self, name, hgnc):
        if not name:
            return
        k = name.strip().upper()
        if not k:
            return
        if k in self.alias and self.alias[k] != hgnc:
            self.alias[k] = None       # ambiguous -> unusable
        elif k not in self.alias:
            self.alias[k] = hgnc

    def load_hgnc(self, path):
        if not path.exists():
            print(f"  WARNING: {path} missing; gene crosswalk will be weak")
            return
        df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
        cols = {c.lower(): c for c in df.columns}
        hc = cols["hgnc_id"]; sc = cols["symbol"]
        ac = cols.get("alias_symbol"); pc = cols.get("prev_symbol")
        ec = cols.get("entrez_id"); gc = cols.get("ensembl_gene_id")
        uc = cols.get("uniprot_ids")
        for _, r in df.iterrows():
            hgnc = r[hc]
            if not isinstance(hgnc, str) or not hgnc.startswith("HGNC:"):
                continue
            s = r[sc]
            if isinstance(s, str) and s.strip():
                self.sym[s.strip().upper()] = hgnc
            if ac and isinstance(r[ac], str):
                for a in r[ac].split("|"):
                    self._add_alias(a, hgnc)
            if pc and isinstance(r[pc], str):
                for a in r[pc].split("|"):
                    self._add_alias(a, hgnc)
            if ec and isinstance(r[ec], str) and r[ec].strip():
                self.entrez[r[ec].strip()] = hgnc
            if gc and isinstance(r[gc], str) and r[gc].strip():
                self.ens[r[gc].strip()] = hgnc
            if uc and isinstance(r[uc], str):
                for u in r[uc].split("|"):
                    u = u.strip()
                    if u:
                        self.uniprot.setdefault(u, hgnc)
        print(f"  HGNC crosswalk: {len(self.sym)} symbols, {len(self.entrez)} entrez, "
              f"{len(self.ens)} ensembl, {len(self.uniprot)} uniprot")

    def load_entity_map(self, path):
        """Add Entrez->HGNC links the project already established (broadens recall)."""
        if not path.exists():
            return
        n = 0
        for chunk in pd.read_csv(path, dtype=str, chunksize=200_000):
            sub = chunk[(chunk["canonical_id"].str.startswith("HGNC:", na=False)) &
                        (chunk["raw_id"].str.startswith("NCBIGene:", na=False))]
            for raw, can in zip(sub["raw_id"], sub["canonical_id"]):
                acc = raw.split(":", 1)[1]
                self.entrez.setdefault(acc, can)
                n += 1
        print(f"  +{n} Entrez->HGNC links from entity_id_map")

    def by_symbol(self, s):
        if not isinstance(s, str):
            return None
        k = s.strip().upper()
        if not k:
            return None
        return self.sym.get(k) or self.alias.get(k)

    def by_entrez(self, e):
        if not isinstance(e, str):
            return None
        e = e.strip()
        if e.endswith(".0"):
            e = e[:-2]
        return self.entrez.get(e) or (f"NCBIGene:{e}" if e.isdigit() else None)

    def by_ensembl(self, g):
        if not isinstance(g, str):
            return None
        return self.ens.get(g.strip().split(".")[0])

    def by_uniprot(self, u):
        if not isinstance(u, str):
            return None
        return self.uniprot.get(u.strip())


# --------------------------------------------------------------------------- #
# Disease crosswalk:  DOID / EFO / Orphanet / OMIM code  ->  MONDO:n
# --------------------------------------------------------------------------- #
def load_disease_maps():
    """Return {prefix_lower: {code: MONDO:n}} from mondo_xrefs.tsv."""
    maps = defaultdict(dict)
    if not MONDO_XREFS.exists():
        return maps
    df = pd.read_csv(MONDO_XREFS, sep="\t", dtype=str)
    for mondo, xref in zip(df["cls"], df["xref"]):
        if not isinstance(xref, str):
            continue
        x = xref.strip().strip('"')
        if ":" not in x:
            continue
        pre, code = x.split(":", 1)
        maps[pre.lower()].setdefault(code.strip(), mondo)
    print("  disease xref maps:", {k: len(v) for k, v in maps.items()
                                    if k in ("doid", "efo", "orphanet", "omim")})
    return maps


# --------------------------------------------------------------------------- #
# Source extractors -- each yields (subject, predicate, object, source)
# --------------------------------------------------------------------------- #
def extract_gene2phenotype(rows):
    """G2P CSV: 'hgnc id' (bare number or HGNC:n) + 'disease MONDO' (MONDO:n)."""
    out = []
    folders = sorted(RAW_DIR.glob("Gene2Phenotype*"))
    for fol in folders:
        for f in fol.glob("data.raw.csv*"):
            try:
                df = pd.read_csv(f, dtype=str, on_bad_lines="skip")
            except Exception as e:
                print(f"    {f}: {e}"); continue
            cols = {c.lower().strip(): c for c in df.columns}
            hc = cols.get("hgnc id"); dc = cols.get("disease mondo")
            if not hc or not dc:
                continue
            for h, d in zip(df[hc], df[dc]):
                if not isinstance(h, str) or not isinstance(d, str):
                    continue
                h = h.strip(); d = d.strip()
                if not h or not d or d.lower() == "nan":
                    continue
                hid = h if h.upper().startswith("HGNC:") else f"HGNC:{h}"
                if not re.match(r"^HGNC:\d+$", hid):
                    continue
                if not re.match(r"^MONDO:\d+$", d):
                    continue
                out.append((hid, "BIOLINK:GENE_ASSOCIATED_WITH_CONDITION", d, "Gene2Phenotype"))
    return out


DRUG_PREFIX = {"chembl": "CHEMBL", "drugbank": "DRUGBANK", "rxcui": "RXCUI"}
GENE_PREFIX = {"hgnc": "HGNC", "ncbigene": "NCBIGene"}


def _norm_concept(val, mapping):
    """'hgnc:2625' -> 'HGNC:2625' (only for prefixes we trust)."""
    if not isinstance(val, str):
        return None
    val = val.strip()
    if ":" not in val:
        return None
    pre, acc = val.split(":", 1)
    std = mapping.get(pre.lower())
    if not std or not acc:
        return None
    return f"{std}:{acc}"


def extract_dgidb(rows):
    """DGIdb interaction TSVs: gene_concept_id + drug_concept_id."""
    out = []
    for fol in sorted(RAW_DIR.glob("DGIDB*")):
        for f in fol.glob("data.raw.tsv*"):
            try:
                df = pd.read_csv(f, sep="\t", dtype=str, on_bad_lines="skip")
            except Exception as e:
                print(f"    {f}: {e}"); continue
            cols = {c.lower().strip(): c for c in df.columns}
            gc = cols.get("gene_concept_id"); dc = cols.get("drug_concept_id")
            if not gc or not dc:
                continue
            for g, d in zip(df[gc], df[dc]):
                gid = _norm_concept(g, GENE_PREFIX)
                did = _norm_concept(d, DRUG_PREFIX)
                if gid and did:
                    out.append((gid, "BIOLINK:INTERACTS_WITH", did, "DGIdb"))
    return out


ORPHA_CAUSAL = "Disease-causing"


def extract_orphadata(orpha2mondo):
    """Orphadata gene-association XML (Orphadata36): OrphaCode <-> Gene(HGNC xref)."""
    out = []
    targets = []
    for fol in sorted(RAW_DIR.glob("Orphadata*")):
        for f in fol.glob("data.raw.xml*"):
            try:
                head = f.open("r", encoding="utf-8", errors="ignore").read(20000)
            except Exception:
                continue
            if "DisorderGeneAssociation" in head:
                targets.append(f)
    for f in targets:
        try:
            tree = ET.parse(f); root = tree.getroot()
        except Exception as e:
            print(f"    {f}: {e}"); continue
        for disorder in root.iter("Disorder"):
            oc = disorder.findtext("OrphaCode")
            if not oc:
                continue
            dis = orpha2mondo.get(str(oc), f"ORPHA:{oc}")
            for assoc in disorder.iter("DisorderGeneAssociation"):
                gene = assoc.find("Gene")
                if gene is None:
                    continue
                hgnc = None
                for ref in gene.iter("ExternalReference"):
                    if (ref.findtext("Source") or "").strip().upper() == "HGNC":
                        r = (ref.findtext("Reference") or "").strip()
                        if r:
                            hgnc = r if r.upper().startswith("HGNC:") else f"HGNC:{r}"
                        break
                if not hgnc or not re.match(r"^HGNC:\d+$", hgnc):
                    continue
                atype = ""
                at = assoc.find("DisorderGeneAssociationType")
                if at is not None:
                    atype = (at.findtext("Name") or "")
                pred = ("BIOLINK:CAUSES" if ORPHA_CAUSAL in atype
                        else "BIOLINK:GENE_ASSOCIATED_WITH_CONDITION")
                out.append((hgnc, pred, dis, "Orphadata"))
    return out


def _find_civic(*name_substrs):
    """Locate a CiVIC dump whose header contains all given column substrings."""
    for fol in sorted(RAW_DIR.glob("CiVIC*")):
        for f in fol.glob("data.raw.tsv*"):
            try:
                hdr = f.open("r", encoding="utf-8", errors="ignore").readline().lower()
            except Exception:
                continue
            if all(s in hdr for s in name_substrs):
                return f
    return None


def extract_civic(xwalk, dmaps):
    """CiVIC: join variants (molecular_profile_id -> entrez gene) with evidence
    (molecular_profile_id -> DOID disease) to emit gene--disease associations."""
    out = []
    # 1) molecular_profile_id -> gene (entrez)  from the variants dump
    var_f = _find_civic("single_variant_molecular_profile_id", "entrez_id")
    mp2gene = {}
    if var_f:
        df = pd.read_csv(var_f, sep="\t", dtype=str, on_bad_lines="skip")
        cols = {c.lower(): c for c in df.columns}
        mc = cols.get("single_variant_molecular_profile_id"); ec = cols.get("entrez_id")
        if mc and ec:
            for mp, e in zip(df[mc], df[ec]):
                if isinstance(mp, str) and isinstance(e, str) and mp.strip() and e.strip():
                    g = xwalk.by_entrez(e)
                    if g:
                        mp2gene.setdefault(mp.strip(), g)
    # 2) molecular_profile_id -> disease (DOID) from the evidence dump
    ev_f = _find_civic("molecular_profile_id", "doid", "evidence_type")
    if ev_f and mp2gene:
        df = pd.read_csv(ev_f, sep="\t", dtype=str, on_bad_lines="skip")
        cols = {c.lower(): c for c in df.columns}
        mc = cols.get("molecular_profile_id"); dc = cols.get("doid")
        if mc and dc:
            doid_map = dmaps.get("doid", {})
            for mp, doid in zip(df[mc], df[dc]):
                if not isinstance(mp, str) or not isinstance(doid, str):
                    continue
                mp = mp.strip(); doid = doid.strip()
                gene = mp2gene.get(mp)
                if not gene or not doid or doid.lower() == "nan":
                    continue
                if doid.endswith(".0"):
                    doid = doid[:-2]
                mondo = doid_map.get(doid)
                if not mondo:
                    continue
                out.append((gene, "BIOLINK:GENE_ASSOCIATED_WITH_CONDITION", mondo, "CiVIC"))
    return out


def extract_gwas(xwalk, dmaps):
    """GWAS catalog: SNP_GENE_IDS (Ensembl) x MAPPED_TRAIT_URI (EFO) -> gene/MONDO.
    Streamed (the dump is ~650 MB)."""
    out = []
    f = RAW_DIR / "GWAS" / "data.raw.tsv.tsv"
    if not f.exists():
        return out
    efo_map = dmaps.get("efo", {})
    seen = set()
    with f.open("r", encoding="utf-8", errors="ignore", newline="") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        header = next(rdr)
        idx = {c.strip().upper(): i for i, c in enumerate(header)}
        gi = idx.get("SNP_GENE_IDS"); ti = idx.get("MAPPED_TRAIT_URI")
        if gi is None or ti is None:
            return out
        for row in rdr:
            if len(row) <= max(gi, ti):
                continue
            genes_raw = row[gi]; traits_raw = row[ti]
            if not genes_raw or not traits_raw:
                continue
            genes = [xwalk.by_ensembl(g) for g in re.split(r"[,;]", genes_raw)]
            genes = [g for g in genes if g]
            if not genes:
                continue
            mondos = []
            for uri in re.split(r"[,;]\s*", traits_raw):
                m = re.search(r"EFO_(\d+)", uri)
                if m:
                    mondo = efo_map.get(m.group(1))
                    if mondo:
                        mondos.append(mondo)
            for g in genes:
                for d in mondos:
                    key = (g, d)
                    if key not in seen:
                        seen.add(key)
                        out.append((g, "BIOLINK:GENE_ASSOCIATED_WITH_CONDITION", d, "GWAS"))
    return out


def extract_drugcentral(xwalk):
    """DrugCentral bioactivity: drug (STRUCT_ID) -- target gene (symbol/UniProt)."""
    out = []
    f = RAW_DIR / "DrugCentral" / "data.raw.tsv.tsv"
    if not f.exists():
        return out
    df = pd.read_csv(f, sep="\t", dtype=str, on_bad_lines="skip")
    cols = {c.lower().strip(): c for c in df.columns}
    sid = cols.get("struct_id"); gc = cols.get("gene")
    acc = cols.get("accession") or cols.get("swissprot")
    org = cols.get("organism")
    if not sid or (not gc and not acc):
        return out
    for i in range(len(df)):
        struct = df[sid].iat[i]
        if not isinstance(struct, str) or not struct.strip():
            continue
        if struct.endswith(".0"):
            struct = struct[:-2]
        # human targets only (keep the graph's human-gene focus clean)
        if org and isinstance(df[org].iat[i], str) and df[org].iat[i].strip() and \
           df[org].iat[i].strip().lower() != "homo sapiens":
            continue
        gene = None
        if gc:
            gene = xwalk.by_symbol(df[gc].iat[i])
        if not gene and acc:
            a = df[acc].iat[i]
            if isinstance(a, str):
                for tok in re.split(r"[|,;]", a):
                    gene = xwalk.by_uniprot(tok)
                    if gene:
                        break
        if not gene:
            continue
        drug = f"DrugCentral:{struct.strip()}"
        out.append((gene, "BIOLINK:INTERACTS_WITH", drug, "DrugCentral"))
    return out


def extract_dida(xwalk):
    """DIDA digenic combinations: Gene A -- Gene B (genetic interaction)."""
    out = []
    for fol in sorted(RAW_DIR.glob("DIDA*")):
        for f in fol.glob("data.raw.tsv*"):
            try:
                hdr = f.open("r", encoding="utf-8", errors="ignore").readline().lower()
            except Exception:
                continue
            if "gene a" not in hdr or "gene b" not in hdr:
                continue
            df = pd.read_csv(f, sep="\t", dtype=str, on_bad_lines="skip")
            cols = {c.lower().strip(): c for c in df.columns}
            ga = cols.get("gene a"); gb = cols.get("gene b")
            if not ga or not gb:
                continue
            for a, b in zip(df[ga], df[gb]):
                g1 = xwalk.by_symbol(a); g2 = xwalk.by_symbol(b)
                if g1 and g2 and g1 != g2:
                    out.append((g1, "BIOLINK:GENETICALLY_INTERACTS_WITH", g2, "DIDA"))
    return out


# --------------------------------------------------------------------------- #
def main():
    print("Loading crosswalks...")
    xwalk = GeneXwalk()
    xwalk.load_hgnc(HGNC_TSV)
    xwalk.load_entity_map(ENTITY_ID_MAP)
    dmaps = load_disease_maps()
    orpha2mondo = {code: mondo for code, mondo in dmaps.get("orphanet", {}).items()}

    all_edges = []
    def run(name, fn, *a):
        print(f"Extracting {name}...")
        e = fn(*a); print(f"  {len(e):,} raw edges"); all_edges.extend(e)

    run("Gene2Phenotype", extract_gene2phenotype, None)
    run("DGIdb", extract_dgidb, None)
    run("Orphadata", extract_orphadata, orpha2mondo)
    run("CiVIC", extract_civic, xwalk, dmaps)
    run("GWAS", extract_gwas, xwalk, dmaps)
    run("DrugCentral", extract_drugcentral, xwalk)
    run("DIDA", extract_dida, xwalk)

    # dedup within these sources, dropping self-loops; aggregate sources
    agg = {}
    for s, r, t, src in all_edges:
        if s == t:
            continue
        agg.setdefault((s, r, t), set()).add(src)

    with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["source_id", "relation", "target_id", "weight", "dataset_sources"])
        for (s, r, t), srcs in agg.items():
            w.writerow([s, r, t, 1, ";".join(sorted(srcs))])

    by_src = defaultdict(int)
    for _key, srcs in agg.items():
        for x in srcs:
            by_src[x] += 1
    print(f"\nWrote {len(agg):,} unique non-self-loop edges -> {OUT_CSV}")
    for k, v in sorted(by_src.items(), key=lambda x: -x[1]):
        print(f"  {k}: {v:,} edges")


if __name__ == "__main__":
    main()
