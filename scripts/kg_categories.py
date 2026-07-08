"""Entity category lookup from ID prefix.

Single source of truth for mapping a canonical node ID (e.g. ``HGNC:1100``,
``MONDO:0007254``) to a biological category (gene, disease, phenotype, ...).
Mirrors the map in step5-1.py so type-aware code (type-matched negative sampling,
prediction filtering) stays consistent with the Neo4j node table.
"""
from __future__ import annotations

PREFIX_CATEGORY_MAP = {
    # genes
    "HGNC": "gene", "NCBIGENE": "gene", "ENTREZGENE": "gene", "ENSEMBL": "gene",
    "ENSG": "gene", "MGI": "gene", "RGD": "gene", "ZFIN": "gene", "FBGN": "gene",
    "VGNC": "gene", "WB": "gene", "XENBASE": "gene",
    # diseases
    "MONDO": "disease", "DOID": "disease", "OMIM": "disease", "ORPHA": "disease",
    "ORPHANET": "disease", "SNOMEDCT": "disease", "UMLS": "disease",
    "MESH": "disease", "EFO": "disease",
    # phenotypes
    "HP": "phenotype", "MP": "phenotype", "ZP": "phenotype",
    # therapies
    "MAXO": "therapy",
    # proteins
    "UNIPROTKB": "protein", "UNIPROT": "protein", "PDB": "protein",
    "ENSP": "protein", "IPR": "protein", "PF": "protein",
    # compounds / drugs
    "CHEBI": "compound", "DRUGBANK": "compound", "CHEMBL": "compound",
    "PUBCHEM": "compound", "CID": "compound", "KEGG": "compound", "CAS": "compound",
    "RXCUI": "compound",
    # variants
    "DBSNP": "variant", "RS": "variant", "CLINVAR": "variant", "COSMIC": "variant",
    "HGVS": "variant", "CAID": "variant", "GNOMAD": "variant",
    # pathways
    "REACTOME": "pathway", "WP": "pathway", "WIKIPATHWAYS": "pathway",
    "GO": "pathway", "KEGG.PATHWAY": "pathway",
    # anatomy
    "UBERON": "anatomy", "FMA": "anatomy", "CL": "anatomy",
    # taxonomy
    "NCBITAXON": "taxonomy",
}


def category_of(node_id: str) -> str:
    """Return the category for a node ID, or 'unknown' if the prefix is unmapped."""
    if not node_id:
        return "unknown"
    prefix = node_id.split(":", 1)[0].upper()
    return PREFIX_CATEGORY_MAP.get(prefix, "unknown")
