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
    # genes -- model-organism databases (FlyBase, dictyBase, SGD, PomBase)
    "FB": "gene", "DICTYBASE": "gene", "SGD": "gene", "POMBASE": "gene",
    # diseases
    "MONDO": "disease", "DOID": "disease", "OMIM": "disease", "ORPHA": "disease",
    "ORPHANET": "disease", "SNOMEDCT": "disease", "UMLS": "disease",
    "MESH": "disease", "EFO": "disease",
    # phenotypes (incl. model-organism phenotype ontologies:
    # FYPO=fission yeast, WBPhenotype=worm, XPO=Xenopus, DDPHENO=Dictyostelium)
    "HP": "phenotype", "MP": "phenotype", "ZP": "phenotype",
    "FYPO": "phenotype", "WBPHENOTYPE": "phenotype", "XPO": "phenotype",
    "DDPHENO": "phenotype",
    # therapies
    "MAXO": "therapy",
    # proteins
    "UNIPROTKB": "protein", "UNIPROT": "protein", "PDB": "protein",
    "ENSP": "protein", "IPR": "protein", "PF": "protein",
    # compounds / drugs
    "CHEBI": "compound", "DRUGBANK": "compound", "CHEMBL": "compound",
    "PUBCHEM": "compound", "CID": "compound", "KEGG": "compound", "CAS": "compound",
    "RXCUI": "compound", "DRUGCENTRAL": "compound",
    # variants
    "DBSNP": "variant", "RS": "variant", "CLINVAR": "variant", "COSMIC": "variant",
    "HGVS": "variant", "CAID": "variant", "GNOMAD": "variant",
    # pathways
    "REACTOME": "pathway", "WP": "pathway", "WIKIPATHWAYS": "pathway",
    "GO": "pathway", "KEGG.PATHWAY": "pathway",
    # anatomy (incl. model-organism anatomy ontologies:
    # EMAPA=mouse, FBbt=fly, WBbt=worm, ZFA=zebrafish)
    "UBERON": "anatomy", "FMA": "anatomy", "CL": "anatomy",
    "EMAPA": "anatomy", "FBBT": "anatomy", "WBBT": "anatomy", "ZFA": "anatomy",
    # taxonomy
    "NCBITAXON": "taxonomy",
}


def category_of(node_id: str) -> str:
    """Return the category for a node ID, or 'unknown' if the prefix is unmapped."""
    if not node_id:
        return "unknown"
    prefix = node_id.split(":", 1)[0].upper()
    return PREFIX_CATEGORY_MAP.get(prefix, "unknown")
