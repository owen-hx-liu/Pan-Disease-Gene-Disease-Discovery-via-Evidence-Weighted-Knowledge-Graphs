import csv
import sys
from pathlib import Path

INPUT_PATH = Path("data/processed/canonical_nodes.csv")
OUTPUT_PATH = Path("data/processed/nodes.csv")

if not INPUT_PATH.exists():
    sys.exit(f"❌ Missing input file: {INPUT_PATH}")

print("📂 Loading canonical nodes...")

# =====================================================
# PREFIX → CATEGORY MAP
# Derived from your Step 3.1 detection logic
# =====================================================

PREFIX_CATEGORY_MAP = {

    # ===== GENES =====
    "HGNC":"gene",
    "NCBIGENE":"gene",
    "ENTREZGENE":"gene",
    "ENSEMBL":"gene",
    "ENSG":"gene",
    "MGI":"gene",
    "RGD":"gene",
    "ZFIN":"gene",
    "FBGN":"gene",
    "VGNC":"gene",
    "WB":"gene",
    "XENBASE":"gene",
    # model-organism gene databases (FlyBase, dictyBase, SGD, PomBase)
    "FB":"gene",
    "DICTYBASE":"gene",
    "SGD":"gene",
    "POMBASE":"gene",

    # ===== DISEASES =====
    "MONDO":"disease",
    "DOID":"disease",
    "OMIM":"disease",
    "ORPHA":"disease",
    "ORPHANET":"disease",
    "SNOMEDCT":"disease",
    "UMLS":"disease",
    "MESH":"disease",
    "EFO":"disease",

    # ===== PHENOTYPES =====
    "HP":"phenotype",
    "MP":"phenotype",
    "ZP":"phenotype",
    # model-organism phenotype ontologies
    "FYPO":"phenotype",        # fission yeast
    "WBPHENOTYPE":"phenotype", # worm
    "XPO":"phenotype",         # Xenopus
    "DDPHENO":"phenotype",     # Dictyostelium

    # ===== THERAPIES =====
    "MAXO":"therapy",

    # ===== PROTEINS =====
    "UNIPROTKB":"protein",
    "UNIPROT":"protein",
    "PDB":"protein",
    "ENSP":"protein",
    "IPR":"protein",
    "PF":"protein",

    # ===== COMPOUNDS / DRUGS =====
    "CHEBI":"compound",
    "DRUGBANK":"compound",
    "CHEMBL":"compound",
    "PUBCHEM":"compound",
    "CID":"compound",
    "KEGG":"compound",
    "CAS":"compound",
    "RXCUI":"compound",
    "DRUGCENTRAL":"compound",

    # ===== VARIANTS =====
    "DBSNP":"variant",
    "RS":"variant",
    "CLINVAR":"variant",
    "COSMIC":"variant",
    "HGVS":"variant",
    "CAID":"variant",
    "GNOMAD":"variant",

    # ===== PATHWAYS =====
    "REACTOME":"pathway",
    "WP":"pathway",
    "WIKIPATHWAYS":"pathway",
    "GO":"pathway",
    "KEGG.PATHWAY":"pathway",

    # ===== ANATOMY =====
    "UBERON":"anatomy",
    "FMA":"anatomy",
    "CL":"anatomy",
    # model-organism anatomy ontologies
    "EMAPA":"anatomy", # mouse
    "FBBT":"anatomy",  # fly
    "WBBT":"anatomy",  # worm
    "ZFA":"anatomy",   # zebrafish

    # ===== TAXONOMY =====
    "NCBITAXON":"taxonomy"
}

nodes = []

with open(INPUT_PATH, newline='', encoding='utf-8') as f:
    reader = csv.DictReader(f)

    if "canonical_id" not in reader.fieldnames:
        sys.exit("❌ canonical_nodes.csv must contain 'canonical_id' column")

    for row in reader:
        cid = row["canonical_id"].strip()
        if not cid:
            continue

        prefix = cid.split(":",1)[0].upper()

        # Determine category
        category = PREFIX_CATEGORY_MAP.get(prefix,"unknown")

        nodes.append({
            "node_id": cid,
            "label": category,
            "prefix": prefix
        })

print(f"✓ Loaded {len(nodes):,} nodes")

print("📝 Writing Neo4j node table...")

OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)

with open(OUTPUT_PATH, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["node_id:ID","category:LABEL","prefix"])

    for n in nodes:
        writer.writerow([n["node_id"], n["label"], n["prefix"]])

print(f"✅ nodes.csv written to {OUTPUT_PATH}")
print("🎉 STEP 5.1 COMPLETE")
