# ===================== 3.1 — ENTITY EXTRACTION AND RELATIONSHIP MAPPING =====================

import pandas as pd
import json
import re
from pathlib import Path
from tqdm import tqdm
import xml.etree.ElementTree as ET
from Bio import SeqIO

RAW_DATA_DIR = Path("data/raw")
MONDO_XREF_FILE = Path("scripts/mondo_xrefs.tsv")
OUTPUT_ENTITIES = Path("data/processed/entity_vocab.csv")
OUTPUT_REL_MAP = Path("data/processed/dataset_relation_map.json")

# ===================== EXTENDED PREFIX MAP =====================

ENTITY_PREFIX_MAP = {

    # ===== GENES =====
    "HGNC:":"Gene","NCBIGENE:":"Gene","ENSEMBL:":"Gene",
    "ENSG":"Gene","ENSMUSG":"Gene","ENSDARG":"Gene",

    "MGI:":"Gene","RGD:":"Gene","ZFIN:":"Gene","SGD:":"Gene",
    "WB:":"Gene","WORMBASE:":"Gene",
    "FB:":"Gene","FBGN:":"Gene",
    "POMBASE:":"Gene","DICTYBASE:":"Gene",
    "XENBASE:":"Gene",

    # ===== DISEASE =====
    "MONDO:":"Disease","DOID:":"Disease","OMIM:":"Disease",
    "ORPHANET:":"Disease","MESH:":"Disease","UMLS:":"Disease",

    # ===== PHENOTYPE =====
    "HP:":"Phenotype","MP:":"Phenotype",
    "WBPHENOTYPE:":"Phenotype",

    # ===== CHEMICAL =====
    "CHEBI:":"Chemical","DRUGBANK:":"Chemical","PUBCHEM:":"Chemical",

    # ===== VARIANT =====
    "CLINVAR:":"Variant","CAID:":"Variant",
    "DBSNP:":"Variant","RS":"Variant",

    # ===== PATHWAY =====
    "GO:":"Pathway","REACT:":"Pathway","KEGG:":"Pathway","REACTOME:":"Pathway",

    # ===== ANATOMY =====
    "UBERON:":"Anatomy","ZFA:":"Anatomy",
    "FBBT:":"Anatomy","WBBT:":"Anatomy",

    # ===== THERAPY =====
    "MAXO:":"Therapy"
}

# ===================== TYPE INFERENCE =====================

def infer_entity_type(eid):

    if not isinstance(eid,str):
        return "Unknown"

    eid_up = eid.upper()

    # ---------- ORIGINAL PREFIX LOGIC ----------
    for prefix,etype in ENTITY_PREFIX_MAP.items():
        if eid_up.startswith(prefix):
            return etype

    # ---------- Namespace inference ----------
    if ":" in eid_up:
        ns = eid_up.split(":",1)[0]

        if ns in {
            "HGNC","NCBIGENE","ENSEMBL","ENSG","ENSMUSG","ENSDARG",
            "MGI","RGD","ZFIN","SGD","WB","FB","POMBASE",
            "DICTYBASE","XENBASE"
        }:
            return "Gene"

        if ns in {"MONDO","DOID","OMIM","ORPHANET"}:
            return "Disease"

        if ns in {"HP","MP","WBPHENOTYPE"}:
            return "Phenotype"

        if ns in {"CHEBI","PUBCHEM","DRUGBANK"}:
            return "Chemical"

        if ns in {"CLINVAR","CAID","DBSNP"}:
            return "Variant"

        if ns in {"GO","REACTOME","KEGG"}:
            return "Pathway"

        if ns in {"UBERON","ZFA","FBBT","WBBT"}:
            return "Anatomy"

        if ns in {"MAXO"}:
            return "Therapy"

    # ---------- STRONG HEURISTICS ----------

    # lncRNAs
    if re.match(r"^LINC\d+",eid_up):
        return "Gene"

    # antisense genes
    if re.search(r"-AS\d*$",eid_up):
        return "Gene"

    # ORFs
    if re.search(r"ORF\d+",eid_up):
        return "Gene"

    # Ensembl IDs
    if re.match(r"^ENS[A-Z]*G\d+",eid_up):
        return "Gene"

    # rs variants
    if re.match(r"^RS\d+",eid_up):
        return "Variant"

    # Symbol:number
    if re.match(r"^[A-Z0-9\-]+:\d+$",eid_up):
        return "Gene"

    return "Unknown"

# ===================== TEXT NORMALIZATION =====================

def normalize_text(x):
    if not isinstance(x,str):
        return ""
    x=x.lower()
    x=re.sub(r"[^a-z0-9\s]"," ",x)
    return re.sub(r"\s+"," ",x).strip()

# ===================== MONDO XREFS =====================

print("Loading MONDO xrefs...")
xref_df=pd.read_csv(MONDO_XREF_FILE,sep="\t",names=["mondo_id","xref"])
xref_df["xref"]=xref_df["xref"].str.strip('"')
XREF_MAP=dict(zip(xref_df["xref"],xref_df["mondo_id"]))

def normalize_disease_id(did):
    if not isinstance(did,str):
        return ""
    if did.startswith("MONDO:"):
        return did
    return XREF_MAP.get(did,did)

# ===================== LOADERS (UNCHANGED) =====================

def load_csv_tsv(path):
    try:
        return pd.read_csv(path,sep="\t" if path.suffix==".tsv" else ",",
                           low_memory=False,on_bad_lines="skip")
    except:
        return pd.DataFrame()

def load_json(path):
    try:
        return pd.DataFrame(json.load(open(path,encoding="utf-8")))
    except:
        return pd.DataFrame()

def load_xml(path):
    rows=[]
    try:
        tree=ET.parse(path)
        root=tree.getroot()
        for rec in root.findall(".//record"):
            rows.append({
                "subject_id":rec.findtext("subject"),
                "subject_name":rec.findtext("subject_label"),
                "object_id":rec.findtext("object"),
                "object_name":rec.findtext("object_label")
            })
    except:
        pass
    return pd.DataFrame(rows)

def load_fasta(path):
    rows=[]
    try:
        for rec in SeqIO.parse(path,"fasta"):
            rows.append({"subject_id":rec.id,"subject_name":rec.description})
    except:
        pass
    return pd.DataFrame(rows)

# ===================== MAIN INGEST (UNCHANGED) =====================

all_entities=[]
dataset_relation_map={}

files=list(RAW_DATA_DIR.rglob("data.raw.*"))
print(f"Found {len(files)} raw datasets\n")

for path in tqdm(files,desc="Processing datasets"):

    df=pd.DataFrame()

    if path.suffix in [".tsv",".csv"]:
        df=load_csv_tsv(path)
    elif path.suffix==".json":
        df=load_json(path)
    elif path.suffix==".xml":
        df=load_xml(path)
    elif path.suffix in [".fasta",".fa"]:
        df=load_fasta(path)
    else:
        continue

    if df.empty:
        continue

    cols={c.lower():c for c in df.columns}

    subj_id=cols.get("subject") or cols.get("subject_id")
    subj_name=cols.get("subject_label") or cols.get("subject_name")
    obj_id=cols.get("object") or cols.get("object_id")
    obj_name=cols.get("object_label") or cols.get("object_name")

    rel_cols=[c for c in df.columns if c not in {subj_id,subj_name,obj_id,obj_name}]
    if rel_cols:
        dataset_relation_map[path.name]=rel_cols

    for _,row in df.iterrows():

        if subj_id:
            sid=row.get(subj_id)
            all_entities.append({
                "id":sid,
                "name":normalize_text(row.get(subj_name)),
                "type":infer_entity_type(sid)
            })

        if obj_id:
            oid=normalize_disease_id(row.get(obj_id))
            all_entities.append({
                "id":oid,
                "name":normalize_text(row.get(obj_name)),
                "type":infer_entity_type(oid)
            })

# ===================== SAVE =====================

entity_df=pd.DataFrame(all_entities)
entity_df=entity_df.dropna(subset=["id"])
entity_df=entity_df.drop_duplicates("id")
entity_df=entity_df.sort_values("type")

OUTPUT_ENTITIES.parent.mkdir(parents=True,exist_ok=True)
entity_df.to_csv(OUTPUT_ENTITIES,index=False)

with open(OUTPUT_REL_MAP,"w") as f:
    json.dump(dataset_relation_map,f,indent=2)

print("\n✅ 3.1 COMPLETE")
print(entity_df["type"].value_counts())
