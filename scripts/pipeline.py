"""
Phase 3: Node and Edge Extraction Pipeline (Steps 3.3 - 3.5)
Supports CSV, TSV, JSON, XML, FASTA datasets.
Memory-safe streaming, deduplication, automatic entity and relationship labeling.
"""

import os
import pandas as pd
from tqdm import tqdm
import networkx as nx
import json
import xml.etree.ElementTree as ET

# Optional Biopython for FASTA
try:
    from Bio import SeqIO
    BIOPYTHON_AVAILABLE = True
except ModuleNotFoundError:
    print("Warning: Biopython not installed. FASTA files will be skipped.")
    BIOPYTHON_AVAILABLE = False

# ----------------------
# CONFIG
# ----------------------
RAW_DATA_DIR = "data/raw"
PROCESSED_DIR = "data/processed"

ENTITY_ID_KEYS = ["GENE", "GENE_ID", "ENTREZ_ID", "DISEASE", "DISEASE_ID",
                  "VARIANT", "CHEMICAL", "PATHWAY", "PHENOTYPE"]

# Mapping from column keywords → standard entity types
ENTITY_TYPE_MAP = {
    "GENE": "Gene",
    "GENE_ID": "Gene",
    "ENTREZ_ID": "Gene",
    "DISEASE": "Disease",
    "DISEASE_ID": "Disease",
    "VARIANT": "Variant",
    "CHEMICAL": "Chemical",
    "PATHWAY": "Pathway",
    "PHENOTYPE": "Phenotype"
}

# Standard relationships (Step 3.2)
RELATIONSHIP_MAP = {
    ("Gene", "Gene"): "INTERACTS_WITH",
    ("Gene", "Disease"): "ASSOCIATED_WITH",
    ("Chemical", "Disease"): "TREATS",
    ("Gene", "Phenotype"): "EXPRESSED_IN",
    ("Gene", "Variant"): "CAUSES",
    # Add more heuristics if needed
}

# ----------------------
# HELPER FUNCTIONS
# ----------------------

def ensure_processed_dir():
    if not os.path.exists(PROCESSED_DIR):
        os.makedirs(PROCESSED_DIR)

def safe_read_table(file_path, nrows=None):
    """Read CSV or TSV safely, skip bad lines"""
    try:
        if file_path.lower().endswith(".tsv"):
            return pd.read_csv(file_path, sep="\t", low_memory=False,
                               nrows=nrows, on_bad_lines='skip')
        elif file_path.lower().endswith(".csv"):
            return pd.read_csv(file_path, low_memory=False,
                               nrows=nrows, on_bad_lines='skip')
        else:
            return None
    except Exception as e:
        print(f"Failed to read {file_path}: {e}")
        return None

def extract_candidate_columns(df):
    """Return columns matching known entity IDs"""
    return [c for c in df.columns if any(k in c.upper() for k in ENTITY_ID_KEYS)]

def infer_entity_type(col_name):
    """Map column name to standard entity type"""
    for key, etype in ENTITY_TYPE_MAP.items():
        if key in col_name.upper():
            return etype
    return "Unknown"

def assign_relation(source_type, target_type):
    """Assign standard relationship based on entity types"""
    return RELATIONSHIP_MAP.get((source_type, target_type), "ASSOCIATED_WITH")

# ----------------------
# NODE EXTRACTION
# ----------------------

def extract_nodes_from_dataframe(df):
    """Extract nodes from a DataFrame"""
    nodes = {}
    if df is None:
        return nodes

    id_columns = extract_candidate_columns(df)
    for col in id_columns:
        entity_type = infer_entity_type(col)
        for _, row in df.iterrows():
            raw_id = row.get(col)
            if pd.isna(raw_id):
                continue
            raw_id = str(raw_id).strip()
            if raw_id not in nodes:
                nodes[raw_id] = {
                    "entity_type": entity_type,
                    "raw_id": raw_id,
                    "name": row.get(col + "_NAME") if col + "_NAME" in df.columns else "",
                    "alt_ids": ""
                }
            # Add alternative IDs (other columns in the same row)
            alt_ids = []
            for other_col in id_columns:
                if other_col == col:
                    continue
                other_val = row.get(other_col)
                if pd.notna(other_val):
                    alt_ids.append(str(other_val).strip())
            if alt_ids:
                nodes[raw_id]["alt_ids"] = "|".join(set(alt_ids))
    return nodes

def extract_nodes_from_fasta(file_path, max_nodes=1000):
    """Extract nodes from FASTA headers"""
    nodes = {}
    if not BIOPYTHON_AVAILABLE:
        return nodes
    try:
        for i, record in enumerate(SeqIO.parse(file_path, "fasta")):
            if i >= max_nodes:
                break
            raw_id = record.id.strip()
            if raw_id not in nodes:
                nodes[raw_id] = {
                    "entity_type": "Gene",  # default for FASTA
                    "raw_id": raw_id,
                    "name": record.description.strip(),
                    "alt_ids": ""
                }
    except Exception as e:
        print(f"Failed to parse FASTA {file_path}: {e}")
    return nodes

def extract_nodes_from_json(file_path, max_rows=1000):
    """Extract nodes from JSON file"""
    nodes = {}
    try:
        with open(file_path) as f:
            try:
                data = json.load(f)
                if isinstance(data, dict):
                    data = [data]
            except json.JSONDecodeError:
                f.seek(0)
                data = []
                for i, line in enumerate(f):
                    data.append(json.loads(line))
                    if i+1 >= max_rows:
                        break

        for obj in data:
            for col, val in obj.items():
                if any(k in col.upper() for k in ENTITY_ID_KEYS):
                    raw_id = str(val).strip()
                    if raw_id not in nodes:
                        nodes[raw_id] = {
                            "entity_type": infer_entity_type(col),
                            "raw_id": raw_id,
                            "name": obj.get(col + "_NAME", ""),
                            "alt_ids": ""
                        }
    except Exception as e:
        print(f"Failed to parse JSON {file_path}: {e}")
    return nodes

def extract_nodes_from_xml(file_path):
    """Extract nodes from XML using streaming"""
    nodes = {}
    try:
        for event, elem in ET.iterparse(file_path, events=("end",)):
            for child in elem:
                for col in child.attrib:
                    if any(k in col.upper() for k in ENTITY_ID_KEYS):
                        raw_id = str(child.attrib[col]).strip()
                        if raw_id not in nodes:
                            nodes[raw_id] = {
                                "entity_type": infer_entity_type(col),
                                "raw_id": raw_id,
                                "name": "",
                                "alt_ids": ""
                            }
            elem.clear()
    except Exception as e:
        print(f"Failed to parse XML {file_path}: {e}")
    return nodes

# ----------------------
# EDGE EXTRACTION
# ----------------------

def extract_edges_from_dataframe(df, dataset_name):
    """Extract candidate edges from a DataFrame"""
    edges = []
    if df is None:
        return edges

    id_columns = extract_candidate_columns(df)
    if len(id_columns) < 2:
        return edges

    for i, src_col in enumerate(id_columns):
        src_type = infer_entity_type(src_col)
        for tgt_col in id_columns[i+1:]:
            tgt_type = infer_entity_type(tgt_col)
            relation = assign_relation(src_type, tgt_type)
            for _, row in df.iterrows():
                src_val = row.get(src_col)
                tgt_val = row.get(tgt_col)
                if pd.isna(src_val) or pd.isna(tgt_val):
                    continue
                edges.append((str(src_val).strip(), str(tgt_val).strip(), relation, dataset_name))
    return edges

def extract_edges_from_fasta(file_path, dataset_name, max_nodes=1000):
    """Extract edges from FASTA (limited connections)"""
    edges = []
    if not BIOPYTHON_AVAILABLE:
        return edges
    try:
        ids = []
        for i, record in enumerate(SeqIO.parse(file_path, "fasta")):
            ids.append(record.id.strip())
            if i+1 >= max_nodes:
                break
        # Create fully connected edges between sampled nodes
        for i, src in enumerate(ids):
            for tgt in ids[i+1:]:
                edges.append((src, tgt, "INTERACTS_WITH", dataset_name))
    except Exception as e:
        print(f"Failed to parse FASTA {file_path}: {e}")
    return edges

def extract_edges_from_json(file_path, dataset_name, max_rows=1000):
    """Extract edges from JSON file"""
    edges = []
    try:
        with open(file_path) as f:
            try:
                data = json.load(f)
                if isinstance(data, dict):
                    data = [data]
            except json.JSONDecodeError:
                f.seek(0)
                data = []
                for i, line in enumerate(f):
                    data.append(json.loads(line))
                    if i+1 >= max_rows:
                        break
        for obj in data:
            id_columns = [k for k in obj.keys() if any(k in k.upper() for k in ENTITY_ID_KEYS)]
            if len(id_columns) < 2:
                continue
            for i, src_col in enumerate(id_columns):
                src_type = infer_entity_type(src_col)
                for tgt_col in id_columns[i+1:]:
                    tgt_type = infer_entity_type(tgt_col)
                    relation = assign_relation(src_type, tgt_type)
                    src_val = obj.get(src_col)
                    tgt_val = obj.get(tgt_col)
                    if src_val and tgt_val:
                        edges.append((str(src_val).strip(), str(tgt_val).strip(), relation, dataset_name))
    except Exception as e:
        print(f"Failed to parse JSON {file_path}: {e}")
    return edges

def extract_edges_from_xml(file_path, dataset_name):
    """Extract edges from XML using streaming"""
    edges = []
    try:
        for event, elem in ET.iterparse(file_path, events=("end",)):
            id_columns = [k for k in elem.attrib if any(k in k.upper() for k in ENTITY_ID_KEYS)]
            if len(id_columns) >= 2:
                for i, src_col in enumerate(id_columns):
                    src_type = infer_entity_type(src_col)
                    for tgt_col in id_columns[i+1:]:
                        tgt_type = infer_entity_type(tgt_col)
                        relation = assign_relation(src_type, tgt_type)
                        src_val = elem.attrib.get(src_col)
                        tgt_val = elem.attrib.get(tgt_col)
                        if src_val and tgt_val:
                            edges.append((src_val.strip(), tgt_val.strip(), relation, dataset_name))
            elem.clear()
    except Exception as e:
        print(f"Failed to parse XML {file_path}: {e}")
    return edges

# ----------------------
# MAIN PROCESS
# ----------------------

def main():
    ensure_processed_dir()
    datasets = [d for d in os.listdir(RAW_DATA_DIR) if os.path.isdir(os.path.join(RAW_DATA_DIR, d))]

    print(f"Processing {len(datasets)} datasets...")
    all_nodes_global = {}
    all_edges_global = []

    for dataset_name in tqdm(datasets):
        dataset_path = os.path.join(RAW_DATA_DIR, dataset_name)
        files = [os.path.join(dataset_path, f) for f in os.listdir(dataset_path)
                 if f.lower().endswith((".csv", ".tsv", ".json", ".xml", ".fasta"))]

        nodes_dataset = {}
        edges_dataset = []

        for file_path in files:
            if file_path.lower().endswith((".csv", ".tsv")):
                df = safe_read_table(file_path, nrows=1000)  # chunking
                nodes_dataset.update(extract_nodes_from_dataframe(df))
                edges_dataset.extend(extract_edges_from_dataframe(df, dataset_name))

            elif file_path.lower().endswith(".json"):
                nodes_dataset.update(extract_nodes_from_json(file_path))
                edges_dataset.extend(extract_edges_from_json(file_path, dataset_name))

            elif file_path.lower().endswith(".xml"):
                nodes_dataset.update(extract_nodes_from_xml(file_path))
                edges_dataset.extend(extract_edges_from_xml(file_path, dataset_name))

            elif file_path.lower().endswith(".fasta"):
                nodes_dataset.update(extract_nodes_from_fasta(file_path))
                edges_dataset.extend(extract_edges_from_fasta(file_path, dataset_name))

        # Save per-dataset CSVs
        nodes_df = pd.DataFrame(nodes_dataset.values())
        nodes_file = os.path.join(PROCESSED_DIR, f"{dataset_name}_nodes.csv")
        nodes_df.to_csv(nodes_file, index=False)

        edges_df = pd.DataFrame(edges_dataset, columns=["source_raw_id", "target_raw_id", "relation", "dataset"])
        edges_file = os.path.join(PROCESSED_DIR, f"{dataset_name}_edges.csv")
        edges_df.to_csv(edges_file, index=False)

        # Add to global
        all_nodes_global.update(nodes_dataset)
        all_edges_global.extend(edges_dataset)

    # Optional: global combined CSVs
    global_nodes_file = os.path.join(PROCESSED_DIR, "all_nodes.csv")
    pd.DataFrame(all_nodes_global.values()).to_csv(global_nodes_file, index=False)

    global_edges_file = os.path.join(PROCESSED_DIR, "all_edges.csv")
    pd.DataFrame(all_edges_global, columns=["source_raw_id", "target_raw_id", "relation", "dataset"]).to_csv(global_edges_file, index=False)

    print(f"✓ Nodes and edges extraction completed.")
    print(f"Global nodes: {len(all_nodes_global)}, Global edges: {len(all_edges_global)}")

    # Optional: build NetworkX graph
    G = nx.Graph()
    for e in all_edges_global:
        G.add_edge(e[0], e[1], relation=e[2], dataset=e[3])
    print(f"Graph built with nodes: {G.number_of_nodes()}, edges: {G.number_of_edges()}")

if __name__ == "__main__":
    main()
