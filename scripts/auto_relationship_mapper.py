import os
import pandas as pd
from tqdm import tqdm
import networkx as nx
from Bio import SeqIO  # For FASTA
import json
import xml.etree.ElementTree as ET

RAW_DATA_DIR = "data/raw"
PROCESSED_DIR = "data/processed"
OUTPUT_FILE = os.path.join(PROCESSED_DIR, "edges_inferred.csv")

ENTITY_ID_KEYS = ["GENE", "GENE_ID", "ENTREZ_ID", "DISEASE", "DISEASE_ID", "VARIANT", "CHEMICAL"]

def safe_read_table(file_path, nrows=1000):
    """Read CSV or TSV safely, skip bad lines, and limit to nrows"""
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


def extract_candidate_edges_from_df(df):
    """Return candidate edges as tuples (source, target)"""
    candidate_edges = []
    if df is None:
        return candidate_edges

    id_columns = [c for c in df.columns if any(k in c.upper() for k in ENTITY_ID_KEYS)]
    if len(id_columns) < 2:
        return candidate_edges

    # Stream rows efficiently
    for i, src in enumerate(id_columns):
        for tgt in id_columns[i+1:]:
            for a, b in zip(df[src], df[tgt]):
                if pd.notnull(a) and pd.notnull(b):
                    candidate_edges.append((a, b))
    return candidate_edges

def extract_edges_from_fasta(file_path, dataset_name, max_nodes=1000):
    """Extract candidate edges from FASTA, limited to max_nodes to avoid explosion"""
    candidate_edges = []
    try:
        ids = []
        for i, record in enumerate(SeqIO.parse(file_path, "fasta")):
            ids.append(record.id)
            if i+1 >= max_nodes:
                break

        # Fully connect sampled IDs (or skip edges if you just want nodes)
        for i, src in enumerate(ids):
            for tgt in ids[i+1:]:
                candidate_edges.append((src, tgt, dataset_name))

    except Exception as e:
        print(f"Failed to parse FASTA {file_path}: {e}")
    return candidate_edges


def extract_edges_from_json(file_path, dataset_name, max_rows=1000):
    candidate_edges = []
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

        # Process candidate edges (same as before)
        for obj in data:
            id_fields = [k for k in obj.keys() if any(key in k.upper() for key in ENTITY_ID_KEYS)]
            if len(id_fields) < 2:
                continue
            for i, src in enumerate(id_fields):
                for tgt in id_fields[i+1:]:
                    a = obj.get(src)
                    b = obj.get(tgt)
                    if a and b:
                        candidate_edges.append((a, b, dataset_name))
    except Exception as e:
        print(f"Failed to parse JSON {file_path}: {e}")
    return candidate_edges


def extract_edges_from_xml(file_path, dataset_name):
    """Infer edges from XML elements with candidate ID attributes or child tags"""
    candidate_edges = []
    try:
        for event, elem in ET.iterparse(file_path, events=("end",)):
            # Collect candidate ID values from attributes or children
            values = []
            for k, v in elem.attrib.items():
                if any(key in k.upper() for key in ENTITY_ID_KEYS) and v:
                    values.append(v)
            for child in elem:
                if child.text and any(key in child.tag.upper() for key in ENTITY_ID_KEYS):
                    values.append(child.text)
            # Make all pairs from collected IDs
            for i, src in enumerate(values):
                for tgt in values[i+1:]:
                    candidate_edges.append((src, tgt, dataset_name))
            elem.clear()  # free memory
    except Exception as e:
        print(f"Failed to parse XML {file_path}: {e}")
    return candidate_edges

def main():
    if not os.path.exists(PROCESSED_DIR):
        os.makedirs(PROCESSED_DIR)

    all_edges = []
    datasets = [d for d in os.listdir(RAW_DATA_DIR) if os.path.isdir(os.path.join(RAW_DATA_DIR, d))]

    print(f"Processing {len(datasets)} datasets...")
    for dataset_name in tqdm(datasets):
        dataset_path = os.path.join(RAW_DATA_DIR, dataset_name)
        files = [
            os.path.join(dataset_path, f)
            for f in os.listdir(dataset_path)
            if f.lower().endswith((".csv", ".tsv", ".fasta", ".json", ".xml"))
        ]

        for file_path in files:
            if file_path.lower().endswith((".csv", ".tsv")):
                df = safe_read_table(file_path, nrows=1000)
                candidate_edges = extract_candidate_edges_from_df(df)
                all_edges.extend([(src, tgt, dataset_name) for src, tgt in candidate_edges])
            elif file_path.lower().endswith(".fasta"):
                fasta_edges = extract_edges_from_fasta(file_path, dataset_name)
                all_edges.extend(fasta_edges)
            elif file_path.lower().endswith(".json"):
                json_edges = extract_edges_from_json(file_path, dataset_name)
                all_edges.extend(json_edges)
            elif file_path.lower().endswith(".xml"):
                xml_edges = extract_edges_from_xml(file_path, dataset_name)
                all_edges.extend(xml_edges)
            else:
                print(f"Skipping unsupported file format: {file_path}")

    edges_df = pd.DataFrame(all_edges, columns=["source", "target", "dataset"])

    if len(edges_df) == 0:
        print("No candidate edges discovered.")
        return

    edges_df.to_csv(OUTPUT_FILE, index=False)
    print(f"✓ Candidate edges saved to: {OUTPUT_FILE}")
    print(f"Total candidate edges: {len(edges_df)}")

    G = nx.Graph()
    for _, row in edges_df.iterrows():
        G.add_edge(row["source"], row["target"], dataset=row["dataset"])
    print("Graph built with nodes:", G.number_of_nodes(), "edges:", G.number_of_edges())

if __name__ == "__main__":
    main()
