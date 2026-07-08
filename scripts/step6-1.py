import pandas as pd
import numpy as np

# ================= SETTINGS =================

RELATIONSHIP_FILE = "data/processed/relationships.csv"

MIN_DATASETS_PER_DISEASE = 5
MIN_DISEASES_PER_GENE = 3
MIN_EDGE_SUPPORT = 2

# confidence weights
W_EDGE = 0.5
W_DISEASE = 0.3
W_GENE = 0.2

# ==========================================


def load_edges(path):
    df = pd.read_csv(path)

    required = {
        "subject",
        "subject_type",
        "object",
        "object_type",
        "relationship",
        "source_file",
    }

    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {missing}")

    return df


def get_gene_disease_edges(df):
    """
    Normalize direction so we always have:
    gene -> disease
    """

    gene_as_subject = df[
        (df.subject_type == "GENE") &
        (df.object_type == "DISEASE")
    ].copy()

    gene_as_object = df[
        (df.subject_type == "DISEASE") &
        (df.object_type == "GENE")
    ].copy()

    gene_as_object = gene_as_object.rename(
        columns={
            "subject": "object",
            "object": "subject",
            "subject_type": "object_type",
            "object_type": "subject_type",
        }
    )

    combined = pd.concat([gene_as_subject, gene_as_object], ignore_index=True)

    combined = combined.rename(columns={
        "subject": "gene",
        "object": "disease"
    })

    return combined


# --------------------------------------------------------

def add_support_metrics(edges):
    """
    Adds:
      edge_support
      disease_dataset_count
      gene_disease_count
    """

    # Edge support
    edge_support = (
        edges.groupby(["gene", "disease"])["source_file"]
        .nunique()
        .reset_index(name="edge_support")
    )

    edges = edges.merge(edge_support, on=["gene", "disease"])

    # Disease robustness
    disease_support = (
        edges.groupby("disease")["source_file"]
        .nunique()
        .reset_index(name="disease_dataset_count")
    )

    edges = edges.merge(disease_support, on="disease")

    # Gene connectivity
    gene_support = (
        edges.groupby("gene")["disease"]
        .nunique()
        .reset_index(name="gene_disease_count")
    )

    edges = edges.merge(gene_support, on="gene")

    return edges


def compute_confidence(edges):

    edges["confidence"] = (
        W_EDGE * np.log1p(edges["edge_support"]) +
        W_DISEASE * np.log1p(edges["disease_dataset_count"]) +
        W_GENE * np.log1p(edges["gene_disease_count"])
    )

    return edges


# --------------------------------------------------------

def main():

    print("Loading relationships...")
    df = load_edges(RELATIONSHIP_FILE)

    print("Extracting gene–disease edges...")
    gd = get_gene_disease_edges(df)

    print(f"Initial edges: {len(gd)}")

    print("Adding support metrics...")
    gd = add_support_metrics(gd)

    print("Applying hard filters...")

    gd = gd[gd.disease_dataset_count >= MIN_DATASETS_PER_DISEASE]
    gd = gd[gd.edge_support >= MIN_EDGE_SUPPORT]
    gd = gd[gd.gene_disease_count >= MIN_DISEASES_PER_GENE]

    print(f"After filters: {len(gd)}")

    print("Computing confidence scores...")
    gd = compute_confidence(gd)

    gd = gd.sort_values("confidence", ascending=False)

    gd.to_csv("trusted_gene_disease_edges.csv", index=False)

    print("\nSaved: trusted_gene_disease_edges.csv")

    print("\nSummary:")
    print("Genes:", gd.gene.nunique())
    print("Diseases:", gd.disease.nunique())
    print("Edges:", len(gd))

    print("\nTop 10 edges by confidence:")
    print(gd[["gene", "disease", "confidence"]].head(10))


if __name__ == "__main__":
    main()
