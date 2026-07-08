import pandas as pd
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.sparse.linalg import eigsh

# --- CONFIGURATION ---
# 1. Repo-relative path, resolved at import (was a hard-coded absolute path)
from config import RELATIONSHIPS_CSV
ABSOLUTE_PATH = str(RELATIONSHIPS_CSV)
# 2. Relative path (Standard project structure)
RELATIVE_PATH = "data/processed/relationships.csv"

OUTPUT_REPORT = "spectral_integrity_report.txt"
SPECTRAL_PLOT = "spectral_partition_map.png"

def find_data_file():
    """Attempts to locate the relationships file using multiple methods."""
    if os.path.exists(ABSOLUTE_PATH):
        print(f"File found via Absolute Path: {ABSOLUTE_PATH}")
        return ABSOLUTE_PATH
    
    if os.path.exists(RELATIVE_PATH):
        print(f"File found via Relative Path: {RELATIVE_PATH}")
        return RELATIVE_PATH
    
    print("Paths not found. Starting deep search in current directory...")
    target = "relationships.csv"
    for root, dirs, files in os.walk(os.getcwd()):
        if target in files:
            found = os.path.join(root, target)
            print(f"File found via Deep Search: {found}")
            return found
            
    return None

def run_spectral_analysis():
    rel_path = find_data_file()
    if not rel_path:
        print("ERROR: 'relationships.csv' could not be located.")
        return

    print("Step 1: Ingesting Graph for Spectral Decomposition...")
    df = pd.read_csv(rel_path, low_memory=False)
    
    start_col = ':START_ID' if ':START_ID' in df.columns else 'START_ID'
    end_col = ':END_ID' if ':END_ID' in df.columns else 'END_ID'
    
    if start_col not in df.columns or end_col not in df.columns:
        print(f"Error: Columns {start_col} or {end_col} not found.")
        return

    # Create Graph and extract Largest Connected Component
    G_full = nx.from_pandas_edgelist(df, start_col, end_col)
    G = G_full.subgraph(max(nx.connected_components(G_full), key=len)).copy()
    
    print(f"Analyzing Main Component: {G.number_of_nodes()} nodes.")

    # Step 2: Algebraic Connectivity (The Fiedler Value)
    print("Step 2: Computing Laplacian Eigenvalues (Fiedler Value)...")
    
    # FIX: Updated for compatibility with newer SciPy/NetworkX csr_array
    L = nx.laplacian_matrix(G)
    # Convert to float64 so the eigenvalue solver (eigsh) can process it
    L = L.astype(float) 
    
    # k=2 for the 2nd smallest eigenvalue (The Fiedler Value)
    # 'SM' = Smallest Magnitude
    eigenvalues, eigenvectors = eigsh(L, k=2, which='SM')
    fiedler_value = eigenvalues[1] 
    fiedler_vector = eigenvectors[:, 1] 

    # Step 3: Spectral Partitioning
    partition = {node: (1 if fiedler_vector[i] > 0 else 0) for i, node in enumerate(G.nodes())}
    
    # Step 4: Measuring Global Influence via Eigen-Centrality
    print("Step 3: Measuring Global Influence...")
    centrality = nx.eigenvector_centrality_numpy(G)

    # Step 5: Visualization
    print("Step 4: Visualizing Spectral Rifts...")
    # Sampling top 300 nodes for a clean plot
    top_nodes = sorted(centrality.items(), key=lambda x: x[1], reverse=True)[:300]
    sub_G = G.subgraph([n[0] for n in top_nodes])
    
    plt.figure(figsize=(10, 8), facecolor='#0f172a')
    ax = plt.gca()
    ax.set_facecolor('#0f172a')
    
    pos = nx.spring_layout(sub_G, seed=42)
    colors = [partition[n] for n in sub_G.nodes()]
    sizes = [centrality[n] * 6000 for n in sub_G.nodes()]

    nx.draw_networkx_edges(sub_G, pos, alpha=0.1, edge_color='white')
    nx.draw_networkx_nodes(sub_G, pos, node_size=sizes, node_color=colors, cmap=plt.cm.coolwarm)
    
    plt.title(f"Spectral Partitioning: Fiedler Value = {fiedler_value:.6f}", color='white', fontsize=12)
    plt.axis('off')
    plt.savefig(SPECTRAL_PLOT, dpi=300, facecolor='#0f172a')

    # Step 6: Generate the Report
    with open(OUTPUT_REPORT, "w") as f:
        f.write("=== SPECTRAL BIOLOGICAL ANALYSIS ===\n")
        f.write(f"Source File: {rel_path}\n")
        f.write(f"Fiedler Value (Algebraic Connectivity): {fiedler_value:.8f}\n\n")
        f.write("INTERPRETATION:\n")
        if fiedler_value < 0.001:
            f.write("- The network is 'Brittle'. It is very easy to disconnect segments.\n")
        else:
            f.write("- The network is 'Resilient'. It has high structural redundancy.\n")
        
        f.write("\nTOP SPECTRAL TARGETS:\n")
        for node, val in sorted(centrality.items(), key=lambda x: x[1], reverse=True)[:10]:
            f.write(f"ID: {node} | Centrality: {val:.6f}\n")

    print(f"Analysis complete. Report: {OUTPUT_REPORT}, Image: {SPECTRAL_PLOT}")

if __name__ == "__main__":
    print("=== STARTING SPECTRAL INTELLIGENCE ENGINE ===")
    run_spectral_analysis()