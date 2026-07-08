import pandas as pd
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import os

# --- CONFIGURATION ---
from config import RELATIONSHIPS_CSV
ABSOLUTE_PATH = str(RELATIONSHIPS_CSV)  # repo-relative, resolved at import
RELATIVE_PATH = "data/processed/relationships.csv"

OUTPUT_DATA = "data/processed/functional_influence_map.csv"
OUTPUT_IMAGE = "data/processed/disease_influence_heatmap.png"

def find_data_file():
    """Dynamically locates the relationships.csv file."""
    if os.path.exists(ABSOLUTE_PATH):
        return ABSOLUTE_PATH
    if os.path.exists(RELATIVE_PATH):
        return RELATIVE_PATH
    target = "relationships.csv"
    for root, dirs, files in os.walk(os.getcwd()):
        if target in files:
            return os.path.join(root, target)
    return None

def ensure_output_dirs():
    os.makedirs("data/processed", exist_ok=True)

def calculate_diffusion_influence(G, seed_node, alpha=0.85):
    """Calculates PageRank-based influence starting from a specific seed."""
    if seed_node not in G:
        print(f"Warning: Seed node {seed_node} not found.")
        return None
    
    personalization = {node: 0 for node in G.nodes()}
    personalization[seed_node] = 1
    
    influence_scores = nx.pagerank(G, alpha=alpha, personalization=personalization)
    return influence_scores

def generate_heatmap(G, influence_map, seed_id):
    """Improved visualization with Log-Scaling to prevent 'blackout' issues."""
    print("Generating High-Fidelity Influence Map...")
    
    try:
        neighbors = nx.single_source_shortest_path_length(G, seed_id, cutoff=2)
        neighborhood_nodes = list(neighbors.keys())
        
        if len(neighborhood_nodes) > 400:
            neighborhood_nodes = sorted(neighborhood_nodes, key=lambda x: influence_map.get(x, 0), reverse=True)[:400]
        
        sub_G = G.subgraph(neighborhood_nodes).copy()
    except Exception as e:
        print(f"Neighborhood search failed: {e}. Falling back to top global nodes.")
        top_nodes = sorted(influence_map.items(), key=lambda x: x[1], reverse=True)[:200]
        sub_G = G.subgraph([n[0] for n in top_nodes]).copy()

    # 2. Layout
    pos = nx.spring_layout(sub_G, k=0.6/np.sqrt(len(sub_G.nodes())), iterations=60, seed=42)
    
    scores = np.array([influence_map.get(node, 0) for node in sub_G.nodes()])
    
    log_scores = np.log10(scores + 1e-10)
    
    norm_scores = (log_scores - log_scores.min()) / (log_scores.max() - log_scores.min() + 1e-9)
    
    plt.figure(figsize=(14, 12), facecolor='#0a0a0c')
    ax = plt.gca()
    ax.set_facecolor('#0a0a0c')
    
    nx.draw_networkx_edges(sub_G, pos, alpha=0.2, edge_color='#555555', width=0.6)
    
    node_sizes = 150 + (norm_scores * 3500)
    
    nx.draw_networkx_nodes(
        sub_G, 
        pos, 
        node_size=node_sizes,
        node_color=norm_scores,
        cmap=plt.cm.magma, 
        alpha=0.85
    )
    
    top_targets = sorted(list(sub_G.nodes()), key=lambda x: influence_map.get(x, 0), reverse=True)[:8]
    labels = {node: node for node in top_targets}
    nx.draw_networkx_labels(sub_G, pos, labels=labels, font_size=10, font_color='white', font_weight='bold')

    if seed_id in pos:
        nx.draw_networkx_nodes([seed_id], pos, node_size=1800, node_color='#00ffff', edgecolors='white', linewidths=2.5)
        plt.annotate("DISEASE SEED", xy=pos[seed_id], xytext=(20, 20), 
                     textcoords='offset points', color='#00ffff', fontweight='bold', fontsize=13,
                     arrowprops=dict(arrowstyle='->', color='#00ffff'))

    plt.title(f"Functional Influence Propagation\nLog-Scaled Map: {seed_id}", 
              color='white', fontsize=20, fontweight='bold', pad=30)
    
    sm = plt.cm.ScalarMappable(cmap=plt.cm.magma, norm=plt.Normalize(vmin=0, vmax=1))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.03, pad=0.04)
    cbar.set_label('Log-Scaled Influence Strength', color='white', fontsize=12)
    cbar.ax.yaxis.set_tick_params(color='white', labelcolor='white')

    plt.axis('off')
    plt.tight_layout()
    plt.savefig(OUTPUT_IMAGE, dpi=300, facecolor='#0a0a0c')
    print(f"Heat Map fixed and saved as {OUTPUT_IMAGE}")
    plt.close()

def run_advanced_topology():
    ensure_output_dirs()
    rel_path = find_data_file()
    
    if not rel_path:
        print("ERROR: 'relationships.csv' not found.")
        return

    print("Step 1: Loading Graph...")
    df = pd.read_csv(rel_path, low_memory=False)
    
    start_col = next((c for c in df.columns if 'START_ID' in c), "START_ID")
    end_col = next((c for c in df.columns if 'END_ID' in c), "END_ID")
    
    G = nx.from_pandas_edgelist(df, start_col, end_col)
    seed_id = "MONDO:0018874" 
    
    if seed_id not in G:
        seed_id = str(df[end_col].iloc[0])
        print(f"Notice: Specific seed not found. Using default: {seed_id}")
    
    print(f"Step 2: Calculating Diffusion for {seed_id}...")
    influence_map = calculate_diffusion_influence(G, seed_id)
    
    if influence_map:
        inf_df = pd.DataFrame(list(influence_map.items()), columns=['node_id', 'influence_score'])
        inf_df = inf_df.sort_values(by='influence_score', ascending=False)
        inf_df.to_csv(OUTPUT_DATA, index=False)
        generate_heatmap(G, influence_map, seed_id)

if __name__ == "__main__":
    print("=== STARTING BIOLOGICAL DIFFUSION ANALYSIS (FIXED VISUALS) ===")
    run_advanced_topology()