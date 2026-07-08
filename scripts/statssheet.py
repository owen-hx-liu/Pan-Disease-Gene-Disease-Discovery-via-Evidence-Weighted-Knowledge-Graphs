import pandas as pd
import networkx as nx
from networkx.algorithms import community as nx_community

# 1. Load the data
INPUT_PATH = "data/processed/consensus similar table.csv"
OUTPUT_REPORT = "science_fair_data_report.txt"

print(f"Reading data from {INPUT_PATH}...")
df = pd.read_csv(INPUT_PATH)

# 2. Build the Graph
G = nx.Graph()
for _, row in df.iterrows():
    G.add_edge(row['head'], row['tail'], weight=row['old_score'], relation=row['relation'])

# 3. Comprehensive Analytics
print("Running deep analysis...")

# Centrality (The 'Hub' Analysis)
degree_dict = dict(G.degree(weight='weight'))
betweenness_dict = nx.betweenness_centrality(G, weight='weight')
try:
    eigenvector_dict = nx.eigenvector_centrality(G, weight='weight', max_iter=1000)
except:
    eigenvector_dict = {node: 0 for node in G.nodes()}

# Community Detection
communities = nx_community.louvain_communities(G, weight='weight')
# Sort communities by size (largest first)
communities = sorted(list(communities), key=len, reverse=True)

# Motif Discovery (Triangles)
triangles = nx.triangles(G)
total_triangles = sum(triangles.values()) // 3

# 4. Write the Report File
with open(OUTPUT_REPORT, "w") as f:
    f.write("=== SCIENCE FAIR GENE-DISEASE NETWORK REPORT ===\n")
    f.write(f"Total Nodes: {G.number_of_nodes()}\n")
    f.write(f"Total Edges: {G.number_of_edges()}\n")
    f.write(f"Total Communities Found: {len(communities)}\n")
    f.write(f"Biological Triangles (Motifs): {total_triangles}\n\n")

    f.write("--- TOP 10 INFLUENTIAL HUBS (Eigenvector Centrality) ---\n")
    f.write("These nodes are connected to other highly important nodes.\n")
    top_eigen = sorted(eigenvector_dict.items(), key=lambda x: x[1], reverse=True)[:10]
    for i, (node, score) in enumerate(top_eigen, 1):
        f.write(f"{i}. {node}: {score:.4f}\n")

    f.write("\n--- TOP 10 BRIDGE NODES (Betweenness Centrality) ---\n")
    f.write("These nodes connect different biological modules together.\n")
    top_between = sorted(betweenness_dict.items(), key=lambda x: x[1], reverse=True)[:10]
    for i, (node, score) in enumerate(top_between, 1):
        f.write(f"{i}. {node}: {score:.4f}\n")

    f.write("\n--- COMMUNITY BREAKDOWN ---\n")
    for i, comm in enumerate(communities[:5], 1): # Showing top 5 largest communities
        f.write(f"Community #{i} (Size: {len(comm)} nodes):\n")
        f.write(f"  Nodes: {', '.join(list(comm)[:15])}...") # List first 15 nodes
        f.write("\n\n")

    f.write("--- MOTIF DISCOVERY ---\n")
    top_triangles = sorted(triangles.items(), key=lambda x: x[1], reverse=True)[:5]
    f.write("Nodes most frequently found in biological triangles:\n")
    for node, count in top_triangles:
        f.write(f" - {node}: {count} triangles\n")

print(f"Report successfully generated: {OUTPUT_REPORT}")