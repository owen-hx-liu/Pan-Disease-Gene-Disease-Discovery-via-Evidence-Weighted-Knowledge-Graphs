import pandas as pd
import networkx as nx
from networkx.algorithms import community as nx_community
from pyvis.network import Network

df = pd.read_csv("data/processed/consensus similar table.csv") 

G = nx.Graph()
for _, row in df.iterrows():
    G.add_edge(row['head'], row['tail'], weight=row['old_score'], relation=row['relation'])

print(f"Graph loaded: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

print("Calculating centralities...")
degree_dict = dict(G.degree(weight='weight'))
betweenness_dict = nx.betweenness_centrality(G, weight='weight')
try:
    eigenvector_dict = nx.eigenvector_centrality(G, weight='weight', max_iter=1000)
except:
    eigenvector_dict = {node: 0 for node in G.nodes()}

nx.set_node_attributes(G, degree_dict, 'degree')
nx.set_node_attributes(G, betweenness_dict, 'betweenness')
nx.set_node_attributes(G, eigenvector_dict, 'eigenvector')

top_hubs = sorted(eigenvector_dict.items(), key=lambda x: x[1], reverse=True)[:5]
print("Top 5 Influential Hubs (Eigenvector):", [h[0] for h in top_hubs])

print("Finding communities...")
communities = nx_community.louvain_communities(G, weight='weight')
partition = {}
for community_id, node_set in enumerate(communities):
    for node in node_set:
        partition[node] = community_id
nx.set_node_attributes(G, partition, 'group')

print("Scanning for motifs...")
triangles = nx.triangles(G)
nx.set_node_attributes(G, triangles, 'triangles')

net = Network(height="800px", width="100%", bgcolor="#222222", font_color="white", select_menu=True, filter_menu=True)
net.from_nx(G)

for node in net.nodes:
    node_id = node['id']
    
    influence = eigenvector_dict.get(node_id, 0)
    node['size'] = 15 + (influence * 100) 
    
    node['title'] = (
        f"<b>Node:</b> {node_id}<br>"
        f"<b>Community:</b> {partition.get(node_id, 0)}<br>"
        f"<b>Influence (Eigen):</b> {influence:.4f}<br>"
        f"<b>Bridge Score:</b> {betweenness_dict.get(node_id, 0):.4f}<br>"
        f"<b>Part of {triangles.get(node_id, 0)} biological triangles</b>"
    )
    
    node['group'] = partition.get(node_id, 0)

for edge in net.edges:
    score = edge.get('weight', 0)
    edge['width'] = score * 3
    edge['color'] = '#4a90e2' if score > 0.7 else '#888888'

net.show_buttons(filter_=['physics', 'nodes', 'edges'])

output_file = "science_fair_graph_explorer.html"
net.save_graph(output_file)

print("-" * 30)
print(f"Analysis Complete!")
print(f"Total Communities: {len(communities)}")
print(f"Triangles Found: {sum(triangles.values()) // 3}")
print(f"File Saved: {output_file}")