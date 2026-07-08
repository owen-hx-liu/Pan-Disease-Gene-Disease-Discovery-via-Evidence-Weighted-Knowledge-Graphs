import pandas as pd
import networkx as nx

def find_high_value_pairs():
    """
    Scans the processed relationship data to find pairs of nodes 
    that have strong biological consensus for AI testing.
    """
    # 1. Load the data
    # We use the consensus table as it already has the scores we need
    file_path = "data/processed/consensus similar table.csv"
    
    print(f"Loading relationship data from {file_path}...")
    try:
        df = pd.read_csv(file_path)
    except FileNotFoundError:
        print("Error: Could not find 'consensus similar table.csv'.")
        print("Please check your file path or run your processing scripts first.")
        return

    # 2. Filter for high-confidence connections
    # We want pairs where 'old_score' is high (close to 1.0 or high integer)
    # and we want a variety of relation types.
    print("Filtering for high-consensus pairs...")
    
    # Sort by score descending to get the strongest links first
    top_pairs_df = df.sort_values(by='old_score', ascending=False).head(50)
    
    # 3. Build a mini-graph to check for 'importance' (Centrality)
    G = nx.Graph()
    for _, row in df.head(10000).iterrows(): # Using a subset for speed
        G.add_edge(row['head'], row['tail'], weight=row['old_score'])
    
    print("Calculating node influence to prioritize pairs...")
    pagerank = nx.pagerank(G, weight='weight')

    # 4. Select the best pairs to show the user
    # A 'Best Pair' = High Link Score + High Node Importance
    results = []
    for _, row in top_pairs_df.iterrows():
        u, v = row['head'], row['tail']
        importance = pagerank.get(u, 0) + pagerank.get(v, 0)
        results.append({
            'Node A': u,
            'Node B': v,
            'Relation': row['relation'],
            'Score': row['old_score'],
            'Importance': importance
        })

    # Sort results by importance so the user sees the 'Hub' pairs first
    results = sorted(results, key=lambda x: x['Importance'], reverse=True)

    print("\n" + "="*60)
    print("🚀 TOP 10 RECOMMENDED TEST PAIRS FOR YOUR AI APP")
    print("="*60)
    print(f"{'Node A (ID)':<20} | {'Node B (ID)':<20} | {'Type':<15}")
    print("-" * 60)
    
    for res in results[:10]:
        print(f"{res['Node A']:<20} | {res['Node B']:<20} | {res['Relation']:<15}")
        
    print("="*60)
    print("\nINSTRUCTIONS:")
    print("1. Copy one ID from 'Node A' and one from 'Node B'.")
    print("2. Paste them into your Streamlit App.")
    print("3. If the AI gives a high score (e.g. >80%), it means the model")
    print("   correctly recognizes this high-consensus biological link!")

if __name__ == "__main__":
    find_high_value_pairs()