import pandas as pd
import numpy as np
import networkx as nx
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, roc_auc_score
import joblib
import ast 
import os
import random

def find_file_globally(filename):
    """Search for a file in the current directory and all subdirectories."""
    print(f"--- Searching for {filename} ---")
    for root, dirs, files in os.walk(os.getcwd()):
        if filename in files:
            found_path = os.path.join(root, filename)
            print(f"SUCCESS: Found {filename} at: {found_path}")
            return found_path
    print(f"FAILED: Could not find {filename} anywhere in {os.getcwd()}")
    return None

def load_data():
    print(f"DEBUG: Python's current 'standing point' (CWD): {os.getcwd()}")
    
    # Let's try to find the files by searching the folders
    rel_path = find_file_globally("relationships.csv")
    node_path = find_file_globally("nodes.csv")
    
    if not rel_path:
        print("!!! ERROR: relationships.csv is missing from your project folder.")
        return None
    
    if not node_path:
        print("!!! WARNING: nodes.csv not found, but we can continue with relationships.")
    else:
        print("MATCH: nodes.csv is confirmed and ready.")

    print(f"Step 1: Loading relationship data from {rel_path}...")
    edges_df = pd.read_csv(rel_path, low_memory=False)
    
    # Detect column names (Neo4j often adds extra characters)
    start_col = next((c for c in edges_df.columns if 'START_ID' in c), None)
    end_col = next((c for c in edges_df.columns if 'END_ID' in c), None)
    weight_col = next((c for c in edges_df.columns if 'weight' in c.lower()), None)
    
    if not start_col or not end_col:
        print(f"Error: Could not find ID columns. Found these instead: {list(edges_df.columns)}")
        return None
        
    if weight_col:
        edges_df['is_high_consensus'] = (edges_df[weight_col] >= 7).astype(int)
    else:
        edges_df['is_high_consensus'] = 1
        
    return edges_df, start_col, end_col

def load_embeddings():
    emb_path = find_file_globally("node2vec_embeddings.csv")
    if not emb_path:
        print("!!! ERROR: node2vec_embeddings.csv not found!")
        return None

    try:
        print(f"Step 2: Loading embeddings from {emb_path}...")
        emb_df = pd.read_csv(emb_path)
        embeddings_dict = {}
        
        for _, row in emb_df.iterrows():
            node_id = str(row['node_id'])
            raw_emb = row['embedding']
            
            if isinstance(raw_emb, str):
                try:
                    vector = np.array(ast.literal_eval(raw_emb.strip()))
                    embeddings_dict[node_id] = vector
                except:
                    continue
            else:
                embeddings_dict[node_id] = np.array(raw_emb)
                
        return embeddings_dict
    except Exception as e:
        print(f"Error processing embeddings: {e}")
        return None

def feature_engineering(G, edges_df, embeddings_dict, start_col, end_col):
    print("Step 3: Engineering features (Creating the AI study guide)...")
    features = []
    targets = []

    # Use a safe sample size for training
    sample_size = min(50000, len(edges_df))
    sampled_edges = edges_df.sample(n=sample_size)

    for _, row in sampled_edges.iterrows():
        u, v = str(row[start_col]), str(row[end_col])
        
        if u in embeddings_dict and v in embeddings_dict:
            try:
                # Structural graph math
                aa = list(nx.adamic_adar_index(G, [(u, v)]))[0][2]
                jc = list(nx.jaccard_coefficient(G, [(u, v)]))[0][2]
            except:
                aa, jc = 0, 0
            
            # Combine vectors
            combined_vec = embeddings_dict[u] * embeddings_dict[v]
            features.append(np.hstack([[aa, jc], combined_vec]))
            targets.append(row['is_high_consensus'])

    return np.array(features), np.array(targets)

def train_model():
    data_result = load_data()
    if not data_result: return
    edges_df, start_col, end_col = data_result
    
    embeddings_dict = load_embeddings()
    if not embeddings_dict: return

    print("Step 4: Building Graph Network...")
    edges_df[start_col] = edges_df[start_col].astype(str)
    edges_df[end_col] = edges_df[end_col].astype(str)
    
    G = nx.from_pandas_edgelist(edges_df, start_col, end_col)
    
    X, y = feature_engineering(G, edges_df, embeddings_dict, start_col, end_col)
    
    if len(X) == 0:
        print("Error: No data matched between files. Your IDs might not match.")
        return

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    
    print(f"Step 5: Training AI Model on {len(X_train)} links...")
    model = RandomForestClassifier(n_estimators=100, max_depth=12, n_jobs=-1, random_state=42)
    model.fit(X_train, y_train)
    
    y_pred = model.predict(X_test)
    print("\n" + "="*40)
    print("SCIENCE FAIR AI MODEL RESULTS")
    print("="*40)
    print(classification_report(y_test, y_pred))
    
    joblib.dump(model, "consensus_predictor_model.pkl")
    print("\nAI Model saved! Ready for your science fair presentation.")

if __name__ == "__main__":
    train_model()