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
    for root, dirs, files in os.walk(os.getcwd()):
        if filename in files:
            return os.path.join(root, filename)
    return None

def generate_negative_samples(G, nodes, num_samples):
    """
    Since your data currently has weights of 1 for everything, 
    we MUST create 'Negative Samples' (fake links) so the AI 
    has something to compare against.
    """
    print(f"Generating {num_samples} negative samples (fake links)...")
    neg_edges = []
    nodes_list = list(nodes)
    existing_edges = set(G.edges())
    
    max_attempts = num_samples * 20
    attempts = 0
    
    while len(neg_edges) < num_samples and attempts < max_attempts:
        u, v = random.sample(nodes_list, 2)
        if u != v and (u, v) not in existing_edges and (v, u) not in existing_edges:
            neg_edges.append((u, v))
        attempts += 1
    return neg_edges

def load_data():
    rel_path = find_file_globally("relationships.csv")
    if not rel_path: 
        print("Error: relationships.csv not found.")
        return None
    
    edges_df = pd.read_csv(rel_path, low_memory=False)
    
    start_col = next((c for c in edges_df.columns if 'START_ID' in c), None)
    end_col = next((c for c in edges_df.columns if 'END_ID' in c), None)
    
    if not start_col or not end_col:
        print("Error: Could not identify ID columns in relationships.csv")
        return None

    edges_df['is_high_consensus'] = 1
        
    return edges_df, start_col, end_col

def load_embeddings():
    emb_path = find_file_globally("node2vec_embeddings.csv")
    if not emb_path: 
        print("Error: node2vec_embeddings.csv not found.")
        return None
        
    emb_df = pd.read_csv(emb_path)
    embeddings_dict = {}
    
    print(f"Processing {len(emb_df)} node embeddings...")
    for _, row in emb_df.iterrows():
        node_id = str(row['node_id'])
        raw_emb = row['embedding']
        
        try:
            if isinstance(raw_emb, str):
                # Converts "[0.1, 0.2...]" string to a numpy array
                vector = np.array(ast.literal_eval(raw_emb.strip()))
                embeddings_dict[node_id] = vector
            else:
                embeddings_dict[node_id] = np.array(raw_emb)
        except Exception as e:
            continue
            
    return embeddings_dict

def feature_engineering(G, pair_list, embeddings_dict, label):
    """Turns pairs of nodes into a single row of data for the AI."""
    features = []
    labels = []
    for u, v in pair_list:
        u, v = str(u), str(v)
        if u in embeddings_dict and v in embeddings_dict:
            try:
                aa = list(nx.adamic_adar_index(G, [(u, v)]))[0][2]
                
                jc = list(nx.jaccard_coefficient(G, [(u, v)]))[0][2]
            except:
                aa, jc = 0.0, 0.0
            
            combined_vec = embeddings_dict[u] * embeddings_dict[v]
            
            features.append(np.hstack([[aa, jc], combined_vec]))
            labels.append(label)
            
    return features, labels

def train_model():
    data_res = load_data()
    if not data_res: return
    edges_df, start_col, end_col = data_res
    
    embeddings_dict = load_embeddings()
    if not embeddings_dict: return

    print("Building Graph Network...")
    edges_df[start_col] = edges_df[start_col].astype(str)
    edges_df[end_col] = edges_df[end_col].astype(str)
    G = nx.from_pandas_edgelist(edges_df, start_col, end_col)
    
    sample_size = min(25000, len(edges_df))
    pos_edges = edges_df.sample(n=sample_size, random_state=42)
    pos_pairs = list(zip(pos_edges[start_col], pos_edges[end_col]))
    
    neg_pairs = generate_negative_samples(G, embeddings_dict.keys(), sample_size)
    
    print("Calculating math features for AI training...")
    pos_X, pos_y = feature_engineering(G, pos_pairs, embeddings_dict, 1)
    neg_X, neg_y = feature_engineering(G, neg_pairs, embeddings_dict, 0)
    
    X = np.array(pos_X + neg_X)
    y = np.array(pos_y + neg_y)

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    
    print(f"Training Random Forest on {len(X_train)} samples...")
    model = RandomForestClassifier(n_estimators=100, max_depth=10, n_jobs=-1, random_state=42, class_weight='balanced')
    model.fit(X_train, y_train)
    
    y_pred = model.predict(X_test)
    print("\n" + "="*40)
    print("SCIENCE FAIR AI MODEL RESULTS")
    print("="*40)
    print(classification_report(y_test, y_pred))
    
    importances = model.feature_importances_
    print("\nProject Insights:")
    print(f"1. Importance of 'Local' Neighbors (Jaccard/Adamic): {importances[0] + importances[1]:.4f}")
    print(f"2. Importance of 'Global' Position (Node2Vec): {np.sum(importances[2:]):.4f}")

    joblib.dump(model, "consensus_predictor_model.pkl")
    print("\nModel successfully saved as 'consensus_predictor_model.pkl'.")

if __name__ == "__main__":
    train_model()