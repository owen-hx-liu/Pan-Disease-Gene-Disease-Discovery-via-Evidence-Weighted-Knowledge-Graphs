import streamlit as st
import joblib
import pandas as pd
import numpy as np
import networkx as nx
import ast
import os

# --- PAGE CONFIGURATION ---
st.set_page_config(
    page_title="Bio-Link Analysis Portal",
    layout="wide"
)

# --- UI STYLING ---
st.markdown("""
    <style>
    /* Dark professional background */
    .stApp {
        background-color: #0f172a;
        color: #f8fafc;
    }
    
    /* Metric Score Styling - FORCE TO WHITE */
    [data-testid="stMetricValue"] {
        color: #ffffff !important;
        font-weight: 700;
        font-size: 3rem !important;
    }
    
    /* Label for the metric */
    [data-testid="stMetricLabel"] {
        color: #94a3b8 !important;
    }

    /* Sidebar Styling */
    section[data-testid="stSidebar"] {
        background-color: #1e293b;
        border-right: 1px solid #334155;
    }
    
    /* Input field styling */
    .stTextInput>div>div>input {
        background-color: #1e293b;
        color: white;
        border: 1px solid #334155;
    }
    
    /* Table styling for dark mode */
    .stTable {
        background-color: #1e293b;
        color: white;
    }
    
    /* Button styling */
    .stButton>button {
        width: 100%;
        background-color: #3b82f6;
        color: white;
        border: none;
        padding: 0.75rem;
        font-weight: 600;
    }
    .stButton>button:hover {
        background-color: #2563eb;
        border: none;
        color: white;
    }
    </style>
    """, unsafe_allow_html=True)

# --- CONFIG & PATHS ---
MODEL_PATH = "data/processed/consensus_predictor_model.pkl"
EMBEDDINGS_PATH = "data/processed/node2vec_embeddings.csv"
RELATIONS_PATH = "data/processed/relationships.csv"

# --- DATA LOADING ---
@st.cache_resource
def load_resources(source_col="START_ID", target_col="END_ID"):
    model, embeddings_dict, G = None, None, None
    
    if os.path.exists(MODEL_PATH):
        try:
            model = joblib.load(MODEL_PATH)
        except:
            pass

    if os.path.exists(EMBEDDINGS_PATH):
        try:
            emb_df = pd.read_csv(EMBEDDINGS_PATH)
            embeddings_dict = {
                str(row['node_id']): (np.array(ast.literal_eval(row['embedding'])) if isinstance(row['embedding'], str) else np.array(row['embedding'])) 
                for _, row in emb_df.iterrows()
            }
        except:
            pass

    if os.path.exists(RELATIONS_PATH):
        try:
            edges_df = pd.read_csv(RELATIONS_PATH, low_memory=False)
            edges_df.columns = [c.strip() for c in edges_df.columns]
            actual_src = next((c for c in edges_df.columns if source_col in c), edges_df.columns[0])
            actual_tgt = next((c for c in edges_df.columns if target_col in c), edges_df.columns[1])
            edges_df[actual_src] = edges_df[actual_src].astype(str)
            edges_df[actual_tgt] = edges_df[actual_tgt].astype(str)
            G = nx.from_pandas_edgelist(edges_df, actual_src, actual_tgt)
        except:
            pass
            
    return model, embeddings_dict, G

# --- SIDEBAR ---
with st.sidebar:
    st.title("Settings")
    st.markdown("---")
    src_input = st.text_input("Source ID Column", value="START_ID")
    tgt_input = st.text_input("Target ID Column", value="END_ID")
    st.info("System checking status of CSV and Model files.")

# --- MAIN UI ---
st.title("Bio-Link Prediction Interface")
st.markdown("Analyze biological associations using graph topology and node embeddings.")

model, embeddings, G = load_resources(src_input, tgt_input)

if G is not None:
    st.caption(f"Network synchronized: {G.number_of_nodes()} unique entities indexed.")

st.divider()

# Input Section
col_a, col_b = st.columns(2)
with col_a:
    node_a = st.text_input("Source Identifier (Node A)", placeholder="Enter ID")
with col_b:
    node_b = st.text_input("Target Identifier (Node B)", placeholder="Enter ID")

if st.button("Generate Link Prediction"):
    if not node_a or not node_b:
        st.warning("Please provide both identifiers.")
    elif G is None or embeddings is None or model is None:
        st.error("Resources missing. Verify model and CSV files exist in the script directory.")
    elif node_a not in embeddings or node_b not in embeddings:
        st.error("Node ID Error: One or both nodes were not found in the embedding dataset.")
    else:
        try:
            # 1. Feature Extraction (Adamic-Adar logic preserved)
            if node_a in G and node_b in G:
                aa_gen = nx.adamic_adar_index(G, [(node_a, node_b)])
                aa_score = list(aa_gen)[0][2]
                jc_gen = nx.jaccard_coefficient(G, [(node_a, node_b)])
                jc_score = list(jc_gen)[0][2]
            else:
                aa_score, jc_score = 0.0, 0.0
            
            # 2. Embedding combination (Hadamard product)
            combined_vec = embeddings[node_a] * embeddings[node_b]
            features = np.hstack([[aa_score, jc_score], combined_vec]).reshape(1, -1)
            
            # 3. Model Inference
            prob = model.predict_proba(features)[0][1]
            pred = model.predict(features)[0]
            
            # 4. Results Display
            st.subheader("Analysis Breakdown")
            res_left, res_right = st.columns([1, 2])
            
            with res_left:
                st.metric("Link Probability", f"{prob:.4f}")
                if pred == 1:
                    st.info("High Confidence Association Detected")
                else:
                    st.warning("Low Confidence Association Detected")
            
            with res_right:
                st.markdown("**Calculated Features**")
                stats_df = pd.DataFrame({
                    "Feature Name": ["Adamic-Adar", "Jaccard Index", "Embedding Dim"],
                    "Value": [f"{aa_score:.6f}", f"{jc_score:.6f}", f"{len(combined_vec)}"]
                })
                st.table(stats_df)
                
        except Exception as e:
            st.error(f"Computation error: {str(e)}")

st.markdown("---")
st.caption("Science Fair Research Tool | Analysis System")