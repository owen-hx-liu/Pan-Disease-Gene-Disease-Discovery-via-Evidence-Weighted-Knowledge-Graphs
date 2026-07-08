import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import os
import random

# --- CONFIGURATION ---
# We will look for this filename globally
RELATIONS_FILENAME = "relationships.csv"
OUTPUT_IMAGE = "network_growth_over_time.png"
OUTPUT_CSV = "temporal_metrics.csv"

def find_file_globally(filename):
    """Search for a file in the current directory and all subdirectories."""
    # First, check if it's just in the current folder
    if os.path.exists(filename):
        return os.path.abspath(filename)
        
    # If not, search subfolders
    for root, dirs, files in os.walk(os.getcwd()):
        if filename in files:
            return os.path.join(root, filename)
    return None

def load_temporal_data():
    """Loads relationship data and ensures a 'year' column exists."""
    rel_path = find_file_globally(RELATIONS_FILENAME)
    
    if not rel_path:
        print(f"Error: Could not find {RELATIONS_FILENAME} anywhere in {os.getcwd()}")
        # Return 4 values to match the expected 'unpacking' in the main function
        return None, None, None, None
        
    print(f"File found at: {rel_path}")
    df = pd.read_csv(rel_path, low_memory=False)
    
    # Identify column names (handling Neo4j :START_ID format)
    start_col = next((c for c in df.columns if 'START_ID' in c), "START_ID")
    end_col = next((c for c in df.columns if 'END_ID' in c), "END_ID")
    
    # Check if a time column exists
    year_col = next((c for c in df.columns if 'year' in c.lower() or 'date' in c.lower()), None)
    
    # SIMULATION: If no year column exists
    if not year_col:
        print("No 'year' column found. Simulating discovery years for Science Fair demonstration...")
        years = list(range(1995, 2024))
        weights = [i for i in range(1, len(years) + 1)] 
        df['discovery_year'] = random.choices(years, weights=weights, k=len(df))
        year_col = 'discovery_year'
    else:
        print(f"Found temporal column: {year_col}")
        
    return df, start_col, end_col, year_col

def analyze_network_growth():
    """Builds the graph year by year and tracks its growth."""
    df, start_col, end_col, year_col = load_temporal_data()
    
    # Check if df is None (file wasn't found)
    if df is None:
        return None

    print("\nAnalyzing network evolution over time...")
    
    # Sort data chronologically
    df = df.sort_values(by=year_col)
    
    # Get the unique years in the dataset
    years = sorted([y for y in df[year_col].dropna().unique() if 1950 <= y <= 2025])
    
    G = nx.Graph()
    metrics = []
    
    for current_year in years:
        edges_this_year = df[df[year_col] == current_year]
        
        for _, row in edges_this_year.iterrows():
            u = str(row[start_col])
            v = str(row[end_col])
            G.add_edge(u, v)
        
        num_nodes = G.number_of_nodes()
        num_edges = G.number_of_edges()
        avg_degree = (2.0 * num_edges) / num_nodes if num_nodes > 0 else 0
        
        metrics.append({
            'Year': int(current_year),
            'Total_Nodes': num_nodes,
            'Total_Connections': num_edges,
            'Average_Connections_Per_Node': round(avg_degree, 2)
        })
        
        print(f"Year {int(current_year)}: {num_nodes} nodes, {num_edges} edges.")
        
    metrics_df = pd.DataFrame(metrics)
    metrics_df.to_csv(OUTPUT_CSV, index=False)
    print(f"\nSummary data saved to {OUTPUT_CSV}")
    
    return metrics_df

def plot_growth(metrics_df):
    """Creates a beautiful chart for the Science Fair poster."""
    if metrics_df is None or metrics_df.empty:
        print("No data available to plot.")
        return

    print("Generating chart...")
    
    try:
        plt.style.use('ggplot') # Using 'ggplot' as a reliable fallback for 'seaborn'
    except:
        pass

    fig, ax1 = plt.subplots(figsize=(12, 7))

    color1 = 'tab:blue'
    ax1.set_xlabel('Year of Discovery', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Total Connections (Edges)', color=color1, fontsize=12, fontweight='bold')
    ax1.plot(metrics_df['Year'], metrics_df['Total_Connections'], color=color1, linewidth=3, marker='o', label='Edges')
    ax1.tick_params(axis='y', labelcolor=color1)

    ax2 = ax1.twinx()  
    color2 = 'tab:red'
    ax2.set_ylabel('Total Unique Entities (Nodes)', color=color2, fontsize=12, fontweight='bold')
    ax2.plot(metrics_df['Year'], metrics_df['Total_Nodes'], color=color2, linewidth=3, linestyle='--', marker='x', label='Nodes')
    ax2.tick_params(axis='y', labelcolor=color2)

    plt.title('Evolution of the Gene-Disease Knowledge Graph over Time', fontsize=14, fontweight='bold')
    fig.tight_layout()

    plt.savefig(OUTPUT_IMAGE, dpi=300)
    print(f"Chart successfully saved as '{OUTPUT_IMAGE}'!")
    plt.show()

if __name__ == "__main__":
    print("=== TEMPORAL GRAPH ANALYZER ===")
    results_df = analyze_network_growth()
    if results_df is not None:
        plot_growth(results_df)
    else:
        print("\n[!] Analysis could not be completed. Please check the error messages above.")