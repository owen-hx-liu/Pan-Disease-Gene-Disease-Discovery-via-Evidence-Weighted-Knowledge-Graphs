from neo4j import GraphDatabase
import pandas as pd

# 1. Update these with your Neo4j info
URI = "bolt://localhost:7687"
AUTH = ("neo4j", "La1nos#b") # Change to your actual password

def export_to_csv():
    driver = GraphDatabase.driver(URI, auth=AUTH)
    query = "MATCH (n) WHERE n.embedding IS NOT NULL RETURN n.node_id AS node_id, n.embedding AS embedding"
    
    print("Connecting to Neo4j and fetching embeddings... this may take a minute.")
    
    with driver.session() as session:
        result = session.run(query)
        # Convert the result directly into a Pandas DataFrame
        df = pd.DataFrame([dict(record) for record in result])
        
    print(f"Successfully fetched {len(df)} nodes.")
    
    # Save to CSV
    df.to_csv("node2vec_embeddings.csv", index=False)
    print("File saved: node2vec_embeddings.csv")
    
    driver.close()

if __name__ == "__main__":
    export_to_csv()