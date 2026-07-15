import os

from neo4j import GraphDatabase
import pandas as pd

# Neo4j connection read from the environment -- no secrets in source.
#   NEO4J_URI       (default bolt://localhost:7687)
#   NEO4J_USER      (default neo4j)
#   NEO4J_PASSWORD  (required -- export before running)
URI = os.environ.get("NEO4J_URI", "bolt://localhost:7687")
_PASSWORD = os.environ.get("NEO4J_PASSWORD")
if not _PASSWORD:
    raise RuntimeError(
        "NEO4J_PASSWORD is not set. Export your Neo4j password before running, e.g. "
        "setx NEO4J_PASSWORD <password> (Windows) or export NEO4J_PASSWORD=<password> (Linux/macOS)."
    )
AUTH = (os.environ.get("NEO4J_USER", "neo4j"), _PASSWORD)

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