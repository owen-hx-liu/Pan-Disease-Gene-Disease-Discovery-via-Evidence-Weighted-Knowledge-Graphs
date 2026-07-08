import pandas as pd
from neo4j import GraphDatabase
import os

from config import neo4j_credentials  # reads NEO4J_* from the environment

# --- CONFIGURATION ---
URI, AUTH = neo4j_credentials()  # no secrets in source; set NEO4J_PASSWORD env var
TSV_PATH = "data/drug.target.interaction.tsv" # Ensure this file is in your folder

class DrugInjector:
    def __init__(self, uri, auth):
        self.driver = GraphDatabase.driver(uri, auth=auth)

    def close(self):
        self.driver.close()

    def inject_drug_central_data(self, file_path):
        """Reads the DrugCentral TSV and injects it into Neo4j."""
        if not os.path.exists(file_path):
            print(f"Error: Could not find {file_path}. Please place the TSV in the project folder.")
            return

        print(f"Reading DrugCentral data from {file_path}...")
        # DrugCentral TSVs are tab-separated
        df = pd.read_csv(file_path, sep='\t')
        
        # We need to map DrugCentral columns to your Graph nodes
        # Usually: 'DRUG_NAME' and 'GENE_SYMBOL' or 'ENTREZ_ID'
        # Let's use GENE_SYMBOL as it's common in bio-graphs
        drug_col = 'DRUG_NAME'
        gene_col = 'GENE_SYMBOL' 
        
        if drug_col not in df.columns or gene_col not in df.columns:
            print("Columns not matching. Found columns:", df.columns.tolist())
            print("Attempting to guess columns based on content...")
            # Fallback guessing logic
            drug_col = [c for c in df.columns if 'drug' in c.lower()][0]
            gene_col = [c for c in df.columns if 'gene' in c.lower() or 'symbol' in c.lower()][0]

        # Filter out rows with missing data
        df = df[[drug_col, gene_col]].dropna().drop_duplicates()
        print(f"Processing {len(df)} drug-target interactions...")

        with self.driver.session() as session:
            # 1. Inject the Drug-to-Gene links
            print("Injecting links into Neo4j (Step 1 of 2)...")
            count = 0
            for _, row in df.head(5000).iterrows(): # Testing with first 5000 for speed
                success = session.execute_write(self._create_drug_connection, row[drug_col], str(row[gene_col]))
                if success: count += 1
            
            print(f"Successfully linked {count} drugs to genes in your graph.")

            # 2. Discover the 'Discovery Delta' (Drug -> Disease)
            print("Identifying Potential Treatments (Step 2 of 2)...")
            session.execute_write(self._discover_potential_treatments)

    @staticmethod
    def _create_drug_connection(tx, drug_name, gene_id):
        # We try to match the Gene by node_id OR prefix/name since bio-data varies
        query = """
        MERGE (d:Drug {name: $drug_name})
        WITH d
        MATCH (g:Gene) 
        WHERE g.node_id = $gene_id OR g.name = $gene_id
        MERGE (d)-[r:TARGETS]->(g)
        RETURN count(r) as created
        """
        result = tx.run(query, drug_name=drug_name, gene_id=gene_id)
        record = result.single()
        return record['created'] > 0 if record else False

    @staticmethod
    def _discover_potential_treatments(tx):
        """
        Calculates which drugs might treat which diseases based on 
        the strength of the gene-disease connection in your graph.
        """
        query = """
        MATCH (d:Drug)-[:TARGETS]->(g:Gene)-[r]->(dis:Disease)
        WHERE NOT type(r) = 'TARGETS'
        WITH d, dis, SUM(toInteger(r.weight)) as repurpose_score
        WHERE repurpose_score > 2  // Only keep significant connections
        MERGE (d)-[pt:POTENTIAL_TREATMENT]->(dis)
        SET pt.confidence = repurpose_score,
            pt.method = 'AI Inferred via Knowledge Graph'
        RETURN count(pt) as new_links
        """
        result = tx.run(query)
        record = result.single()
        print(f"Discovery complete: {record['new_links'] if record else 0} POTENTIAL_TREATMENT links created.")

if __name__ == "__main__":
    print("=== DRUGCENTRAL KNOWLEDGE INJECTOR ===")
    injector = DrugInjector(URI, AUTH)
    try:
        injector.inject_drug_central_data(TSV_PATH)
    finally:
        injector.close()