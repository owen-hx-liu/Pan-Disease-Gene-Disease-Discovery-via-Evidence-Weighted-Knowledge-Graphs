import pandas as pd
from neo4j import GraphDatabase

from config import neo4j_credentials, NODES_CSV, RELATIONSHIPS_CSV

NEO4J_URI, _AUTH = neo4j_credentials()  # no secrets in source; set NEO4J_PASSWORD env var

NODES_PATH = str(NODES_CSV)
RELS_PATH = str(RELATIONSHIPS_CSV)

driver = GraphDatabase.driver(NEO4J_URI, auth=_AUTH)

nodes_df = pd.read_csv(NODES_PATH)
rels_df = pd.read_csv(RELS_PATH)

print("NODES columns:", nodes_df.columns.tolist())
print("RELATIONSHIPS columns:", rels_df.columns.tolist())

with driver.session() as session:
    # Upload nodes
    for _, row in nodes_df.iterrows():
        props = row.to_dict()
        label = props.pop("label", "Node")
        query = f"MERGE (n:{label} {{name: $name}}) SET n += $props"
        session.run(query, name=props.get("name", str(props)), props=props)
    print(f"✅ Uploaded {len(nodes_df)} nodes")

    # Upload relationships
    for _, row in rels_df.iterrows():
        props = row.to_dict()
        source = props.pop("source", None)
        target = props.pop("target", None)
        rel_type = props.pop("type", "INTERACTS_WITH")
        query = f"""
            MATCH (a {{name: $source}})
            MATCH (b {{name: $target}})
            MERGE (a)-[r:{rel_type}]->(b)
            SET r += $props
        """
        session.run(query, source=source, target=target, props=props)
    print(f"✅ Uploaded {len(rels_df)} relationships")

driver.close()
print("🎉 Database restored!")