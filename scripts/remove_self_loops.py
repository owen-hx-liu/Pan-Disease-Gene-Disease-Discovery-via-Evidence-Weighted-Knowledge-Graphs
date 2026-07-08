import pandas as pd

EDGES_FILE = "data/processed/edges_deduplicated.csv"
OUTPUT_FILE = "data/processed/edges_clean.csv"

INVALID_SELF_LOOP_RELATIONS = {
    "biolink:orthologous_to",
    "biolink:interacts_with",
    "biolink:genetically_interacts_with",
}

print("Loading edges...")
df = pd.read_csv(EDGES_FILE)

print("Initial edges:", len(df))

# Normalize relation case
df["relation_lower"] = df["relation"].str.lower()

# Filter
mask = ~(
    (df["source_id"] == df["target_id"]) &
    (df["relation_lower"].isin(INVALID_SELF_LOOP_RELATIONS))
)

clean_df = df[mask].drop(columns=["relation_lower"])

print("Removed self-loops:", len(df) - len(clean_df))
print("Remaining edges:", len(clean_df))

print("Saving cleaned edges...")
clean_df.to_csv(OUTPUT_FILE, index=False)

print("✅ Done")
