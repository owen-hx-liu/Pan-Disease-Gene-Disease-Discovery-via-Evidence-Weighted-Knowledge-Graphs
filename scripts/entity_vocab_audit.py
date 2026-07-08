import pandas as pd

df = pd.read_csv("data/processed/entity_vocab.csv")

print("\n=== BASIC INFO ===")
print(df.shape)

print("\n=== TYPE COUNTS ===")
print(df["type"].value_counts().head(20))

print("\n=== IDS WITHOUT NAMESPACE (no colon) ===")
no_ns = df[~df["id"].astype(str).str.contains(":")]
print(no_ns.head())
print("Count:", len(no_ns))

print("\n=== LOWERCASE TYPES (bad) ===")
print(df[df["type"].str.lower() != df["type"]].head())

print("\n=== SAMPLE IDS ===")
print(df["id"].sample(20).tolist())
