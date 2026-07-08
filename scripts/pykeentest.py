import pandas as pd

df = pd.read_csv("data/processed/transetableuseforpykeen.csv")

print(df.head())
print(df.columns)
print(df.shape)

# 🔥 CRITICAL CHECK
raw = pd.read_csv("data/processed/transetableuseforpykeen.csv", header=None)
print(raw.head())
print(raw.shape)