import os, time, numpy as np, pandas as pd
SRC = "data/processed/edges_clean.csv"
OUT = "data/processed/benchmark/cache"
os.makedirs(OUT, exist_ok=True)
t0 = time.time()
df = pd.read_csv(SRC, usecols=["source_id","relation","target_id"], dtype=str)
print(f"read {len(df)} rows t={time.time()-t0:.1f}s", flush=True)
# factorize nodes jointly over head+tail
both = pd.concat([df["source_id"], df["target_id"]], ignore_index=True)
codes, uniques = pd.factorize(both, sort=False)
n = len(df)
head = codes[:n].astype(np.int32); tail = codes[n:].astype(np.int32)
del both, codes
rcodes, runiques = pd.factorize(df["relation"], sort=False)
rel = rcodes.astype(np.int16)
np.save(f"{OUT}/head.npy",head); np.save(f"{OUT}/tail.npy",tail); np.save(f"{OUT}/rel.npy",rel)
np.save(f"{OUT}/node_labels.npy", np.asarray(uniques, dtype=object))
np.save(f"{OUT}/rel_labels.npy", np.asarray(runiques, dtype=object))
print(f"DONE n_edges={n} n_nodes={len(uniques)} n_rels={len(runiques)} t={time.time()-t0:.1f}s", flush=True)
