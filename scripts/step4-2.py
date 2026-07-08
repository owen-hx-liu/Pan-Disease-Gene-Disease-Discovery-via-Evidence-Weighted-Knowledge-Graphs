#!/usr/bin/env python3

import pandas as pd
import json
from pathlib import Path
from collections import defaultdict
from tqdm import tqdm

ENTITY_VOCAB_FILE = Path("data/processed/entity_vocab.csv")
ID_PRIORITY_FILE = Path("scripts/id_priority.json")
OUTPUT_MAP_FILE = Path("data/processed/entity_id_map.csv")

# ================= LOAD =================

print("Loading entity vocab...")
entity_df = pd.read_csv(ENTITY_VOCAB_FILE)

print("Loading ID priority rules...")
with open(ID_PRIORITY_FILE) as f:
    ID_PRIORITY = json.load(f)

# ================= HELPERS =================

def namespace(x):
    if isinstance(x,str) and ":" in x:
        return x.split(":",1)[0].upper()
    return "UNKNOWN"

# Clean columns
entity_df = entity_df.dropna(subset=["id"])
entity_df["id"] = entity_df["id"].astype(str).str.strip()
entity_df["name"] = entity_df["name"].fillna("").astype(str).str.lower().str.strip()
entity_df["type"] = entity_df["type"].fillna("Unknown")

# ================= UNION FIND =================

class UF:
    def __init__(self):
        self.p={}
    def find(self,x):
        if self.p.get(x,x)!=x:
            self.p[x]=self.find(self.p[x])
        return self.p.get(x,x)
    def union(self,a,b):
        ra,rb=self.find(a),self.find(b)
        if ra!=rb:
            self.p[rb]=ra

uf=UF()

# ================= GROUP BY (TYPE + NAME) =================

print("\nGrouping by (type,name)...")

name_groups=defaultdict(list)

for _,r in tqdm(entity_df.iterrows(),
                total=len(entity_df),
                desc="Grouping"):
    key=(r["type"], r["name"])
    name_groups[key].append(r["id"])

# ================= MERGE WITHIN NAME GROUP =================

print("\nMerging IDs within same-name groups...")

for (etype,name),ids in tqdm(name_groups.items(),
                             desc="Processing groups"):

    if len(ids)<=1:
        continue

    priorities=ID_PRIORITY.get(etype,[])

    # group by namespace
    ns_map=defaultdict(list)
    for i in ids:
        ns_map[namespace(i)].append(i)

    # build priority-ordered list
    ordered=[]
    for p in priorities:
        ordered+=ns_map.get(p,[])

    # include remaining namespaces
    for ns,vals in ns_map.items():
        if ns not in priorities:
            ordered+=vals

    base=ordered[0]
    for other in ordered[1:]:
        uf.union(base,other)

# ================= BUILD CLUSTERS =================

print("\nBuilding clusters...")

clusters=defaultdict(list)

for rid in tqdm(entity_df["id"], desc="Clustering"):
    clusters[uf.find(rid)].append(rid)

print("Total clusters:",len(clusters))

# ================= CANONICAL SELECTION =================

def pick_canonical(ids,etype):
    priorities=ID_PRIORITY.get(etype,[])
    ns_map=defaultdict(list)

    for i in ids:
        ns_map[namespace(i)].append(i)

    for p in priorities:
        if ns_map.get(p):
            return ns_map[p][0]

    return min(ids,key=len)

print("\nSelecting canonical IDs...")

rows=[]

for root,members in tqdm(clusters.items(),
                         desc="Canonical mapping"):

    etype = entity_df.loc[
        entity_df["id"]==members[0],"type"
    ].values[0]

    canon=pick_canonical(members,etype)

    for m in members:
        rows.append({
            "raw_id":m,
            "canonical_id":canon,
            "entity_type":etype
        })

# ================= SAVE =================

print("\nSaving...")

out=pd.DataFrame(rows)

OUTPUT_MAP_FILE.parent.mkdir(parents=True,exist_ok=True)
out.to_csv(OUTPUT_MAP_FILE,index=False)

print("\n✅ DONE")
print("Total mappings:",len(out))
