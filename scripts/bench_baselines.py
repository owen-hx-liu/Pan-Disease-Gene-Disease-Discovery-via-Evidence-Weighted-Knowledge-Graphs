import time, json, numpy as np
C="data/processed/benchmark/cache"
t0=time.time(); rng=np.random.default_rng(7)
indptr=np.load(f"{C}/indptr.npy"); indices=np.load(f"{C}/indices.npy"); deg=np.load(f"{C}/deg.npy")
test_u=np.load(f"{C}/test_u.npy"); test_v=np.load(f"{C}/test_v.npy")
neg_u=np.load(f"{C}/neg_u.npy"); neg_v=np.load(f"{C}/neg_v.npy")
ukey=np.load(f"{C}/ukey.npy"); N=int(np.load(f"{C}/Nnodes.npy")[0])
# AA weight per node from train degree
with np.errstate(divide='ignore'):
    aaw=np.where(deg>=2, 1.0/np.log(deg.astype(np.float64)), 0.0)

def nbrs(u): return indices[indptr[u]:indptr[u+1]]
def cn_aa(u,v):
    c=np.intersect1d(nbrs(u),nbrs(v),assume_unique=True)
    return len(c), float(aaw[c].sum())

def score_arr(us,vs):
    cn=np.empty(len(us)); aa=np.empty(len(us))
    for i in range(len(us)):
        cn[i],aa[i]=cn_aa(int(us[i]),int(vs[i]))
    return cn,aa

# ---- classification scores (AUROC/AUPRC): pos vs neg ----
pos_cn,pos_aa=score_arr(test_u,test_v)
neg_cn,neg_aa=score_arr(neg_u,neg_v)
print(f"scored class sets t={time.time()-t0:.1f}s",flush=True)

def rankdata_avg(x):
    order=np.argsort(x,kind='stable'); r=np.empty(len(x)); sx=x[order]
    i=0
    while i<len(x):
        j=i
        while j+1<len(x) and sx[j+1]==sx[i]: j+=1
        r[order[i:j+1]]=(i+j)/2.0+1.0; i=j+1
    return r
def auroc(pos,neg):
    s=np.concatenate([pos,neg]); y=np.concatenate([np.ones(len(pos)),np.zeros(len(neg))])
    r=rankdata_avg(s); np_,nn=len(pos),len(neg)
    return (r[y==1].sum()-np_*(np_+1)/2.0)/(np_*nn)
def auprc(pos,neg):
    s=np.concatenate([pos,neg]); y=np.concatenate([np.ones(len(pos)),np.zeros(len(neg))])
    order=np.argsort(-s,kind='stable'); y=y[order]
    tp=np.cumsum(y); fp=np.cumsum(1-y)
    prec=tp/(tp+fp); rec=tp/y.sum()
    rec=np.concatenate([[0],rec]); prec=np.concatenate([[1],prec])
    return float(np.sum((rec[1:]-rec[:-1])*prec[1:]))

# ---- sanity check metrics on a known toy case ----
# perfect separation -> AUROC 1, AUPRC 1 ; identical -> AUROC 0.5
assert abs(auroc(np.array([3.,4,5]),np.array([0.,1,2]))-1.0)<1e-9
assert abs(auroc(np.array([1.,1,1]),np.array([1.,1,1]))-0.5)<1e-9
assert abs(auprc(np.array([3.,4,5]),np.array([0.,1,2]))-1.0)<1e-9

rand_pos=rng.random(len(test_u)); rand_neg=rng.random(len(neg_u))
res={"protocol":{"test_edges":int(len(test_u)),"neg_per_pos_class":1,
  "ranking_negatives_per_pos":50,"split":"random held-out (no real temporal data available)","seed":42},
  "classification":{
    "Random":{"AUROC":auroc(rand_pos,rand_neg),"AUPRC":auprc(rand_pos,rand_neg)},
    "CommonNeighbors":{"AUROC":auroc(pos_cn,neg_cn),"AUPRC":auprc(pos_cn,neg_cn)},
    "AdamicAdar":{"AUROC":auroc(pos_aa,neg_aa),"AUPRC":auprc(pos_aa,neg_aa)}}}

# ---- ranking: tail corruption with K sampled filtered negatives ----
K=50
def rank_metrics(score_fn):
    rr=[]; h1=h3=h10=0
    for i in range(len(test_u)):
        u=int(test_u[i]); vt=int(test_v[i])
        # sample K negative tails not real neighbors of u and != u
        cand=rng.integers(0,N,size=K*2)
        a=np.minimum(u,cand).astype(np.int64); b=np.maximum(u,cand).astype(np.int64)
        k=a*N+b; pos=np.clip(np.searchsorted(ukey,k),0,len(ukey)-1)
        good=(ukey[pos]!=k)&(cand!=u)
        cand=cand[good][:K]
        ts=score_fn(u,vt)
        ns=np.array([score_fn(u,int(c)) for c in cand])
        rank=1+int(np.sum(ns>ts))+0.5*int(np.sum(ns==ts))
        rr.append(1.0/rank); h1+=rank<=1; h3+=rank<=3; h10+=rank<=10
    n=len(test_u)
    return {"MRR":float(np.mean(rr)),"Hits@1":h1/n,"Hits@3":h3/n,"Hits@10":h10/n}

def cn_fn(u,v): return float(len(np.intersect1d(nbrs(u),nbrs(v),assume_unique=True)))
def aa_fn(u,v): return float(aaw[np.intersect1d(nbrs(u),nbrs(v),assume_unique=True)].sum())
def rand_fn(u,v): return rng.random()

res["ranking"]={"Random":rank_metrics(rand_fn)}
print(f"rand rank done t={time.time()-t0:.1f}s",flush=True)
res["ranking"]["CommonNeighbors"]=rank_metrics(cn_fn)
print(f"cn rank done t={time.time()-t0:.1f}s",flush=True)
res["ranking"]["AdamicAdar"]=rank_metrics(aa_fn)
print(f"aa rank done t={time.time()-t0:.1f}s",flush=True)
json.dump(res,open(f"{C}/baseline_results.json","w"),indent=2)
print(json.dumps(res,indent=2))
print(f"ALL DONE t={time.time()-t0:.1f}s")
