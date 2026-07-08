import time, numpy as np
C="data/processed/benchmark/cache"
t0=time.time(); SEED=42; rng=np.random.default_rng(SEED)
head=np.load(f"{C}/head.npy"); tail=np.load(f"{C}/tail.npy")
N=int(max(head.max(),tail.max()))+1
a=np.minimum(head,tail).astype(np.int64); b=np.maximum(head,tail).astype(np.int64)
loop=a==b; n_loops=int(loop.sum())
a=a[~loop]; b=b[~loop]
key=a*N+b
ukey=np.unique(key)                      # unique undirected simple edges
print(f"N={N} raw_dir={len(head)} self_loops={n_loops} unique_undirected={len(ukey)} t={time.time()-t0:.1f}s",flush=True)
# split
TEST=4000
perm=rng.permutation(len(ukey))
test_idx=perm[:TEST]; train_idx=perm[TEST:]
test_key=ukey[test_idx]; train_key=ukey[train_idx]
test_u=(test_key//N).astype(np.int32); test_v=(test_key%N).astype(np.int32)
tr_u=(train_key//N).astype(np.int32); tr_v=(train_key%N).astype(np.int32)
# CSR from TRAIN edges (both directions), built by sorting
src=np.concatenate([tr_u,tr_v]); dst=np.concatenate([tr_v,tr_u])
order=np.argsort(src,kind='stable'); src=src[order]; dst=dst[order]
indptr=np.zeros(N+1,dtype=np.int64)
np.add.at(indptr, src+1, 1)
indptr=np.cumsum(indptr)
indices=dst.astype(np.int32)
deg=np.diff(indptr)
assert int(deg.sum())==len(src), "CSR degree mismatch"
# sort neighbors within each row for fast intersection
for i in range(N):
    s,e=indptr[i],indptr[i+1]
    if e>s+1: indices[s:e]=np.sort(indices[s:e])
np.save(f"{C}/indptr.npy",indptr); np.save(f"{C}/indices.npy",indices); np.save(f"{C}/deg.npy",deg)
np.save(f"{C}/test_u.npy",test_u); np.save(f"{C}/test_v.npy",test_v)
np.save(f"{C}/ukey.npy",ukey); np.save(f"{C}/Nnodes.npy",np.array([N]))
# negatives for AUROC/AUPRC: same count as test, filtered (not in full undirected set, not self loop)
def sample_neg(m):
    out_u=np.empty(m,dtype=np.int32); out_v=np.empty(m,dtype=np.int32); got=0
    while got<m:
        cu=rng.integers(0,N,size=2*(m-got)); cv=rng.integers(0,N,size=2*(m-got))
        aa=np.minimum(cu,cv).astype(np.int64); bb=np.maximum(cu,cv).astype(np.int64)
        ok=aa!=bb
        k=aa*N+bb
        pos=np.searchsorted(ukey,k)
        pos=np.clip(pos,0,len(ukey)-1)
        exists=ukey[pos]==k
        ok&=~exists
        cu=cu[ok][:m-got]; cv=cv[ok][:m-got]
        out_u[got:got+len(cu)]=cu; out_v[got:got+len(cv)]=cv; got+=len(cu)
    return out_u,out_v
neg_u,neg_v=sample_neg(TEST)
np.save(f"{C}/neg_u.npy",neg_u); np.save(f"{C}/neg_v.npy",neg_v)
print(f"split done test={TEST} train_edges={len(train_key)} mean_deg={deg.mean():.2f} max_deg={deg.max()} t={time.time()-t0:.1f}s",flush=True)
