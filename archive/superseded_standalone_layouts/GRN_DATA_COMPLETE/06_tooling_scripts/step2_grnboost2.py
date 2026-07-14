#!/usr/bin/env python
"""
STEP 2 - GRNBoost2 co-expression network (complementary observational view), per lineage.
regulators = perturbed targets (present as genes) + canonical T-cell TFs.
targets(labels) = effector program genes + costim receptors.
Sign of each edge = sign of within-lineage Pearson corr(regulator, target) * importance.
LocalCluster capped at 40 workers (shared box). Runs CD8 then CD4 on ONE cluster (sequential).

NOTE: all heavy work (h5ad load, cluster, grnboost2) lives under __main__ so dask worker
subprocesses that RE-IMPORT this module only pick up imports + constants + function defs,
never the 3GB load. (Previous runaway = 40 workers each loading the h5ad on import.)
"""
import os, sys, json, time
os.environ["NUMEXPR_MAX_THREADS"]="64"; os.environ["NUMEXPR_NUM_THREADS"]="8"
os.environ["OMP_NUM_THREADS"]="1"; os.environ["OPENBLAS_NUM_THREADS"]="1"; os.environ["MKL_NUM_THREADS"]="1"
import numpy as np, scipy.sparse as sp, pandas as pd, h5py

def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

BASE="/media/balthasar-lab/RAID4/costim_engager_counterscreen"
F=f"{BASE}/data/schmidt2022_tcell_perturbseq/cd8_grn_input.h5ad"
WORK=f"{BASE}/work/cd8cd4_gof_grn"
SEED=0; NWORKERS=40

EFFECTOR=["IFNG","GZMB","GZMA","GZMK","PRF1","TNF","IL2","FASLG","NKG7","KLRG1","GNLY","CCL4","CCL5","XCL1"]
COSTIM=["CD28","CD27","TNFRSF9","CD2","LTBR"]                 # 5 arms present in screen (drive-scored)
COSTIM12=["CD28","CD27","TNFRSF9","TNFRSF4","TNFRSF18","TNFRSF25",
          "TNFRSF14","TNFRSF8","ICOS","CD226","CD2","CD40"]  # parent's 12-arm plotting panel
COSTIM_PANEL=sorted(set(COSTIM12)|set(COSTIM))               # all costim NODES for the network (13)
TFS=["TBX21","EOMES","PRDM1","TCF7","RUNX3","ID2","ZEB2","BATF","IRF4","FOXO1"]

def read_col(h, path):
    o=h[path]
    if hasattr(o,"keys"):
        cats=o["categories"][:]; codes=o["codes"][:]
        cats=np.array([c.decode() if isinstance(c,bytes) else c for c in cats])
        return np.where(codes>=0, cats[codes.clip(min=0)], None)
    arr=o[:]
    if arr.dtype.kind in ("S","O"):
        arr=np.array([a.decode() if isinstance(a,bytes) else a for a in arr])
    return arr

def signed_corr(Xdense, cols, reg_names, tgt_names):
    name2c={g:i for i,g in enumerate(cols)}
    R=Xdense[:, [name2c[g] for g in reg_names]]
    T=Xdense[:, [name2c[g] for g in tgt_names]]
    Rz=(R-R.mean(0)); Tz=(T-T.mean(0))
    Rs=np.sqrt((Rz**2).sum(0)); Ts=np.sqrt((Tz**2).sum(0))
    num=Rz.T@Tz
    den=np.outer(Rs,Ts); den[den==0]=np.nan
    return num/den

def main():
    os.makedirs(WORK, exist_ok=True)
    log("loading", F)
    with h5py.File(F,"r") as h:
        shape=tuple(h["X"].attrs["shape"])
        X=sp.csr_matrix((h["X/data"][:],h["X/indices"][:],h["X/indptr"][:]),shape=shape)
        sym=read_col(h,"var/sym")
        guide=read_col(h,"obs/guide_target").astype(str)
    gene2col={g:i for i,g in enumerate(sym)}
    log(f"X {X.shape}")

    perturbed=sorted(set(guide)-{"NO-TARGET"})
    # regulators = perturbed targets + canonical TFs + all 12+ costim panel arms (present as genes)
    regs=[]
    for g in perturbed+TFS+COSTIM_PANEL:
        if g in gene2col and g not in regs: regs.append(g)
    missing_reg=[str(g) for g in set(perturbed+TFS+COSTIM_PANEL) if g not in gene2col]
    log(f"regulators: {len(regs)}; missing-from-var: {missing_reg}")
    # expression columns = union(regs, effector program) -> effector genes are modelled as targets
    cols=[]
    for g in regs+EFFECTOR:
        if g in gene2col and g not in cols: cols.append(g)
    col_idx=[gene2col[g] for g in cols]
    log(f"expression columns (union reg+effector): {len(cols)}")

    def col_pos(g): return (X[:,gene2col[g]].toarray().ravel() > 0)
    cd8a=col_pos("CD8A"); cd8b=col_pos("CD8B"); cd4=col_pos("CD4")
    masks={"CD8": cd8a|cd8b, "CD4": cd4 & ~cd8a & ~cd8b}

    # --- monkeypatch: arboreto 0.1.6 calls from_delayed([]) for the unused meta frame
    #     when include_meta=False; dask>=2024 raises on empty list. Return empty frame. ---
    import arboreto.core as ac
    import dask.dataframe as dd
    _orig_fd = ac.from_delayed
    def _patched_fd(dfs, *a, **k):
        if isinstance(dfs,(list,tuple)) and len(dfs)==0:
            meta=k.get("meta")
            if meta is None and a: meta=a[0]
            return dd.from_pandas(meta.copy(), npartitions=1)
        return _orig_fd(dfs,*a,**k)
    ac.from_delayed=_patched_fd
    log("applied from_delayed empty-list patch")

    from arboreto.algo import grnboost2
    from distributed import Client, LocalCluster
    log(f"starting LocalCluster n_workers={NWORKERS} (workers hold only the ~{len(cols)}-gene matrix)")
    cluster=LocalCluster(n_workers=NWORKERS, threads_per_worker=1, memory_limit="3GB",
                         processes=True, dashboard_address=None)
    client=Client(cluster)
    log("cluster up:", client.dashboard_link if client.dashboard_link else "(no dashboard)")

    results={}
    for lin, mask in masks.items():
        t0=time.time()
        idx=np.where(mask)[0]
        Xd=X[idx][:,col_idx].toarray().astype(np.float32)
        ex=pd.DataFrame(Xd, columns=cols)
        log(f"{lin}: expression {ex.shape}; running grnboost2...")
        net=grnboost2(expression_data=ex, tf_names=regs, client_or_address=client,
                      seed=SEED, verbose=False)   # native cols: TF, target, importance
        log(f"{lin}: grnboost2 returned {net.shape[0]} edges in {time.time()-t0:.1f}s")
        # sign EVERY edge by within-lineage Pearson corr(TF,target)
        all_tgt=[g for g in cols]
        corr=signed_corr(Xd, cols, regs, all_tgt)
        rpos={g:i for i,g in enumerate(regs)}; tpos={g:i for i,g in enumerate(all_tgt)}
        signs=[]
        for tf,tg in zip(net["TF"].values, net["target"].values):
            c=corr[rpos[tf], tpos[tg]]
            signs.append(0.0 if np.isnan(c) else float(np.sign(c)))
        net["corr_sign"]=signs
        net["signed_importance"]=net["importance"]*net["corr_sign"]
        # FULL native adjacency (all edges) with EXACT [TF,target,importance] first
        net=net.sort_values("importance", ascending=False).reset_index(drop=True)
        net.to_csv(f"{WORK}/grnboost2_{lin}.csv", index=False)
        results[lin]=dict(n_edges_total=int(net.shape[0]), t_sec=round(time.time()-t0,1))
        log(f"{lin}: saved FULL {net.shape[0]} edges -> grnboost2_{lin}.csv")

    client.close(); cluster.close()
    with open(f"{WORK}/step2_summary.json","w") as f:
        json.dump(dict(regulators=regs, missing_reg=missing_reg, cols=cols,
                       EFFECTOR=EFFECTOR, COSTIM=COSTIM, COSTIM12=COSTIM12, COSTIM_PANEL=COSTIM_PANEL,
                       results=results, seed=SEED, nworkers=NWORKERS), f, indent=2)
    log("STEP2 DONE")

if __name__=="__main__":
    main()
