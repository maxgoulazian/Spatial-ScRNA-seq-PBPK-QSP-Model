#!/usr/bin/env python
"""
STEP 1 - Perturbation-response network (primary causal GOF readout) + lineage split.
Dataset: GSE190604 CRISPRa Perturb-seq = CD4+CD8 POOL (not CD8-only).
Split: CD8 lineage = CD8A>0 or CD8B>0 ; CD4 lineage = CD4>0 & CD8A==0 & CD8B==0.
For EACH lineage: DE signature per perturbed target = mean(log1p target) - mean(log1p NO-TARGET).
Directed target->gene edges (agonism = direct CRISPRa direction, GOF).
"""
import os, sys, json, time
os.environ["NUMEXPR_MAX_THREADS"]="64"; os.environ["NUMEXPR_NUM_THREADS"]="16"
os.environ["OMP_NUM_THREADS"]="8"; os.environ["OPENBLAS_NUM_THREADS"]="8"; os.environ["MKL_NUM_THREADS"]="8"
import numpy as np, scipy.sparse as sp, h5py

def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

BASE="/media/balthasar-lab/RAID4/costim_engager_counterscreen"
F=f"{BASE}/data/schmidt2022_tcell_perturbseq/cd8_grn_input.h5ad"
WORK=f"{BASE}/work/cd8cd4_gof_grn"; os.makedirs(WORK, exist_ok=True)

EFFECTOR=["IFNG","GZMB","GZMA","GZMK","PRF1","TNF","IL2","FASLG","NKG7","KLRG1","GNLY","CCL4","CCL5","XCL1"]
COSTIM=["CD28","CD27","TNFRSF9","CD2","LTBR"]
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

log("loading", F)
with h5py.File(F,"r") as h:
    shape=tuple(h["X"].attrs["shape"])
    X=sp.csr_matrix((h["X/data"][:],h["X/indices"][:],h["X/indptr"][:]),shape=shape)
    sym=read_col(h,"var/sym"); ens=read_col(h,"var/_index")
    guide=read_col(h,"obs/guide_target").astype(str)
gene2col={g:i for i,g in enumerate(sym)}
G=X.shape[1]
log(f"X {X.shape} nnz={X.nnz:,}")
np.save(f"{WORK}/var_sym.npy", sym); np.save(f"{WORK}/var_ens.npy", ens)

targets=sorted(set(guide)-{"NO-TARGET"})
assert all(g in gene2col for g in EFFECTOR+COSTIM+TFS), "missing gene"

def col_pos(g): return (X[:,gene2col[g]].toarray().ravel() > 0)
cd8a=col_pos("CD8A"); cd8b=col_pos("CD8B"); cd4=col_pos("CD4")
masks={"CD8": cd8a|cd8b, "CD4": cd4 & ~cd8a & ~cd8b}
log(f"CD8A+ {cd8a.mean():.3f} CD8B+ {cd8b.mean():.3f} CD4+ {cd4.mean():.3f}")
log(f"CD8 {masks['CD8'].sum():,} | CD4 {masks['CD4'].sum():,} | neither {(~masks['CD8']&~masks['CD4']).sum():,}")

summary={"lineages":{}, "genes":{"n_var":int(G)},
         "EFFECTOR":EFFECTOR,"COSTIM":COSTIM,"TFS":TFS,"targets":targets}

for lin, mask in masks.items():
    t0=time.time()
    idx=np.where(mask)[0]; g_lin=guide[idx]; Xl=X[idx].tocsr()
    ctrl_sel=(g_lin=="NO-TARGET"); n_ctrl=int(ctrl_sel.sum())
    Xc=Xl[ctrl_sel]
    ctrl_mean=np.asarray(Xc.mean(axis=0)).ravel()
    ctrl_ex2=np.asarray(Xc.multiply(Xc).mean(axis=0)).ravel()
    ctrl_var=np.maximum(ctrl_ex2-np.square(ctrl_mean),0.0)
    tgs=list(targets); tpos={t:i for i,t in enumerate(tgs)}
    ncell={t:int((g_lin==t).sum()) for t in tgs}
    rows=[];cols=[];vals=[]
    for ci,tg in enumerate(g_lin):
        if tg in tpos and ncell[tg]>0:
            rows.append(tpos[tg]); cols.append(ci); vals.append(1.0/ncell[tg])
    Ind=sp.csr_matrix((vals,(rows,cols)),shape=(len(tgs),Xl.shape[0]))
    tgt_mean=(Ind@Xl).toarray()              # n_tgt x G dense (sparse@sparse -> sparse, densify)
    Xsq=Xl.multiply(Xl)
    tgt_ex2=(Ind@Xsq).toarray()
    tgt_var=np.maximum(tgt_ex2-np.square(tgt_mean),0.0)
    lfc=tgt_mean-ctrl_mean[None,:]
    n_t=np.array([ncell[t] for t in tgs],float)[:,None]
    se=np.sqrt(tgt_var/np.maximum(n_t,1)+ctrl_var[None,:]/max(n_ctrl,1))
    tstat=np.divide(lfc,se,out=np.zeros_like(lfc),where=se>0)
    np.savez_compressed(f"{WORK}/pertresp_{lin}.npz",
        lfc=lfc.astype(np.float32), tstat=tstat.astype(np.float32),
        tgt_mean=tgt_mean.astype(np.float32), ctrl_mean=ctrl_mean.astype(np.float32),
        targets=np.array(tgs), ncell=np.array([ncell[t] for t in tgs]))
    ifc=gene2col["IFNG"]
    ifng_lfc={c: round(float(lfc[tpos[c],ifc]),5) for c in COSTIM}
    ifng_t={c: round(float(tstat[tpos[c],ifc]),3) for c in COSTIM}
    summary["lineages"][lin]=dict(n_cells=int(mask.sum()), n_ctrl=n_ctrl,
        costim_ncell={c:ncell[c] for c in COSTIM}, ifng_lfc=ifng_lfc, ifng_t=ifng_t)
    log(f"{lin}: {mask.sum():,} cells, {n_ctrl} ctrl, lfc {lfc.shape} in {time.time()-t0:.1f}s")
    log(f"{lin} IFNG dlog1p per arm: {ifng_lfc}")

with open(f"{WORK}/step1_summary.json","w") as f: json.dump(summary,f,indent=2)
log("STEP1 DONE ->", f"{WORK}/step1_summary.json")
