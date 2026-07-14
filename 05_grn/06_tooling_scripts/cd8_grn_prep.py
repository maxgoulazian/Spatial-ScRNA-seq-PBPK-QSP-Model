"""Preprocess Schmidt CD8 CRISPRa Perturb-seq for GRN inference.
Loads GSE190604 mtx (103,805 cells x 36,601 genes + 154 guides),
assigns each cell its guide target, QC-filters, saves clean AnnData.
"""
import scanpy as sc, anndata as ad, numpy as np, pandas as pd, scipy.io as sio, scipy.sparse as sp, gzip, os, sys, time
t0=time.time()
D="/media/balthasar-lab/RAID4/costim_engager_counterscreen/data/schmidt2022_tcell_perturbseq/GSE190604_CRISPRa_Perturbseq"
OUT="/media/balthasar-lab/RAID4/costim_engager_counterscreen/data/schmidt2022_tcell_perturbseq/cd8_grn_input.h5ad"

def log(m): print(f"[+{time.time()-t0:5.0f}s] {m}", flush=True)

log("reading features/barcodes")
feats=pd.read_csv(f"{D}/GSE190604_features.tsv.gz",sep="\t",header=None,names=["ens","sym","ftype"])
bcs=pd.read_csv(f"{D}/GSE190604_barcodes.tsv.gz",header=None,names=["bc"])
is_gex=(feats.ftype=="Gene Expression").values
is_grna=(feats.ftype=="CRISPR Guide Capture").values
log(f"{is_gex.sum()} genes, {is_grna.sum()} guides, {len(bcs)} cells")

log("reading mtx (264M nonzeros — ~2-4 min)")
M=sio.mmread(f"{D}/GSE190604_matrix.mtx.gz").tocsr()   # features x cells
log(f"mtx shape {M.shape}")
X=M[is_gex,:].T.tocsr()   # cells x genes
log(f"GEX cells x genes = {X.shape}")

A=ad.AnnData(X=X, obs=pd.DataFrame(index=bcs.bc.values),
             var=pd.DataFrame({"sym":feats.sym[is_gex].values}, index=feats.ens[is_gex].values))
A.var_names_make_unique()

# assign guide from guidecalls (singlets only)
log("assigning guides")
gc=pd.read_csv(f"{D}/GSE190604_cellranger-guidecalls-aggregated-unfiltered.txt.gz",sep="\t")
gc=gc[gc.num_features==1].copy()                       # singlet guide only
gc["target"]=gc.feature_call.str.replace(r"-\d+$","",regex=True)  # ABCB10-1 -> ABCB10
g=gc.set_index("cell_barcode")["target"]
A.obs["guide_target"]=A.obs_names.map(g)
A.obs["is_control"]=A.obs.guide_target.eq("NO-TARGET")
n_assigned=A.obs.guide_target.notna().sum()
log(f"{n_assigned} cells have singlet guide ({100*n_assigned/A.n_obs:.0f}%)")

# keep only guide-assigned cells, basic QC
A=A[A.obs.guide_target.notna()].copy()
sc.pp.filter_cells(A, min_genes=500)
sc.pp.filter_genes(A, min_cells=20)
A.obs["n_counts"]=np.asarray(A.X.sum(1)).ravel()
log(f"after QC: {A.shape}")

# normalize + log for GRN
A.layers["counts"]=A.X.copy()
sc.pp.normalize_total(A, target_sum=1e4)
sc.pp.log1p(A)
log(f"guide targets: {A.obs.guide_target.nunique()} unique")
print(A.obs.guide_target.value_counts().head(20).to_string(), flush=True)

A.write(OUT)
log(f"saved {OUT} ({os.path.getsize(OUT)/1e6:.0f} MB)")
