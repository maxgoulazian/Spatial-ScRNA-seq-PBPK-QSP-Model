"""wholebody_percell.py — production whole-body per-cell engine (BEC-entry / LEC-exit structure).

Physiological structure (per organ, per user design):
  systemic plasma --Q blood flow--> ORGAN VASCULAR COMPARTMENT (Vv) --extravasation at BEC cells-->
  per-cell interstitial graph (diffuse, ECM-hindered; each cell its own Rhoden binder) --drain at LEC cells-->
  organ lymph --> systemic lymph pool --> plasma (recirculation).

Every barcode is (1) a node in the kNN spatial diffusion graph AND (2) its own Rhoden binder. Drug ENTERS
only at blood-endothelial (BEC) cells and LEAVES only at lymphatic-endothelial (LEC) cells — the real
vasculature entry / lymphatic exit points from the Xenium markers (LYVE1/PROX1/PDPN). TMDD emerges as
Σ(per-cell internalization). Transport: implicit-Euler + prefactored sparse LU (row-sum-zero Laplacian).
"""
import numpy as np, os
from scipy import sparse
from scipy.sparse.linalg import splu
from sklearn.neighbors import NearestNeighbors

# ---- Tissue-graph cache directory --------------------------------------------
# Where to cache/reuse the expensive per-organ diffusion graphs. Override with env
# WB_GRAPH_CACHE (set to "" or "0"/"off" to DISABLE caching and always rebuild).
_env_cache=os.environ.get("WB_GRAPH_CACHE", None)
if _env_cache is None:
    _GRAPH_CACHE_DIR="/media/balthasar-lab/RAID4/atlas_spatial_omics/organ_spatial/graph_cache"
elif _env_cache.strip().lower() in ("","0","off","none","false"):
    _GRAPH_CACHE_DIR=None
else:
    _GRAPH_CACHE_DIR=_env_cache

AVO=6.02214076e23
ECM_GENES=["COL11A1","FAP","VCAN","THBS1","ACTA2","TAGLN","PDGFRA","THBS4","COL17A1","LAMC2"]
STROMAL_LABELS={"hepatic stellate cell","fibroblast","myofibroblast cell","fibroblast_stromal",
                "Fibroblast_stromal","stromal cell","smooth muscle cell","pericyte"}
R_CELL_UM=8.0; R_AB_UM_DEFAULT=0.0125

def percell_ageff_nM(R_copies,r_cell_um=R_CELL_UM,r_Ab_um=R_AB_UM_DEFAULT):
    N=1e9; Ag_bulk=R_copies*N/AVO*1e9
    SA_cell=4.0*np.pi*r_cell_um**2; SA_Ab=np.pi*r_Ab_um**2; V_Ab=(2.0/3.0)*np.pi*r_Ab_um**3
    return (Ag_bulk/N/SA_cell)*SA_Ab/V_Ab*1e15

import kinetic_rhoden_percell as _KRP
_NM_PER_COPY=6.0/257000.0   # synapse reaction-volume basis; IDENTICAL to wholebody_pd.NM_PER_COPY
def percell_rhoden_qss_Cvec(C_nM_arr,R_copies,KD_nM,n_arm):
    C=np.maximum(np.asarray(C_nM_arr,float),0.0); KD=max(float(KD_nM),1e-9); Rt=np.asarray(R_copies,float)
    if n_arm<2:
        S=Rt*(C/(C+KD)); return S,np.zeros_like(S)
    Ag=percell_ageff_nM(Rt)
    a=2.0*Ag*C; b=KD*Rt*(2.0*C+KD); cc=-(KD*KD)*(Rt*Rt)
    with np.errstate(over='ignore',invalid='ignore',divide='ignore'):
        disc=b*b-4.0*a*cc
        Rfree=np.where(a>1e-30,(-b+np.sqrt(disc))/(2.0*np.maximum(a,1e-30)),Rt/np.maximum(1.0+2.0*C/KD,1e-30))
    Rfree=np.clip(Rfree,0.0,Rt); S=(2.0*C/KD)*Rfree; D=Ag*Rfree*S/(2.0*np.maximum(Rt,1e-30)*KD)
    return S,D


# ---- Tissue-graph disk cache -------------------------------------------------
# The per-organ graph (kNN Laplacian L, LEC-distance drainage weights, BEC/LEC masks)
# is a PURE FUNCTION of (coords, ECM, bec/lec indices, k, D, alpha_D). It is the
# expensive part of the build (kNN over up to ~400k cells). Cache it keyed on a hash
# of those inputs: first run builds+saves, later runs load in ~ms. dt-dependent splu
# is always rebuilt from cached L (SuperLU is not picklable, and dt can vary per run).
import hashlib as _hashlib
def _graph_cache_key(organ, x, y, ecm, bec_idx, lec_idx, k, D_um2s, alpha_D):
    h=_hashlib.sha1()
    h.update(organ.encode())
    for a in (np.ascontiguousarray(x,dtype=np.float64),
              np.ascontiguousarray(y,dtype=np.float64),
              np.ascontiguousarray(ecm,dtype=np.float64),
              np.ascontiguousarray(np.asarray(bec_idx,dtype=np.int64)),
              np.ascontiguousarray(np.asarray(lec_idx,dtype=np.int64))):
        h.update(a.tobytes())
    h.update(np.array([k,D_um2s,alpha_D],dtype=np.float64).tobytes())
    return h.hexdigest()[:20]

class TissueGraph:
    """One tissue: per-cell kNN diffusion graph + explicit vascular compartment; BEC entry, LEC exit."""
    def __init__(self, organ, x, y, ecm, labels, bec_idx, lec_idx, Vis_L, Vv_L,
                 dt=0.01, k=6, D_um2s=10.0, alpha_D=3.0):
        self.organ=organ; self.x=x; self.y=y; self.n=len(x); self.dt=dt
        self.Vis=Vis_L; self.Vv=Vv_L
        self.ecm = ecm if ecm is not None else np.zeros(self.n)
        self.phi_D = 1.0/(1.0+alpha_D*self.ecm)
        # ---- CACHE CHECK: the graph (L, BEC/LEC masks, LEC-distance drainage) is a pure
        # function of (coords, ECM, bec/lec idx, k, D, alpha_D). Load if present, else build+save.
        cached=None
        if _GRAPH_CACHE_DIR is not None:
            key=_graph_cache_key(organ,x,y,self.ecm,bec_idx,lec_idx,k,D_um2s,alpha_D)
            cpath=os.path.join(_GRAPH_CACHE_DIR,f"graph_{organ}_{key}.npz")
            if os.path.exists(cpath):
                try:
                    z=np.load(cpath)
                    self.L=sparse.csr_matrix((z['L_data'],z['L_indices'],z['L_indptr']),shape=tuple(z['L_shape']))
                    self.bec=z['bec']; self.lec=z['lec']; self.drain_w=z['drain_w']
                    self.n_bec=int(self.bec.sum()); self.n_lec=int(self.lec.sum())
                    cached=cpath
                except Exception:
                    cached=None
        if cached is None:
            # kNN graph Laplacian (1/day), ECM-hindered edge conductance, row-sum-zero (mass-conserving)
            P=np.column_stack([x,y]); kk=min(k,max(2,self.n-1))
            nn=NearestNeighbors(n_neighbors=kk+1).fit(P); dist,idx=nn.kneighbors(P)
            Dcell=D_um2s*self.phi_D
            rows=[];cols=[];vals=[]
            for j in range(1,kk+1):
                nb=idx[:,j]; d=np.maximum(dist[:,j],0.003)  # 3 nm membrane-to-membrane physical floor
                ge=0.5*(Dcell+Dcell[nb])*86400.0/(d*d)
                rows.append(np.arange(self.n));cols.append(nb);vals.append(ge)
            W=sparse.csr_matrix((np.concatenate(vals),(np.concatenate(rows),np.concatenate(cols))),shape=(self.n,self.n))
            W=W.maximum(W.T); deg=np.asarray(W.sum(axis=1)).ravel()
            self.L=(W - sparse.diags(deg)).tocsr()              # row sums EXACTLY 0
            # BEC (entry) and LEC (exit) masks
            self.bec=np.zeros(self.n,bool); self.bec[np.asarray(bec_idx,int)]=True
            self.lec=np.zeros(self.n,bool); self.lec[np.asarray(lec_idx,int)]=True
            if self.bec.sum()==0: self.bec[np.argsort(x+y)[:max(1,self.n//50)]]=True
            if self.lec.sum()==0: self.lec[np.argsort(x+y)[-max(1,self.n//100):]]=True
            self.n_bec=int(self.bec.sum()); self.n_lec=int(self.lec.sum())
            # distance-graded lymphatic drainage weight: whole tissue drains, but rate decays with distance from
            # the nearest LEC (initial lymphatics are distributed; interstitial fluid convects toward them, fastest
            # near the vessels). lam = lymphatic catchment length (~100 um, physiological inter-lymphatic spacing).
            lam=100.0  # um
            if self.n_lec>0:
                lec_xy=np.column_stack([x[self.lec],y[self.lec]])
                nnl=NearestNeighbors(n_neighbors=1).fit(lec_xy)
                dl,_=nnl.kneighbors(np.column_stack([x,y]))   # distance to nearest LEC (um)
                self.drain_w = np.exp(-dl.ravel()/lam) + 0.05  # graded weight + small floor so ENTIRE tissue drains
            else:
                self.drain_w = np.ones(self.n)
            if _GRAPH_CACHE_DIR is not None:
                Lc=self.L.tocsr()
                try:
                    os.makedirs(_GRAPH_CACHE_DIR,exist_ok=True)
                    np.savez(cpath, L_data=Lc.data, L_indices=Lc.indices, L_indptr=Lc.indptr,
                             L_shape=np.array(Lc.shape), bec=self.bec, lec=self.lec, drain_w=self.drain_w)
                except Exception:
                    pass
        # prefactor pure-diffusion implicit operator (dt-dependent -> always rebuilt from cached L)
        I=sparse.identity(self.n,format="csc")
        self.lu=splu((I-dt*self.L).tocsc())
        self.C=np.zeros(self.n)          # per-cell interstitial free-drug (nM)
        self.Cvasc=0.0                   # organ vascular-compartment free-drug (nM)

    count_scale=1.0   # sampled cells -> physiological pool (set by builder; heme/blood pattern)
    def set_antigen(self,R_copies,kdeg_perday=None,kon1=None,koff1=None,span_coeng_nm=12.5):
        """Init per-cell receptor STATE for kinetic TMDD. R_copies -> Ag0 set-point (nM on pericellular
        volume); live pools (free Ag, singly-bound BAg1, avidity Bdbl) carried across steps. If kdeg/kon/koff
        are None the step falls back to the legacy QSS sink (backward compatible)."""
        self.R=np.asarray(R_copies,float)
        # synapse reaction-volume basis — IDENTICAL to PD (wholebody_pd.NM_PER_COPY = 6.0/257000 nM/copy)
        self.Ag0_nM=self.R*_NM_PER_COPY
        self.Ag_nM=self.Ag0_nM.copy(); self.BAg1=np.zeros(self.n); self.Bdbl=np.zeros(self.n)
        self.kin_kdeg=kdeg_perday; self.kin_kon1=kon1; self.kin_koff1=koff1
        self.kin_span_coeng_nm=float(span_coeng_nm)   # per-molecule co-engagement span (nm) for bivalent avidity
        self.kinetic_tmdd=(kdeg_perday is not None and kon1 is not None and koff1 is not None)

    def step(self, C_plasma_nM, Q_Lday, PS_extrav_Lday, PS_lymph_Lday, KD_nM, n_arm, kint_perday):
        """One dt. Physiological BEC-in/LEC-out transport.
        C_plasma_nM : systemic plasma conc feeding this organ's vascular compartment.
        Q_Lday      : organ blood flow (vascular exchange with systemic plasma).
        PS_extrav   : permeability-surface product for BEC extravasation (L/day) = k_dist*L*(1-sigV).
        PS_lymph    : lymphatic drainage conductance (L/day) = k_dist*L*(1-sigL).
        Returns (organ_sink nmol/day into TMDD, drain_to_lymph nmol/day, S, D).
        """
        dt=self.dt
        # (1) organ vascular compartment — QUASI-STEADY (blood flow Q >> extravasation PS, so the vascular
        #     pool equilibrates near-instantly each step; explicit integration would be catastrophically
        #     stiff at Q/Vv ~ 1e5/day). Solve 0 = Q*(Cpl - Cvasc) - PS_extrav*(Cvasc - Cis_bec):
        #        Cvasc = (Q*Cpl + PS_extrav*Cis_bec) / (Q + PS_extrav)
        Cis_bec = self.C[self.bec].mean() if self.n_bec else 0.0
        self.Cvasc = (Q_Lday*C_plasma_nM + PS_extrav_Lday*Cis_bec)/max(Q_Lday+PS_extrav_Lday,1e-12)
        J_extrav = PS_extrav_Lday*(self.Cvasc - Cis_bec)          # nmol/day vascular->interstitium (at BEC)
        # (2) interstitial per-cell graph: diffuse (implicit) + BEC source + LEC drain + per-cell binding
        vcell=max(self.Vis/self.n,1e-18)
        if getattr(self,'kinetic_tmdd',False):
            # KINETIC Rhoden per-cell binding with receptor turnover (ONE solve; feeds PK sink here + PD).
            # SAME-antigen valency: n_arm>=2 -> bivalent avidity (Bdbl forms via BAg1 crosslink);
            # n_arm<2 -> monovalent (AgEFF=0, Bdbl stays 0). Backward-Euler: stable + census-exact + fast
            # (nsub from slow scales, NOT pinned by the stiff crosslink). Replaces the old two-pool
            # rhoden_bivalent_step call whose Bdbl was INERT for single-antigen (fed from BAg2=0).
            Ag1EFF = _KRP.geo_ageff_nM(self.R, r_cell_um=8.0, span_nm=getattr(self,'kin_span_coeng_nm',12.5)) if n_arm>=2 else np.zeros(self.n)
            self.Ag_nM, self.BAg1, self.Bdbl, intern_flux_nM_day = _KRP.rhoden_samecell_bivalent_step(
                self.C, self.Ag_nM, self.BAg1, self.Bdbl, self.Ag0_nM,
                self.kin_kon1, self.kin_koff1, Ag1EFF, self.kin_kdeg, kint_perday, dt)
            # intern_flux is nM/day of internalized drug on the SYNAPSE volume; convert to copies/cell/day
            # via NM_PER_COPY (copies = nM / NM_PER_COPY), matching how PD counts bound receptors.
            intern_copies_cell = (intern_flux_nM_day/_NM_PER_COPY)*self.count_scale   # copies/cell/day
            loss_nM_day = intern_copies_cell/AVO/vcell*1e9                            # -> interstitial nM/day (drug removal)
            S = (self.BAg1 + self.Bdbl)/_NM_PER_COPY   # bound copies/cell (PD readout/report)
            D = self.Bdbl/_NM_PER_COPY
        else:
            S,D=percell_rhoden_qss_Cvec(self.C,self.R,KD_nM,n_arm)
            intern_copies_cell = kint_perday*(S+2.0*D)*self.count_scale
            loss_nM_day = intern_copies_cell/AVO/vcell*1e9
        loss_nM_day = np.minimum(loss_nM_day, self.C/dt*0.9)
        src=np.zeros(self.n)
        # BEC extravasation source: J_extrav (nmol/day) spread over BEC cells' pericellular volume -> nM/day
        if self.n_bec:
            src[self.bec] += (J_extrav/self.Vis)*(self.n/self.n_bec)
        # Lymphatic drainage: organ-total flux = PS_lymph * Cis_mean (2-pore mass-correct), realized as a
        # distance-GRADED spatial field over the WHOLE tissue (drain_w decays with distance from nearest LEC),
        # normalized to integrate EXACTLY to that flux. Whole tissue drains -> no diffusion bottleneck / no
        # accumulation; gradient -> drug preferentially exits near the lymphatics (physiological convection).
        Cis_mean = self.C.mean()
        drain_to_lymph = PS_lymph_Lday*Cis_mean                        # nmol/day organ -> lymph
        w = self.drain_w*self.C                                        # weight * local conc
        wint = w.mean()*self.Vis                                       # ~ integral of w over interstitium (nM*L)
        drain_field = w*(drain_to_lymph/wint) if wint>1e-18 else np.zeros(self.n)  # nM/day per cell
        local = src - drain_field - loss_nM_day
        self.C = np.maximum(self.lu.solve(self.C + dt*local), 0.0)
        organ_sink=float(intern_copies_cell.sum())/AVO*1e9            # nmol/day (TMDD)
        return organ_sink, drain_to_lymph, J_extrav, S, D
