"""coupled_percell_pk.py — whole-body PK with per-organ vascular compartment + per-cell interstitial graph.
Structure (user design): systemic plasma -Q-> organ vascular compartment -extravasation@BEC-> per-cell graph
(diffuse+bind) -drain@LEC-> systemic lymph pool -> plasma. Shah-Betts fixed physiology (Q, L, sigV, sigL,
Vv, Vis, k_dist, k_cat, k_lymph_return); only KD/n_arm/kint/mw/fFcRn/pools are per-molecule.
"""
import numpy as np, os, time, sys
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from wholebody_percell import TissueGraph, STROMAL_LABELS, AVO

def load_tissue(organ, target_col, agents_dir, heme_npz="heme_abm_cells.npz"):
    if organ=="tumor":
        f=os.path.join(agents_dir,"..","tumor_abm_cells.npz")
    elif organ=="heme_tumor":
        # heme malignancy compartment (AML/myeloma/lymphoma) — per-target npz from routing map
        f=os.path.join(agents_dir,heme_npz)
        if not os.path.exists(f): f=os.path.join(agents_dir,"..",heme_npz)
    elif organ=="blood":
        # normal circulating blood cells (Tabula Sapiens) — always-present compartment (TMDD sink + PD + CRS)
        f=os.path.join(agents_dir,"agents_blood.npz")
        if not os.path.exists(f): f=os.path.join(agents_dir,"..","agents_blood.npz")
    else:
        f=os.path.join(agents_dir, f"agents_{organ}.npz")
    z=np.load(f,allow_pickle=True)
    # FIX-3: HGNC-symbol <-> common-name aliases (SAME gene — safe identities, not affinity/epitope claims),
    # tried ONLY as a fallback so a resolving lookup is never changed; rescues the symbol-schism that otherwise
    # binds at silent ZERO density. On a true miss, WARN LOUDLY (not silent) — the run continues (zero density)
    # but the invalid target is visible in the log rather than producing quiet garbage.
    _COLNAME_ALIAS={'MS4A1':'CD20','CD20':'MS4A1','ERBB2':'HER2','HER2':'ERBB2','FOLH1':'PSMA','PSMA':'FOLH1',
                    'TNFRSF17':'BCMA','BCMA':'TNFRSF17','FCRL5':'FcRH5','FcRH5':'FCRL5'}
    def _resolve_density_col(cands_primary):
        cands=list(cands_primary)
        _al=_COLNAME_ALIAS.get(target_col)
        if _al: cands += ["dens_"+_al, _al+"_copies", _al]
        col=next((c for c in cands if c in z.files), None)
        if col is None:
            import sys as _sys
            _sys.stderr.write(f"[FIX-3] {organ}: NO density column for target '{target_col}' (tried {cands}); "
                              f"ZERO density -> results for this target are INVALID. "
                              f"npz density-like cols: {[c for c in z.files if ('dens' in c or 'copies' in c)][:10]}\n")
        return col
    if organ in ("tumor","heme_tumor","blood"):
        x,y=z["X"].astype(float),z["Y"].astype(float)
        # decode per-cell labels: heme npz stores per-cell 'labels' (len==N); tumor npz stores a
        # category-name list in 'labels' + per-cell integer codes in 'lab' -> index names by codes.
        if "labels" in z.files and len(z["labels"])==len(x):
            labs=np.asarray(z["labels"])
        elif "lab" in z.files:
            names=z["labels"] if "labels" in z.files else z["labels_map"]
            labs=np.asarray(names)[np.asarray(z["lab"]).astype(int)]
        else:
            labs=z["labels"]
        # tumor uses 'dens_<gene>'; heme uses '<gene>_copies'. Also accept raw '<gene>'.
        col=_resolve_density_col(["dens_"+target_col, target_col+"_copies", target_col])
        ag=z[col].astype(float) if col else np.zeros(len(x))
    else:
        x,y=z["x"].astype(float),z["y"].astype(float); labs=z["fine"] if "fine" in z.files else z["coarse"]
        _c2=_resolve_density_col([target_col])
        ag=z[_c2].astype(float) if _c2 else np.zeros(len(x))
    good=np.isfinite(x)&np.isfinite(y); return x[good],y[good],ag[good],np.asarray(labs)[good],good

class CoupledPerCellPK:
    def __init__(self, organs, target, target_col, KD_nM, n_arm, kint_perday, mw_kda, fFcRn,
                 agents_dir, pools, Q, L, sigV, Vis, Vv, bec_lec, wellmixed=None, sigL=0.85, k_dist=3.0,
                 k_lymph_return=24.0, CLup=0.3503, k_renal=0.0, V_pl=3.1, V_ly=2.6, dt=0.01, D_um2s=10.0,
                 kon1_perM_pers=None, koff1_pers=None, kdeg_perday=None, span_coeng_nm=12.5):
        self.organs=organs; self.KD=KD_nM; self.n_arm=n_arm; self.kint=kint_perday
        # kinetic-TMDD params: kon [/M/s]->[/nM/day], koff [/s]->[/day]; None -> QSS fallback in the step
        _SPD=86400.0
        self.kon1 = (kon1_perM_pers/1e9*_SPD) if kon1_perM_pers is not None else None
        self.koff1 = (koff1_pers*_SPD) if koff1_pers is not None else None
        self.kdeg = kdeg_perday
        self.span_coeng_nm = span_coeng_nm
        self.mw=mw_kda; self.V_pl=V_pl; self.V_ly=V_ly; self.dt=float(__import__('os').environ.get('PD_DT', dt))
        self.Q=Q; self.L=L; self.sigV=sigV; self.Vis=Vis; self.Vv=Vv; self.sigL=sigL; self.k_dist=k_dist
        self.k_lymph_return=k_lymph_return; self.k_cat=CLup*(1.0-fFcRn)+k_renal
        self.pools=pools; self.graphs={}; self._agents_dir=agents_dir
        self.wm = dict(wellmixed) if wellmixed else {}
        self.A_wm = {o:0.0 for o in self.wm}
        t0=time.time()
        for o in organs:
            x,y,ag,labs,good=load_tissue(o,target_col,agents_dir)
            lab_s=np.array([str(v) for v in labs]); stro=np.array([v in STROMAL_LABELS for v in lab_s],float)
            pool_nM=pools.get(o,0.0); ag_sum=ag.sum()
            # UNIFIED with PD + heme/blood: keep the REAL per-cell receptor copies (Glassman/IHC, per
            # RECEPTORS_PER_CELL_METHOD.md) so binding is computed at physiological receptor number
            # (correct avidity for bivalent) and per-cell heterogeneity is preserved. Scale only the CELL
            # COUNT via count_scale = (physiological pool copies)/(sum of sampled real copies), so the summed
            # organ sink still integrates to the physiological receptor pool. For monovalent this is numerically
            # identical to the old pool-redistribution; for bivalent (ternary) it is the physically correct path.
            if pool_nM>0 and ag_sum>0:
                tot_copies=pool_nM*Vis[o]*AVO/1e9
                Rpercell=np.asarray(ag,float)          # REAL per-cell copies (NOT redistributed)
                cscale=tot_copies/ag_sum               # count_scale: sampled cells -> physiological pool
            else:
                Rpercell=np.zeros(len(ag)); cscale=0.0
            # map stored bec/lec indices (built on full table) through the finite-coord mask
            bl=bec_lec[o]; remap=-np.ones(len(good),int); remap[np.where(good)[0]]=np.arange(good.sum())
            bec=remap[np.asarray(bl['bec'],int)]; bec=bec[bec>=0]
            lec=remap[np.asarray(bl['lec'],int)]; lec=lec[lec>=0]
            g=TissueGraph(o,x,y,stro,labs,bec,lec,Vis[o],Vv[o],dt=dt,D_um2s=D_um2s)
            g.set_antigen(Rpercell, kdeg_perday=self.kdeg, kon1=self.kon1, koff1=self.koff1, span_coeng_nm=self.span_coeng_nm); g.count_scale=cscale
            g.labs=lab_s          # per-cell type labels (for snapshots)
            self.graphs[o]=g
        self.build_s=time.time()-t0
        self.pool=ThreadPoolExecutor(max_workers=min(len(organs),9))

    def simulate(self, schedule, tsim, route="IV", F_sc=0.6, ka_sc=0.25, inf_rate=0.0, inf_dur=0.0, rec=400, snap_times=None):
        dt=self.dt; nstep=int(tsim/dt)+1
        A_pl=0.0; A_ly=0.0; A_sc=0.0
        sched=sorted(schedule); si=0
        t_rec=[]; C_rec=[]; sink_rec={o:[] for o in self.organs}; bound_rec={o:[] for o in self.organs}
        snap_times = sorted(snap_times) if snap_times else []
        self.snaps={}; _snap_i=0
        for k in range(nstep):
            t=k*dt
            while si<len(sched) and sched[si][0]<=t+1e-9:
                mg=sched[si][1]
                if route=="SC": A_sc+=mg/self.mw*1e3
                else: A_pl+=mg/self.mw*1e3
                si+=1
            C_pl=A_pl/self.V_pl
            J_sc=ka_sc*A_sc; A_sc-=dt*J_sc
            infn=(inf_rate/self.mw*1e3) if (inf_rate>0 and t<=inf_dur) else 0.0
            def _one(o):
                g=self.graphs[o]
                PS_ex=self.k_dist*self.L[o]*(1.0-self.sigV[o]); PS_ly=self.k_dist*self.L[o]*(1.0-self.sigL)
                sink,drain,J_extrav,S,D=g.step(C_pl,self.Q[o],PS_ex,PS_ly,self.KD,self.n_arm,self.kint)
                return o,sink,drain,J_extrav,S
            tot_drain=0.0; tot_extrav=0.0
            results=list(self.pool.map(_one,self.organs)) if self.pool else [ _one(o) for o in self.organs]
            do_rec = (k % max(1,nstep//rec)==0)
            for o,sink,drain,J_extrav,S in results:
                tot_extrav += J_extrav; tot_drain += drain
                if do_rec:
                    sink_rec[o].append(sink)
                    bound_rec[o].append(float(S.sum()/AVO*1e9/self.Vis[o]) if self.pools.get(o,0)>0 else 0.0)
            # spatial snapshot: full per-cell fields at requested times (for overlay/per-celltype figures)
            if _snap_i < len(snap_times) and t >= snap_times[_snap_i]-1e-9:
                tsnap=snap_times[_snap_i]; self.snaps[tsnap]={}
                Sd={o:S for o,_,_,_,S in results}
                for o in self.organs:
                    g=self.graphs[o]
                    self.snaps[tsnap][o]=dict(x=g.x.copy(), y=g.y.copy(), C=g.C.copy(),
                        S=Sd[o].copy(), R=g.R.copy(), labs=g.labs, bec=g.bec.copy(), lec=g.lec.copy())
                _snap_i+=1
            for o,(Lo,sVo,Viso) in self.wm.items():
                Cis = self.A_wm[o]/max(Viso,1e-9)
                Jex = self.k_dist*Lo*(1.0-sVo)*C_pl
                Jre = self.k_dist*Lo*(1.0-self.sigL)*Cis
                self.A_wm[o] = max(self.A_wm[o] + dt*(Jex - Jre), 0.0)
                tot_extrav += Jex; tot_drain += Jre
            dA_pl = infn + F_sc*J_sc + self.k_lymph_return*A_ly - self.k_cat*A_pl - tot_extrav
            dA_ly = tot_drain - self.k_lymph_return*A_ly
            A_pl=max(A_pl+dt*dA_pl,0.0); A_ly=max(A_ly+dt*dA_ly,0.0)
            if k % max(1,nstep//rec)==0:
                t_rec.append(t); C_rec.append(A_pl/self.V_pl*self.mw/1e3)
        return dict(t=np.array(t_rec), Cplasma_ugml=np.array(C_rec),
                    sink=sink_rec, bound_nM=bound_rec, build_s=self.build_s)
