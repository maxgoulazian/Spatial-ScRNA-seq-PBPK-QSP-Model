"""wholebody_pd.py — per-cell PD everywhere (killing + cytokine + Treg suppression) for the whole-body
single-cell PBPK-PD model. Generalizes the VALIDATED tumor ABM kill law (unified_pbpk_pd.py tumor_percell:
Schropp ternary trimer -> kill_cd8 = Cb/(1+0.25*n_treg) -> apportioned to neighbor targets) to ALL organs,
per barcode. One agent = one barcode; each T cell forms its OWN trimer from its OWN local drug (from the
per-cell transport grid), its OWN CD3 density, and its neighbouring target cells' OWN antigen density.

Design (matches user constraints):
 - individual cells, not pools; per-cell Rhoden/Schropp binding; PD tracked in EVERY organ.
 - effector arm R_CD3  <- per-cell CD3E expression (agent table column 'CD3E')
 - target arm  R_TAA   <- per-cell antigen density (the same per-cell R used by the transport binding)
 - trimer Cb           <- Schropp ternary_equilibrium(C_cell, R_CD3, R_TAA, KD_CD3, KD_TAA)  (prozone emergent)
 - kill per T cell     <- Cb / (1 + 0.25*n_treg_local)   (Treg suppression, validated tumor constant)
 - cytokine per T cell <- proportional to engaged trimer, hierarchy IL6:IFN:TNF:IL2 (project CYTO_HIER)
 - suppression         <- local Treg density (spatial), damps BOTH kill and cytokine
"""
import numpy as np
from scipy.spatial import cKDTree
import scipy.sparse as sp

AVO = 6.02214076e23
# cytokine hierarchy weights (project-standard CYTO_HIER, mosunetuzumab-anchored)
CYTO_HIER = {"IL6":1.0, "IFN":0.36, "TNF":0.31, "IL2":0.18}
CYTO_IL6_CLINICAL_ANCHOR_PGML = 570.0   # mosunetuzumab peak IL-6 ~570 pg/mL (Hosseini 2020 Fig5A). Anchored by
#   CONSTRUCTION (il6_eng_scale=570/peak_raw) so the scale is dose-invariant; run uses the true C1D1 step-up
#   FIRST dose = 1 mg (CRS worst at first exposure). Same clinical anchor as the well-mixed QSP.
# Acute cytokine DESENSITIZATION (tachyphylaxis) — ported from qsp_costim_window_v2 cytokine calib so the
# whole-body first-dose IL6 spike SELF-LIMITS (~2-3 d) instead of tracking bound drug. Dcyto in [0,1) builds
# under engagement, recovers slowly; cytokine production is gated by resp=(1-Dcyto).
K_CYTO_DESENS = 30.0     # /day desensitization build rate per unit (normalized) engagement
K_CYTO_RECOV  = 0.003    # /day recovery of cytokine responsiveness
DESENS_ENG_REF = 1.0e4   # engagement normalizer (raw eng-sum scale; ~ mosun 1mg peak) so build is O(1)
def cytokine_to_pgml(raw_cyto, il6_eng_scale):
    """Convert a raw cytokine dict (engagement-sum units) to physical pg/mL using the mosunetuzumab-anchored
    il6_eng_scale (pg/mL per raw-IL6-unit). Hierarchy ratios are already baked into raw_cyto via CYTO_HIER."""
    return {k: v*il6_eng_scale for k,v in raw_cyto.items()}
R_SYN_UM = 30.0       # synapse reach (validated tumor value)
# ---- TISSUE MYELOID count_scale (physiological macrophage+monocyte count / sampled myeloid count) ----------
# Citation-gated organ myeloid census, loaded from handoff/organ_myeloid_counts.json if present:
#   {"spleen": 1234.5, "bone": ..., ...}   (drug-INDEPENDENT tissue property; NOT the antigen count_scale)
# Absent -> {} -> every organ falls back to 1.0 (sampled counts). Plasma IL-6 is EXACTLY LINEAR in this scale,
# so a run at 1.0 can be re-scaled analytically from the per-organ il6_prod trace (no re-run required).
def _load_myeloid_scales():
    """Locate organ_myeloid_counts.json ROBUSTLY. This module is imported from EITHER engine/ OR handoff/
    (run_tce_pd_reval.py puts handoff/ at sys.path[0]), so a path relative to __file__ that assumes engine/
    silently resolves to rundir/rundir/handoff/ and finds NOTHING -> every organ falls back to count_scale=1.0
    -> IL-6 comes out ~290,000x too low (measured: 2.38 pg/hr instead of ~5e5). Search several roots and
    ANNOUNCE what was loaded, so a silent miss can never happen again."""
    import json as _json, os as _os, sys as _sys
    _here = _os.path.dirname(_os.path.abspath(__file__))
    _kws = _os.environ.get("KWS_ROOT", "")
    cands = [
        _os.path.join(_here, "organ_myeloid_counts.json"),                                  # handoff/ (import site)
        _os.path.join(_here, "..", "handoff", "organ_myeloid_counts.json"),                 # engine/ -> ../handoff (FINAL_MODEL layout)
        _os.path.join(_here, "..", "rundir", "handoff", "organ_myeloid_counts.json"),       # engine/ -> rundir/handoff (source layout)
        _os.path.join(_here, "..", "..", "rundir", "handoff", "organ_myeloid_counts.json"),
        (_os.path.join(_kws, "handoff", "organ_myeloid_counts.json") if _kws else ""),      # KWS_ROOT/handoff
        "organ_myeloid_counts.json",                                                        # cwd
    ]
    cands = [c for c in cands if c]
    for _c in cands:
        try:
            if _os.path.exists(_c):
                _d = _json.load(open(_c))
                out = {k: float(v) for k, v in _d.items() if isinstance(v, (int, float))}
                _sys.stderr.write(f"[myeloid-IL6] organ census LOADED from {_os.path.abspath(_c)} "
                                  f"({len(out)} organs, e.g. spleen={out.get('spleen')})\n")
                return out
        except Exception:
            pass
    _sys.stderr.write("[myeloid-IL6] *** WARNING: organ_myeloid_counts.json NOT FOUND -> every organ "
                      "count_scale=1.0 -> IL-6 will be ~1e5x TOO LOW. Searched: "
                      f"{[_os.path.abspath(c) for c in cands]}\n")
    return {}
_MYELOID_COUNT_SCALE = _load_myeloid_scales()
R_TREG_UM = 50.0      # Treg suppression neighborhood (validated tumor value)
TREG_K = 0.25         # per-Treg suppression constant (validated tumor value)
RCAP_CD3 = 2.0        # validated CD3 receptor-capacity scale (unified_pbpk_pd Params.pbpk.Rcap_CD3)
RCAP_TAA = 6.0        # validated TAA receptor-capacity scale (unified_pbpk_pd Params.pbpk.Rcap_TAA)
# LITERAL absolute copies->local synapse concentration (nM). ONE physical constant: the synapse reaction
# volume. Pinned by the validated tumor kill anchor: Rcap_TAA=6.0 nM at CEACAM5 257,000 copies/cell =>
# 6.0/257000 = 2.335e-5 nM/copy  (equivalently a ~71 pL synapse reaction volume, V=1e9/(AVO*f)). This
# REPLACES the per-organ-mean normalization (which erased absolute abundance so low-copy BCMA saturated
# like high-copy CD20). Now each cell's arm concentration = copies * NM_PER_COPY, absolute; the ternary
# trimer, prozone, kill and cytokine are then fully emergent from real per-cell binding. Same constant for
# CD3 and TAA arms (CD3 ~92k copies -> ~2.1 nM, recovers the tumor's uniform Rcap_CD3=2.0; in linear regime
# below KD_CD3=40 as CD3 physically is). k_death is re-anchored ONCE against the tumor 28.7% on this basis.
NM_PER_COPY = 6.0/257000.0   # 2.335e-5 nM per receptor copy (71 pL synapse); FIXED for all molecules/targets/organs

def ternary_equilibrium(C, R_A, R_B, KD1, KD2, alpha=1.0):
    """FULL Schropp (2019) closed-form bridged-trimer equilibrium — the SAME implementation the validated
    tumor ABM uses (qsp_costim_window_v2, Eqs 26-30 & 33): free receptors R_A,R_B are solved from the QE
    quadratic BEFORE RC_AB, so the prozone/hook is emergent (at high C both arms saturate as binary
    complexes -> trimer collapses). Vectorized over PER-CELL C, R_A, R_B (each cell its own arms).
    Replaces the earlier reduced-linear form (C*R_A*R_B/KD1KD2) which had NO free-receptor depletion and
    thus NO prozone — that shortcut overstated kill/cytokine at high drug. KDs in nM; returns nM."""
    C=np.maximum(np.asarray(C,float),0.0)
    RA=np.asarray(R_A,float)*np.ones_like(C); RB=np.asarray(R_B,float)*np.ones_like(C)
    aKK=alpha*KD1*KD2
    out=np.zeros_like(C)
    m=(C>0)&(RA>0)&(RB>0)
    c=C[m]; rA=RA[m]; rB=RB[m]
    a=(1.0+c/KD2)*c/aKK                                   # Eq.28
    b=c*(rA-rB)/aKK + (1.0+c/KD2)*(1.0+c/KD1)             # Eq.29
    d=-rB*(1.0+c/KD1)                                      # Eq.30
    disc=b*b-4.0*a*d
    RBf=(-b+np.sqrt(np.maximum(disc,0.0)))/(2.0*a)        # Eq.27 (free R_B)
    RAf=rA/(1.0+c/KD1+RBf*c/aKK)                          # Eq.26 (free R_A)
    out[m]=c*RAf*RBf/aKK                                  # Eq.33
    return out

class OrganPD:
    """Per-cell PD for one organ, built on the SAME cell population as the transport graph.
    Consumes per-cell local drug C (nM) each step; emits per-cell kill hazard, organ cytokine, tracked."""
    def __init__(self, organ, x, y, labs, R_CD3_percell, R_TAA_percell,
                 KD_CD3_nM=40.0, KD_TAA_nM=1.45, costim_boost=0.0,
                 costim_arm=None, R_costim_percell=None, KD_costim_nM=1.0,
                 pd_kinetics=False, kon_CD3_perM_s=1e5, kon_TAA_perM_s=1e5,
                 koff_CD3_pers=None, koff_TAA_pers=None, kint_bridge_perday=0.9,
                 span_bridge_nm=12.5, span_cis_nm=12.5, cis_avidity=0.0,
                 n_CD3=1, n_TAA=1, n_costim=0, span_coeng_T_nm=None, span_coeng_tumor_nm=None,
                 k_hit_perday=12.0):
        self.organ=organ
        self.x=np.asarray(x,float); self.y=np.asarray(y,float)
        self.labs=np.array([str(v) for v in labs])
        # ---- MECHANISTIC IL-6 (replaces the fitted mosun-anchored IL6_SCALE) --------------------------
        # IL-6 is MYELOID-derived (Giavridis PMID 29808005; Norelli PMID 29808007), contact-activated by
        # ENGAGED T cells. Each monocyte/macrophage agent is its own emitter at its MEASURED per-cell rate
        # (0.0196 pg/hr/cell, PMID 37533643). Saturation + per-drug differences EMERGE from the finite,
        # spatially-distributed myeloid pool — no Emax, no EC50, no fitted scale.
        # CONTACT distance (14.1 um = r_macrophage 10.6 + r_Tcell 3.5; PMID 9400735 / PMID 30571054), NOT the
        # 30 um T:target synapse REACH. CD40L-CD40 is membrane-bound -> the cells must physically touch.
        from myeloid_il6 import MyeloidIL6 as _MyeIL6, R_CONTACT_UM as _RC
        self.myeloid = _MyeIL6(self.x, self.y, self.labs, r_contact_um=_RC)
        self.il6_prod_pg_hr = 0.0          # this organ's instantaneous myeloid IL-6 production (pg/hr)
        # MYELOID count_scale: sampled myeloid -> physiological myeloid. MUST be a TISSUE property
        # (drug-INDEPENDENT). Do NOT reuse graphs[o].count_scale: that is antigen-derived and VERIFIED
        # drug-dependent (0.18x-5.34x between MS4A1 and TNFRSF17 in the same organ, measured 2026-07-13) ->
        # it would scale the identical monocyte population differently per drug = counterscreen-corrupting.
        # BLOOD is handled separately by blood_count_scale (real circulating count / sampled count = correct,
        # cell-based, per-lineage). Solid-tissue myeloid cellularity is the ONE remaining DATA GATE:
        # default 1.0 (sampled counts) until an organ macrophage/monocyte census is sourced. This affects the
        # ABSOLUTE pg/mL only; the per-drug ORDERING and all relative comparisons are already correct.
        # Loaded from handoff/organ_myeloid_counts.json when present (citation-gated census: physiological
        # macrophage+monocyte count per organ / sampled myeloid count). Absent -> 1.0 (sampled counts), which
        # affects the ABSOLUTE pg/mL only; per-drug ORDERING is already correct and drug-independent.
        # Because plasma IL-6 is EXACTLY LINEAR in this scale, a run made at 1.0 can be re-scaled analytically
        # from the per-organ il6_prod trace -- no re-run needed when the census lands.
 # A COUNT IS NOT A SCALE. (bug found + fixed 2026-07-13 by a physical-ceiling check.)
        # organ_myeloid_counts.json holds RAW PHYSIOLOGICAL COUNTS (bone = 33,633,218 macrophages). The code
        # was using those counts DIRECTLY as a multiplicative scale, i.e. multiplying this organ's production
        # by the whole macrophage population instead of by (population / sampled agents). That inflates IL-6
        # by exactly n_sampled_myeloid -- measured 942x above the model's own physical ceiling (every secretor
        # macrophage in the census, fully activated, at the measured per-cell secretion rate, can make at most
        # 1,821 pg/hr; the model was reporting 1,716,379 pg/hr).
        #
        # The scale is: how many PHYSIOLOGICAL myeloid cells does ONE SAMPLED AGENT stand for?
        #     scale = census_count(organ) / n_sampled_myeloid_agents(organ)
        _cen = _MYELOID_COUNT_SCALE.get(organ)
        _n_sampled = max(int(getattr(self.myeloid, "n_myeloid", 0)), 1)
        if _cen is None or _n_sampled == 0:
            self.myeloid_count_scale = 1.0            # no census for this organ -> sampled counts
        else:
            self.myeloid_count_scale = float(_cen) / float(_n_sampled)
        self.myeloid_census_count = float(_cen) if _cen is not None else None
        self.myeloid_n_sampled = _n_sampled
        self.KD_CD3=KD_CD3_nM; self.KD_TAA=KD_TAA_nM; self.costim_boost=costim_boost
        # kinetic-synapse config (literal engage/hit/detach ODE, user 2026-07-12). QSS path is default.
        self.pd_kinetics=bool(pd_kinetics)
        self._kin_cfg=dict(kon_CD3_perM_s=kon_CD3_perM_s, kon_TAA_perM_s=kon_TAA_perM_s,
                           koff_CD3_pers=koff_CD3_pers, koff_TAA_pers=koff_TAA_pers,
                           kint_bridge_perday=kint_bridge_perday, span_bridge_nm=span_bridge_nm,
                           span_cis_nm=span_cis_nm, cis_avidity=cis_avidity, k_hit_perday=k_hit_perday)
        self.kin=None                                   # built after neighborhoods (needs Wt_norm)
        n=len(self.x)
        # classify cells
        low=np.char.lower(self.labs)
        self.is_T   = np.array([('t cell' in s or 't_cell' in s or 'cd8' in s or 'cd4' in s or
                                  'regulatory t' in s or 'nk t' in s) for s in low])
        self.is_CD8 = np.array([('cd8' in s) for s in low])
        self.is_treg= np.array([('regulatory t' in s or 'treg' in s) for s in low])
        self.is_target = (np.asarray(R_TAA_percell,float) > 0)   # any antigen-bearing cell is killable
        R_CD3_raw = np.asarray(R_CD3_percell,float)
        R_TAA_raw = np.asarray(R_TAA_percell,float)
        # LITERAL ABSOLUTE conversion: each cell's arm concentration (nM) = its own receptor copies x the
        # fixed synapse conversion NM_PER_COPY. NO per-organ-mean normalization — absolute abundance is
        # preserved, so genuinely low-copy targets (BCMA ~11k -> 0.26 nM, below KD_TAA) form a weaker bridge
        # than high-copy targets (CD20 ~95k -> 2.2 nM) and their depletion becomes exposure-dependent rather
        # than saturating identically. The trimer/prozone/kill are then emergent from real per-cell binding.
        self.R_CD3 = R_CD3_raw * NM_PER_COPY
        self.R_TAA = R_TAA_raw * NM_PER_COPY
        self.R_CD3_raw = R_CD3_raw; self.R_TAA_raw = R_TAA_raw
        # cumulative per-cell kill hazard (integral of kill rate * dt); survival = exp(-hazard)
        self.kill_hazard = np.zeros(n)
        self.alive = np.ones(n, bool)
        # precompute synapse neighborhoods: each T cell -> target cells within R_SYN; T cell -> Tregs within R_TREG
        self._build_neighborhoods()
        # cumulative cytokine (organ-summed, per cytokine species)
        self.cyto = {k:0.0 for k in CYTO_HIER}          # cumulative integral (AUC-like)
        self.cyto_rate = {k:0.0 for k in CYTO_HIER}     # INSTANTANEOUS rate (unit-comparable to clinical peak)
        self.Dcyto = 0.0                                # acute desensitization state [0,1)
        self._killfrac = 0.0
        # ---- PER-CELL COSTIM SIGNALING (replaces the scalar costim_boost) --------------------------
        # Each T cell / Treg carries its OWN costim-receptor copies (R_costim_percell, HPA/Glassman).
        # The costim arm engages at occupancy occ = C/(C+KD_costim) scaled by that cell's receptor copies
        # (relative to the population anchor). Occupancy drives a per-cell, per-program signaling state
        # (signaling_dynamics.PerCellSignaling): magnitude from the lane drive table, timescale from the
        # hero Rest/8hr/48hr kinetics, per-cell heterogeneity from receptor-copy spread. The CD8/effector
        # cells read effector/cytokine/exhaustion; Tregs read suppression. This is the QSP efficacy+tox core.
        self.costim_arm=costim_arm
        self.sig=None; self.R_costim=None; self.p_cis=0.0; self.n_costim=int(n_costim)
        self.costim_ind=None            # activation-induced receptor density (None => static resting, legacy)
        if costim_arm is not None and R_costim_percell is not None:
            import signaling_dynamics as _sigmod
            self._sigmod=_sigmod
            Rc=np.asarray(R_costim_percell,float)
            # .copy() is REQUIRED: activation-induction WRITES into self.R_costim each step. np.asarray does
            # not copy an existing float array, so without this every organ would mutate the SHARED caller
            # array and leak its induced densities into the others.
            self.R_costim=Rc.copy(); self.KD_costim=KD_costim_nM
            # per-cell receptor scale for the T-cell population (occupancy amplitude multiplier)
            Tpos=Rc[self.Tidx] if len(self.Tidx)>0 else np.array([1.0])
            self._costim_anchor=Tpos[Tpos>0].mean() if (Tpos>0).any() else 1.0
            # --- FORMAT geometry for the costim arm (multiarm_binding): cis coincidence fraction ---
            # p_cis in [0,1] from the T-side co-engagement span (height-matched->1 cis; 60A mismatch->0 trans).
            # DEFAULT (span_coeng_T_nm=None) -> p_cis=0 => costim drive stays CELL-AUTONOMOUS (trans / legacy).
            self.n_costim=int(n_costim); self.n_CD3=int(n_CD3); self.n_TAA=int(n_TAA)
            self.span_coeng_T_nm=span_coeng_T_nm; self.span_coeng_tumor_nm=span_coeng_tumor_nm
            try:
                from multiarm_binding import _cis_feasibility as _cisf
                self.p_cis = _cisf(span_coeng_T_nm) if (self.n_costim>0 and span_coeng_T_nm is not None) else 0.0
            except Exception:
                self.p_cis = 0.0
            # one signaling integrator over the T-cell population (vectorized); Treg mask within it
            self.sig=_sigmod.PerCellSignaling(len(self.Tidx), costim_arm, drive_source="lane")
            self._is_treg_T = self.is_treg[self.Tidx]      # which of the T-cell agents are Tregs
            self._is_cd8_T  = self.is_CD8[self.Tidx]
            # ---- ACTIVATION-INDUCED costim receptor density (4-1BB/OX40/ICOS/GITR are induced, not resting) ----
            # ON for every arm and every organ. Constitutive arms carry fold=1.0 -> R(t)==R_rest -> the old
            # static path exactly. Inducible arms with an UNSOURCED fold RAISE rather than silently run at the
            # resting density (which would under-rate them by the induction factor — the exact bug this fixes).
            try:
                from costim_induction import CostimInduction as _CI
                self.costim_ind = _CI(costim_arm, len(self.Tidx), Rc[self.Tidx], strict=True)
            except ImportError:
                self.costim_ind = None
        # ---- build the KINETIC synapse engine (if enabled) — needs Wt_norm from _build_neighborhoods ----
        self.per_target_surv=None
        if self.pd_kinetics and len(self.Tidx)>0 and self.tgtidx.size>0:
            import kinetic_synapse as _ks
            c=self._kin_cfg
            SPD=86400.0
            kon_CD3=c['kon_CD3_perM_s']*SPD/1e9        # /M/s -> /nM/day
            kon_TAA=c['kon_TAA_perM_s']*SPD/1e9
            # koff from literature if given (per s -> per day); else koff = kon*KD (KD-consistent split)
            koff_CD3=(c['koff_CD3_pers']*SPD) if c['koff_CD3_pers'] else kon_CD3*self.KD_CD3
            koff_TAA=(c['koff_TAA_pers']*SPD) if c['koff_TAA_pers'] else kon_TAA*self.KD_TAA
            # per-T CD3 capacity (nM) and per-target TAA capacity (nM): reuse the SAME absolute R arrays
            R_CD3_T = self.R_CD3[self.Tidx]
            R_TAA_tgt = self.R_TAA[self.tgtidx]
            # per-T mean neighbour TAA density (copies) for the geometric ceff: recover copies from nM
            dens_copies_T = (self.syn_TAA_mean/NM_PER_COPY)
            self.kin=_ks.KineticSynapse(R_CD3_T, R_TAA_tgt, self.Wt_norm, dens_copies_T,
                kon_CD3, koff_CD3, kon_TAA, koff_TAA, c['kint_bridge_perday'],
                span_bridge_nm=c['span_bridge_nm'], span_cis_nm=c['span_cis_nm'],
                cis_avidity=c['cis_avidity'], k_hit_perday=c['k_hit_perday'],
                W_incidence=self.W)
            # free-CD3 receptor turnover (KSYN=RC0*kdeg): kdeg from cfg if provided, else 0 (static, back-compat)
            self.kin.set_turnover(float(c.get('kdeg_CD3_perday', 0.0)))
            self.per_target_surv=np.ones(self.tgtidx.size)

    def _build_neighborhoods(self):
        Tidx=np.where(self.is_T)[0]; tgtidx=np.where(self.is_target)[0]; tregidx=np.where(self.is_treg)[0]
        self.Tidx=Tidx; self.tgtidx=tgtidx
        if len(Tidx)==0:
            self.syn_nb=[]; self.n_treg=np.zeros(0); return
        xyT=np.c_[self.x[Tidx], self.y[Tidx]]
        if len(tgtidx)>0:
            tree_t=cKDTree(np.c_[self.x[tgtidx], self.y[tgtidx]])
            self.syn_nb=tree_t.query_ball_point(xyT, r=R_SYN_UM)     # per T cell: list of local target (row in tgtidx)
            # sparse synapse incidence T x n_target (1 where target in reach); used for vectorized apportionment
            rows=[]; cols=[]
            for j,nb in enumerate(self.syn_nb):
                for c in nb: rows.append(j); cols.append(c)
            self.W = sp.csr_matrix((np.ones(len(rows)), (rows, cols)),
                                   shape=(len(Tidx), len(tgtidx))) if rows else sp.csr_matrix((len(Tidx),len(tgtidx)))
            # per-T synapse TAA weight vector prep: each edge weighted by target antigen; row-normalized at use
            self.tgt_taa = self.R_TAA[tgtidx]                        # antigen per target cell (col-aligned)
            Wt = self.W.multiply(self.tgt_taa[None,:]).tocsr()       # edge weight = target antigen
            rs = np.asarray(Wt.sum(1)).ravel(); rs[rs==0]=1.0
            self.Wt_norm = sp.diags(1.0/rs) @ Wt                     # row-normalized apportionment matrix (T x target)
            self.syn_TAA_mean = np.asarray(self.W.multiply(self.tgt_taa[None,:]).sum(1)).ravel() / np.maximum(np.asarray(self.W.sum(1)).ravel(),1.0)
        else:
            self.syn_nb=[[] for _ in Tidx]
            self.W=sp.csr_matrix((len(Tidx),0)); self.Wt_norm=sp.csr_matrix((len(Tidx),0))
            self.tgt_taa=np.zeros(0); self.syn_TAA_mean=np.zeros(len(Tidx))
        if len(tregidx)>0:
            tree_r=cKDTree(np.c_[self.x[tregidx], self.y[tregidx]])
            self.n_treg=np.array([len(v) for v in tree_r.query_ball_point(xyT, r=R_TREG_UM)],float)
        else:
            self.n_treg=np.zeros(len(Tidx))


    def _apply_cis_coincidence(self, occ, cd3_engagement):
        """CIS/TRANS coincidence gating of costim occupancy (format-driven, per-cell).
        occ_eff = occ * [(1-p_cis) + p_cis * f_cd3], where f_cd3 is the per-T normalized CD3-engagement
        (QSS trimer Cb, or kinetic bridged B2). p_cis=0 (trans / no-costim / legacy default) -> occ
        UNCHANGED (cell-autonomous). p_cis=1 (cis, height-matched) -> costim fully gated on the SAME
        cell's CD3 engagement: the coincident-signal design. Coincidence is EMERGENT from span geometry
        (p_cis, Rhoden reach) x real per-cell CD3 binding, not a phenomenological constant."""
        p = getattr(self, "p_cis", 0.0)
        if p <= 0.0:
            return occ                                   # trans / no costim arm -> autonomous (byte-identical)
        e = np.maximum(np.asarray(cd3_engagement, float), 0.0)
        ref = np.median(e[e > 0]) if np.any(e > 0) else 1.0   # per-cell normalizer (self-scaling)
        f_cd3 = e / (e + max(ref, 1e-12))                # normalized CD3-engagement [0,1)
        return occ * ((1.0 - p) + p * f_cd3)

    def step(self, C_percell, dt, k_death=1.0):
        """One PD step. C_percell = per-cell local free drug (nM) from the transport grid (length n).
        k_death (1/day) = trimer->death rate constant, calibrated to the validated tumor 28.7% (calibrate_kdeath).
        Returns dict(kill_frac, cyto_rate, n_engaged). Accumulates DRUG-GRADED hazard + cytokine."""
        Tidx=self.Tidx
        if len(Tidx)==0:
            return dict(kill_frac=0.0, cyto_rate={k:0.0 for k in CYTO_HIER})
        # ---- KINETIC-SYNAPSE PATH (literal engage/hit/detach ODE) --------------------------------------
        if self.pd_kinetics and self.kin is not None:
            return self._step_kinetic(C_percell, dt, k_death)
        Cd=np.maximum(C_percell[Tidx],0.0)
        RA=self.R_CD3[Tidx]
        RB=self.syn_TAA_mean                                        # precomputed per-T synapse antigen mean
        # per-T-cell trimer (Schropp ternary, prozone emergent). R_CD3/R_TAA are brought onto the VALIDATED
        # tumor Rcap basis (self.R_CD3_cap, self.R_TAA_cap) so the whole-body kill scale matches the
        # tumor's validated 28.7% snapshot rather than raw copy magnitudes.
        Cb=ternary_equilibrium(Cd, RA, RB, self.KD_CD3, self.KD_TAA)
        # ---- per-cell costim signaling (replaces scalar costim_boost) --------------------------------
        # advance each T cell's signaling state by its OWN costim-receptor occupancy at local drug, then
        # read per-cell program modifiers: effector -> kill gain, suppression -> extra Treg-like damping,
        # exhaustion -> kill decay. Legacy scalar path preserved when no arm wired.
        g_eff=1.0; supp_extra=0.0
        if self.sig is None and self.costim_boost>0:
            g_eff=(1.0+self.costim_boost)         # LEGACY scalar path preserved when no per-cell arm wired
        if self.sig is not None:
            # ---- ACTIVATION-INDUCED COSTIM RECEPTOR DENSITY -------------------------------------------
            # 4-1BB / OX40 / ICOS / GITR are ABSENT on resting T cells and appear only AFTER TCR engagement.
            # That conditionality is the entire reason the field targets them. Reading a STATIC resting
            # density (the old behaviour: self.R_costim set once in __init__, never updated) sees them at
            # ~zero and systematically UNDER-RATES exactly the arms that matter, while OVER-RATING the
            # constitutive ones (CD28/CD2/CD27) -> a confidently-wrong ranking ("CD2 beats 4-1BB").
            #
            # Induction is driven by THIS T cell's OWN engaged trimer -> a T cell that never engages never
            # upregulates. The tumour-conditionality is therefore EMERGENT, not imposed.
            # Constitutive arms have fold=1.0 -> R(t) == R_rest -> byte-identical to the old path.
            if self.costim_ind is not None:
                p_eng = Cb/(Cb + RA + 1e-12)                       # per-T-cell engaged fraction in [0,1)
                self.R_costim[Tidx] = self.costim_ind.step(dt, p_eng)
            Rc_T=self.R_costim[Tidx]
            occ = (Cd/(Cd+self.KD_costim)) * (Rc_T/max(self._costim_anchor,1e-9))   # per-cell occupancy [0,~]
            occ = self._apply_cis_coincidence(occ, Cb)   # cis: gate costim on per-cell CD3 trimer; trans: unchanged
            occ=np.clip(occ,0.0,5.0)
            self.sig.step(occ, dt*24.0)                              # signaling kinetics are per-HOUR; dt is days
            eff_p=self.sig.program("effector"); supp_p=self.sig.program("suppression")
            exh_p=self.sig.program("exhaustion")
            # effector program -> multiplicative kill gain on CD8/effector T cells (bounded, agonism>0 => >1)
            kE_gain=0.55                                             # effector->kill sensitivity (locked calib)
            g_eff=np.ones(len(Tidx))
            g_eff[self._is_cd8_T] = np.exp(kE_gain*eff_p[self._is_cd8_T])
            # exhaustion program -> attenuates kill (agonism that raises exhaustion loses durable killing)
            g_eff = g_eff*np.exp(-0.30*np.maximum(exh_p,0.0))
            # suppression program on Tregs -> extra suppression of the whole synapse neighborhood
            supp_extra = float(np.mean(np.maximum(supp_p[self._is_treg_T],0.0))) if self._is_treg_T.any() else 0.0
        Cb = Cb*g_eff                                               # per-cell effector-gained trimer
        kill_T = Cb/(1.0+TREG_K*self.n_treg*(1.0+supp_extra))       # Treg-damped + signaling-suppression per-T kill
        # VECTORIZED apportionment: dkill_target = Wt_norm^T @ kill_T  (antigen-weighted, row-normalized)
        n=len(self.x); dkill=np.zeros(n)
        dkill_tgt = self.Wt_norm.T @ kill_T                         # length n_target
        dkill[self.tgtidx] = dkill_tgt
        # accumulate hazard as a DRUG-GRADED RATE: the apportioned trimer (dkill) is the instantaneous
        # per-target kill propensity (the SAME species the validated tumor ABM reads as 1-exp(-Cb)).
        # Integrating dt*kkill*dkill makes the cumulative killed fraction depend on BOTH drug level
        # (via dkill) AND time — so low exposure -> low rate -> incomplete kill within the window
        # (drug-graded plateau), unlike a saturating fixed hazard. kkill is calibrated against the
        # validated tumor 28.7% (see calibrate_kkill), NOT the tumor's snapshot constant.
        self.kill_hazard += dt*k_death*dkill
        # cytokine: each engaged T cell emits proportional to its trimer, hierarchy-weighted, Treg-damped
        eng = float((Cb/(1.0+TREG_K*self.n_treg)).sum())
        # signaling cytokine drive: costim engagement scales cytokine emission per its IFNG/TNF/IL2 programs
        cyto_sig_gain=1.0
        if self.sig is not None:
            ifn=self.sig.program("cyto_IFNG"); tnf=self.sig.program("cyto_TNF"); il2=self.sig.program("cyto_IL2")
            # mean over engaged T cells; agonism-positive raises the storm, cold arms lower it (bounded >=0.2)
            cs=float(np.mean(0.45*ifn+0.32*tnf+0.18*il2))
            cyto_sig_gain=max(0.2, 1.0+cs)
        # DESENSITIZATION gate: cytokine responsiveness resp=(1-Dcyto); the storm self-limits as Dcyto builds
        resp = max(1.0 - self.Dcyto, 0.0)
        cyto_rate={k: w*eng*resp*cyto_sig_gain for k,w in CYTO_HIER.items()}   # gated + signaling-scaled rate
        for k in self.cyto: self.cyto[k]+= dt*cyto_rate[k]
        self.cyto_rate = dict(cyto_rate)                # keep the (gated) instantaneous rate for peak readout
        # MECHANISTIC IL-6 (QSS mirror of the kinetic path): per-T engagement -> per-cell myeloid activation.
        # QSS has no B2; the per-T trimer Cb saturates to its own receptor pool -> engaged fraction Cb/(Cb+RA).
        p_eng = np.clip(Cb/np.maximum(RA, 1e-30), 0.0, 1.0)
        self.il6_prod_pg_hr = self.myeloid.step(dt, self.x[Tidx], self.y[Tidx], p_eng)
        # evolve Dcyto: builds under engagement (normalized), recovers slowly (both bounded to [0,1))
        eng_norm = eng/DESENS_ENG_REF
        self.Dcyto += dt*(K_CYTO_DESENS*eng_norm*resp - K_CYTO_RECOV*self.Dcyto)
        self.Dcyto = min(max(self.Dcyto,0.0),0.999)
        killfrac=float((1-np.exp(-self.kill_hazard[self.is_target])).mean()) if self.is_target.any() else 0.0
        self._killfrac=killfrac   # cache for O(1) summary()
        return dict(kill_frac=killfrac, cyto_rate=cyto_rate, n_engaged=int((Cb>1e-6).sum()))

    def move_immune(self, dt_days, speed_um_per_min=5.0, chemotax=0.4, rebuild=True):
        """Move motile immune cells (T cells) over dt with a persistent random walk + weak chemotactic bias
        toward antigen-dense (target) neighbourhoods, then rebuild the synapse graph so NEW contacts form.
        T cells migrate ~2-10 um/min in tissue (Miller 2002/2004 two-photon); default 5 um/min. Bounded to
        the tissue extent. chemotax in [0,1] = fraction of the step biased up the local target gradient.
        rebuild=True re-runs _build_neighborhoods (KDTree) so W, syn_nb, Wt_norm reflect the new positions;
        the kinetic synapse engine is rebuilt to match. Non-motile cells (targets, stroma) stay put."""
        Tidx=self.Tidx
        if len(Tidx)==0: return
        step_um = speed_um_per_min*60.0*24.0*dt_days   # um per this move interval
        rng=getattr(self,'_rng',None)
        if rng is None: self._rng=rng=np.random.RandomState(0)
        # persistent random-walk heading (retain a per-T heading, nudge it) + chemotactic bias
        th=getattr(self,'_theta',None)
        if th is None or len(th)!=len(Tidx): self._theta=th=rng.uniform(0,2*np.pi,len(Tidx))
        th = th + rng.normal(0,0.5,len(Tidx))            # heading persistence with angular diffusion
        self._theta=th
        dx=np.cos(th)*step_um; dy=np.sin(th)*step_um
        # chemotactic bias: pull toward local target centroid within reach (if any targets)
        if chemotax>0 and self.tgtidx.size>0:
            xT=self.x[Tidx]; yT=self.y[Tidx]
            # vector to mean of currently-reachable targets (from existing syn_nb); fall back to global centroid
            gx=np.mean(self.x[self.tgtidx]); gy=np.mean(self.y[self.tgtidx])
            vx=gx-xT; vy=gy-yT; vn=np.hypot(vx,vy); vn[vn==0]=1.0
            dx=(1-chemotax)*dx + chemotax*step_um*vx/vn
            dy=(1-chemotax)*dy + chemotax*step_um*vy/vn
        # apply + clamp to tissue bounds
        self.x[Tidx]=np.clip(self.x[Tidx]+dx, self.x.min(), self.x.max())
        self.y[Tidx]=np.clip(self.y[Tidx]+dy, self.y.min(), self.y.max())
        if rebuild:
            _B1=getattr(self.kin,'B1',None) if getattr(self,'kin',None) is not None else None
            self._build_neighborhoods()                  # W, syn_nb, Wt_norm, syn_TAA_mean, n_treg refreshed
            # rebuild kinetic synapse against the new neighbourhoods (preserve bound state where shapes match)
            if getattr(self,'pd_kinetics',False) and getattr(self,'kin',None) is not None and self.tgtidx.size>0:
                import kinetic_synapse as _ks
                c=self._kin_cfg; SPD=86400.0
                kon_CD3=c['kon_CD3_perM_s']*SPD/1e9; kon_TAA=c['kon_TAA_perM_s']*SPD/1e9
                koff_CD3=(c['koff_CD3_pers']*SPD) if c['koff_CD3_pers'] else kon_CD3*self.KD_CD3
                koff_TAA=(c['koff_TAA_pers']*SPD) if c['koff_TAA_pers'] else kon_TAA*self.KD_TAA
                dens=(self.syn_TAA_mean/NM_PER_COPY)
                newkin=_ks.KineticSynapse(self.R_CD3[Tidx], self.R_TAA[self.tgtidx], self.Wt_norm, dens,
                    kon_CD3,koff_CD3,kon_TAA,koff_TAA,c['kint_bridge_perday'],
                    span_bridge_nm=c['span_bridge_nm'],span_cis_nm=c['span_cis_nm'],
                    cis_avidity=c['cis_avidity'],k_hit_perday=c['k_hit_perday'],W_incidence=self.W)
                newkin.set_turnover(float(c.get('kdeg_CD3_perday',0.0)))
                if _B1 is not None and len(_B1)==len(Tidx):      # carry CD3-side bound state (T-indexed)
                    newkin.B1=_B1; newkin.B2=self.kin.B2
                self.kin=newkin

    def _step_kinetic(self, C_percell, dt, k_death):
        """Kinetic-synapse PD step: route through the per-cell engage/hit/detach ODE (kinetic_synapse).
        Per-target survival evolves so serial killing + avidity dwell-time are emergent. Cytokine is
        driven by ENGAGED DWELL (synapse-stability), signaling-scaled the same way as the QSS path."""
        Tidx=self.Tidx
        Cd_T=np.maximum(C_percell[Tidx],0.0)
        # per-cell costim signaling (identical to QSS path): occupancy -> program modifiers
        g_eff=1.0; supp_extra=0.0
        if self.sig is not None:
            Rc_T=self.R_costim[Tidx]
            occ=(Cd_T/(Cd_T+self.KD_costim))*(Rc_T/max(self._costim_anchor,1e-9))
            occ=self._apply_cis_coincidence(occ, self.kin.B2)   # cis: gate costim on per-cell CD3 engagement (B2); trans: unchanged
            occ=np.clip(occ,0.0,5.0)
            self.sig.step(occ, dt*24.0)
            eff_p=self.sig.program("effector"); supp_p=self.sig.program("suppression"); exh_p=self.sig.program("exhaustion")
            kE_gain=0.55
            g_eff=np.ones(len(Tidx)); g_eff[self._is_cd8_T]=np.exp(kE_gain*eff_p[self._is_cd8_T])
            g_eff=g_eff*np.exp(-0.30*np.maximum(exh_p,0.0))
            supp_extra=float(np.mean(np.maximum(supp_p[self._is_treg_T],0.0))) if self._is_treg_T.any() else 0.0
        # advance the kinetic synapse: returns per-target hazard increment (already dt*k_death*serial_rate)
        # effector gain scales k_death (more effector -> more lethal per engaged synapse); Treg damping
        # applied to the per-T serial output via a global suppression factor (matches QSS TREG_K structure).
        k_death_eff=k_death
        dH_tgt=self.kin.step(Cd_T, dt, k_death_eff, self.per_target_surv)
        # apply effector/Treg modifiers as a scalar on the hazard increment (per-target already apportioned)
        gscale=float(np.mean(g_eff)) if np.ndim(g_eff) else g_eff
        treg_damp=1.0/(1.0+TREG_K*float(np.mean(self.n_treg))*(1.0+supp_extra)) if len(self.n_treg) else 1.0
        dH_tgt=dH_tgt*gscale*treg_damp
        # accumulate hazard on the TARGET cells (map tgtidx -> full cell index)
        self.kill_hazard[self.tgtidx]+=dH_tgt
        self.per_target_surv=np.exp(-self.kill_hazard[self.tgtidx])
        # cytokine from ENGAGED DWELL (synapse-stability): more sustained engagement -> more cytokine.
        # eng ~ total engaged synapse level (sum p_eng), Treg-damped; signaling-scaled like QSS path.
        eng=self.kin.engaged_dwell_rate()*treg_damp
        cyto_sig_gain=1.0
        if self.sig is not None:
            ifn=self.sig.program("cyto_IFNG"); tnf=self.sig.program("cyto_TNF"); il2=self.sig.program("cyto_IL2")
            cs=float(np.mean(0.45*ifn+0.32*tnf+0.18*il2)); cyto_sig_gain=max(0.2,1.0+cs)
        resp=max(1.0-self.Dcyto,0.0)
        cyto_rate={k:w*eng*resp*cyto_sig_gain for k,w in CYTO_HIER.items()}
        for k in self.cyto: self.cyto[k]+=dt*cyto_rate[k]
        self.cyto_rate=dict(cyto_rate)
        # ---- MECHANISTIC IL-6: per-cell myeloid emitters, contact-activated by THIS organ's engaged T cells.
        # p_eng = B2/RC = each T cell's own engaged-synapse fraction (the CD40L-bearing, activated T cells).
        # Every myeloid agent integrates ITS OWN local engaged-T contact -> its own activation -> its own
        # secretion at the measured per-cell rate. Organ IL-6 production (pg/hr) is summed into plasma by
        # CoupledPerCellPD against the measured IL-6 clearance. NOTHING here is fitted.
        p_eng = np.clip(self.kin.B2/np.maximum(self.kin.RC,1e-30), 0.0, 1.0)
        self.il6_prod_pg_hr = self.myeloid.step(dt, self.x[Tidx], self.y[Tidx], p_eng)
        eng_norm=eng/DESENS_ENG_REF
        self.Dcyto=min(max(self.Dcyto+dt*(K_CYTO_DESENS*eng_norm*resp-K_CYTO_RECOV*self.Dcyto),0.0),0.999)
        killfrac=float((1-np.exp(-self.kill_hazard[self.is_target])).mean()) if self.is_target.any() else 0.0
        self._killfrac=killfrac
        return dict(kill_frac=killfrac, cyto_rate=cyto_rate, n_engaged=int((self.kin.B2>1e-6).sum()))

    def summary(self):
        # O(1) read of the values step() already computed — no per-record recompute over all cells
        return dict(organ=self.organ, kill_frac=self._killfrac, cyto=dict(self.cyto),
                    cyto_rate=dict(self.cyto_rate))    # instantaneous (for clinical-peak anchoring)

    def summary_full(self):
        surv=np.exp(-self.kill_hazard)
        return dict(organ=self.organ, n_cells=len(self.x), n_T=int(self.is_T.sum()),
                    n_target=int(self.is_target.sum()), n_treg=int(self.is_treg.sum()),
                    kill_frac=float((1-surv[self.is_target]).mean()) if self.is_target.any() else 0.0,
                    cyto=dict(self.cyto))
