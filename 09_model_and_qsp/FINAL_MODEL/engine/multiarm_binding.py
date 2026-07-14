"""multiarm_binding.py — generalized per-cell multivalent engager binding, FORMAT-DRIVEN geometry.

Three arms, each valency n in {0,1,2}: CD3 (T cell), COSTIM (T cell), TAA (tumor cell).
THREE independent geometric spans encode the CONSTRUCT FORMAT (BiTE / IgG / 2+1 / tetravalent C-term
fusion / ...). Each is settable per construct (AF3-derived override) or defaulted by format:

  span_bridge_nm       CROSS-CELL: T-cell-binding module <-> tumor-binding module. Sets synapse-cleft
                       feasibility + cross-cell bridge occupancy -> KILLING (CD3xTAA) and signal-2
                       delivery (costimxTAA).
  span_coeng_T_nm      T-SIDE co-engagement: between two binders on the T-cell end. Governs BOTH
                       (a) CD3<->costim CIS coincidence (one molecule co-engaging both on one T cell)
                       and (b) bivalent-CD3 / bivalent-costim avidity on one T cell.
  span_coeng_tumor_nm  TUMOR-SIDE co-engagement: between two binders on the tumor end. Governs
                       bivalent-TAA avidity AND dual-TAA co-engagement (two TAAs on one tumor cell).

Bound species returned (these DRIVE SIGNALS downstream):
  Cb_kill    CD3.drug.TAA bridged trimer  -> kill + CRS cytokine (signal 1)
  Cb_costim  costim-engaged drug          -> costim program (effector/supp/exh)
  Cb_cis     costim COINCIDENT with CD3 on the same T cell (cis subset) -> coincident signal 1+2
  Cb_costTAA costim.drug.TAA bridge (signal-2, TAA-anchored, CD3-independent)
One physics: the SAME Rhoden-2016 geometric ageff used by percell_binding / kinetic_synapse.
"""
import numpy as np
AVO=6.02214076e23

# arm reach when a span is unset (compact within-module reach, nm)
DEFAULT_ARM_REACH_NM = 12.5

def geo_ageff_nM(rec_pc, r_cell_um, span_nm):
    """Rhoden geometric effective 2nd-arm conc (nM); IDENTICAL to qsp _geo_ageff / percell ageff_nM.
    Larger span -> larger explored shell -> LOWER c_eff (dilution); shorter -> higher."""
    _cc=1e9
    SA_cell=4.0*np.pi*r_cell_um**2; r_Ab_um=max(span_nm,1e-3)*1e-3
    SA_Ab=np.pi*r_Ab_um**2; V_Ab=(2.0/3.0)*np.pi*r_Ab_um**3
    Ag_bulk=np.asarray(rec_pc,float)*_cc/AVO*1e9
    Am_cell=Ag_bulk/_cc; Am_SA=Am_cell/SA_cell; Ag_SA=Am_SA*SA_Ab
    return Ag_SA/V_Ab*1e15

def _cis_feasibility(span_coeng_nm, gap_match_nm=12.5, tol_nm=8.0):
    """CIS co-engagement feasibility [0,1]: can ONE molecule co-engage two binders on the SAME cell
    surface? Feasible when the construct co-engagement span MATCHES the inter-epitope gap (height-
    matched). A deliberate epitope-HEIGHT mismatch (e.g. tall CRD1 4-1BB, ~60A) forces span far from
    gap_match -> arms cannot co-reach -> ->0 (TRANS design). Gaussian tolerance in span."""
    if span_coeng_nm is None: return 0.0
    d=(float(span_coeng_nm)-gap_match_nm)/max(tol_nm,1e-6)
    return float(np.exp(-0.5*d*d))

def _bridge_feasibility(span_bridge_nm, cleft_min=13.0, cleft_max=40.0):
    """Cross-cell bridge feasibility [0,1]: the T-arm<->tumor-arm span must hold the ~13-40nm immune
    synapse cleft. Same convention as kinetic_synapse.cleft_feasibility (span/cleft ramps to 1)."""
    cleft=float(np.clip(span_bridge_nm, cleft_min, cleft_max))
    g=span_bridge_nm/max(cleft,1e-6)
    return float(np.clip(g,0.0,1.0))

def bound_arm(n_arm, R_local, C, rec_pc, r_cell_um, KD, coeng_span_nm=None):
    """Rhoden per-arm bound conc. n=0 -> arm ABSENT (0). n=1 -> monovalent 1:1.
    n>=2 -> bivalent saturating avidity using the CELL-SIDE co-engagement span."""
    R_local=np.maximum(np.asarray(R_local,float),0.0); C=np.maximum(np.asarray(C,float),0.0)
    if n_arm<=0: return np.zeros_like(R_local*1.0)
    if n_arm==1: return R_local*C/(KD+C+1e-12)
    span=coeng_span_nm if coeng_span_nm is not None else DEFAULT_ARM_REACH_NM
    ageff=geo_ageff_nM(rec_pc, r_cell_um, span)
    a=2.0*ageff*C/KD**2; b=1.0+2.0*C/KD
    x=np.where(a>1e-30,(-b+np.sqrt(b*b+4.0*a))/(2.0*a+1e-300),1.0/np.maximum(b,1e-30))
    x=np.clip(x,0.0,1.0)
    return R_local*(2.0*x*C/KD + ageff*x*x*C/KD**2)

def multiarm_bound(Cfree, R_CD3, R_TAA, R_cos, *,
                   n_CD3=1, n_TAA=1, n_cos=0,
                   KD_CD3=3.0, KD_TAA=0.3, KD_cos=1.0,
                   rec_CD3=1e4, rec_TAA=1e5, rec_cos=1e4,
                   r_cell_T=4.0, r_cell_tum=8.0,
                   span_bridge_nm=12.5, span_coeng_T_nm=None, span_coeng_tumor_nm=None):
    """Per-cell bound species (nM), vectorized over per-cell arrays. See module docstring for the
    three spans and returned species. span_coeng_* None -> that co-engagement absent/compact-default."""
    Cfree=np.maximum(np.asarray(Cfree,float),0.0); z=np.zeros_like(Cfree*np.asarray(R_CD3,float)*1.0)
    # ---- single-arm bound pools on home cells (bivalent uses the cell-side co-engagement span) ----
    Cb3 = bound_arm(n_CD3, R_CD3, Cfree, rec_CD3, r_cell_T,   KD_CD3, span_coeng_T_nm)      # CD3 (T)
    CbC = bound_arm(n_cos, R_cos, Cfree, rec_cos, r_cell_T,   KD_cos, span_coeng_T_nm)      # costim (T)
    CbT = bound_arm(n_TAA, R_TAA, Cfree, rec_TAA, r_cell_tum, KD_TAA, span_coeng_tumor_nm)  # TAA (tumor)
    # ---- cross-cell bridges (feasibility from span_bridge vs cleft) ----
    f_bridge=_bridge_feasibility(span_bridge_nm)
    occ_b = Cfree/(KD_TAA/max(f_bridge,1e-3) + Cfree + 1e-12)
    Cb_kill    = np.minimum(Cb3, CbT)*occ_b if (n_CD3>0 and n_TAA>0) else z
    Cb_costTAA = np.minimum(CbC, CbT)*occ_b if (n_cos>0 and n_TAA>0) else z
    # ---- CIS: costim co-engaged with CD3 on the SAME T cell (T-side span) ----
    if n_CD3>0 and n_cos>0:
        p_cis = _cis_feasibility(span_coeng_T_nm, gap_match_nm=12.5, tol_nm=8.0)  # span-match -> co-engage
        Cb_cis = np.minimum(Cb3, CbC)*p_cis
    else:
        Cb_cis=z; p_cis=0.0
    return dict(Cb_kill=Cb_kill, Cb_costim=CbC, Cb_cis=Cb_cis, Cb_costTAA=Cb_costTAA,
                Cb3=Cb3, CbC=CbC, CbT=CbT, p_cis=float(p_cis))

# --------------------------------------------------------------------------------------------------
# FORMAT LIBRARY — commonly reported TCE construct formats -> (arm counts + characteristic geometry).
# Spans are literature/architecture geometry estimates (Fab-Fab ~13nm; tandem-scFv flexible ~13nm;
# C-term fusion end ~5-6nm; Fc top-to-bottom ~10-15nm); AF3-derived spans OVERRIDE per construct.
# span_coeng_* = None means that cell-side has no co-engagement (single arm there).
# --------------------------------------------------------------------------------------------------
FORMATS = {
    # 1+1 CD3xTAA
    "BiTE":            dict(n_CD3=1, n_TAA=1, n_cos=0, span_bridge_nm=13.0,
                            span_coeng_T_nm=None, span_coeng_tumor_nm=None, mw_kda=55),
    "IgG_1x1":         dict(n_CD3=1, n_TAA=1, n_cos=0, span_bridge_nm=13.0,
                            span_coeng_T_nm=None, span_coeng_tumor_nm=None, mw_kda=150),
    "DART_Fc":         dict(n_CD3=1, n_TAA=1, n_cos=0, span_bridge_nm=6.5,
                            span_coeng_T_nm=None, span_coeng_tumor_nm=None, mw_kda=105),
    # 2+1 (bivalent TAA + mono CD3) — glofitamab-like tumor-side avidity
    "IgG_2TAA_1CD3":   dict(n_CD3=1, n_TAA=2, n_cos=0, span_bridge_nm=13.0,
                            span_coeng_T_nm=None, span_coeng_tumor_nm=13.0, mw_kda=195),
    # trispecific CD3 x TAA x costim — tetravalent homodimer w/ C-term fusion
    #   CIS (height-matched costim epitope): T-side co-engagement span ~ compact -> cis co-fires
    "tetravalent_Cterm_cis":  dict(n_CD3=1, n_TAA=2, n_cos=1, span_bridge_nm=13.0,
                            span_coeng_T_nm=12.5, span_coeng_tumor_nm=13.0, mw_kda=200),
    #   TRANS (tall CRD1 costim epitope, ~60A height mismatch): T-side span large -> cis fails
    "tetravalent_Cterm_trans":dict(n_CD3=1, n_TAA=2, n_cos=1, span_bridge_nm=13.0,
                            span_coeng_T_nm=60.0, span_coeng_tumor_nm=13.0, mw_kda=200),
    # trispecific with bivalent TAA + bivalent CD3/costim variants can be added per construct CSV.
}

def bound_for_format(Cfree, R_CD3, R_TAA, R_cos, fmt, **kd_rec_over):
    """Convenience: look up a named format's arm counts + spans, then call multiarm_bound.
    kd_rec_over lets the caller supply per-molecule KD_*/rec_* (target-specific), overriding defaults."""
    if fmt not in FORMATS: raise KeyError(f"unknown format {fmt}; have {list(FORMATS)}")
    g={k:v for k,v in FORMATS[fmt].items() if k!='mw_kda'}
    g.update(kd_rec_over)
    return multiarm_bound(Cfree, R_CD3, R_TAA, R_cos, **g)
