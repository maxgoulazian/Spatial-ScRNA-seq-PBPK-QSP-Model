"""Full whole-body TCE PD re-validation — drives the signaling-wired per-cell PD engine per engager over
the FULL gridded organ set, reads systemic IL-6 (physically calibrated) + target-cell depletion, compares
to the clinical PD curves in mab_tce_pkpd.sqlite. Build-once per engager (arm-independent transport)."""
import numpy as np, json, sys, time, sqlite3
KWS=sys.argv[1]; sys.path.insert(0,f"{KWS}/handoff"); sys.path.insert(0,KWS)
import qsp_costim_window_v2 as q
from coupled_percell_pd import CoupledPerCellPD
pb=q._PBPKArrays()
# per-antigen internalization rate kint (1/day) from the curated antigen kinetics table (CD20 0.02, BCMA 2.0,
# GPRC5D 0.2, ...). The kint drives the per-cell Rhoden TMDD sink; it is a TARGET property, never a constant.
_AKT=json.load(open(f"{KWS}/handoff/antigen_kinetics_table.json"))["membrane"]
# model targets are HGNC gene symbols; the kinetics table keys some antigens by common name.
_KINT_ALIAS={"MS4A1":"CD20","FOLH1":"PSMA","IL3RA":"CD123","TNFRSF17":"TNFRSF17","ERBB2":"HER2"}
def antigen_kint(tgt, default=0.15):
    v=_AKT.get(tgt) or _AKT.get(_KINT_ALIAS.get(tgt,""))
    if isinstance(v,dict): return float(v.get("kint_perday", v.get("kint", default)))
    if isinstance(v,(list,tuple)) and v: return float(v[0])   # [central, alt, rationale]
    return float(default)
def antigen_kdeg(tgt, default=0.5):
    """free-receptor turnover kdeg (/day) = element [1] of the membrane kinetics table; ksyn=kdeg*Rtot set-point."""
    v=_AKT.get(tgt) or _AKT.get(_KINT_ALIAS.get(tgt,""))
    if isinstance(v,dict): return float(v.get("kdeg_perday", v.get("kdeg", default)))
    if isinstance(v,(list,tuple)) and len(v)>1: return float(v[1])
    return float(default)
def kinetic_binding_params(cfg, tgt):
    """Per-molecule (kon1 /M/s, koff1 /s, kdeg /day) for the TMDD receptor arm. kon/koff from the molecule
    spec when the param agents supplied them; else standard mAb kon=1e5 /M/s and koff=kon*KD (KD in M)."""
    # ONE MOLECULE, ONE AFFINITY: prefer the MEASURED KD (KD_norm, merged from eng_params_normalized.json)
    # over the ENG dict's value, exactly as attach_pd does. Using cfg['KD'] here while the PD path used
    # KD_norm gave the SAME antigen two different affinities in the same run (talquetamab 2.0 vs 11.0 nM).
    KD_nM=float(cfg.get('KD_norm', cfg['KD'])); kdeg=cfg.get('kdeg_perday', antigen_kdeg(tgt))
    kon=cfg.get('kon_TAA_perM_pers', cfg.get('kon_perM_pers'))
    koff=cfg.get('koff_TAA_pers', cfg.get('koff_pers'))
    # NOTE ON PROVENANCE: kon=1e5 is a GENERIC PLACEHOLDER, not a measurement. The SQLite audit (2026-07-13)
    # shows 15 of 22 molecules carry kon=1e5 with koff BACK-DERIVED as kon*KD. For those molecules the
    # KINETICS ARE ASSUMED, and any conclusion that depends on kon/koff (rather than on KD alone) is
    # ASSUMPTION-DRIVEN for them. Measured kinetics exist only for: teclistamab, elranatamab, catumaxomab,
    # and the rituximab-lineage CD20 kon (4.3e5). Do not present the rest as kinetically validated.
    if kon is None: kon=1e5                              # GENERIC PLACEHOLDER — see note above
    if koff is None: koff=kon*(KD_nM*1e-9)               # DERIVED, not measured: koff = kon*KD (nM->M)
    return float(kon), float(koff), float(kdeg)
CAL=json.load(open(f"{KWS}/handoff/wholebody_cyto_calib.json")); IL6_SCALE=CAL['il6_eng_scale_pgml_per_raw']
import os
import pd_model_config as PDC          # SINGLE SOURCE OF TRUTH for engine + params
PD_KINETICS = PDC.is_kinetic()         # canonical engine choice (env override for A/B: PD_ENGINE_OVERRIDE)
_ov=os.environ.get("PD_ENGINE_OVERRIDE","")
if _ov in ("kinetic","qss"): PD_KINETICS=(_ov=="kinetic")
if PD_KINETICS:
    IL6_SCALE=PDC.IL6_SCALE_KIN; K_DEATH_KIN=PDC.K_DEATH; KIN_PARAMS=PDC.kin_params()
else:
    K_DEATH_KIN=1.0; KIN_PARAMS=None   # QSS uses the QSS IL6_SCALE already loaded above
REG=json.load(open(f"{KWS}/handoff/regimen_schedules_final.json"))
RT=json.load(open(f"{KWS}/handoff/Rtot_wholebody_final.json"))['Rtot_nM']
bl=json.load(open(f"{KWS}/handoff/bec_lec_masks.json"))
# full gridded organ set (target-bearing organs matter most; include all for systemic cytokine sum)
ORGANS=['spleen','bone','large_int','liver','lung','small_int','pancreas','kidney','skin','heart','adipose']
def arr(a): return {o:a[i] for i,o in enumerate(pb.names)}  # ALL pbpk tissues incl tumor (idx 14)
#   built over full pb.names so run_organs=ORGANS+['tumor'] (solid runs) resolves Vis/Vv/Q/L/sigV['tumor']
Q,L,sigV,Vis,Vv=arr(pb.Q),arr(pb.L),arr(pb.sigV),arr(pb.Vis),arr(pb.Vv)
# engager -> target gene col, KD_TAA, mw, fFcRn, n_arm(TAA valency), IL6-observed peak (clinical anchor)
ENG={
 # route per clinical label: mosun/glofit IV; epcoritamab/BCMA/GPRC5D SC (F_sc bioavail, ka slow absorption)
 # il6_obs = clinical IL-6 anchor (pg/mL). DISPLAY/SCORING ONLY -- it never enters the mechanism (the model's
 # IL-6 is emergent from per-cell myeloid emitters), so it cannot bias a simulated value.
 #
 # ANCHOR PROVENANCE AUDIT, 2026-07-13. THE PREVIOUS VALUES HERE WERE CONTAMINATED.
 # EVERY anchor must be the SAME STATISTIC. The old table silently mixed a population MEAN (teclistamab)
 # with individual PEAKS and with numbers that had no source at all. That manufactured a fake ~26x clinical
 # dynamic range which the mechanism was then repeatedly "fixed" to chase. What the old values actually were:
 #   mosunetuzumab 570  -> NO SOURCE EXISTS. Deleted. The real population MEAN is 152.
 #   elranatamab   340  -> cited to "MagnetisMM-3 Fig 6". That paper contains ZERO mentions of IL-6 and has
 #                         only four figures. The citation is to a figure that does not exist.
 #   elranatamab   191  -> A PAGE NUMBER. A dot-leader page reference from the Table of Figures of FDA BLA
 #                         761345. It is not a concentration. It was being used as one.
 #   elranatamab   230 / 366.88 -> NO SOURCE EXISTS.
 #   teclistamab   288  -> REAL, but it is the HIGHEST INDIVIDUAL PATIENT Cmax -- an ORDER STATISTIC that
 #                         scales with cohort size N, so it is NOT comparable across trials and NOT a
 #                         central-tendency anchor. Not usable here.
 # ELRANATAMAB HAS NO CLINICAL IL-6 VALUE IN EXISTENCE (no paper, no FDA review, no label) -> il6_obs=None.
 #
 # SURVIVING ANCHORS = POPULATION MEAN IL-6 peak Cmax. One statistic. Verbatim-sourced. Nothing else.
 #   teclistamab 21  -- PMID 38831634: "The mean IL-6 peak concentration (Cmax) was 21 pg/mL"
 #   mosunetuzumab 152
 # Individual severe-CRS outliers (elran 23600 pleural/case-report, tecli 7180 pleural) EXCLUDED.
 'mosunetuzumab':dict(tgt='MS4A1',KD=5.0,mw=146.0,fFcRn=0.89,narm=1,il6_obs=152.0,route='IV'),   # population MEAN
 'glofitamab':   dict(tgt='MS4A1',KD=5.0,mw=195.0,fFcRn=0.89,narm=2,il6_obs=None,route='IV'),   # 2:1 bivalent CD20
 'epcoritamab':  dict(tgt='MS4A1',KD=5.0,mw=148.0,fFcRn=0.89,narm=1,il6_obs=None,route='SC',F_sc=0.6,ka=0.25),
 'elranatamab':  dict(tgt='TNFRSF17',KD=0.15,mw=148.5,fFcRn=0.70,narm=1,il6_obs=None,route='SC',F_sc=0.6,ka=0.25),  # mw 148.5 (Elrexfio label), IgG2 fFcRn 0.70. il6_obs=None: NO clinical IL-6 value exists for this molecule. measured kon/koff/KD_TAA=0.04 from normalized merge
 'teclistamab':  dict(tgt='TNFRSF17',KD=0.15,mw=146.0,fFcRn=0.89,narm=1,il6_obs=21.0,route='SC',F_sc=0.6,ka=0.25),   # population MEAN (PMID 38831634)
 'talquetamab':  dict(tgt='GPRC5D',KD=2.0,mw=146.0,fFcRn=0.89,narm=1,il6_obs=None,route='SC',F_sc=0.6,ka=0.25),      # 18.2 was a MEDIAN, not a mean -> not comparable to the mean pair. Excluded pending a mean.
 # ROUTE-MATCHED PK VALIDATION VARIANTS (2026-07-13):
 # PK route-matched variants. il6_obs corrected to the SAME sourced population MEANS (was 570 / 21.5 — the
 # fabricated and mis-rounded values). These variants are for PK validation; the IL-6 anchor is display-only.
 'mosunetuzumab_sc':dict(tgt='MS4A1',KD=5.0,mw=146.0,fFcRn=0.89,narm=1,il6_obs=152.0,route='SC',F_sc=0.6,ka=0.25),
 'teclistamab_iv': dict(tgt='TNFRSF17',KD=0.15,mw=146.0,fFcRn=0.89,narm=1,il6_obs=21.0,route='IV'),
 'teclistamab_iv_low':dict(tgt='TNFRSF17',KD=0.15,mw=146.0,fFcRn=0.89,narm=1,il6_obs=21.0,route='IV'),
 # ---- AUTO-GENERATED from eng_params_normalized.json (SQLite param audit, 2026-07-13) ----
 # Every field READ from the normalized file. Nothing invented. il6_obs=None: no sourced clinical
 # IL-6 anchor exists for these, and a fabricated one is worse than none.
 'odronextamab':dict(tgt='MS4A1',KD=5.0,mw=146.0,fFcRn=1.0,narm=1,il6_obs=None,route='IV'),   # AUTO from eng_params_normalized (audit 2026-07-13)
 'blinatumomab':dict(tgt='CD19',KD=1.49,mw=54.0,fFcRn=0.0,narm=1,il6_obs=None,route='IV'),   # AUTO from eng_params_normalized (audit 2026-07-13)
 'linvoseltamab':dict(tgt='TNFRSF17',KD=1.0,mw=150.0,fFcRn=0.7,narm=1,il6_obs=None,route='IV'),   # AUTO from eng_params_normalized (audit 2026-07-13)
 'alnuctamab':dict(tgt='TNFRSF17',KD=0.3,mw=150.0,fFcRn=0.7,narm=2,il6_obs=None,route='SC',F_sc=0.6,ka=0.25),   # AUTO from eng_params_normalized (audit 2026-07-13)
 'cevostamab':dict(tgt='FCRL5',KD=0.016,mw=145.2,fFcRn=0.7,narm=1,il6_obs=None,route='IV'),   # AUTO from eng_params_normalized (audit 2026-07-13)
 'tarlatamab':dict(tgt='DLL3',KD=0.64,mw=100.0,fFcRn=0.6,narm=1,il6_obs=None,route='IV'),   # AUTO from eng_params_normalized (audit 2026-07-13)
 'cibisatamab':dict(tgt='CEACAM5',KD=16.0,mw=194.0,fFcRn=0.75,narm=2,il6_obs=None,route='IV'),   # AUTO from eng_params_normalized (audit 2026-07-13)
 'tebentafusp':dict(tgt='PMEL',KD=0.015,mw=77.0,fFcRn=0.0,narm=1,il6_obs=None,route='IV'),   # AUTO from eng_params_normalized (audit 2026-07-13)
 'catumaxomab':dict(tgt='EPCAM',KD=0.56,mw=150.0,fFcRn=0.3,narm=1,il6_obs=None,route='IP'),   # AUTO from eng_params_normalized (audit 2026-07-13)
 'solitomab':dict(tgt='EPCAM',KD=5.0,mw=55.0,fFcRn=0.0,narm=1,il6_obs=None,route='IV'),   # AUTO from eng_params_normalized (audit 2026-07-13)
 'pasotuxizumab':dict(tgt='FOLH1',KD=47.0,mw=55.0,fFcRn=0.0,narm=1,il6_obs=None,route='SC',F_sc=0.6,ka=0.25),   # AUTO from eng_params_normalized (audit 2026-07-13)
 'acapatamab':dict(tgt='FOLH1',KD=1.0,mw=100.0,fFcRn=0.6,narm=1,il6_obs=None,route='IV'),   # AUTO from eng_params_normalized (audit 2026-07-13)
 'runimotamab':dict(tgt='ERBB2',KD=0.3,mw=150.0,fFcRn=0.75,narm=1,il6_obs=None,route='IV'),   # AUTO from eng_params_normalized (audit 2026-07-13)
 'cinrebafusp_alfa':dict(tgt='ERBB2',KD=0.3,mw=170.0,fFcRn=0.75,narm=2,il6_obs=None,route='IV'),   # AUTO from eng_params_normalized (audit 2026-07-13)
}
# ---- OPTIONAL CONSTRUCT SCREEN --------------------------------------------------------------------------
# ENG holds CLINICAL molecules. A designed CONSTRUCT (e.g. scr_bi_MS4A1_TNFRSF9) is not in it, so
# `cfg = ENG[name]` KeyErrors and the whole screen is unrunnable. costim_screen.py already emits construct
# rows in EXACTLY the ENG schema (tgt/KD/mw/fFcRn/narm/route + costim_arm/KD_costim_nM/KD_CD3_nM), so the
# only thing missing was a way to hand them to the runner. Merge them in from a JSON when asked.
# SCREEN_JSON unset -> ENG is untouched -> clinical runs are byte-identical.
_scr = os.environ.get("SCREEN_JSON")
if _scr and os.path.exists(_scr):
    _sj = json.load(open(_scr))
    _rows = _sj.get("screen", _sj)
    ENG.update(_rows)
    print(f"[screen] merged {len(_rows)} construct definitions from {_scr}", flush=True)
    # DOSING FOR SCREEN CONSTRUCTS.  A counterscreen must hold the regimen CONSTANT across every
    # construct -- otherwise the ranking is confounded by the schedule and we would be ranking dosing,
    # not biology.  Each construct therefore inherits ONE common reference regimen (default: the
    # published teclistamab step-up, from regimen_schedules_final.json), unless the row carries its own
    # 'sched'.  This is a stated DESIGN CHOICE, not a fitted parameter; it is recorded on every result.
    _ref = os.environ.get("SCREEN_REGIMEN", "teclistamab")
    if _ref not in REG:
        raise RuntimeError(f"[screen] SCREEN_REGIMEN={_ref!r} is not in regimen_schedules_final.json "
                           f"(have: {sorted(REG)[:8]}...). Refusing to invent a dose schedule.")
    _n_inherit = 0
    for _k, _v in _rows.items():
        if isinstance(_v, dict) and _v.get("sched"):
            REG[_k] = _v["sched"]
        else:
            REG[_k] = REG[_ref]
            _n_inherit += 1
    print(f"[screen] {_n_inherit}/{len(_rows)} constructs inherit the common reference regimen "
          f"'{_ref}' ({len(REG[_ref])} doses; held CONSTANT across the screen so the arm/target is the "
          f"only variable)", flush=True)

TSIM={'mosunetuzumab':49.0,'glofitamab':24.0,'epcoritamab':24.0,'elranatamab':24.0,'teclistamab':200.0,'talquetamab':24.0,'mosunetuzumab_sc':160.0,'teclistamab_iv':10.0,'teclistamab_iv_low':10.0,'odronextamab':24.0,'blinatumomab':24.0,'linvoseltamab':24.0,'alnuctamab':24.0,'cevostamab':24.0,'tarlatamab':24.0,'cibisatamab':24.0,'tebentafusp':24.0,'catumaxomab':24.0,'solitomab':24.0,'pasotuxizumab':24.0,'acapatamab':24.0,'runimotamab':24.0,'cinrebafusp_alfa':24.0}
def run(name, tsim=None):
    cfg=ENG[name]; tgt=cfg['tgt']; tsim=tsim or TSIM.get(name,24.0)
    # TSIM_DAYS env override. The long TSIMs (teclistamab 200 d) exist to reach PK steady state; the CRS/IL-6
    # PEAK is an EARLY event (priming/step-up doses, first ~1-2 weeks). So an IL-6 validation run does NOT need
    # the full PK window -- TSIM_DAYS=24 captures the peak at ~1/8 the cost. PK Cmax/AUC from such a run are
    # TRUNCATED and must NOT be quoted as steady-state values.
    _ts=os.environ.get("TSIM_DAYS")
    if _ts:
        tsim=float(_ts)
        print(f"[TSIM override] {name}: tsim={tsim} d (IL-6/CRS peak window; PK is TRUNCATED - do not quote steady-state PK)", flush=True)
    # ===== MERGE MEASURED PER-MOLECULE KINETICS from eng_params_normalized.json =====
    # The ENG dict above carries route/variant specifics (SC vs IV, il6_obs); the normalized file carries
    # the MEASURED kon/koff (both arms), KD, kint, kdeg, mw, fFcRn the param agents gathered. Merge so every
    # molecule runs on its real kinetics, not the generic kon=1e5 fallback. cfg keys win only where the
    # normalized file lacks them (route/variant). Variant names (teclistamab_iv) map to base (teclistamab).
    try:
        _NORM=json.load(open(f"{KWS}/handoff/eng_params_normalized.json"))
    except Exception:
        _NORM={}
    _base=name
    for _suf in ("_iv_low","_iv","_sc"):
        if _base.endswith(_suf): _base=_base[:-len(_suf)]; break
    _np=_NORM.get(_base) or _NORM.get(name)
    if _np:
        # map normalized schema -> cfg fields consumed by kinetic_binding_params + build
        _map={'kon_TAA_perM_pers':'kon_TAA_perM_pers','koff_TAA_pers':'koff_TAA_pers',
              'kon_CD3_perM_pers':'kon_CD3_perM_pers','koff_CD3_pers':'koff_CD3_pers',
              'kdeg_perday':'kdeg_perday','kint_perday':'kint_perday'}
        for _k,_dst in _map.items():
            if _np.get(_k) not in (None,0): cfg.setdefault(_dst,_np[_k])
        # KD: prefer measured TAA KD if the ENG entry didn't pin one deliberately
        if _np.get('KD_TAA_nM') not in (None,0): cfg.setdefault('KD_norm',_np['KD_TAA_nM'])
        cfg.setdefault('KD_CD3_nM', _np.get('KD_CD3_nM',40.0))
        cfg.setdefault('_measured_vs_derived', _np.get('measured_vs_derived',''))
    HEME_TARGETS={'MS4A1','TNFRSF17','GPRC5D','IL3RA','CD33','FLT3','CD19','FCRL5'}
    cancer_type=cfg.get('cancer_type') or os.environ.get("CANCER_TYPE") or ("heme" if tgt in HEME_TARGETS else "solid")
    # SOLID tumor is a SPATIAL ABM organ (ECM + BEC/LEC + per-cell Rhoden): include it in the organ list so
    # attach_pd builds its per-cell kill and it participates in transport. HEME/BLOOD are plasma-driven (separate).
    run_organs = (ORGANS + ['tumor']) if cancer_type=="solid" else list(ORGANS)
    pools={o:RT.get(tgt,{}).get(o,0.0) for o in run_organs}
    _kon1,_koff1,_kdeg = kinetic_binding_params(cfg, tgt)
    # FIX-1-UPGRADE: thread the MEASURED per-molecule kon/koff (BOTH arms) into the PD synapse so PK and PD
    # bind on ONE identical scheme (was: PD synapse rebuilt from generic kon=1e5, koff=kon*KD_registry).
    # eng_params stores kon in /M/s and koff in /s -- exactly the units OrganPD expects (kon_*_perM_s, koff_*_pers).
    # Fallback to the generic KIN_PARAMS defaults where a molecule has no measured value.
    _kpm = dict(KIN_PARAMS) if KIN_PARAMS else None
    if _kpm is not None:
        if cfg.get('kon_TAA_perM_pers'): _kpm['kon_TAA_perM_s']=cfg['kon_TAA_perM_pers']
        if cfg.get('koff_TAA_pers'):     _kpm['koff_TAA_pers']=cfg['koff_TAA_pers']
        if cfg.get('kon_CD3_perM_pers'): _kpm['kon_CD3_perM_s']=cfg['kon_CD3_perM_pers']
        if cfg.get('koff_CD3_pers'):     _kpm['koff_CD3_pers']=cfg['koff_CD3_pers']
 # ONE MOLECULE, ONE AFFINITY. (bug found by the SQLite param audit, 2026-07-13.)
    # This line used cfg['KD'] (the ENG dict's value) while attach_pd/attach_heme_pd/attach_blood_pd below
    # all used cfg.get('KD_norm', cfg['KD']) -- the MEASURED value merged from eng_params_normalized.json.
    # So the PK antigen sink and the PD synapse were binding the SAME antigen with DIFFERENT affinities:
    #     talquetamab  PK 2.0 nM  vs  PD 11.0 nM   (5.5x apart)
    #     elranatamab  PK 0.15 nM vs  PD 0.04 nM   (3.75x apart)
    # These are physically the same binding event. The measured KD wins everywhere.
    _KD_used = cfg.get('KD_norm', cfg['KD'])
    if abs(_KD_used - cfg['KD']) > 1e-12:
        print(f"[KD] {name}: using MEASURED KD_TAA={_KD_used:g} nM everywhere "
              f"(ENG carried {cfg['KD']:g} nM; PK and PD are now consistent)", flush=True)
    m=CoupledPerCellPD(run_organs,tgt,tgt,_KD_used,cfg['narm'],antigen_kint(tgt),cfg['mw'],cfg['fFcRn'],
                       f"{KWS}/handoff/agents",pools,Q,L,sigV,Vis,Vv,bl,dt=0.02,
                       kon1_perM_pers=_kon1, koff1_pers=_koff1, kdeg_perday=_kdeg,
                       span_coeng_nm=cfg.get('span_coeng_nm', cfg.get('span_bridge_nm', 12.5)))
    # COSTIM ARM: attach_pd has always supported it (costim_arm -> per-cell '<ARM>_copies' + the per-cell
    # signaling engine), but the runner never PASSED it -- so every costim construct silently ran as a plain
    # TCE. Thread it from cfg so the costim screen is actually possible. cfg['costim_arm']=None (the clinical
    # validation molecules) reproduces the previous plain-TCE behaviour EXACTLY -- this is additive.
    # 11 costim arms have per-cell density in all 11 organs (verified 2026-07-13): TNFRSF9 (4-1BB), TNFRSF4
    # (OX40), TNFRSF18 (GITR), CD28, ICOS, CD27, CD40, TNFRSF25, CD2, CD226 (DNAM-1), TNFRSF14 (HVEM).
    m.attach_pd(tgt, KD_CD3_nM=cfg.get('KD_CD3_nM',40.0), KD_TAA_nM=cfg.get('KD_norm', cfg['KD']), cd3_col='CD3E',
                costim_arm=cfg.get('costim_arm'), KD_costim_nM=cfg.get('KD_costim_nM', 1.0),
                pd_kinetics=PD_KINETICS, kin_params=(_kpm if PD_KINETICS else None))
    # ===== CANCER-TYPE SWITCH (cancer_type computed above): heme subtype XOR solid, never both. =====
    # The 11 healthy organs + normal blood are ALWAYS present (body: PK sinks, CRS, tox). SOLID tumor is
    # already in the organ list (built with attach_pd); HEME malignancy is a separate plasma-driven compartment.
    m.cancer_type=cancer_type
    is_heme = (cancer_type in ("AML","DLBCL","myeloma","heme")) or (cancer_type not in ("solid",) and tgt in HEME_TARGETS)
    if is_heme:
        # heme malignancy compartment (routing map picks AML/DLBCL/myeloma npz by target antigen)
        m.attach_heme_pd(tgt, KD_CD3_nM=cfg.get('KD_CD3_nM',40.0), KD_TAA_nM=cfg.get('KD_norm', cfg['KD']), cd3_col='CD3E',
                         pd_kinetics=PD_KINETICS, kin_params=(_kpm if PD_KINETICS else None))
    # solid: tumor already included as a spatial organ above (ECM-throttled per-cell ABM kill via attach_pd)
    # NORMAL BLOOD compartment — ALWAYS attached (TMDD sink + normal-cell depletion tox + circulating-T CRS)
    if os.environ.get("USE_BLOOD_COMPARTMENT","1")=="1":
        m.attach_blood_pd(tgt, KD_CD3_nM=cfg.get('KD_CD3_nM',40.0), KD_TAA_nM=cfg.get('KD_norm', cfg['KD']), cd3_col='CD3E',
                          pd_kinetics=PD_KINETICS, kin_params=(_kpm if PD_KINETICS else None))
    sched=[(t,d) for t,d in REG[name]]
    rt=cfg.get('route','IV')
    _snap=[float(x) for x in os.environ.get("SNAP_DAYS","1,7").split(",") if x.strip()]
    _prog=int(os.environ.get('PROGRESS_EVERY','0'))
    kw=dict(route=rt, k_death=(K_DEATH_KIN if PD_KINETICS else 1.0), rec=200, pd_every=3, snap_times=_snap,
            progress_every=_prog, progress_tag=name)
    if rt=='SC': kw.update(F_sc=cfg.get('F_sc',0.6), ka_sc=cfg.get('ka',0.25))
    else: kw.update(iv_inf_h=cfg.get('iv_inf_h',2.0))   # real IV infusion duration (clinical TCE ~2h), not bolus
    r=m.simulate_pd(sched, tsim=tsim, **kw)
    # --- SPATIAL SNAPSHOTS: multi-timepoint per-cell fields (render-ready: x,y,C,bound_nM,R,labs,surv) ---
    # m.snaps[day][organ] built by simulate_pd(snap_times=...). Save as pickle (nested dict, ragged arrays).
    import pickle as _pkl
    _tag2=os.environ.get("PD_OUT_TAG","")
    _snap_out={}
    for _day,_od in getattr(m,"snaps",{}).items():
        _snap_out[_day]={}
        for _o,_d in _od.items():
            _snap_out[_day][_o]={k:(v if v is None else np.asarray(v)) for k,v in _d.items()}
    # heme compartment final survival field (circulating blasts, no transport graph)
    if hasattr(m,"heme_pd"):
        _h=m.heme_pd
        _snap_out["heme_final"]={"heme":dict(x=np.asarray(_h.x),y=np.asarray(_h.y),
            surv=np.exp(-_h.kill_hazard), is_target=np.asarray(_h.is_target),
            labs=np.asarray(_h.labs) if hasattr(_h,"labs") else np.array([]))}
    if hasattr(m,"blood_pd"):
        _b=m.blood_pd
        _snap_out["blood_final"]={"blood":dict(x=np.asarray(_b.x),y=np.asarray(_b.y),
            surv=np.exp(-_b.kill_hazard), is_target=np.asarray(_b.is_target),
            labs=np.asarray(_b.labs) if hasattr(_b,"labs") else np.array([]))}
    _pkl.dump(_snap_out, open(f"{KWS}/handoff/tce_spatial_{_tag2}{name}.pkl","wb"))
    # ===== IL-6: MECHANISTIC ONLY. Per-cell myeloid emitters (PMID 29808005/29808007) -> plasma ODE against
    # measured clearance. Saturation and per-molecule differences EMERGE from the finite, spatially-distributed
    # myeloid pool. No Emax, no EC50, no fitted scale.
    #
    # THE FITTED LEGACY PATH IS DELETED (2026-07-13, Max: "just delete the legacy value from the code").
    # It computed  il6 = sys_cyto_rate['IL6'] * IL6_SCALE, where IL6_SCALE was ONE constant fitted so that
    # mosunetuzumab -> 570 pg/mL. That 570 has since been shown to have NO SOURCE -- it is a fabricated
    # anchor. Worse, the old code did:
    #       if il6.size==0: il6 = il6_legacy      # "safety: never crash a run"
    # i.e. on ANY failure of the mechanistic recorder it SILENTLY substituted the fitted number and emitted
    # it under the same field name (il6_peak). A fit to a fabricated anchor, wearing a mechanistic label,
    # with no warning. That is exactly the failure mode that let a PAGE NUMBER (elranatamab "191", a
    # dot-leader from an FDA Table of Figures) circulate as a clinical IL-6 value.
    #
    # It is now a HARD ERROR. An empty mechanistic array means the model did not compute IL-6 -- that is a
    # bug to fix, not a number to substitute. Crashing loudly beats reporting a fabricated constant quietly.
    il6=np.array(r.get('il6_plasma_pgml') or [])
    if il6.size==0:
        raise RuntimeError(
            f"[{name}] MECHANISTIC IL-6 was not recorded (il6_plasma_pgml is empty). Refusing to fall back to "
            f"the retired fitted IL6_SCALE path (it was fitted to mosunetuzumab=570, a FABRICATED anchor). "
            f"Fix the IL-6 recorder in coupled_percell_pd.py; do not substitute a constant.")
    # SAVE EVERYTHING — PK curve, per-organ kill trajectories, ALL cytokine time-courses, heme depletion.
    # Never discard a computed time-course: re-runs are expensive and the curves are the primary output.
    cyto_sys={c:[float(v) for v in r['sys_cyto_rate'][c]] for c in r.get('sys_cyto_rate',{})}   # systemic cytokine RATE (clinical-peak comparable)
    cyto_sys_cum={c:[float(v) for v in r['sys_cyto'][c]] for c in r.get('sys_cyto',{})}          # systemic cytokine (cumulative)
    return dict(name=name, tgt=tgt,
                # --- time bases ---
                t=r['pd_t'].tolist(), t_pk=r['t'].tolist(),
                # --- PK (plasma) ---
                Cplasma_ugml=r['Cplasma_ugml'].tolist(), Cmax=float(r['Cplasma_ugml'].max()),
                # --- PD: cytokines (all species, systemic; IL-6 in pg/mL via scale) ---
                il6_pgml=il6.tolist(), il6_peak=float(il6.max()),
                # il6_legacy_pgml / il6_legacy_peak / IL6_SCALE: DELETED 2026-07-13. The fitted path is gone;
                # emitting it alongside the mechanistic value invited exactly the mix-up it caused.
                il6_prod_pg_hr=r.get('il6_prod_pg_hr'),      # myeloid production (pg/hr) that drove the plasma ODE
                il6_method="mechanistic_myeloid_percell",    # provenance: no fitted scale, no Emax, no EC50
                cyto_sys_rate=cyto_sys, cyto_sys_cum=cyto_sys_cum,
                # --- PD: PER-ORGAN cytokine time-courses (off-tumor / off-site tox signal, all species) ---
                cyto_organ={o:{c:[float(v) for v in r['cyto'][o][c]] for c in r['cyto'][o]} for o in m.organs},
                # --- PD: per-organ kill TRAJECTORIES (not just final) ---
                kill_traj={o:[float(v) for v in r['kill_frac'][o]] for o in m.organs},
                kill_final={o:float(r['kill_frac'][o][-1]) for o in m.organs},
                target_organ_kill=float(max(r['kill_frac'][o][-1] for o in m.organs)),
                kill_ntgt={o:int(r.get('n_target',{}).get(o,0)) for o in m.organs},
                depletion_weighted=float(
                    sum(r['kill_frac'][o][-1]*r.get('n_target',{}).get(o,0) for o in m.organs)
                    / max(sum(r.get('n_target',{}).get(o,0) for o in m.organs),1)),
                # --- PD: heme malignant-blast depletion trajectory ---
                heme_depletion=(float(r['heme_kill'][-1]) if r.get('heme_kill') else None),
                heme_ntgt=int(r.get('heme_ntgt',0)),
                heme_depletion_t=([float(v) for v in r['heme_kill']] if r.get('heme_kill') else None),
                # --- NORMAL blood-cell depletion (on-target tox readout, e.g. circulating B-cell depletion) ---
                cancer_type=getattr(m,"cancer_type","solid"),
                blood_depletion=(float(r['blood_kill'][-1]) if r.get('blood_kill') else None),
                blood_ntgt=int(r.get('blood_ntgt',0)),
                blood_depletion_t=([float(v) for v in r['blood_kill']] if r.get('blood_kill') else None))
# single-engager mode: argv[2] = engager name -> writes handoff/tce_pd_<name>.json (parallel-safe, own file)
import os as _os
if _os.environ.get("CALIB_IMPORT_ONLY"):
    pass  # inert import for calibration harness
elif len(sys.argv)>2:
    name=sys.argv[2]; t0=time.time(); o=run(name); o['secs']=time.time()-t0
    obs=ENG[name]['il6_obs']; tag=f" (clinical IL6 {obs})" if obs else ""
    hd=o.get('heme_depletion'); hds=f" heme_blast_depl={hd:.3f}" if hd is not None else ""
    print(f"{name:15s} IL6peak={o['il6_peak']:7.1f} pg/mL{tag}  depl(wt)={o.get('depletion_weighted',float('nan')):.3f}{hds}  Cmax={o['Cmax']:.3f} ({o['secs']:.0f}s)",flush=True)
    _tag=os.environ.get("PD_OUT_TAG","")
    json.dump(o, open(f"{KWS}/handoff/tce_pd_{_tag}{name}.json","w"), indent=1)
    print(f"DONE {name}")
else:
    res={}
    for name in ENG:
        t0=time.time(); o=run(name); o['secs']=time.time()-t0; res[name]=o
        print(f"{name:15s} IL6peak={o['il6_peak']:7.1f} depl(wt)={o.get('depletion_weighted',float('nan')):.3f} ({o['secs']:.0f}s)",flush=True)
        json.dump(res, open(f"{KWS}/handoff/tce_pd_reval_result.json","w"), indent=1)
    print("DONE")
