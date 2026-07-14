"""signaling_dynamics.py — TIME-RESOLVED per-cell costim signaling for the whole-body QSP.
Couples three real data layers:
  1. per-cell costim RECEPTOR occupancy  (rna_to_receptor HPA/Glassman copies -> occupancy by arm at
     local drug conc)  -> the per-cell perturbation amplitude
  2. GRN signal PROPAGATION (signaling_engine, CellOracle-style) -> per-cell direction/shape of the
     program response, heterogeneous across cells
  3. program KINETICS from the hero time axis (Rest/8hr/48hr) -> real first-order rates k_on and the
     acute-vs-sustained shape; per-arm signed effect on each program

Per-cell signaling state S_p(t) for each program p follows first-order kinetics driven by engagement:
    dS_p/dt = k_on_p * ( drive_p(arm,cell) * occ_cell(t) - S_p )        [rise toward driven setpoint]
    (+ acute overshoot term for cytokines: fast rise k_on, then decay k_off to sustained fraction)
S_p(t) then modifies the PD ODE rates (prolif expansion, exhaustion decay, cytokine secretion, Treg
suppression) PER CELL, complementing (not replacing) the binary PROLIF/EXH veto gate.

drive_p(arm,cell) = arm's signed program effect (agonism = -KD_effect from kinetics) x GRN-propagated
per-cell shape (so two cells with the same arm but different network state respond differently)."""
import numpy as np, json, os
KWS="/media/balthasar-lab/RAID4/claude-science/orgs/3761bdf5-5d07-49ad-9607-2e0cc29f640d/workspaces/ad4cbbf0-ef0e-43dd-8e36-52612f369b40"
PROGRAMS=["proliferation","exhaustion","cyto_IL6","cyto_IFNG","cyto_TNF","cyto_IL2","suppression","effector"]
_K=None
def kinetics():
    global _K
    if _K is None:
        import os as _o
        _here=_o.path.dirname(_o.path.abspath(__file__))
        _c=[_o.path.join(_o.environ.get("SIGNAL_SURROGATE_DIR",""),"program_kinetics.json"),
            _o.path.join(_here,"signaling_surrogate","program_kinetics.json"),
            _o.path.join(_o.getcwd(),"rundir","handoff","signaling_surrogate","program_kinetics.json"),
            _o.path.join(_o.getcwd(),"handoff","signaling_surrogate","program_kinetics.json"),
            f"{KWS}/handoff/program_kinetics.json"]
        _f=next((x for x in _c if x and _o.path.exists(x)), None)
        if _f is None:
            raise RuntimeError("[signaling_dynamics] program_kinetics.json not found. Searched:\n  " + "\n  ".join(c for c in _c if c))
        _K=json.load(open(_f))
    return _K
# resolvable-rate ceiling: cytokines overshoot by 8hr, so cap k_on and add explicit decay to 48hr level
def _rates():
    K=kinetics(); R={}
    for p in PROGRAMS:
        d=K["programs"].get(p,{}); t=d.get("WT_traj",{})
        kon=d.get("k_on_perhr"); kon=kon if (kon and np.isfinite(kon)) else 0.1
        L0=t.get("Rest"); L8=t.get("Stim8hr"); L48=t.get("Stim48hr")
        # acute overshoot flag: 8hr level > 48hr level (cytokine burst then contraction)
        acute = (L8 is not None and L48 is not None and L0 is not None and L8>L48 and L8>L0)
        # k_off (decay from 8hr peak to 48hr sustained) if acute
        koff=np.nan
        if acute and (L8-L0)>1e-6:
            frac_remaining=max((L48-L0)/(L8-L0),1e-3)   # fraction of peak still present at 48hr
            koff=-np.log(frac_remaining)/40.0            # 8->48hr = 40hr window
        R[p]={"k_on_perhr":min(kon,1.5),"k_off_perhr":(koff if np.isfinite(koff) else 0.0),
              "acute":acute,"sustained_frac":(max((L48-L0)/(L8-L0),0.0) if (acute and (L8-L0)>1e-6) else 1.0)}
    return R
RATES=_rates()
def arm_drive(arm):
    """Per-arm signed program drive = agonism direction = -(KD effect at 48hr) from hero kinetics."""
    K=kinetics(); out={}
    for p in PROGRAMS:
        e=K["programs"].get(p,{}).get("arm_effect",{}).get(arm)
        out[p]= -e if (e is not None and np.isfinite(e)) else 0.0   # agonism = opposite of knockdown
    return out

# ---- lane drive-table adapter: MAGNITUDE from the signaling-surrogate lane's calibrated z-scores ----
# (per README: magnitude comes from deliverable-4 drive table; propagation/baseline give per-cell shape;
#  my hero kinetics give the timescale). Arm label alias: drive table uses receptor common names.
_ARM_ALIAS={"TNFRSF9":"4-1BB","TNFRSF4":"OX40","TNFRSF18":"GITR","TNFRSF8":"CD30","TNFRSF14":"HVEM",
            "TNFRSF25":"DR3","CD226":"DNAM1"}  # agent-column gene -> drive-table receptor label

# ---------------------------------------------------------------------------
# DRIVE-TABLE RESOLUTION  [added 2026-07-13]
# The two drive tables were loaded by ABSOLUTE Xeon path with an os.path.exists()
# guard.  Off the Xeon (RunPod, laptop, any container) the path does not exist, the
# guard returned None, and arm_drive_lane() then returned drive=0.0 FOR EVERY ARM --
# silently.  A 132-construct counterscreen would have completed, ranked nothing, and
# looked fine.  These tables ARE the ranking.  Missing table => hard stop, never zero.
# ---------------------------------------------------------------------------
def _surrogate_path(fname):
    here = os.path.dirname(os.path.abspath(__file__))
    cands = []
    env = os.environ.get("SIGNAL_SURROGATE_DIR")
    if env:
        cands.append(os.path.join(env, fname))
    cands += [
        os.path.join(here, "signaling_surrogate", fname),
        os.path.join(here, "..", "handoff", "signaling_surrogate", fname),
        os.path.join(os.getcwd(), "handoff", "signaling_surrogate", fname),
        os.path.join(os.getcwd(), "rundir", "handoff", "signaling_surrogate", fname),
        "/media/balthasar-lab/RAID4/costim_engager_counterscreen/handoff/signaling_surrogate/" + fname,
        "/Volumes/RAID4/costim_engager_counterscreen/handoff/signaling_surrogate/" + fname,
    ]
    for c in cands:
        if os.path.exists(c):
            return c
    raise RuntimeError(
        "[signaling_dynamics] REQUIRED drive table not found: %s\n"
        "  Searched:\n    %s\n"
        "  This table IS the per-arm costim ranking. Without it every arm scores drive=0.0\n"
        "  and the counterscreen silently ranks nothing. Refusing to run.\n"
        "  Set SIGNAL_SURROGATE_DIR to the directory holding the signaling_surrogate CSVs."
        % (fname, "\n    ".join(cands)))

_DRIVE_TBL=None
def _load_drive_table():
    global _DRIVE_TBL
    if _DRIVE_TBL is None:
        import pandas as pd
        f=_surrogate_path("per_arm_drive_magnitude_uncertainty.csv")
        _DRIVE_TBL=pd.read_csv(f).set_index("receptor")
    return _DRIVE_TBL
def arm_drive_lane(arm, kE=0.11, kC=0.20, kS=0.20):
    """Per-program signed DRIVE setpoint from the lane's calibrated z-scores (magnitude source of truth).
    effector_z on raw Schmidt CRISPRa scale -> x kE (locked). crs_z/supp_z on empirical-null +-3 -> x kC/kS.
    Cytokines split by CYTO_HIER-consistent share of crs_z. Returns {program: signed setpoint}."""
    T=_load_drive_table(); lbl=_ARM_ALIAS.get(arm,arm)
    if T is None:
        raise RuntimeError("[signaling_dynamics] drive table is None -- refusing to score arm %r with zero drive" % (arm,))
    if lbl not in T.index:
        raise RuntimeError("[signaling_dynamics] arm %r (label %r) is ABSENT from the drive table.\n"
                           "  Known arms: %s\n"
                           "  Scoring it as drive=0.0 would silently rank it as inert." % (arm, lbl, list(T.index)))
    r=T.loc[lbl]
    eff=kE*float(r["effector_z"]); crs=kC*float(r["crs_z"]); supp=kS*float(r["supp_z"])
    # cytokine split: IFNG/TNF/IL2 share the CRS drive by hero-derived hierarchy; IL6 is myeloid-amplified (small direct)
    sh={"cyto_IFNG":0.45,"cyto_TNF":0.32,"cyto_IL2":0.18,"cyto_IL6":0.05}
    d={"effector":eff,"exhaustion":0.0,"suppression":supp,"proliferation":0.0}
    for k,f in sh.items(): d[k]=crs*f
    # uncertainty (z_sd) carried per axis for the CI band
    sd={"effector":float(r["effector_z_sd"])*kE,"suppression":float(r["supp_z_sd"])*kS,
        "proliferation":0.0,"exhaustion":0.0}
    for k,f in sh.items(): sd[k]=float(r["crs_z_sd"])*kC*f
    # (3) prolif/exh magnitude from the lane's AUTHORITATIVE per_arm_prolif_exh_drive.csv (arm-specific,
    #     agonism-direction, with sd/q). kP scales onto the program rate-modifier basis. Timescale stays
    #     my hero k_on/k_off; veto gate unchanged -> no double-count.
    PE=_load_prolif_exh(); kP=1.0
    if PE is not None and lbl in PE.index:
        pr=PE.loc[lbl]
        d["proliferation"]=kP*float(pr["prolif_agon"]); d["exhaustion"]=kP*float(pr["exh_agon"])
        sd["proliferation"]=kP*float(pr["prolif_sd"]);  sd["exhaustion"]=kP*float(pr["exh_sd"])
    return d, sd
_PE_TBL=None
def _load_prolif_exh():
    global _PE_TBL
    if _PE_TBL is None:
        import pandas as pd
        f=_surrogate_path("per_arm_prolif_exh_drive.csv")
        if True:
            t=pd.read_csv(f)
            # prolif/exh table uses gene names (TNFRSF9 etc) directly; alias to receptor labels for lookup
            t["receptor"]=t["receptor"].map(lambda g:_ARM_ALIAS.get(g,g))
            _PE_TBL=t.set_index("receptor")
    return _PE_TBL

class PerCellSignaling:
    """Per-cell, per-program signaling state integrated over the PD time axis. Vectorized over a cell
    population. State S: (n_programs, n_cells). Driven by per-cell arm occupancy occ(t) in [0,1]."""
    def __init__(self, n_cells, arm, grn_shape=None, drive_source="lane"):
        self.arm=arm; self.n=n_cells; self.progs=PROGRAMS
        self.S=np.zeros((len(PROGRAMS),n_cells),float)
        self.peak=np.zeros((len(PROGRAMS),n_cells),float)   # track acute peak for decay phase
        if drive_source=="lane":
            d,sd=arm_drive_lane(arm); self.drive_sd={p:sd.get(p,0.0) for p in PROGRAMS}
        else:
            d=arm_drive(arm); self.drive_sd={p:0.0 for p in PROGRAMS}
        self.drive=np.array([[d[p]] for p in PROGRAMS])      # (n_prog,1) signed setpoint per program
        # GRN per-cell shape multiplier (heterogeneity); default 1.0 if not supplied
        self.shape=np.ones((len(PROGRAMS),n_cells)) if grn_shape is None else grn_shape
        self.kon=np.array([[RATES[p]["k_on_perhr"]] for p in PROGRAMS])
        self.koff=np.array([[RATES[p]["k_off_perhr"]] for p in PROGRAMS])
        self.acute=np.array([[1.0 if RATES[p]["acute"] else 0.0] for p in PROGRAMS])
        self.sfrac=np.array([[RATES[p]["sustained_frac"] for p in PROGRAMS]]).T
    def step(self, occ, dt_hr):
        """Advance signaling dt_hr with per-cell occupancy occ (len n, in [0,1]).
        Rise: dS = k_on*(setpoint - S)*dt ; setpoint = drive*shape*occ.
        Acute programs additionally decay toward sustained_frac*peak once occupancy-driven rise slows."""
        occ=np.asarray(occ,float)[None,:]                    # (1,n)
        setpoint=self.drive*self.shape*occ                   # (n_prog,n)
        # EXACT exponential integrator for first-order relaxation dS/dt=k_on*(setpoint-S):
        #   S(t+dt) = setpoint + (S-setpoint)*exp(-k_on*dt)   -- unconditionally stable for any dt.
        # (explicit Euler diverges here because k_on*dt can exceed 1 at the PD step size.)
        decay_on=np.exp(-self.kon*dt_hr)
        self.S = setpoint + (self.S-setpoint)*decay_on
        self.peak=np.maximum(self.peak,np.abs(self.S))
        # acute contraction: cytokines relax from peak toward sustained_frac*peak, also exact-exponential
        target_acute=np.sign(self.S)*self.sfrac*self.peak
        decay_off=np.exp(-self.koff*dt_hr)
        self.S = np.where(self.acute>0, target_acute+(self.S-target_acute)*decay_off, self.S)
        return self.S
    def program(self, name):
        return self.S[self.progs.index(name)]

if __name__=="__main__":
    print("=== data-derived program rates (from hero Rest/8hr/48hr) ===")
    for p in PROGRAMS:
        r=RATES[p]; print(f"  {p:14s} k_on={r['k_on_perhr']:.3f}/hr acute={r['acute']} k_off={r['k_off_perhr']:.4f}/hr sust_frac={r['sustained_frac']:.2f}")
    # simulate 4-1BB vs CD28 over 14 days at full occupancy, watch cytokine (acute) vs exhaustion (late)
    print("\n=== per-cell signaling trajectory (occ=1.0), 4-1BB vs CD28 ===")
    for arm in ["TNFRSF9","CD28"]:
        sig=PerCellSignaling(1,arm); tr=[]
        for step in range(int(14*24/0.5)):  # 14 days, 0.5hr steps
            sig.step(np.array([1.0]),0.5); 
            if step in [16,96,672]: tr.append((step*0.5/24, sig.program("cyto_IFNG")[0], sig.program("exhaustion")[0], sig.program("effector")[0]))
        print(f"  {arm}: "+" ".join(f"d{d:.0f}[IFNG={i:+.2f} exh={e:+.2f} eff={f:+.2f}]" for d,i,e,f in tr))
