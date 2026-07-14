"""coupled_percell_pd.py — whole-body per-cell PK WITH live per-cell PD (killing + cytokine + Treg
suppression) in EVERY organ. Extends CoupledPerCellPK (validated transport core, untouched) by attaching an
OrganPD to each organ and advancing it with that organ's per-cell interstitial drug g.C each step.

The PD is TCE-specific (CD3 x TAA trimer): effector arm = per-cell CD3E, target arm = per-cell antigen.
PD is OFF by default (pd_target=None) so transport-only runs are byte-identical to CoupledPerCellPK.
"""
import numpy as np, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from coupled_percell_pk import CoupledPerCellPK, load_tissue, AVO
import kinetic_rhoden_percell as _KRP
_NM_PER_COPY_PD=6.0/257000.0
from wholebody_percell import percell_rhoden_qss_Cvec
from wholebody_pd import OrganPD, CYTO_HIER

class CoupledPerCellPD(CoupledPerCellPK):
    def attach_pd(self, pd_target_col, KD_CD3_nM=40.0, KD_TAA_nM=1.45, costim_boost=0.0, cd3_col="CD3E",
                  costim_arm=None, KD_costim_nM=1.0, agents_dir=None, pd_kinetics=False, kin_params=None):
        """Build one OrganPD per organ. pd_target_col = the antigen column (per-cell TAA density);
        cd3_col = per-cell CD3 column (effector arm). costim_arm = costim receptor column (e.g. 'TNFRSF9')
        whose per-cell HPA/Glassman copies drive the per-cell signaling. Reads the SAME agent tables."""
        agents_dir = agents_dir or self._agents_dir
        self.pd={}; self.pd_target=pd_target_col; self.costim_arm=costim_arm
        # costim receptor copies live in a '<ARM>_copies' column (converted per-cell); orig-7 have no suffix
        costim_col=None
        if costim_arm is not None:
            costim_col = costim_arm+"_copies"
        t0=time.time()
        for o in self.organs:
            g=self.graphs[o]
            # reload per-cell CD3 + TAA aligned to this organ's finite-coord cells (same mask as transport)
            x,y,taa,labs,good=load_tissue(o, pd_target_col, agents_dir)
            _,_,cd3,_,_=load_tissue(o, cd3_col, agents_dir)
            n=len(g.x)
            if len(taa)!=n:   # safety: reload can differ if target/cd3 missing; pad/truncate to graph
                taa=np.resize(taa,n); cd3=np.resize(cd3,n); labs=np.resize(labs,n)
            Rcos=None
            if costim_col is not None:
                try:
                    _,_,Rcos,_,_=load_tissue(o, costim_col, agents_dir)
                    if len(Rcos)!=n: Rcos=np.resize(Rcos,n)
                except Exception:
                    Rcos=None   # organ lacks this costim column (e.g. brain overlay gap) -> legacy path
            _kp=kin_params or {}
            self.pd[o]=OrganPD(o, g.x, g.y, g.labs, cd3, taa,
                               KD_CD3_nM=KD_CD3_nM, KD_TAA_nM=KD_TAA_nM, costim_boost=costim_boost,
                               costim_arm=(costim_arm if Rcos is not None else None),
                               R_costim_percell=Rcos, KD_costim_nM=KD_costim_nM,
                               pd_kinetics=pd_kinetics, **_kp)
        self.pd_build_s=time.time()-t0
        return self

    def attach_heme_pd(self, pd_target_col, KD_CD3_nM=40.0, KD_TAA_nM=1.45, cd3_col="CD3E",
                       agents_dir=None, pd_kinetics=False, kin_params=None):
        """Additive HEME-malignancy PD compartment (circulating + marrow blasts, van Galen AML/MM/lymphoma).
        Plasma-driven (no extravasation barrier for circulating blasts) — the physical heme-vs-solid difference
        is the dense, no-ECM synapse geometry baked into heme_abm_cells.npz. k_death identical to organs."""
        agents_dir = agents_dir or self._agents_dir
        # route to the correct malignancy compartment by target antigen
        heme_npz="heme_abm_cells.npz"
        try:
            import json as _json
            for cand in (os.path.join(agents_dir,"..","params","heme_compartment_routing.json"),
                         os.path.join(agents_dir,"heme_compartment_routing.json"),
                         os.path.join(os.path.dirname(agents_dir),"heme_compartment_routing.json")):
                if os.path.exists(cand):
                    rt=_json.load(open(cand)).get("by_target",{})
                    if pd_target_col in rt: heme_npz=rt[pd_target_col]["npz"]
                    break
        except Exception: pass
        self.heme_npz=heme_npz
        x,y,taa,labs,good=load_tissue("heme_tumor", pd_target_col, agents_dir, heme_npz=heme_npz)
        _,_,cd3,_,_       =load_tissue("heme_tumor", cd3_col,      agents_dir, heme_npz=heme_npz)
        _kp=kin_params or {}
        self.heme_pd=OrganPD("heme_tumor", x, y, labs, cd3, taa,
                             KD_CD3_nM=KD_CD3_nM, KD_TAA_nM=KD_TAA_nM,
                             pd_kinetics=pd_kinetics, **_kp)
        self.heme_target=pd_target_col
        # ---- HEME TMDD SINK: emergent from per-cell Rhoden binding + internalization of complex with kint ----
        # Malignant blasts express the target, bind circulating drug, and internalize the complex -> a real
        # target-mediated clearance the plasma balance must feel, additive to the organ sinks. Magnitude is
        # anchored to a PHYSIOLOGICAL malignant burden pool (heme_pool_nM), distributed over the sampled blasts
        # by their own per-cell target copies -- identical convention to the organ antigen pools
        # (tot_copies = pool_nM * V * AVO/1e9). Without this anchor the sink would track the arbitrary ABM
        # sample size. Only malignant (target+) cells contribute; T-effector cells do not.
        self.heme_sink_on = bool(int(os.environ.get("HEME_TMDD_SINK","1")))
        _malig_name=None
        try:
            import json as _json2
            for cand in (os.path.join(agents_dir,"..","params","heme_compartment_routing.json"),
                         os.path.join(agents_dir,"heme_compartment_routing.json"),
                         os.path.join(os.path.dirname(agents_dir),"heme_compartment_routing.json")):
                if os.path.exists(cand):
                    _malig_name=_json2.load(open(cand)).get("by_target",{}).get(pd_target_col,{}).get("malignancy")
                    break
        except Exception: pass
        # Malignant burden as a TOTAL ANTIGEN AMOUNT (nmol), anchored to the homeostatic organ target pool for
        # this antigen scaled by a disease-burden multiple. Homeostatic amount = sum_organ(pool_nM * Vis_L),
        # the physiological antigen load already driving PK; a heavy myeloma/lymphoma carries a few-fold more.
        # This amount-based anchor keeps the sink physiological (sink ~ linear clearance at therapeutic C) and
        # units-consistent with the organ sinks, unlike a flat nM pool placed in V_pl. HEME_BURDEN_MULT overrides.
        _homeo_nmol=0.0
        try:
            import json as _json3
            for cand in (os.path.join(agents_dir,"Rtot_wholebody_final.json"),
                         os.path.join(os.path.dirname(agents_dir),"Rtot_wholebody_final.json"),
                         os.path.join(agents_dir,"..","params","Rtot_wholebody_final.json")):
                if os.path.exists(cand):
                    _rt=_json3.load(open(cand)).get("Rtot_nM",{}).get(pd_target_col,{})
                    _Vis={"spleen":0.0777,"bone":1.0,"large_int":0.226,"lung":0.10,"kidney":0.047,
                          "liver":0.30,"pancreas":0.016,"small_int":0.10,"heart":0.0088,"skin":0.10,
                          "adipose":0.28,"brain":0.10}
                    _homeo_nmol=sum(_rt.get(o,0)*_Vis.get(o,0.05) for o in _rt if isinstance(_rt.get(o),(int,float)))
                    break
        except Exception: pass
        _burden_mult=float(os.environ.get("HEME_BURDEN_MULT","3.0"))     # disease burden vs homeostatic antigen
        _burden_nmol=_burden_mult*_homeo_nmol if _homeo_nmol>0 else 0.0
        # KEEP the REAL per-cell IHC/Glassman copies (taa) in the Rhoden binding -- do NOT redistribute a pool,
        # which would flatten the heterogeneity. Instead scale only the CELL COUNT: the sampled malignant cells
        # are a representative draw; the true systemic burden carries _burden_nmol of target antigen total.
        # count_scale = (burden antigen copies) / (sum of sampled malignant per-cell copies) -> multiplies the
        # summed sampled sink up to the real burden while every cell binds at its own real receptor number.
        self.heme_R_sink=np.asarray(taa,float)                    # REAL per-cell IHC copies (malignant = target+)
        _malig_mask=self.heme_R_sink>0
        _sampled_copies=float(self.heme_R_sink[_malig_mask].sum())
        _burden_copies=_burden_nmol*1e-9*AVO
        self.heme_count_scale=(_burden_copies/_sampled_copies) if (_burden_copies>0 and _sampled_copies>0) else 0.0
        self.heme_sink_KD=float(KD_TAA_nM); self.heme_sink_narm=int(self.n_arm)
        self.heme_sink_kint=float(self.kint); self.heme_burden_nmol=_burden_nmol
        # kinetic receptor state for the heme sink (synapse-basis nM; turnover + avidity like the organ engine)
        self._heme_Ag0=self.heme_R_sink*_NM_PER_COPY_PD; self._heme_Ag=self._heme_Ag0.copy()
        self._heme_BAg1=np.zeros_like(self._heme_Ag0); self._heme_Bdbl=np.zeros_like(self._heme_Ag0)
        self._heme_kdeg=(self.kdeg if getattr(self,'kdeg',None) is not None else 0.5)
        _kon_hd=(self.kon1 if getattr(self,'kon1',None) is not None else 8.64e-3)
        self._heme_kon=_kon_hd; self._heme_koff=(self.koff1 if getattr(self,'koff1',None) is not None else _kon_hd*self.heme_sink_KD)
        self.heme_homeo_nmol=_homeo_nmol; self.heme_malignancy=_malig_name
        return self

    def attach_blood_pd(self, pd_target_col, KD_CD3_nM=40.0, KD_TAA_nM=1.45, cd3_col="CD3E",
                        agents_dir=None, pd_kinetics=False, kin_params=None):
        """Always-present NORMAL circulating blood compartment (Tabula Sapiens blood). Plasma-driven, like
        heme. Contributes: (1) TMDD sink from normal B/T/myeloid binding+internalization (real per-cell IHC
        copies, count-scaled to physiological blood lineage counts); (2) on-target NORMAL-cell depletion
        (tox readout, e.g. circulating B-cell depletion); (3) CRS from engaged circulating T cells."""
        agents_dir = agents_dir or self._agents_dir
        x,y,taa,labs,good=load_tissue("blood", pd_target_col, agents_dir)
        _,_,cd3,_,_       =load_tissue("blood", cd3_col,      agents_dir)
        _kp=kin_params or {}
        self.blood_pd=OrganPD("blood", x, y, labs, cd3, taa,
                              KD_CD3_nM=KD_CD3_nM, KD_TAA_nM=KD_TAA_nM,
                              pd_kinetics=pd_kinetics, **_kp)
        self.blood_target=pd_target_col
        # ---- BLOOD TMDD SINK (emergent Rhoden + kint on real per-cell copies, count-scaled to real blood) ----
        # Normal circulating target+ cells bind circulating drug and internalize the complex. Burden here is
        # the REAL physiological cell count (5 L blood x lineage densities), not a disease burden. The ABM is a
        # sample; count_scale lifts the sampled per-lineage sink up to the true circulating count for that lineage.
        self.blood_sink_on = bool(int(os.environ.get("BLOOD_TMDD_SINK","1")))
        _bl=np.load(os.path.join(agents_dir,"agents_blood.npz") if os.path.exists(os.path.join(agents_dir,"agents_blood.npz"))
                    else os.path.join(agents_dir,"..","agents_blood.npz"), allow_pickle=True)
        self.blood_R_sink=np.asarray(taa,float)                 # REAL per-cell IHC copies (target+ normal cells)
        # per-lineage count_scale = real_circulating_count / sampled_count  (respects that ABM is a subsample)
        _labmap=list(_bl["real_lineage_names"]); _rc=np.asarray(_bl["real_lineage_counts"],float)
        _real={_labmap[i]:_rc[i] for i in range(len(_labmap))}
        _labs_good=np.asarray(_bl["labels"])[good]
        _scale=np.ones(len(taa))
        for _l in np.unique(_labs_good):
            _m=_labs_good==_l; _ns=int(_m.sum())
            if _ns>0 and str(_l) in _real: _scale[_m]=_real[str(_l)]/_ns
        self.blood_count_scale=_scale                            # per-cell multiplier to real circulating count
        # MYELOID IL-6: hand each blood myeloid agent ITS OWN per-lineage real-count scale, so monocytes are
        # lifted by the REAL circulating monocyte count (not a population average, which would apply a
        # lymphocyte-weighted factor to monocytes). This scale is CELL-COUNT based and DRUG-INDEPENDENT --
        # unlike the antigen-derived graphs[o].count_scale (measured drug-dependent 0.18x-5.34x). Blood
        # monocytes are the dominant CRS IL-6 source (Norelli PMID 29808007), so this is the key anchor.
        # ---- BLOOD myeloid do NOT sustained-contact-activate (MECHANISTIC, not an approximation) ----------
        # IL-6 induction requires SUSTAINED CD40L-CD40 contact ("IL-6 induction and myeloid activation require
        # proximity of CAR T cells and myeloid cells", Giavridis PMID 29808005). In FLOWING blood leukocytes are
        # ~50 um apart AND in motion -> contacts are transient collisions, not synapses. Sustained contact only
        # occurs where myeloid are ADHERENT or RESIDENT. Per PMID 3944542, 60% of blood monocytes are
        # MARGINATING (adherent to endothelium) and only 40% freely circulating; the marginating/extravasated
        # monocytes ARE the myeloid already resident in the organ ABMs (spleen 54,373 macrophages, marrow 3,179).
        # So the contact-competent myeloid are counted in TISSUE, not in the flowing-blood pool -- adding blood
        # here would DOUBLE-COUNT them and, worse, the blood ABM's coordinates are a SYNTHETIC 2D grid (mean
        # spacing ~15 um) whose spatial contact is non-physical: it activates 98.6% of 2.0e9 monocytes and
        # yields 61,874 pg/mL (measured this session) vs a severe-CRS ceiling of ~10-20k. Setting the blood
        # myeloid scale to 0 removes a GEOMETRY ARTIFACT; it does not remove biology.
        try:
            self.blood_pd.myeloid.set_count_scale(0.0)
        except Exception as _e:
            print(f"[myeloid-IL6] blood myeloid gate not applied ({_e})", flush=True)
        self.blood_sink_KD=float(KD_TAA_nM); self.blood_sink_narm=int(self.n_arm); self.blood_sink_kint=float(self.kint)
        self._blood_Ag0=self.blood_R_sink*_NM_PER_COPY_PD; self._blood_Ag=self._blood_Ag0.copy()
        self._blood_BAg1=np.zeros_like(self._blood_Ag0); self._blood_Bdbl=np.zeros_like(self._blood_Ag0)
        self._blood_kdeg=(self.kdeg if getattr(self,'kdeg',None) is not None else 0.5)
        _kon_bd=(self.kon1 if getattr(self,'kon1',None) is not None else 8.64e-3)
        self._blood_kon=_kon_bd; self._blood_koff=(self.koff1 if getattr(self,'koff1',None) is not None else _kon_bd*self.blood_sink_KD)
        return self

    def simulate_pd(self, schedule, tsim, route="IV", F_sc=0.6, ka_sc=0.25, inf_rate=0.0, inf_dur=0.0,
                    iv_inf_h=2.0, rec=400, k_death=1.0, pd_every=1, progress_tag="", progress_every=0, snap_times=None):
        """Same transport as simulate(), but advances per-cell PD each step using g.C. Returns the transport
        dict PLUS pd time-courses: kill_frac[organ], cyto[organ][species], and systemic cytokine sum."""
        if not hasattr(self,"pd"): raise RuntimeError("call attach_pd(...) first")
        dt=self.dt; nstep=int(tsim/dt)+1
        A_pl=0.0; A_ly=0.0; A_sc=0.0
        sched=sorted(schedule); si=0
        _iv_infusions=[]   # active IV infusions: [t_start_d, t_end_d, rate_nmol_per_day]
        t_rec=[]; C_rec=[]
        snap_times=sorted(snap_times) if snap_times else []; self.snaps={}; _snap_i=0   # spatial snapshots (x/y/C/bound/R/labs/surv) at requested days
        import os as _os
        self._motility_on = bool(int(_os.environ.get("IMMUNE_MOTILITY","0")))
        self._motility_stride = int(_os.environ.get("MOTILITY_STRIDE","100"))   # every ~1 day at dt=0.01
        self._motility_speed = float(_os.environ.get("MOTILITY_SPEED_UM_MIN","5.0"))
        self._motility_chemotax = float(_os.environ.get("MOTILITY_CHEMOTAX","0.4"))
        pd_t=[]; kill_rec={o:[] for o in self.organs}; heme_kill_rec=[]; blood_kill_rec=[]; cyto_rec={o:{c:[] for c in CYTO_HIER} for o in self.organs}
        sys_cyto_rec={c:[] for c in CYTO_HIER}
        sys_cyto_rate_rec={c:[] for c in CYTO_HIER}   # INSTANTANEOUS systemic rate (clinical-peak comparable)
        # MECHANISTIC plasma IL-6 (pg/mL, physical units — NOT a scaled rate). Literature-only constants.
        from myeloid_il6 import PlasmaIL6 as _PlasmaIL6
        self._plasma_il6 = _PlasmaIL6()
        self._il6_prod_pg_hr = 0.0
        il6_plasma_rec=[]; il6_prod_rec=[]
        il6_prod_organ_rec={o:[] for o in self.organs}   # UNSCALED per-organ production -> post-hoc census rescale
        import time as _time, sys as _sys
        _t0=_time.time(); _last_day=-1e9
        for k in range(nstep):
            t=k*dt
            # ---- optional tagged progress: construct + sim-day + per-organ kill snapshot ----
            if progress_every and (k % progress_every == 0):
                _ks={o:self.pd[o].summary()['kill_frac'] for o in self.organs} if hasattr(self,'pd') else {}
                _kstr=" ".join(f"{o[:5]}={_ks[o]:.2f}" for o in self.organs) if _ks else ""
                _mk=(sum(_ks.values())/len(_ks)) if _ks else 0.0
                print(f"[{progress_tag or 'PD'}] day {t:5.2f} ({k:5d}/{nstep}) Cpl={A_pl/self.V_pl:7.3f}nM meanKill={_mk:.3f} | {_kstr} | {_time.time()-_t0:5.0f}s",flush=True)
            while si<len(sched) and sched[si][0]<=t+1e-9:
                mg=sched[si][1]
                if route=="SC":
                    A_sc+=mg/self.mw*1e3
                elif iv_inf_h and iv_inf_h>0:
                    # real IV infusion: spread this dose over iv_inf_h hours as a constant-rate input.
                    # MASS-EXACT: track remaining nmol, clamp per-step delivery so total = exactly
                    # (mg/mw*1e3) regardless of dt/dur_d alignment (fixes ~1.2x overage when dur_d
                    # is not an integer multiple of dt).
                    dur_d=iv_inf_h/24.0
                    tot_nmol=mg/self.mw*1e3
                    _iv_infusions.append([t, t+dur_d, tot_nmol/dur_d, tot_nmol])  # [t0,t1,rate,remaining]
                else:
                    A_pl+=mg/self.mw*1e3
                si+=1
            # deliver active IV infusions into plasma this step (clamp to remaining -> mass-exact)
            iv_in=0.0
            for _inf in _iv_infusions:
                if _inf[0]<=t<_inf[1] and _inf[3]>0.0:
                    step_dose=min(_inf[2]*dt, _inf[3])
                    iv_in+=step_dose/dt
                    _inf[3]-=step_dose
            A_pl+=iv_in*dt
            C_pl=A_pl/self.V_pl
            J_sc=ka_sc*A_sc; A_sc-=dt*J_sc
            infn=(inf_rate/self.mw*1e3) if (inf_rate>0 and t<=inf_dur) else 0.0
            tot_drain=0.0; tot_extrav=0.0
            def _one(o):
                g=self.graphs[o]
                PS_ex=self.k_dist*self.L[o]*(1.0-self.sigV[o]); PS_ly=self.k_dist*self.L[o]*(1.0-self.sigL)
                sink,drain,J_extrav,S,D=g.step(C_pl,self.Q[o],PS_ex,PS_ly,self.KD,self.n_arm,self.kint)
                return o,drain,J_extrav,S
            results=list(self.pool.map(_one,self.organs)) if self.pool else [ _one(o) for o in self.organs]
            Sd={}
            for o,drain,J_extrav,S in results:
                tot_extrav+=J_extrav; tot_drain+=drain; Sd[o]=S
            # ---- live per-cell PD, every pd_every steps, driven by each organ's per-cell g.C ----
            if k % pd_every == 0:
                for o in self.organs:
                    self.pd[o].step(self.graphs[o].C, dt*pd_every, k_death=k_death)
                # ---- MECHANISTIC PLASMA IL-6 (replaces the fitted mosun-anchored IL6_SCALE) -------------
                # Sum every compartment's per-cell MYELOID production (pg/hr), count_scale-lifted to the
                # physiological population, then integrate the MEASURED first-order IL-6 clearance
                # (kdeg 0.20/hr, t1/2 ~3.5 h; PMID 31268236). This converts a production RATE into a
                # CONCENTRATION -- the step the engine never had, which the fitted scale was standing in for.
                # count_scale MUST be a MYELOID CELL-COUNT scale (sampled myeloid -> physiological myeloid),
                # NOT graphs[o].count_scale (= tot_antigen_copies/sampled_antigen_copies, a TARGET-cell scale
                # that DIFFERS PER DRUG -> would scale the same monocytes differently for mosun vs elran = a
                # per-drug artifact that would corrupt the counterscreen). Myeloid scale is a TISSUE property:
                # drug-independent, identical across all molecules. Blood uses its cell-count-based
                # blood_count_scale (real circulating count / sampled count) which is already correct.
                _il6_prod = 0.0
                for o in self.organs:
                    _cs = float(getattr(self.pd[o], "myeloid_count_scale", 1.0) or 1.0)   # tissue cellularity
                    _il6_prod += _cs * float(getattr(self.pd[o], "il6_prod_pg_hr", 0.0))
                _b = getattr(self, "blood_pd", None)
                if _b is not None:   # blood myeloid uses blood_pd's SCALAR myeloid_count_scale (tissue property),
                    # NOT self.blood_count_scale (a per-CELL antigen-sink array -> float() crashes; wrong scale anyway)
                    _il6_prod += float(getattr(_b, "myeloid_count_scale", 1.0) or 1.0) * float(getattr(_b, "il6_prod_pg_hr", 0.0))
                _h = getattr(self, "heme_pd", None)
                if _h is not None:   # heme_count_scale is ANTIGEN-derived (malignant burden) -> NOT valid for
                    _il6_prod += float(getattr(_h, "myeloid_count_scale", 1.0) or 1.0) * float(getattr(_h, "il6_prod_pg_hr", 0.0))
                # ---- PHYSICAL-CEILING GUARD (added 2026-07-13 after this exact bug shipped) --------------
                # IL-6 production CANNOT exceed (every secretor macrophage in the census, fully activated)
                # x (the measured per-cell secretion rate). If it does, a scale is being misapplied -- which
                # is precisely what happened: raw census COUNTS were used as a multiplicative SCALE, putting
                # production 942x above what the cells can physically secrete. The model was reporting IL-6
                # that its own cells could not make, and nothing noticed for a full day of runs.
                #
                # This is a CONSERVATION CHECK. It is cheap, it is unambiguous, and it would have caught the
                # bug on the first run. Warn loudly (once) rather than crash mid-run.
                if not getattr(self, "_il6_ceiling_warned", False):
                    _ceil = getattr(self, "_il6_ceiling_pg_hr", None)
                    if _ceil is None:
                        try:
                            from myeloid_il6 import S_MAX_PG_PER_HR as _SM, SECRETOR_FRACTION as _SF
                            _tot_mye = 0.0
                            for _o in self.organs:
                                _c = getattr(self.pd[_o], "myeloid_census_count", None)
                                _tot_mye += float(_c) if _c else float(getattr(self.pd[_o].myeloid, "n_myeloid", 0))
                            _ceil = _tot_mye * _SF * _SM
                            self._il6_ceiling_pg_hr = _ceil
                        except Exception:
                            _ceil = None
                    if _ceil and _il6_prod > _ceil * 1.001:
                        import sys as _s
                        _s.stderr.write(
                            f"[IL6-CEILING] *** PRODUCTION EXCEEDS THE PHYSICAL MAXIMUM ***\n"
                            f"    reported {_il6_prod:,.0f} pg/hr  vs  ceiling {_ceil:,.0f} pg/hr "
                            f"({_il6_prod/_ceil:.1f}x)\n"
                            f"    The census's secretor macrophages, ALL fully activated, at the measured\n"
                            f"    per-cell secretion rate, cannot make this much IL-6. A scale is being\n"
                            f"    misapplied (a COUNT used as a SCALE?). DO NOT TRUST THIS RUN'S IL-6.\n")
                        self._il6_ceiling_warned = True
                self._il6_prod_pg_hr = _il6_prod
                self._plasma_il6.step(dt*pd_every, _il6_prod)
                # ---- immune motility: T cells migrate + synapse graph rebuilt on a stride (opt-in) ----
                if self._motility_on and (k % self._motility_stride == 0):
                    _mv_dt=dt*self._motility_stride
                    for o in self.organs:
                        self.pd[o].move_immune(_mv_dt, speed_um_per_min=self._motility_speed, chemotax=self._motility_chemotax)
                    if hasattr(self,"heme_pd"): self.heme_pd.move_immune(_mv_dt, speed_um_per_min=self._motility_speed, chemotax=0.0)  # heme well-mixed: no chemotaxis
                    if hasattr(self,"blood_pd"): self.blood_pd.move_immune(_mv_dt, speed_um_per_min=self._motility_speed, chemotax=0.0)
                if hasattr(self,"heme_pd"):
                    # circulating heme blasts see plasma drug directly (nM = A_pl/V_pl)
                    C_heme=np.full(len(self.heme_pd.x), C_pl, float)
                    self.heme_pd.step(C_heme, dt*pd_every, k_death=k_death)
                    # HEME TMDD SINK (emergent): bound complex on malignant burden internalized at kint.
                    # S = bound monovalent, D = bivalent-bridged; internalized copies/day = kint*(S+2D).
                    if getattr(self,"heme_sink_on",False) and getattr(self,"heme_count_scale",0.0)>0:
                        # bind at REAL per-cell copies, internalize at kint, then scale sampled sum up to burden.
                        # Weight by LIVE-blast survival: dead blasts no longer internalize drug, so the TMDD sink
                        # shrinks as the malignant burden is depleted (prevents a static full-burden sink from
                        # over-draining plasma into a non-physiological terminal cliff at low dose).
                        surv_h=np.exp(-self.heme_pd.kill_hazard) if hasattr(self.heme_pd,"kill_hazard") else np.ones(len(self.heme_R_sink))
                        # KINETIC heme sink: per-cell Rhoden binding with receptor turnover, state carried across steps.
                        AgEFF_h=_KRP.geo_ageff_nM(self.heme_R_sink,8.0,12.5) if self.heme_sink_narm>=2 else 0.0
                        self._heme_Ag,self._heme_BAg1,self._heme_Bdbl,heme_flux=_KRP.rhoden_samecell_bivalent_step(
                            np.full(len(self.heme_R_sink),C_pl,float),self._heme_Ag,self._heme_BAg1,self._heme_Bdbl,
                            self._heme_Ag0,self._heme_kon,self._heme_koff,AgEFF_h,self._heme_kdeg,self.heme_sink_kint,dt*pd_every)
                        # flux (nM/day, synapse basis) -> copies/cell/day -> nmol/day, weighted by live-blast survival
                        self._heme_sink_nmol_day=self.heme_count_scale*float(((heme_flux/_NM_PER_COPY_PD)*surv_h).sum())/AVO*1e9
                    else:
                        self._heme_sink_nmol_day=0.0
                if hasattr(self,"blood_pd"):
                    # normal circulating blood cells see plasma drug directly (always-present compartment)
                    C_blood=np.full(len(self.blood_pd.x), C_pl, float)
                    self.blood_pd.step(C_blood, dt*pd_every, k_death=k_death)
                    # BLOOD TMDD SINK (emergent): per-cell Rhoden bind + kint, per-lineage count-scaled to real blood
                    if getattr(self,"blood_sink_on",False):
                        AgEFF_b=_KRP.geo_ageff_nM(self.blood_R_sink,8.0,12.5) if self.blood_sink_narm>=2 else 0.0
                        self._blood_Ag,self._blood_BAg1,self._blood_Bdbl,blood_flux=_KRP.rhoden_samecell_bivalent_step(
                            np.full(len(self.blood_R_sink),C_pl,float),self._blood_Ag,self._blood_BAg1,self._blood_Bdbl,
                            self._blood_Ag0,self._blood_kon,self._blood_koff,AgEFF_b,self._blood_kdeg,self.blood_sink_kint,dt*pd_every)
                        self._blood_sink_nmol_day=float((self.blood_count_scale*(blood_flux/_NM_PER_COPY_PD)).sum())/AVO*1e9
                    else:
                        self._blood_sink_nmol_day=0.0
            # ---- spatial snapshot: full per-cell fields at requested days (for 4-panel overlay + kill maps) ----
            if _snap_i < len(snap_times) and t >= snap_times[_snap_i]-1e-9:
                tsnap=snap_times[_snap_i]; self.snaps[tsnap]={}
                # bound per-cell (nM): transport S is bound copies/cell -> nM via pericellular volume (Vis/n)
                for o in self.organs:
                    g=self.graphs[o]; pdo=self.pd.get(o)
                    vcell=max(self.Vis[o]/g.n,1e-18)
                    bound_nM=(np.asarray(Sd.get(o,np.zeros(g.n)))/AVO/vcell*1e9)
                    surv=np.exp(-pdo.kill_hazard) if pdo is not None else None
                    self.snaps[tsnap][o]=dict(x=g.x.copy(), y=g.y.copy(), C=g.C.copy(),
                        bound_nM=bound_nM, R=g.R.copy(), labs=np.asarray(g.labs),
                        surv=(surv.copy() if surv is not None else None),
                        is_target=(pdo.is_target.copy() if pdo is not None else None))
                _snap_i+=1
            for o,(Lo,sVo,Viso) in self.wm.items():
                Cis=self.A_wm[o]/max(Viso,1e-9)
                Jex=self.k_dist*Lo*(1.0-sVo)*C_pl; Jre=self.k_dist*Lo*(1.0-self.sigL)*Cis
                self.A_wm[o]=max(self.A_wm[o]+dt*(Jex-Jre),0.0); tot_extrav+=Jex; tot_drain+=Jre
            dA_pl=infn+F_sc*J_sc+self.k_lymph_return*A_ly-self.k_cat*A_pl-tot_extrav-getattr(self,'_heme_sink_nmol_day',0.0)-getattr(self,'_blood_sink_nmol_day',0.0)  # + heme blast + normal-blood TMDD sinks
            dA_ly=tot_drain-self.k_lymph_return*A_ly
            A_pl=max(A_pl+dt*dA_pl,0.0); A_ly=max(A_ly+dt*dA_ly,0.0)
            if k % max(1,nstep//rec)==0:
                t_rec.append(t); C_rec.append(A_pl/self.V_pl*self.mw/1e3)
                pd_t.append(t)
                sys_c={c:0.0 for c in CYTO_HIER}; sys_cr={c:0.0 for c in CYTO_HIER}
                if hasattr(self,"heme_pd"):
                    hs=self.heme_pd.summary(); heme_kill_rec.append(hs['kill_frac'])
                if hasattr(self,"blood_pd"):
                    bs=self.blood_pd.summary(); blood_kill_rec.append(bs['kill_frac'])
                for o in self.organs:
                    s=self.pd[o].summary(); kill_rec[o].append(s['kill_frac'])
                    for c in CYTO_HIER:
                        cyto_rec[o][c].append(s['cyto'][c]); sys_c[c]+=s['cyto'][c]
                        sys_cr[c]+=s['cyto_rate'][c]
                if hasattr(self,"blood_pd"):
                    _bcy=self.blood_pd.summary()
                    for c in CYTO_HIER:
                        sys_c[c]+=_bcy['cyto'][c]; sys_cr[c]+=_bcy['cyto_rate'][c]   # circulating-T CRS
                for c in CYTO_HIER:
                    sys_cyto_rec[c].append(sys_c[c]); sys_cyto_rate_rec[c].append(sys_cr[c])
                # MECHANISTIC IL-6: record the PHYSICAL plasma concentration (pg/mL) + production (pg/hr)
                il6_plasma_rec.append(float(self._plasma_il6.C))
                il6_prod_rec.append(float(self._il6_prod_pg_hr))
                # PER-ORGAN production at count_scale=1 (UNSCALED). Plasma IL-6 is EXACTLY LINEAR in the
                # per-organ myeloid count_scale, so recording this lets the citation-gated organ myeloid census
                # be applied ANALYTICALLY afterwards -- the run does NOT have to be repeated when it lands.
                for _o in self.organs:
                    il6_prod_organ_rec[_o].append(float(getattr(self.pd[_o], "il6_prod_pg_hr", 0.0)))
        out=dict(t=np.array(t_rec), Cplasma_ugml=np.array(C_rec), pd_t=np.array(pd_t),
                 il6_plasma_pgml=il6_plasma_rec, il6_prod_pg_hr=il6_prod_rec,   # MECHANISTIC (no fitted scale)
                 il6_prod_organ_pg_hr=il6_prod_organ_rec,                        # UNSCALED per-organ -> census rescale
                 myeloid_count_scale={o: float(getattr(self.pd[o],"myeloid_count_scale",1.0)) for o in self.organs},
                 myeloid_n={o: int(getattr(self.pd[o].myeloid,"n_myeloid",0)) for o in self.organs},
                 kill_frac=kill_rec, cyto=cyto_rec, sys_cyto=sys_cyto_rec,
                 sys_cyto_rate=sys_cyto_rate_rec, build_s=self.build_s,
                 n_target={o:int(self.pd[o].tgtidx.size) for o in self.organs}, heme_kill=heme_kill_rec, heme_ntgt=(int(self.heme_pd.tgtidx.size) if hasattr(self,'heme_pd') else 0), blood_kill=blood_kill_rec, blood_ntgt=(int(self.blood_pd.tgtidx.size) if hasattr(self,'blood_pd') else 0))  # target-cell count per organ (for count-weighted depletion)
        return out
