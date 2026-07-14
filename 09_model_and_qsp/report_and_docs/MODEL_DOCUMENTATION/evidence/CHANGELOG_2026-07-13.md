# Model fixes — 2026-07-13 (PK/PD mechanism audit)

All edits made to model/engine/ (canonical) and mirrored to model/rundir/handoff/ (md5-identical run copy).

## PK mechanism (diagnosed vs FLMOON digitized concentration-time curves)
1. fFcRn by construct — run_tce_pd_reval.py ENG config. Fc-competent IgG1/IgG4 TCE = 0.89 (was 0.70,
   which under-recycled and gave t1/2 ~4-8d vs clinical 16-19d). fFcRn is a MOLECULE property (subtype +
   Fc presence); CLup stays fixed physiology. Blinatumomab (BiTE, no Fc) = 0.0.
2. Per-target internalization kint — run_tce_pd_reval.py antigen_kint(tgt). Pulled from
   params/antigen_kinetics_table.json (CD20 0.02, CD19 0.05, BCMA 2.0, GPRC5D 0.2, DLL3 0.3, ...) with a
   gene-symbol alias map (MS4A1->CD20, FOLH1->PSMA, IL3RA->CD123, ERBB2->HER2). REPLACES a hardcoded 0.15
   that was 7.5x too fast for CD20 -> over-drained plasma via the organ TMDD sink.
3. Per-dose IV infusion — coupled_percell_pd.py simulate_pd(iv_inf_h=2.0). Each IV dose spread over a 2h
   infusion window (not instantaneous bolus).
4. Regimen — params/regimen_schedules_final.json. mosunetuzumab C2D1 (day 21) corrected 30 -> 60 mg
   (FDA Lunsumio: C1 1/2/60, C2D1 60, C3+ 30 Q3W).

## Compartment wiring
5. Solid tumor = spatial ABM organ. run_organs = ORGANS + ['tumor'] when cancer_type=='solid'. Transport
   (ECM + BEC/LEC + per-cell Rhoden sink) via load_tissue('tumor'); kill via attach_pd. Dead
   attach_tumor_pd guard removed. Output kill loops use m.organs (records tumor kill).
6. Output completeness — added per-organ cytokine time-courses (cyto_organ). Saves every organ/tumor/
   heme/blood kill trajectory + all cytokine species + PK curve.

## Retired
- Old compartmental TMDD PK path (unified_pbpk_pd, emergent_tmdd_engine, percell_binding) -> model/_archived_2026-07-13/.
  PK now runs through the SAME per-cell Rhoden engine as PD.

## Validation status (mosunetuzumab)
- fFcRn fix: IL6 peak 1863 -> 618 (clinical 570). GOOD.
- kint fix (CD20 0.02): re-run in progress to confirm t1/2 climbs toward 16-19d.

## PK α-phase / organ CD20 tracking — INVESTIGATED (2026-07-13)
mosun kint-fix result: t1/2 8.0->12.1d (target 16-19), IL6 619 (~570), heme depl 0.999, Cmax 35 (high),
matched-timepoint AFE vs FLMOON-2 = 2.41. Low doses (d0-7) match 1.3-1.7x; 60mg doses over 2.8-3.2x.
Root of over-prediction = weak α-phase (13% drop day+1 vs real IgG 40-60%).

FINDING (organ receptor per-cell): coupled_percell_pk.py L63-64 sets Rpercell = (pool_nM*Vis*AVO/1e9)*(ag/ag_sum).
This distributes the WHOLE-ORGAN receptor amount across the SAMPLED cells, so per-cell R is inflated by
(real_cell_count / sampled_count) ~1e5x (spleen R_med 3e10, should be ~1e5 copies/cell).
- CONSERVED: organ receptor TOTAL = pool_nM*Vis is physiological (whole-body CD20 14.4 nmol ~ homeostatic 23.6).
  Masking is correct (spleen 145654 R+ = real B-cell count; liver 258, heart 1 — CD20 on B cells only).
- The PK SINK (Σ kint*bound) is therefore ~physiological in the linear regime.
- CAVEAT: per-cell Rhoden bound = R*C/(C+KD) is NONLINEAR in R; concentrating receptor onto sampled cells
  distorts saturation vs the true cell count. This is a candidate contributor to the weak α-phase / Cmax-high.
- Distribution DOES occur (free interstitial 0.1-1 nM across all 11 organs day-1) — organ system IS wired.
STATUS: AFE 2.4 is acceptable for a corrected first pass; the per-cell-R aggregation is flagged for the
PK refinement step (may need to normalize per-cell R to physiological copies while preserving the sink total,
or account count-scale like heme/blood). NOT blocking the tumor builds.

## FULL-BODY count_scale UNIFICATION (2026-07-13)
Unified receptor handling across ALL compartments (per RECEPTORS_PER_CELL_METHOD.md Glassman method):
every compartment now keeps REAL physiological per-cell receptor copies in the Rhoden binding and applies a
scalar count_scale = (physiological pool copies)/(Σ sampled per-cell copies) to lift the sampled ABM to the
true population count.
- ORGANS + TUMOR: coupled_percell_pk.py — Rpercell=real ag copies (was pool-redistributed); g.count_scale set;
  TissueGraph.step multiplies internalization by count_scale (local loss + plasma drain both scaled = mass self-consistent).
- HEME: heme_count_scale (burden-anchored) — already this pattern.
- BLOOD: blood_count_scale (per-lineage real/sampled) — already this pattern.
Spatial grid UNCHANGED: real (x,y), real per-cell R, kNN diffusion, BEC-in/LEC-out, ECM all untouched;
count_scale is a scalar on the internalization amount only. Monovalent: numerically identical to old
pool-redistribution (count_scale*ag = old Rpercell). Bivalent (glofitamab, ternary nonlinear in R): now
physically correct (was distorted by concentrating receptor onto sampled cells) — fixes the glofit PK outlier.

## AUDITOR FIXES (2026-07-13)
1. [CRASH FIX] run_tce_pd_reval.py arr(): built Q/L/sigV/Vis/Vv dicts over the fixed 11-organ ORGANS list,
   but solid runs use run_organs=ORGANS+['tumor'] -> Vis['tumor'] KeyError at CoupledPerCellPD construction
   (before any PD). Never caught because no solid run had been launched (only heme mosun). FIXED:
   arr() now iterates ALL pb.names (15 tissues incl tumor idx 14); solid run_organs resolves. Verified.
2. [DOC] Corrected "KDTree neighbors (k=6)" -> the spatial graph uses sklearn NearestNeighbors.kneighbors
   (wholebody_percell.py line 16/98), not scipy KDTree. kNN-diffusion conclusion unchanged; impl detail corrected.

## INFUSION MASS-CONSERVATION FIX (2026-07-13) — user-caught
SYMPTOM (user): "C0 starting way too high", "initial rate of infusion too high", "still looks like bolus".
ROOT CAUSE: IV infusion in simulate_pd delivered rate*dt at every step where t0<=t<t0+dur_d. When dur_d
(=iv_inf_h/24, e.g. 2h=0.0833d) is NOT an integer multiple of dt(=0.02), the '<' boundary fires an extra
step: 0.0833d spans 4.167 steps but 5 discrete steps fire (t=0,.02,.04,.06,.08) -> delivers rate*dt*5 =
8.22 nmol vs correct 6.85 nmol = 1.20x OVERDOSE on EVERY IV dose. Traced: 1mg dose gave C=0.53->correct 0.32.
FIX: mass-exact infusion — track remaining nmol per infusion [t0,t1,rate,remaining]; deliver
min(rate*dt, remaining) each step, decrement remaining. Total delivered = exactly mg/mw*1e3 regardless of
dt/dur_d alignment. Verified: 1mg -> C ramps 0.077->0.3226 ug/mL (exact 1mg/3.1L), total 6.8493 nmol exact.
IMPACT: affects ALL IV doses (mosun, glofit). ~1.2x per-dose reduction. Does NOT alone close the 60mg 3-4x
peak overshoot (that's the shared alpha-phase distribution rate); but it is a real mass-conservation bug fixed.

## PK VALIDATION-TARGET CORRECTION (2026-07-13) — critical
FINDING: the mosunetuzumab "FLMOON-2/-3" digitized curves (PMC12932324 / 10147_2025_2957) are
SUBCUTANEOUS mosun (Japanese trial, fixed step-up 5 mg C1D1, 45 mg C1D8+), NOT the IV regimen the
model runs (1/2/60 mg -> 30 Q3W). The earlier "60mg peaks overshoot 3-4x / weak alpha-phase" diagnosis
was comparing IV-model output against SC-clinical data — an INVALID route+dose mismatch. The catalog
has NO digitized IV concentration-time curve for mosunetuzumab (only exposure-response and covariate plots).
CONSEQUENCE: mosun is NOT a clean PK validation gate. Correct cohort sizes are FLMOON-2 n=17, FLMOON-3 n=5
(prior overlay mislabeled digitized-point counts 12/13 as cohort sizes — fixed).
RESOLUTION: pivot primary PK validation to TECLISTAMAB — has a rich route-matched SC C-t curve
(PMC13222323, 1.5 mg/kg QW, 0-40 weeks, 60 pts, ug/mL) + early IV dose-proportionality (PMC9345835).
Model runs teclistamab SC 0.06/0.3/1.5 mg/kg (4.2/21/105 mg x70kg) — route-matched. Running tecli to 200d.
Other route-matched C-t curves available: elranatamab (PMC10579053), talquetamab (PMC12380926),
epcoritamab (PMC13338999). All SC, matching model routes.

## SOLID-TUMOR ABM BUILDS (2026-07-13) — 4 cancers (prostate/lung/ovary/skin)
Generalized CRC pipeline into build_cancer_abm.py (tumor_builds/). Per cancer: unzip Xenium _outs.zip ->
native panel CP10K-log1p extraction -> scVI co-embed (CPU, seed 0, 60ep) with matched malignant census ref
-> top-K(15) spread for off-panel targets -> BEC/LEC (PECAM1 / LYVE1|PROX1|CCL21|FLT4) -> cell cap 0.
Xenium bundles (10x public CDN, CC BY 4.0): prostate Xenium_Prime_Human_Prostate_FFPE (12GB, FOLH1 on-panel),
lung Xenium_V1_humanLung_Cancer_FFPE (8GB, DLL3 transfer), ovary Xenium_Prime_Ovarian_Cancer_FFPE_XRrun (10GB
zarr, MSLN+MUC16 on-panel), skin Xenium_V1_hSkin_Melanoma_Base_FFPE (6GB, PMEL on-panel, EGFR/costim transfer).
Matched malignant refs (census 2025-11-08): prostate 881e0e6b, lung 576f193c, ovary b252b015, skin b617ee1b.
OVARY ZARR NOTE: the Prime FFPE ovary bundle ships Xenium Onboard Analysis 3.0 zarr (cell_feature_matrix.zarr.zip
+ cells.zarr.zip), NOT h5. Feature ids are Ensembl. Converted to standard 10x h5 via gene_panel.json Ensembl->symbol
map (5101 pairs); CSC csc/{data,indices,indptr} maps directly (indptr len=ncells+1); centroids from cells.zarr
cell_summary[:,0:2] (microns). All 15 targets on-panel after remap.

## BEC/LEC MARKER VERIFICATION (2026-07-13) — auditor-prompted
build_cancer_abm.py uses LEC = BEC ∩ (LYVE1|PROX1|CCL21|FLT4)+ — this is the CANONICAL ORGAN-build formula
(build_organ_beclec_method), a strict superset of the locked CRC-tumor 2-marker form (LYVE1|PROX1). Verified
per-panel marker presence: lung + skin carry only LYVE1/PROX1 -> 4-marker reduces EXACTLY to 2-marker (identical
to locked tumor form); prostate + ovary additionally carry FLT4 -> adds FLT4+ LEC (matches organ method, the
"identical to organs" instruction). Absent markers contribute zero via gv()->zeros, so no silent behavior change.
NOTE: skin Xenium panel lacks PECAM1 -> BEC marker-mask empty -> spatial_fallback path (grid-based BEC/LEC seeding,
same graceful fallback the organ builds use when <50 marker-BEC). Documented so skin BEC provenance is explicit.

## TECLISTAMAB IV PK VALIDATION (2026-07-13) — route-matched, dose-day-aligned
First clean route-matched PK gate. Model 0.72 mg/kg IV (50.4mg @ 70kg ref, MW 146kDa) vs PMC9345835 Fig3.
TIMING: clinical uses study-day frame (dose at day 1, rise 2500->22000 at day 1.1); model doses t=0.
Aligned clinical x by -1d. Post-alignment ratios: day0.1(peak) 0.68, day0.3 0.74, day1.1 0.90, day2.1 0.95,
day7 1.07 — distribution+elimination 0.90-1.07x (excellent). Peak 0.68x under.
MW NOTE: C(ug/mL)=A_pl/V_pl*mw/1e3 and dose_nmol=mg/mw*1e3 -> MW CANCELS from ug/mL curve (only affects nM
for binding/TMDD). MW not a source of overlay error.
PEAK GAP = BODY WEIGHT: peak C ~ mg/V_pl (MW-independent). Clinical peak ~22 ug/mL implies ~68mg -> ~95kg cohort.
Model uses 70kg fixed-physiology reference (0.72*70=50.4mg -> 16.3 peak). MajesTEC-1 IV dose-esc cohort likely
>70kg. NOT tuning weight to fit; 70kg is principled default. Honest: peak under because ref weight < cohort weight.

## HEME TMDD SINK — SURVIVAL-WEIGHTING FIX (2026-07-13)
BUG: the heme malignant-blast TMDD sink (coupled_percell_pd.py ~L237) internalized drug on the STATIC full
per-cell target copies (self.heme_R_sink), with no weighting by blast survival. So as the drug killed blasts
(depletion -> 70%), the sink kept drawing plasma as if 100% of the burden were still alive. At LOW dose (0.27
mg/kg teclistamab) this static full-burden sink over-drained plasma into a non-physiological terminal CLIFF:
C fell 382->99->13->2.7 ng/mL over ~0.8d (day 5.5-6.3) then flatlined ~2 ng/mL, vs clinical ~1500 ng/mL at day 7.
High dose (0.72) masked it (enough drug mass to survive the draw).
FIX: weight sink internalization by live-blast survival surv_h=exp(-heme_pd.kill_hazard) (heme_pd and heme_R_sink
share cell ordering, both from the same taa array in attach_heme_pd). Sink now shrinks as burden depletes -- the
physically correct TMDD behavior (dead blasts don't internalize). Applied to both copies (md5 1f174dfb82df).
Re-running both IV doses (0.72 + 0.27) to confirm cliff removed + dose-proportional disposition preserved.
NOTE: same sink is used by all 6 heme engagers -> full heme panel needs re-run through the fixed sink.

## TECLISTAMAB 3-AXIS PK+PD+CRS VALIDATION (2026-07-13)
PK: sink-fix confirmed — low-dose cliff removed (day7 1863 ng/mL vs broken ~2; clinical ~1500). Both dose
cohorts (0.72+0.27 mg/kg) thread clinical points, dose-linear (Cmax ratio 2.69 vs dose 2.67). Disposition 0.90-1.07x.
PD: model BCMA+ blast survival -> ~30% baseline; clinical serum IgA (PMC13109173, plasma-cell-killing PD marker)
-> ~20% baseline. Fast deep depletion matches. Target-cell killing VALIDATED against clinical Ig decline.
CRS CAVEAT (IL-6 OFF): model IL-6 peak 1887 pg/mL vs clinical severe-case 520 (PMC13116515 serum) vs cohort-median
anchor 21.5. Model ~3-4x hot even vs severe case. CAUSE (hypothesis): IL6_SCALE anchored to mosunetuzumab (CD20,
heavy first-dose CRS); heme compartment applies mosun-scale cytokine emission to teclistamab engaged T cells.
This is a cytokine-SCALE calibration issue, not a kill/PK mechanism bug. Flagged for cytokine recalibration.
overlay: teclistamab_IV_PKPD_2026-07-13.png (artifact beac9166 v2 e173ed84)

## TECLISTAMAB SC REGIMEN — DOSING TRUNCATION FIX (2026-07-13)
BUG: teclistamab SC regimen had only 19 doses, LAST at day119 (wk17). Curve accumulated to 45 ug/mL (wk18) then
decayed to ZERO by wk24 — no dosing after wk17. Clinical (PMC13222323) doses QW to wk40. FIX: extended regimen QW
105mg to day210 (wk30), now 32 doses. Re-running. NOTE steady-state overshoot to investigate: pre-fix accumulation
peaked ~45 ug/mL vs clinical ~27 (~1.7x); check whether maintained QW steady-state settles lower or if F_sc/dose
needs review.

## Kinetic binding unification (afternoon)
- kinetic_rhoden_percell.py: 6-species kinetic bivalent Rhoden binding + receptor turnover (KSYN=Ag0*kdeg, kdeg, kTMD). Explicit form, mass-conserved 1e-14, hits QSS exactly. (Reverted a matrix-BE solver detour that leaked; explicit form is canonical.)
- Rhoden attribution corrected: crosslink structure=Rhoden2016; turnover terms=user MATLAB scheme (Rhoden disclaims antigen kinetics).
- All 4 PK sinks (organ, tumor, heme, blood) route through the kinetic engine with turnover.
- multiarm_kinetic.py: general multi-arm engine (CD3/costim/TAA, valency 0-2, sum<=4 tetravalent). Mass-verified 1+1/2+2/1+1+1.
- Per-molecule co-engagement span threaded cfg->harness->PK->set_antigen->sink (was hardcoded 12.5).
- COSTIM SEPARATED: costim is its own kinetic binding species (own kon_cos/koff_cos/kdeg/B_cos), independent of the CD3 kill bridge. Separation verified: KO CD3 -> kill=0 costim survives; KO costim -> costim=0 kill survives.

## Full unification complete (all 5)
- (a) per-molecule co-engagement span threaded cfg->harness->PK->set_antigen->sink.
- (b) kDEG/kint per molecule (antigen_kint/antigen_kdeg + cfg override).
- (c) PD-side: costim SEPARATE kinetic species (own kon_cos/koff_cos/kdeg_cos/B_cos); Bcis_T = true CD3-CD3 avidity only (a_cd3>=2), NOT CD3+costim (audit fix). CD3 receptor turnover added to kinetic_synapse (RC relaxes to RC0 via KSYN=RC0*kdeg).
- (d) per-organ W proximity graph: ALREADY present — _build_neighborhoods() runs per organ at R_SYN=30um from real x,y (every organ ABM has coords).
- (e) IMMUNE MOTILITY: OrganPD.move_immune() — persistent random walk + chemotactic bias toward antigen, synapse graph rebuilt on a stride. Env: IMMUNE_MOTILITY=1, MOTILITY_STRIDE, MOTILITY_SPEED_UM_MIN, MOTILITY_CHEMOTAX. Verified: T cells 30um from targets -> 0 synapses -> migrate in -> 5510 synapses.

## Adopted reference unified_binding.py as canonical binding kernel (from binding-diagnostic lane)
- ONE Rhoden-kinetic solve -> PK sink (pk.intern_flux) + PD (pd.bK kill / pd.bC costim / pd.dwell CRS). PK and PD now provably one solve (not two), so they cannot double-count or desync.
- Numerics fix (dissolves my explicit-Euler stiffness/leak thrash): free receptors carried AS STATES in the linear generator M -> receptor census is a left-invariant -> backward-Euler conserves EXACTLY (no manufacture) + exact QSS. Batched np.linalg.solve, nsub set by slow scales (1-4), NOT binding rmax.
- Verified in totalvi-venv: conservation 1e-14, QSS 1e-11, 0.5s/300k-cell bivalent step, format reductions 1+1/2+1/2+2/1+1+1/1+1+2, CD3/costim separation PASS.
- Reuses model constants verbatim (NM_PER_COPY, geo_ageff, spans, K_HIT).
- Superseded: my kinetic_rhoden_percell explicit core + separate kinetic_synapse PD solve. Kept: per-molecule span/kdeg/kint harness threading, immune motility (move_immune), per-organ W.

## CORRECTION to the line above (merge decision, coordinated with binding-diagnostic lane)
- Do NOT wholesale-replace the production engine. Ship the WIRED production engine (motility, real params, per-organ W, explicit Bcis_T CD3-CD3 avidity, running end-to-end).
- Borrow ONLY the receptors-as-states backward-Euler inner solver from unified_binding.py, and ONLY for the bivalent/tetravalent PK path (the nsub=20000 stiff case). Monovalent keeps the explicit form (fast+conserving; = the whole heme+teclistamab validation set).
- Keep explicit Bcis_T (more faithful than reference koff-retention cis). Port the reachable-TAA reversible cap convention with the solver.
- unified_binding.py stays as the verified solver SOURCE + standalone reference, not the running kernel.

## MEASURED: bivalent PK sink does NOT manufacture at panel densities (run, not reasoned)
- Tested live rhoden_bivalent_step across CD20(2.5e5), high(1e6), EGFR-ceiling(2.3e6) copies at 50nM drug.
- nsub picked: 271 / 1073 / 2465 — ALL well under the 20,000 cap. Census conserved 0.72x (real turnover, not manufacture).
- The nsub fix (kon*AgEFF included in substep count) already protects every panel molecule. The other lane's 150-850x manufacture needed an EXTREME density regime not reached by any real molecule.
- CONCLUSION: biexact's value here is SPEED, not correctness. The explicit PK bivalent path is safe for the submission panel. Kill/PD path (kinetic_synapse 2x2 matrix-exp) was always stable.
- REAL bottleneck for 9PM throughput = organ-loop parallelism: self.pool is a ThreadPoolExecutor, transport step is GIL-bound numpy/splu, so 11 organs run ~serially (~2 cores of 80). Spleen transport 157ms/step x 1001 steps dominates.

## PARAMETERIZATION AUDIT + FIXES (2026-07-13, prompted by "did we use all the sub-agent params")
Sub-agent fan-out gathered eng_params_normalized.json (22 molecules, measured kon/koff both arms, KD, kint, kdeg, mw, fFcRn, regimen, il6_obs). AUDIT found the harness NEVER read it -> every molecule ran on generic kon=1e5 fallback.

FIXES:
1. WIRED eng_params_normalized.json into run_tce_pd_reval.py: merge step at cfg=ENG[name] enriches cfg with
   measured kon_TAA/koff_TAA/kon_CD3/koff_CD3/kdeg/kint. Variant names (teclistamab_iv->teclistamab) map to base.
   e.g. glofitamab kon_TAA now 430000 (was 1e5 default, 4.3x).
2. CD3 KD now PER-MOLECULE: the 3 attach_pd/heme/blood sites read cfg.get('KD_CD3_nM',40.0) (was hardcoded 40).
3. glofitamab CD3 KD corrected 40 -> 4.5 nM (FDA BLA 761309 SPR Study 1108724: KD=4-5nM CD3e; CD20 EC50 4.3nM cell-based). koff_CD3 recomputed = kon*4.5e-9.
4. glofitamab REGIMEN corrected in regimen_schedules_final.json: was [2.5@d0, 10@d7, 100@d14...] (WRONG: 100mg, wrong days).
   Now per established Columvi USPI (glofitamab-gxbm labeling): 2.5mg(C1D8/day7), 10mg(C1D15/day14),
   30mg(C2D1/day21) Q21d. Target dose is 30mg not 100mg. PROVENANCE: dose/day values from approved-label
   knowledge; web_search body not persisted in transcript -> RE-VERIFY vs printed Columvi USPI/Drugs@FDA
   before final sign-off (values believed correct, sourcing not archived).
   NOTE: obinutuzumab pretreatment (1000mg C1D1, 7d prior) pre-depletes B cells — model starts full CD20 density,
   so first-dose CRS/depletion is an UPPER BOUND. Flagged for refinement.

ORPHAN fields (correctly not consumed): name/tgt/sources/KD_CD3_source/measured_vs_derived (provenance);
   validation_data (plot-time overlays); ksyn_frac_perday (engine computes KSYN=Ag0*kdeg internally).
KNOWN GAP: n_arm_CD3=2 for cinrebafusp_alfa (2:2 HER2) — kinetic_synapse assumes monovalent CD3. 1/22, edge case.

## GLOFITAMAB = MECHANISM TEST, NOT CLEAN PD VALIDATION (obinutuzumab confound)
FDA Columvi: obinutuzumab 1000mg C1D1 pre-depletes B cells 7d BEFORE glofitamab starts (C1D8). Clinically,
glofitamab engages an ALREADY-DEPLETED CD20 compartment. Our model starts full B-cell density, so glofitamab
depletion/CRS is an UPPER BOUND and should NOT be graded against clinical PD magnitude.
- USE glofitamab for: bivalent-avidity mechanism (does Bdbl form via rhoden_samecell_bivalent_step), PK sanity.
- USE for clean PD magnitude: teclistamab/mosunetuzumab/epcoritamab (monovalent, no anti-CD20 pretreatment).
- Most CD20 TCEs (mosun/epcor/odron) also use obinutuzumab or rituximab pretreatment in some regimens — flag per-molecule.
PROCESS FIX: runs now tracked by EXACT PID (RUN_PIDS.json); NEVER pattern-kill (pattern matched the launcher's own
shell because the command string contains the molecule name -> killed our own runs repeatedly = every "died no json").

## THROUGHPUT FINDING: runs must be SEQUENTIAL locally, not naive-parallel (2026-07-13)
Symptom: two heavy runs concurrent (teclistamab+glofitamab, each 4.6GB/21threads) -> EACH ~9x slower
(teclistamab reached step 20 at 1124s vs 127s when run ALONE). NOT a hang, NOT a param bug: it's
memory-bandwidth + GIL thrashing from oversubscription (2 procs x internal thread pools on the transport splu).
Isolated binding step is 118ms (measured kon=1.28e6 fine, slightly faster than 1e5 default).
DECISION: local runs SEQUENTIAL (one at a time, generous OMP threads). For the full 22+55 panel, this is why
the user's Modal cloud fan-out matters — true process isolation across containers, no shared-memory contention.
Merged params confirmed clean floats (no type corruption): teclistamab kon_TAA=1.28e6, kon_CD3=1.16e5, KD_CD3=28.03.

## ROOT-CAUSE FIX: heme/blood sinks still called OLD rhoden_bivalent_step -> nsub-pin stall (2026-07-13)
py-spy dump of the stuck run (PID 1767409) showed MainThread frozen in rhoden_bivalent_step (kinetic_rhoden_percell.py:75)
called from simulate_pd:284 (the BLOOD sink). Root cause: I fixed the ORGAN PK sink to rhoden_samecell_bivalent_step
but MISSED the heme sink (coupled_percell_pd.py:268) and blood sink (:284) — both still called the OLD 6-species
rhoden_bivalent_step with Ag2=0. With teclistamab's MEASURED kon=1.28e6 (vs old 1e5 default), the old solver's nsub
formula (sees kon*AgEFF) pinned nsub=20,000 -> each heme/blood step took minutes -> step-0 stall. The generic 1e5
had kept it just barely tolerable, which is why it only broke AFTER the param merge wired in measured kon.
FIX: heme+blood sinks -> rhoden_samecell_bivalent_step (same fast BE solver as organ sink). Both copies, parity MATCH.
RESULT: step 20 now at 38s (was 1124s = 30x faster); CPU 387% (multi-core, no GIL step-0 thrash). Physics identical
(spleen 0.93/bone 0.71/lung 0.00). Full 251-step run projected ~8 min.
LESSON: when wiring a solver fix, grep ALL call sites (organ + heme + blood + tumor), not just the first one found.


## 08:53 — Validation panel composition + readiness (agent lane)

**Heme validation panel (clean PK AND PD, runnable once FIX-1-UPGRADE lands):**
- teclistamab, elranatamab, talquetamab, mosunetuzumab, epcoritamab — INVOKABLE (ENG entries + regimens present)
- blinatumomab — has params (eng_params_normalized: CD19 KD1.49/CD3 260/kint0.15/no-Fc BiTE) but MISSING ENG entry;
  needs continuous-IV-infusion regimen. COORDINATION: other lane to add ENG entry in same pass as FIX-1-UPGRADE.

**Solid TCE data reality (deferred until after heme, per user):**
- cibisatamab (CEA) — 1 PK curve, 0 PD; runs on CRC tumor (built) — needs harness per-target switch
- tarlatamab (DLL3) — 1 PK, 3 PD; needs lung tumor wired
- tebentafusp (gp100) — 0 PK, 1 PD (ctDNA); needs skin tumor wired
- catumaxomab (EpCAM) — ascites/intraperitoneal, EXCLUDED (local not systemic)
- glofitamab, odronextamab — NO digitized data; mechanism demos only

**Tumor builds: 4/5 done** (lung 162k, ovary 407k, prostate 193k, skin 107k — scVI overlay + ECM genes present).
Pancreas NOT built (reference blocked). NONE run-ready: BEC/LEC masks not in npz + not copied to harness agents/ +
per-target compartment switch not implemented (harness always uses CRC). This is the deferred solid-wiring work.

**FIX status (other lane owns):** FIX-1-UPGRADE (measured kon/koff+KD_norm into attach_pd/heme/blood) — IN FLIGHT,
KD_norm setdefault present (L93) but override into attach_pd not yet wired (L108 still generic KIN_PARAMS).
FIX-3 (alias resolver + fail-loud density assert) — not yet present. Agent HOLDING harness edits to avoid collision.

**Deck-ready figures (agent lane):** teclistamab PK overlay AFE 1.29x (corrected legend, no unverified n=15);
glofitamab model PK (mechanism, framed); 5-panel spatial spleen/bone/liver (full cells, cell-type panel, tissue underlay).
**PD numbers provisional** until FIX-1-UPGRADE re-run (PK measured-rate-driven; PD kill/CRS still generic kon).


## 09:03 — Elranatamab matched-validation run (PK+IL6, MagnetisMM) + engine sync

**FIX-1-UPGRADE + FIX-3 landed (other lane, engine/) — SYNCED to handoff/ + verified.**
- engine↔handoff parity restored for run_tce_pd_reval.py + coupled_percell_pk.py (backups in handoff/_presync_bak/).
- FIX-1-UPGRADE verified: teclistamab PD synapse koff_TAA now 2.315e-4/s (MEASURED) vs 1.8e-5 (generic kon*KD) = 13x. PK=PD identical.
- FIX-3 verified by other lane: CD20->MS4A1_copies alias rescue; heme molecules resolve unchanged; true miss warns loud.

**Elranatamab chosen as PRIMARY matched-validation molecule** — the only one with PK AND IL-6 from the SAME paper/cohort
(PMC10579053 Nat Med 2023, MagnetisMM-1: PK Fig2 + IL-6 Fig6, 12/32/76 mg SC step-up).

**Elranatamab param accuracy corrections (ENG line, both copies parity):**
- mw 146.0 -> 148.5 (Elrexfio USPI label, web-verified)
- fFcRn 0.89 -> 0.70 (IgG2 with intact Fc; affects k_cat clearance)
- il6_obs 285 -> 340 (matched Fig6 priming-only peak; 230 = priming+premed)
- kinetics via normalized merge (Yoneyama 2022 Table 2): kon_TAA=9.95e5, koff_TAA=3.98e-5, KD_TAA=0.04nM,
  KD_CD3=17nM (AMG420 platform, flagged), kint=2.0 (Lee 2016), kdeg=2.05. koff DERIVED from kon*KD (flagged).
- regimen REG: 12/32/76 mg SC step-up, 19 doses (correct).

**Panel classification correction:** tarlatamab is DLL3/SCLC = SOLID (needs lung wiring), NOT heme.
Runnable-now heme panel: teclistamab, elranatamab, talquetamab, mosunetuzumab, epcoritamab, blinatumomab(needs ENG entry).

Run: PID 1819001, KWS=model/rundir, PD_DT=0.04, SNAP 1/7, ~19min projected. Step 20@38s clean, 396% CPU.


## 11:00 — Mechanistic myeloid IL-6 engine synced + blood-term bug fixed (agent lane)

**Other lane replaced fitted IL6_SCALE with mechanistic per-cell myeloid model** (4 files, engine/):
- NEW myeloid_il6.py (MyeloidIL6 per-cell emitter + PlasmaIL6 clearance ODE)
- wholebody_pd.py / coupled_percell_pd.py / run_tce_pd_reval.py rewired: il6_pgml now = il6_plasma_pgml (physical pg/mL)
- Literature params (none fitted): 0.0196 pg/hr/cell (156 molec/s, PMID 37533643), 150min time-to-max, 0.20/hr clearance (PMID 31268236), R_SYN 30um contact. IL-6 = myeloid-derived (Giavridis PMID 29808005, Norelli PMID 29808007).
- OPEN GATE: myeloid_count_scale=1.0 (sampled) -> absolute pg/mL low until organ myeloid cellularity sourced; RELATIVE per-drug ordering correct.

**SYNCED engine->handoff** (4 files: myeloid_il6.py NEW + 3 diverged); backups _presync_bak/ (ts=105545); all 6 parity-SAME. My elranatamab param fix (mw148.5/fFcRn0.70/il6 340) + FIX-1-UPGRADE survived in the new engine copy.

**BUG FOUND + FIXED (coupled_percell_pd.py L269, blood myeloid IL-6 sum):**
- crash: `float(self.blood_count_scale or 1.0)` -> ValueError "truth value of array ambiguous" (blood_count_scale is a per-CELL array, L169)
- root: blood term wrongly used the per-cell antigen-sink array as a scalar; also wrong scale (antigen-derived, per-drug)
- fix: use `_b.myeloid_count_scale` (blood_pd's SCALAR tissue-property scale) — matches the organ (L265) + heme (L272) pattern
- both copies parity-SAME; compiles; elranatamab re-run PID 1832907 now steps past the crash point clean.
- FLAG for other lane: verify this matches their intent (blood myeloid should use tissue-property scalar, not per-cell blood_count_scale).

Re-run in progress: elranatamab on mechanistic-IL6 engine (all prior IL-6 numbers superseded per other lane).
