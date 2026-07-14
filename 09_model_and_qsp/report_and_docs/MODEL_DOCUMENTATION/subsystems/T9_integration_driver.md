---
title: "T9 — Integration driver: run assembly, PK/PD cadence, organ loop, compartment wiring"
subsystem: T9 (integration / dataflow spine)
model: costim_engager_counterscreen (whole-body per-cell PBPK-PD TCE counterscreen)
live_source_files:
  - engine/coupled_percell_pd.py   (395 lines, read in full this task)
  - engine/run_tce_pd_reval.py     (274 lines, read in full this task)
note_on_paths: >
  engine/*.py and rundir/handoff/*.py are byte-identical copies (verified this task for
  coupled_percell_pd, run_tce_pd_reval, wholebody_pd, wholebody_percell, coupled_percell_pk,
  myeloid_il6, kinetic_rhoden_percell, pd_model_config). The driver puts {KWS}/handoff at the
  front of sys.path, so the handoff copy is the one actually imported — the line numbers below
  are valid for both. The JSON data files are NOT all duplicated safely: see the duplicated-regimen-file warning in §3.3.
date: 2026-07-13
generated_by: workflow-subagent T9
---

up:: [[00_INDEX]]
tags:: #model-doc #integration #dataflow
sibling docs:: [[T1_shah_betts_pbpk_backbone]] · [[T2_whole_body_per_cell_pk]] · [[T4_multi_arm_format_geometry]] · [[T5_kinetic_immune_synapse]] · [[T6_per_cell_pd_killing]] · [[T7_costim_signaling_activation_induced_receptor_density]] · [[T8_mechanistic_crs_il_6]]

> [!abstract]+ At a glance
> **T9 is the spine.** It is the only code that knows how a *run* is assembled: which molecule, which target, which
> compartments exist, what the clock is, in what order the physics fires, and what gets written out. Every other
> T-subsystem is a callee of the loop documented here.
> - **Files owned:** `engine/coupled_percell_pd.py` (build + step loop), `engine/run_tce_pd_reval.py` (molecule registry + CLI).
> - **Live import graph (verified this task):** `run_tce_pd_reval` → `coupled_percell_pd`, `pd_model_config`, `qsp_costim_window_v2`; `coupled_percell_pd` → `coupled_percell_pk`, `kinetic_rhoden_percell`, `myeloid_il6`, `wholebody_percell`, `wholebody_pd`; `wholebody_pd` → `costim_induction`, `kinetic_synapse`, `multiarm_binding`, `myeloid_il6`. **`signaling_dynamics` is NOT in this graph**: `wholebody_pd.py:209` imports it, but only inside the `costim_arm is not None` branch, and **the module does not exist anywhere in the repo** (RUN-verified: `importlib.util.find_spec('signaling_dynamics')` is `None` on the driver's own `sys.path`). Every clinical molecule runs with `costim_arm=None`, so the branch never fires. It is a **crash landmine**, not a live import — see §5.10.
> - **NOT in the execution path (do not cite):** `cytokine_pbpk.py`, `il6_pbpk.py`, `unified_binding.py`, `multiarm_kinetic.py`, `biexact_solver.py`, `rna_to_receptor.py`, `convert_copies_ALL.py`, `calib_kdeath.py`.
> - **Every parameter below is tagged** `[MEASURED: …]` / `[DERIVED: …]` / `[FITTED: …]` / `[ASSUMED: …]` / `[UNSOURCED — TBD]`, plus `[UNVERIFIED CITATION]` wherever an in-code comment names a source I could not confirm from the code itself.

> [!danger]+ Three traps a reader must know before quoting any number out of this engine
> 1. **The simulation horizon can silently amputate the dose that causes the effect.** `TSIM_DAYS=7` gives
>    mosunetuzumab **only its 1 mg priming dose** — its 60 mg dose is at **t = 14 d** in the live schedule
>    (`handoff/regimen_schedules_final.json`, loaded at `run_tce_pd_reval.py:44`), while teclistamab has already
>    received its full **105 mg** by t = 7 d. Comparing molecules at 7 days compares **different rungs of different
>    step-up ladders**. See §2 EQ-1 / §3.3 / §5.2.
> 2. **Only IL-6 has physical units.** `il6_plasma_pgml` is pg/mL. `cyto_sys_rate` / `cyto_sys_cum` (IL6/IFN/TNF/IL2)
>    are **raw engagement-sum units** — the pg/mL converter (`wholebody_pd.cytokine_to_pgml`) is **never called** in
>    this driver. Do not present them as concentrations. See §5.5.
> 3. **The IL-6 clinical anchors in the registry are display-only, and most historical values were contaminated.**
>    The only surviving anchors in live code are **mosunetuzumab 152** and **teclistamab 21** (both *population means*);
>    everything else is `il6_obs=None` (`run_tce_pd_reval.py:53–89`). The retired values (570 / 340 / 230 / 191 / 288 /
>    366.88) are documented in §3.8 **as contamination, not as data** — 191 was a *page number*; 288 is the *highest
>    individual patient*; elranatamab has **no** clinical IL-6 value in existence.

---

## 1. PURPOSE & DATAFLOW POSITION

### 1.1 What T9 does

T9 answers one question: **given a molecule name on the command line, what exactly gets built, in what order does it
step, and what comes out?** It owns:

1. **The molecule registry** (`ENG`, `run_tce_pd_reval.py:53–89`) — target gene symbol, KD, MW, FcRn-protected fraction,
   TAA valency, route, and the (display-only) clinical IL-6 anchor for each of 9 engager entries.
2. **Per-molecule kinetics merge** (`run_tce_pd_reval.py:106–124`) — the measured kon/koff/kint/kdeg from
   `handoff/eng_params_normalized.json` are merged over the registry so PK and PD bind on **one identical scheme**.
3. **Compartment assembly** — 11 healthy organs (+ solid tumor as a 12th spatial ABM organ), a plasma-driven **heme
   malignancy** compartment, and an always-present **normal blood** compartment (`run_tce_pd_reval.py:142–168`;
   `coupled_percell_pd.py:17–197`).
4. **The step loop** (`coupled_percell_pd.py:226–386`) — dosing → plasma → per-organ transport → per-cell PD →
   myeloid IL-6 → TMDD sinks → plasma/lymph mass balance → recording.
5. **The PK/PD cadence split** — transport every step (`dt = 0.02 d`), PD every 3rd step (`pd_every = 3` → 0.06 d).
6. **The output contract** (`coupled_percell_pd.py:387–394`, `run_tce_pd_reval.py:225–255`).

### 1.2 Where it sits in the life of the molecule

```
  mg on a clinical schedule  (regimen_schedules_final.json)
        │  EQ-1  mg → nmol
        ▼
  A_pl  (plasma amount, nmol)  ──EQ-2 IV zero-order infusion / EQ-4 SC depot──►  C_pl = A_pl/V_pl   [EQ-3]
        │
        ├── per organ o ∈ ORGANS(+tumor):  PS_ex, PS_ly  [EQ-5]  ─►  TissueGraph.step   [EQ-6 → T2]
        │        vascular QSS → BEC extravasation → per-cell interstitial diffusion+binding → LEC drain
        │        returns:  organ TMDD sink, lymph drain, J_extrav, bound S, bridged D
        │                        │
        │                        └── g.C  (per-cell free drug, nM)  ──►  OrganPD.step   [EQ-7 → T5/T6/T7]
        │                                     per-cell trimer → kill hazard → engaged-T fraction p_eng
        │                                                                            │
        │                                                                            ▼
        │                                                          MyeloidIL6.step (per-organ)  [→ T8]
        │                                                          il6_prod_pg_hr (unscaled, sampled cells)
        │                                                                            │
        │                        Σ_organs  count_scale_o × prod_o  [EQ-8]  ──────────┘
        │                                     │
        │                                     ▼
        │                        PlasmaIL6.step   [EQ-9]  →  il6_plasma_pgml  (the ONE physical PD readout)
        │
        ├── heme blasts (plasma-driven): PD kill + kinetic-Rhoden TMDD sink  [EQ-11..EQ-13]
        ├── normal blood (plasma-driven): PD kill (tox) + TMDD sink; myeloid GATED OFF  [EQ-10, EQ-14]
        │
        ▼
  dA_pl/dt = infusion + F_sc·J_sc + k_lymph_return·A_ly − k_cat·A_pl − Σ J_extrav − heme_sink − blood_sink   [EQ-16]
  dA_ly/dt = Σ drain − k_lymph_return·A_ly                                                                  [EQ-17]
```

### 1.3 Build phase (once per molecule) — what is constructed and from what

| Step | Code | What it builds |
|---|---|---|
| PBPK arrays | `run_tce_pd_reval.py:8, 49–51` | `Q, L, sigV, Vis, Vv` for **all 15 PBPK tissues** (incl. `tumor`, idx 14) from `qsp_costim_window_v2._PBPKArrays` → T1 |
| target kinetics | `run_tce_pd_reval.py:11–33` | `kint`, `kdeg`, `kon`, `koff` per target/molecule |
| engine choice | `run_tce_pd_reval.py:36–43` | `PD_KINETICS` from `pd_model_config` (canonical = **kinetic**), env A/B override `PD_ENGINE_OVERRIDE` |
| organ set | `run_tce_pd_reval.py:48, 126–129` | 11 organs; **+`tumor`** iff `cancer_type=="solid"` |
| transport core | `run_tce_pd_reval.py:142–145` → `coupled_percell_pk.CoupledPerCellPK.__init__` | one `TissueGraph` per organ: KDTree graph Laplacian, BEC/LEC masks, per-cell antigen state, `count_scale` |
| organ PD | `run_tce_pd_reval.py:152–154` → `coupled_percell_pd.attach_pd:17–51` | one `OrganPD` per organ (CD3 arm, TAA arm, optional costim arm, myeloid emitters) |
| heme PD (+sink) | `run_tce_pd_reval.py:160–163` → `attach_heme_pd:53–137` | malignancy ABM chosen by target via `heme_compartment_routing.json` |
| blood PD (+sink) | `run_tce_pd_reval.py:166–168` → `attach_blood_pd:139–197` | Tabula-Sapiens blood ABM; **blood myeloid gated OFF** (`:188`) |

**Cancer-type switch is exclusive** (`run_tce_pd_reval.py:125–129, 155–164`): heme malignancy **XOR** solid tumor,
never both. The 11 healthy organs + normal blood are **always** present (they are the PK sinks, the CRS source and the
on-target tox readout).

### 1.4 The molecule/ENG registry and its clinical anchors

`ENG` (`run_tce_pd_reval.py:53–89`) is 9 entries (6 clinical molecules + 3 route-matched PK-validation variants).
Fields: `tgt` (HGNC symbol → per-cell density column), `KD` (nM, the registry KD used for the **PK** TMDD arm),
`mw` (kDa), `fFcRn` (FcRn-protected fraction → sets catabolic clearance, EQ-16), `narm` (TAA valency → bivalent
avidity in the Rhoden binder), `il6_obs` (clinical IL-6 anchor, **display/scoring only** — the header states
explicitly it "never enters the mechanism", `:55–56`), `route`, and for SC `F_sc`/`ka`.

> **Anchor integrity (the code's own audit block, `run_tce_pd_reval.py:58–76`).** Every anchor must be the *same
> statistic*. The surviving anchors are **population-mean IL-6 peak Cmax**: teclistamab **21** (verbatim quote in code:
> *"The mean IL-6 peak concentration (Cmax) was 21 pg/mL"*, attributed to **PMID 38831634** — `[UNVERIFIED CITATION]`,
> I read the attribution in code, I did not open the paper) and mosunetuzumab **152** (population mean; **no PMID is
> given in the code** → `[UNSOURCED — TBD]` for the citation, the *value* is the one the code carries). Everything else
> is `None`. **Elranatamab has no clinical IL-6 value in existence** (`:71`). See §3.8 for the retired values and why
> each is unusable.

### 1.5 Step phase — exact order of operations inside one `k`-step

`coupled_percell_pd.simulate_pd`, loop at `:226`:

| # | Line(s) | Operation | Cadence |
|---|---|---|---|
| 0 | `:227` | `t = k·dt` | every step |
| 1 | `:229–233` | optional tagged progress print | `PROGRESS_EVERY` |
| 2 | `:234–248` | consume all scheduled doses with `t_dose ≤ t + 1e-9` | every step |
| 3 | `:249–256` | deliver active IV infusions into `A_pl` (mass-exact clamp) | every step |
| 4 | `:257–259` | `C_pl`; SC depot flux `J_sc`; constant-rate infusion `infn` | every step |
| 5 | `:261–269` | **organ loop** (ThreadPool over organs): `TissueGraph.step` → `drain`, `J_extrav`, `S` | every step |
| 6 | `:271–273` | **PD step** per organ with that organ's per-cell `g.C` and `dt·pd_every` | every `pd_every`(=3) |
| 7 | `:285–297` | myeloid IL-6 census-weighted sum → `PlasmaIL6.step` | every `pd_every` |
| 8 | `:299–304` | immune motility (opt-in, default OFF) | `k % MOTILITY_STRIDE` **within** the PD block |
| 9 | `:305–325` | heme PD step + heme kinetic-Rhoden TMDD sink | every `pd_every` |
| 10 | `:326–338` | blood PD step + blood kinetic-Rhoden TMDD sink | every `pd_every` |
| 11 | `:339–352` | spatial snapshot (x,y,C,bound,R,labs,surv,is_target) | at `SNAP_DAYS` |
| 12 | `:353–356` | well-mixed compartments (`self.wm`) — **EMPTY in this driver** (`wellmixed=None`) | every step |
| 13 | `:357–359` | plasma + lymph mass balance, non-negativity clamp | every step |
| 14 | `:360–386` | recording (PK, kill, cytokines, IL-6, per-organ unscaled IL-6 production) | `nstep//rec` |

**Zero-order hold on the sinks.** `_heme_sink_nmol_day` / `_blood_sink_nmol_day` are computed **only** inside the
PD block (every 3rd step) but are subtracted from plasma **every** step (`:357`, via `getattr(...,0.0)`). Between PD
updates the sink is held constant. Before the first PD step the `getattr` default supplies 0.0 — no NameError, no
spurious sink.

---

## 2. GOVERNING EQUATIONS

Notation: `A_pl`, `A_ly`, `A_sc` = amounts (nmol) in plasma / lymph / SC depot; `C_pl` = plasma conc (nM);
`mw` = molecular weight (kDa); `dt` = transport step (day); `Vis[o]`, `Vv[o]` = interstitial / vascular volume (L);
`Q[o]` = organ plasma flow (L/day); `L[o]` = organ lymph flow (L/day); `AVO = 6.02214076e23`.

---

### EQ-1 — Dose schedule → plasma amount (mg → nmol) (`coupled_percell_pd.py:234–248`)

```
while sched[si].t ≤ t + 1e-9:
    n_dose[nmol] = mg / mw[kDa] · 1e3
    route == "SC"                    →  A_sc  += n_dose                       (:237)
    route == "IV" and iv_inf_h > 0   →  push infusion [t, t+dur_d, n_dose/dur_d, n_dose]   (:245)
    route == "IV" and iv_inf_h == 0  →  A_pl  += n_dose                       (:247)
```
- **Biology.** A clinical mg dose becomes molar drug in the vascular space. Dimensionally: 1 mg of a 146 kDa IgG =
  1/146 µmol = **6.849 nmol** (RUN-verified in Python this task); into `V_pl = 3.1 L` that is **2.209 nM** — the whole
  system's initial forcing.
- **Mechanistic rationale.** The schedule is a *list of (day, mg)* read verbatim from
  `handoff/regimen_schedules_final.json` (`run_tce_pd_reval.py:44, 169`), i.e. the **real step-up regimens**, not a
  single bolus. Rejected alternative: dosing at a fixed nominal concentration — that would erase the step-up ladder,
  which is precisely the structure that determines when CRS occurs.
- **⚠ THE HORIZON TRAP.** The `while` at `:234` only ever admits doses with `t_dose ≤ tsim`. Live mosunetuzumab
  schedule is `[(0,1),(7,2),(14,60),(21,60),(42,30),(63,30),(84,30),(105,30)]` mg. Therefore:
  - `TSIM_DAYS=7` → mosunetuzumab has received **1 mg** (its 2 mg dose lands on the *final* step, t = 7.0, with no time
    to act; its **60 mg** dose at t = 14 d never happens);
  - the same 7-day window gives teclistamab `[(0,4.2),(3,21),(7,105),…]` its **full 105 mg**.
  The code's own guidance (`run_tce_pd_reval.py:93–96`) is `TSIM_DAYS=24` as the safe truncation for an IL-6 peak run,
  with the explicit warning that PK Cmax/AUC from a truncated run are **not** steady-state values. Default `TSIM`
  (`:90`) is mosunetuzumab **49 d**, teclistamab **200 d**, most others 24 d.
- **Units.** `mg` [mg]; `mw` [kDa]; `n_dose` [nmol]; `t`, `dur_d` [day].

---

### EQ-2 — Mass-exact IV zero-order infusion (`coupled_percell_pd.py:243–256`)

```
dur_d      = iv_inf_h / 24                                            (:243)
tot_nmol   = mg / mw · 1e3                                            (:244)
rate       = tot_nmol / dur_d                                         (:245)   [nmol/day]
per step, for each active infusion with t0 ≤ t < t1 and remaining > 0:
    step_dose = min(rate · dt, remaining)                             (:253)
    iv_in    += step_dose / dt ;   remaining -= step_dose             (:254–255)
A_pl += iv_in · dt                                                    (:256)
```
- **Biology.** A clinical TCE is a ~2 h IV infusion, not a bolus. Infusing rather than bolusing changes the **peak**
  `C_pl` the first cohort of engaged T cells ever sees — and CRS is a peak-driven event.
- **Mechanistic rationale.** The `min(rate·dt, remaining)` clamp makes the delivered mass **exactly** `tot_nmol`
  regardless of whether `dur_d` is an integer multiple of `dt`. The in-code comment records that the unclamped version
  over-delivered by ~1.2× (`:240–242`). Rejected alternative: bolus (`iv_inf_h=0`, still reachable at `:247`).
- **Units.** `iv_inf_h` [h]; `rate` [nmol/day]; `iv_in` [nmol/day]; `A_pl` [nmol].

---

### EQ-3 — Plasma concentration (`coupled_percell_pd.py:257`)

```
C_pl = A_pl / V_pl                                   [nM] = [nmol] / [L]
```
- **Biology.** The systemic free-drug concentration every organ's vascular compartment sees, and the concentration the
  plasma-driven compartments (heme blasts, normal blood cells) see **directly** (`:307`, `:328`).
- **Rationale.** One well-mixed central compartment; spatial structure lives *inside* organs, not in plasma.
- Reported as µg/mL at `:361`: `C[µg/mL] = A_pl/V_pl · mw/1e3` (1 nM of 146 kDa = 0.146 µg/mL; RUN-verified).

---

### EQ-4 — SC depot: first-order absorption with bioavailability (`coupled_percell_pd.py:258`, `:357`)

```
J_sc  = ka_sc · A_sc                    [nmol/day]           (:258)
A_sc ← A_sc − dt · J_sc                                       (:258)
plasma gains  F_sc · J_sc                                     (:357)
```
- **Biology.** SC-administered TCEs (epcoritamab, teclistamab, elranatamab, talquetamab) are absorbed slowly from the
  injection site; a fraction `(1 − F_sc)` is destroyed presystemically (lymphatic/interstitial catabolism).
- **Mechanistic rationale.** The depot **loses the full `J_sc`** but plasma **gains only `F_sc·J_sc`** — bioavailability
  is implemented as a *loss*, not as a dose scaling. This is the correct structure: it preserves the absorption
  time-course while removing mass. Rejected alternative: `A_sc += F_sc·dose` at input, which would make the lost
  fraction disappear instantaneously rather than over the absorption profile.
- **Units.** `ka_sc` [1/day]; `F_sc` [—].

---

### EQ-5 — Per-organ permeability-surface products (`coupled_percell_pd.py:263`)

```
PS_ex[o] = k_dist · L[o] · (1 − sigV[o])          [L/day]   vascular → interstitium (at BEC)
PS_ly[o] = k_dist · L[o] · (1 − sigL)             [L/day]   interstitium → lymph (at LEC)
```
- **Biology.** Two-pore-style convective transport: an IgG crosses the vascular wall at a rate set by the organ's
  lymph flow `L[o]` and its vascular reflection coefficient `sigV[o]`, and leaves with lymph against the lymphatic
  reflection coefficient `sigL`.
- **Mechanistic rationale.** `k_dist` (=3.0, inherited default) is a lumped distribution multiplier on the
  lymph-flow-driven convective term. **This is the Shah–Betts-style formulation the PK core claims**
  (`coupled_percell_pk.py:1–5`) — `[UNVERIFIED CITATION]`: the docstring names "Shah-Betts fixed physiology" but the
  code contains **no PMID**; I did not verify the parameter set against the paper.
- **Units.** `k_dist` [—]; `L[o]` [L/day]; `sigV`, `sigL` [—].

---

### EQ-6 — Organ transport step (interface) (`coupled_percell_pd.py:261–269`)

```
sink, drain, J_extrav, S, D = g.step(C_pl, Q[o], PS_ex, PS_ly, KD, n_arm, kint)      (:264)
tot_extrav = Σ_o J_extrav ;  tot_drain = Σ_o drain                                    (:268–269)
```
- **Biology.** Each organ is a **per-cell spatial graph**: drug enters at BEC cells, diffuses through the interstitium
  on a KDTree graph Laplacian, binds/internalises on **each individual cell's own receptor number**, and drains at LEC
  cells. The returned `S` (bound copies/cell) and `D` (bivalently bridged) are what the snapshots and the TMDD sink read.
- **Mechanistic rationale for T9's part.** T9 **parallelises** the organ loop over a `ThreadPoolExecutor`
  (`coupled_percell_pk.py:106`, `min(len(organs), 9)` workers; used at `coupled_percell_pd.py:266`). The organs are
  **independent within a step** — they couple only through the scalar `C_pl` (read at the top of the step) and through
  the summed fluxes (written after the map). This is what makes the parallel map exact rather than an approximation.
- **Internals are owned by [[T2_whole_body_per_cell_pk]]** (`wholebody_percell.py:154–209`): vascular QSS
  `Cvasc = (Q·C_pl + PS_ex·C_bec)/(Q + PS_ex)` (`:168`), graded LEC drainage, per-cell Rhoden TMDD.
- **Units.** `sink`, `drain`, `J_extrav` [nmol/day]; `S`, `D` [copies/cell].

---

### EQ-7 — PD cadence: the multi-rate split (`coupled_percell_pd.py:271–273`)

```
if k % pd_every == 0:
    for o in organs:  pd[o].step(graphs[o].C, dt · pd_every, k_death=k_death)
```
- **Biology.** Transport (diffusion, extravasation, binding) is *fast*; the immune-synapse/kill/cytokine layer changes
  on a slower timescale. The engine therefore advances PD on a **coarser clock**, handing it the elapsed time
  `dt·pd_every` so that all PD rate constants (per-day) integrate the correct amount of time.
- **Mechanistic rationale.** `dt = 0.02 d` (28.8 min) and `pd_every = 3` → **PD step = 0.06 d (86.4 min)**. This is a
  multi-rate (not a sub-cycled) scheme: PD sees the *instantaneous* per-cell `g.C` at the PD tick, not its average over
  the 3 transport steps. That is a **first-order-accurate sampling approximation** and is the single clearest
  numerical simplification in the driver (see §5.3). The stiff parts (receptor turnover, crosslink) are handled inside
  the callees by their own sub-stepping (`kinetic_rhoden_percell.rhoden_samecell_bivalent_step`, backward Euler).
- **Note (RUN-verified):** a `tsim = 49 d` mosunetuzumab run is `nstep = int(49/0.02)+1 = 2451` transport steps and
  **817** PD steps.
- **Units.** `dt` [day]; `pd_every` [—]; `k_death` [1/day].

---

### EQ-8 — Myeloid IL-6 production aggregation (the census weighting) (`coupled_percell_pd.py:285–296`)

```
il6_prod = Σ_{o ∈ organs}  cs_o · P_o                                   (:286–288)
           + cs_blood · P_blood                                          (:289–292)
           + cs_heme  · P_heme                                           (:293–295)

where   P_o  = pd[o].il6_prod_pg_hr           [pg/hr]  (sampled myeloid, unscaled)
        cs_o = float( pd[o].myeloid_count_scale  or  1.0 )               (:287)
```
- **Biology.** IL-6 in CRS is **myeloid-derived** and **contact-gated**: each macrophage/monocyte agent integrates its
  own local contact with engaged T cells and secretes at its own measured per-cell rate (owned by
  [[T8_mechanistic_crs_il_6]], `myeloid_il6.py`). The agent tables are a **sample** of the tissue, so the sampled
  production must be lifted to the physiological myeloid population before it can be summed into plasma.
- **Mechanistic rationale — why a *separate* scale.** The code is explicit (`:280–284`): the myeloid scale **must not**
  be `graphs[o].count_scale`, which is *antigen*-derived (`tot_antigen_copies / sampled_antigen_copies`,
  `coupled_percell_pk.py:91–94`) and therefore **drug-dependent** (the comment records a measured 0.18×–5.34× spread
  between MS4A1 and TNFRSF17 in the *same* organ). Using it would scale the **identical monocyte population**
  differently for mosunetuzumab and elranatamab — a per-drug artifact that would corrupt a counterscreen whose entire
  job is to *rank* molecules. Myeloid cellularity is a **tissue** property: drug-independent, identical across molecules.
- **The census** is loaded by `wholebody_pd._load_myeloid_scales` (`wholebody_pd.py:42–70`) from
  `handoff/organ_myeloid_counts.json` and attached per organ at `wholebody_pd.py:144`. **Plasma IL-6 is exactly linear
  in this scale**, which is why the driver also records the **unscaled** per-organ production
  (`coupled_percell_pd.py:382–386`) — a run made at scale 1.0 can be re-scaled analytically without re-running.
- **⚠ `or 1.0` DEFECT (found this task).** `float(x or 1.0)` maps `0.0 → 1.0` in Python. The census file sets liver,
  large_int, brain, adipose, heart, kidney, pancreas and skin to **0.0** (intended: *excluded*), but `:287` converts
  that 0.0 into **1.0**, so those organs contribute their **sampled** (unscaled) myeloid production instead of zero.
  Magnitude (RUN-verified this task, myeloid-token scan of every `agents_<organ>.npz`): the eight zero-census organs
  contribute **28,811 sampled agents** in total (liver 18,501 · skin 5,840 · kidney 2,272 · adipose 1,294 · heart 458 ·
  pancreas 436 · brain 6 · large_int 4) at scale 1.0, against **9.53e10** physiological myeloid in the four
  census-scaled organs (spleen 1.70e10 + bone 4.86e10 + lung 2.86e10 + small_int 1.05e9). The leak is therefore
  **≈3.0e-7 of the total — about 6–7 orders of magnitude down**, and does **not** change any reported number
  materially. But the recorded diagnostic `myeloid_count_scale` in the output dict (`:390`) has **no** `or 1.0` and
  therefore reports **0.0** for exactly those organs — the output and the accumulator disagree. This is a latent trap,
  not a live error. See §5.4.
- **Units.** `P` [pg/hr]; `cs` [physiological cells per sampled cell]; `il6_prod` [pg/hr].

---

### EQ-9 — Plasma IL-6 ODE (production rate → concentration) (`coupled_percell_pd.py:297`; body `myeloid_il6.py:206–215`)

```
dC/dt = R − kdeg · C ,      R = il6_prod / V_dist          [pg/mL/hr]
exact step:   C_ss = R/kdeg ;   C ← C_ss + (C − C_ss)·exp(−kdeg·Δt_hr)          (myeloid_il6.py:213–214)
```
- **Biology.** This is the step the engine historically **did not have**: it converts a *production rate* (pg/hr) into a
  *concentration* (pg/mL), which is what a clinical assay measures. Without it, a rate was being multiplied by a fitted
  scalar and reported as a concentration.
- **Mechanistic rationale.** The analytic (exponential) update rather than explicit Euler → **unconditionally stable**
  at the 0.06 d (1.44 h) PD step, where `kdeg·Δt ≈ 0.29` and Euler would already be losing accuracy.
- **⚠ `kdeg = 0.20/hr` is NOT a measured constant** — see §3.6. The in-code comment (`coupled_percell_pd.py:277`,
  `myeloid_il6.py:64`) calls it "the MEASURED first-order IL-6 clearance" and cites **PMID 31268236** with a range
  "0.18–0.25 /hr". Per the repo's own 2026-07-13 provenance audit
  (`reference_unified_binding/IL6_MECHANISM_PROVENANCE.md:30`, `docs/PROVENANCE_AND_VALIDATION.md:191`,
  `IL6_ANCHORS_VERIFIED_2026-07-13.md:31`, sibling doc `T8:455`), that PMID is **Chen 2019, *Clin Transl Sci*** — a
  **semi-mechanistic PK/PD modelling** paper, i.e. the number is a **model-fitted estimate borrowed from another
  model, not a measurement**. Correct tag: **`[FITTED-ELSEWHERE / NOT MEASURED]`**. *Precision note (this audit):
  nothing in **this** engine tunes `kdeg` — it is not fitted here, and the earlier draft's claim that it was "fitted to
  land the plasma peak" was not supported by the code.* The word "MEASURED" in the `PlasmaIL6` docstring
  (`myeloid_il6.py:201`) is nonetheless wrong and should be removed from the source. Because
  `C_peak = production_peak / CL`, production and clearance are **not separately identifiable** from a plasma peak.
  **Between-molecule ratios are clearance-invariant** (CL is a property of IL-6, not of the drug), so the *ranking*
  deliverable survives; the *absolute pg/mL* does not.
- **Units.** `C` [pg/mL]; `R` [pg/mL/hr]; `kdeg` [1/hr]; `V_dist` [mL].

---

### EQ-10 — Blood myeloid gate (`coupled_percell_pd.py:187–190`)

```
self.blood_pd.myeloid.set_count_scale(0.0)          →  cs_i = 0 ∀ i  →  P_blood ≡ 0
```
- **Biology (this is a mechanism decision, not an approximation).** IL-6 induction requires **sustained CD40L–CD40
  contact** ("IL-6 induction and myeloid activation require proximity of CAR T cells and myeloid cells", cited in code
  to **PMID 29808005** — `[UNVERIFIED CITATION]`). In **flowing** blood, leukocytes are ~50 µm apart *and in motion*:
  encounters are transient collisions, not synapses. Sustained contact happens only where myeloid cells are **adherent
  or resident**.
- **Mechanistic rationale (three independent reasons, all in `:175–186`):**
  1. **Double-counting.** The code cites **PMID 3944542** (`[UNVERIFIED CITATION]`) for 60% of blood monocytes being
     **marginating** (adherent to endothelium) — and those marginating/extravasated monocytes **are** the myeloid cells
     already counted in the tissue ABMs. Adding the blood pool would count them twice.
  2. **Geometry artifact.** The blood ABM's coordinates are a **synthetic 2D grid** (mean spacing ~15 µm). Contact in
     that frame is non-physical: the code records a measured run in which it activated **98.6% of 2.0e9 monocytes** and
     produced **61,874 pg/mL** — against a severe-CRS ceiling of ~10–20k pg/mL.
  3. Therefore setting the blood myeloid scale to 0 **removes a geometry artifact; it does not remove biology.**
- **Implementation note.** The gate is applied *inside* `MyeloidIL6.cs` (the per-cell scale array), not at the
  aggregation site, so `P_blood ≡ 0` regardless of what `cs_blood` at `:292` is (it is in fact 1.0 — "blood" is not a
  census key). The gate is wrapped in `try/except` with a printed warning (`:189–190`) — it fails **loud**, not silent.

---

### EQ-11 — Heme malignant burden anchor (`coupled_percell_pd.py:102–127`)

```
homeo_nmol   = Σ_o  Rtot_nM[tgt][o] · Vis_hardcoded[o]                      (:113)   [nmol]
burden_nmol  = HEME_BURDEN_MULT · homeo_nmol            (default 3.0)        (:116–117)
burden_copies= burden_nmol · 1e-9 · AVO                                       (:126)
count_scale  = burden_copies / Σ_{malignant i} R_i                           (:127)
```
- **Biology.** The sampled heme ABM (~27–41k cells) is **not** the patient's tumour burden. The malignant blasts bind
  and internalise circulating drug — a real target-mediated clearance the plasma balance must feel. Its magnitude must
  be anchored to a **physiological antigen amount**, not to the arbitrary ABM sample size.
- **Mechanistic rationale.** The anchor is an **amount (nmol)**, homologous to the organ antigen pools, scaled by a
  disease-burden multiple over the homeostatic whole-body antigen load. Crucially the code scales the **cell count**,
  not the per-cell receptor number (`:118–122`): every blast keeps its **real IHC/Glassman copies** (preserving the
  heterogeneity that drives per-cell binding and avidity), and only the *number of such cells* is lifted. Rejected
  alternative (stated in code, `:101`): a flat nM pool placed in `V_pl`, which is neither units-consistent with the
  organ sinks nor linear-in-C at therapeutic concentrations.
- **⚠ The `_Vis` dict at `:110–112` is HARDCODED and disagrees with the engine's real `Vis`.** Compare (hardcoded vs
  `qsp_costim_window_v2._PBPKArrays.Vis`, RUN-verified this task): spleen 0.0777 **vs 0.038**; large_int 0.226 **vs
  0.074**; heart 0.0088 **vs 0.1056**; skin 0.10 **vs 1.0268**; adipose 0.28 **vs 1.755**; liver 0.30 **vs 0.36**;
  small_int 0.10 **vs 0.13**; kidney 0.047 **vs 0.062**; pancreas 0.016 **vs 0.0324**; lung 0.10 **vs 0.094**; bone 1.0
  **vs 1.0** (only bone agrees). This affects **only** `homeo_nmol` → `burden_nmol` → the heme TMDD sink magnitude. It
  is a **duplicated-constant hazard**: the same physiological quantity exists twice with different values. `[UNSOURCED
  — TBD]` for the hardcoded dict.
- **Units.** `Rtot_nM` [nM]; `Vis` [L]; `homeo_nmol`, `burden_nmol` [nmol]; `count_scale` [—].

---

### EQ-12 — Heme kinetic-Rhoden TMDD sink (`coupled_percell_pd.py:311–323`)

```
surv_h   = exp(−heme_pd.kill_hazard)                                                        (:316)
AgEFF_h  = geo_ageff_nM(R_heme, r_cell=8.0 µm, span=12.5 nm)   if n_arm ≥ 2 else 0          (:318)
(Ag, BAg1, Bdbl, flux) = rhoden_samecell_bivalent_step(C_pl, Ag, BAg1, Bdbl, Ag0,
                                                       kon, koff, AgEFF_h, kdeg, kint, Δt)  (:319–321)
sink[nmol/day] = count_scale · Σ_i [ (flux_i / NM_PER_COPY) · surv_h,i ] / AVO · 1e9        (:323)
```
- **Biology.** Blasts express the target, bind circulating drug at their own receptor number (with **bivalent avidity**
  when `n_arm ≥ 2`), internalise the complex at the antigen's own `kint`, and thereby **clear drug from plasma**. As
  the blasts are killed, the sink must **shrink** — a dead blast does not internalise.
- **Mechanistic rationale.** Three deliberate choices:
  1. **Kinetic, not QSS.** Receptor state (`Ag`, `BAg1`, `Bdbl`) is **carried across steps** (`:131–132, 319–321`) with
     turnover — so TMDD has memory (receptor depletion, resynthesis), which a QSS sink cannot represent.
  2. **Survival-weighted** (`:313–316`): weighting the per-cell flux by `exp(−hazard)` prevents a static full-burden
     sink from over-draining plasma into "a non-physiological terminal cliff at low dose" (code's own words, `:315`).
  3. **Count-scaled, not pool-redistributed** — see EQ-11.
- **`NM_PER_COPY = 6.0/257000 = 2.335e-5 nM/copy`** (`coupled_percell_pd.py:12`) is the *synapse reaction-volume basis*
  shared with the PD engine (`wholebody_pd.py:83`) and the transport binder (`wholebody_percell.py:148`) — **one**
  conversion constant across all three, which is what makes the sink and the PD trimer count the *same* bound receptors.
- **Units.** `flux` [nM/day, synapse basis]; `sink` [nmol/day]; `kint`, `kdeg` [1/day]; `kon` [1/nM/day]; `koff` [1/day].

---

### EQ-13 — Heme is plasma-driven, not extravasation-limited (`coupled_percell_pd.py:305–308`)

```
C_heme = full(len(heme_pd.x), C_pl)                                (:307)
heme_pd.step(C_heme, dt·pd_every, k_death)                          (:308)
```
- **Biology.** Circulating blasts (AML/DLBCL/myeloma in marrow sinusoids and blood) are **not behind a vascular
  barrier**. Every blast sees the plasma concentration directly — this is the physical heme-vs-solid difference, and it
  is why heme TCEs work at doses at which solid TCEs do not. The other heme-vs-solid difference (dense, no-ECM synapse
  geometry) is baked into the heme ABM coordinates themselves (`:56–57`), **not** into a different kill constant:
  `k_death` is *identical* to the organs.
- **Rationale.** Rejected alternative: give heme its own `k_death`. The code is explicit that the depletion gap must
  **emerge** from contact geometry + effector:target ratio + antigen density/KD, not from a per-setting fitted constant.

---

### EQ-14 — Blood: per-lineage count scale + TMDD sink (`coupled_percell_pd.py:160–169`, `:331–336`)

```
scale_i = real_lineage_count[lab_i] / n_sampled[lab_i]          per lineage           (:166–169)
AgEFF_b = geo_ageff_nM(R_blood, 8.0 µm, 12.5 nm)   if n_arm ≥ 2                        (:332)
(…, blood_flux) = rhoden_samecell_bivalent_step(C_pl, …)                               (:333–335)
sink[nmol/day] = Σ_i [ scale_i · blood_flux_i / NM_PER_COPY ] / AVO · 1e9              (:336)
```
- **Biology.** Normal circulating target-positive cells (B cells for CD20; plasma cells for BCMA) are a **real,
  always-present** TMDD sink and the **on-target tox readout** (circulating B-cell depletion). Unlike heme, the burden
  here is the **real physiological cell count** (5 L of blood × lineage densities), not a disease burden.
- **Mechanistic rationale — per-**lineage**, not population-average.** A single mean scale would apply a
  lymphocyte-weighted factor to monocytes (wrong cell type, `:170–174`). The scale is therefore an **array**, one entry
  per cell, computed lineage-by-lineage. This matters for IL-6 too (monocytes are the dominant CRS source), which is why
  `MyeloidIL6.set_count_scale` accepts a full-length per-cell array (`myeloid_il6.py:148–158`) — even though for blood
  it is then overwritten with 0 by EQ-10.
- **⚠ Type hazard, correctly handled.** `self.blood_count_scale` is a **per-cell array**, while the IL-6 aggregator at
  `:292` needs a **scalar** — the code notes that passing the array there would crash `float()` *and* be the wrong scale
  anyway, and instead reads `blood_pd.myeloid_count_scale` (a scalar tissue property). Good defensive wiring.

---

### EQ-15 — Snapshot: bound copies → pericellular nM (`coupled_percell_pd.py:344–351`)

```
v_cell   = max(Vis[o] / g.n, 1e-18)                        [L per cell]         (:345)
bound_nM = S / AVO / v_cell · 1e9                          [nM]                 (:346)
surv     = exp(−pd[o].kill_hazard)                                              (:347)
```
- **Biology.** Converts a per-cell bound **copy number** into the local **concentration** of bound drug in that cell's
  pericellular interstitial volume, so a spatial overlay of "drug bound" is on the same axis as "drug free" (`g.C`).
- **Rationale.** `Vis[o]/n` is the **mean** pericellular volume — the ABM is a 2D sample of a 3D tissue, so a per-cell
  Voronoi volume is not available. `[ASSUMED: uniform pericellular volume]`.

---

### EQ-16 — Plasma mass balance (`coupled_percell_pd.py:357`, `:359`)

```
dA_pl/dt = infn + F_sc·J_sc + k_lymph_return·A_ly
                 − k_cat·A_pl − Σ_o J_extrav[o] − heme_sink − blood_sink
A_pl ← max(A_pl + dt·dA_pl, 0)
```
- **Biology.** Plasma gains: constant-rate infusion, SC absorption, lymph return. Plasma loses: **catabolic clearance**
  (`k_cat = CLup·(1−fFcRn) + k_renal`, `coupled_percell_pk.py:76` — endosomal uptake minus the FcRn-rescued fraction),
  **extravasation** into every organ, and the two **target-mediated** sinks (malignant blasts + normal blood cells).
- **Mechanistic rationale.** FcRn protection enters as a *fraction of uptake that is rescued*, so a molecule's half-life
  is an **emergent** consequence of its `fFcRn` rather than a fitted CL. RUN-verified: mosunetuzumab `fFcRn = 0.89` →
  `k_cat = 0.3503·0.11 = 0.0385 /day` → catabolic t½ ≈ **18.0 d**; elranatamab `fFcRn = 0.70` → `k_cat = 0.105 /day` →
  t½ ≈ **6.6 d** (before extravasation and TMDD, which shorten both).
- **Non-negativity clamp** (`:359`) is a hard floor: it can, in principle, silently absorb mass if the sinks over-drain
  in one step. That is the numerical reason the heme sink is survival-weighted (EQ-12) and the interstitial loss term is
  capped at `0.9·C/dt` inside the organ step (`wholebody_percell.py:192`).
- **Units.** all amounts [nmol]; all rate constants [1/day]; all fluxes [nmol/day].

---

### EQ-17 — Lymph compartment (`coupled_percell_pd.py:358–359`)

```
dA_ly/dt = Σ_o drain[o] − k_lymph_return · A_ly ;   A_ly ← max(A_ly + dt·dA_ly, 0)
```
- **Biology.** All organ lymph drainage collects in one systemic lymph pool that returns to plasma at
  `k_lymph_return = 24 /day` (i.e. mean residence 1 h). This closes the two-pore recirculation loop.
- **Rationale.** A single lumped lymph node compartment; no node-level structure. `[ASSUMED]`.

---

### EQ-18 — Well-mixed side compartments (`coupled_percell_pd.py:353–356`) — **INERT in this driver**

```
Cis   = A_wm[o] / max(Vis[o], 1e-9)
J_ex  = k_dist·L[o]·(1 − sigV[o])·C_pl ;   J_re = k_dist·L[o]·(1 − sigL)·Cis
A_wm[o] ← max(A_wm[o] + dt·(J_ex − J_re), 0)
```
- **Status.** `self.wm` is built from `wellmixed=None` (`coupled_percell_pk.py:78`) and `run_tce_pd_reval.py:142–145`
  never passes it → **`self.wm == {}`** → this block is a no-op in every run of this driver.
- **Consequence (§5.1):** the PBPK tissues **not** in `ORGANS` — **muscle (Vis 3.6 L), brain (0.2175 L), stomach
  (0.0255 L)** — are represented **nowhere**. RUN-verified: the 11 run organs carry **4.678 L** of the model's
  **8.551 L** total interstitial volume = **54.7%**; **45.3% of the interstitium (chiefly muscle) is absent from the PK
  mass balance.** The hook to fix this exists and is unused.

---

### EQ-19 — Recording stride and systemic cytokine sums (`coupled_percell_pd.py:360–386`)

```
record when   k % max(1, nstep // rec) == 0                (rec = 200)                    (:360)
sys_cyto[c]      = Σ_o pd[o].cyto[c]      + blood_pd.cyto[c]          (cumulative)         (:371, :376)
sys_cyto_rate[c] = Σ_o pd[o].cyto_rate[c] + blood_pd.cyto_rate[c]     (instantaneous)      (:372, :376)
il6_plasma_rec  ← PlasmaIL6.C          [pg/mL]                                             (:380)
il6_prod_organ_rec[o] ← pd[o].il6_prod_pg_hr    (UNSCALED, count_scale = 1)                (:385–386)
```
- **Rationale for recording *both* cumulative and instantaneous cytokine:** a clinical IL-6 **peak** is an
  instantaneous concentration; an AUC is the integral. Recording only one makes the other unrecoverable.
- **Rationale for recording UNSCALED per-organ IL-6 production:** plasma IL-6 is **exactly linear** in the myeloid
  census, so an expensive run can be **re-scaled analytically** when the census changes (`:382–384`). This is a
  deliberate, and rather good, design decision: it decouples the expensive simulation from a data gate that is still open.
- **⚠ Asymmetry:** the systemic cytokine sums include organs + **blood** but **not heme** (`:364–376`), whereas the
  mechanistic IL-6 production sum (EQ-8) **does** include heme (`:293–295`). The two cytokine readouts therefore have
  different compartment coverage. See §5.5.

---

### EQ-20 — Derived readouts (`run_tce_pd_reval.py:229–246`)

```
Cmax               = max(Cplasma_ugml)                                            (:229)
il6_peak           = max(il6_plasma_pgml)                                         (:231)
target_organ_kill  = max_o  kill_frac[o][-1]                                      (:242)
depletion_weighted = Σ_o kill_frac[o][-1] · n_target[o]  /  max(Σ_o n_target[o], 1)  (:244–246)
```
- **Biology.** `depletion_weighted` is a **target-cell-count-weighted** whole-body depletion — the right summary for
  "did this molecule deplete the B-cell/plasma-cell compartment?", because an organ with 4 target cells must not carry
  the same weight as marrow.
- **⚠ `Cmax` from a `TSIM_DAYS`-truncated run is a truncated Cmax** and the code says so (`:96`). It is not steady-state.

---

### EQ-21 — Per-molecule kinetics resolution (`run_tce_pd_reval.py:14–33`, `:106–141`)

```
kint(tgt)  = antigen_kinetics_table["membrane"][tgt][0]      (default 0.15 /day)         (:14–18)
kdeg(tgt)  = antigen_kinetics_table["membrane"][tgt][1]      (default 0.5 /day)          (:19–24)
kon        = cfg.kon_TAA_perM_pers   else 1e5 /M/s                                       (:31)
koff       = cfg.koff_TAA_pers       else kon · KD·1e-9                                  (:32)
_kpm       = KINETIC ⊕ {kon/koff for BOTH arms from eng_params_normalized}                (:136–141)
```
- **Biology.** `kint` is a **target property** (CD20 = 0.02/day, non-internalising; BCMA = 2.0/day; GPRC5D = 0.2/day —
  RUN-verified from `handoff/antigen_kinetics_table.json` this task), never a global constant. That single number is
  what makes BCMA TCEs behave differently from CD20 TCEs in the TMDD sink.
- **Mechanistic rationale for the merge (`:101–105`, "FIX-1-UPGRADE").** Before it, the **PK** arm used the measured
  kon/koff while the **PD** synapse was rebuilt from a generic `kon = 1e5` and `koff = kon·KD_registry` — i.e. the same
  antibody bound its target with **two different rate constants** in two parts of the same model. The merge threads the
  measured values into **both**. Unit handling: `eng_params` stores /M/s and /s, exactly what `OrganPD` expects, and the
  PK core converts (`coupled_percell_pk.py:70–71`: `kon[/nM/day] = kon[/M/s]/1e9·86400`; RUN-verified: 1e5 /M/s = **8.64
  /nM/day**; teclistamab 1.2847e6 /M/s = **111.0 /nM/day**, koff 2.315e-4 /s = **20.0 /day**).
- **`koff = kon·KD` is a DERIVED split, not a measurement**, wherever the normalized file lacks a measured koff — and the
  normalized file's own `measured_vs_derived` strings say so explicitly per molecule.
- **⚠ Fallback hazard (dead in this driver).** `coupled_percell_pd.py:134` / `:195` fall back to `kon = 8.64e-3
  /nM/day` if `self.kon1` is None — that is equivalent to **100 /M/s**, i.e. **1000× below** the "standard mAb 1e5 /M/s"
  the code claims elsewhere. It is **unreachable** in this driver (`run_tce_pd_reval.py:144` always passes `kon1_perM_pers`),
  but it is a live landmine for any other caller. `[UNSOURCED — TBD]`.

---

### EQ-22 — Immune motility (opt-in; default OFF) (`coupled_percell_pd.py:211–214`, `:299–304`)

```
if IMMUNE_MOTILITY and (k % MOTILITY_STRIDE == 0):        # inside the  k % pd_every == 0  block
    pd[o].move_immune(dt·MOTILITY_STRIDE, speed=MOTILITY_SPEED_UM_MIN, chemotax=MOTILITY_CHEMOTAX)
```
- **Biology.** T cells migrate (~2–10 µm/min in tissue) and form **new** synapses; a static synapse graph under-counts
  serial killing and under-counts myeloid contact. Default OFF → the synapse graph is frozen at t=0 positions.
- **⚠ Two defects (found this task).**
  1. **Stale comment:** `:212` says `MOTILITY_STRIDE = 100` is "every ~1 day at dt=0.01" — but the driver runs at
     **dt = 0.02** (`run_tce_pd_reval.py:143`), so 100 steps = **2 days**.
  2. **Compounded cadence:** the stride test is *nested inside* the `k % pd_every == 0` block, so motility fires only
     when `k` is divisible by **both** 3 and 100 → every **300 steps = 6.0 days**, not 1 day and not 2 days.
  Both are inert while the feature is default-OFF, but the moment anyone flips `IMMUNE_MOTILITY=1` the motility
  interval is 6× longer than the comment claims.

---

## 3. PARAMETERS OWNED

Everything below is read from the two live files (or from the JSON they load, as marked). Provenance tags are
`[MEASURED]` / `[DERIVED]` / `[FITTED]` / `[ASSUMED]` / `[UNSOURCED — TBD]`, plus `[UNVERIFIED CITATION]` where a comment
names a source I could not confirm.

### 3.1 Driver numerics (the clock and the cadence)

| Symbol | Value | Units | Provenance | Source (file:line) | Mechanistic rationale |
|---|---|---|---|---|---|
| `dt` | **0.02** | day (28.8 min) | `[ASSUMED: numerical]` | `run_tce_pd_reval.py:143` (constructor arg); env `PD_DT` override at `coupled_percell_pk.py:74` | Transport clock. Must resolve the fastest *explicit* process; the stiff binding/diffusion pieces are handled implicitly or sub-stepped inside callees. No convergence study is present in the code. |
| `pd_every` | **3** | steps | `[ASSUMED: numerical]` | `run_tce_pd_reval.py:173` | PD clock = 0.06 d. Multi-rate split (EQ-7). No stated justification in code. |
| `rec` | **200** | records | `[ASSUMED]` | `run_tce_pd_reval.py:173` | Output resolution; stride = `nstep//rec` (`coupled_percell_pd.py:360`). |
| `snap_times` | **[1, 7]** | day | `[ASSUMED]` | `run_tce_pd_reval.py:171` (`SNAP_DAYS` env) | Spatial snapshot days. |
| `k_death` | **1.0** | 1/day | `[ASSUMED: "LOCKED"]` — provenance string in `handoff/kinetic_calib.json` claims a *primary* anchor of the engaged-CTL serial ceiling `k_hit·koff/(k_hit+koff) = 11.6/day` and "NO fitting", but its *secondary* check is stated against the **retired 570 pg/mL** IL-6 anchor | `run_tce_pd_reval.py:41, 173`; `pd_model_config.py:59, 63` | Trimer→death rate constant. One shared value across all engagers (deliberately: a per-setting `k_death` would hide mechanism in a fitted constant). |
| `k_hit` | **12.0** | 1/day | `[LITERATURE-CLAIMED — NO PMID IN CODE]` (the comment asserts "FIXED from serial-killing literature, not fitted", but names no paper) → `[UNVERIFIED CITATION]`. **Not** `[MEASURED]`: nothing in the code substantiates it. | `pd_model_config.py:41`; live value loaded from `handoff/kinetic_calib.json` (`k_hit: 12.0`) at `pd_model_config.py:59`, fallback literal `:63` | ~1 lethal hit / 2 h. Consumed by T5. |
| `iv_inf_h` | **2.0** | h | `[ASSUMED: clinical TCE infusion ~2 h]` | `run_tce_pd_reval.py:176` | Zero-order infusion duration (EQ-2). |
| `F_sc` | **0.6** | — | `[ASSUMED]` (no source in code) | `run_tce_pd_reval.py:79–82, 175`; default `coupled_percell_pd.py:199` | SC bioavailability. Applied as a *loss* (EQ-4). |
| `ka_sc` | **0.25** | 1/day | `[ASSUMED]` (no source in code) | same as above | SC absorption rate (t½ absorption ≈ 2.8 d). |
| `span_coeng_nm` | **12.5** | nm | `[ASSUMED: default]` — molecules may override via `cfg['span_coeng_nm']`/`span_bridge_nm` | `run_tce_pd_reval.py:145` | Bivalent same-antigen co-engagement span → `geo_ageff_nM` (T4). |
| `KD_CD3` | **40.0** | nM | `[ASSUMED: schema default]`; overridden per molecule by `eng_params_normalized.json` (e.g. teclistamab 28.03, elranatamab 17.0, glofitamab 4.5) | `run_tce_pd_reval.py:152, 162, 167`; `coupled_percell_pd.py:17` | CD3 arm affinity. |
| `KD_costim` | **1.0** | nM | `[ASSUMED: placeholder]` | `run_tce_pd_reval.py:153`; `coupled_percell_pd.py:18` | Costim arm affinity; only used when `costim_arm` is set. |
| `_NM_PER_COPY_PD` | **6.0/257000 = 2.335e-5** | nM/copy | `[DERIVED: pinned by the validated tumour kill anchor Rcap_TAA = 6.0 nM at CEACAM5 257,000 copies/cell]` (rationale in `wholebody_pd.py:75–83`) → equivalently a ~71 pL synapse reaction volume | `coupled_percell_pd.py:12` | The ONE conversion linking copies ↔ synapse nM. Shared identically by transport, PD and both TMDD sinks — that identity is what stops the same bound receptor being counted twice with two different constants. |
| `geo_ageff` cell radius | **8.0** | µm | `[ASSUMED]` (hardcoded at both sink call sites) | `coupled_percell_pd.py:318, 332` | Blast/leukocyte radius for the Rhoden geometric 2nd-arm concentration. |
| `geo_ageff` span | **12.5** | nm | `[ASSUMED]` (hardcoded — **does not** read `self.span_coeng_nm`) | `coupled_percell_pd.py:318, 332` | ⚠ The organ transport path uses the per-molecule `span_coeng_nm` (`wholebody_percell.py:178`), but the heme/blood sinks hardcode 12.5 — a per-molecule span override silently does **not** reach the heme/blood sinks. |

### 3.2 The ENG molecule registry (`run_tce_pd_reval.py:53–89`)

| Molecule | tgt | KD (nM) | mw (kDa) | fFcRn | n_arm | route | `il6_obs` (pg/mL) | Provenance of `il6_obs` |
|---|---|---|---|---|---|---|---|---|
| mosunetuzumab | MS4A1 | 5.0 | 146.0 | 0.89 | 1 | IV | **152.0** | **`[UNSOURCED — TBD]`.** The code calls it a population MEAN (`:75`) but gives **NO PMID and no citation of any kind**. It must **not** be tagged `[MEASURED]` — nothing in the code substantiates the value. Display/scoring only. |
| glofitamab | MS4A1 | 5.0 | 195.0 | 0.89 | **2** (2:1 bivalent CD20) | IV | **None** | no valid anchor exists in code |
| epcoritamab | MS4A1 | 5.0 | 148.0 | 0.89 | 1 | SC (F=0.6, ka=0.25) | **None** | — |
| elranatamab | TNFRSF17 | 0.15 | 148.5 | **0.70** (IgG2) | 1 | SC | **None** | **No clinical IL-6 value exists for this molecule** (code `:71`) |
| teclistamab | TNFRSF17 | 0.15 | 146.0 | 0.89 | 1 | SC | **21.0** | `[MEASURED: PMID 38831634 — "The mean IL-6 peak concentration (Cmax) was 21 pg/mL", quoted verbatim at :74]` `[UNVERIFIED CITATION]` |
| talquetamab | GPRC5D | 2.0 | 146.0 | 0.89 | 1 | SC | **None** | the available 18.2 is a **MEDIAN**, not a mean → not comparable to the mean pair; excluded (`:82`) |
| mosunetuzumab_sc | MS4A1 | 5.0 | 146.0 | 0.89 | 1 | SC | 152.0 | PK route-matched variant; anchor display-only |
| teclistamab_iv | TNFRSF17 | 0.15 | 146.0 | 0.89 | 1 | IV | 21.0 | PK route-matched variant |
| teclistamab_iv_low | TNFRSF17 | 0.15 | 146.0 | 0.89 | 1 | IV | 21.0 | PK route-matched variant |

- **Registry `KD` is the PK-arm KD.** The **PD** arm uses `cfg['KD_norm']` — the **measured** TAA KD merged in from
  `eng_params_normalized.json` when present (`run_tce_pd_reval.py:122, 152`). These differ (e.g. elranatamab registry
  0.15 vs measured **0.04**; teclistamab registry 0.15 vs measured **0.18**; mosunetuzumab registry 5.0 vs measured
  **5.45**). `[DERIVED: registry values are class/rounded; KD_norm is the measured one]`. **A reviewer will ask why the
  PK sink and the PD synapse see different KDs — they do, by construction, and the merge only fixes kon/koff, not the
  registry `KD` passed to `CoupledPerCellPD(...)` at `:142`.** See §5.6.
- `fFcRn` values: `[ASSUMED: isotype class values]` (0.89 IgG1-class, 0.70 IgG2 per code comment at `:80`). No PMID in code.
- `mw` elranatamab 148.5 "(Elrexfio label)" `[UNVERIFIED CITATION]`; others `[ASSUMED: class values]`.

### 3.3 Simulation horizons and step-up regimens

`TSIM` (`run_tce_pd_reval.py:90`), days: mosunetuzumab **49**, glofitamab **24**, epcoritamab **24**, elranatamab **24**,
teclistamab **200**, talquetamab **24**, mosunetuzumab_sc **160**, teclistamab_iv **10**, teclistamab_iv_low **10**.
All overridable by `TSIM_DAYS` (`:97–100`). `[ASSUMED: per-molecule, to reach PK steady state]`.

Live schedules — the file the driver **actually loads** is `{KWS}/handoff/regimen_schedules_final.json`
(`run_tce_pd_reval.py:44, 169`), day → mg; **0-based days** in the file, i.e. day 0 = clinical C1D1.

> **⚠ Duplicated-file hazard (found this task).** There is a second copy at `engine/regimen_schedules_final.json`
> that is **never loaded** and is **not the same file**: the handoff copy is a strict superset (it additionally
> carries acapatamab, alnuctamab, blinatumomab, catumaxomab, cevostamab, cibisatamab, linvoseltamab, odronextamab,
> pasotuxizumab, runimotamab, solitomab, tebentafusp, …). RUN-verified this task: for the 9 molecules in `ENG` the
> two copies agree exactly (the diff is pure addition), so every number in the table below is correct — but **quote
> the handoff path, not the engine path**, and do not edit the engine copy expecting a run to change.

| Molecule | Schedule (day, mg) | The trap |
|---|---|---|
| **mosunetuzumab** | (0, **1**) · (7, **2**) · (14, **60**) · (21, 60) · (42, 30) · (63, 30) · (84, 30) · (105, 30) | **A 7-day sim gives ONLY the 1 mg priming dose.** The 2 mg dose lands on the last step (t = 7.0, no time to act); the **60 mg** dose is at **t = 14 d** (clinical C1D15) and never occurs. |
| **teclistamab** | (0, 4.2) · (3, 21) · (7, **105**) · then 105 q7d to day 210 | At 7 days it has its **full 105 mg** RP2D. |
| elranatamab | (0, 12) · (3, 32) · (7, **76**) · then 76 q7d | full dose by day 7 |
| epcoritamab | (0, 0.16) · (7, 0.8) · (14, **48**) · then 48 q7d | full dose only from day 14 |
| glofitamab | (7, 2.5) · (14, 10) · (21, **30**) · then 30 q21d | **no dose at all before day 7** |
| talquetamab | (0, 0.7) · (3, 4.2) · (7, **28**) · then 28 q7d | full dose by day 7 |

> **Consequence.** Any cross-molecule comparison at a horizon shorter than ~24 d is comparing **different rungs of
> different ladders**. Mosunetuzumab and epcoritamab reach their therapeutic dose at day 14; glofitamab at day 21;
> teclistamab/elranatamab/talquetamab at day 7. The code's own recommendation is `TSIM_DAYS=24` (`:93–96`).

### 3.4 Per-organ myeloid census (`handoff/organ_myeloid_counts.json` → `wholebody_pd.py:42–70, 144`)

The file's values are **SCALES** (physiological myeloid count ÷ sampled myeloid agents), despite the file name saying
"counts". RUN-verified this task by counting myeloid-labelled agents in each `agents_<organ>.npz` with the module's own
token list (`macrophage|monocyte|myeloid|kupffer|microglia`, `myeloid_il6.py:122`):

| Organ | census scale (file) | sampled myeloid agents (RUN-verified) | ⇒ implied physiological macrophages | Provenance |
|---|---|---|---|---|
| **spleen** | **290,206.39** | 58,579 | **1.70e10** | **`[LITERATURE-CLAIMED — NO SOURCE IN THE CODE OR THE JSON]`.** The only provenance anywhere is an external repo doc, `reference_unified_binding/IL6_MECHANISM_PROVENANCE.md:86`, which attributes the census to **PMID 37871201 (Sender 2023 PNAS)** — verified this task to exist at exactly that file:line, and its four values (spleen 1.70e10 · bone 4.86e10 · lung 2.86e10 · small_int 1.05e9) reproduce the scale×agents products below exactly. **`[UNVERIFIED CITATION]`: I confirmed the arithmetic and the in-repo attribution, NOT the paper.** Do **not** call this `[MEASURED]` — neither live file nor the JSON carries a provenance field. |
| **bone** | **33,633,218.0** | 1,445 | **4.86e10** | same |
| **lung** | **1,083,826.0** | 26,388 | **2.86e10** | same |
| **small_int** | **84,910.24** | 12,366 | **1.05e9** | same |
| liver | **0.0** | 18,501 | — | **EXCLUDED: Kupffer-cell census could not be sourced.** |
| large_int | 0.0 | 4 | — | EXCLUDED: only **4** myeloid agents sampled — unrepresentative (the intestine is macrophage-rich in reality) |
| brain | 0.0 | 6 | — | EXCLUDED: 6 agents; microglia not sampled |
| adipose / heart / kidney / pancreas / skin | 0.0 | 1,294 / 458 / 2,272 / 436 / 5,840 | — | not in the sourced census |
| blood | (not a census key) | 74,207 | — | **GATED OFF by mechanism** (EQ-10) |
| heme_tumor | (not a census key ⇒ 1.0) | myeloma ABM **0**; DLBCL ABM **567**; AML ABM **27,931** | — | contributes at **sampled** counts (scale 1.0) |

> **⚠ THE LIVER EXCLUSION BIASES IL-6 DOWNWARD, AND WE SAY SO.** The liver's Kupffer cells are the body's largest
> resident-macrophage population. They are **absent from the census** (scale 0.0) because no citation could be sourced
> for their count. Every IL-6 number this engine produces is therefore a **LOWER BOUND** with respect to hepatic
> myeloid contribution. This is a **stated data gap**, not a modelling choice — and it must never be "fixed" by tuning
> another constant upward to compensate.
>
> Note the interaction with the `or 1.0` defect (EQ-8): the excluded organs do not contribute **zero**; they contribute
> their **sampled** production at scale 1.0. All eight zero-census organs together contribute 28,811 sampled myeloid
> agents' worth of production against 9.53e10 physiological myeloid in the four scaled organs — **≈3.0e-7 of the total,
> ~6–7 orders of magnitude down**, hence numerically negligible, but conceptually the "exclusion" is not what the file
> says it is.

### 3.5 Heme / blood compartment parameters

| Symbol | Value | Units | Provenance | Source | Rationale |
|---|---|---|---|---|---|
| `HEME_TMDD_SINK` | **1** (on) | bool | `[ASSUMED]` | `coupled_percell_pd.py:86` | env kill-switch |
| `HEME_BURDEN_MULT` | **3.0** | × homeostatic antigen | **`[ASSUMED: "a heavy myeloma/lymphoma carries a few-fold more" — no source]`** | `coupled_percell_pd.py:116` | Sets the whole magnitude of the heme TMDD sink. A reviewer will attack this: it is the single least-grounded number in the heme path. |
| `_Vis` (hardcoded) | spleen .0777, bone 1.0, large_int .226, lung .10, kidney .047, liver .30, pancreas .016, small_int .10, heart .0088, skin .10, adipose .28, brain .10 | L | **`[UNSOURCED — TBD]` and INCONSISTENT with the engine's real `Vis`** (see EQ-11) | `coupled_percell_pd.py:110–112` | Used only for `homeo_nmol`. |
| `heme` fallback `kdeg` | 0.5 | 1/day | `[ASSUMED: fallback]` | `coupled_percell_pd.py:133` | reached only if `self.kdeg` is None |
| `heme`/`blood` fallback `kon` | **8.64e-3** | 1/nM/day | `[UNSOURCED — TBD]` — equals **100 /M/s**, 1000× below the code's own "standard mAb" 1e5 /M/s | `coupled_percell_pd.py:134, 195` | **Dead in this driver** (`kon1` always supplied). Landmine for other callers. |
| `BLOOD_TMDD_SINK` | **1** (on) | bool | `[ASSUMED]` | `coupled_percell_pd.py:157` | env kill-switch |
| `USE_BLOOD_COMPARTMENT` | **1** (on) | bool | `[ASSUMED]` | `run_tce_pd_reval.py:166` | blood is always attached |
| blood myeloid scale | **0.0** | — | `[DERIVED: mechanism — sustained CD40L–CD40 contact impossible in flowing blood; marginating monocytes already counted in tissue]` (EQ-10); cited to **PMID 29808005** and **PMID 3944542** `[UNVERIFIED CITATION]` | `coupled_percell_pd.py:188` | Removes a measured 61,874 pg/mL synthetic-grid artifact |
| blood `real_lineage_counts` | from `agents_blood.npz` | cells | `[UNSOURCED — TBD]` in this file (the npz carries them; no PMID in code) | `coupled_percell_pd.py:162–163` | per-lineage count scale (EQ-14) |

### 3.6 IL-6 chain constants consumed by this driver (owned by [[T8_mechanistic_crs_il_6]])

| Symbol | Value | Units | Provenance | Source |
|---|---|---|---|---|
| `KDEG_IL6_PER_HR` | **0.20** | 1/hr | **`[FITTED-ELSEWHERE / NOT MEASURED]`** — in-code (`myeloid_il6.py:64`, `coupled_percell_pd.py:277`) it is called "MEASURED" and cited to **PMID 31268236** (range "0.18–0.25 /hr"). Per the repo's own audit (`reference_unified_binding/IL6_MECHANISM_PROVENANCE.md:30`; `docs/PROVENANCE_AND_VALIDATION.md:191`; T8:455) that PMID is **Chen 2019, *Clin Transl Sci*** — a **semi-mechanistic MODELLING** paper, so the value is a **model-fitted estimate borrowed from another model**. `[UNVERIFIED CITATION]` — I read the attribution in the repo, I did not open the paper. **Not fitted inside this engine.** Since `C_peak = production/CL`, production and clearance are not separately identifiable from a peak. | `myeloid_il6.py:64`, used at `coupled_percell_pd.py:297` |
| `V_PLASMA_ML` | **11,650** | mL | `[DERIVED: the model's OWN PBPK volumes — Σ Vis (8.5508 L) + V_pl (3.10 L) = 11.6508 L; RUN-verified this task from `qsp_costim_window_v2._PBPKArrays`]` | `myeloid_il6.py:65` |
| `S_MAX_MOLEC_PER_S` | **10.6** | molec/s/cell | `[MEASURED: PMID 20376398 — mean over actively-secreting monocytes]` `[UNVERIFIED CITATION]` | `myeloid_il6.py:43` |
| `SECRETOR_FRACTION` | **0.039** | — | `[MEASURED: PMID 37533643 — cell-intrinsic, measured on fully LPS-stimulated cells]` `[UNVERIFIED CITATION]` | `myeloid_il6.py:53` |
| `R_CONTACT_UM` | **14.1** | µm | `[DERIVED: r_macrophage 10.6 (PMID 9400735) + r_Tcell 3.5 (PMID 30571054)]` `[UNVERIFIED CITATION]` | `myeloid_il6.py:101` |

> **⚠ Stale comment (cross-file, flagged for the T8 owner):** `myeloid_il6.py:111` and `wholebody_pd.py:124` both say the
> per-cell secretion is "**~0.0196 pg/hr/cell**". At the **live** `S_MAX = 10.6 molec/s` the value is **0.00133
> pg/hr/cell** (RUN-verified: 10.6 × 21000 / 6.022e23 × 1e12 × 3600). 0.0196 pg/hr/cell corresponds to the **retired
> 156 molec/s** peak-tail rate. The **code is right; the comments are stale by 14.7×.** Do not quote 0.0196.

### 3.7 Inherited PK defaults active in every run (owned by [[T1_shah_betts_pbpk_backbone]] / [[T2_whole_body_per_cell_pk]])

Never overridden by `run_tce_pd_reval.py` — they are the *silent* physiology of every result.

| Symbol | Value | Units | Provenance | Source |
|---|---|---|---|---|
| `V_pl` | 3.1 | L | `[ASSUMED: standard human plasma volume]`; `[UNVERIFIED CITATION]` "Shah-Betts" | `coupled_percell_pk.py:65` |
| `V_ly` | 2.6 | L | `[ASSUMED]` (declared but the lymph ODE uses amounts, not `V_ly`) | `coupled_percell_pk.py:65` |
| `CLup` | 0.3503 | L/day | `[UNVERIFIED CITATION: "Shah-Betts"; no PMID in code]` | `coupled_percell_pk.py:65, 76` |
| `k_renal` | 0.0 | L/day | `[ASSUMED: IgG is not renally cleared]` | `coupled_percell_pk.py:65` |
| `sigL` | 0.85 | — | `[UNVERIFIED CITATION: "Shah-Betts"]` | `coupled_percell_pk.py:64` |
| `k_dist` | 3.0 | — | `[UNSOURCED — TBD]` | `coupled_percell_pk.py:64` |
| `k_lymph_return` | 24.0 | 1/day | `[ASSUMED: 1 h lymph residence]` | `coupled_percell_pk.py:65` |
| `D_um2s` | 10.0 | µm²/s | `[ASSUMED]` | `coupled_percell_pk.py:65` |
| `ORGANS` | spleen, bone, large_int, liver, lung, small_int, pancreas, kidney, skin, heart, adipose (+`tumor` if solid) | — | `[ASSUMED: selection]` — **muscle, brain, stomach omitted** | `run_tce_pd_reval.py:48, 129` |

### 3.8 Vestigial / contaminated constants still present in live code (**no consumer in the live import graph** — verified by grep this task)

| Constant | Value | Where | Status |
|---|---|---|---|
| `IL6_SCALE` | 0.10678 (QSS, from `handoff/wholebody_cyto_calib.json:il6_eng_scale_pgml_per_raw` = 0.10677852891678402 — itself constructed as 570/peak_raw) / 0.05473 (kinetic, from `handoff/kinetic_calib.json:IL6_SCALE_kin` = 0.054730085) | `run_tce_pd_reval.py:34, 41`; `pd_model_config.py:45, 59, 63` | **LOADED BUT NEVER READ ON THE LIVE PATH.** The fitted legacy IL-6 path was deleted 2026-07-13 (`run_tce_pd_reval.py:203–214`) and replaced by a **hard error** if the mechanistic recorder is empty (`:216–220`) — the old code silently substituted the fitted value under the same field name. **Correction to an earlier draft of this doc: it is *not* "consumed by nothing"** — one module outside the live import graph still reads it (`H.IL6_SCALE`), so deleting the assignment is a *breaking* change for that caller, not a free one. Its QSS value is *by construction* an anchor to the fabricated 570 (`anchor_pgml: 570.0`). **Recommend deletion + fixing that one caller**, so it cannot be re-wired. |
| `CYTO_IL6_CLINICAL_ANCHOR_PGML` | **570.0** | `wholebody_pd.py:23` | **DEFINED, NEVER READ.** 570 is one of the **fabricated** anchors — **it has no source.** Its comment cites "Hosseini 2020 Fig5A" → `[UNVERIFIED CITATION]`; treat as contamination. |
| `cytokine_to_pgml()` | — | `wholebody_pd.py:32–35` | **DEFINED, NEVER CALLED.** Hence `cyto_sys_rate`/`cyto_sys_cum` in the output are **raw units**, not pg/mL. |
| `clin_il6: 570.0`, `il6_raw_peak`, `mosun_kill{}` | — | `handoff/kinetic_calib.json` | Calibration record whose **only** live consumers are `k_hit` and `k_death`. The 570 anchor inside it is **retired**; `k_death`'s own provenance string still cites it as a *secondary check*. |

**Values that must NEVER be re-introduced** (provenance audit 2026-07-13, encoded in `run_tce_pd_reval.py:58–76`):
`570` (mosunetuzumab — **no source exists**), `340` (elranatamab — cited to a MagnetisMM-3 figure that **does not exist**;
that paper contains **zero** mentions of IL-6 and has only four figures), `191` (elranatamab — **a PAGE NUMBER**, a
dot-leader from the Table of Figures of FDA BLA 761345, used as a concentration), `230` and `366.88` (elranatamab —
**no source**), `288` (teclistamab — **real** but the **highest individual patient Cmax**, an *order statistic* that
scales with cohort N and is therefore not comparable across trials, and not a central-tendency anchor).

---

## 4. WHAT IS EMERGENT vs IMPOSED

### 4.1 EMERGENT (computed from mechanism inside the loop)

| Quantity | Emerges from |
|---|---|
| **Plasma PK curve** (Cmax, half-life, distribution phase) | dose schedule (EQ-1/2/4) + FcRn-set catabolism `k_cat = CLup(1−fFcRn)` + Σ organ extravasation + **target-mediated** heme/blood sinks (EQ-16). *No fitted CL, no fitted Vd.* Elranatamab's shorter catabolic t½ (6.6 d vs mosunetuzumab's 18.0 d) is a **consequence** of its IgG2 `fFcRn = 0.70`, not an input. |
| **Per-organ drug exposure** | vascular QSS → BEC extravasation → per-cell graph diffusion → LEC drainage (T2). Nothing per-organ is dialled. |
| **TMDD nonlinearity** | per-cell Rhoden binding at **real** receptor numbers + `kint` (a *target* property: CD20 0.02/day vs BCMA 2.0/day) + receptor turnover. The heme sink **shrinks as blasts die** (survival-weighted, EQ-12). |
| **Killing** (per organ, per cell) | per-cell trimer / kinetic synapse (T5/T6) driven by that cell's own local `g.C`, its own CD3/TAA copies, its own Treg neighbourhood. |
| **IL-6 magnitude, saturation, and per-molecule differences** | **where** each drug engages: the finite, spatially-distributed myeloid pool + contact gating. There is **no Emax, no EC50, no Hill, no fitted scale** in the production term. Saturation emerges from the per-cell `(1−a)` ceiling and the finite pool. |
| **Depletion ordering across organs** | count-weighted per-organ kill (EQ-20) — an organ with 4 target cells cannot outvote marrow. |
| **Heme-vs-solid depletion gap** | contact geometry (dense no-ECM ABM) + effector:target ratio + antigen density/KD. Explicitly **not** a per-setting `k_death` (`coupled_percell_pd.py:56–58`; `kinetic_calib.json:k_death_provenance`). |

### 4.2 IMPOSED (handed to the loop as constants)

| Quantity | Why it is imposed | Honest severity |
|---|---|---|
| **IL-6 clearance `kdeg = 0.20/hr`** | Nothing in the model computes cytokine disposition. | **SEVERE — a MODEL-FITTED constant (fitted in the cited paper, not measured) carrying a PMID and the word "MEASURED" in the source.** It sets the *absolute* pg/mL scale entirely (`C_peak = production/CL`). Between-molecule **ratios are invariant to it**; absolute values are not. |
| **IL-6 distribution volume 11.65 L** | Lumped ECF compartment. | The mechanistically complete version (route IL-6 through the antibody's own PBPK transport with IL-6's reflection coefficient) is **built in this engine and currently bypassed** (`myeloid_il6.py:77–81`). A named, open uncertainty. |
| **Myeloid census scales** | External data (`organ_myeloid_counts.json`). | Data gate, **exactly linear** → analytically re-scalable. Liver = 0 ⇒ IL-6 is a lower bound. |
| **`k_death = 1.0 /day`** | One shared constant for all engagers. | Locked, not per-molecule → it cannot manufacture a per-molecule ordering. But it does set the absolute depletion scale. |
| **`HEME_BURDEN_MULT = 3.0`** | Disease burden vs homeostatic antigen. | **Unsourced.** Directly scales the heme TMDD sink → directly scales heme-target PK. |
| **`F_sc = 0.6`, `ka_sc = 0.25/day`** | SC absorption. | Unsourced generic values applied to **all four SC molecules identically** — so SC-vs-SC comparisons are internally consistent, but SC-vs-IV comparisons inherit them. |
| **Regimens, MW, fFcRn, KD, kon/koff, kint** | Read from data files. | Provenance is per-molecule and *documented in the data file itself* (`measured_vs_derived` strings, read this task) — the honest answer is "mixed", and **thinner than the earlier draft implied**: only **teclistamab** has both arms MEASURED ("BCMA & CD3 kon/koff: MEASURED-in-QSP", Janssen translational QSP). **Elranatamab's kon is measured (Yoneyama Table 2) but its `koff_BCMA` is explicitly DERIVED = kon·KD.** Mosunetuzumab/epcoritamab/glofitamab CD20 kon are *class*-measured (rituximab/anti-CD20 SPR 4.3e5) with **DERIVED** koff and **class-estimate** KDs; talquetamab's affinities are a flagged cross-reference from a "same range as forimtamig" statement. Several CD3 KDs are schema-default 40 nM. |
| **`dt = 0.02 d`, `pd_every = 3`** | Numerics. | No convergence study exists in code. |

### 4.3 Where emergence stops — the one-sentence version

**Everything upstream of the plasma IL-6 ODE is emergent; the ODE's clearance constant is not.** The engine computes,
from per-cell mechanism, a *production rate* in pg/hr. Converting that to the pg/mL a clinician measures requires a
clearance the literature does not appear to contain. The defensible deliverable is therefore a **relative / rank-order
CRS-risk counterscreen** (clearance cancels), **not an absolute pg/mL predictor**.

---

## 5. KNOWN LIMITATIONS / OPEN QUESTIONS

### 5.1 45% of the interstitium is not in the model
`ORGANS` (`run_tce_pd_reval.py:48`) omits **muscle (Vis = 3.6 L)**, **brain (0.2175 L)** and **stomach (0.0255 L)**, and
the well-mixed hook that exists to carry them (`coupled_percell_pd.py:353–356`) is **never populated** (`wellmixed=None`).
RUN-verified: the 11 run organs hold **4.678 L of the model's 8.551 L** interstitial volume — **54.7% coverage**. The
missing extravascular volume is a **missing drug sink**, so plasma concentrations are biased **upward** and the
distribution phase is too shallow. *A reviewer will spot this immediately: an IgG PBPK without muscle.* Fixing it is
mechanical (pass `wellmixed={o:(L[o],sigV[o],Vis[o])}` for the three missing tissues).

### 5.2 The horizon trap is a live foot-gun, not a historical one
`TSIM_DAYS` is an env var with no guard rail. Nothing in the code prevents a 7-day run and nothing annotates the output
JSON with "this molecule had only reached rung 1 of its step-up". The **only** protection is a print statement (`:100`).
**Recommendation:** the driver should refuse (or loudly annotate) any run whose horizon does not include the molecule's
maximum scheduled dose. Every cross-molecule claim must state the horizon and the dose rung reached.

### 5.3 The PD cadence is a sampling approximation, unvalidated
PD is advanced with `dt·pd_every` but is handed the **instantaneous** `g.C` at the tick, not its 3-step average. During
an IV infusion — precisely when `C` is changing fastest and CRS is being set — this is where the error is largest. There
is **no convergence study** in the repo (no `pd_every=1` vs `3` comparison, no `dt` halving). This is cheap to do and
should be done before any number is defended.

### 5.4 The `or 1.0` census defect
`float(getattr(...,"myeloid_count_scale",1.0) or 1.0)` (`coupled_percell_pd.py:287`) maps a **legitimate 0.0** to
**1.0**. Organs the census deliberately excludes (liver, large_int, brain, …) therefore contribute their **sampled**
myeloid production rather than none — while the output diagnostic (`:390`, no `or`) reports **0.0** for them. The
numerical impact is negligible (**≈3.0e-7 of the total — ~6–7 orders of magnitude below the scaled organs**;
RUN-verified), but **the code does not do what the data file says**, and the two places disagree with each other.
Fix: `_cs = float(x) if x is not None else 1.0`.

### 5.5 Two cytokine readouts with different units and different compartment coverage
- `il6_plasma_pgml` — **pg/mL, physical**, mechanistic, includes organs + heme (blood contributes 0 by EQ-10).
- `cyto_sys_rate` / `cyto_sys_cum` (IL6/IFN/TNF/IL2) — **raw engagement-sum units**, includes organs + **blood**, excludes
  **heme**. The pg/mL converter is defined and never called (§3.8).
So the two "IL-6" numbers in the output JSON are **not the same quantity, not in the same units, and not summed over the
same compartments**. Only `il6_plasma_pgml` may be compared to a clinical value. IFN/TNF/IL2 have **no** physical
calibration at all and must not be presented as concentrations.

### 5.6 The PK arm and the PD arm can see different KDs
`CoupledPerCellPD(..., cfg['KD'], ...)` (`run_tce_pd_reval.py:142`) passes the **registry** KD to the transport TMDD
binder, while `attach_pd(..., KD_TAA_nM=cfg.get('KD_norm', cfg['KD']), ...)` (`:152`) passes the **measured** KD to the
PD synapse. For elranatamab that is **0.15 nM (PK) vs 0.04 nM (PD)** — a 3.75× discrepancy in the same molecule in the
same run. The `kon/koff` merge (FIX-1-UPGRADE) unified the *rate constants* but left this KD asymmetry in place. Either
it is deliberate (registry KD is a deliberately-pinned PK value) and should be documented as such, or it is a residue of
the pre-merge world and should be closed.

### 5.7 Heme burden magnitude rests on one unsourced multiplier
`HEME_BURDEN_MULT = 3.0` (`:116`) × a `homeo_nmol` computed with a **hardcoded interstitial-volume dict that disagrees
with the engine's own `Vis`** (EQ-11). The heme TMDD sink — and therefore the plasma PK of every BCMA/CD20 heme molecule
— is linear in this product. **Two independent unsourced/duplicated quantities multiply here.** This is the most
attackable number in the heme path.

### 5.8 The heme/blood sinks hardcode the co-engagement span
`geo_ageff_nM(R, 8.0, 12.5)` at `:318` and `:332` ignores `self.span_coeng_nm` (which the organ transport path *does*
honour, `wholebody_percell.py:178`). A per-molecule format/span override (the whole point of T4) therefore reaches the
tissue avidity but **not** the circulating-compartment avidity.

### 5.9 Heme myeloid are unscaled; blood myeloid are zeroed; solid tumour myeloid are unscaled
Only 4 organs have a sourced myeloid census. The heme compartment (which for the AML ABM contains **27,931** myeloid
agents) contributes IL-6 at **scale 1.0** (`:295`, no census key), i.e. at *sampled* counts — an arbitrary number tied
to the ABM sample size, exactly the failure mode the count-scaling was invented to avoid. For the two clinically
validated heme molecules this happens to be near-harmless (myeloma ABM has **0** myeloid agents; DLBCL ABM has 567), but
an AML/CD123 run would silently take an unscaled, sample-size-dependent myeloid IL-6 contribution.

### 5.10 The costim path this driver now threads CANNOT RUN — `signaling_dynamics` does not exist
**(found by this audit; the most actionable defect in T9.)**
`run_tce_pd_reval.py:146–153` (added 2026-07-13) fixed a real bug: the runner "never PASSED" `costim_arm`, so "every
costim construct silently ran as a plain TCE". It now threads `cfg['costim_arm']` into `attach_pd`. But the code that
consumes it, `wholebody_pd.py:208–209`, does:
```python
if costim_arm is not None and R_costim_percell is not None:
    import signaling_dynamics as _sigmod          # <-- unguarded
```
and **`signaling_dynamics.py` does not exist anywhere in the repository.** RUN-verified this task:
`importlib.util.find_spec('signaling_dynamics')` is `None` on the driver's own `sys.path`
(`{KWS}`, `{KWS}/handoff`, `engine/`), while `costim_induction`, `kinetic_synapse`, `multiarm_binding` and
`myeloid_il6` all resolve. A repo-wide `find` returns nothing.

Consequences:
- **All 9 `ENG` entries have `costim_arm = None`**, so the branch never fires and every clinical validation run is
  unaffected. Nothing already reported is wrong because of this.
- **The moment anyone adds a costim arm to `ENG` — the entire purpose of a "costim engager counterscreen" — the run
  dies with `ModuleNotFoundError` inside `attach_pd`, at build time, for the first organ.**
- The doc's own earlier draft listed `signaling_dynamics` as part of the "live import graph (verified this task)".
  It is not live, and it is not present. That claim has been corrected.

**Recommendation:** either restore/write `signaling_dynamics.py`, or guard the import
(`try: import signaling_dynamics ... except ImportError: raise RuntimeError("costim_arm requested but the
signaling engine is missing")`) so the failure names its cause. Until then, **T9 cannot run a costim construct**, and
no costim result may be claimed.

### 5.11 Open questions a committee will ask
0. **Where is `signaling_dynamics.py`?** (§5.10 — blocks every costim run.)
1. **What is human IL-6 clearance?** The one in use is a *modelling* paper's fitted estimate (§3.6). Until a measured
   value is found (or the IL-6-through-PBPK path is un-bypassed), the engine can rank molecules but cannot predict
   pg/mL. Is there a published IL-6 IV PK study at all?
2. **What is the liver Kupffer census?** Its absence makes every IL-6 number a lower bound.
3. **Does the ranking survive `TSIM_DAYS=24` with like-for-like (same-statistic) anchors?** With only two valid anchors
   left (mosunetuzumab 152, teclistamab 21 — both population means, different targets, different routes), the
   validation set is **n = 2**. A Spearman on n = 2 is not a validation.
4. **Is `pd_every = 3` converged?** Unknown.
5. **What is the SC bioavailability of each molecule?** One shared `F_sc = 0.6` is applied to four different molecules.
6. **Why do the PK and PD arms carry different KDs (§5.6)?**

---

## Appendix — RUN-verifications performed for this document (Python, conda `claude-skills`)

| Claim | Verification |
|---|---|
| 1 mg of 146 kDa = 6.849 nmol → 2.209 nM in `V_pl = 3.1 L` | `1/146*1e3 = 6.8493`; `/3.1 = 2.2095` |
| mosunetuzumab 60 mg → 410.96 nmol → 132.6 nM | `60/146*1e3/3.1` |
| `k_cat` mosun (fFcRn 0.89) = 0.03853/day, t½ 18.0 d; elran (0.70) = 0.10509/day, t½ 6.6 d | `0.3503*(1−f)`, `ln2/k` |
| kon unit conversion 1e5 /M/s = **8.64** /nM/day (⇒ the 8.64e-3 fallback = 100 /M/s) | `1e5/1e9*86400` |
| teclistamab kon 1.2847e6 /M/s = 111.0 /nM/day; koff 2.315e-4/s = 20.0 /day | `x/1e9*86400`, `x*86400` |
| IL-6 distribution volume 11.65 L = Σ Vis (8.5508) + V_pl (3.10) | `qsp_costim_window_v2._PBPKArrays().Vis.sum()`, `q.V_PLASMA` |
| Organ-set interstitial coverage 4.678 / 8.551 L = **54.7%** | summed `Vis` over `ORGANS` |
| Census values are **scales**: spleen 290,206 × 58,579 sampled myeloid = **1.70e10** physiological; bone 33,633,218 × 1,445 = **4.86e10**; lung 1,083,826 × 26,388 = **2.86e10**; small_int 84,910 × 12,366 = **1.05e9** | myeloid-token count over every `agents_<organ>.npz` + the census JSON |
| Per-cell secretion at the live `S_MAX = 10.6 molec/s` = **0.00133 pg/hr/cell** (NOT the 0.0196 the comments claim) | `10.6*21000/6.022e23*1e12*3600` |
| `tsim = 49 d, dt = 0.02` → 2451 transport steps, 817 PD steps | `int(49/0.02)+1`, `//3` |
| Myeloid agent counts per compartment (used for the census table) | `agents_*.npz` label scan with `myeloid_il6.MYELOID_TOKENS` |
| `CYTO_IL6_CLINICAL_ANCHOR_PGML`, `cytokine_to_pgml` have **no consumer at all**; `IL6_SCALE` has **no consumer in the live import graph** (but one dead-module caller still reads it) | grep over `engine/*.py` |

### Appendix B — corrections made by the adversarial re-verification (2026-07-13)

| # | Claim in the previous draft | Verdict | Correction |
|---|---|---|---|
| 1 | `wholebody_pd` → … → `signaling_dynamics` listed as a **verified live import** | **FALSE** | The module **does not exist in the repo** (`find_spec` → `None`). Removed from the graph; escalated to **§5.10** as the defect that blocks every costim run. |
| 2 | Live schedules live in `engine/regimen_schedules_final.json` | **WRONG FILE** | The driver loads `{KWS}/handoff/regimen_schedules_final.json` (`:44`). The `engine/` copy is a never-loaded, non-identical duplicate (handoff is a strict superset). The 6 tabulated schedules agree in both, so **no schedule value changed**. |
| 3 | `or 1.0` census leak is "~5 orders of magnitude" below the scaled organs | **WRONG VALUE** | RUN-verified: 28,811 excluded sampled agents vs 9.53e10 scaled physiological myeloid = **3.0e-7 → ~6–7 orders**. (Conclusion — negligible — unchanged.) |
| 4 | `kdeg=0.20/hr` is `[FITTED: to land the plasma peak]` | **OVERCLAIM** | Nothing in this engine tunes it. Correct tag: **`[FITTED-ELSEWHERE / NOT MEASURED]`** — a fitted estimate from a semi-mechanistic *modelling* paper (Chen 2019, per the repo's own audit). Conclusion (absolute pg/mL not defensible; ranking is) unchanged. |
| 5 | `[MEASURED]` on the myeloid census, `k_hit`, and mosunetuzumab `il6_obs=152` | **MIS-TAGGED** | None has a source **in the code**. Retagged `[LITERATURE-CLAIMED — NO SOURCE IN CODE]` / `[UNSOURCED — TBD]`. (The census's external attribution, **PMID 37871201**, was *confirmed* to exist at `reference_unified_binding/IL6_MECHANISM_PROVENANCE.md:86`.) |
| 6 | `sigL`, `k_dist` at `coupled_percell_pk.py:65`; IL-6 exact step at `myeloid_il6.py:212–214` | **OFF-BY-ONE** | `:64` and `:213–214` respectively. |
| 7 | `IL6_SCALE` "consumed by nothing" | **OVERCLAIM** | One out-of-graph module still reads it; deleting the assignment is a breaking change for that caller. |
| 8 | File lengths 396 / 275 lines | **WRONG** | 395 / 274. |
| 9 | Elranatamab kon/koff "MEASURED from published QSP" | **OVERSTATED** | kon measured; **koff is DERIVED = kon·KD** per the file's own `measured_vs_derived` string. |
| — | **All PMIDs cited in this doc** (20376398, 29808005, 30571054, 31268236, 37533643, 37871201, 38831634, 3944542, 9400735) | **NO FABRICATIONS** | Each traces either to the two live files or (37871201 only) to `reference_unified_binding/IL6_MECHANISM_PROVENANCE.md:86` at the exact file:line cited. |
| — | Vis mismatch table, 54.7% coverage, agent counts, census scales, regimen values, KD asymmetry, kon/koff conversions, k_cat/t½, the `or 1.0` defect, the 6-day motility cadence, the 0.0196 stale comment, 71 pL synapse volume, nstep 2451/817, every `file:line` in the two live files | **CONFIRMED** | Re-verified independently; no change. |
