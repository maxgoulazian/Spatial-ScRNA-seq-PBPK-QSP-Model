---
title: "T2 — Whole-body per-cell PK (vascular QSS → BEC extravasation → per-cell interstitial graph → LEC drainage → lymph → plasma)"
subsystem: T2
model: costim_engager_counterscreen
live_source_files:
  - engine/coupled_percell_pk.py   (159 lines, read in full)
  - engine/wholebody_percell.py    (209 lines, read in full)
date: 2026-07-13
revised: 2026-07-13 (adversarial verification pass — see §6 CORRECTION LOG)
generated_by: workflow-subagent T2
provenance_rule: "Every number below was READ from the two live source files (or from the live module that injects it, cited by file:line). Nothing is inferred from memory. Where the code gives no source, the row says so."
---

# T2 — Whole-Body Per-Cell PK

> **Scope discipline.** This doc documents ONLY `engine/coupled_percell_pk.py` and `engine/wholebody_percell.py`.
> Where a value is *injected* by a live caller (`run_tce_pd_reval.py`, `qsp_costim_window_v2.py`) I cite that
> file:line and mark the row **consumed, not owned**. The dead modules (`cytokine_pbpk.py`, `il6_pbpk.py`,
> `unified_binding.py`, `multiarm_kinetic.py`, `biexact_solver.py`, `rna_to_receptor.py`, `convert_copies_ALL.py`,
> `calib_kdeath.py`) are **not** cited anywhere below, even though two of them contain a near-verbatim copy of
> the quasi-steady vascular algebra (EQ-9). They are not in the execution path and are not evidence for anything.

---

## 1. PURPOSE & DATAFLOW POSITION

### 1.1 What this subsystem is

T2 is the **transport core**: it carries an antibody molecule from the syringe to the surface of a named cell in a
named organ, and back out again. It is the only place in the model where "plasma concentration" and "the drug
concentration this particular Xenium barcode sees" are connected. Everything downstream (synapse formation,
killing, IL-6) reads `TissueGraph.C[i]` — the free interstitial drug concentration at cell *i* — and everything
upstream (dose, route, schedule) writes into `A_pl` / `A_sc`.

The physiological chain, exactly as the code implements it (`wholebody_percell.py:4-11`):

```
        dose (IV bolus | IV infusion | SC depot)
                 │
          ┌──────▼──────┐   k_cat = CLup·(1−fFcRn) + k_renal      (pk:76)
          │  A_pl       │◄──────────── k_lymph_return · A_ly       (pk:153)
          │  plasma     │──── catabolism ─────────────────────────► ✝
          └──────┬──────┘
                 │ Q (organ blood flow)          ── PER ORGAN ──
        ┌────────▼─────────┐
        │ organ VASCULAR   │  ***QUASI-STEADY*** — no state variable  (wbp:163-168)
        │ compartment Cvasc│
        └────────┬─────────┘
                 │ PS_ex = k_dist·L·(1−σ_V)     extravasation ONLY at BEC cells   (pk:127, wbp:169,196)
        ┌────────▼──────────────────────────────────────────────┐
        │ per-cell INTERSTITIAL GRAPH  C[i]  (nM, one node/cell)│
        │  · kNN Laplacian diffusion, ECM-hindered   (wbp:99-109)│
        │  · per-cell ANTIGEN SINK (TMDD)            (wbp:172-192)│
        └────────┬──────────────────────────────────────────────┘
                 │ PS_ly = k_dist·L·(1−σ_L), spatially graded toward LEC cells  (wbp:201-205)
          ┌──────▼──────┐
          │  A_ly lymph │──── k_lymph_return ──► back to plasma     (pk:154)
          └─────────────┘
```

### 1.2 Feeds in

| Input | Where from | Notes |
|---|---|---|
| `organs`, `agents_dir`, per-cell `x,y`, cell-type labels, per-cell target copies | `load_tissue()` (`coupled_percell_pk.py:11-60`) reading the Xenium/Tabula-Sapiens `.npz` agent files | tumor npz uses `dens_<gene>`; heme npz uses `<gene>_copies`; normal organs use the bare gene column (pk:53-59) |
| `bec_lec` (per-organ BEC / LEC cell indices) | caller-supplied JSON (`run_tce_pd_reval.py:46`, built from LYVE1/PROX1/PDPN markers per `wholebody_percell.py:9-10`) | remapped through the finite-coordinate mask at pk:98-100 |
| `Q, L, sigV, Vis, Vv` | `qsp_costim_window_v2.py:126-143` frozen physiology → `run_tce_pd_reval.py:51` | **consumed, not owned** (see §3.4) |
| `pools` (per-organ physiological antigen pool, nM) | `run_tce_pd_reval.py:130` | sets `count_scale` (EQ-1) |
| `KD_nM, n_arm, kint, mw, fFcRn, kon1, koff1, kdeg, span_coeng_nm` | per-molecule config, `run_tce_pd_reval.py:142-145` | the ONLY per-molecule inputs (module docstring, pk:3-5) |

### 1.3 Feeds out

| Output | Consumer |
|---|---|
| `TissueGraph.C[i]` — free interstitial drug at every cell | the PD layer (`coupled_percell_pd.py`) → synapse/kill/IL-6 |
| `S` (bound copies/cell), `D` (avidity-bound) | PD readout; `bound_nM` record (pk:137) |
| `organ_sink` (nmol/day internalized) | **emergent TMDD** — recorded per organ by `simulate()` (pk:136). ⚠ **The production PD driver discards it** (`coupled_percell_pd.py:264-265`; no `sink_rec` in that file) — see EQ-22. |
| `Cplasma_ugml(t)` | the clinical PK observable (pk:157) |
| `self.snaps[t][organ]` | spatial overlay figures (pk:139-146) |

### 1.4 Live-path caveat you must know before reading §2

`CoupledPerCellPK.simulate()` (pk:108-159) is the **PK-only** driver. The production runner calls
`CoupledPerCellPD.simulate_pd()` (`coupled_percell_pd.py:199`), which **re-implements the same loop** with PD
attached — the systemic ODE at `coupled_percell_pd.py:357-358` is line-for-line the same as pk:153-154 plus two
extra sink terms (`_heme_sink_nmol_day`, `_blood_sink_nmol_day`). The per-organ physics (`TissueGraph.step`) is
**shared, not duplicated** (`coupled_percell_pd.py:263-264` calls the same `g.step`). So:
- `wholebody_percell.py` is 100% live.
- `coupled_percell_pk.py`'s `__init__` + `load_tissue` are 100% live (PD inherits them).
- `coupled_percell_pk.py`'s `simulate()` body is the **reference/transport-only** driver; its equations are live
  *as physics* (they are mirrored exactly in the PD driver) but the PD runs execute the copy, not this function.
  I flag this rather than pretend the whole file is on the hot path.

---

## 2. GOVERNING EQUATIONS

Units convention throughout: **time = days**, **concentration = nM**, **amount = nmol**, **volume = L**,
**distance = µm** (graph) or **nm** (antibody span), **receptor = copies/cell**.
`AVO = 6.02214076e23` (`wholebody_percell.py:29`).

---

### EQ-1 — Physiological antigen pool → per-cell copies + `count_scale` (coupled_percell_pk.py:91-96)

```
tot_copies = pool_nM · Vis[o] · AVO / 1e9          # copies in the whole organ interstitium
R_percell  = ag                                     # REAL per-cell copies from the npz — NOT redistributed
count_scale = tot_copies / Σ_i ag_i                 # sampled cells → physiological cell population
```
- **(a) Math.** `tot_copies` [copies] = (nM = 1e-9 mol/L) × (L) × (copies/mol) / 1e9 — the `1e9` divides out the
  nano prefix. `count_scale` [dimensionless] = physiological receptor pool ÷ receptor pool actually present in the
  sampled ABM cells.
- **(b) Biology.** The spatial atlas is a *subsample* of a real organ (a Xenium section is not 30 L of muscle).
  Each sampled cell therefore stands for `count_scale` real cells. The code deliberately **keeps the measured
  per-cell copy number** (so avidity, which is nonlinear in receptor density, is computed at the physiological
  receptor number) and scales the **cell count** instead.
- **(c) Rationale / rejected alternative.** The comment at pk:85-90 states the rejected alternative explicitly:
  the *old* path redistributed the organ pool evenly over sampled cells ("pool-redistribution"). That is
  numerically identical for a monovalent binder (sink linear in R) but **wrong for a bivalent binder**, because
  the avidity term is ∝ R² (EQ-17) — flattening the receptor distribution destroys the avidity heterogeneity.
  This is the single most important structural choice in the antigen sink.
- **(d) Units.** `pool_nM` [nM]; `Vis[o]` [L]; `ag_i` [copies/cell]; `count_scale` [cells/cell, dimensionless].
- **Guard:** if `pool_nM ≤ 0` or `Σag = 0` → `R=0, count_scale=0` (pk:95-96) → the organ is a pure
  transport compartment with no sink.

---

### EQ-2 — Per-molecule kinetic-rate unit conversion (coupled_percell_pk.py:69-71)

```
_SPD  = 86400                                       # s/day
kon1  [1/nM/day] = kon1[1/M/s] / 1e9 · 86400
koff1 [1/day]    = koff1[1/s]  · 86400
```
- **Biology.** SPR reports on-rates in M⁻¹s⁻¹; the whole-body engine runs in nM and days. Nothing physical
  happens here — but a wrong 1e9 here is a 9-order-of-magnitude affinity error, so it is documented.
- **Rationale.** `None` → the kinetic branch is disabled and `TissueGraph.step` falls back to the QSS binding
  (EQ-17). This is the switch between "measured kon/koff exist for this molecule" and "only KD is known".

---

### EQ-3 — Systemic elimination: FcRn salvage as a lumped catabolic rate (coupled_percell_pk.py:76)

```
k_cat = CLup · (1 − fFcRn) + k_renal                [1/day]
```
- **(a) Math.** First-order plasma elimination rate constant. Default `CLup = 0.3503 /day`, `k_renal = 0.0 /day`
  (pk:65); `fFcRn` is per-molecule (0.89 for the IgG1 engagers, 0.70 for IgG2 elranatamab —
  `run_tce_pd_reval.py:77-81`).
- **(b) Biology.** The endothelium pinocytoses plasma at a rate `CLup`. Inside the acidified endosome, FcRn binds
  the Fc and **recycles the fraction `fFcRn` back to the surface**; the remaining `(1−fFcRn)` is delivered to the
  lysosome and catabolised. FcRn salvage is the reason an IgG lives ~3 weeks instead of ~1 day.
- **(c) Rationale / rejected alternative.** The **rejected** alternative is an explicit endosomal compartment
  (pinocytosis → endosome with FcRn binding at pH 6 → recycling/degradation split), which is what a full
  Shah–Betts PBPK does. This model **collapses that to one algebraic product**. Consequences a reviewer will
  probe: (i) FcRn is not saturable here, so no IVIG-style competition and no nonlinearity at high dose; (ii) there
  is no membrane-bound FcRn pool, so no dose-dependent half-life; (iii) the salvage term is applied only to the
  *plasma* pool, not to drug that has extravasated. In exchange, elimination costs one multiplication and cannot
  be stiff.
  **Numerically (computed this session from the code's own defaults):** `k_cat = 0.3503 × (1−0.89) = 0.03853/day`
  → catabolic t½ = ln2/k_cat = **17.99 d**; at `fFcRn = 0.70` (IgG2) → `k_cat = 0.1051/day` → t½ = **6.60 d**.
  So Fc-isotype is the *only* lever that moves systemic half-life in this engine.
- **(d) Units.** `CLup` [1/day]; `fFcRn` [dimensionless, 0–1]; `k_renal` [1/day]; `k_cat` [1/day].
- **⚠ Honest flag.** `k_renal` defaults to **0.0** in this module (pk:65) and the live runner never overrides it
  (`run_tce_pd_reval.py:142-145` passes no `k_renal`). The size-gated glomerular-sieving law
  (`k_renal = k_renal_max / (1 + (mw/mw50)^hill)`, `qsp_costim_window_v2.py:216-218`) exists in the *lumped* QSP
  model but is **NOT connected to this per-cell PK engine**. Therefore a no-Fc, sub-glomerular construct (BiTE)
  simulated through T2 would get FcRn loss but **no renal filtration** — its half-life would be wrong (too long).
  All current live engagers are Fc-bearing IgGs (≈146–195 kDa), so this is latent, not active. It is still a real
  gap and I am not going to hide it.

---

### EQ-4 — Dose administration: IV bolus vs SC depot (coupled_percell_pk.py:117-120)

```
while sched[si].t ≤ t + 1e-9:
    n_dose[nmol] = mg / mw_kda · 1e3
    route == "SC":  A_sc += n_dose        # subcutaneous depot
    else:           A_pl += n_dose        # IV bolus straight into plasma
```
- **(a) Math.** mg ÷ kDa × 1e3 = nmol. Checked against a real molecule in the live table: mosunetuzumab
  `mw = 146.0` kDa (`run_tce_pd_reval.py:77`) ⇒ 1 mg = 1/146.0 × 1e3 = **6.849 nmol**. The live `mw` range is
  146.0–195.0 kDa (`run_tce_pd_reval.py:77-88`).
- **(b) Biology.** An IV bolus is instantaneously mixed into the central plasma pool; an SC injection forms a
  depot in the subcutis from which drug must first be absorbed (EQ-5).
- **(c) Rationale.** The bolus idealisation (instantaneous mixing) is standard for mAbs because mixing (~1 min) is
  orders of magnitude faster than the ~18-day elimination.
- **⚠ Live-path caveat.** This bolus branch is the behaviour of `coupled_percell_pk.simulate()`. The **production
  PD driver does not run an IV bolus**: for every IV row the runner passes `iv_inf_h = 2.0`
  (`run_tce_pd_reval.py:176` — "real IV infusion duration (clinical TCE ~2h), not bolus") into
  `CoupledPerCellPD.simulate_pd()` (`coupled_percell_pd.py:199`). So the bolus documented here is the
  *reference-driver* behaviour, not what the live TCE runs execute. IV **infusion** in *this* module is EQ-6.
- **(d) Units.** `mg` [mg]; `mw` [kDa]; `A_pl`, `A_sc` [nmol].

---

### EQ-5 — SC absorption: first-order depot with bioavailability (coupled_percell_pk.py:123, 153)

```
J_sc  = ka_sc · A_sc                       [nmol/day]     (pk:123)
A_sc ← A_sc − dt · J_sc                                    (pk:123, explicit Euler)
dA_pl += F_sc · J_sc                                       (pk:153)
```
- **(a) Math.** The depot empties with rate constant `ka_sc` = **0.25 /day** (pk:108) → absorption t½ =
  ln2/0.25 = **2.77 d** (computed). Only the fraction `F_sc` = **0.6** (pk:108) reaches plasma; the balance
  `(1−F_sc)·J_sc` leaves the depot and **is not tracked anywhere** — it is destroyed.
- **(b) Biology.** SC-injected IgG is absorbed mainly by convective uptake into the lymphatics of the subcutis;
  transit through the draining lymph node and pre-systemic proteolysis/FcRn-independent catabolism cause an
  absolute bioavailability well below 1. The slow absorption is what gives SC engagers (teclistamab,
  epcoritamab, talquetamab, elranatamab in `run_tce_pd_reval.py:79-82`) their **flat-topped, blunted Cmax** — which
  is precisely the mechanism claimed to reduce CRS.
- **(c) Rationale / rejected alternative.** The rejected alternative is routing the SC depot through the *model's
  own lymph pool* (`A_ly`), which would make bioavailability emergent from `k_lymph_return` and lymph-node
  catabolism instead of imposed as `F_sc`. The code chose the imposed-`F_sc` shortcut. This is an **imposed, not
  emergent** feature (see §4) and `F_sc = 0.6` carries **no citation anywhere in the code**.
- **(d) Units.** `ka_sc` [1/day]; `F_sc` [dimensionless]; `A_sc` [nmol]; `J_sc` [nmol/day].
- **⚠ Note on the loss term.** Because the depot is debited the **full** `J_sc` but plasma is credited only
  `F_sc·J_sc`, mass is *deliberately* not conserved across the SC interface. That is correct pharmacokinetics
  (the lost 40% really is destroyed pre-systemically) but it means the model cannot report where it went.

---

### EQ-6 — IV infusion (coupled_percell_pk.py:124, 153)

```
infn = (inf_rate / mw · 1e3)   if (inf_rate > 0 and t ≤ inf_dur) else 0     [nmol/day]
dA_pl += infn
```
- **Biology.** Zero-order infusion into plasma for a finite duration. `inf_rate` is in mg/day; the same
  mg→nmol conversion as EQ-4. Defaults `inf_rate = 0, inf_dur = 0` (pk:108) → OFF.
- **Rationale.** Rectangular (not ramped) infusion: the clinical step-up dosing of TCEs is delivered over a fixed
  window at a fixed rate, so a rectangle is exact, not an approximation.

---

### EQ-7 — Plasma concentration and the reported observable (coupled_percell_pk.py:122, 157)

```
C_pl [nM]      = A_pl / V_pl                                   (pk:122)
C_rec [µg/mL]  = (A_pl / V_pl) · mw_kda / 1e3                  (pk:157)
```
- **(a) Math.** `V_pl` = 3.1 L (pk:65). Checked on a live molecule: 1 nM of the 146.0-kDa mosunetuzumab
  (`run_tce_pd_reval.py:77`) = 1e-9 mol/L × 146 000 g/mol = 0.1460 mg/L = **0.1460 µg/mL** — which is exactly
  `1 × 146.0/1e3`. The conversion is checked and correct. (Per-molecule: glofitamab `mw = 195.0` ⇒ 0.1950 µg/mL
  per nM, `run_tce_pd_reval.py:78`.)
- **(b) Biology.** `C_pl` is what every organ's vascular compartment sees (EQ-9), and `C_rec` is the only quantity
  in this subsystem that can be compared to a clinical label.
- **(c) Rationale.** The central volume is set to the **physiological plasma volume** (3.1 L), *not* to a fitted
  `Vc`. That is a real commitment: in a classical 2-compartment PopPK fit, `Vc` is a free parameter that absorbs
  early-distribution error. Here it cannot — so any Vss/Vc discrepancy is forced into `σ_L`/`k_dist` (EQ-8),
  which is where those two knobs got their calibrated values.

---

### EQ-8 — Two-pore permeability–surface products (coupled_percell_pk.py:127)

```
PS_ex[o] = k_dist · L[o] · (1 − σ_V[o])        [L/day]    vascular → interstitium
PS_ly[o] = k_dist · L[o] · (1 − σ_L)           [L/day]    interstitium → lymph
```
- **(a) Math.** `L[o]` = organ **lymph flow** = `Q[o]/500` (`_LYMPH_RATIO`, `qsp_costim_window_v2.py:101`; applied
  at `qsp_costim_window_v2.py:141`); `σ_V[o]` = organ vascular reflection coefficient (0.75 tumor … 0.99 brain,
  table rows `qsp_costim_window_v2.py:84-98`); `σ_L` = 0.85 lymphatic reflection (pk:64); `k_dist` = 3.0 (pk:64).
- **(b) Biology.** This is the **convection-dominated (two-pore) picture** of macromolecular transport: an IgG
  does not diffuse across the endothelium, it is *dragged* through large pores by the fluid flux, and the pore
  partially rejects it. `(1 − σ)` is the fraction of the fluid flux that carries the antibody. Leaky tumor
  endothelium ⇒ low σ_V (0.75) ⇒ 5× the extravasation of brain (σ_V = 0.99) at the same lymph flow.
- **(c) Rationale / rejected alternative.** The rejected alternative is a diffusive (Fickian) `PS` term. It is
  rejected because for a 150-kDa protein the Péclet number across the endothelium is ≫1 — convection wins. The
  **steady-state interstitial:plasma concentration ratio** of this pair is `(1−σ_V)/(1−σ_L)`, which is
  **independent of `k_dist`** — that is the entire reason the code splits distribution into an *extent* knob
  (σ_L, sets Vss) and a *rate* knob (k_dist, sets the α-phase). See §3.2.
  Computed at the code's values: lung `(1−0.95)/(1−0.85) = 0.333`, tumor `(1−0.75)/(1−0.85) = 1.667` — i.e. at
  equilibrium the tumor interstitium is *above* plasma concentration, the classic leaky-vessel result.
- **(d) Units.** `L` [L/day]; σ [dimensionless]; `k_dist` [dimensionless]; `PS` [L/day].

---

### EQ-9 — ★ THE VASCULAR COMPARTMENT IS QUASI-STEADY ★ (wholebody_percell.py:163-168)

```
Cis_bec = mean( C[i] : i ∈ BEC )                                       (wbp:167)

0 = Q·(C_pl − C_vasc) − PS_ex·(C_vasc − Cis_bec)          ← the QSS balance (wbp:165)

⇒  C_vasc = ( Q·C_pl + PS_ex·Cis_bec ) / ( Q + PS_ex )                 (wbp:168)
```
- **(a) Math.** A flow-limited well-stirred organ vascular pool: blood enters at plasma concentration `C_pl` and
  leaves at `C_vasc` (net exchange `Q·(C_pl − C_vasc)`), while `PS_ex·(C_vasc − Cis_bec)` leaks into the tissue.
  Setting the time derivative to zero and solving gives a **flow-weighted average of plasma and BEC-local
  interstitium**, weighted by the two conductances. `max(Q+PS_ex, 1e-12)` guards a zero-flow organ.
- **(b) Biology.** Blood in an organ's capillary bed is exchanged with the systemic circulation on the timescale
  of a **single pass** (seconds). Antibody leaking out of that capillary bed happens on the timescale of **hours
  to days**. The vascular pool therefore has no memory: whatever the plasma is doing, the capillary content has
  already caught up. `C_vasc` is a slaved variable, not a state.
- **(c) MECHANISTIC RATIONALE — why not integrate it explicitly (this is the question the committee will ask).**
  The vascular compartment's own relaxation rate is `(Q + PS_ex)/Vv`. Computed **this session from the live
  physiology** (table `qsp_costim_window_v2.py:82-101`, built at `:127-143`; `Q = Qfrac/Σ Qfrac_nonlung × 5000`
  at `:139-140`, `L = Q/500` at `:141`, `Vv = fV·V` at `:142`; `k_dist=3, σ_L=0.85`):

  | organ | Q (L/day) | Vv (L) | **Q/Vv (1/day)** | PS_ex (L/day) | Q/PS_ex |
  |---|---|---|---|---|---|
  | lung | 5000.0 | 0.0525 | **95 238** | 1.500 | 3333 |
  | kidney | 1005.3 | 0.0325 | **30 885** | 0.603 | 1667 |
  | small_int | 529.1 | 0.0247 | **21 421** | 0.318 | 1667 |
  | tumor | 105.8 | 0.0070 | **15 117** | 0.159 | **667** |
  | liver | 343.9 | 0.2070 | 1 661 | 0.310 | 1111 |
  | bone | 264.6 | 0.4100 | 645 | 0.079 | 3333 |

  - **Separation of timescales is 3–4 orders of magnitude:** since `PS_ex = k_dist·(Q/500)·(1−σ_V)`, the ratio is
    purely a function of σ_V: `Q/PS_ex = 500/(3·(1−σ_V)) = 166.7/(1−σ_V)`. Over the 15 tabulated organs it ranges
    from **667 (tumor, σ_V = 0.75)** to **16 667 (brain, σ_V = 0.99)**. The premise "Q ≫ PS" is therefore satisfied
    *everywhere*, worst case by a factor of 667. The code comment's claim `Q/Vv ~ 1e5/day` (wbp:165) is
    **verified** (lung: 95 238/day).
    *Scope note:* brain is in the physiology table but **not** in the live organ set —
    `ORGANS = [spleen, bone, large_int, liver, lung, small_int, pancreas, kidney, skin, heart, adipose]`
    (+`tumor` for solid runs), `run_tce_pd_reval.py:48, 129`. The worst case (tumor, 667) **is** live; the 16 667
    max is a table statement.
  - **Explicit integration would be catastrophically stiff.** Forward-Euler stability on the vascular pool needs
    `dt < 2·Vv/(Q + PS_ex)`. The tightest organ (lung, and lung *is* live) gives `dt < 2.10e-5 day ≈ 1.8 s`. The
    production step is `dt = 0.02 day` (`run_tce_pd_reval.py:143`) — i.e. **953× larger than the stability limit**
    (computed). An explicit vascular ODE at the production dt would not be inaccurate, it would **blow up**.
    Keeping it would force either dt ≈ 2e-5 d (a 953× increase in cost, over the ~400k-cell graphs the code itself
    sizes at wbp:58 × 11–12 organs × 21+ days — completely infeasible) or an implicit/stiff solver for a variable
    that carries no information.
  - **The approximation costs almost nothing.** The QSS error is O(dt·(Q+PS)/Vv)⁻¹ relative — i.e. the vascular
    pool is wrong only during the first ~1e-5 day after a bolus, which is below the recording resolution.
  - **Rejected alternatives:** (i) explicit vascular ODE — stiff, see above; (ii) implicit/stiff integrator on the
    full coupled system — would have to include the 400k-node interstitial graph and the binding ODEs, destroying
    the operator-split structure that makes the per-cell graph tractable; (iii) dropping the vascular compartment
    entirely and extravasating straight from plasma (`J = PS·(C_pl − C_is)`) — this is what the *well-mixed*
    fallback organs do (EQ-21), and it is exactly the `Q → ∞` limit of EQ-10. Keeping the finite-Q form retains
    **flow-limitation**, which matters for the low-`Q/PS_ex` organs (tumor, 667) and is free.
- **(d) Units.** `Q` [L/day]; `PS_ex` [L/day]; `C_pl`, `C_vasc`, `Cis_bec` [nM]; `Vv` [L, **and see the flag below**].
- **⚠ Consequence: `Vv` is a dead parameter inside `TissueGraph`.** It is passed in (`wholebody_percell.py:76`) and
  stored (`:79`) but appears **nowhere** in `step()`. The QSS eliminates it by construction. It is retained only so
  the caller's physiology tuple stays intact. A reviewer noticing "you list a vascular volume but never use it" is
  right — the honest statement is: *the vascular volume sets only the relaxation rate we deliberately took to
  infinity; it cannot affect any slow observable.*

---

### EQ-10 — BEC extravasation flux, and its series-conductance identity (wholebody_percell.py:169)

```
J_extrav = PS_ex · ( C_vasc − Cis_bec )                     [nmol/day]
```
Substituting EQ-9:
```
J_extrav = [ PS_ex·Q / (Q + PS_ex) ] · ( C_pl − Cis_bec )   ≡  G_ex · (C_pl − Cis_bec)
```
- **(a) Math.** The extravasation conductance is the **series (harmonic) combination** of blood delivery `Q` and
  endothelial permeability `PS_ex` — exactly Ohm's law for two resistors in series. In the permeability-limited
  regime (`PS_ex ≪ Q`, which holds everywhere, §EQ-9 table) `G_ex → PS_ex` and the vascular compartment becomes
  transparent; in a hypothetical flow-limited organ (`Q ≪ PS_ex`) it would clamp to `Q`.
- **(b) Biology.** Drug crosses the endothelium **only at blood-endothelial cells** — the model does not let
  antibody appear in the middle of a tumor nest by fiat, it must enter at a vessel and diffuse. This is the
  structural reason the model can produce a *penetration gradient* at all.
- **(c) Rationale.** `Cis_bec` is the **mean over BEC nodes**, not the organ-mean interstitium. That is the correct
  local driving force (the endothelium sees its own neighbourhood, not the organ average) — and it is what makes
  the extravasation flux *fall* when the perivascular region saturates, i.e. **back-pressure**. Note the sign is
  free: if `Cis_bec > C_vasc`, `J_extrav < 0` and drug **reabsorbs** into blood. The code allows this
  (no clamp) — correct, and necessary for the terminal phase.
- **(d) Units.** `PS_ex` [L/day]; concentrations [nM]; `J_extrav` [nmol/day].
- **⚠ Asymmetry worth flagging:** extravasation is driven by the **BEC-local** mean `Cis_bec` (wbp:167) but
  lymph drainage (EQ-12) is driven by the **whole-tissue** mean `C.mean()` (wbp:201). The two fluxes therefore use
  different definitions of "the interstitial concentration". Defensible (they happen at different places) but it is
  an inconsistency a reviewer can pick at.

---

### EQ-11 — BEC source field: distributing the flux onto the entry cells (wholebody_percell.py:193-196)

```
vcell = Vis / n                                              # pericellular volume per sampled cell (wbp:171)
src[i∈BEC] += ( J_extrav / Vis ) · ( n / n_bec )             [nM/day]
```
- **(a) Math.** Mass check (verified algebraically): `Σ_{i∈BEC} src_i · vcell = n_bec · (J_extrav/Vis)·(n/n_bec) ·
  (Vis/n) = J_extrav`. **Exactly mass-conserving** — the `(n/n_bec)` factor is precisely what concentrates an
  organ-average concentration rate onto the BEC subset without creating or destroying drug.
- **(b) Biology.** The extravasated antibody appears in the interstitial fluid immediately abluminal to the
  vessel wall, i.e. in the pericellular volume of the BEC barcodes. From there it must **diffuse** (EQ-13) to reach
  a tumor cell three cell-diameters away.
- **(c) Rationale / rejected alternative.** The rejected alternative — spreading `J_extrav` uniformly over all
  cells — is exactly what a compartmental PBPK does, and it **eliminates the penetration gradient**, which is the
  entire scientific point of a per-cell model. Concentrating the source at BEC nodes is what makes the
  perivascular-vs-core exposure difference emergent rather than imposed.
- **(d) Units.** `J_extrav` [nmol/day]; `Vis` [L]; `src` [nM/day]; `vcell` [L].

---

### EQ-12 — Lymphatic drainage: organ-exact flux, spatially graded field (wholebody_percell.py:197-205)

```
Cis_mean       = mean(C)                                                    (wbp:201)
drain_to_lymph = PS_ly · Cis_mean                             [nmol/day]     (wbp:202)
w_i            = drain_w_i · C_i                                            (wbp:203)
wint           = mean(w) · Vis                                [nM·L]        (wbp:204)
drain_field_i  = w_i · ( drain_to_lymph / wint )              [nM/day]      (wbp:205)
```
- **(a) Math.** Mass check (verified algebraically): `Σ_i drain_field_i · vcell = (drain_to_lymph/wint) ·
  Σ_i w_i · (Vis/n) = (drain_to_lymph/wint) · mean(w) · Vis = drain_to_lymph`. The normalisation makes the
  *spatial* field integrate **exactly** to the *compartmental* 2-pore flux. Guard: `wint > 1e-18` else zero field.
- **(b) Biology.** Interstitial fluid is continuously produced by capillary filtration and removed by the initial
  lymphatics; the antibody is swept along by that convective flow (`(1−σ_L)` of it, EQ-8). Because initial
  lymphatics are *distributed* through the tissue, the whole tissue drains — but faster near a lymphatic.
- **(c) Rationale / rejected alternative — this one is subtle and the code comment says so (wbp:197-200).** The
  "obvious" implementation is to remove drug **only at LEC nodes**. That was rejected because it makes lymphatic
  exit **diffusion-limited**: drug far from an LEC cannot leave, so it accumulates, and the organ's terminal
  washout becomes an artifact of the LEC sampling density in the Xenium section. The chosen construction keeps
  the *total* flux exactly at the 2-pore compartmental value (so whole-body PK is unchanged and mass-correct)
  while imposing a **physiological spatial gradient** on where it leaves. This decouples "how much leaves" (physiology,
  imposed) from "where it leaves" (geometry, emergent).
- **(d) Units.** `PS_ly` [L/day]; `Cis_mean`,`C_i` [nM]; `drain_w_i` [dimensionless]; `wint` [nM·L];
  `drain_field` [nM/day].

---

### EQ-13 — Interstitial diffusion: ECM-hindered kNN graph Laplacian (wholebody_percell.py:81, 98-109)

```
φ_D,i   = 1 / ( 1 + α_D · ecm_i )                                          (wbp:81)
D_i     = D_um2s · φ_D,i                                    [µm²/s]        (wbp:101)
d_ij    = max( |x_i − x_j| , 0.003 )                        [µm]           (wbp:104)
g_ij    = ½·(D_i + D_j) · 86400 / d_ij²                     [1/day]        (wbp:105)
W       = max(W, Wᵀ)                                        (symmetrise)   (wbp:108)
L       = W − diag(rowsum W)                                → row sums ≡ 0  (wbp:109)
```
- **(a) Math.** A finite-volume discretisation of ∇·(D∇C) on the k-nearest-neighbour graph of the real cell
  centroids (`k = 6`, wbp:77). Edge conductance is the **arithmetic mean diffusivity ÷ distance²** — the standard
  two-point flux approximation. `86400` converts µm²/s → µm²/day... and note the length units cancel:
  `[µm²/day]/[µm²] = [1/day]`, so `L` is a **rate matrix**, not a length-scaled operator. **Row sums are exactly
  zero** (wbp:109 comment) ⇒ the pure-diffusion operator conserves mass to machine precision.
- **(b) Biology.** An antibody in the interstitium random-walks through a collagen/hyaluronan mesh. Where the
  matrix is dense, its effective diffusivity drops. `k = 6` neighbours ≈ the coordination number of a packed 2D
  cell sheet: each cell exchanges with its immediate physical neighbours.
- **(c) Rationale / rejected alternative.** A regular finite-difference grid was rejected implicitly by the
  design: the model's whole premise is that the **cells are the mesh** — every barcode is a node (wbp:8). A graph
  Laplacian on real centroids preserves the actual tissue architecture (vessel spacing, nest geometry) that a
  regular grid would smear. The 3 nm floor on `d_ij` (wbp:104) is a physical floor (two membranes cannot be closer)
  that prevents `1/d²` blowing up on duplicated coordinates.
- **(d) Units.** `ecm_i` [dimensionless, here a 0/1 indicator]; `α_D` = 3.0 [dimensionless]; `D_um2s` = 10.0
  [µm²/s]; `d_ij` [µm]; `g_ij`, `L` [1/day].
- **⚠ MAJOR HONESTY FLAG — "ECM" is not ECM.** `phi_D` is computed from the argument `ecm`, and the only live
  caller passes `stro` — a **binary stromal-cell-type indicator** built by string-matching cell labels against
  `STROMAL_LABELS` (`coupled_percell_pk.py:83`; `STROMAL_LABELS` at `wholebody_percell.py:31-32`). So
  `ecm_i ∈ {0,1}` and `φ_D ∈ {1, 0.25}` — **a stromal cell diffuses 4× slower, every other cell is unhindered.**
  The module also defines an `ECM_GENES` list (`wholebody_percell.py:30`: COL11A1, FAP, VCAN, THBS1, ACTA2,
  TAGLN, PDGFRA, THBS4, COL17A1, LAMC2) — **this list is referenced nowhere in the live code** (grepped: its
  definition is the only hit). The *intended* continuous, expression-derived ECM score is **not wired**. What is
  live is a coarse cell-type proxy.

---

### EQ-14 — Lymphatic catchment kernel (wholebody_percell.py:116-126)

```
dl_i     = distance from cell i to the NEAREST LEC cell         [µm]
drain_w_i = exp( −dl_i / λ ) + 0.05,    λ = 100 µm                          (wbp:124)
(if no LEC in the organ: drain_w ≡ 1)                                        (wbp:126)
```
- **(a) Math.** An exponential decay kernel in distance-to-nearest-lymphatic, plus a **floor of 0.05** so that no
  cell has literally zero drainage.
- **(b) Biology.** Interstitial fluid convects toward the initial lymphatics; the convective velocity (and hence
  the antibody flux carried with it) is highest near the vessel. λ = 100 µm is stated in the code comment as the
  physiological inter-lymphatic spacing / catchment length.
- **(c) Rationale.** The floor of 0.05 is a **numerical safety device with a physical excuse**: without it, a cell
  1 mm from any LEC would have `w ≈ 0`, all of its drainage weight would be redistributed to near-LEC cells by the
  normalisation in EQ-12, and drug would accumulate irreversibly in the far field. With the floor, the entire
  tissue drains and there is no accumulation artifact (wbp:199-200).
- **(d) Units.** `dl` [µm]; λ [µm]; `drain_w` [dimensionless].
- **Provenance:** λ = 100 µm and the 0.05 floor have **no citation in the code**.

---

### EQ-15 — Time integration: operator-split implicit-Euler diffusion (wholebody_percell.py:135-137, 206-207)

```
LU = splu( I − dt·L )                       # prefactored ONCE at build           (wbp:136-137)
local_i = src_i − drain_field_i − loss_i    # all point sources/sinks [nM/day]     (wbp:206)
C ← max( LU.solve( C + dt·local ) , 0 )                                            (wbp:207)
```
- **(a) Math.** Backward-Euler on the diffusion operator, **forward** (explicit) on the local source/sink terms —
  a standard IMEX split. Solving `(I − dt·L) C^{n+1} = C^n + dt·local` is unconditionally stable for the diffusion
  part regardless of `dt`, `k`, or the cell spacing. The sparse LU is computed **once** at construction and reused
  every step, so each step costs one triangular solve (O(nnz)), not a factorisation.
- **(b) Biology.** Nothing biological — this is what makes a tissue graph of the size the code itself anticipates
  ("kNN over up to ~400k cells", wbp:58) integrable over 21+ days.
- **(c) Rationale / rejected alternative.** Explicit diffusion would require `dt < d²/(2·D)`. With `d` down to
  the 3 nm floor and `D` = 10 µm²/s = 864 000 µm²/day, that is a sub-microsecond step — as unusable as the
  explicit vascular pool (EQ-9). **The same stiffness argument that forces vascular QSS forces implicit diffusion.**
  The two design choices are the same choice, applied to the two fastest processes in the model.
- **(d) Guard.** The `max(·, 0)` clamp prevents any negative concentration surviving into the next step. With the
  loss clamp (EQ-18) it should never bind; it is a hard safety net.

---

### EQ-16 — Per-cell antigen sink, KINETIC branch (wholebody_percell.py:172-187)

Active when `kdeg`, `kon1`, `koff1` are all supplied (`set_antigen`, wbp:152).

```
Ag1EFF_i = geo_ageff_nM( R_i , r_cell = 8 µm , span = span_coeng_nm )   if n_arm ≥ 2
         = 0                                                            if n_arm < 2       (wbp:178)

(Ag, BAg1, Bdbl, intern_flux) = rhoden_samecell_bivalent_step(
      C, Ag, BAg1, Bdbl, Ag0, kon1, koff1, Ag1EFF, kdeg, kint, dt )                        (wbp:179-181)

intern_copies_i = ( intern_flux_i / NM_PER_COPY ) · count_scale     [copies/cell/day]      (wbp:184)
loss_i          = intern_copies_i / AVO / vcell · 1e9               [nM/day]               (wbp:185)
S_i = (BAg1_i + Bdbl_i)/NM_PER_COPY      D_i = Bdbl_i/NM_PER_COPY    [copies/cell]         (wbp:186-187)
```
with `NM_PER_COPY = 6.0/257000 = 2.3346e-5 nM/copy` (wbp:41) — which implies a reaction volume of
`1/(NM_PER_COPY · 1e-9 · AVO)` = **71.1 pL** (computed; matches the "71 pL synapse" comment in the PD module).

- **(a) Math.** Receptor state is carried **across steps** (`Ag`, `BAg1`, `Bdbl` are members, wbp:149) — this is a
  true kinetic ODE, not an equilibrium assumption. The `Ag0` set-point drives receptor synthesis/turnover at
  `kdeg`. `Ag1EFF` is the Rhoden **geometric effective second-arm concentration**: the local concentration a
  tethered second arm experiences, ∝ (surface density) × (SA_arm/V_arm) = (surface density) × 3/(2·span).
- **(b) Biology.** An antibody arm binds a receptor (`kon1`) → `BAg1`. If the molecule has ≥2 arms against the
  **same** antigen, the second arm now searches a small hemisphere of radius `span_coeng_nm` on the *same cell
  membrane*, at effective concentration `Ag1EFF` — if it finds a partner it crosslinks into `Bdbl` (avidity).
  Bound complexes are internalised at `kint`; free receptor turns over at `kdeg`. **This is where TMDD is born**:
  the drug is consumed by being eaten, cell by cell.
- **(c) Rationale / rejected alternative — three, all recorded in the code:**
  1. **Backward-Euler sub-cycling inside the binding step**, not the outer step (wbp:175-177): "stable +
     census-exact + fast (nsub from slow scales, NOT pinned by the stiff crosslink)". The crosslink reaction is
     the stiffest thing in the model (`Ag1EFF` reaches 1e3–1e4 nM — computed: 2478 nM at 1e4 copies, 12 388 nM at
     5e4 copies, at span 12.5 nm / r_cell 8 µm). Letting it set the global `dt` would be ruinous.
  2. The **rejected predecessor** is named in the comment (wbp:177): the old `rhoden_bivalent_step` two-pool call,
     whose `Bdbl` was **INERT for single-antigen** because it was fed from `BAg2 = 0`. i.e. the previous code
     silently had **no avidity at all** for a same-antigen bivalent binder. `rhoden_samecell_bivalent_step` fixes it.
  3. `n_arm < 2 ⇒ Ag1EFF = 0` — monovalent binders get **no avidity by construction**, not by a tuned-to-zero
     parameter.
- **(d) Units.** `Ag1EFF`, `Ag`, `BAg1`, `Bdbl`, `intern_flux` [nM or nM/day on the 71 pL synapse volume];
  `kon1` [1/nM/day]; `koff1`, `kdeg`, `kint` [1/day]; `loss` [nM/day on the *interstitial pericellular* volume];
  `S`,`D` [copies/cell].
- **⚠ The `(S + 2D)`-style flux is antibody flux, not receptor flux.** The internalisation flux out of the **live**
  routine `rhoden_samecell_bivalent_step` is `kTMD·(BAg1 + 2·Bdbl)` (docstring, `kinetic_rhoden_percell.py:98`;
  the function is defined at `:83`). *(Do not confuse this with `kinetic_rhoden_percell.py:46`, which is the
  docstring of the **rejected predecessor** `rhoden_bivalent_step` — a different, 6-species flux
  `kTMD·(BAg1 + BAg2 + 2·Bdbl)`. That function is not called from `wholebody_percell.py`.)*
  The factor 2 on the
  doubly-bound pool is **correct as an antibody flux** under independent-receptor internalisation: a doubly-engaged
  antibody is dragged in when *either* of its two receptors internalises, so its removal rate is `2·kint`. It is
  *not* a double-count of the drug. (Under the alternative reading — "count occupied receptors" — the numbers
  coincide, which is why this is easy to get wrong.)
- **⚠ Binding does not deplete free drug.** `C` is an **input only** to the binding step; the call at wbp:179-181
  returns four values and `C` is not among them. Free interstitial drug is reduced **only** by internalisation
  (`loss`), never by the act of binding. So the bound pool is not subtracted from the free pool: the model assumes
  the surface-bound reservoir is small compared with the interstitial free pool, or equivalently that binding is
  instantaneously re-supplied. At very high receptor density and low dose (the antigen-sink-dominated regime) this
  will **over-estimate** free interstitial drug. This is a genuine, unhedged approximation.

---

### EQ-17 — Per-cell antigen sink, QSS fallback branch (wholebody_percell.py:35-52, 188-191)

Active when any of `kon1/koff1/kdeg` is `None` (i.e. only a KD is known for the molecule).

```
# effective 2nd-arm concentration (wbp:35-38), r_cell = 8 µm, r_Ab = 0.0125 µm = 12.5 nm:
Ag_eff(R) = ( R/AVO·1e9 / SA_cell ) · SA_Ab / V_Ab · 1e15,
            SA_cell = 4πr_cell², SA_Ab = πr_Ab², V_Ab = (2/3)πr_Ab³      ⇒ SA_Ab/V_Ab = 3/(2·r_Ab)

# monovalent (n_arm < 2), wbp:45:
S = R · C/(C + KD) ;   D = 0

# bivalent (n_arm ≥ 2), wbp:46-52 — closed-form root of the avidity quadratic:
a = 2·Ag_eff·C ;   b = KD·R·(2C + KD) ;   c = −KD²·R²
R_free = ( −b + √(b² − 4ac) ) / (2a)            [clipped to 0 ≤ R_free ≤ R]
S      = (2C/KD)·R_free
D      = Ag_eff·R_free·S / (2·R·KD)
```
- **(a) Math.** The quadratic `a·R_free² + b·R_free + c = 0` is **exactly** the receptor-conservation identity
  `R = R_free + S + 2D` (I verified this algebraically: substituting `S` and `D` and dividing through by `R²`
  reproduces the quadratic term-for-term). `S` = singly-engaged antibodies (factor 2 = two arms, so the monovalent
  on-rate statistical factor is explicit); `D` = doubly-engaged (avidity) antibodies, each occupying **two**
  receptors. The `a > 1e-30` guard falls back to the linear (no-avidity) solution `R/(1 + 2C/KD)`.
- **(b) Biology.** Same physics as EQ-16, but assuming binding is **instantaneously equilibrated** with the local
  free-drug concentration. Valid when `koff` ≫ the transport rates — true for most engagers (koff ~1e-4/s ≈ 9/day,
  vs. transport ~0.1–1/day), which is why this fallback is defensible rather than merely convenient.
- **(c) Rationale / rejected alternative.** The rejected alternative is refusing to simulate molecules for which
  only a `KD` was published. The QSS path lets the counterscreen run on partially-characterised constructs. The
  cost is that receptor **turnover** (`kdeg`) and receptor **depletion** by internalisation are absent: `R` is
  frozen at its initial value in this branch, so it cannot report antigen down-modulation.
- **(d) Units.** `R` [copies/cell]; `C`, `KD`, `Ag_eff` [nM]; `S`, `D` [copies/cell — see the flag below];
  `r_cell` [µm]; `r_Ab` [µm].
- **⚠ `S` means different things in the two branches (readout inconsistency).** In the kinetic branch,
  `S = (BAg1 + Bdbl)/NM_PER_COPY` = **all bound antibody** (wbp:186). In the QSS branch, `S` = **singly-bound
  only**, with the doubly-bound pool held separately in `D` (wbp:51). Anything that consumes `S` alone — e.g. the
  `bound_nM` record (EQ-22) — therefore reports a *different quantity* depending on which branch ran. The sink
  itself is not affected (each branch computes its own internalisation flux consistently: kinetic via
  `intern_flux`, QSS via `kint·(S + 2D)`, wbp:190).

---

### EQ-18 — Sink clamp (wholebody_percell.py:192)

```
loss_i ← min( loss_i , 0.9 · C_i / dt )
```
- **(a) Math.** No cell may lose more than 90% of its free drug in one step.
- **(b) Rationale.** With an *explicit* treatment of the local sink (EQ-15), a large `kint·R` at small `C` could
  drive `C` negative. This is a numerical guard, and the 0.9 (rather than 1.0) leaves headroom for the
  simultaneous drainage term.
- **⚠ Mass-balance flag.** The clamp is applied to `loss` (which enters the interstitial update, wbp:206) but
  **`organ_sink` is computed from the *unclamped* `intern_copies_cell`** (wbp:208). So in any step where the clamp
  binds, the **reported TMDD sink exceeds the drug actually removed from the interstitium**. In the production
  regime (dt = 0.02 d, kint ~0.1–0.5/day) the clamp should essentially never bind — but the code does not check,
  and does not warn. **Recommended fix (not applied here — this doc does not modify code): compute `organ_sink`
  from the clamped `loss` field.**

---

### EQ-19 — Organ sink aggregation → emergent TMDD (wholebody_percell.py:208)

```
organ_sink = Σ_i intern_copies_i / AVO · 1e9              [nmol/day]
```
- **(a) Math.** copies/day → mol/day → nmol/day. Consistency check (verified algebraically):
  `Σ_i loss_i · vcell = Σ_i intern_copies_i/AVO·1e9 = organ_sink` — the local concentration sink and the reported
  organ flux are the **same number**, up to the clamp caveat in EQ-18.
- **(b) Biology.** **This is the target-mediated drug disposition, and it is not a parameter.** No `k_TMDD` was
  fitted; no Michaelis–Menten term was imposed on the plasma. The nonlinear, dose-dependent, target-dependent
  clearance of the engager emerges from summing what individual cells actually ate. Change the target from CD20 to
  BCMA and the TMDD changes by itself, because the per-cell copy numbers changed.
- **(c) Rationale.** This is the entire justification for a per-cell PK engine. A compartmental model would need
  `Vmax/Km` per target per organ.

---

### EQ-20 — Systemic plasma & lymph balance (coupled_percell_pk.py:153-155)

```
dA_pl/dt = infn + F_sc·J_sc + k_lymph_return·A_ly − k_cat·A_pl − Σ_o J_extrav[o]
dA_ly/dt = Σ_o drain_to_lymph[o] − k_lymph_return·A_ly

A_pl ← max(A_pl + dt·dA_pl, 0)      A_ly ← max(A_ly + dt·dA_ly, 0)
```
- **(a) Math.** Explicit Euler on the two systemic pools. This is legitimate *here* (and only here) because every
  rate in these two equations is **slow**: `k_cat ≈ 0.0385/day`; the whole-body extravasation conductance
  `Σ_o PS_ex[o] / V_pl` = **3.38 L/day ÷ 3.1 L = 1.09/day** over the 11 base organs (**3.54 L/day ⇒ 1.14/day**
  including tumor) — computed this session from the live physiology (`run_tce_pd_reval.py:48, 129` organ set;
  `PS_ex = 3·(Q/500)·(1−σ_V)`); and `k_lymph_return = 24/day`. The fastest is the lymph return at 24/day, giving an
  explicit stability limit `dt < 2/24 = 0.083 day`, above the production `dt = 0.02 day` by ~4×.
  **The vascular pool, at 95 238/day, is ~4000× faster than anything else in this equation (95 238/24 = 3968) —
  which is exactly why it had to be removed (EQ-9).**
- **(b) Biology.** Drug leaves plasma by (i) FcRn-unsalvaged catabolism and (ii) extravasation into all organs.
  It returns from the tissues **only via the lymph** (there is no direct interstitium→plasma back-flux term other
  than a negative `J_extrav`), with the whole-body lymph pool acting as a fast-turnover transit compartment.
- **(c) Rationale / rejected alternative.** The lymph pool could have been eliminated by QSS too
  (`k_lymph_return = 24/day` is fast). It was kept explicit — legitimately, since `dt·k_lymph_return =
  0.02 × 24 = 0.48 < 2` (stable, ~4× margin) and it introduces a real (if short, ~1 h) transit delay. **No
  double-counting:** plasma is debited
  `J_extrav` and *not* debited the organ TMDD sink, because internalised drug was already extravasated. The QSS
  vascular pool makes this exact — since `C_vasc` holds no mass, net plasma→organ transfer *is* `J_extrav`.
- **(d) Units.** `A_pl`, `A_ly` [nmol]; all rate constants [1/day]; fluxes [nmol/day].

---

### EQ-21 — Well-mixed organ fallback (coupled_percell_pk.py:147-152)

```
C_is  = A_wm[o] / max(Vis[o], 1e-9)
J_ex  = k_dist·L[o]·(1 − σ_V[o]) · C_pl                    # NOTE: driven by C_pl, no back-flux
J_re  = k_dist·L[o]·(1 − σ_L)    · C_is
A_wm[o] ← max( A_wm[o] + dt·(J_ex − J_re), 0 )
tot_extrav += J_ex ;  tot_drain += J_re
```
- **(a) Math.** A classical 2-pore compartment for organs that have **no spatial atlas**. Note this is the
  `Q → ∞` limit of EQ-9/EQ-10 **with the back-pressure term dropped**: `J_ex` is proportional to `C_pl` alone,
  not to `(C_pl − C_is)`. Uptake therefore never saturates on the interstitial side; the only thing balancing it is
  the separate `J_re` return. At steady state `C_is/C_pl = (1−σ_V)/(1−σ_L)`, i.e. the intended Vss ratio — so the
  *equilibrium* is right even though the *driving force* is not the true gradient.
- **(b) Biology.** Same 2-pore transport, no geometry.
- **(c) Status.** **The live runner never passes `wellmixed`** (grepped `run_tce_pd_reval.py`: zero hits; in the
  engine the only occurrences of the identifier are the parameter and its assignment, pk:64, 78) → `self.wm = {}`
  → this loop body does not execute in the production TCE runs. It is live, reachable code with no live caller.
  The **same** inert loop is duplicated in the production PD driver (`coupled_percell_pd.py:353-356`), so switching
  `wellmixed` on would activate it there too. Documented for completeness, flagged as currently inert.

---

### EQ-22 — Recorded observables (coupled_percell_pk.py:136-137, 157)

```
sink_rec[o]  ← organ_sink                                       [nmol/day]     (pk:136)
bound_rec[o] ← Σ_i S_i / AVO · 1e9 / Vis[o]   if pools[o] > 0   [nM]           (pk:137)
             ← 0.0                            otherwise                        (pk:137)
C_rec        ← A_pl/V_pl · mw/1e3                               [µg/mL]        (pk:157)
```
- **⚠ `sink_rec` does not exist in the production driver.** `simulate()` records `organ_sink` (pk:136), but
  `CoupledPerCellPD.simulate_pd()` — the driver the live runs execute (§1.4) — **discards it**: the inner worker at
  `coupled_percell_pd.py:264-265` unpacks `sink` from `g.step` and then returns only `(o, drain, J_extrav, S)`;
  there is no `sink_rec` anywhere in `coupled_percell_pd.py` (grepped). The TMDD sink still **acts** (it is the
  `loss` term inside `TissueGraph.step`, wbp:185/191) — it is simply not **reported** per organ in the PD runs.
  Any claim of the form "the model reports X nmol/day of organ TMDD" cannot be sourced from a production PD run
  without adding the record back.
- **⚠ `bound_nM` omits `count_scale`.** `S` is bound copies per **sampled** cell; `organ_sink` (EQ-16/17/19) is
  multiplied by `count_scale`, but `bound_rec` (pk:137) is **not**. The reported bound concentration is therefore
  the *sampled-cell* aggregate — it is off from the physiological bound pool by exactly a factor of `count_scale`
  (`bound_nM_reported = bound_nM_physiological / count_scale`). This is a **readout** defect only: `bound_rec` feeds
  no dynamics. It must not be plotted against a measured receptor-occupancy without correcting it.
- Recording cadence: every `nstep//rec` steps, `rec = 400` (pk:108; the PD runner uses `rec = 200`,
  `run_tce_pd_reval.py:173`).

---

## 3. PARAMETERS

### 3.1 Owned by `coupled_percell_pk.py` (defaults in the constructor / `simulate` signature)

| Symbol | Value | Units | Provenance | Source (as found in code) | Mechanistic rationale |
|---|---|---|---|---|---|
| `sigL` (σ_L) | 0.85 | — | **[FITTED: Vss/Vc ≈ 2.1]** | pk:64 (no comment); the *same* value in the live QSP module carries the comment "CALIBRATED so the steady-state interstitial:plasma amount ratio gives Vss/Vc ≈ 2.1 (class-typical mAb; pembrolizumab 2.17, trastuzumab 2.7, mosunetuzumab 2.1)" (`qsp_costim_window_v2.py:170-175`). **[UNVERIFIED CITATION]** — the three literature Vss values are stated in a comment with no PMID; I could not verify them from code. | Sets **distribution EXTENT**. SS ratio `(1−σ_V)/(1−σ_L)` is independent of `k_dist`, so σ_L is the clean Vss knob. |
| `k_dist` | 3.0 | — | **[FITTED: α-phase depth]** | pk:64; comment in `qsp_costim_window_v2.py:176-180`: "2-pore convective distribution-RATE multiplier … calibrated vs pembrolizumab day-1/3/7 fall". **[UNVERIFIED CITATION]** | Sets **distribution RATE**; cancels from the SS ratio, so it is orthogonal to σ_L. Multiplies BOTH `PS_ex` and `PS_ly`. |
| `k_lymph_return` | 24.0 | 1/day | **[ASSUMED: "fast turnover"]** | pk:65; `qsp_costim_window_v2.py:190` comment "/day lymph pool → plasma return (fast turnover)". No source. | Whole-body lymph turns over ~hourly; the exact value is immaterial as long as it is ≫ the elimination rate (it only sets a short transit delay). |
| `CLup` | 0.3503 | 1/day | **[FITTED: mosunetuzumab t½ 16.1 d]** | pk:65; `qsp_costim_window_v2.py:182-188` comment: "CALIBRATED (PK/PD validation track) so the backbone IgG-engager terminal t½ matches the mosunetuzumab clinical anchor 16.1 d exactly (FDA LUNSUMIO label)". **[UNVERIFIED CITATION]** — I did not verify the 16.1 d label value; it is a code comment. | Pinocytic uptake rate of plasma by endothelium. Only `CLup·(1−fFcRn)` is observable, so this is one fitted number, not two. |
| `k_renal` | **0.0** | 1/day | **[ASSUMED: inactive]** | pk:65. The runner never overrides it. | Renal filtration OFF. See EQ-3 flag — the size-gated law exists in the QSP module but is **not wired to this engine**. |
| `V_pl` | 3.1 | L | **[ASSUMED: reference human]** | pk:65; the same constant `V_PLASMA = 3.1` at `qsp_costim_window_v2.py:362` with the comment "plasma & lymph physiological volumes (L)". Module header says physiology is "reference-human (ICRP/Brown 71 kg)" (`qsp_costim_window_v2.py:75-76`) — **[UNVERIFIED CITATION]**, no PMID/DOI in code. | Physiological plasma volume, **not** a fitted central volume. |
| `V_ly` | 2.6 | L | **[ASSUMED: reference human]** | pk:65; `V_LYMPH = 2.6`, `qsp_costim_window_v2.py:363`. Same [UNVERIFIED CITATION] status. | **Note:** `V_ly` is stored (pk:74) but **never used** — the lymph balance (EQ-20) is written in *amounts*, so the lymph volume never enters. Dead parameter. |
| `dt` | 0.01 (runner: **0.02**) | day | [CODE-INTERNAL] | pk:65; overridable by env `PD_DT` (pk:74); runner passes 0.02 (`run_tce_pd_reval.py:143`) | Step size. See the `PD_DT` hazard in §5. |
| `D_um2s` | 10.0 | µm²/s | **[UNSOURCED — TBD]** | pk:65. No comment, no citation anywhere in either file. | Antibody diffusivity in the interstitium. This is a load-bearing number (it sets the penetration length) with **no provenance in the code**. |
| `span_coeng_nm` | 12.5 | nm | **[ASSUMED — the *value* is unsourced; the *functional form* is cited]** | pk:66; runner default 12.5 (`run_tce_pd_reval.py:145`); also `SPAN_BRIDGE_DEFAULT_NM = SPAN_CIS_DEFAULT_NM = 12.5` (`kinetic_synapse.py:39-40`). The consuming function `geo_ageff_nM` (`kinetic_rhoden_percell.py:37-38`) is documented "Rhoden geometric effective 2nd-arm conc". **Correction (2026-07-13): the code DOES carry a citation for the scheme** — `kinetic_rhoden_percell.py:3-5`: *"the bivalent CROSSLINK structure (BAg1/BAg2/Bdbl with the geometric effective 2nd-arm concentration AgEFF) is Rhoden et al. 2016 (kinetic form: bioRxiv 10.1101/2022.09.12.507653)"*. **[UNVERIFIED CITATION]** — I did not open that DOI in this task. What has **no** source anywhere in the code is the **numeric span = 12.5 nm** itself. | Co-engagement reach of the second arm. `Ag_eff ∝ 3/(2·span)` — halving the span **doubles** the avidity concentration, so this is a high-leverage number whose *value* is unsourced. |
| `F_sc` | 0.6 | — | **[UNSOURCED — TBD]** | pk:108; runner passes `cfg.get('F_sc', 0.6)` (`run_tce_pd_reval.py:175`) and every SC engager row hard-codes `F_sc=0.6` (`run_tce_pd_reval.py:79-82`) | SC absolute bioavailability. Same value for four different molecules — an assumption, not four measurements. |
| `ka_sc` | 0.25 | 1/day | **[UNSOURCED — TBD]** | pk:108; runner `cfg.get('ka', 0.25)`, hard-coded 0.25 for all SC rows (`run_tce_pd_reval.py:79-82`) | SC absorption rate → t½,abs = 2.77 d (computed). Again identical across four molecules. |
| `inf_rate`, `inf_dur` | 0.0, 0.0 | mg/day, day | [CODE-INTERNAL] | pk:108 | IV-infusion OFF by default. |
| `rec` | 400 | — | [CODE-INTERNAL] | pk:108 | Number of recorded time points. |
| `_SPD` | 86400 | s/day | [EXACT] | pk:69 | Unit conversion. |
| `AVO` | 6.02214076e23 | 1/mol | [EXACT] (SI defined) | wbp:29, imported at pk:9 | Avogadro. |
| `max_workers` | min(n_organs, 9) | — | [CODE-INTERNAL] | pk:106 | Thread pool over organs (the per-organ step is embarrassingly parallel — each organ only reads `C_pl`). |
| `_COLNAME_ALIAS` | MS4A1↔CD20, ERBB2↔HER2, FOLH1↔PSMA, TNFRSF17↔BCMA, FCRL5↔FcRH5 | — | [CODE-INTERNAL, identity-safe] | pk:29-30 | HGNC-symbol ↔ common-name aliases, tried **only as a fallback** so a resolving lookup is never changed (pk:26-28). On a true miss: **loud stderr warning + zero density** (pk:37-40) — the run continues but is explicitly marked INVALID for that target. This is the anti-silent-failure guard. |

### 3.2 Owned by `wholebody_percell.py`

| Symbol | Value | Units | Provenance | Source (as found in code) | Mechanistic rationale |
|---|---|---|---|---|---|
| `AVO` | 6.02214076e23 | 1/mol | [EXACT] | wbp:29 | — |
| `_NM_PER_COPY` | 6.0/257000 = 2.3346e-5 | nM/copy | **[DERIVED: 71 pL synapse reaction volume]** | wbp:41; "IDENTICAL to wholebody_pd.NM_PER_COPY"; the PD module comment says "2.335e-5 nM per receptor copy (71 pL synapse); FIXED for all molecules/targets/organs" (`wholebody_pd.py:83`). I verified the implied volume: 1/(2.3346e-5·1e-9·AVO) = **7.11e-11 L = 71.1 pL** (computed). | Converts receptor copies → the concentration a *synaptic* reaction volume sees. Fixing it across all molecules/targets is a deliberate choice: it means differences between targets come from **copy number**, not from a per-target volume knob. |
| `R_CELL_UM` | 8.0 | µm | **[ASSUMED]** | wbp:33; the same 8.0 appears as "target-cell radius (µm), Rhoden default" in `kinetic_synapse.py:38` — **[UNVERIFIED CITATION]**, no PMID. | Cell radius → surface area for the density→`Ag_eff` conversion. |
| `R_AB_UM_DEFAULT` | 0.0125 µm = **12.5 nm** | µm | **[ASSUMED]** | wbp:33. No comment. Numerically equals the `span_coeng_nm` default. | Antibody arm reach in the QSS `Ag_eff` (EQ-17). |
| `ECM_GENES` | COL11A1, FAP, VCAN, THBS1, ACTA2, TAGLN, PDGFRA, THBS4, COL17A1, LAMC2 | — | **DEAD — referenced nowhere** | wbp:30 | Intended continuous ECM score. **Not wired.** See EQ-13 flag. |
| `STROMAL_LABELS` | {hepatic stellate cell, fibroblast, myofibroblast cell, fibroblast_stromal, Fibroblast_stromal, stromal cell, smooth muscle cell, pericyte} | — | [CODE-INTERNAL] | wbp:31-32; consumed at pk:83 | The **actual** live "ECM" definition: a cell-type string set. |
| `alpha_D` | 3.0 | — | **[UNSOURCED — TBD]** | wbp:77. No comment, no citation. | ECM hindrance strength: `φ_D = 1/(1+3·ecm)` ⇒ stromal cells diffuse **4× slower**. Load-bearing and unsourced. |
| `k` (kNN) | 6 | neighbours | **[ASSUMED: 2D packing coordination]** | wbp:77 | Graph connectivity. |
| `D_um2s` | 10.0 | µm²/s | **[UNSOURCED — TBD]** | wbp:77 (default; also pk:65) | See §3.1. |
| distance floor | 0.003 | µm | **[ASSUMED: 3 nm membrane-to-membrane physical floor]** | wbp:104 (comment states exactly this) | Prevents `1/d²` divergence on coincident/duplicate coordinates. |
| `lam` (λ) | 100 | µm | **[ASSUMED]** | wbp:119 + comment wbp:117-118: "lymphatic catchment length (~100 um, physiological inter-lymphatic spacing)". **No citation.** | Sets how steeply drainage falls off with distance from a lymphatic. |
| drain floor | 0.05 | — | **[ASSUMED: numerical]** | wbp:124: "graded weight + small floor so ENTIRE tissue drains" | Prevents far-field accumulation (EQ-14). |
| loss clamp | 0.9 | — | **[ASSUMED: numerical]** | wbp:192 | Max fractional depletion per step. |
| BEC fallback | `n//50` cells (2%) | — | [CODE-INTERNAL fallback] | wbp:113 | If an organ has **no** BEC cells, the 2% with lowest `x+y` are declared BEC. **A silent geometric fiction** — see §5. |
| LEC fallback | `n//100` cells (1%) | — | [CODE-INTERNAL fallback] | wbp:114 | Same, highest `x+y`. |
| `86400` | | s/day | [EXACT] | wbp:105 | µm²/s → µm²/day in the Laplacian. |
| `_GRAPH_CACHE_DIR` | `/media/balthasar-lab/RAID4/atlas_spatial_omics/organ_spatial/graph_cache` (env `WB_GRAPH_CACHE`) | path | [CODE-INTERNAL] | wbp:21-27 | The kNN graph is a **pure function** of (coords, ecm, bec/lec idx, k, D, α_D) (wbp:56-60) so it is SHA1-keyed and cached (wbp:62-72). The `dt`-dependent LU is **always rebuilt** (SuperLU is unpicklable, and dt varies) — wbp:135-137. |
| ε guards | 1e-12 (Q+PS), 1e-18 (vcell, wint), 1e-30 (a, R) | — | [CODE-INTERNAL] | wbp:168, 171, 205; wbp:50-51 | Division guards. |

### 3.3 Per-molecule inputs (consumed, not owned — supplied per construct)

| Symbol | Units | Where set | Provenance |
|---|---|---|---|
| `KD_nM` | nM | `run_tce_pd_reval.py:77-88` per engager | per-molecule; **not owned by T2** |
| `n_arm` | — | same (1 for most; **2** for glofitamab 2:1) | structural (format), not fitted |
| `kint_perday` | 1/day | `antigen_kint(tgt)` (`run_tce_pd_reval.py:142`) | per-target; owned elsewhere |
| `mw_kda` | kDa | per engager (146.0–195.0) | label-derived per code comments |
| `fFcRn` | — | per engager (0.89 IgG1; 0.70 IgG2 elranatamab) | **[ASSUMED / UNSOURCED]** — the runner comment says "IgG2 fFcRn 0.70"; no citation |
| `kon1_perM_pers`, `koff1_pers`, `kdeg_perday` | 1/M/s, 1/s, 1/day | "normalized merge" of measured values (`run_tce_pd_reval.py:102-105`) | claimed MEASURED upstream; **not verified in this task** |
| `pools[o]` | nM | `run_tce_pd_reval.py:130` from a receptor table | drives `count_scale` (EQ-1) |

### 3.4 Physiology (consumed, not owned — injected from the live QSP module)

`Q, L, sigV, Vis, Vv` are built once at import by `_PBPKArrays.__init__` (`qsp_costim_window_v2.py:127-143`) from
the `_PBPK_TISSUES` table (`:82-99`, 15 organs) plus `_PLASMA_CO` (`:100`) and `_LYMPH_RATIO` (`:101`), and handed
to T2 by `run_tce_pd_reval.py:51, 143`.

| Quantity | Construction | Provenance |
|---|---|---|
| `Q[o]` | **non-lung `Qfrac` are first RENORMALISED to sum to 1** (`qsp_costim_window_v2.py:139`: `Qf[mask] /= Qf[mask].sum()`, Σ raw non-lung Qfrac = 0.945), then × `_PLASMA_CO = 5000 L/day` (`:100, :140`); lung is then overwritten to 5000 (full CO). ⇒ kidney = 0.190/0.945 × 5000 = **1005.3**, *not* 0.190 × 5000 = 950. This renormalisation is what the EQ-9 table reflects. | **[UNVERIFIED CITATION]** — header claims "reference-human (ICRP/Brown 71 kg)" (`qsp_costim_window_v2.py:75-76`), no PMID/DOI in code |
| `L[o]` | `Q[o]/500` (`_LYMPH_RATIO = 1.0/500.0`, `qsp_costim_window_v2.py:101`; applied `:141`) | code comment attributes to "Shah-Betts" (`qsp_costim_window_v2.py:101`) — **[UNVERIFIED CITATION]**, no PMID/DOI |
| `sigV[o]` | table rows `qsp_costim_window_v2.py:84-98`, 0.75 (tumor) … 0.99 (brain) | **[UNVERIFIED CITATION]** — table has no per-row source |
| `Vv[o]`, `Vis[o]` | `fV·V` (`:142`), `fIS·V` (`:143`) | same |

> **T2 does not own these and must not be cited as their source.** They are listed only because the QSS
> justification (EQ-9) is *quantitatively* a statement about `Q`, `PS_ex` and `Vv`, so the reader needs them.

---

## 4. WHAT IS EMERGENT vs IMPOSED

### Genuinely emergent (computed from mechanism, not handed over)

| Emergent quantity | From what |
|---|---|
| **TMDD / nonlinear clearance** | `Σ` per-cell internalisation (EQ-19). No `Vmax`, no `Km`, no `k_TMDD` anywhere in either file. Change the target and the clearance changes because the per-cell copy numbers changed. |
| **Bivalent avidity and its receptor-density nonlinearity** | The `Ag_eff ∝ R` term inside the binding quadratic (EQ-17) / kinetic crosslink (EQ-16). Avidity is ∝R², so high-copy cells are eaten disproportionately — *without* a "high-expressor" flag. |
| **Tissue penetration gradient** | BEC-only entry (EQ-11) + graph diffusion (EQ-13) + per-cell sink. Perivascular cells see more drug than core cells because of geometry, not a fitted "tumor penetration factor". |
| **The extravasation back-pressure** | `J_extrav ∝ (C_vasc − Cis_bec)` (EQ-10) — extravasation self-limits as the perivascular interstitium fills. |
| **Organ-to-organ exposure differences** | Emerge from `Q`, `σ_V`, `Vis` and per-organ cell composition; no per-organ tuning parameter exists in T2. |
| **`count_scale`** | Computed from the atlas vs the physiological pool (EQ-1); not chosen. |
| **Where drug exits the tissue** | The `drain_w` field (EQ-14) makes exit preferentially perilymphatic — but see below, *how much* exits is imposed. |
| **Flow- vs permeability-limitation** | The series conductance `PS·Q/(Q+PS)` (EQ-10) decides this per organ from the physiology; the model was not told which regime it is in. |

### Imposed (a constant, handed to the subsystem)

| Imposed quantity | What it is, honestly |
|---|---|
| **`k_cat` = CLup·(1−fFcRn) + k_renal** | FcRn salvage is a **lumped algebraic fraction**, not an endosomal mechanism (EQ-3). No saturation, no competition, no dose-dependence. `CLup` is **fitted to a clinical half-life**. |
| **`F_sc` = 0.6** | SC bioavailability is a **number**, not an outcome. The model does not derive it from lymphatic uptake and pre-systemic catabolism (which it structurally could). |
| **`ka_sc` = 0.25/day** | SC absorption rate is imposed, identical for four different molecules. |
| **`σ_V`, `σ_L`, `k_dist`** | The two-pore *transport extent and rate* are imposed/fitted, not derived from pore radii and Stokes radius. σ_L and k_dist are explicitly **calibration knobs**. |
| **Total lymph flux out of an organ** | `PS_ly·mean(C)` is compartmental (EQ-12). Only its *spatial distribution* is emergent. |
| **`D_um2s`, `alpha_D`, `λ`, `span`, `r_cell`** | Geometry/transport constants, **unsourced in code** (§3). |
| **The vascular compartment** | Not integrated at all — algebraically slaved (EQ-9). Correct, but it means no vascular transit-time phenomenon can emerge. |
| **The `ECM` field** | A binary stromal cell-type indicator, not an expression-derived matrix density (EQ-13 flag). |

### The honest boundary

The **drug's fate inside a tissue** is emergent. The **drug's fate in the body's plumbing** (elimination, SC
absorption, transport extent) is imposed and, in the case of `CLup`/`σ_L`/`k_dist`, **fitted to clinical PK**.
The model's claim is *not* "we predicted the half-life from first principles" — it is "given a calibrated systemic
backbone, everything that happens at the cell is computed, not assumed". Any statement stronger than that is not
supported by this code.

---

## 5. KNOWN LIMITATIONS / OPEN QUESTIONS — where a reviewer will attack

1. **`D_um2s = 10 µm²/s` and `alpha_D = 3.0` have no source in the code.** The penetration length of the drug into
   a tumor nest scales as √(D/k_sink) — so the single most visible output of the spatial model (the perivascular
   gradient) is set by an **unsourced constant**. *This must be sourced or a sensitivity analysis must be shown.*
   Same for `λ = 100 µm`, `span = 12.5 nm`, `r_cell = 8 µm`. Marked **[UNSOURCED — TBD]** in §3, not dressed up.

2. **`CLup` is fitted, and the anchor it is fitted to is a code comment.** `CLup = 0.3503` is stated (in the live
   QSP module) to be calibrated so the IgG backbone matches "the mosunetuzumab clinical anchor 16.1 d (FDA LUNSUMIO
   label)". I did **not** verify that 16.1 d value in this task — it is tagged **[UNVERIFIED CITATION]**. If it is
   wrong, every plasma curve in the model is wrong by the same factor. A committee will ask for the label page.

3. **FcRn is not a mechanism here, it is a multiplication.** No endosomal compartment, no pH-dependent binding, no
   receptor saturation. Consequences: (a) the model **cannot** produce dose-dependent half-life; (b) it cannot
   simulate FcRn competition (IVIG, Fc-engineering, anti-FcRn); (c) `fFcRn` per molecule (0.89 / 0.70) is
   **[ASSUMED]** with no citation. Given that the group's other work is *literally about FcRn*, this is the
   likeliest single line of attack.

4. **The renal-clearance law is not wired in** (EQ-3 flag). `k_renal = 0` in every live run. The engine would give
   a badly wrong half-life for a sub-glomerular no-Fc construct — which is exactly the kind of molecule a
   counterscreen might want to include. Latent, but a real correctness hole.

5. **Binding does not deplete free drug** (EQ-16 flag). `C` is read-only in the binding step; only internalisation
   removes drug. In the antigen-excess / low-dose regime this **over-estimates free interstitial drug**, which
   propagates into an over-estimate of downstream synapse formation. The magnitude has not been bounded.

6. **`organ_sink` is reported from the *unclamped* flux while the interstitium is debited the *clamped* one**
   (EQ-18 flag). If the 0.9·C/dt clamp ever binds, the TMDD ledger and the interstitial mass balance disagree, and
   nothing warns. Suggested fix stated in EQ-18.

7. **`bound_nM` omits `count_scale`** (EQ-22 flag) → it is a *sampled-cell* quantity being reported next to
   *physiological* fluxes. Readout-only, but it must not be compared to a measured receptor occupancy as-is.

8. **`S` is a different quantity in the kinetic vs QSS branch** (EQ-17 flag: all-bound vs singly-bound-only). Any
   downstream consumer of `S` alone silently changes meaning depending on whether the molecule had measured
   kon/koff. This is exactly the kind of thing that produces a "why did the number move when we added kinetics"
   ghost.

9. **The `ECM` hindrance field is a cell-type string match, not ECM.** `ECM_GENES` is defined and never used
   (EQ-13 flag). The claim "ECM-hindered diffusion" in the module docstring (wbp:10) is, as implemented,
   "stromal-cell-hindered diffusion" with a 4× penalty. Either wire the gene score or change the docstring.

10. **`PD_DT` env override desynchronises the two clocks.** `self.dt = float(os.environ.get('PD_DT', dt))` (pk:74).
    `TissueGraph` is then constructed at pk:101 with **`dt=dt` — the raw kwarg, *not* `self.dt`**. So the override
    is applied to the driver's clock and **not** propagated to the tissue: `TissueGraph.step` uses its own
    `self.dt` (wbp:162, = the kwarg) and prefactors `splu(I − dt·L)` for that kwarg (wbp:137), while the driver
    loop advances plasma with the class's `self.dt` (pk:109).
    If `PD_DT` is set to anything other than the `dt` kwarg, **the tissue and the plasma advance at different time
    steps and the implicit-Euler LU is factorised for the wrong dt.** Not triggered in the live runs (no `PD_DT`
    anywhere in `run_tce_pd_reval.py` — grepped; `dt=0.02` passed explicitly at `run_tce_pd_reval.py:143`), but it
    is a live footgun. *(Corrected 2026-07-13: an earlier draft of this doc claimed `self.dt` was read from the env
    **after** `TissueGraph` was built. That ordering is wrong — pk:74 precedes pk:101. The defect is the kwarg-vs-
    attribute mismatch, not the ordering.)*

11. **BEC/LEC fallback is a geometric fiction.** If an organ has no annotated BEC (or LEC) cells, the code silently
    declares the 2% of cells with the smallest `x+y` to be blood endothelium and the 1% with the largest `x+y` to be
    lymphatics (wbp:113-114). This produces a **corner-to-corner** drug gradient that has nothing to do with the
    tissue's real vasculature — and it does so **without a warning**. Every organ used in production must be audited
    for real BEC/LEC annotation, or this must be made loud.

12. **`Cis_bec` (BEC-local mean) vs `C.mean()` (organ mean) asymmetry** between the extravasation and drainage
    driving forces (EQ-10 flag). Defensible, but inconsistent.

13. **Dead parameters that suggest unfinished work:** `Vv` (eliminated by the QSS, EQ-9) and `V_ly` (the lymph
    balance is in amounts, EQ-20) are both passed, stored, and never read. Harmless, but they make the interface
    look richer than the physics is.

14. **The QSS itself.** It is well-justified (worst-case separation `Q/PS_ex = 667`, and the explicit step would
    need to be **953× smaller**; both computed in this session from the code's own numbers — EQ-9). But it is a
    *singular* perturbation: the vascular pool is assumed to have no mass. For an antibody, the vascular sub-volume
    holds real drug (`Vv·C_vasc`); the model's total-body mass balance therefore **omits the vascular drug mass**
    (~0.2 L × C_pl across all organs, i.e. a few percent of the central pool). This is a small, systematic,
    quantifiable bias in Vss — it should be stated, not discovered.

15. **The well-mixed fallback (EQ-21) has no back-pressure** (`J_ex ∝ C_pl`, not `C_pl − C_is`), unlike the
    per-cell organs. It is currently inert (no caller passes `wellmixed`), but if it were ever switched on, the
    two organ classes would obey **different transport laws**.

---

## Appendix — every arithmetic claim in this doc, and where it came from

All values below were **computed in this session** from constants read out of the live files; none is quoted from
memory or from literature.

| Claim | Inputs (file:line) | Result |
|---|---|---|
| lung `Q/Vv` | `_PLASMA_CO = 5000` (`qsp:100`), lung Q overwritten to full CO (`qsp:140`), `Vv = fV·V = 0.105×0.50` (`qsp:84, 142`) | 5000/0.0525 = 95 238 /day — confirms the wbp:165 comment "~1e5/day" |
| worst-case explicit-vascular stability limit | `2·Vv/(Q+PS_ex)` over all 15 tabulated organs | min = lung, **2.10e-5 day**; production dt = 0.02 d (`run_tce_pd_reval.py:143`) ⇒ **953× over the limit** |
| worst-case timescale separation | `Q/PS_ex = 500/(3·(1−σ_V)) = 166.7/(1−σ_V)` | min = **667** (tumor, σ_V=0.75 — **live**), max = 16 667 (brain, σ_V=0.99 — **tabulated, not in the live organ set**) |
| `k_cat` at fFcRn = 0.89 (IgG1 engagers) | pk:76, CLup=0.3503 (pk:65), fFcRn (`run_tce_pd_reval.py:77-79, 81-82`) | 0.03853 /day → t½ = **17.99 d** |
| `k_cat` at fFcRn = 0.70 (IgG2 elranatamab) | pk:76, `run_tce_pd_reval.py:80` | 0.1051 /day → t½ = **6.60 d** |
| SS interstitium:plasma ratio | `(1−σ_V)/(1−σ_L)`, σ_V from `qsp:84-98`, σ_L=0.85 (pk:64) | lung 0.333; tumor **1.667** |
| SC absorption t½ | ka=0.25/day (pk:108) | **2.77 d** |
| Σ organ extravasation conductance | `Σ_o PS_ex[o]`, live organ set (`run_tce_pd_reval.py:48, 129`) | **3.38 L/day** (11 base organs); **3.54 L/day** incl. tumor ⇒ /V_pl = 1.09–1.14 /day |
| 1 mg → nmol | mw = 146.0 kDa, mosunetuzumab (`run_tce_pd_reval.py:77`) | **6.849 nmol** — confirms pk:119-120 |
| 1 nM IgG in µg/mL | mw = 146.0 (`run_tce_pd_reval.py:77`); live mw range 146.0–195.0 (`:77-88`) | **0.1460 µg/mL** — confirms pk:157. *(An earlier draft of this doc used mw = 146.9 → 0.1469 µg/mL. **No molecule in the live table has mw = 146.9.** Corrected 2026-07-13.)* |
| synapse volume implied by `NM_PER_COPY` | 6.0/257000 (wbp:41) | **71.1 pL** — confirms the `wholebody_pd.py:83` comment |
| `Ag_eff` magnitude (avidity driver) | `geo_ageff_nM` (`kinetic_rhoden_percell.py:37-44`), span 12.5 nm, r_cell 8 µm ⇒ `Ag_eff = 0.2478 · R` nM | 248 nM @ 1e3 copies; 2478 nM @ 1e4; **12 388 nM @ 5e4** — this is why the crosslink is the stiffest reaction in the model |
| BEC source mass conservation | wbp:171, 196 | `Σ src·vcell = J_extrav` **exactly** (algebra) |
| drainage field mass conservation | wbp:202-205 | `Σ drain_field·vcell = drain_to_lymph` **exactly** (algebra) |
| QSS binding quadratic ≡ receptor conservation | wbp:47-51 | `R = R_free + S + 2D` reproduces the quadratic term-for-term (algebra) |

---

## 6. CORRECTION LOG — adversarial verification pass, 2026-07-13

Every claim in §1–§5 was re-checked line-by-line against the live sources. Defects found and fixed **in this file**:

| # | Defect | Was | Now |
|---|---|---|---|
| 1 | **Wrong parameter value + false citation.** `mw = 146.9 kDa` cited to `run_tce_pd_reval.py:77`. That line reads `mw=146.0`. **No engager in the live table has mw = 146.9** (values are 146.0, 195.0, 148.0, 148.5, 146.0, 146.0). The two derived numbers were therefore wrong. | 146.9 kDa; 6.807 nmol/mg; 0.1469 µg/mL per nM | 146.0 kDa; **6.849 nmol/mg**; **0.1460 µg/mL** per nM (EQ-4, EQ-7, Appendix) |
| 2 | **Mis-cited function.** The internalisation flux was cited to `kinetic_rhoden_percell.py:46` — the docstring of `rhoden_bivalent_step`, the **explicitly rejected predecessor** (flux `kTMD·(BAg1+BAg2+2·Bdbl)`), not the live `rhoden_samecell_bivalent_step`. | `kinetic_rhoden_percell.py:46` | `kinetic_rhoden_percell.py:98` (def at `:83`); the wrong-function trap is now called out (EQ-16) |
| 3 | **False absence-of-citation claim.** §3.1 asserted "no PMID/DOI appears in the code" for the Rhoden geometric term. The code **does** carry one: `kinetic_rhoden_percell.py:3-5` cites *Rhoden et al. 2016, bioRxiv 10.1101/2022.09.12.507653*. | "no PMID/DOI appears in the code" | scheme is cited (DOI, unverified); the **12.5 nm value** is what is unsourced (§3.1) |
| 4 | **Wrong `file:line` for the physiology.** `_LYMPH_RATIO` + "Shah-Betts" comment is at `qsp:101` (doc said 102); `self.L = Q·_LYMPH_RATIO` is at `qsp:141` (doc said 143 — which is `Vis`); `_PLASMA_CO` is at `qsp:100`; the σ_V table rows are `qsp:84-98`. | 102 / 143 / 83-100 / 83-101 | 101 / 141 / 84-98 / 82-101 (EQ-8, EQ-9, §3.4, Appendix) |
| 5 | **Wrong causal story in Limitation #10.** Claimed `self.dt` is read from `PD_DT` *after* `TissueGraph` is built. It is read at pk:74, **before** the build at pk:101. The real defect is that pk:101 passes the raw `dt` **kwarg**, not `self.dt`. | "after ... already been constructed" | kwarg-vs-attribute mismatch (§5 #10) |
| 6 | **Unsupported number.** "24/day is only 480× slower than the production step" — nothing in the code yields 480. | 480× | `dt·k_lymph_return = 0.02×24 = 0.48 < 2` ⇒ stable, ~4× margin (EQ-20c) |
| 7 | **Under-stated aggregate.** `Σ PS_ex ≈ 3.3 L/day` given with no derivation. | ≈3.3 | **3.38** L/day (11 base organs), **3.54** incl. tumor; inputs cited (EQ-20a, Appendix) |
| 8 | **Incomplete construction of `Q[o]`.** §3.4 said "flow fraction × 5000", which yields kidney = 950 and **contradicts the doc's own EQ-9 table (1005.3)**. The non-lung fractions are renormalised to Σ=1 first (`qsp:139`; Σ raw = 0.945). | "flow fraction × 5000" | renormalisation made explicit (§3.4) |
| 9 | **Missing live-path fact.** The production PD driver **discards `organ_sink`** (`coupled_percell_pd.py:264-265`; no `sink_rec` in that file). The doc advertised it as "recorded per organ". | "recorded per organ" | flagged: it acts but is not reported in PD runs (§1.3, EQ-22) |
| 10 | **Missing live-path fact.** Production IV runs are **not boluses** — the runner passes `iv_inf_h = 2.0` (`run_tce_pd_reval.py:176`). | EQ-4 implied bolus is what runs | caveat added (EQ-4) |
| 11 | **Scope error.** "Q/PS_ex ranges 667–16 667" presented as the live range; **brain is not in the live organ set** (`run_tce_pd_reval.py:48`). | live range | tabulated vs live distinguished (EQ-9) |
| 12 | Line-count and line-range slips: `coupled_percell_pk.py` is **159** lines (frontmatter said 160); EQ-14's comment block starts at wbp:**116**; the `~400k cells` figure is the code's own (wbp:58) and is now cited rather than asserted. | — | fixed |

**Claims re-verified as CORRECT (no change):** every equation in EQ-1…EQ-22 against its cited `pk:`/`wbp:` lines;
the EQ-9 organ table (Q, Vv, PS_ex, Q/Vv, Q/PS_ex — all reproduce from `qsp:84-98` + the renormalisation);
`k_cat` = 0.03853 /day ⇒ t½ 17.99 d and 0.1051 /day ⇒ 6.60 d; SC t½ 2.77 d; 71.1 pL synapse volume;
`Ag_eff = 0.2478·R` (248 / 2478 / 12 388 nM); the QSS-quadratic ≡ `R = R_free + S + 2D` algebra; the BEC-source and
drainage-field mass-conservation identities; `ECM_GENES` is defined at wbp:30 and referenced **nowhere** in the
engine (grepped — hits only in this doc and in archived copies); `wellmixed` is never passed by the live runner;
`Vv` and `V_ly` are stored and never read.

**Citation audit (the thing that matters most).** Every literature-flavoured string in this doc was traced to a
code comment and quoted verbatim from it. **No fabricated citation was found.** The four quoted anchors —
"pembrolizumab 2.17, trastuzumab 2.7, mosunetuzumab 2.1" (`qsp:172-173`), "mosunetuzumab clinical anchor 16.1 d
(FDA LUNSUMIO label)" (`qsp:185-186`), "reference-human (ICRP/Brown 71 kg)" (`qsp:75`), "Shah-Betts" (`qsp:101`) —
all exist in the code exactly as quoted, and all remain tagged **[UNVERIFIED CITATION]** because the code gives no
PMID/DOI for any of them. Nothing in this doc is tagged **[MEASURED]**. The dead modules (`cytokine_pbpk`,
`il6_pbpk`, `unified_binding`, `multiarm_kinetic`, `biexact_solver`, `rna_to_receptor`, `convert_copies_ALL`,
`calib_kdeath`) are cited nowhere as evidence. Nothing in this doc is claimed "validated".
