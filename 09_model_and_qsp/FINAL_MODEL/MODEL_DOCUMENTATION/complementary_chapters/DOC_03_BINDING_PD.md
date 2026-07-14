> **CITATION CAVEAT (added on integration):** this narrative chapter's factual claims are correct,
> but its own reviewer found several `wholebody_pd.py` line-number citations drifted from live source
> (notably `p_eng = B2/R_C` cited at :485, which is actually a survival line; §7.3/§7.4 also drift).
> For authoritative file:line equation citations, use the adversarially-verified **subsystems/T5,T6** docs.
> The engine changed mid-authoring; T-docs cite against a frozen checksummed snapshot, this chapter did not.

# DOC_03 — The Binding + PD Arm

**Costim-engager counter-screen QSP model — mechanistic, biological, and mathematical walkthrough of the binding and pharmacodynamic engine.**

Scope: the *in-use* execution path only. Every equation, constant, and parameter below is quoted against live engine source at
`model/engine/` with `file:line`. Where the model's design-brief shorthand differs from what the source actually computes, the **source is ground truth** and the discrepancy is flagged explicitly. Provisional / unvalidated behaviour is either omitted or labelled **PROVISIONAL**.

---

## 0. Provenance, scope, and verification method

### 0.1 Files documented (all in the frozen `LIVE_FINAL_PATH`)

This chapter covers the binding + PD arm. The engine files it traces, and their live roles:

| File | Live role in the binding + PD arm |
|---|---|
| `wholebody_pd.py` | Per-cell PD engine per organ: builds the synapse-incidence graph, runs the kill core, applies costim signaling + Treg damping, accumulates cytokine. Contains **two** kill laws (QSS + kinetic). |
| `kinetic_synapse.py` | The **canonical** kill core: per-T-cell engage → lethal-hit → detach bond ODE (exact 2×2 matrix-exponential), serial killing, avidity dwell-time, two-sided conservation. |
| `multiarm_binding.py` | Multi-arm engager geometry (CD3/costim/TAA, valency 0–2, three spans) + the `p_cis` cis-feasibility function. **Live use is narrow** (see §5.4). |
| `kinetic_rhoden_percell.py` | The bivalent Rhoden TMDD binding kernel that produces the per-cell local free drug `Cd` the PD consumes (transport binding + heme/blood sinks). |
| `coupled_percell_pd.py` | Couples PK ↔ PD in one time loop; drives each `OrganPD.step` with that organ's per-cell drug; integrates the plasma IL-6 ODE. |
| `pd_model_config.py` | Single source of truth for engine selection (`PD_ENGINE="kinetic"`) and the locked calibration (`k_hit`, `k_death`, IL-6 scale). |
| `costim_induction.py` | Activation-induced costim receptor density scaffold (default holds constitutive arms static; **refuses** inducible arms without an assumed fold). |
| `myeloid_il6.py` | The CRS output the PD arm feeds: per-cell myeloid IL-6 emitters + plasma IL-6 ODE. |
| `wholebody_percell.py` | Spatial substrate: per-cell interstitial diffusion graph + per-cell Rhoden binding that yields `Cd` (context for the PD input). |

### 0.2 Verification method

Source was read directly at `/media/balthasar-lab/RAID4/costim_engager_counterscreen/model/engine/`. Every numeric constant and equation below was located and quoted from the live `.py`. Live-vs-dead classification was cross-checked three independent ways:

1. **Import trace** — `run_tce_pd_reval.py` (the runner) imports `qsp_costim_window_v2`, `coupled_percell_pd`, `pd_model_config`; the PD chain resolves to `wholebody_pd → {kinetic_synapse, multiarm_binding, kinetic_rhoden_percell, costim_induction, myeloid_il6}`.
2. **`__pycache__` bytecode** — compiled modules exist for exactly the live files; **no `signaling_dynamics.pyc` exists anywhere** in the tree (a decisive fact, see §7).
3. **Runtime path** — `run_tce_pd_reval.py:5` prepends `{KWS}/handoff` to `sys.path`; the runtime handoff dir carries the live `.py` copies and, again, **no `signaling_dynamics.py`**.

### 0.2.1 Version anchor (line numbers are reproducible against this snapshot)

The live `model/engine/` is an active working directory; two files (`wholebody_pd.py`, `coupled_percell_pd.py`) were edited *during* authoring. To keep every `file:line` reference reproducible, all citations in this document are resolved against a frozen, checksummed snapshot taken at authoring time. Line numbers below are exact for these versions:

| File | Lines | sha256 (first 12) | mtime |
|---|---|---|---|
| `wholebody_pd.py` | 520 | `00b3fa641456` | 2026-07-13 16:34 |
| `kinetic_synapse.py` | 242 | `57510344c626` | 2026-07-13 05:43 |
| `multiarm_binding.py` | 128 | `a0875052b7a5` | 2026-07-12 18:01 |
| `kinetic_rhoden_percell.py` | 162 | `e9ad227a5cfa` | 2026-07-13 14:44 |
| `coupled_percell_pd.py` | 427 | `f72b622f2c2d` | 2026-07-13 16:34 |
| `pd_model_config.py` | 75 | `a4fba1917e36` | 2026-07-12 11:06 |
| `costim_induction.py` | 135 | `eebd80298085` | 2026-07-13 16:24 |
| `myeloid_il6.py` | 215 | `a5fb1c93e45f` | 2026-07-13 13:55 |
| `wholebody_percell.py` | 209 | `2de45274ac31` | 2026-07-13 06:28 |

The snapshot is saved alongside this document as `engine_snapshot/`. If a line number ever fails to resolve in the current source, diff against the snapshot — the *equation and constant content* was verified identical; only line positions can drift as the working directory evolves.

### 0.3 Three truths that override the design-brief shorthand

Because the brief and older docs paraphrase the code, three points must be stated up front; each is proven in the section noted.

- **T1 — The canonical kill core is the *kinetic* synapse, not the QSS trimer.** `pd_model_config.py:29` sets `PD_ENGINE = "kinetic"`. The manifest's cited lines (`wholebody_pd.py:353/360/362`) are the **QSS** copy of the shared signaling block; the *canonical* kill executes in `_step_kinetic` (`wholebody_pd.py:455`) → `kinetic_synapse.KineticSynapse.step` (`kinetic_synapse.py:148`). Both are documented; the kinetic one is what runs. (§6, §8)
- **T2 — The costim signaling layer is present but *dormant* in the clinical-TCE re-validation runs.** None of the six clinical engagers in the runner's `ENG` table set a `costim_arm`, so `OrganPD.sig = None`, the `occ → programs → g_eff` block is never entered, and `g_eff = 1.0`. It is a wired hook for the (separate) costim design sweep, not part of the validated CD3×TAA kill. (§7)
- **T3 — The costim receptor density is read *static from resting* copies.** Even when the costim layer is exercised, constitutive arms are held at resting density by construction and inducible arms (4-1BB / OX40 / ICOS / GITR) **refuse to run** without an explicitly-assumed fold. A resting-copy ranking therefore **under-rates 4-1BB / OX40 / ICOS** and produces a spurious "CD2 wins" ordering. This is a mandatory disclosure. (§7.3, §11)

---

## 1. The binding → PD signal chain (the trace)

At each PD sub-step, for one organ, for every T cell *i* and its neighbouring target cells *j*:

```
              PK/transport                     synapse graph                 kinetic bond ODE
 plasma C  ──────────────▶  per-cell local    ──────────────▶  armed B1  ──────────────▶  bridged B2
 (A_pl/V)   Rhoden TMDD      free drug  Cd_i     W (R_SYN=30µm)   (drug·CD3)   2×2 expm      (CD3·drug·TAA)
            (kinetic_rhoden)  + per-cell CD3     Wt_norm apport.  kinetic_synapse.py           │
                              + per-cell TAA                                                    │  engage→hit→detach
                                                                                                ▼  serial killing
   costim arm (DORMANT in clinical runs) ───────────────────────────────────────────▶  serial_rate (k_hit=12/d)
   occ=(Cd/(Cd+KD_cos))·(Rc/anchor)                                                             │
   → cis-coincidence gate (p_cis)                                                               ▼
   → programs (effector/supp/exh) ─▶ g_eff=exp(0.55·eff_p)·exp(−0.30·exh_p)          dkill = k_death·serial_rate
                                                                                                │
                              Treg damping  treg_damp = 1/(1+0.25·n_treg·(1+supp_extra))        ▼
                                                                      hazard  += g_eff·treg_damp·dkill  (apportioned by Wtᵀ)
                                                                      survival  = exp(−hazard)
                                                                      engaged dwell ─▶ myeloid IL-6 (CRS)
```

The remainder of this document walks each stage with full equations. Read §2 (units) first; the single conversion constant `NM_PER_COPY` threads through every downstream expression.

---

## 2. Units and the one physical constant

The entire binding + PD arm runs on **one** physical constant: the synapse reaction-volume conversion from receptor copies to a local concentration.

```
NM_PER_COPY = 6.0/257000.0 = 2.335e-5 nM per receptor copy      # wholebody_pd.py:83
                                                                # kinetic_synapse.py:37 (identical)
                                                                # coupled_percell_pd.py:12 (_NM_PER_COPY_PD)
```

**Derivation and meaning** (`wholebody_pd.py:75-83`). The value is pinned by the validated tumor kill anchor: the TAA receptor-capacity scale `Rcap_TAA = 6.0 nM` was reached at CEACAM5 = 257,000 copies/cell, so one copy ≙ `6.0/257000 = 2.335e-5 nM`. Equivalently this is a **~71 pL synapse reaction volume**: a single copy (1/N_A mol) dissolved in V gives `(1/N_A)/V` mol/L; setting that equal to `2.335e-14 mol/L` yields `V = 71.1 pL`. The same constant is used for **both** arms and **all** molecules/targets/organs.

Why it matters biologically: it replaces an earlier per-organ-mean normalization that erased absolute antigen abundance (so a low-copy target like BCMA saturated identically to high-copy CD20). With an absolute copies→nM map, a genuinely low-copy target (BCMA ~11k copies → 0.26 nM, below `KD_TAA`) forms a **weaker** bridge than a high-copy one (CD20 ~95k → 2.2 nM), and its depletion becomes exposure-dependent rather than saturating (`wholebody_pd.py:178-186`).

Per-cell arm concentrations (`wholebody_pd.py:185-186`):

```
R_CD3 = R_CD3_raw * NM_PER_COPY      # per-cell CD3 capacity (nM)
R_TAA = R_TAA_raw * NM_PER_COPY      # per-cell TAA capacity (nM)
```

Rate/time conventions: PD rates are **per day**; the kinetic-synapse step is per day; the signaling step (dormant) is per hour (`dt*24`); the myeloid/IL-6 sub-model is per hour. On-rates arrive as `/M/s` and are converted to `/nM/day` with `SPD = 86400` via `kon[/nM/day] = kon[/M/s] · 86400 / 1e9` (`wholebody_pd.py:247-249`).

---

## 3. The spatial substrate: per-cell inputs and the synapse-incidence graph

### 3.1 Per-cell inputs

Each organ is a population of real single cells with real coordinates `(x, y)`, cell-type labels, and per-cell receptor copies loaded from the agent tables (`coupled_percell_pd.attach_pd`, `coupled_percell_pd.py:17-51`):

- **T side (effector arm):** per-cell CD3 copies from the `CD3E` column → `R_CD3` (nM).
- **Target side:** per-cell antigen copies from the TAA column → `R_TAA` (nM); any cell with `R_TAA > 0` is a killable target (`wholebody_pd.py:177`, `is_target`).
- **Costim (optional):** per-cell costim-receptor copies from a `<ARM>_copies` column → `R_costim` (used only when a costim arm is wired; dormant in clinical runs).

Cell classification (`wholebody_pd.py:173-177`): T cells (`cd8`/`cd4`/`regulatory t`/`nk t`/`t cell`), CD8 subset, Treg subset (`regulatory t`/`treg`), targets (`R_TAA>0`).

### 3.2 The synapse-incidence graph W at R_SYN = 30 µm

Killing is spatial: a T cell can only engage targets within its synapse reach. The neighbourhood builder (`wholebody_pd._build_neighborhoods`, `wholebody_pd.py:267-296`) constructs, per organ:

```
R_SYN_UM  = 30.0        # T:target synapse reach (validated tumor value)   wholebody_pd.py:36
R_TREG_UM = 50.0        # Treg suppression neighbourhood                   wholebody_pd.py:71
```

Using a KD-tree over target coordinates, `query_ball_point(xy_T, r=R_SYN)` gives, for each T cell, the list of targets within 30 µm. This builds the **incidence matrix** `W` (T × target, entry 1 where a target is reachable) (`wholebody_pd.py:280-280`).

The apportionment matrix is antigen-weighted and row-normalized (`wholebody_pd.py:284-286`):

$$
W^{\text{taa}}_{ij} = W_{ij}\, R^{\text{TAA}}_j, \qquad
W^{\text{norm}}_{ij} = \frac{W^{\text{taa}}_{ij}}{\sum_k W^{\text{taa}}_{ik}}
$$

so a T cell's killing is distributed across its reachable targets **in proportion to each target's antigen density** — high-antigen neighbours absorb more of the kill. `Wt_norm` (rows sum to 1) is the operator used to map per-T kill onto per-target hazard everywhere downstream.

Two derived quantities:

- **`syn_TAA_mean`** (`wholebody_pd.py:287`): each T cell's antigen-weighted mean neighbour TAA density — the QSS path's target-arm concentration `R_B`.
- **`W` (raw 0/1 incidence)** is retained as `W_incidence` for the kinetic path's **two-sided conservation cap** (§6.5).

Treg neighbourhood (`wholebody_pd.py:294`): `n_treg[i]` = number of Tregs within `R_TREG_UM = 50 µm` of T cell *i* — the spatial driver of Treg suppression (§8.3). This is computed **always**, independent of the costim layer.

---

## 4. Binding layer A — per-cell Rhoden TMDD (produces the local drug `Cd`)

The PD consumes per-cell local free drug `Cd`. That field is produced upstream by the bivalent Rhoden TMDD binding kernel `kinetic_rhoden_percell.py`, on the interstitial diffusion graph (`wholebody_percell.py`) and as the heme/blood clearance sinks (`coupled_percell_pd.py`). This is the "binding" half of the arm; it is documented here because its output is the PD input, and because it shares one geometric physics with the synapse.

### 4.1 The geometric effective second-arm concentration

The avidity of any bivalent interaction is set by how concentrated the second epitope appears to a tethered first-bound arm. All three engines (transport, synapse, multi-arm) use **one** convention (`kinetic_rhoden_percell.geo_ageff_nM`, `kinetic_rhoden_percell.py:37-45`; identical in `kinetic_synapse.ageff_nM:53` and `multiarm_binding.geo_ageff_nM:29`):

$$
c_{\text{eff}} = \frac{A_{m}/\mathrm{SA_{cell}} \cdot \mathrm{SA_{Ab}}}{V_{Ab}}\times 10^{15}
\quad\text{with}\quad
\mathrm{SA_{cell}}=4\pi r_{\text{cell}}^2,\;
\mathrm{SA_{Ab}}=\pi r_{Ab}^2,\;
V_{Ab}=\tfrac{2}{3}\pi r_{Ab}^3
$$

where `r_Ab = span_nm/1000` (µm). **Larger span → larger explored shell → lower `c_eff` (dilution); shorter span → higher.** This is the single lever by which construct geometry (arm reach) tunes avidity, and it is byte-identical across PK and PD.

### 4.2 The 6-species bivalent TMDD scheme

`kinetic_rhoden_percell.py` implements the full Rhoden 2016 bivalent crosslink (states `Ag, BAg1, BAg2, Bdbl`) with **added TMDD receptor turnover** (synthesis `KSYN = Ag₀·kDEG`, degradation `−kDEG·Ag`, internalization `−kTMD` of every bound species) taken verbatim from the user-supplied MATLAB scheme (`kinetic_rhoden_percell.py:21-27`). At `C=0` the free receptor relaxes exactly to `Ag₀` (no QSS assumption).

The live calls use the **same-antigen reduced form** `rhoden_samecell_bivalent_step` (`kinetic_rhoden_percell.py:83-162`) — the correct 3-state system (`Ag, BAg1, Bdbl`) when both arms bind one antigen (e.g. 2×CD20 glofitamab, 2×BCMA). It is integrated **backward-Euler with receptors-as-states**, which is unconditionally stable and conserves the receptor census exactly even at the stiff crosslink rate:

$$
(\,I - h M\,)\,y_{t+h} = y_t + h\,s,\qquad y=[\,Ag,\,BAg1,\,Bdbl\,]^\top,\quad s=[\,KSYN,0,0\,]^\top
$$

with generator (`kinetic_rhoden_percell.py:123-125`):

```
M00 = −(kb + kDEG)        M01 = (koff − kx)      M02 = 2·koff
M10 = kb                  M11 = −(koff+kx+kTMD)  M12 = 2·koff
M20 = 0                   M21 = kx               M22 = −(2·koff+2·kTMD)
```

where `kb = kon·C` (first-arm on-rate), `kx = kon·x` (crosslink on-rate), `x = AgEFF·(Ag/Ag₀)` is the frozen effective second-arm concentration. A singular-row guard freezes any non-finite cell rather than crash or propagate NaN (`kinetic_rhoden_percell.py:134-158`). The internalized flux returned each step, `kTMD·(BAg1 + 2·Bdbl)`, is the TMDD contribution to plasma clearance.

**Live use of this kernel:**
- **Transport binding** (`wholebody_percell.py:178-179`): per-cell interstitial drug binds each cell's receptors; the per-cell free drug remaining is the field `g.C = Cd` handed to the PD.
- **Heme/blood TMDD sinks** (`coupled_percell_pd.py:350-368`): circulating blasts and normal blood cells bind + internalize drug, subtracted from plasma (`coupled_percell_pd.py:389`).

The PD engine (`OrganPD`) therefore receives `Cd = C_percell[Tidx]` as an input; it does not itself run the transport binding.

---

## 5. Binding layer B — multi-arm engager geometry (`multiarm_binding.py`)

This module encodes **construct format** as geometry: three arms — CD3, costim, TAA — each with valency `n ∈ {0,1,2}`, and three independent geometric spans that map to real formats (BiTE / IgG / 2+1 / tetravalent C-term fusion).

### 5.1 The three spans (`multiarm_binding.py:7-14`)

| Span | Physical meaning | Governs |
|---|---|---|
| `span_bridge_nm` | cross-cell: T-module ↔ tumor-module | synapse-cleft feasibility + bridge occupancy → **killing** and signal-2 delivery |
| `span_coeng_T_nm` | T-side co-engagement (two binders on the T end) | CD3↔costim **cis** coincidence + bivalent-CD3/costim avidity |
| `span_coeng_tumor_nm` | tumor-side co-engagement | bivalent-TAA avidity + dual-TAA co-engagement |

Default arm reach when a span is unset: `DEFAULT_ARM_REACH_NM = 12.5` (`multiarm_binding.py:27`).

### 5.2 The cis-feasibility Gaussian — `p_cis`

This is the mathematical heart of *subset selectivity by geometry*. Can one molecule co-engage CD3 **and** costim on the **same** T cell? Only if the construct's T-side co-engagement span matches the inter-epitope height gap (`multiarm_binding._cis_feasibility`, `multiarm_binding.py:39-46`):

```
d = (span_coeng_nm − gap_match_nm) / max(tol_nm, 1e-6)     # gap_match=12.5, tol=8.0
p_cis = exp(−0.5 · d²)
```

$$
\boxed{\,p_{\text{cis}}(\text{span}) = \exp\!\Big(-\tfrac{1}{2}\Big(\tfrac{\text{span}-12.5}{8}\Big)^{2}\Big)\,}
$$

> **⚠ Source-vs-brief discrepancy (flagged per truth policy).** The project brief writes this as `exp(−((span−12.5)/8)²)` — **without** the factor ½. The **live source has the ½** (`multiarm_binding.py:46`: `np.exp(-0.5*d*d)`), i.e. a Gaussian kernel with standard deviation `σ = tol = 8 nm`. Document and use the source form.

Worked values (σ = 8 nm):

| `span_coeng_T` (nm) | `d` | `p_cis` | interpretation |
|---|---|---|---|
| 12.5 (matched) | 0 | **1.00** | full cis — costim co-fires on the CD3-engaged cell (gated design) |
| 20.5 | 1 | 0.607 | partial |
| 28.5 | 2 | 0.135 | weak |
| 60.0 (tall CRD1 epitope) | 5.94 | **≈ 2.3e-8 ≈ 0** | trans / decoupled |

**Key architectural point (matches truth policy):** there is **no fixed "trans-60 nm" constant** in the physics. The `60.0` that appears is a *per-format geometry example* in the format library (`tetravalent_Cterm_trans`, `multiarm_binding.py:117-118`), representing a deliberately tall (~60 Å) costim epitope. Decoupling is **emergent** from the span's distance to the 12.5 nm gap via the Gaussian — not an imposed switch.

### 5.3 The cross-cell bridge feasibility

`multiarm_binding._bridge_feasibility` (`multiarm_binding.py:48-53`) gates the CD3×TAA (and costim×TAA) cross-cell bridge against the immune-synapse cleft: `cleft = clip(span_bridge, 13, 40)`, `feas = clip(span_bridge/cleft, 0, 1)`. Within the 13–40 nm window a matched span gives full reach. (Note: the *kinetic* synapse uses its own slightly different cleft ramp — §6.2 — because that is the engine that actually kills.)

### 5.4 Bound species and the format library — and what is actually live

`multiarm_bound` (`multiarm_binding.py:68-93`) returns the per-cell bound pools that *would* drive signals: `Cb_kill` (CD3·drug·TAA bridged trimer), `Cb_costim`, `Cb_cis` (costim coincident with CD3), `Cb_costTAA`. The monovalent arm is `R·C/(KD+C)` (`multiarm_binding.py:60`); the bivalent arm uses the saturating avidity closed form with `c_eff` (`multiarm_binding.py:63-66`). The `FORMATS` dictionary (`multiarm_binding.py:101-120`) maps named constructs (BiTE, IgG_1x1, DART_Fc, IgG_2TAA_1CD3, tetravalent_Cterm_{cis,trans}) to arm counts + spans.

> **Live-scope note (verified by import trace).** In the whole-body kinetic PD loop, `wholebody_pd.py` imports **only** `_cis_feasibility` from this module (`wholebody_pd.py:225`) to compute `p_cis`. The full `multiarm_bound` species calculator and the `FORMATS` library are the **format-geometry design API** (for the costim construct sweep); they are **not** called inside the kinetic kill loop. The kill itself comes from `kinetic_synapse` (§6). This module's live contribution to the validated CD3×TAA runs is therefore the single scalar `p_cis` — and even that is only consulted when a costim arm is wired (dormant in clinical runs, §7).

---

## 6. The kinetic synapse — the canonical kill core (`kinetic_synapse.py`)

This is the engine that kills. Per `pd_model_config.py:29` (`PD_ENGINE = "kinetic"`), every clinical PK-PD validation and design-sweep run routes through `OrganPD._step_kinetic` → `KineticSynapse.step`. It replaces the QSS equilibrium trimer (§8.4) with a **literal, time-integrated engage → lethal-hit → detach bond ODE**, so synapse-lifetime effects (serial killing, kinetic proofreading, dwell-time cytokine) are captured — the axes a CD3/TAA/costim affinity + format sweep tunes (`kinetic_synapse.py:1-30`).

### 6.1 States

Per T cell *i* (nM, synapse reaction-volume basis) — `kinetic_synapse.py:14-17`:

- **`B1_i`** = [drug·CD3] — an **armed** T cell (drug singly bound to its CD3), not yet bridged.
- **`B2_i`** = [CD3·drug·TAA] — the **bridged trimer**, summed over the T cell's alive target neighbours.

Per target cell *j*: **`surv_j ∈ [0,1]`** — survival. As targets die, `surv_j` falls, freed CD3 re-bridges the survivors → **serial killing emerges** (`kinetic_synapse.py:17`).

### 6.2 Geometry: emergent cleft and the trans effective TAA concentration

The cleft relaxes toward the bound-complex span, clamped to the physical window (`kinetic_synapse.py:41-42, 128`):

```
CLEFT_MIN_NM = 13.0        # ~TCR–pMHC dimension (floor)
CLEFT_MAX_NM = 40.0        # ceiling before the bond is mechanically unfavourable
cleft_nm = clip(span_bridge_nm, 13, 40)
```

Cleft feasibility (`kinetic_synapse.cleft_feasibility`, `kinetic_synapse.py:67-72`): `g = span_bridge/cleft`, `feas = clip((g−0.6)/0.4, 0, 1)` — below 0.6× the cleft the arm cannot bridge; at ≥1.0× it bridges fully. The armed CD3 arm then samples an **effective alive-neighbour TAA concentration** (`kinetic_synapse.py:131`):

$$
c_{\text{eff,trans}} = \text{ageff}(\text{dens}, r_{\text{cell}}{=}8, \text{span\_bridge}) \times \text{feas}
$$

### 6.3 Avidity retention: costim slows CD3 detachment

A cis-co-engaged costim arm (span-matched) buys avidity by slowing the **effective** CD3 off-rate (`kinetic_synapse.py:143`):

$$
k_{\text{off,CD3}}^{\text{eff}} = k_{\text{off,CD3}}\,\big(1 - \mathrm{clip}(\text{cis\_avidity}, 0, 0.95)\big)
$$

For the CD3×TAA clinical validation engagers `cis_avidity = 0` (`pd_model_config.py:40`), so `k_off,CD3^eff = k_off,CD3`. This is the hook the trispecific sweep uses, inactive in clinical runs.

### 6.4 The bond ODE and its exact 2×2 matrix-exponential integrator

Each step, the alive bridgeable TAA fraction per T cell is read from the survivor-weighted apportionment (`kinetic_synapse.py:163-166`): `alive_frac_i = Σ_j Wt_norm[i,j]·surv_j`, `Tfree_i = c_eff,trans,i · alive_frac_i`. With `rate_on = kon_CD3·C` and `kf = kon_TAA·Tfree`, the two-state reduction of the CD3→armed→bridged system is (`kinetic_synapse.py:169-171` documented; `:207-212` coded):

$$
\begin{aligned}
\dot B_1 &= \text{rate\_on}\,(R_C - B_1 - B_2) - k_{\text{off,CD3}}^{\text{eff}} B_1 - k_f B_1 + k_{\text{off,TAA}} B_2\\
\dot B_2 &= k_f B_1 - k_{\text{off,TAA}} B_2 - k_{\text{int}} B_2
\end{aligned}
$$

**Key avidity physics** (`kinetic_synapse.py:203-208`): from the bridged trimer `B2`, drug can leave the synapse **only** by the TAA arm releasing (`k_off,TAA` → back to `B1`) or by internalization (`k_int`). The CD3 arm releasing does **not** dissolve the trimer (drug is still held by TAA) — so the trimer lifetime is set by the *slower* (TAA) arm. This is why `B2` must not decay at `k_off,CD3`.

Because `rate_on·(R_C − B1 − B2)` makes this linear in `(B1, B2)` with a constant source `rate_on·R_C`, it is a 2×2 affine system solved **exactly** each step by `_expm2x2_apply` (`kinetic_synapse.py:74-102`):

$$
X(t{+}\Delta t)=e^{M\Delta t}X + M^{-1}\!\big(e^{M\Delta t}-I\big)b,\qquad b=[\,\text{rate\_on}\cdot R_C,\;0\,]^\top
$$

using the closed-form 2×2 matrix exponential (Sylvester: `expm(MΔt) = a₀I + a₁M`, eigenvalues `λ = ½(tr ± √(tr²−4det))`, with a coincident-eigenvalue series fallback). Unconditionally stable at the hour-scale PD step; the QSS limit emerges automatically when `koff` is fast versus the step (`kinetic_synapse.py:3-7`). The matrix entries (`kinetic_synapse.py:207-212`):

```
m11 = −(rate_on + koff_CD3_eff + kf)     m12 = koff_TAA − rate_on
m21 = kf                                 m22 = −(koff_TAA + kint)
b1  = rate_on · RC
```

### 6.5 Two-sided conservation (restores antigen-abundance grading)

Two clamps run after the bond solve:

1. **CD3-side** (`kinetic_synapse.py:213-217`): `B1 + B2 ≤ R_C` per T cell (scale down if exceeded).
2. **TAA-side** (`kinetic_synapse.py:219-227`): the bridged trimer cannot exceed the alive TAA reachable in the synapse, `taa_cap_i = Σ_j W_inc[i,j]·R_TAA_j·surv_j`. Excess `B2` is returned to `B1` (drug still held by CD3).

Biological consequence (`kinetic_synapse.py:119-123`): with the TAA cap, a **high-abundance** target is CD3-limited (full kill) and a **low-abundance** target is TAA-limited (partial) — restoring the linear-in-`R_TAA` abundance grading that the Schropp QSS engine had, which a CD3-only cap would destroy (low-copy BCMA would be killed like high-copy CD20).

### 6.6 Serial killing: the engage → hit → detach race

Once bridged, two clocks race per engaged synapse (`kinetic_synapse.py:209-227` documented):

- **`k_hit`** — deliver the lethal hit (needs sustained engagement).
- **`k_off,CD3^eff`** — detach and move on.

$$
P(\text{hit before detach}) = \frac{k_{\text{hit}}}{k_{\text{hit}}+k_{\text{off,CD3}}^{\text{eff}}},\qquad
\text{cycling rate}\approx k_{\text{off,CD3}}^{\text{eff}}
$$

The engaged fraction is `p_eng = clip(B2/R_C, 0, 1)` (`kinetic_synapse.py:216`). The **productive serial rate** (targets/day per T cell) is the harmonic product (`kinetic_synapse.py:228`):

$$
\boxed{\;\text{serial\_rate} = \frac{k_{\text{hit}}\,k_{\text{off,CD3}}^{\text{eff}}}{k_{\text{hit}}+k_{\text{off,CD3}}^{\text{eff}}}\;\cdot\;p_{\text{eng}}\;\cdot\;\text{has\_live}\;}
$$

```
K_HIT_DEFAULT = 12.0   # /day (~1 lethal hit per 2 h engaged)   kinetic_synapse.py:48
```

`k_hit = 12/day` is **FIXED from serial-killing literature, not fitted** (`kinetic_synapse.py:48-51`, `pd_model_config.py:41`). The harmonic form makes CD3 affinity a **window optimum, not a monotone knob**: slow `koff` → cannot cycle (rate → `koff`, low); fast `koff` → detaches before hitting (rate → `k_hit`, but hit probability → 0). An empirical sweep (this build) puts the endpoint kill peak at `koff ≈ 346/day` (≈ KD_CD3 40 nM × kon) — i.e. clinical CD3 affinity sits on the high-`koff` plateau, and the *window* metric (kill per engaged-dwell) improves monotonically with `koff` (`kinetic_synapse.py:218-227`).

### 6.7 Serial re-normalization onto survivors and the returned hazard

As neighbours die, a T cell's bridged capacity concentrates on the remaining live targets. The apportionment is re-normalized onto survivors each step (`kinetic_synapse.py:200-207`): `Ws[i,j] = Wt_norm[i,j]·surv_j`, rows re-normalized to sum 1 over live neighbours. The per-target hazard increment returned (`kinetic_synapse.py:233-234`):

$$
\Delta H^{\text{tgt}} = \Delta t \cdot k_{\text{death}} \cdot \big(W_s^\top \cdot \text{serial\_rate}\big)
$$

and the cumulative **engaged dwell** `dwell_engaged += Δt·p_eng` (`kinetic_synapse.py:231`) is the synapse-stability cytokine driver (§9) — a slow-`koff`, over-stable synapse emits more cytokine per kill.

---

## 7. Signaling layer C — costim occupancy → programs → effector gain (WIRED BUT DORMANT)

This section documents the costim-signaling block **and proves it is inactive in the clinical-TCE re-validation runs** (truth T2). It is the intended efficacy+toxicity core of the *costim design sweep*; it is not part of the validated CD3×TAA kill.

### 7.1 The occupancy law and cis-coincidence gate

When a costim arm is wired, each T cell's costim occupancy is (`wholebody_pd.py:353` QSS / `:465` kinetic):

$$
\text{occ} = \frac{C_d}{C_d + K_{D,\text{costim}}}\cdot\frac{R_{c,T}}{\text{anchor}},\qquad
\text{anchor} = \overline{R_c}\big|_{R_c>0}
$$

`KD_costim` default = 1.0 nM (`pd_model_config`/`OrganPD` default `KD_costim_nM=1.0`; manifest `KD_costim_nM_default=1.0`). The anchor is the mean costim-receptor copy over receptor-positive T cells (`wholebody_pd.py:218`). Occupancy is then passed through the **cis-coincidence gate** (`_apply_cis_coincidence`, `wholebody_pd.py:299-312`):

$$
\text{occ}_{\text{eff}} = \text{occ}\cdot\big[(1-p_{\text{cis}}) + p_{\text{cis}}\, f_{\text{CD3}}\big],\qquad
f_{\text{CD3}} = \frac{e}{e + \text{median}(e)}
$$

where `e` is the per-T CD3 engagement (QSS trimer `Cb`, or kinetic bridged `B2`). `p_cis = 0` (trans / no costim) → occ unchanged (cell-autonomous). `p_cis = 1` (cis, height-matched) → costim fully gated on the **same cell's** CD3 engagement — the coincident-signal design. Coincidence is emergent from span geometry × real per-cell CD3 binding. Occupancy is clipped to `[0, 5]` (`wholebody_pd.py:355/467`).

### 7.2 Programs and the effector gain

Occupancy drives a per-cell, per-program signaling integrator (`PerCellSignaling`), read as five programs. The effector gain applied to killing (`wholebody_pd.py:360-366` QSS / `:470-473` kinetic):

```
kE_gain = 0.55                                         # effector→kill sensitivity (LOCKED calib)   :360 / :470
g_eff[cd8] = exp(kE_gain · eff_p[cd8])                 # CD8/effector T only                        :362 / :471
g_eff      = g_eff · exp(−0.30 · max(exh_p, 0))        # exhaustion attenuates ALL T                :364 / :472
supp_extra = mean_over_Tregs( max(supp_p, 0) )         # extra suppression from Treg program        :366 / :473
```

$$
g_{\text{eff}} = \underbrace{\exp(0.55\,\text{eff\_p})}_{\text{CD8 only}}\cdot\underbrace{\exp(-0.30\,\max(\text{exh\_p},0))}_{\text{all T}}
$$

So agonism that raises the effector program amplifies CD8 killing (`kE_gain = 0.55`, locked); agonism that raises exhaustion attenuates it (`−0.30`); the Treg suppression program adds to the synapse-neighbourhood damping via `supp_extra`. Cytokine emission is likewise scaled by the IFN-γ/TNF/IL-2 programs (`cyto_sig_gain = max(0.2, 1 + mean(0.45·ifn + 0.32·tnf + 0.18·il2))`, `wholebody_pd.py:384-388/490-492`).

### 7.3 PROOF that this layer is dormant in the clinical runs

Three independent confirmations (all verified against source):

1. **No clinical engager wires a costim arm.** The runner's `ENG` table (`run_tce_pd_reval.py:77-88`) defines mosunetuzumab, glofitamab, epcoritamab, elranatamab, teclistamab, talquetamab (+ route variants). **None** sets a `costim_arm` key. The runner passes `costim_arm=cfg.get('costim_arm')` → `None` (`run_tce_pd_reval.py:153`), and its own comment states this "reproduces the previous plain-TCE behaviour EXACTLY" (`run_tce_pd_reval.py:146-149`).
2. **`None` costim arm ⇒ `sig = None`.** `OrganPD.__init__` enters the signaling block only `if costim_arm is not None and R_costim_percell is not None` (`wholebody_pd.py:208`); otherwise `self.sig = None` (`wholebody_pd.py:206`). In `step`/`_step_kinetic`, the `occ → programs → g_eff` block is guarded by `if self.sig is not None` (`wholebody_pd.py:338`, `:463`). With `sig = None`, **`g_eff = 1.0` and `supp_extra = 0.0`** (the initialized values, `wholebody_pd.py:335/462`).
3. **The signaling module does not exist on disk.** `OrganPD` would `import signaling_dynamics` (`wholebody_pd.py:209`) to build the integrator — but `signaling_dynamics.py` is **absent from the entire tree** (confirmed by `find`), and **no `signaling_dynamics.pyc`** was ever produced in `__pycache__` (whereas every live module has one). If the signaling branch were ever taken in these runs, the import would raise `ImportError`. Its absence is direct evidence the branch is never entered.

**Consequence for the validated CD3×TAA kill:** `g_eff = 1`, `supp_extra = 0`. Killing is driven purely by the kinetic synapse (serial rate) and the *spatial* Treg damping (§8.3); the costim occupancy/program/effector-gain machinery contributes nothing to the six clinical engagers. It is a hook for the costim sweep, documented here for completeness and honesty.

### 7.4 The static-R_costim scaffold and its refusal (`costim_induction.py`)

When the costim layer *is* exercised, receptor density is read **static from resting** copies (truth T3). `costim_induction.py` provides an activation-induced density model `R(t) = R_rest·(1 + (fold−1)·a(t))` with per-cell activation memory `da/dt = k_on·p_eng·(1−a) − k_off·a` (`costim_induction.py:14-15, 119-127`). But:

- **Constitutive arms** (CD28, CD2, CD27) carry `fold = 1.0` → `R(t) ≡ R_rest` — exactly the static resting density (`costim_induction.py:44-46`).
- **Inducible arms** (4-1BB/TNFRSF9, OX40/TNFRSF4, ICOS, GITR/TNFRSF18) have `fold = None` because the surface-density fold-change is **NOT_FOUND in the literature** (audited; only %-positive and kinetics exist) (`costim_induction.py:33-57`). With `strict=True` (how `OrganPD` instantiates it, `wholebody_pd.py:239`), an inducible arm **raises `ValueError`** rather than silently run at resting density (`costim_induction.py:107-111`). A sensitivity sweep requires an explicit `COSTIM_FOLD` env var, and the assumed fold is recorded on every result (`fold_is_assumed=True`).

This is the mechanism behind the mandatory limitation in §11: a resting-copy ranking under-rates exactly the inducible arms (4-1BB/OX40/ICOS), and either under-rates them (if run at fold=1) or refuses to run them (strict) — so a naive static ranking yields a spurious "CD2 wins."

---

## 8. The kill core: hazard accumulation and apportionment

### 8.1 Canonical (kinetic) kill law — `_step_kinetic` (`wholebody_pd.py:455-508`)

This is what runs. Per PD step:

1. Signaling modifiers (dormant in clinical runs → `g_eff=1`, `supp_extra=0`): §7.
2. Advance the synapse and get the per-target hazard increment (already `Δt·k_death·serial_rate`, apportioned by survivor-renormalized `Wsᵀ`) — `wholebody_pd.py:478`:
   ```
   dH_tgt = self.kin.step(Cd_T, dt, k_death_eff=k_death, per_target_surv)
   ```
3. Apply the effector and Treg modifiers as scalar multipliers (`wholebody_pd.py:480-482`):

$$
\text{gscale} = \overline{g_{\text{eff}}},\qquad
\text{treg\_damp} = \frac{1}{1 + \text{TREG\_K}\cdot \overline{n_{\text{treg}}}\,(1+\text{supp\_extra})},\qquad
\Delta H^{\text{tgt}} \leftarrow \Delta H^{\text{tgt}}\cdot\text{gscale}\cdot\text{treg\_damp}
$$

4. Accumulate and update survival (`wholebody_pd.py:484-485`):

$$
H_j \leftarrow H_j + \Delta H^{\text{tgt}}_j,\qquad \text{surv}_j = e^{-H_j}
$$

> **Note on the kinetic path's effector gain.** In `_step_kinetic`, `k_death_eff = k_death` exactly (`wholebody_pd.py:477`); the effector gain enters as a **scalar mean** `gscale = mean(g_eff)` multiplied onto the apportioned hazard (`wholebody_pd.py:480-482`), not per-cell. This differs from the QSS path (§8.2), which applies `g_eff` per-cell to `Cb` *before* apportionment. Because `g_eff = 1` in the clinical runs, the two are identical there; the distinction matters only for the (dormant) costim sweep.

### 8.2 The kill-core region cited in the brief (QSS `step`, `wholebody_pd.py:314-404`)

The manifest points at `wholebody_pd.py` L330–360 as "the kill core." That region is the **QSS** (Schropp equilibrium) path — the comparator engine, *not* canonical. It is documented for completeness; the canonical kill is §8.1. Its per-cell kill law (`wholebody_pd.py:367-368`):

```
Cb = ternary_equilibrium(Cd, R_CD3, syn_TAA_mean, KD_CD3, KD_TAA)   # per-cell trimer (prozone emergent)
Cb = Cb · g_eff                                                      # per-cell effector-gained trimer   :367
kill_T = Cb / (1 + TREG_K · n_treg · (1 + supp_extra))               # Treg-damped per-T kill            :368
```

$$
\boxed{\;\text{kill\_T}_i = \frac{C_{b,i}\, g_{\text{eff},i}}{1 + \text{TREG\_K}\cdot n_{\text{treg},i}\,(1+\text{supp\_extra})}\;}\qquad \text{TREG\_K}=0.25
$$

apportioned to targets and integrated as a **drug-graded hazard rate** (`wholebody_pd.py:371-379`):

$$
\Delta H = \Delta t\cdot k_{\text{death}}\cdot\big(W^{\text{norm}\,\top}\cdot \text{kill\_T}\big),\qquad H \mathrel{+}= \Delta H
$$

Here `g_eff` is per-cell (CD8 mask), applied to `Cb` before apportionment. The killed fraction is `mean(1 − e^{−H})` over targets (`wholebody_pd.py:402/506`).

### 8.3 Treg suppression — spatial, always live

```
TREG_K = 0.25          # per-Treg suppression constant (validated tumor value)   wholebody_pd.py:72
```

$$
\text{damping} = \frac{1}{1 + 0.25\cdot n_{\text{treg}}\,(1+\text{supp\_extra})}
$$

`n_treg` is the count of Tregs within 50 µm of the T cell (§3.2), computed every step regardless of the costim layer — so **Treg suppression via spatial density is always active** in the clinical runs. Only the *signaling-driven* extra term `supp_extra` (from the costim suppression program) is dormant (`supp_extra=0` when `sig=None`). The functional form (kill divided by `1 + 0.25·n_treg`) is the validated tumor ABM constant (`wholebody_pd.py:12`).

### 8.4 The Schropp ternary equilibrium (`ternary_equilibrium`, `wholebody_pd.py:85-104`)

The QSS trimer `Cb` used above is the full Schropp (2019) closed form (free receptors solved from the QE quadratic *before* the trimer, so the prozone/hook is emergent). Eqs. 26–33 as coded (`wholebody_pd.py:98-104`):

```
a = (1 + C/KD2)·C/(αKD1KD2)                          # Eq.28
b = C·(R_A − R_B)/(αKD1KD2) + (1+C/KD2)(1+C/KD1)     # Eq.29
d = −R_B·(1 + C/KD1)                                  # Eq.30
R_Bf = (−b + √(b²−4ad)) / (2a)                        # Eq.27
R_Af = R_A / (1 + C/KD1 + R_Bf·C/(αKD1KD2))           # Eq.26
Cb   = C·R_Af·R_Bf/(αKD1KD2)                          # Eq.33
```

At high `C` both arms saturate as binary complexes and the trimer collapses (prozone), which the earlier reduced-linear form `C·R_A·R_B/(KD1·KD2)` could not reproduce (`wholebody_pd.py:89-92`).

### 8.5 Apportionment summary

In **both** engines, per-T kill is mapped to per-target hazard through the antigen-weighted, row-normalized `Wt_norm` (kinetic path additionally re-normalizes onto survivors each step, §6.7): `dkill_tgt = Wt_normᵀ · (per-T kill)`. Hazard accumulates additively; survival is `exp(−hazard)`; the organ killed fraction is the mean of `1 − exp(−hazard)` over target cells.

### 8.6 The shared, locked kill constants

```
k_death = 1.0     # trimer/serial → death rate; LOCKED, ONE shared value all engagers   pd_model_config.py:63
k_hit   = 12.0    # /day serial lethal-hit ceiling; FIXED from literature, not fitted    kinetic_synapse.py:48
```

`k_death = 1.0` is passed from the runner (`run_tce_pd_reval.py:173`, `k_death = K_DEATH_KIN = PDC.K_DEATH`) into `simulate_pd` → `OrganPD.step(..., k_death)` → `_step_kinetic(..., k_death)` → `KineticSynapse.step(..., k_death)`, where it multiplies the serial-kill hazard (`kinetic_synapse.py:234`). It is the **only** PD parameter that was calibrated (anchored to the engaged-CTL serial ceiling ≈ 11.6/day at KD_CD3 = 40 nM, within the Halle band 2–16/day); the tumor 28.7% depletion is an **output** of the calibrated engine, not a fitting target (`kinetic_synapse.py:48-51`).

---

## 9. PK ↔ PD coupling in one time loop (`coupled_percell_pd.py`)

The binding + PD arm does not run in isolation: it is advanced inside the whole-body PK loop so that drug depletion (TMDD), extravasation, and PD all feel each other each step. `CoupledPerCellPD.simulate_pd` (`coupled_percell_pd.py:199-427`) is the driver.

### 9.1 The per-step sequence

At each transport step `k` (dt from config, default 0.02 d):

1. **Dose input** — SC first-order absorption `J_sc = ka·A_sc` or mass-exact IV infusion over `iv_inf_h` hours (`coupled_percell_pd.py:236-255`).
2. **Plasma concentration** `C_pl = A_pl / V_pl` (`coupled_percell_pd.py:257`).
3. **Per-organ transport + TMDD** — each organ graph steps its per-cell Rhoden binding at `C_pl`, returning extravasation, drainage, and per-cell bound drug (`coupled_percell_pd.py:259-269`).
4. **Live PD, every `pd_every` steps** (`coupled_percell_pd.py:271-273`): for each organ,
   ```
   self.pd[o].step(self.graphs[o].C, dt·pd_every, k_death=k_death)
   ```
   i.e. the per-cell local interstitial drug `g.C = Cd` (from step 3) drives the synapse/kill of §6/§8. Heme and blood compartments step against plasma drug directly (`coupled_percell_pd.py:339-340, 328-329`).
5. **Plasma IL-6** — sum every compartment's myeloid production and integrate the IL-6 ODE (§10).
6. **Plasma mass balance** (`coupled_percell_pd.py:389`):
   $$
   \dot A_{\text{pl}} = \text{inf} + F_{\text{sc}}J_{\text{sc}} + k_{\text{ly}}A_{\text{ly}} - k_{\text{cat}}A_{\text{pl}} - J_{\text{extrav}} - J_{\text{heme sink}} - J_{\text{blood sink}}
   $$
   where `k_cat = CLup·(1−fFcRn) + k_renal` (`coupled_percell_pk.py:76`) is the catabolic clearance modulated by FcRn recycling, and the heme/blood TMDD sinks are the internalized-drug fluxes from the Rhoden kernel (§4.2).

The PD is **off by default** (`pd_target=None`) so transport-only runs stay byte-identical to the validated PK core (`coupled_percell_pd.py:5-6`).

### 9.2 What the PD arm reads and writes each step

- **Reads:** per-cell interstitial drug `Cd` (from the transport/TMDD binding), per-cell CD3/TAA/(costim) copies, the synapse graph.
- **Writes:** per-target kill hazard → survival; engaged dwell → myeloid IL-6 production; cumulative + instantaneous cytokine. The instantaneous myeloid production feeds back into the plasma IL-6 concentration (the CRS readout), and the TMDD internalization feeds back into plasma drug clearance.

---

## 10. The CRS output: engaged dwell → myeloid IL-6 (`myeloid_il6.py`)

The trace ends at the toxicity readout. In the live path CRS is **mechanistic and myeloid-derived** (Giavridis/Norelli), not a fitted scale on bound drug. The PD arm's contribution is the **engaged-T signal** it hands to the myeloid emitters.

### 10.1 The bridge from PD to myeloid activation

Each PD step, the kinetic engine exposes each T cell's engaged-synapse fraction `p_eng = B2/R_C` (`wholebody_pd.py:502`). This is the CD40L-bearing, activated-T signal. Each myeloid agent integrates its **own** local engaged-T contact within a **contact** radius (not the 30 µm synapse reach):

```
R_CONTACT_UM = 14.1     # r_macrophage 10.6 + r_Tcell 3.5 (CD40L–CD40 is membrane-bound)   myeloid_il6.py:47
```

Per-cell myeloid activation ODE (`myeloid_il6.py:183-184`): `da/dt = K_ON·contact·(1−a) − K_OFF·a`, with `contact_i = Σ` engaged-T `p_eng` within 14.1 µm of myeloid cell *i*.

### 10.2 Emission and the measured per-cell constants

Only the intrinsic IL-6-secretor subset emits, at its measured per-cell maximum (`myeloid_il6.py:190`):

$$
\text{prod}_{\text{pg/hr}} = \sum_i a_i \cdot S_{\max}\cdot \text{cs}_i\cdot \mathbb{1}[\text{secretor}_i]
$$

```
S_MAX_PG_PER_HR   ≈ 0.0196 pg/hr/cell   (10.6 molec/s · 21 kDa)   myeloid_il6.py:43,111
SECRETOR_FRACTION = 0.039               (~3.9% of stimulated monocytes are IL-6 secretors)   myeloid_il6.py:53
```

No Emax, no EC50, no fitted scale — saturation emerges from the finite, spatially-distributed myeloid pool and the `(1−a)` structural cap.

### 10.3 The plasma IL-6 ODE

Summed production over all compartments (organ + blood + heme, each count-scale-lifted, `coupled_percell_pd.py:285-329`) drives the plasma IL-6 concentration, integrated with an **exact exponential step** (`myeloid_il6.py:206-214`):

$$
\frac{dC}{dt} = \frac{\text{prod}}{V} - k_{\text{deg}}C
\;\Rightarrow\;
C(t{+}\Delta t) = C_{ss} + (C - C_{ss})e^{-k_{\text{deg}}\Delta t_{\text{hr}}},\quad C_{ss}=\frac{\text{prod}/V}{k_{\text{deg}}}
$$

```
KDEG_IL6_PER_HR = 0.20        # IL-6 first-order elimination, t½ ≈ 3.5 h (PMID 31268236)   myeloid_il6.py:64
V_PLASMA_ML     = 11650.0     # IL-6 distribution space = interstitium 8.55 L + plasma 3.10 L = 11.65 L   myeloid_il6.py:65
```

The interstitial+plasma (ECF) distribution space is used because 21 kDa IL-6 is made in the interstitium and has no FcRn recycling — the physiologically correct denominator (committed before results).

### 10.4 The signaling-driven cytokine scaling is also dormant

The per-organ cytokine hierarchy (`CYTO_HIER = {IL6:1.0, IFN:0.36, TNF:0.31, IL2:0.18}`, `wholebody_pd.py:22`) accumulates from engaged dwell, and in the costim sweep would be scaled by `cyto_sig_gain` (from the IFN/TNF/IL2 programs). With `sig = None` in the clinical runs, `cyto_sig_gain = 1.0` (`wholebody_pd.py:490-492`); the IL-6 CRS readout is driven purely by the mechanistic myeloid path above.

---

## 11. Mandatory disclosures and truth-policy boundaries

### 11.1 The static-R_costim limitation (REQUIRED disclosure)

`R_costim` is set once at init from **resting** copy numbers and read unchanged every step (`wholebody_pd.py:207-215`; induction defaults to static per §7.4). The model captures costim *conditionality* through binding **geometry** (the cis-coincidence gate `p_cis`, §5.2/§7.1) but **not** through activation-induced receptor upregulation. Consequences:

- A resting-copy ranking **under-rates 4-1BB (TNFRSF9), OX40 (TNFRSF4), and ICOS** — precisely the arms that are near-absent on resting T cells and appear only after TCR engagement (which is *why* the field targets them).
- It over-rates constitutive arms (CD28, CD2, CD27).
- A naive static ranking therefore yields a **spurious "CD2 wins"** ordering that an immunologist would reject on sight (`costim_induction.py:9-11`).
- The code's honest response: constitutive arms run static (correct for them); inducible arms **refuse to run** under `strict=True` without an explicitly-assumed, recorded fold (§7.4).

**Any costim ranking produced by this engine on resting densities must be reported with this caveat.**

### 11.2 Provisional items — OMITTED or labelled

- **OX40 / GITR net-negative kill in Treg-rich settings — PROVISIONAL, not asserted.** This claim is in-chat only, not artifact-backed, and is structurally impossible under static `R_costim` (the signaling layer that would produce it is dormant, §7). It is **omitted** from the model's validated conclusions.
- **Activation-induced `R_costim`** is canonical only if wired with literature-sourced induction kinetics + a version tag; absent a sourced fold, the **static** version is canonical, and inducible arms are refused rather than guessed (§7.4).

### 11.3 Source-vs-brief discrepancies flagged in this document

1. **`p_cis` Gaussian** has a factor ½ in source (`exp(−½·((span−12.5)/8)²)`, `multiarm_binding.py:46`); the brief drops it. Source form documented (§5.2).
2. **Canonical kill core** is the *kinetic* synapse (`_step_kinetic` → `kinetic_synapse.step`), not the QSS `step` region (L297–374) the manifest cites for the occ/g_eff lines. Both documented; kinetic is canonical (§6, §8, T1).
3. **`signaling_dynamics.py` does not exist** — the occ→programs→g_eff layer cannot execute in these runs; it is a dormant hook (§7.3, T2).
4. **No fixed "trans-60 nm" span constant** exists; decoupling is emergent from span distance to 12.5 nm (§5.2).

---

## 12. Parameter and constant reference (all verified against live source)

| Symbol | Value | Meaning | Source (file:line) | Status |
|---|---|---|---|---|
| `NM_PER_COPY` | 2.335e-5 nM/copy | copies→synapse nM (71 pL) | `wholebody_pd.py:83`; `kinetic_synapse.py:37` | FIXED (tumor anchor) |
| `R_SYN_UM` | 30.0 µm | T:target synapse reach (W graph) | `wholebody_pd.py:36` | validated tumor value |
| `R_TREG_UM` | 50.0 µm | Treg suppression neighbourhood | `wholebody_pd.py:71` | validated tumor value |
| `TREG_K` | 0.25 | per-Treg suppression constant | `wholebody_pd.py:72` | validated tumor value |
| `R_CELL_UM` | 8.0 µm | cell radius (Rhoden geometry) | `kinetic_synapse.py:38` | Rhoden default |
| `CLEFT_MIN/MAX` | 13 / 40 nm | synapse cleft window | `kinetic_synapse.py:41-42` | physical bounds |
| `SPAN_BRIDGE/CIS` default | 12.5 nm | arm spans (AF3/format override) | `kinetic_synapse.py:39-40`; `pd_model_config.py:38-39` | default |
| `p_cis` gap_match / tol | 12.5 / 8 nm | cis-feasibility Gaussian center/σ | `multiarm_binding.py:39,46` | geometry |
| `k_hit` | 12.0 /day | serial lethal-hit ceiling | `kinetic_synapse.py:48`; `pd_model_config.py:41` | FIXED, not fitted |
| `k_death` | 1.0 /day | trimer/serial→death rate | `pd_model_config.py:63`; runner `:173` | LOCKED, calibrated (only fit) |
| `kE_gain` | 0.55 | effector→kill sensitivity | `wholebody_pd.py:360,470` | locked calib (dormant in clinical) |
| exhaustion coef | 0.30 | exhaustion→kill attenuation | `wholebody_pd.py:364,472` | locked (dormant) |
| `KD_costim` | 1.0 nM | costim arm affinity default | `wholebody_pd.py` `OrganPD` default | default (dormant) |
| `KD_CD3` | 40.0 nM | CD3 arm affinity default | `wholebody_pd.py:111` (`OrganPD` default) | per-construct override |
| `kint_bridge` | 0.9 /day | trimer internalization | `pd_model_config.py:37` | literature |
| `CYTO_HIER` | IL6:1, IFN:0.36, TNF:0.31, IL2:0.18 | cytokine hierarchy weights | `wholebody_pd.py:22` | mosun-anchored |
| `R_CONTACT_UM` | 14.1 µm | myeloid↔engaged-T contact | `myeloid_il6.py:47` | PMID 9400735/30571054 |
| `S_MAX` | 0.0196 pg/hr/cell | per-cell IL-6 max secretion | `myeloid_il6.py:43,111` | PMID 37533643 |
| `SECRETOR_FRACTION` | 0.039 | IL-6-secretor monocyte fraction | `myeloid_il6.py:53` | PMID 37533643 |
| `KDEG_IL6` | 0.20 /hr | IL-6 clearance (t½≈3.5h) | `myeloid_il6.py:64` | PMID 31268236 |
| `V_IL6` (ECF) | 11.65 L | IL-6 distribution space | `myeloid_il6.py:65` | interstitium+plasma |

Engine selection: `PD_ENGINE = "kinetic"` (`pd_model_config.py:29`).

---

## 13. One-paragraph summary

Local free drug `Cd` (from the per-cell bivalent Rhoden TMDD kernel, `kinetic_rhoden_percell.py`) meets per-cell CD3 on the T side and per-cell TAA density on the target side across a 30 µm synapse-incidence graph `W` (antigen-weighted, row-normalized to `Wt_norm`). Construct format enters as geometry — three arms (CD3/costim/TAA, valency 0–2) and three spans — with cis-feasibility `p_cis = exp(−½·((span−12.5)/8)²)` deciding whether costim co-engages the *same* T cell; decoupling is emergent from span distance, with no fixed trans constant. The **canonical** kill core is the kinetic synapse (`kinetic_synapse.py`): armed `B1` and bridged trimer `B2` evolve by an exact 2×2 matrix-exponential bond ODE with two-sided (CD3 and alive-TAA) conservation, and killing is the serial engage→hit→detach race `serial_rate = k_hit·koff/(k_hit+koff)·p_eng` at fixed `k_hit = 12/day`, re-bridging survivors as they die. Hazard accumulates as `Δt·k_death·serial_rate` (shared `k_death = 1`), damped by spatial Tregs `1/(1+0.25·n_treg)`, apportioned to targets by `Wt_normᵀ`; survival is `exp(−hazard)`, and engaged dwell drives mechanistic myeloid IL-6 (the CRS readout). The costim occupancy→program→effector-gain layer (`occ = (Cd/(Cd+KD))·(Rc/anchor)`, `g_eff = exp(0.55·eff_p)·exp(−0.30·exh_p)`) is fully wired but **dormant** in the six clinical CD3×TAA re-validation engagers (no `costim_arm` set; `signaling_dynamics.py` absent; `g_eff = 1`), and is the hook the separate costim design sweep exercises. Because `R_costim` is read static from resting copies, any resting-copy costim ranking under-rates the activation-induced arms (4-1BB/OX40/ICOS) and must carry that caveat.

---

*All equations, constants, and line references above were read and verified directly from live engine source under `model/engine/` at freeze (2026-07-13). Where a number is not present in source or manifest, it is not stated. Dead files (`unified_binding.py`, `biexact_solver.py`, `multiarm_kinetic.py`, `il6_pbpk.py`, `cytokine_pbpk.py`, `convert_copies_ALL.py`, `rna_to_receptor.py`, `calib_kdeath.py`) are excluded from the live path and are not documented here.*