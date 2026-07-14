# DOC_02 — PK / Distribution: The Life of the Molecule

**QSP costim-engager counter-screen — mechanistic model documentation, chapter 2**

**Scope:** the pharmacokinetic (PK) and distribution arm only — how a dosed T-cell-engager
molecule enters the body, distributes to tissue, is bound and internalized by its targets, is
salvaged or catabolized, and is cleared. The pharmacodynamic (killing / cytokine / suppression)
arm is chapter 3.

---

## 0. What this chapter documents, and what it does not

This chapter documents **only the files on the verified live execution path** — the runtime
import trace of the production runner `run_tce_pd_reval.py`. Every equation, constant, and
parameter below is quoted against live engine source at `model/engine/` with a `file:line`
citation. Nothing is taken from memory or from superseded design notes; where a number is not in
the source or the manifest, this document says so rather than inventing it.

**Live PK files documented here (and nothing else):**

| File | Role in the PK arm |
|---|---|
| `qsp_costim_window_v2.py` (L63–934) | Layer-1 **well-mixed full-body PBPK** (Shah–Betts 2-pore + FcRn); frozen physiology table + plasma/interstitium/lymph ODEs |
| `coupled_percell_pk.py` | **Per-cell PK core**: systemic plasma → organ vascular → BEC extravasation → per-cell interstitial graph → LEC drainage → lymph → plasma recirculation |
| `wholebody_percell.py` | `TissueGraph` spatial substrate: kNN diffusion graph, BEC/LEC entry/exit, per-cell Rhoden binder wiring, QSS conversion fallback |
| `kinetic_rhoden_percell.py` | Core **per-cell bivalent kinetic-TMDD sink**: receptors-as-states backward-Euler solver, free-receptor turnover, avidity crosslink |
| `multiarm_binding.py` | Multi-arm valency (CD3 / costim / TAA, each n∈{0,1,2}) + cis/trans geometric feasibility |
| `pd_model_config.py` | Single source of truth for engine flags + shared kinetic constants (read by PK too, because binding is shared) |
| `run_tce_pd_reval.py` | Harness that assembles physiology + per-molecule PK and drives the simulation |
| `coupled_percell_pd.py` | Inherits the PK core unchanged; the PK↔PD coupling loop lives here |

> **Explicitly NOT documented (dead code, excluded from the live path):**
> `unified_binding.py`, `biexact_solver.py`, `multiarm_kinetic.py`, `il6_pbpk.py`,
> `cytokine_pbpk.py`, `convert_copies_ALL.py`, `rna_to_receptor.py`, `calib_kdeath.py`.
> Older documentation labelled `unified_binding.py` the canonical binding kernel. **It is not
> loaded at runtime.** The real binding kernel on the live path is
> `kinetic_rhoden_percell.py` + `multiarm_binding.py` (+ `kinetic_synapse.py` for the PD synapse).

### 0.1 Two PK layers, one physiology

There are **two** PK implementations on the live path, and it matters which does what:

- **Layer-1, well-mixed** (`qsp_costim_window_v2.py`): a compact ODE system with one lumped
  interstitial compartment per organ, integrated with `scipy.solve_ivp`. This is the abstract
  window-scoring layer used for the nomination. Its PBPK block is also the *physiological
  reference*: the frozen tissue table (volumes, flows, reflection coefficients, cell densities)
  that every other layer reads.
- **Per-cell, spatial** (`coupled_percell_pk.py` + `wholebody_percell.py`): the same whole-body
  transport topology, but each organ's interstitium is resolved into a spatial graph of individual
  cells (from the Xenium/ABM tissue maps), each carrying its own receptor state and Rhoden binder.
  This is what the production runner `run_tce_pd_reval.py` actually instantiates
  (`run_tce_pd_reval.py:142`, class `CoupledPerCellPD`, which inherits the PK core).

Both layers draw Q, L, σ_V, V_is, V_v from the **same** `_PBPKArrays` table
(`run_tce_pd_reval.py:8,51`), so the physiology is identical; only the interstitial spatial
resolution differs. This chapter traces the molecule through the **per-cell** path (the one that
runs in production), using the well-mixed layer as the source of the physiological constants and
the closed-form statement of the same balances.

### 0.2 The load-bearing architectural claim: PK binding **is** PD binding

The single most important structural fact of this model — the one the whole three-axis nomination
leans on — is that **the drug-binding reaction is computed exactly once and used identically by PK
and PD.**

- `CoupledPerCellPD` **inherits** from `CoupledPerCellPK`
  (`coupled_percell_pd.py:16`), and its docstring states it "Extends CoupledPerCellPK (validated
  transport core, **untouched**)" and that "PD is OFF by default … so transport-only runs are
  byte-identical to CoupledPerCellPK" (`coupled_percell_pd.py:2,6`).
- The PK-only driver `simulate()` and the PK+PD driver `simulate_pd()` call the **same**
  `TissueGraph.step(...)` — compare `coupled_percell_pk.py:128` with
  `coupled_percell_pd.py:264`. That one call runs the per-cell Rhoden solve; its internalization
  flux is the TMDD sink for PK, and its bound-receptor census `S`,`D` is the engagement readout for
  PD (`kinetic_rhoden_percell.py:13–15`).

Biologically this is not a convenience — it is a correctness requirement. The drug cannot bind a
receptor "for clearance" at one affinity and "for signaling" at another; there is one occupancy
per cell per instant, and it simultaneously (a) removes drug from the interstitium (TMDD) and
(b) determines the signal delivered to that cell. Splitting them would let the model report
efficacy the pharmacokinetics never paid for. **The nomination's efficacy and toxicity axes are
therefore read from the same molecular event that sets clearance.**

---

## 1. The molecule and its units

A T-cell engager here is a multi-arm antibody construct. The PK layer needs only a handful of
**format** properties, because the entire class — whose plasma half-lives span two orders of
magnitude — is reproduced from format alone, not from free per-molecule knobs
(`qsp_costim_window_v2.py:191–199`):

- **Molar mass** `mw_kda` (default 146.9 kDa for IgG; `qsp_costim_window_v2.py:198`). Sets the
  µg/mL ↔ nM conversion and the renal sieving cutoff.
- **Fc status** `has_fc` (`qsp_costim_window_v2.py:199`). If `False`, FcRn salvage is forced to 0.
- **FcRn salvage fraction** `fFcRn` (see §7).
- **Valency** `n_arm` per side (see §8).
- **Target-arm affinity** `KD`, on/off rates `kon1`/`koff1`, and per-antigen internalization and
  turnover rates `kint`/`kdeg` (see §9).

**Amount vs. concentration.** All PK state variables are **amounts** (nmol). A dose in mg enters as
`A += mg/mw*1e3` nmol (`coupled_percell_pk.py:119–120`; `coupled_percell_pd.py:237,247`). Plasma
concentration is `C_pl = A_pl / V_pl` in **nM** (`coupled_percell_pk.py:122`). For a 60 mg dose of a
146 kDa IgG this is `60/146×10³ = 411 nmol`. Because the dose enters linearly as the initial plasma
amount, the whole trajectory scales with dose, so real-unit readouts (µg/mL, nM) are recovered by
`A_mg = A × (dose_mg/dose)`, `C_ugml = A_mg/V_pl`, `C_nM = C_ugml/mw_kda×1000`
(well-mixed readout, `qsp_costim_window_v2.py:1272–1281`).

**Rate-unit convention.** External kinetic constants are quoted in SI (`kon` /M/s, `koff` /s) and
converted once at construction to the model's internal per-day, per-nM basis using seconds-per-day
`_SPD = 86400`: `kon1 = kon1_perM_pers/1e9*86400` (→ /nM/day) and `koff1 = koff1_pers*86400`
(→ /day) — `coupled_percell_pk.py:69–71`. A standard mAb `kon = 10⁵` /M/s therefore becomes
`8.64` /nM/day; with `KD = 5` nM, `koff = kon·KD = 43.2`/day.

---

## 2. Dose: the molecule enters the body

The molecule can enter by two clinically used routes, both handled in the per-cell driver's dosing
loop (`coupled_percell_pk.py:117–124`; PD driver `coupled_percell_pd.py:234–259`).

### 2.1 Intravenous (IV)

The scheduled dose (mg at time t) is converted to nmol and deposited into the **plasma amount**:

```
A_pl += mg/mw*1e3           # IV, coupled_percell_pk.py:120
```

The production runner refines a pure bolus into a **real short infusion** — the clinical TCE
standard ~2 h — rather than an instantaneous spike (`run_tce_pd_reval.py:176`, `iv_inf_h=2.0`).
In `simulate_pd` this is realized mass-exactly: each scheduled IV dose registers an active infusion
`[t₀, t₀+dur, rate, remaining]` with `dur = iv_inf_h/24` days and `rate = (mg/mw·10³)/dur`, and
each step delivers `min(rate·dt, remaining)` so the delivered total equals the intended nmol
regardless of how `dt` aligns with the duration (`coupled_percell_pd.py:238–256`). This fixes an
otherwise ~1.2× overage when the duration is not an integer multiple of `dt`.

A separate **continuous infusion** channel (e.g. blinatumomab) is available via `inf_rate` (mg/day)
delivered for `t ≤ inf_dur` as `infn = inf_rate/mw*1e3` nmol/day (`coupled_percell_pk.py:124`;
`coupled_percell_pd.py:259`). Default 0 ⇒ pure bolus/short-infusion.

### 2.2 Subcutaneous (SC)

SC dosing deposits the dose into an **absorption depot** `A_sc` instead of plasma
(`coupled_percell_pk.py:119`):

```
A_sc += mg/mw*1e3                     # SC depot
J_sc  = ka_sc * A_sc ;  A_sc -= dt*J_sc   # first-order absorption, coupled_percell_pk.py:123
```

The depot empties by first-order absorption with rate constant `ka_sc` (default 0.25 /day), and
only a **bioavailable fraction** `F_sc` (default 0.6) reaches plasma — the plasma balance receives
`F_sc·J_sc`, not `J_sc` (see §6, `coupled_percell_pk.py:153`; `coupled_percell_pd.py:357`). The
lost `(1−F_sc)` represents pre-systemic catabolism/lymphatic loss at the injection site. These
defaults are applied for every SC-labelled molecule in the runner
(`run_tce_pd_reval.py:175,79–82`), e.g. epcoritamab, teclistamab, talquetamab, elranatamab.
Biologically SC entry is slower and lossier than IV, giving a lower, later, flatter C_max — which
is exactly why SC dosing is used clinically to blunt the first-dose cytokine spike.

---

## 3. Whole-body physiology: the frozen tissue table

Distribution runs on a **reference-human (ICRP/Brown, 71 kg) PBPK skeleton** built once at import
into `_PBPKArrays` (`qsp_costim_window_v2.py:82–163`). The runner instantiates this same table and
hands the per-organ arrays to the per-cell PK core (`run_tce_pd_reval.py:8,49–51`). The model runs
**15 tissues** — lung, heart, kidney, brain, muscle, skin, adipose, bone, stomach, small_int,
large_int, pancreas, spleen, liver — **plus tumor** (`qsp_costim_window_v2.py:82–99`).

### 3.1 Per-tissue physiological constants

For each tissue the frozen table gives (`qsp_costim_window_v2.py:82–99`):

| Symbol | Meaning | How used |
|---|---|---|
| `V` (L) | total tissue volume | split into vascular + interstitial sub-volumes |
| `Qfrac` | fraction of cardiac plasma output | organ blood flow `Q` |
| `sigV` (σ_V) | **vascular reflection coefficient** | 2-pore extravasation gate |
| `fV` | vascular volume fraction | `V_v = fV·V` |
| `fIS` | interstitial volume fraction | `V_is = fIS·V` |
| `portal` | drains to liver? | portal→liver routing (well-mixed layer) |

Selected values (`qsp_costim_window_v2.py:84–98`): plasma flow scale `_PLASMA_CO = 5000` L/day
(cardiac output × plasma fraction, L100); lung σ_V = 0.95, brain σ_V = **0.99** (near-impermeable,
the blood–brain barrier), spleen σ_V = **0.80** and liver σ_V = **0.85** (fenestrated sinusoids),
tumor σ_V = **0.75** (leakiest — the EPR-like tumor vasculature). Interstitial fraction is largest
in tumor (fIS = 0.30) and skin (0.302).

### 3.2 Blood flow, lymph flow, and the L = Q/500 rule

`_PBPKArrays.__init__` (`qsp_costim_window_v2.py:127–163`) computes:

- **Blood flow** `Q = Qfrac · 5000` L/day, with the parallel organs' fractions renormalized to
  conserve the circulatory loop and lung set to the full plasma output (it is in series)
  (`qsp_costim_window_v2.py:138–140`).
- **Lymph flow** `L = Q · (1/500)` (`qsp_costim_window_v2.py:141`, `_LYMPH_RATIO = 1/500`, L101).
  This is the **Shah & Betts (2012)** platform rule: antibody tissue uptake is *lymph-flow-limited*,
  and across tissues lymph flow is ~1/500 of blood flow. This single ratio is why mAb distribution
  is slow and restricted rather than perfusion-limited like a small molecule.
- **Sub-volumes** `V_v = fV·V` (vascular) and `V_is = fIS·V` (interstitial)
  (`qsp_costim_window_v2.py:142–143`).

> **Why lymph flow, not blood flow, governs antibody distribution.** A 150 kDa IgG cannot diffuse
> across the tight vascular endothelium; it crosses by **convection** through inter-endothelial pores,
> carried by the fluid flux that ultimately returns to blood as lymph. So the *rate* of tissue uptake
> is set by that trans-endothelial fluid flux (∝ lymph flow), and the *extent* is set by how much the
> pores sieve out the antibody (the reflection coefficients, §5). This is the mechanistic content of
> the Shah–Betts platform (`qsp_costim_window_v2.py:890–897`, citing Shah DK & Betts AM 2012,
> *J Pharmacokinet Pharmacodyn* 39:67–86, DOI 10.1007/s10928-011-9232-2, PMID 22143261).

### 3.3 Cell-density fields (drive where TMDD and tox live)

Two per-tissue density fields, both **data-grounded from the Tabula Sapiens single-cell atlas**
(~1.14M cells, CZ CELLxGENE Census 2025-11-08; normalized to spleen = 1), scale target capacity
and toxicity substrate (`qsp_costim_window_v2.py:103–123`):

- **T-cell density** `PB.tcell` (`qsp_costim_window_v2.py:111–115`) → CD3-arm receptor capacity per
  tissue, i.e. where CD3-mediated TMDD is strongest. Largest in large_int (2.36), skin (1.51),
  small_int (1.24); smallest in brain (0.02).
- **Myeloid density** `PB.myeloid` (`qsp_costim_window_v2.py:119–123`) → local CRS amplification and
  the **liver-tox substrate**. Liver myeloid = 0.84 (Kupffer cells) — the mechanistic substrate for
  4-1BB / costim hepatotoxicity, which an imposed prior would understate (L116–118). This is a PD
  quantity but it is *sited* by the same physiology table, so it is noted here.

Reference: Tabula Sapiens Consortium, *Science* 2022, DOI 10.1126/science.abl4896
(`qsp_costim_window_v2.py:110`).

### 3.4 Plasma and lymph pool volumes

The systemic pools are physiological (`qsp_costim_window_v2.py:362–363`;
`coupled_percell_pk.py:65`):

- **Plasma volume** `V_PLASMA = V_pl = 3.1 L`.
- **Lymph pool volume** `V_LYMPH = V_ly = 2.6 L`.

A deliberate modelling choice (`qsp_costim_window_v2.py:898–903`): the organ *vascular* spaces hold
drug at plasma concentration and are part of the same circulating plasma, so they are **not** added
as a separate central volume. The central volume V_c is therefore the physiological plasma volume
3.1 L — the value 2-compartment clinical fits recover (~3.0–3.5 L). Adding the vascular spaces
separately pushed V_c to 5.16 L and depressed the predicted C₀ ~40% below clinical, so it is not
done.

---

## 4. The transport topology: plasma → organ vascular → BEC → interstitium → LEC → lymph → plasma

This is the heart of the distribution model. The per-cell PK core
(`coupled_percell_pk.py`, module docstring L1–4; `wholebody_percell.py:1–11`) implements the
physiological structure the user designed:

```
                     Q (blood flow)
  systemic plasma  ───────────────►  ORGAN VASCULAR COMPARTMENT (V_v)
       ▲                                     │
       │ recirculation                       │ extravasation at BEC cells (2-pore convection, σ_V)
       │ (k_lymph_return)                     ▼
  systemic lymph pool ◄──────────  PER-CELL INTERSTITIAL GRAPH  ── diffuse (ECM-hindered) + bind (Rhoden)
                       drain at            each cell = 1 graph node AND 1 Rhoden binder
                       LEC cells (σ_L)
```

Drug **enters** tissue only at **blood-endothelial (BEC)** cells and **leaves** only at
**lymphatic-endothelial (LEC)** cells — the real vasculature-entry / lymphatic-exit points, read
from Xenium markers LYVE1/PROX1/PDPN (`wholebody_percell.py:8–11`). TMDD emerges as the sum of
per-cell internalization (`wholebody_percell.py:10–11`). The BEC/LEC index masks per organ are
supplied to the runner from `bec_lec_masks.json` (`run_tce_pd_reval.py:46`).

### 4.1 The per-cell interstitial substrate: `TissueGraph`

Each organ's interstitium is a **k-nearest-neighbour spatial diffusion graph** over its sampled
cells (`wholebody_percell.py:74–139`), built once and disk-cached (the kNN over up to ~400k cells is
the expensive step; `wholebody_percell.py:55–72,82–96`). Construction:

1. **kNN graph** (`k = 6`, `wholebody_percell.py:77,99–100`) over cell (x,y) coordinates.
2. **ECM-hindered edge conductance.** Each cell's local diffusivity is throttled by an ECM score
   from stromal/matrix genes (COL11A1, FAP, VCAN, … ; `wholebody_percell.py:30`):
   `phi_D = 1/(1 + alpha_D·ecm)` with `alpha_D = 3.0` (`wholebody_percell.py:81,77`). This is why
   a dense-stroma tumor penetrates poorly. Base diffusivity `D_um2s = 10 µm²/s`
   (`wholebody_percell.py:77`).
3. **Edge weight** between neighbours i, nb at distance d (with a 3 nm membrane-to-membrane floor):
   `ge = 0.5·(D_i + D_nb)·86400/d²` (µm²/day → /day), symmetrized `W = max(W, Wᵀ)`
   (`wholebody_percell.py:104–108`).
4. **Graph Laplacian** `L = W − diag(deg)`, row sums **exactly 0** ⇒ pure diffusion conserves mass
   (`wholebody_percell.py:108–109`).
5. **BEC/LEC masks** from the supplied indices, with safe fallbacks if a mask is empty
   (`wholebody_percell.py:110–115`).
6. **Distance-graded lymphatic drainage weight** `drain_w = exp(−d_LEC/λ) + 0.05`, λ = 100 µm
   (physiological inter-lymphatic spacing) — the whole tissue drains, fastest near lymphatics
   (`wholebody_percell.py:116–126`).
7. **Prefactored implicit-diffusion operator** `(I − dt·L)` LU-factored once with SuperLU
   (`wholebody_percell.py:135–137`); each step solves a sparse linear system (unconditionally
   stable, mass-conserving).

### 4.2 One step of transport (the `TissueGraph.step` method)

`TissueGraph.step` (`wholebody_percell.py:154–209`) advances one organ by `dt`. It is called
identically by PK's `simulate` (`coupled_percell_pk.py:128`) and PD's `simulate_pd`
(`coupled_percell_pd.py:264`). The stages:

**(1) Organ vascular compartment — quasi-steady.** Blood flow Q vastly exceeds the extravasation
conductance (Q/V_v up to ~10⁵/day), so the vascular pool equilibrates near-instantly each step and
is solved algebraically rather than integrated (which would be catastrophically stiff)
(`wholebody_percell.py:163–168`):

```
C_vasc = (Q·C_plasma + PS_extrav·C_is,BEC) / (Q + PS_extrav)       # wholebody_percell.py:168
```

where `C_is,BEC` is the mean interstitial concentration at the BEC cells. The extravasation flux
into the interstitium is then

```
J_extrav = PS_extrav · (C_vasc − C_is,BEC)          # nmol/day, wholebody_percell.py:169
```

**(2) Interstitial per-cell update.** Over the graph: implicit diffusion + a BEC source spreading
`J_extrav` over the BEC cells' pericellular volume + a distance-graded LEC drain + the per-cell
binding sink (§4.3). The concentration field update is the implicit-Euler solve
`C ← max(LU.solve(C + dt·local), 0)` where `local = source − drain_field − binding_loss`
(`wholebody_percell.py:193–207`). The drain field is normalized to integrate **exactly** to the
organ-total lymph flux `PS_lymph · C_is,mean` (`wholebody_percell.py:201–205`), so drainage is
mass-correct even though it is spatially graded.

**(3) Returns.** `(organ_sink, drain_to_lymph, J_extrav, S, D)` — the TMDD sink (nmol/day into
internalization), the lymph drainage (nmol/day), the extravasation flux, and the bound-receptor
census `S` (singly + doubly bound copies/cell) and `D` (doubly/avidity-bound)
(`wholebody_percell.py:160,208–209`).

### 4.3 Where PS_extrav and PS_lymph come from (2-pore convection)

The permeability–surface products are assembled in the driver's per-organ inner function and passed
into `step` (`coupled_percell_pk.py:127`; identically `coupled_percell_pd.py:263`):

```
PS_extrav = k_dist · L · (1 − sigV)        # plasma → interstitium
PS_lymph  = k_dist · L · (1 − sigL)        # interstitium → lymph
```

This is the 2-pore convective statement: convective flux = (fluid flux L) × (fraction of antibody
**not** sieved out, `1−σ`), scaled by the global distribution-rate knob `k_dist` (§5). σ_V gates
entry; σ_L gates return.

---

## 5. The two-pore reflection coefficients and the distribution-rate knob

Three PBPK system parameters set distribution *extent* and *rate*; all are on the `PBPK` dataclass
(`qsp_costim_window_v2.py:166–190`) and threaded to the per-cell core by the runner.

### 5.1 σ_V — vascular reflection (per organ, entry gate)

`sigV` is the per-organ vascular reflection coefficient from the frozen table (§3.1). A high σ_V
(brain 0.99) means the endothelium sieves out almost all antibody — little extravasation; a low
σ_V (tumor 0.75) means leaky vasculature and more entry. It enters as `(1 − σ_V)` in PS_extrav.

### 5.2 σ_L — lymphatic (return) reflection, one global value

`sigL = 0.85` (`qsp_costim_window_v2.py:170`). **Calibrated** so the steady-state
interstitial:plasma amount ratio gives V_ss/V_c ≈ 2.1, the class-typical mAb value (pembrolizumab
2.17, trastuzumab 2.7, mosunetuzumab 2.1) — `qsp_costim_window_v2.py:170–175`. The steady-state
ratio is `(1−σ_V)/(1−σ_L)`, **independent of k_dist**, so this knob sets the distribution *extent*
(V_ss) without touching the *rate*.

> **Worked check (verified in-kernel).** For a typical tissue σ_V = 0.95 with σ_L = 0.85:
> `(1−0.95)/(1−0.85) = 0.33` — interstitial amount ~⅓ of the vascular-equivalent at steady state
> (restricted distribution, as expected for IgG). For tumor σ_V = 0.75: `0.25/0.15 = 1.67`, so the
> leaky tumor accumulates antibody in its interstitium relative to a tight tissue.

### 5.3 k_dist — distribution-rate multiplier, one global value

`k_dist = 3.0` (`qsp_costim_window_v2.py:176–180`). A pure **rate** knob on both extravasation and
lymph return. It **cancels** from the steady-state ratio (appears identically in PS_extrav and
PS_lymph), so it sets the α-phase (distribution) *duration/depth* — calibrated to the class-typical
~1.3-day biexponential distribution (vs. the pembrolizumab day-1/3/7 fall), without changing V_ss.

The clean separation is deliberate: **σ_L sets *where* the drug ends up (extent), k_dist sets *how
fast* it gets there (rate).**

---

## 6. Recirculation and the closing plasma/lymph balance

Each step, the per-organ contributions are summed (`tot_extrav = ΣJ_extrav`,
`tot_drain = Σdrain_to_lymph`) and the systemic pools are advanced by explicit Euler
(`coupled_percell_pk.py:130–155`; identical form in `coupled_percell_pd.py:357–359`):

```
dA_pl = infusion + F_sc·J_sc + k_lymph_return·A_ly − k_cat·A_pl − tot_extrav      # ck_pk.py:153
dA_ly = tot_drain − k_lymph_return·A_ly                                            # ck_pk.py:154
A_pl  = max(A_pl + dt·dA_pl, 0);  A_ly = max(A_ly + dt·dA_ly, 0)                    # ck_pk.py:155
```

Reading the plasma balance term by term — this **is** the life of the molecule in one line:

| Term | Sign | Meaning |
|---|---|---|
| `infusion` | + | continuous IV infusion input (nmol/day), if any (§2.1) |
| `F_sc·J_sc` | + | bioavailable SC absorption from the depot (§2.2) |
| `k_lymph_return·A_ly` | + | **recirculation**: extravasated drug returns from the lymph pool |
| `−k_cat·A_pl` | − | **linear catabolic clearance** (FcRn-modulated, §7) |
| `−tot_extrav` | − | net loss to tissue interstitia (2-pore extravasation) |

The lymph pool is a fast-turnover buffer: it fills from tissue drainage (`tot_drain`) and empties to
plasma at `k_lymph_return = 24 /day` (`qsp_costim_window_v2.py:190`; `coupled_percell_pk.py:65`) —
a residence time of ~1 hour, so lymph is a transit compartment, not a storage depot. The recirculation
loop (extravasate → interstitium → drain → lymph → back to plasma) is what makes the peripheral
distribution *reversible* and produces the biexponential plasma profile: a fast α-phase (distribution
into tissue) followed by a slow β-phase (elimination, rate-limited by k_cat).

In the production PD driver the plasma balance additionally subtracts the **heme-blast** and
**normal-blood** TMDD sinks (`coupled_percell_pd.py:357`, `_heme_sink_nmol_day`,
`_blood_sink_nmol_day`) — circulating target cells that bind and internalize drug directly from
plasma (no extravasation barrier). These are additive plasma sinks with the same internalization
mechanism (§9); they are PD-attached compartments but they act on the PK plasma pool.

---

## 7. FcRn salvage and linear catabolic clearance

### 7.1 The QSS FcRn recycling model

An IgG's long half-life comes from **FcRn salvage**: pinocytosed antibody that binds FcRn in the
acidic endosome is recycled to the surface rather than routed to the lysosome. This model does **not**
integrate an explicit endosome ODE; it uses the **quasi-steady-state reduction** — a recycled
fraction `fFcRn` — for numerical robustness while preserving the platform's long IgG half-life
(`qsp_costim_window_v2.py:69–71,181`).

The net linear catabolic clearance rate is (`coupled_percell_pk.py:76`; well-mixed identical
`qsp_costim_window_v2.py:930`):

```
k_cat = CLup·(1 − fFcRn) + k_renal          # /day
```

- **`CLup = 0.3503` /day** — plasma pinocytic uptake rate (`qsp_costim_window_v2.py:182`;
  `coupled_percell_pk.py:65`). The fraction `(1 − fFcRn)` of pinocytosed drug is **not** salvaged and
  is catabolized. **Calibrated** (PK/PD validation track) so the backbone IgG-engager terminal t₁/₂
  matches the mosunetuzumab clinical anchor 16.1 d (FDA LUNSUMIO label), while the comparative window
  ranking is unchanged — tumor/plasma AUC ratio 0.234 and the prozone bell-shape are preserved
  (`qsp_costim_window_v2.py:182–188`).
- **`k_renal`** — size-gated renal/catabolic clearance (§7.3), ~0 for an intact IgG.

### 7.2 fFcRn by construct — the values actually used

`fFcRn` is a per-molecule input threaded from the runner's engager table
(`run_tce_pd_reval.py:77–88`, passed at `run_tce_pd_reval.py:142`). The **verified live values**:

| Construct | fFcRn | MW (kDa) | Rationale (from source comment) |
|---|---|---|---|
| default `PBPK` dataclass | 0.90 | 146.9 | platform default → k_cat = 0.0350/day, t₁/₂ ≈ 20 d (`qsp_costim_window_v2.py:181,188`) |
| mosunetuzumab, glofitamab, epcoritamab, teclistamab, talquetamab | 0.89 | 146–195 | IgG1/IgG4-class Fc → k_cat = 0.0385/day, t₁/₂ ≈ 18 d (`run_tce_pd_reval.py:77–82,86–88`) |
| **elranatamab** | **0.70** | 148.5 | explicitly annotated **"IgG2 fFcRn 0.70"** in source (`run_tce_pd_reval.py:80`) → k_cat = 0.105/day, t₁/₂ ≈ 6.6 d |
| no-Fc BiTE (`has_fc=False`) | forced 0 | ~54 | FcRn salvage OFF (`eff_fFcRn()=0`) → catabolism + renal dominate |

> **Verified half-life arithmetic (in-kernel).** With `CLup = 0.3503` and k_renal ≈ 0 (large IgG):
> fFcRn 0.90 → k_cat 0.0350/day → t₁/₂ **19.8 d**; fFcRn 0.89 → 0.0385/day → **18.0 d**; fFcRn 0.70
> → 0.105/day → **6.6 d**; fFcRn 0 (no-Fc) → 0.3503/day → **2.0 d** from catabolism alone (renal
> then takes it to hours). These reproduce the source comment's "t₁/₂ ~18–20 d" for IgG
> (`qsp_costim_window_v2.py:188`) and the "losing FcRn drops t₁/₂ ~18 d → ~2 d" claim
> (`qsp_costim_window_v2.py:203`).

The elranatamab IgG2 value is the one case where the model departs from the IgG1/4 platform default,
and it is source-annotated — worth flagging because IgG2 has measurably weaker FcRn recycling, which
the shorter modelled half-life reflects.

### 7.3 The no-Fc / BiTE size gate (why the class spans two orders of magnitude in t₁/₂)

The `eff_fFcRn()` accessor returns 0 when `has_fc` is False (`qsp_costim_window_v2.py:212–214`), and
a **size-gated renal clearance** turns on for small scaffolds (`qsp_costim_window_v2.py:216–218`):

```
k_renal = k_renal_max / (1 + (mw_kda / mw50_renal)^hill_renal)
```

with `k_renal_max = 8.70` /day, `mw50_renal = 69` kDa (albumin glomerular cutoff), `hill_renal = 10`
(sharp threshold) — `qsp_costim_window_v2.py:205–210`. This is a Hill sieve on molecular weight:
~0/day for a 150 kDa IgG (not filtered), dominant for a ~55 kDa BiTE (filtered + proteolysed).
**Calibrated** so a no-Fc BiTE at 54 kDa (blinatumomab) reaches the clinical ~2.1 h terminal t₁/₂
(`qsp_costim_window_v2.py:205–208`). This is the mechanistic reason one PK layer reproduces both a
weeks-long IgG-TCE and an hours-long BiTE **from format alone**: Fc status flips FcRn salvage, and
molecular weight flips renal filtration.

> **Note on layer availability.** The MW-Hill renal term and `has_fc` accessor live on the
> **well-mixed** `PBPK` dataclass (`qsp_costim_window_v2.py:200–218`). The per-cell core receives an
> explicit `k_renal` argument (default 0) and computes `k_cat = CLup·(1−fFcRn) + k_renal`
> (`coupled_percell_pk.py:65,76`); the production runner instantiates the per-cell core for
> Fc-bearing IgG-class TCEs (fFcRn 0.70–0.89) and does not pass a nonzero `k_renal`, so the BiTE
> renal branch is exercised in the well-mixed layer, not in the per-cell production runs documented
> here.

---

## 8. Valency: `n_arm` per side and why binding must be identical to PD

### 8.1 What `n_arm` controls

Every binding arm targets a cell-surface antigen with a valency n ∈ {0, 1, 2}
(`qsp_costim_window_v2.py:241–254`; `multiarm_binding.py:3,55–66`). In the per-cell PK core the
single scalar `n_arm` selects the binding regime for the target arm
(`wholebody_percell.py:44,174–178`):

- **n_arm = 0** — arm absent, no binding (`multiarm_binding.py:59`).
- **n_arm = 1** — **monovalent** 1:1 Langmuir. QSS occupancy `S = Rtot·C/(C+KD)`
  (`wholebody_percell.py:44–45`); kinetically, AgEFF = 0 so the avidity crosslink is inert and only
  singly-bound species form (`wholebody_percell.py:178`).
- **n_arm ≥ 2** — **bivalent avidity**: the second arm can crosslink a neighbouring copy of the same
  antigen, forming the doubly-bound species `Bdbl` via the Rhoden geometric effective 2nd-arm
  concentration AgEFF (`wholebody_percell.py:178`; `kinetic_rhoden_percell.py:83–93`). This is the
  glofitamab-like 2:1 (2×CD20) case (`run_tce_pd_reval.py:78`, narm = 2).

The full three-arm format space (CD3 / costim / TAA, each 0–2, up to tetravalent, three geometric
spans) lives in `multiarm_binding.py` (`multiarm_binding.py:68–93`) and is used by the PD synapse;
the PK sink uses the single target-arm valency `n_arm` that the runner threads
(`run_tce_pd_reval.py:142`, `cfg['narm']`).

### 8.2 The geometric effective 2nd-arm concentration (AgEFF)

Bivalent avidity is not a fitted enhancement factor — it is the **Rhoden 2016 geometry**
(`kinetic_rhoden_percell.py:37–45`; identical formula `multiarm_binding.py:29–37`,
`wholebody_percell.py:35–38`). Once one arm is bound, the second arm explores a small shell of
radius = the antibody arm reach; the local concentration of a second antigen copy in that shell is

```
geo_ageff_nM(rec_pc, r_cell_um=8.0, span_nm=12.5):
    SA_cell = 4π·r_cell²                      # cell surface area
    r_Ab    = span_nm (→µm)                    # arm reach as a sphere radius
    SA_Ab   = π·r_Ab² ;  V_Ab = (2/3)π·r_Ab³   # explored cap area / volume
    Ag_bulk = rec_pc·1e9/AVO·1e9               # receptor copies → bulk nM
    return (Ag_bulk / SA_cell)·SA_Ab / V_Ab · 1e15   # surface density → local nM
```

Biologically: more receptors per cell (higher `rec_pc`) → higher local 2nd-arm concentration → more
avidity; a **longer arm reach** (`span_nm`) dilutes it (larger explored shell → lower c_eff,
`multiarm_binding.py:31`). Default cell radius 8 µm and arm reach 12.5 nm
(`kinetic_rhoden_percell.py:37`; `multiarm_binding.py:27`). This same AgEFF is what couples avidity
to *conditionality* (a well-matched cis span co-engages, a mismatched one decouples — §8.4).

### 8.3 Why PK binding **must** be identical to PD binding

This is the reason the binding kernel is shared, restated at the level of the equations. The Rhoden
solve returns a bound-receptor census and an internalization flux **from the same states**
(`kinetic_rhoden_percell.py:13–15`):

- **PK uses the internalization flux** as the TMDD sink: `intern_flux = kTMD·(BAg1 + 2·Bdbl)`
  (`kinetic_rhoden_percell.py:98,161`) — every bound drug species is internalized, removing drug
  from the interstitium.
- **PD uses the bound census** `S = (BAg1 + Bdbl)/NM_PER_COPY`, `D = Bdbl/NM_PER_COPY`
  (`wholebody_percell.py:186–187`) as the engagement that drives killing/cytokine/costim signal.

If PK and PD used different binding, a molecule could report signaling that its own clearance never
accounted for (or vice versa): the receptor occupancy that removes the drug and the occupancy that
signals would disagree. Because they are the **same `S`, `D`** from the **same solve at the same
affinity/valency/geometry**, occupancy is conserved. Concretely, `simulate` (PK) and `simulate_pd`
(PD) issue the identical `g.step(C_pl, Q, PS_ex, PS_ly, KD, n_arm, kint)` call
(`coupled_percell_pk.py:128` ≡ `coupled_percell_pd.py:264`), and `CoupledPerCellPD` inherits the
transport core "untouched" (`coupled_percell_pd.py:2,16`). **This is the invariant that lets the
three-axis nomination trust that efficacy and clearance are read from one physical event.**

### 8.4 Conditionality is emergent from geometry, not a hardcoded "trans" switch

For the trispecific costim constructs, whether the costim arm co-engages the *same* T cell as CD3
(cis, gated, safe) or acts decoupled (trans-like) is set by the **cis-feasibility Gaussian**
(`multiarm_binding.py:39–46`):

```
p_cis = exp(−½·((span_coeng − 12.5)/8)²)      # gap_match = 12.5 nm, tol = 8 nm
```

A well-matched cis span (~12.5 nm) → p_cis ≈ 1 (costim co-fires on the CD3-engaged T cell); a large
epitope-height mismatch (e.g. a tall CRD1 4-1BB epitope, ~60 Å) forces the span far from 12.5 nm →
p_cis → 0 (decoupled / trans design). **There is no fixed "trans" span constant** — decoupling is an
emergent function of geometric distance to the 12.5 nm gap match (this document's manifest confirms:
"NO fixed 'trans' span constant exists in code"). The bridge feasibility across the immune-synapse
cleft is the companion gate: span clipped to 13–40 nm, `feas = clip(span/cleft, 0, 1)`
(`multiarm_binding.py:48–53`; cleft constants `kinetic_synapse.py:41–42`).

---

## 9. The per-cell Rhoden bivalent kinetic-TMDD sink (`kinetic_rhoden_percell.py`)

This is the molecular engine of clearance: the reaction that binds drug to receptor on each cell,
internalizes the bound complex, and thereby removes drug from the interstitium. It is a **kinetic**
(not equilibrium) scheme, with receptor turnover added on top of the Rhoden bivalent-binding
structure.

### 9.1 Receptors as states, with explicit turnover

The scheme carries the **free antigen (receptor) as a live state** with synthesis and degradation —
the piece a quasi-steady-state sink lacks, and the piece that makes terminal TMDD behave correctly
(`kinetic_rhoden_percell.py:17–19`):

```
dAg/dt   has  +KSYN (= Ag0·kDEG set-point synthesis)  and  −kDEG·Ag (turnover)
every bound species internalizes at kTMD
At C=0 steady state, Ag → Ag0 exactly.  No QSS assumption anywhere.
```

- **`KSYN = Ag0·kDEG`** — zero-order receptor synthesis pinned to the resting set-point Ag0, so free
  receptor returns exactly to Ag0 at zero drug (`kinetic_rhoden_percell.py:57,104`). This is the
  free-receptor turnover the task specifies as `KSYN = Rtot·kdeg`.
- **`kDEG`** (per-antigen, /day) — free-receptor degradation/turnover.
- **`kTMD = kint`** (per-antigen, /day) — internalization of every bound drug species; **this flux
  is the TMDD clearance.**

The receptor set-point on each cell is `Ag0_nM = R_copies · NM_PER_COPY`
(`wholebody_percell.py:148`), with **`NM_PER_COPY = 6.0/257000 = 2.335×10⁻⁵` nM/copy**
(`wholebody_percell.py:41`) — the synapse reaction-volume basis, held **identical to the PD engine's
`NM_PER_COPY`** (`wholebody_percell.py:41,147`). This is another face of the shared-binding
invariant: one copy of receptor is worth the same concentration to the clearance sink and to the
signaling readout.

### 9.2 The reference 6-species scheme (two co-engageable antigens)

For two distinct antigens Ag1, Ag2 (e.g. CD3 + costim on the T cell, or a CD3–TAA trans bridge), the
full kinetic bivalent system is (`kinetic_rhoden_percell.py:21–27,47–80`; species Ag1, Ag2 free +
BAg1, BAg2 singly-bound + Bdbl doubly/crosslinked):

```
dAg1  = KSYN1 − kon1·C·Ag1 + koff1·BAg1 − kon1·x1·BAg2 + koff1·Bdbl − kDEG·Ag1
dAg2  = KSYN2 − kon2·C·Ag2 + koff2·BAg2 − kon2·x2·BAg1 + koff2·Bdbl − kDEG·Ag2
dBAg1 = kon1·C·Ag1 − koff1·BAg1 − kon2·x2·BAg1 + koff2·Bdbl − kTMD·BAg1
dBAg2 = kon2·C·Ag2 − koff2·BAg2 − kon1·x1·BAg2 + koff1·Bdbl − kTMD·BAg2
dBdbl = kon1·x1·BAg2 + kon2·x2·BAg1 − (koff1+koff2)·Bdbl − 2·kTMD·Bdbl
   with  x1 = Ag1EFF·(Ag1/Ag1_0),  x2 = Ag2EFF·(Ag2/Ag2_0)   # geometric crosslink drive
```

Every bound state loses drug to internalization at kTMD (Bdbl at 2·kTMD, it carries two drug
molecules). The crosslink terms `kon·xᵢ·Bⱼ` are the avidity: a singly-bound arm reaches a second
antigen copy at the geometric effective concentration AgEFF (§8.2), scaled by current free-antigen
availability `(Ag/Ag0)`. This reference form is integrated by an adaptively sub-cycled explicit Euler
(`kinetic_rhoden_percell.py:47–80`), sub-stepped so the stiff crosslink rate cannot overshoot.

### 9.3 The 3-species same-antigen scheme — what the live PK path actually calls

For the single-antigen bivalent case (both arms bind the **same** antigen — glofitamab 2×CD20,
alnuctamab 2×BCMA — and, importantly, the general monovalent case where the crosslink is simply
inert), the model uses the **reduced 3-state** system `rhoden_samecell_bivalent_step`
(`kinetic_rhoden_percell.py:83–161`). This is the function the `TissueGraph.step` actually invokes
in production (`wholebody_percell.py:179–181`); it replaced the 2-pool call whose `Bdbl` was inert
when Ag2 = 0 (`wholebody_percell.py:176–177`). States: Ag (free), BAg1 (one arm bound), Bdbl (both
arms → two neighbour copies) (`kinetic_rhoden_percell.py:89–93`):

```
dAg   = KSYN − kon·C·Ag + koff·BAg1 − kon·x·BAg1 + 2·koff·Bdbl − kDEG·Ag     (x = AgEFF·Ag/Ag0)
dBAg1 = kon·C·Ag − koff·BAg1        − kon·x·BAg1 + 2·koff·Bdbl − kTMD·BAg1
dBdbl = kon·x·BAg1                  − 2·koff·Bdbl              − 2·kTMD·Bdbl
KSYN  = Ag0·kDEG ;  census Ag_tot = Ag + BAg1 + 2·Bdbl conserved (turnover/internalization aside)
```

**Backward-Euler, receptors-as-states.** The step is solved implicitly: a 3×3 linear generator M is
built per cell, and each sub-step solves `(I − h·M)·y = y₀ + h·s`
(`kinetic_rhoden_percell.py:119–159`). This is **unconditionally stable and census-exact** — no mass
is manufactured even at very large AgEFF — so the number of sub-steps is set by the *slow* scales
(kDEG, kTMD), not pinned by the stiff crosslink: `nsub = clip(⌈(kDEG + 2·kTMD)·dt/0.25⌉, 1, 16)`
(`kinetic_rhoden_percell.py:107–108`). Biologically: receptors are conserved (they turn over and
internalize, but they are not conjured), and the solver stays stable across the full physiological
range of receptor densities.

**Singular-row guard (robustness, counted not hidden).** Because `A = I − h·M` has diagonal
`1 + h·(non-negative rates) ≥ 1`, it cannot be singular for finite inputs; a singular row therefore
means a non-finite rate leaked in (NaN/Inf drug conc, or Ag0 = 0). Such a cell is **frozen**
(y = y₀) rather than crashing the whole run — the physically-correct no-op for a fully depleted cell
— and the event is counted and warned, not silently propagated
(`kinetic_rhoden_percell.py:134–158`). This is documented because it was a real failure mode: one
bad cell out of ~10⁵ (elranatamab, TSIM = 24) once killed a 40-minute simulation with a LinAlgError
(`kinetic_rhoden_percell.py:136–138`).

### 9.4 From internalization flux to interstitial drug loss

The solve returns the step-averaged internalization flux `intern_flux = kTMD·(BAg1 + 2·Bdbl)`
(nM/day, on the synapse reaction volume) (`kinetic_rhoden_percell.py:98,161`). The graph step then
converts it to an actual interstitial drug removal (`wholebody_percell.py:182–187`):

```
intern_copies_cell = (intern_flux / NM_PER_COPY) · count_scale      # copies/cell/day
loss_nM_day        = intern_copies_cell / AVO / v_cell · 1e9        # interstitial nM/day removed
organ_sink         = Σ intern_copies_cell / AVO · 1e9               # nmol/day (the TMDD sink)
```

- **`count_scale`** lifts the *sampled* cells to the *physiological* population: it is
  `tot_copies/ag_sum` where `tot_copies = pool_nM · V_is · AVO/1e9`
  (`coupled_percell_pk.py:91–94,102`). Real per-cell receptor copies are kept unchanged (correct
  avidity for the bivalent case); only the **cell count** is scaled so the summed organ sink
  integrates to the physiological receptor pool (`coupled_percell_pk.py:85–96`). The receptor pools
  per organ come from `Rtot_wholebody_final.json` (`run_tce_pd_reval.py:45,130`).
- The loss is clamped to at most 90% of the local drug per step for positivity
  (`wholebody_percell.py:192`).

So the organ TMDD sink is genuinely **emergent**: it is the sum over cells of receptor-bound drug
being internalized, largest where target density and drug concentration are both high (tumor, and
T-cell-dense tissues for the CD3 arm). Nothing about TMDD is imposed
(`qsp_costim_window_v2.py:219–222`).

### 9.5 Per-antigen kint and kdeg — the target-property table

`kint` and `kdeg` are **target properties, never constants** (`run_tce_pd_reval.py:9–10`), read from
the curated `antigen_kinetics_table.json` (`run_tce_pd_reval.py:11`; the runtime handoff copy is
byte-identical to `params/antigen_kinetics_table.json`, verified this session). Each entry is
`[kint_perday, kdeg_perday, rationale]`. The **verified live values** for the antigens the screen
uses:

| Antigen | kint (/day) | kdeg (/day) | Rationale (source) |
|---|---|---|---|
| CD20 (MS4A1) | 0.02 | 0.05 | non-internalizing type-I marker |
| CD19 | 0.05 | 0.10 | stable B-lineage marker, slow internalization |
| BCMA (TNFRSF17) | 2.0 | 1.5 | rapid internalization; Lee 2016 *Br J Haematol* 174:911, PMID 27313079 |
| GPRC5D | 0.2 | 0.3 | internalizing orphan GPCR |
| FcRH5 (FCRL5) | 0.2 | 0.25 | internalizing surface Ig-family receptor |
| HER2 (ERBB2) | 0.17 | 0.10 | moderate RTK internalization/recycling |
| EGFR | 0.95 | 0.30 | receptor-mediated internalization; Lobet 2023, DOI 10.1007/s40262-023-01270-2 |
| CD38 | 0.10 | 0.20 | internalizing ectoenzyme |
| DLL3 | 0.3 | 0.5 | internalizing ADC/TCE target |
| PD-1 (PDCD1) | 0.05 | 0.5 | slow internalization, activated-T turnover |

The runner's accessors default to `kint = 0.15`, `kdeg = 0.5` /day for an antigen absent from the
table (`run_tce_pd_reval.py:14,19`). This is the single largest driver of *cross-target* TMDD
differences: a BCMA engager (kint 2.0/day) has a far stronger target sink than a CD20 engager
(kint 0.02/day) at the same affinity and density — the pharmacokinetic signature of a rapidly
internalizing versus a non-internalizing antigen.

**Per-molecule on/off rates.** `kon1`/`koff1` for the target arm come from the molecule's measured
kinetics when available (`eng_params_normalized.json`, merged at `run_tce_pd_reval.py:101–124`),
else the standard mAb fallback `kon = 10⁵` /M/s with `koff = kon·KD` (KD-consistent)
(`run_tce_pd_reval.py:25–33`). These are converted to per-day/per-nM at construction
(`coupled_percell_pk.py:69–71`, §1). If the kinetic params are None, the graph step falls back to the
**QSS conversion** `percell_rhoden_qss_Cvec` — monovalent `S = Rtot·C/(C+KD)`, or the closed-form
bivalent ternary root for n_arm ≥ 2 (`wholebody_percell.py:42–52,188–191`) — so a run is always
well-defined even without measured rates.

---

## 10. The well-mixed layer's emergent TMDD (for completeness)

The Layer-1 well-mixed PBPK (`qsp_costim_window_v2.py`) implements the same physics with one lumped
interstitial compartment per organ. Its plasma/interstitium/lymph ODEs
(`qsp_costim_window_v2.py:904–934`):

```
C_pl     = A_pl / V_pl                                            # nM (V_pl = 3.1 L)
Cis      = A_is / PB.Vis                                          # per-organ interstitial conc
J_extrav = k_dist·PB.L·(1 − PB.sigV)·C_pl                         # plasma → each interstitium
J_return = k_dist·PB.L·(1 − sigL)·Cis                             # each interstitium → lymph
tmdd_sink = (kint_mono·Cb_mono + kint_cplx·Cb_bridge)·PB.Vis      # saturable target-mediated elim.
k_cat    = CLup·(1 − eff_fFcRn()) + k_renal()                     # linear catabolism
dA_pl = infusion − k_cat·A_pl − ΣJ_extrav + k_lymph_return·A_ly
dA_is = J_extrav − J_return − tmdd_sink
dA_ly = ΣJ_return − k_lymph_return·A_ly
```

The topology is identical to the per-cell core; the difference is that the interstitial binding is a
lumped bivalent-binding call `bivalent_binding(Cis, R_CD3, R_TAA, p)` returning monovalent and
bridged complex pools (`qsp_costim_window_v2.py:918`), internalized at `kint_mono = 0.25` /day and
`kint_cplx = 0.90` /day (`qsp_costim_window_v2.py:287–288`). CD3 capacity scales with tissue T-cell
density `R_CD3 = Rcap_CD3·PB.tcell` (`qsp_costim_window_v2.py:916`, `Rcap_CD3 = 2.0`), TAA capacity
is tumor-only `R_TAA = Rcap_TAA·PB.taa` (`Rcap_TAA = 6.0`). Because the bound complex saturates at
the finite receptor pool, fractional clearance falls as dose rises — the emergent TMDD nonlinearity
(`qsp_costim_window_v2.py:919–922`).

This layer is the byte-frozen abstract window-scoring engine (`pd_model_config.py:9–11`); its PK
numbers are comparative, not clinical predictions (§12). The per-cell core (§4–9) is the one that
produces the spatially-resolved, clinically-anchored PK the validation track uses.

---

## 11. Numerical integration and time-stepping

- **Per-cell core:** fixed-step explicit Euler on the systemic pools, implicit (LU) diffusion on the
  graph, backward-Euler on the per-cell binding. Default `dt = 0.01` day
  (`coupled_percell_pk.py:65`), overridable by env `PD_DT` (`coupled_percell_pk.py:74`); the
  production runner sets `dt = 0.02` (`run_tce_pd_reval.py:143`). Organs are stepped in parallel
  across a thread pool (`coupled_percell_pk.py:106,131`). PD is advanced every `pd_every = 3` steps
  (`run_tce_pd_reval.py:173`).
- **Well-mixed layer:** `scipy.integrate.solve_ivp` with method **LSODA**, `rtol = 1e-6`,
  `atol = 1e-8`, `max_step = 0.5` day, dense output (`qsp_costim_window_v2.py:1229–1232`). IV bolus
  enters as the initial plasma amount `y0[I_PL] = dose` (`qsp_costim_window_v2.py:1228`).
- **Simulation horizon** is per-run; for IL-6/CRS-peak validation the runner truncates the horizon
  and explicitly warns that steady-state PK must not be quoted from a truncated run
  (`run_tce_pd_reval.py:100`).

---

## 12. Calibration philosophy and honest scope

The PBPK/PD layer is **shape-validated** against literature PK/PD (Betts/Schropp/Hosseini) so
absolute time-courses are credible, but **per-arm window numbers remain comparative, not clinical
predictions** (`qsp_costim_window_v2.py:48–55`). The PK constants that are calibrated (and to what)
are stated at their definitions:

- `CLup = 0.3503` → mosunetuzumab terminal t₁/₂ 16.1 d (LUNSUMIO label)
  (`qsp_costim_window_v2.py:182–188`).
- `sigL = 0.85` → V_ss/V_c ≈ 2.1 (class-typical mAb) (`qsp_costim_window_v2.py:170–175`).
- `k_dist = 3.0` → ~1.3-day α-phase vs. pembrolizumab (`qsp_costim_window_v2.py:176–180`).
- `k_renal_max = 8.70` → blinatumomab ~2.1 h terminal t₁/₂ (Blincyto label)
  (`qsp_costim_window_v2.py:205–208`).

Everything else (Q, L = Q/500, σ_V, V, fV, fIS, cell densities) is **fixed reference physiology**,
not fitted (`qsp_costim_window_v2.py:75–78`).

---

## 13. Limitations and disclosures

### 13.1 Static R_costim — mandatory disclosure

**The costim-arm receptor number is set once at initialization from resting copy numbers and read
unchanged at every step.** The model captures conditionality through binding *geometry* (the cis
p_cis gate, §8.4) but **not** activation-induced upregulation of the costim receptor. Consequently a
resting-copy ranking **under-rates 4-1BB / OX40 / ICOS** — receptors that are strongly upregulated on
activated T cells — and yields a spurious ordering in which CD2 (constitutively high at rest) appears
to win on receptor availability. This is a known limitation of the current model state and must be
read alongside any receptor-availability comparison (per the submission manifest's
`static_R_costim_LIMITATION`).

The engine does contain a `costim_induction.py` module, but it is **default OFF** and active only for
costim-armed sweep constructs; the canonical PK/PD path documented here runs with static R_costim.
Activation-induced R_costim is only to be asserted if wired with literature-sourced induction
kinetics and a version tag; otherwise the static version is canonical.

### 13.2 Explicitly PROVISIONAL — not asserted here

The claim that OX40/GITR arms are **net-negative on killing in Treg-rich settings** is **PROVISIONAL**
— it is not artifact-backed and is structurally impossible under static R_costim (a resting-copy model
cannot represent the activation-dependent reversal). It is omitted from this document's conclusions and
flagged here only so the omission is deliberate.

### 13.3 Structural scope notes

- **Renal/BiTE branch:** the MW-Hill renal sieve and `has_fc` accessor are implemented on the
  well-mixed `PBPK` dataclass; the per-cell production runs documented here are Fc-bearing IgG-class
  TCEs (fFcRn 0.70–0.89) and do not exercise a nonzero `k_renal` (§7.3).
- **FcRn as QSS fraction:** salvage is a recycled-fraction reduction of catabolism, not an explicit
  endosomal on/off ODE — a deliberate robustness choice (`qsp_costim_window_v2.py:69–71`); it
  reproduces the long IgG half-life but does not resolve endosomal dynamics or FcRn-competition
  (Fc-engineering, high-dose FcRn saturation) explicitly.
- **Vascular QSS:** organ vascular pools are algebraic each step (Q/V_v ≫ pharmacology rate),
  physically exact for perfusion-limited spaces (`wholebody_percell.py:163–168`); this removes
  artificial stiffness but means intra-vascular transients faster than a step are not resolved.
- **Comparative window numbers:** per-arm window magnitudes are relative rankings anchored to
  reference molecules, not absolute clinical PK predictions (§12).

---

## 14. Parameter reference (verified live values)

All values quoted from live source with `file:line`; SI→internal conversions per §1.

| Parameter | Value | Units | Source (`file:line`) |
|---|---|---|---|
| `V_PLASMA` (V_pl, V_c) | 3.1 | L | `qsp…:362`; `ck_pk:65` |
| `V_LYMPH` (V_ly) | 2.6 | L | `qsp…:363`; `ck_pk:65` |
| `_PLASMA_CO` | 5000 | L/day | `qsp…:100` |
| `_LYMPH_RATIO` (L = Q/500) | 1/500 | — | `qsp…:101,141` |
| `sigL` (σ_L) | 0.85 | — | `qsp…:170`; `ck_pk:64` |
| `k_dist` | 3.0 | — | `qsp…:176`; `ck_pk:64` |
| `fFcRn` (platform default) | 0.90 | — | `qsp…:181` |
| `fFcRn` (IgG1/4 clinical set) | 0.89 | — | `run…:77–82,86–88` |
| `fFcRn` (elranatamab, IgG2) | 0.70 | — | `run…:80` |
| `CLup` | 0.3503 | /day | `qsp…:182`; `ck_pk:65` |
| `k_cat` = CLup·(1−fFcRn)+k_renal | derived | /day | `ck_pk:76`; `qsp…:930` |
| `Kdeg` (unsalvaged endosomal) | 26.0 | /day | `qsp…:189` |
| `k_lymph_return` | 24.0 | /day | `qsp…:190`; `ck_pk:65` |
| `k_renal_max` | 8.70 | /day | `qsp…:205` |
| `mw50_renal` | 69 | kDa | `qsp…:209` |
| `hill_renal` | 10 | — | `qsp…:210` |
| `mw_kda` (IgG default) | 146.9 | kDa | `qsp…:198` |
| SC bioavailability `F_sc` | 0.6 | — | `ck_pk:108`; `run…:175` |
| SC absorption `ka_sc` | 0.25 | /day | `ck_pk:108`; `run…:175` |
| IV infusion duration `iv_inf_h` | 2.0 | h | `run…:176` |
| `NM_PER_COPY` | 2.335×10⁻⁵ | nM/copy | `wb_pc:41` |
| cell radius (Rhoden geo) | 8.0 | µm | `krp:37`; `wb_pc:33` |
| arm reach (span default) | 12.5 | nm | `krp:37`; `mab:27` |
| cis gap-match / tol | 12.5 / 8.0 | nm | `mab:39` |
| synapse cleft min / max | 13 / 40 | nm | `ks:41–42`; `mab:48` |
| diffusivity `D_um2s` | 10.0 | µm²/s | `wb_pc:77` |
| ECM hindrance `alpha_D` | 3.0 | — | `wb_pc:77` |
| kNN neighbours `k` | 6 | — | `wb_pc:77` |
| lymphatic catchment `lam` | 100 | µm | `wb_pc:119` |
| well-mixed `kint_mono` / `kint_cplx` | 0.25 / 0.90 | /day | `qsp…:287–288` |
| default per-cell dt | 0.01 (prod 0.02) | day | `ck_pk:65`; `run…:143` |

**Per-antigen kint / kdeg (/day):** CD20 0.02/0.05 · CD19 0.05/0.10 · BCMA 2.0/1.5 · GPRC5D 0.2/0.3 ·
FcRH5 0.2/0.25 · HER2 0.17/0.10 · EGFR 0.95/0.30 · CD38 0.10/0.20 · DLL3 0.3/0.5 · PD-1 0.05/0.5;
default 0.15/0.5 (`params/antigen_kinetics_table.json`; accessors `run…:14,19`).

File abbreviations: `qsp…` = `qsp_costim_window_v2.py`, `ck_pk` = `coupled_percell_pk.py`,
`ck_pd` = `coupled_percell_pd.py`, `wb_pc` = `wholebody_percell.py`,
`krp` = `kinetic_rhoden_percell.py`, `mab` = `multiarm_binding.py`, `ks` = `kinetic_synapse.py`,
`run…` = `run_tce_pd_reval.py`.

---

## 15. Summary — the molecule's life in one paragraph

A dose enters plasma (IV, as a ~2 h mass-exact infusion) or a SC depot (first-order absorption `ka`,
bioavailability `F_sc`). From the 3.1 L plasma pool it perfuses 15 organs + tumor; each organ's
vascular space equilibrates instantly (QSS) and antibody crosses the blood-endothelial cells by
2-pore **convection**, gated by the vascular reflection σ_V and carried at the lymph-flow rate
`L = Q/500`. In the interstitium it diffuses on a per-cell spatial graph (ECM-hindered) and **binds
its target on each cell** via the Rhoden bivalent kinetic scheme — receptors as live states with
turnover `KSYN = Ag0·kDEG`, avidity crosslink when the arm is bivalent, and internalization `kTMD` of
every bound species. That internalization **is** the target-mediated clearance (TMDD), emergent and
saturable, strongest where target density and drug are both high. Unbound interstitial drug drains at
the lymphatic-endothelial cells (reflection σ_L), collects in a fast lymph pool, and **recirculates**
to plasma (`k_lymph_return`). In plasma the molecule is salvaged by FcRn (recycled fraction `fFcRn`)
and the unsalvaged remainder is catabolized linearly at `k_cat = CLup·(1−fFcRn) + k_renal` — the term
that sets the terminal half-life (≈18–20 d for an IgG1/4 TCE, ≈6.6 d for the IgG2 elranatamab, ≈2 d
and falling for a no-Fc BiTE). **The very same per-cell binding solve that clears the drug is the
engagement that drives killing and cytokine release in the PD arm** — one occupancy, two consequences
— which is the structural guarantee that the counter-screen's efficacy and toxicity axes are read
from the same physical event that governs the pharmacokinetics.

