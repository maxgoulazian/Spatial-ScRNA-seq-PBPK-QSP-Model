---
title: "T6 — Per-cell PD: trimer formation, killing, Treg damping & cytokine output"
subsystem: T6 (per-cell pharmacodynamics; whole-body, every organ)
model: costim_engager_counterscreen (whole-body single-cell PBPK-PD)
live_source_files:
  - engine/wholebody_pd.py   (503 lines, read in full this task)
  - engine/pd_model_config.py (75 lines, read in full this task)
supporting_files_read_for_provenance:
  - params/kinetic_calib.json  (the file that SETS k_death, k_hit, IL6_SCALE_kin)
  - engine/kinetic_synapse.py  (k_hit consumer; documented in its own subsystem doc)
  - engine/run_tce_pd_reval.py (top-level driver — the caller that supplies k_death)
  - engine/coupled_percell_pd.py (the caller that constructs OrganPD)
date: 2026-07-13
generated_by: workflow-subagent T6
verified: 2026-07-13 adversarial re-verification against live code (every equation, line number, parameter
  value and citation re-checked; prozone table, V_syn, serial ceiling and IL6_SCALE arithmetic re-RUN).
  Corrections applied in-place — see the changelog at the end of §6.
repo_state: NOT a git repository (confirmed: `git rev-parse` → "not a git repository"; no HEAD to pin) —
  line numbers are against the working tree as read 2026-07-13 (`wholebody_pd.py` = 503 lines,
  `pd_model_config.py` = 75 lines, both confirmed by `wc -l`)
---

up:: [[00_INVENTORY_AND_MAP]]
tags:: #model-doc #pd #tce #cytokine

> [!danger]+ LIVE CODE DEFECT — read §5.12 before quoting any COSTIM result from this subsystem
> `costim_ind.step()` (activation-induced costim receptor density) is called **only in the QSS path**
> (`wholebody_pd.py:334`). It is **never called in `_step_kinetic`**, which is the **CANONICAL** engine
> (`pd_model_config.PD_ENGINE = "kinetic"`). In every canonical run, `R_costim` therefore stays frozen at the
> **resting** density, and the inducible arms (4-1BB / OX40 / ICOS / GITR) are read at ≈0 — exactly the
> "confidently-wrong ranking ('CD2 beats 4-1BB')" failure the induction code was written to prevent
> (`:322–331`). Clinical validation molecules (`costim_arm=None`) are unaffected. **§5.12 has the one-line
> fix. Until it lands, no costim-arm ranking from the kinetic engine is trustworthy.**

> [!warning]+ PROVENANCE HEALTH WARNING — read §3.5 before quoting any number from this subsystem
> This module contains **three constants that are fitted to a clinical anchor which the 2026-07-13 provenance
> audit found to have NO SOURCE** (`CYTO_IL6_CLINICAL_ANCHOR_PGML = 570.0`, `IL6_SCALE_kin = 0.05473`, and the
> `clin_il6: 570.0` entry inside `params/kinetic_calib.json`). Two of those three are now **dead code / dead
> data** — the fitted IL-6 pathway was deleted from the driver on 2026-07-13 and replaced by the mechanistic
> myeloid emitter — but **the constants are still physically present in the live module** and a reader who
> greps for "570" will find them. They are documented here as CONTAMINATED-AND-ORPHANED so that nobody
> re-wires them. See §3.5 and §5.2.

---

## 1. PURPOSE & DATAFLOW POSITION

**T6 is the per-cell pharmacodynamic layer.** It answers one question, once per PD step, for every T cell in
every organ: *given the drug concentration at THIS T cell, the CD3 copies on THIS T cell, and the antigen
copies on the target cells physically within reach of THIS T cell — how much bridged trimer forms, how much
lethal signal does that trimer deliver to those specific neighbours, and how much cytokine does that T cell
drive?*

Everything in the module is **per agent**. There is no tumour compartment, no effector pool, no "T-cell
density × antigen density" product. One agent = one barcode from the scRNA-seq-derived cell table; each T
cell forms its own trimer from its own local free drug (delivered by the transport grid), its own `CD3E`
expression, and its own spatial neighbourhood of antigen-bearing cells (module docstring,
`engine/wholebody_pd.py:1–15`).

### Where it sits in the life of the molecule

```
   PK (plasma)                     TRANSPORT (per-cell grid)                 T6 — PER-CELL PD
   ───────────                     ─────────────────────────                 ────────────────
   dose → plasma C  ──────────►  per-cell free drug C_i (nM)  ──────────►  Schropp ternary trimer  Cb_i
                                  wholebody_percell.py                       (EQ-2)
   scRNA-seq agent table                                                          │
   ───────────────────                                                            ├──► kill hazard on the
   CD3E copies      ─────────────────────────────────────────────────────────►    │    T cell's NEIGHBOURS
   TAA copies       ─────────────────────────────────────────────────────────►    │    (EQ-8, EQ-9, EQ-10)
   costim copies    ──────────► costim_induction (R(t))  ────────────────────►    │
   x,y coordinates  ──────────► cKDTree synapse graph (EQ-4/5/6)  ──────────►     ├──► cytokine rate
   cell-type labels ──────────► T / CD8 / Treg / target masks                     │    (EQ-13/14/15)
                                                                                  └──► engaged fraction p_eng
                                                                                       ──► myeloid_il6 (IL-6)
```

**Consumes (inputs):**
- `C_percell` — per-cell free drug (nM), handed in every step by `CoupledPerCellPD` from the per-cell
  transport grid (`step(self, C_percell, dt, k_death=1.0)`, `engine/wholebody_pd.py:297`).
- `R_CD3_percell`, `R_TAA_percell` — per-cell receptor **copy numbers** from the agent table
  (`__init__`, `:110`; converted to nM at `:168–169`).
- `R_costim_percell` — per-cell costim-receptor copies (optional; `:110–111`).
- `x`, `y`, `labs` — agent coordinates and cell-type labels (`:110`).
- `k_death` — the trimer→death potency constant, passed **in from the caller** every step (`:297`), sourced
  from `pd_model_config.K_DEATH` (`engine/pd_model_config.py:47,59,63`).
- `k_hit_perday` — serial-killing ceiling, passed at construction (`:117`), sourced from
  `pd_model_config.KINETIC['k_hit_perday']` (`engine/pd_model_config.py:41,60`).

**Feeds (outputs):**
- `kill_hazard[i]` — cumulative per-cell hazard; survival = `exp(-hazard)` (`:362`, `:385`, `:499`).
- `cyto` / `cyto_rate` — cumulative and instantaneous organ cytokine, per species in `CYTO_HIER` (`:374–376`).
- `il6_prod_pg_hr` — this organ's myeloid IL-6 production, computed by the `myeloid_il6` subsystem from the
  per-T **engaged fraction** T6 supplies (`:380`, `:486`). This is the ONLY live IL-6 channel.
- `summary()` — O(1) organ readout of values `step()` already cached (`:493–496`). `summary_full()` (`:498–503`) is
  **not** O(1): it recomputes `exp(−H)` over every cell (`:499`).

**Two mutually-exclusive kill laws live in this module**, selected by `pd_kinetics`:
- `pd_kinetics=False` → **QSS path** (`step`, `:297–387`): Schropp ternary equilibrium, evaluated fresh each
  step. Equilibrium magnitude only — no dwell time, no serial-killing cycle.
- `pd_kinetics=True` → **KINETIC path** (`_step_kinetic`, `:438–491`): T6 becomes a thin wrapper that routes
  the per-cell drug into `kinetic_synapse.KineticSynapse` (its own subsystem) and applies T6's Treg/costim
  modifiers to the hazard that comes back. **This is the CANONICAL engine** (`pd_model_config.py:18,25,29`).

The selection is made in exactly one place — `pd_model_config.PD_ENGINE = "kinetic"`
(`engine/pd_model_config.py:29`) — which is the entire reason `pd_model_config.py` exists (its docstring,
`:3–4`: *"lets make sure to separate out and make a 1 final model so we dont get confused on other runs later"*).

> **Read §5.12 before trusting any costim ranking.** The two paths are **not** the same code. The
> activation-induced costim receptor density (EQ-11, `:332–334`) is stepped **only in the QSS path**. In the
> **canonical kinetic path it is never stepped**, so `R_costim` stays at its resting value for the whole run.

---

## 2. GOVERNING EQUATIONS

Notation: `C` = free drug at the cell (nM); `R_A` = CD3 arm capacity (nM); `R_B` = TAA arm capacity (nM);
`KD1` = `KD_CD3` (nM); `KD2` = `KD_TAA` (nM); `Cb` = bridged ternary trimer (nM); `dt` = step (days);
`N_A` = 6.02214076e23 (`:20`).

---

### EQ-1 — Absolute copies → local synapse concentration (`engine/wholebody_pd.py:83`, applied `:168–169`)

```
NM_PER_COPY = 6.0 / 257000                       = 2.3346303e-5  nM per receptor copy
R_CD3_i [nM] = R_CD3_raw_i [copies] · NM_PER_COPY        (:168)
R_TAA_i [nM] = R_TAA_raw_i [copies] · NM_PER_COPY        (:169)
```

**Biological meaning.** A receptor is not a concentration; it is a count on a membrane. To put it into a
mass-action equilibrium you must declare a **reaction volume** — the volume of the immunological synapse in
which the drug, the CD3 and the TAA actually meet. `NM_PER_COPY` *is* that declaration: it is the reciprocal
of a synapse reaction volume. Inverting it (RUN-verified this task):

```
V_syn = 1e9 / (N_A · NM_PER_COPY) = 7.113e-11 L = 71.1 pL          [RUN-verified, this task]
```

**Mechanistic rationale, and the alternative that was rejected.** The rejected alternative is spelled out in
the code comment (`:75–83`, `:164–167`): the module previously **normalised each cell's receptor count by the
per-organ mean**. That erases absolute abundance — a low-copy target (BCMA ~11k copies) would present the same
normalised "1.0" as a high-copy target (CD20 ~95k copies), and would therefore saturate the trimer identically.
For a **counterscreen whose entire job is to rank constructs against targets of different abundance, that is
fatal**. The literal-absolute conversion preserves abundance, so (RUN-verified this task):

| receptor | copies | → nM | vs KD |
|---|---|---|---|
| CD3 (T cell) | ~92,000 | **2.148 nM** | ≪ KD_CD3 = 40 nM → **linear regime**, as CD3 physically is |
| CEACAM5 (the anchor) | 257,000 | **6.000 nM** | — (this is the pin, by construction) |
| CD20 (mosunetuzumab) | ~95,000 | **2.218 nM** | > KD_TAA → strong bridge |
| BCMA (elranatamab) | ~11,000 | **0.257 nM** | ≲ KD_TAA → weak bridge, exposure-dependent depletion |

That BCMA/CD20 spread is exactly the counterscreen signal, and it exists **only** because the conversion is
absolute. **Units:** copies [dimensionless] × nM·copy⁻¹ → nM.

> **The pin is circular and must be stated as such.** `NM_PER_COPY = 6.0/257000` is defined by forcing the
> *old tumour model's* `Rcap_TAA = 6.0 nM` to land on CEACAM5's 257,000 copies/cell (`:76–78`). `Rcap_TAA=6.0`
> was itself a **receptor-capacity scale in a compartmental QSP model** (`engine/qsp_costim_window_v2.py:227`),
> not a measured synapse volume. So the 71 pL synapse is **[DERIVED from a fitted compartmental scale]**, not
> measured. See §3.5 and §5.1.

---

### EQ-2 — Schropp ternary equilibrium (the bridged trimer) — `ternary_equilibrium`, `engine/wholebody_pd.py:85–105`

Written exactly as coded (equation numbers are the module's own references to Schropp 2019, `:86–87`):

```
aKK  = α · KD1 · KD2                                                                    (:94)

a    = (1 + C/KD2) · C / aKK                                             [Schropp Eq.28]  (:98)
b    = C·(R_A − R_B)/aKK + (1 + C/KD2)·(1 + C/KD1)                       [Schropp Eq.29]  (:99)
d    = −R_B · (1 + C/KD1)                                                [Schropp Eq.30]  (:100)

R_Bf = ( −b + sqrt( max(b² − 4·a·d, 0) ) ) / (2a)                        [Schropp Eq.27]  (:102)
R_Af = R_A / ( 1 + C/KD1 + R_Bf·C/aKK )                                  [Schropp Eq.26]  (:103)

Cb   = C · R_Af · R_Bf / aKK                                             [Schropp Eq.33]  (:104)
```

Guards: `C` clamped ≥0 (`:92`); the solve runs only on the mask `m = (C>0) & (R_A>0) & (R_B>0)` (`:96`), all
other cells get `Cb = 0` (`:95`); the discriminant is floored at 0 (`:102`).

**Biological meaning.** A T-cell engager must bridge: one arm on CD3 on the T cell, one arm on the TAA on the
target cell, **simultaneously**. `Cb` is the concentration of that bridged trimer — the only species that
delivers a lethal signal. The two binary complexes (drug·CD3 and drug·TAA) are *dead ends*: they consume
receptor and consume drug but do not kill.

**Mechanistic rationale — and this is the load-bearing design decision of the whole subsystem.** The quadratic
is a **quasi-equilibrium receptor-conservation solve**: it computes the **free** receptor concentrations
`R_Af`, `R_Bf` from the *total* receptor concentrations `R_A`, `R_B` **before** it forms the trimer. I verified
the algebra by hand this task — the quadratic in `R_Bf` is precisely what falls out of the two conservation
statements

```
R_A,tot = R_Af·(1 + C/KD1) + C·R_Af·R_Bf/aKK        (free + binary + ternary)
R_B,tot = R_Bf·(1 + C/KD2) + C·R_Af·R_Bf/aKK
```
after eliminating `R_Af`. The coefficients `a`, `b`, `d` in the code are algebraically identical to that
elimination. **This is what makes the prozone/hook effect EMERGENT.** At high drug, both arms of the drug get
saturated as *separate binary* complexes on *separate* molecules; there is no free receptor left for a *bridged*
molecule to find, and the trimer **collapses**.

The **rejected alternative is named in the code** (`:90–91`): the reduced-linear form
`Cb = C·R_A·R_B/(KD1·KD2)`, which has **no free-receptor depletion** and therefore **no prozone** — it rises
monotonically with drug forever and *"overstated kill/cytokine at high drug."* That shortcut is gone.

**RUN-VERIFIED PROZONE (this task, live code, CD3=92k copies, TAA=257k copies, KD_CD3=40, KD_TAA=1.45):**

| C (nM) | 0.01 | 0.1 | 1 | 3.16 | **10** | 31.6 | 100 | 1000 | 1e4 |
|---|---|---|---|---|---|---|---|---|---|
| Cb (nM) | 2.20e-3 | 2.05e-2 | 1.19e-1 | 1.82e-1 | **1.977e-1 (grid max)** | 1.55e-1 | 8.59e-2 | 1.23e-2 | 1.28e-3 |

`C = 10 nM` is the maximum **of this sampled grid**, not the true maximum. A 400-point log scan of the live
function (this task) puts the **true peak at C = 7.58 nM, Cb = 0.1994 nM**. The hook is real, it peaks
between 7 and 8 nM, and it decays as ~1/C. **Nobody imposed it. It is a consequence of receptor
conservation.** This is the single strongest "emergence" claim T6 can make (§4).

**Units.** `C, R_A, R_B, KD1, KD2, Cb` all in nM. `α` dimensionless. `aKK` in nM².

**Two honest caveats on this equation (see §5):**
1. **α (cooperativity) defaults to 1.0** (`:85`) and is **never passed** by any caller — the model assumes no
   cooperativity between the two arms.
2. **Free drug `C` is NOT depleted by binding here.** The conservation solve conserves *receptors* but treats
   `C` as a fixed free-ligand reservoir. Drug consumption is handled upstream (the transport/antigen-sink
   layer); within a PD step the trimer does not draw down `C`.

---

### EQ-3 — Cell classification masks (`engine/wholebody_pd.py:155–160`)

```
low      = lowercase(labels)                                                        (:155)
is_T     = ('t cell' ∨ 't_cell' ∨ 'cd8' ∨ 'cd4' ∨ 'regulatory t' ∨ 'nk t') ∈ low    (:156–157)
is_CD8   = 'cd8' ∈ low                                                              (:158)
is_treg  = ('regulatory t' ∨ 'treg') ∈ low                                          (:159)
is_target = (R_TAA_percell > 0)                                                     (:160)
```

**Biological meaning.** Who can kill (`is_T`), who gets the effector boost (`is_CD8`), who suppresses
(`is_treg`), and who is killable (`is_target`).

**Mechanistic rationale.** `is_target` is defined **by antigen, not by cell type** (`:160`, comment: *"any
antigen-bearing cell is killable"*). This is the correct and important choice for a **counterscreen**: on-target
off-tumour killing of healthy antigen-positive cells is not a special case bolted on — it is the *same*
equation. A healthy hepatocyte expressing the TAA is killed by exactly the mechanism that kills a tumour cell.

**Units.** Boolean masks.

**Two structural consequences a reviewer will notice (both flagged in §5):**
- `is_treg ⊂ is_T` — **Tregs are effectors too**. They enter `Tidx`, form trimers (EQ-2), and contribute to
  `kill_T` (EQ-8) at baseline gain, while *simultaneously* suppressing their neighbourhood (EQ-7/EQ-8).
- Classification is **substring matching on free-text labels**. `'cd8'` will match any label containing that
  substring. This is fragile against label-vocabulary changes in the agent table.

---

### EQ-4 — Synapse neighbourhood (who can I kill?) — `engine/wholebody_pd.py:257–258`

```
tree_t   = cKDTree( (x,y) of target cells )                                (:257)
syn_nb_j = { targets t : ‖ (x,y)_Tj − (x,y)_t ‖₂  ≤  R_SYN_UM = 30 µm }    (:258)
```

**Biological meaning.** A T cell can only kill a cell it can physically reach. `R_SYN_UM = 30 µm` is the
synapse **reach** — the radius within which a T cell is deemed able to form an immunological synapse with a
target during the step.

**Mechanistic rationale.** The alternative — a well-mixed organ where every T cell sees every target — is
exactly what this whole model exists to replace. A finite reach is what makes killing depend on **where the T
cells are** (bone-marrow vs spleen vs skin), which is the source of the model's organ-to-organ selectivity.
The 30 µm reach is the "validated tumor value" (`:36`) — it is a **spatial-model parameter, not a measurement**
(§3).

**Units.** µm.

---

### EQ-5 — Antigen-weighted apportionment matrix (how do I split the kill?) — `engine/wholebody_pd.py:263–269`

```
W[j,t]      = 1  if target t ∈ syn_nb_j, else 0            (sparse T×n_target incidence)    (:263–264)
Wt[j,t]     = W[j,t] · R_TAA[t]                            (edge weighted by target antigen) (:267)
rs[j]       = Σ_t Wt[j,t]   (with rs[j]=0 → 1 to avoid /0)                                   (:268)
Wt_norm     = diag(1/rs) · Wt                              (row-stochastic)                  (:269)
```

**Biological meaning.** A T cell in reach of five targets does not kill all five equally: it preferentially
engages the ones presenting the most antigen. `Wt_norm[j,t]` is the fraction of T-cell *j*'s lethal output
that lands on target *t*, proportional to *t*'s antigen density.

**Mechanistic rationale.** Row-normalisation enforces a **conservation of lethal output**: each T cell
distributes exactly 1.0 of its kill signal across its neighbours (never more), so a T cell in a dense
neighbourhood spreads its output thin rather than multiplying it. The alternative (unnormalised weights)
would let a single T cell in a dense field deliver unbounded total kill — a well-known artefact of naive ABM
kill rules.

**Units.** `W` dimensionless; `Wt` in nM; `Wt_norm` dimensionless (row-stochastic).

---

### EQ-6 — Per-T mean neighbourhood antigen (the TAA arm of the trimer) — `engine/wholebody_pd.py:270`

```
syn_TAA_mean[j] = ( Σ_t W[j,t]·R_TAA[t] ) / max( Σ_t W[j,t], 1 )       [nM]     (:270)
```

**Biological meaning.** The TAA arm concentration that T cell *j* actually experiences is the **mean antigen
density of the target cells within its synapse reach** — not the organ mean, not the tumour mean. A T cell
sitting in an antigen-poor pocket forms a weak bridge even at saturating drug.

**Mechanistic rationale.** This is what makes the trimer (EQ-2) **spatially heterogeneous**. `R_B` in the
Schropp solve is `syn_TAA_mean`, i.e. a *per-T-cell*, *neighbourhood-derived* quantity (used at `:309`, and
re-used as the ceff density basis for the kinetic engine at `:240`). Spatial antigen heterogeneity therefore
propagates directly into trimer, kill and cytokine, with no extra machinery.

**Units.** nM.

---

### EQ-7 — Treg neighbourhood count — `engine/wholebody_pd.py:275–277`

```
n_treg[j] = | { Tregs r : ‖ (x,y)_Tj − (x,y)_r ‖₂ ≤ R_TREG_UM = 50 µm } |     (:277)
```

**Biological meaning.** Regulatory T cells suppress effector function over a local radius (contact + secreted
IL-10/TGF-β/IL-2 sink). `n_treg[j]` is how many Tregs are close enough to suppress T cell *j*.

**Mechanistic rationale.** Suppression is made **spatial and countable** rather than a global immunosuppression
scalar. The suppression a T cell feels emerges from the local Treg census of the scRNA-seq-derived agent table
— so a Treg-rich organ (or a Treg-rich tumour region) damps killing there and nowhere else.

**Units.** dimensionless count; `R_TREG_UM` in µm.

---

### EQ-8 — Treg-damped per-T-cell kill (QSS) — `engine/wholebody_pd.py:350–351`

```
Cb'_j    = Cb_j · g_eff_j                                                             (:350)
kill_T_j = Cb'_j / ( 1 + TREG_K · n_treg_j · (1 + supp_extra) )      [nM]              (:351)

with TREG_K = 0.25                                                                     (:72)
```

**Biological meaning.** The lethal drive a T cell delivers is its bridged trimer, boosted by its costimulatory
effector program (`g_eff`, EQ-11), and **divided down** by the suppressive Tregs in its neighbourhood.

**Mechanistic rationale.** The **hyperbolic** (`1/(1+K·n)`) form, rather than exponential or linear-subtractive,
is the standard saturable-inhibition form: one Treg cuts the kill to 80% (`1/(1+0.25)`), four Tregs to 50%,
and it never goes negative (a linear `Cb·(1 − K·n_treg)` form would). `supp_extra` (EQ-12) makes the *strength*
of each Treg itself a function of that Treg's costim-driven suppression program — so an agonist that
inadvertently boosts Tregs raises the denominator, which is the tox axis the counterscreen is looking for.

**Units.** `Cb` nM; `TREG_K` per-Treg (dimensionless); `n_treg` count; `kill_T` nM.

---

### EQ-9 — Vectorised apportionment to targets — `engine/wholebody_pd.py:354–355`

```
dkill_tgt = Wt_normᵀ · kill_T                  [nM, length n_target]     (:354)
dkill[tgtidx] = dkill_tgt                                                 (:355)
```

**Biological meaning.** Each target cell's instantaneous kill propensity is the antigen-weighted sum of the
lethal drive of every T cell that can reach it. A target surrounded by three engaged T cells accumulates
hazard three times as fast.

**Mechanistic rationale.** The transpose of the row-stochastic matrix is what makes killing **many-to-many and
conservative**: T cells share their output across targets (EQ-5), targets sum the output they receive. This is
the ABM kill rule, expressed as one sparse matrix-vector product (which is also why it is fast enough to run
in every organ).

**Units.** nM.

---

### EQ-10 — Hazard accumulation and survival — `engine/wholebody_pd.py:362`, `:385`

```
H_t  ←  H_t + dt · k_death · dkill_t                                     (:362)
S_t  =  exp( −H_t )                                                      (:385, :468, :499)
kill_frac = mean over targets of ( 1 − exp(−H_t) )                       (:385)
```

**Biological meaning.** Cell death is a **stochastic hazard**, not a threshold. A target accumulates hazard at
a rate proportional to the lethal drive it is receiving, and its survival probability decays exponentially in
the accumulated hazard. The reported organ `kill_frac` is the mean death probability across that organ's
antigen-positive cells.

**Mechanistic rationale (stated explicitly in the code, `:356–361`).** The rejected alternative is a
**fixed/saturating hazard** — "if engaged, die at rate k". That makes depletion depend only on *time*, so any
dose above a threshold gives the same depletion, and the model loses all dose-response. Integrating
`dt·k_death·dkill` instead makes cumulative depletion depend on **both drug level (via `dkill` → `Cb`) and
exposure time**: a low dose gives a low rate and an *incomplete* plateau within the treatment window. That
drug-graded plateau is the observable a PK-PD model exists to reproduce.

**Units — AND A UNIT INCONSISTENCY THAT MUST BE STATED (see §5.4).** `dt` [day], `dkill` [nM], `H` must be
dimensionless ⇒ **`k_death` in the QSS path carries units of nM⁻¹·day⁻¹**, *not* the "1/day" its own docstring
claims (`:299`). In the **kinetic** path (EQ-16) the same numeric constant multiplies a quantity already in
day⁻¹, so there it is **dimensionless**. The same number 1.0 is used in both. This is defensible only because
1.0 is the identity, but it means **`k_death` is not transferable between the two engines** and its stated
units are wrong for one of them.

---

### EQ-11 — Costim occupancy → per-cell signalling → effector gain / exhaustion decay — `engine/wholebody_pd.py:333–347`

```
--- QSS PATH ONLY (:332–334). NOT executed in the canonical kinetic path — see §5.12 ---
p_eng_j   = Cb_j / (Cb_j + R_A_j + 1e-12)                        engaged fraction, drives induction   (:333)
R_cos_j(t)= CostimInduction.step(dt, p_eng)   → writes self.R_costim[Tidx]                            (:334)

--- BOTH PATHS (QSS :336–347;  kinetic :448–455) ---
A         = mean over T cells of the POSITIVE RESTING costim copies,  fixed once in __init__          (:200–201)
occ_j     = [ C_j / (C_j + KD_costim) ] · [ R_cos_j / max(A, 1e-9) ]                                  (:336, :448)
occ_j     ← _apply_cis_coincidence(occ_j, e_j)   e = Cb (QSS) | kin.B2 (kinetic)   (EQ-12b)           (:337, :449)
occ_j     ← clip(occ_j, 0, 5)                                                                         (:338, :450)

PerCellSignaling.step(occ, dt·24)                                signalling kinetics are per-HOUR     (:339, :451)

g_eff_j   = exp( kE_gain · eff_p_j )   for CD8 cells only,  kE_gain = 0.55                            (:343–345, :453–454)
g_eff_j   ← g_eff_j · exp( −0.30 · max(exh_p_j, 0) )             exhaustion attenuation               (:347, :455)
```

> **The normalizer `A` (`_costim_anchor`) is a CONSTANT, not a running mean.** It is computed **once, in
> `__init__`** (`:200–201`), from the **resting** costim copies of the T-cell population
> (`Tpos[Tpos>0].mean()`). It is *not* re-derived from the induced `R_cos(t)`. This distinction is
> load-bearing: had the denominator been the *running* mean of `R_cos`, induction would **cancel out of
> `occ` entirely** and EQ-11's whole induction mechanism would be inert. Any statement of this equation that
> writes the denominator as "mean(R_cos over T cells)" is **wrong**.

**Biological meaning.** The costimulatory arm (4-1BB / OX40 / ICOS / GITR / CD28 / CD2 / CD27) engages its
receptor at Langmuir occupancy, scaled by *this* T cell's receptor copies relative to the population. That
occupancy drives a per-cell signalling state (an ODE integrator in `signaling_dynamics`, a separate subsystem)
which reads out as an **effector program** (more killing) and an **exhaustion program** (less durable killing).
Only CD8/effector cells get the effector gain (`:345`).

**Mechanistic rationale.** Two decisions matter here.
1. **`p_eng` drives induction, so tumour-conditionality is emergent — IN THE QSS PATH ONLY** (`:322–334`).
   4-1BB/OX40/ICOS/GITR are *absent on resting T cells* and appear only after TCR engagement. Reading a static
   resting density would see them at ≈0 and systematically **under-rate exactly the arms the field cares
   about**, while over-rating the constitutive arms (CD28/CD2/CD27) — producing, in the code's own words,
   *"a confidently-wrong ranking ('CD2 beats 4-1BB')"*. Because induction is driven by **this T cell's own
   trimer**, a T cell that never engages never upregulates: conditionality is computed, not imposed.
   Constitutive arms carry fold = 1.0 → the induced density collapses to the resting density → byte-identical
   to the old static path.
   > 🔴 **BUT `costim_ind.step()` IS CALLED IN EXACTLY ONE PLACE — `:334`, inside `step()` (the QSS path).**
   > `_step_kinetic` never calls it (grep-verified this task: `costim_ind` appears only at `:190`, `:222`,
   > `:332`, `:334`). So in the **CANONICAL** engine `R_costim` is **frozen at the resting density for the
   > entire run** and the bug this block was written to fix is **still live in the engine that actually
   > runs.** See §5.12. Do not quote this row as an emergent property of the canonical model.
2. **The gain is exponential in the program level, not linear.** `exp(0.55·eff_p)` is bounded below by 0 and is
   multiplicative with the exhaustion term `exp(−0.30·exh_p)`, so the two programs compose additively in the
   exponent — an agonist that raises both effector *and* exhaustion nets out, which is the actual pharmacology
   of 4-1BB agonism.

**Units.** `occ` dimensionless (Langmuir occupancy × relative copy number, clipped to [0,5]); `eff_p`, `exh_p`,
`supp_p` are dimensionless program levels from `signaling_dynamics`; `kE_gain`, `0.30` dimensionless
sensitivities; `dt·24` converts days → hours because the signalling integrator runs per-hour (`:339`).

> **`kE_gain = 0.55` and the exhaustion coefficient `0.30` are FITTED** — the code calls `kE_gain` a *"locked
> calib"* (`:343`) and gives no source for either. See §3.4.

---

### EQ-12 — Treg suppression program (extra damping) — `engine/wholebody_pd.py:349`

```
supp_extra = mean over Treg T-cells of  max( supp_p_j , 0 )        (scalar, organ-level)     (:349)
```

**Biological meaning.** Costim agonism does not only boost effectors — it can boost **Tregs** (GITR and CD28
agonism are the canonical worries). `supp_extra` is the average suppression program level across this organ's
Tregs, and it multiplies the per-Treg suppression strength in EQ-8.

**Mechanistic rationale.** This is the counterscreen's tox arm: an agonist whose signalling drives the
`suppression` program in Tregs raises `supp_extra`, which raises the EQ-8 denominator, which *lowers* the kill
of every T cell in the neighbourhood — so a nominally "more potent" costim can come out with **worse** net
efficacy. That trade-off is computed, not asserted.

**Units.** dimensionless.

---

### EQ-12b — CIS/TRANS coincidence gating of costim occupancy — `_apply_cis_coincidence`, `engine/wholebody_pd.py:282–295` (body `:289–295`)

```
if p_cis ≤ 0:  return occ                                     (trans / no costim → unchanged)  (:290–291)

e       = max(cd3_engagement, 0)                              (Cb in QSS; B2 in kinetic)       (:292)
ref     = median( e[e>0] )                                    (self-scaling normalizer)        (:293)
f_cd3   = e / (e + max(ref, 1e-12))                           ∈ [0,1)                          (:294)
occ_eff = occ · [ (1 − p_cis) + p_cis · f_cd3 ]                                                (:295)
```

`p_cis` itself comes from the T-side co-engagement span geometry via
`multiarm_binding._cis_feasibility(span_coeng_T_nm)` (`:207–211`); default `span_coeng_T_nm=None` → `p_cis = 0`.

**Biological meaning.** In a trispecific (CD3 × TAA × costim), whether the costim arm can engage its receptor
*on the same T cell, at the same time* as the CD3 arm depends on the **height-matching of the two epitopes**.
Height-matched → the two arms co-engage in *cis* → the costim signal fires **only** on T cells that are already
CD3-engaged. Mismatched (≈60 Å) → the costim arm engages in *trans*, i.e. on any T cell, engaged or not.

**Mechanistic rationale.** This is the coincident-signal design in a single line: at `p_cis=1` the costim
occupancy is multiplied by the T cell's own normalised CD3 engagement, so costim is gated on target engagement
— which is the *entire safety premise* of a tumour-conditional costim. At `p_cis=0` the costim is
cell-autonomous (the systemic-agonism failure mode). Coincidence is **emergent from span geometry × real
per-cell CD3 binding**, not a phenomenological "tumour-conditional" flag. Default `p_cis=0` → returns `occ`
byte-identically, so the clinical validation molecules (no costim arm) are unaffected.

**Units.** `p_cis`, `f_cd3`, `occ` dimensionless. `ref` in the units of the engagement measure (nM).

---

### EQ-13 — Organ engagement sum (the cytokine driver) — `engine/wholebody_pd.py:364`

```
eng = Σ_j  Cb'_j / ( 1 + TREG_K · n_treg_j )         [nM, summed over T cells]      (:364)
```

**Biological meaning.** Total engaged trimer across all T cells in the organ, Treg-damped. This is the quantity
cytokine output is proportional to.

**Mechanistic rationale.** Cytokine release is engagement-driven, so the driver is the **sum** of per-cell
trimer (an extensive quantity — more engaged T cells, more cytokine), whereas killing is the *apportioned*
per-target trimer (EQ-9). They are the same underlying species read two different ways.

> **Inconsistency to note (§5.3): the cytokine Treg damping OMITS `supp_extra`.** EQ-8 damps kill by
> `1/(1 + TREG_K·n_treg·(1+supp_extra))`; EQ-13 damps cytokine by `1/(1 + TREG_K·n_treg)` only. So a costim
> that drives the Treg suppression program suppresses **killing** but not **cytokine**. This looks unintended.

**Units.** nM (summed across cells — the code calls these *"engagement-sum units"*, `:33`).

---

### EQ-14 — Cytokine hierarchy and per-T-cell cytokine output — `engine/wholebody_pd.py:22`, `:370–376`

```
CYTO_HIER = { IL6: 1.00,  IFN: 0.36,  TNF: 0.31,  IL2: 0.18 }                             (:22)

cs            = mean_j ( 0.45·IFNG_p_j + 0.32·TNF_p_j + 0.18·IL2_p_j )                     (:370)
cyto_sig_gain = max( 0.2, 1 + cs )                                                         (:371)
resp          = max( 1 − Dcyto, 0 )                                                        (:373)

cyto_rate[k]  = CYTO_HIER[k] · eng · resp · cyto_sig_gain      for k in {IL6,IFN,TNF,IL2}  (:374)
cyto[k]      ← cyto[k] + dt · cyto_rate[k]                     (cumulative, AUC-like)      (:375)
```

**Biological meaning.** Every engaged T cell emits cytokine in proportion to its bridged trimer. The **four
species are not modelled independently** — they are emitted in a **fixed stoichiometric ratio** set by
`CYTO_HIER`, with IL-6 as the numeraire (1.00) and IFN-γ, TNF-α, IL-2 at 0.36, 0.31, 0.18 of it. The ratio is
then modulated *uniformly* by two gains: the costim signalling drive (`cyto_sig_gain`) and acute
desensitisation (`resp`, EQ-15).

**Mechanistic rationale — and its honest limit.** The hierarchy encodes the clinical observation that in TCE
CRS, IL-6 is the dominant and most-measured species and the others track below it. Modelling it as a fixed
ratio means **the model has exactly one cytokine degree of freedom** — all four species are the same time
course scaled by a constant. That is a deliberate simplification (there is no independent IFN-γ production
mechanism anywhere in the module), and it means the model **cannot** reproduce a dissociation between IL-6 and
IFN-γ kinetics. The `cyto_sig_gain` weights (0.45/0.32/0.18) are the *costim programs'* contribution to the
cytokine gain and are a **different, unrelated set of weights** from `CYTO_HIER` — they are not normalised and
have no source.

**Units.** `CYTO_HIER` weights dimensionless (ratios to IL-6). `cyto_rate` in **raw engagement-sum units**
(nM-summed-over-cells), **NOT pg/mL**. `cyto` = raw units × days.

> **CRITICAL: `cyto_rate` IS NOT A CONCENTRATION.** The only function that ever converted it to pg/mL is
> `cytokine_to_pgml` (`:32–35`), and that function is **ORPHANED — it has zero callers anywhere in the live
> tree** (verified by repo-wide grep this task). The live IL-6 **concentration** comes from a completely
> different, mechanistic path (`myeloid_il6` → plasma ODE, `:377–380`, `:480–486`). `cyto_rate` survives as a
> *relative* readout (`sys_cyto_rate` in the driver output) and must never be quoted in pg/mL.

> **THREE DIFFERENT "ENGAGED FRACTIONS" LIVE IN THIS MODULE — do not conflate them (§5.13).**
> | name | formula | line | consumer |
> |---|---|---|---|
> | induction driver (QSS only) | `Cb/(Cb + R_A + 1e-12)` | `:333` | `CostimInduction.step` |
> | myeloid driver (QSS) | `clip(Cb / max(R_A, 1e-30), 0, 1)` | `:379` | `MyeloidIL6.step` |
> | myeloid driver (kinetic) | `clip(B2 / max(RC, 1e-30), 0, 1)` | `:485` | `MyeloidIL6.step` |
> The first two are **different functions of the same two quantities** (`x/(x+R)` vs `x/R`) and they are not
> equal (they agree only as `Cb/R_A → 0`). Both are called `p_eng` in the source.

---

### EQ-15 — Acute cytokine desensitisation (tachyphylaxis) — `engine/wholebody_pd.py:29–31`, `:373`, `:382–384`

```
eng_norm = eng / DESENS_ENG_REF                        DESENS_ENG_REF = 1.0e4          (:31, :382)
resp     = max( 1 − Dcyto, 0 )                                                          (:373)

dDcyto/dt = K_CYTO_DESENS · eng_norm · resp  −  K_CYTO_RECOV · Dcyto                    (:383)
            K_CYTO_DESENS = 30.0 /day        K_CYTO_RECOV = 0.003 /day                  (:29, :30)
Dcyto    ← clip( Dcyto, 0, 0.999 )                                                      (:384)
```

**Biological meaning.** Clinically, CRS is **worst at the first exposure** and attenuates on subsequent doses
even as drug exposure rises — that is why every approved TCE uses a step-up dosing schedule. `Dcyto` is a
lumped "cytokine responsiveness has been used up" state: it builds while T cells are engaged and recovers over
days. Production is gated by `resp = 1 − Dcyto`.

**Mechanistic rationale.** Without this, cytokine output simply **tracks bound drug** (`:26–28`), so the model
predicts the *worst* CRS at the *highest* dose — the opposite of the clinic. The `(1 − Dcyto)` factor inside the
build term (`:383`) makes the state saturate structurally at 1 rather than needing a cap (the cap at 0.999 is
a numerical guard). The asymmetry `K_CYTO_DESENS = 30 /day` vs `K_CYTO_RECOV = 0.003 /day` (a **10,000×**
separation) is what produces a fast spike-and-collapse (~2–3 days, `:27`) followed by very slow recovery.

**Units.** `K_CYTO_DESENS`, `K_CYTO_RECOV` in day⁻¹; `DESENS_ENG_REF` in raw engagement-sum units;
`Dcyto`, `resp` dimensionless ∈ [0,1).

> **All three constants are FITTED** — *"ported from qsp_costim_window_v2 cytokine calib"* (`:26`), with
> `DESENS_ENG_REF` explicitly *"~ mosun 1mg peak"* (`:31`), i.e. **normalised against the same mosunetuzumab
> anchor whose IL-6 value has no source**. See §3.5.

---

### EQ-16 — KINETIC path (the CANONICAL engine): T6 as a modifier wrapper — `engine/wholebody_pd.py:438–491`

When `pd_kinetics=True`, T6 does **not** compute a trimer. It computes the costim modifiers (`:446–456`), then
delegates the bond kinetics to `kinetic_synapse.KineticSynapse` and applies its modifiers to the returned
hazard:

> ⚠️ **The costim block here is NOT identical to the QSS block.** `:446–456` reproduces the occupancy →
> signalling → `g_eff` / `supp_extra` chain, but it **silently omits the activation-induction step**
> (`:332–334`): there is no `costim_ind.step()` call anywhere in `_step_kinetic`, so `self.R_costim` is read
> (`:447`) but never written. The two paths differ in exactly one thing, and it is the thing the costim
> screen depends on. **§5.12.**

```
dH_tgt   = KineticSynapse.step( C_T, dt, k_death, per_target_surv )    [per-target hazard incr]  (:461)
gscale   = mean( g_eff )                                               [SCALAR over T cells]     (:463)
treg_damp= 1 / ( 1 + TREG_K · mean(n_treg) · (1 + supp_extra) )        [SCALAR over T cells]     (:464)
dH_tgt  ← dH_tgt · gscale · treg_damp                                                            (:465)

H[tgtidx]        ← H[tgtidx] + dH_tgt                                                            (:467)
per_target_surv  = exp( −H[tgtidx] )      → fed BACK into the next step (serial killing)         (:468)

eng      = KineticSynapse.engaged_dwell_rate() · treg_damp             [engaged DWELL, not Cb]   (:471)
p_eng_j  = clip( B2_j / RC_j , 0, 1 )                                  → myeloid IL-6 driver     (:485)
```

Kinetic-engine construction (`:227–248`), including the unit conversions T6 owns:

```
kon [/nM/day]  = kon [/M/s] · 86400 / 1e9                                                (:231–232)
koff[/day]     = koff_lit [/s] · 86400   , or, if no literature koff:  koff = kon · KD    (:234–235)
dens_copies_T  = syn_TAA_mean / NM_PER_COPY          (recover copies from nM for the ceff) (:240)
```

**Biological meaning.** The equilibrium trimer (EQ-2) knows *how much* bridge exists but nothing about *how
long it lasts*. Killing is a **cycle**: engage → deliver a lethal hit (takes time) → detach → find the next
target. The kinetic engine models that cycle literally, so **synapse dwell time and serial killing are
emergent**, and `per_target_surv` feeding back (`:468`) means a T cell stops wasting hits on cells it has
already killed.

**Mechanistic rationale for making this canonical** (`pd_model_config.py:18–25`): dwell time and serial
throughput are *precisely the axes an affinity + format design sweep tunes*. A QSS equilibrium law is
structurally blind to them — it would rank a slow-koff high-affinity CD3 binder as strictly better, when in
fact such a T cell **cannot detach** and is throughput-limited (`kinetic_synapse.py:46–47`). A counterscreen
built on the QSS law would therefore give confidently wrong format guidance.

**The structural losses of the wrapper (flagged in §5.5, §5.12).** `gscale` (`:463`) and `treg_damp` (`:464`)
are **organ-mean scalars**, whereas the QSS path applies both **per T cell** (`:345`, `:351`). So in the
canonical engine:
- a T cell's own effector gain does **not** attach to the targets *it* kills — every target gets the organ-mean gain;
- a target in a Treg-dense pocket is damped by the **organ-mean** Treg count, not its local one;
- **activation-induced costim receptor density is never updated** (no `costim_ind.step()` in `_step_kinetic`) —
  every inducible arm runs at its **resting** density for the whole simulation (§5.12).

The spatial resolution of the Treg and costim effects, and the induction of the costim receptor, are therefore
**lost in the canonical path** and retained only in the non-canonical QSS path.

**Units.** `dH_tgt` dimensionless; `gscale`, `treg_damp` dimensionless; `k_death` dimensionless *here* (see EQ-10).

---

### EQ-17 — T-cell migration (persistent random walk + chemotaxis) — `move_immune`, `engine/wholebody_pd.py:389–436`

```
step_um  = speed_um_per_min · 60 · 24 · dt_days                           speed = 5 µm/min      (:398)
θ        ← θ + N(0, 0.5)                                                  heading persistence   (:404)
(dx,dy)  = ( cos θ, sin θ ) · step_um                                                            (:406)

if chemotax > 0 and targets exist:
    (vx,vy) = (centroid of targets) − (x_T, y_T);   v̂ = (vx,vy)/‖·‖                              (:411–412)
    dx = (1−χ)·dx + χ·step_um·v̂ₓ ;   dy = (1−χ)·dy + χ·step_um·v̂_y      χ = chemotax = 0.4      (:413–414)

x,y  ← clip( x+dx, y+dy,  tissue bounds )                                                        (:416–417)
rebuild → _build_neighborhoods()  (W, syn_nb, Wt_norm, syn_TAA_mean, n_treg all refreshed)       (:420)
```

**Biological meaning.** T cells are motile. Over a multi-day treatment window they migrate, so the synapse
graph is **not static** — new T:target contacts form, which is how a T cell that started in an antigen-poor
region eventually reaches targets. `chemotax = 0.4` biases 40% of each step up the antigen gradient.

**Mechanistic rationale.** Without migration, killing saturates at whatever contacts existed at t=0 and the
model would systematically under-predict late depletion. The persistent random walk (heading + angular
diffusion, `:404`) rather than an isotropic random walk reproduces the observed **directional persistence** of
T-cell tracks. On rebuild, the kinetic-synapse engine is reconstructed against the new neighbourhoods with the
**CD3-side bound state carried over** (`:434–435`) so migration does not silently reset the bond state.

**Units.** `speed_um_per_min` µm·min⁻¹; `step_um` µm; `θ` rad; `chemotax` dimensionless ∈ [0,1].

> The chemotactic bias points at the **global target centroid** (`:411`), not the *local* antigen gradient —
> the code's own comment says "fall back to global centroid" but there is **no local-gradient branch
> implemented**. See §5.6.

---

### EQ-18 — Calibration loading & engine selection — `engine/pd_model_config.py:49–75`

```
PD_ENGINE = "kinetic"                                                                    (:29)

_load_calib():
    read the FIRST of:  <engine>/handoff/kinetic_calib.json  |  <engine>/kinetic_calib.json    (:55–56)
    IL6_SCALE_KIN = d["IL6_SCALE_kin"];  K_HIT = d.get("k_hit", 12.0);  K_DEATH = d.get("k_death", 1.0)  (:59)
    KINETIC["k_hit_perday"] = K_HIT                                                       (:60)

    if NO file found → hard-coded fallback:  IL6_SCALE_KIN = 0.05473; K_HIT = 12.0; K_DEATH = 1.0   (:63)

_load_calib() runs at import time                                                         (:75)
```

**Mechanistic rationale.** One file decides which kill law every harness runs, so a run can never
silently mix the abstract window-scoring law, the QSS law and the kinetic law (`pd_model_config.py:7`:
*"they are NOT interchangeable, do NOT mix them in one run"*).

> **The fallback at `:63` is a silent-failure hazard.** If the calibration JSON is missing, the module does
> **not** raise — it substitutes literals and continues. `IL6_SCALE_KIN = 0.05473` is a fit to the
> contaminated 570 anchor (§3.5), so a missing-file run silently proceeds carrying it. (It is currently
> harmless *only* because `IL6_SCALE_KIN` turns out to have **no live consumer** — see §3.5.)

---

## 3. PARAMETERS OWNED

Every constant defined in `engine/wholebody_pd.py` and `engine/pd_model_config.py`. Provenance tags per the
task rules. **Where a code comment asserts a source I could not confirm, the row is tagged
[UNVERIFIED CITATION] and the comment's claim is quoted rather than endorsed.**

### 3.1 Geometry & spatial constants

| Symbol | Value | Units | Tag | Source | Mechanistic rationale |
|---|---|---|---|---|---|
| `R_SYN_UM` (`:36`) | 30.0 | µm | **[UNSOURCED — TBD]** | code comment: *"synapse reach (validated tumor value)"* — a prior-model value, no citation, no PMID | The radius within which a T cell can form a synapse. Sets how many targets one T cell can reach → directly sets kill density. **A pure spatial-model parameter with no measurement behind it.** |
| `R_TREG_UM` (`:71`) | 50.0 | µm | **[UNSOURCED — TBD]** | code comment: *"Treg suppression neighborhood (validated tumor value)"* — no citation | Radius over which a Treg suppresses. Larger than the synapse reach (suppression is partly paracrine, contact-independent) — a defensible ordering, but the value itself is not sourced. |
| `AVO` (`:20`) | 6.02214076e23 | mol⁻¹ | **[DEFINED — SI exact, not measured]** | SI 2019 redefinition: the mole is *defined* by this exact integer. **Not a measurement; nothing in this doc is [MEASURED].** | Avogadro. |
| `speed_um_per_min` (`:389`, default arg) | 5.0 | µm·min⁻¹ | **[UNVERIFIED CITATION]** | code comment (`:393`): *"T cells migrate ~2-10 um/min in tissue (Miller 2002/2004 two-photon)"* — **no PMID or DOI in code**; I did not verify the Miller papers in this task | Migration speed. 5 is the mid-point of the quoted 2–10 band. The **band** is plausible two-photon physiology; the **point value** is a mid-range pick. |
| `chemotax` (`:389`, default arg) | 0.4 | — | **[UNSOURCED — TBD]** | none in code | Fraction of each migration step biased toward antigen. |
| `R_TREG`/`R_SYN` KD-tree metric | Euclidean 2-D | — | [ASSUMED: the agent tables are 2-D spatial reconstructions] | `:257`, `:276` | Tissue is modelled as a 2-D slice. |

### 3.2 Receptor-capacity scaling — RCAP_CD3, RCAP_TAA, NM_PER_COPY

| Symbol | Value | Units | Tag | Source | Mechanistic rationale |
|---|---|---|---|---|---|
| `RCAP_CD3` (`:73`) | 2.0 | nM | **[FITTED: compartmental QSP receptor-capacity scale] — and DEAD in this module** | Code comment claims *"validated CD3 receptor-capacity scale (unified_pbpk_pd Params.pbpk.Rcap_CD3)"*. `unified_pbpk_pd` is **not in this tree**; the value is present at `engine/qsp_costim_window_v2.py:226` as `Rcap_CD3: float = 2.0  # CD3 receptor capacity scale (x tissue T-cell density)` — i.e. **a multiplier on a tissue T-cell density in a compartmental QSP model, with no citation there either.** **[UNSOURCED — TBD]** at origin. | **NOT USED.** Repo-wide grep (this task): `RCAP_CD3` appears in `wholebody_pd.py` **only at its own definition line 73** and in comments. It has **zero consumers**. It was superseded by `NM_PER_COPY` (EQ-1) and left behind. Its *value* survives only as the sanity-check at `:81` ("CD3 ~92k copies → ~2.1 nM, recovers the tumor's uniform Rcap_CD3=2.0" — RUN-verified this task: **2.148 nM**). |
| `RCAP_TAA` (`:74`) | 6.0 | nM | **[FITTED: same compartmental scale]** — and DEAD as a symbol, but **LOAD-BEARING as a numerator** | Same: `engine/qsp_costim_window_v2.py:227` — `Rcap_TAA: float = 6.0  # TAA receptor capacity scale in tumor (high antigen density)`. **No citation at origin. [UNSOURCED — TBD].** | The **symbol** has zero consumers, but its **value 6.0 is hard-coded into `NM_PER_COPY = 6.0/257000`** (`:83`). So the entire per-cell concentration basis of the model rests on this unsourced compartmental scale. |
| `NM_PER_COPY` (`:83`) | 6.0/257000 = **2.3346e-5** | nM·copy⁻¹ | **[DERIVED: from RCAP_TAA=6.0 (fitted) pinned at CEACAM5 = 257,000 copies/cell]** | `:76–78`. The 257,000 CEACAM5 copies/cell figure is **not sourced in this module** — **[UNSOURCED — TBD in this file]** (it may be sourced in the antigen-copy pipeline; I did not confirm it in this task). | Declares the synapse reaction volume: **71.1 pL** (RUN-verified this task). The one physical constant that turns copies into concentration. **See §5.1 — this is the model's most consequential unsourced number.** |
| implied `V_syn` | 7.113e-11 | L | [DERIVED: 1e9/(AVO·NM_PER_COPY), RUN-verified this task] | — | 71 pL is **enormous** for an immunological synapse (a whole T cell is ~0.2–0.5 pL). It is best read as an *effective* reaction volume absorbing all the unmodelled geometry, **not** as a physical synapse cleft volume. |

### 3.3 Binding & kill-law constants (the ones the committee will ask about)

| Symbol | Value | Units | Tag | Source | Mechanistic rationale |
|---|---|---|---|---|---|
| `k_death` (`:297` arg; set by `pd_model_config.K_DEATH`, `pd_model_config.py:47,59,63`) | **1.0** | nM⁻¹·day⁻¹ (QSS) / dimensionless (kinetic) — **see EQ-10** | **[FITTED: clinical anchoring — see the full honest account in §3.5]** | `params/kinetic_calib.json` → `"k_death": 1.0` with a `k_death_provenance` string (quoted in §3.5). | Converts bridged trimer into a death hazard. **This is the single free potency knob of the kill law.** |
| `k_hit_perday` (`:117`; `pd_model_config.py:41`) | **12.0** | day⁻¹ | **[UNVERIFIED CITATION — no PMID in code]** | Code says *"FIXED from serial-killing literature, not fitted"* (`pd_model_config.py:41`) and *"(~1 lethal hit per 2 h engaged) — FIXED from serial-killing literature, NOT fitted"* (`kinetic_synapse.py:48–49`). **No PMID, no DOI, no author, no year appears anywhere in the code for this value.** The only named source in the whole repo is inside `params/kinetic_calib.json`: *"inside Halle 2-16/day band"* — **an author name with no PMID, which I did not verify in this task.** | Ceiling on the serial-kill rate an engaged synapse can deliver. **Consumed by `kinetic_synapse`, not by T6 itself** — T6 only carries and forwards it. Actual throughput is `k_hit·koff/(k_hit+koff)`, so a slow-koff T cell is koff-limited, not k_hit-limited. |
| `KD_CD3_nM` (`:111`, default) | 40.0 | nM | **[UNSOURCED — TBD]** — no source anywhere in the live tree | Default in `OrganPD.__init__` (`:111`), `CoupledPerCellPD.attach_pd` (`coupled_percell_pd.py:17`) and the driver (`run_tce_pd_reval.py:123,152`). No comment, no citation. | The CD3 arm affinity. At 40 nM and CD3 = 2.15 nM, the CD3 arm is deep in the **linear regime** — a deliberate, physiologically-correct posture (`:81`), but the *value* is unsourced. |
| `KD_TAA_nM` (`:111`, default) | 1.45 | nM | **[UNSOURCED — TBD]; not exercised by the canonical driver** | The default is **always overridden on the driver path**, which passes a per-drug KD (`run_tce_pd_reval.py:152,162,167` → `cfg.get('KD_norm', cfg['KD'])`). It is **not** dead in general: it is still the fallback default of both `OrganPD.__init__` (`:111`) and `CoupledPerCellPD.attach_pd` (`coupled_percell_pd.py:17`), so **any other caller that omits `KD_TAA_nM` silently gets 1.45**. | Per-drug TAA affinity, supplied by the caller. |
| `KD_costim_nM` (`:112`, default) | 1.0 | nM | **[UNSOURCED — TBD]** | Default; driver passes `cfg.get('KD_costim_nM', 1.0)` (`run_tce_pd_reval.py:153`) — i.e. **the 1.0 default is live whenever a config omits it**. | Costim arm affinity (Langmuir occupancy, EQ-11). A generic placeholder. |
| `alpha` (`:85`, default) | 1.0 | — | **[ASSUMED: no cooperativity between the two arms]** | Never passed by any caller. | Schropp's cooperativity factor. `α=1` ⇒ binding of the first arm does not change the affinity of the second. Physically this is an assumption, not a null. |
| `kon_CD3_perM_s`, `kon_TAA_perM_s` (`:113`; `pd_model_config.py:33–34`) | 1e5 | M⁻¹·s⁻¹ | **[ASSUMED: "standard mAb assoc"]** — code comment `pd_model_config.py:33` | No PMID. | Association rates for the kinetic engine. Sets the koff via `koff = kon·KD` when no literature koff is supplied (`:234–235`). |
| `koff_CD3_pers`, `koff_TAA_pers` (`:114`; `pd_model_config.py:35–36`) | `None` | s⁻¹ | [CODE-INTERNAL: None → derive from kon·KD] | — | Literature-koff override hook. **Currently unused (None) → every koff in the canonical engine is DERIVED from the assumed kon=1e5 and the KD.** So the affinity is real (per-drug KD) but the *kinetics* (the thing the kinetic engine exists to resolve!) are set by an assumed kon. **See §5.7.** |
| `kint_bridge_perday` (`:114`; `pd_model_config.py:37`) | 0.9 | day⁻¹ | **[UNSOURCED — TBD]** | no comment source | Trimer internalisation rate (consumed by `kinetic_synapse`). |
| `span_bridge_nm`, `span_cis_nm` (`:115`; `pd_model_config.py:38–39`) | 12.5 | nm | **[ASSUMED / format-overridable]** — comment: *"AF3/format override per construct"* | — | Arm spans; feed the cis-feasibility geometry (EQ-12b) and the kinetic engine's cleft geometry. |
| `cis_avidity` (`:115`; `pd_model_config.py:40`) | 0.0 | — | [CODE-INTERNAL: 0 for plain CD3×TAA clinical validation] | — | Costim co-engagement avidity; off for the validation molecules. |
| `n_CD3`, `n_TAA`, `n_costim` (`:116`) | 1, 1, 0 | valency | [CODE-INTERNAL defaults] | — | Arm valencies. |
| `TREG_K` (`:72`) | **0.25** | per-Treg | **[UNSOURCED — TBD]** | code comment: *"per-Treg suppression constant (validated tumor value)"* — no citation, no PMID | Sets the strength of spatial Treg suppression (EQ-8). One Treg in reach → kill ×0.80; four → ×0.50. **This constant single-handedly sets how much a Treg-rich organ is protected, and it is unsourced.** |

### 3.4 Costim signalling coupling constants (all inside `wholebody_pd.py`)

| Symbol | Value | Units | Tag | Source | Mechanistic rationale |
|---|---|---|---|---|---|
| `kE_gain` (`:343`, `:453`) | 0.55 | — | **[FITTED: code says "locked calib" (`:343`); fit target NOT stated in code]** | none | Effector-program → kill-gain sensitivity, in the exponent: `g_eff = exp(0.55·eff_p)`. |
| exhaustion coefficient (`:347`, `:455`) | 0.30 | — | **[UNSOURCED — TBD]** | bare literal, no comment | Exhaustion attenuation `× exp(−0.30·exh_p)`. |
| cytokine-drive weights (`:370`, `:475`) | IFNG 0.45, TNF 0.32, IL2 0.18 | — | **[UNSOURCED — TBD]** | bare literals | Weights of the three costim cytokine programs into `cyto_sig_gain`. **Note these are NOT `CYTO_HIER` and do not match it.** |
| `cyto_sig_gain` floor (`:371`, `:475`) | 0.2 | — | [ASSUMED: numerical floor so a "cold" arm cannot zero cytokine] | — | `max(0.2, 1+cs)`. |
| occupancy clip (`:338`, `:450`) | [0, 5] | — | [ASSUMED: numerical guard on the receptor-copy amplification] | — | Bounds `occ` so a rare very-high-copy T cell cannot blow up the signalling integrator. |
| signalling time base (`:339`, `:451`) | `dt·24` | h | [CODE-INTERNAL unit bridge] | — | `signaling_dynamics` runs per-hour; T6 runs per-day. |
| `_costim_anchor` (`:200–201`) | mean of the **positive RESTING** costim copies over T cells | copies | [CODE-INTERNAL: population normaliser, **frozen at construction**] | — | Makes `occ` a *relative* receptor amplification, so the absolute copy scale of a costim arm cancels at t=0 — **this deliberately removes cross-arm absolute-abundance differences from the occupancy term**. It is a **constant**, computed once in `__init__` from the resting `Rc`, **not** a running mean of the induced `R_costim`; if it were the running mean, induction (EQ-11) would cancel out of `occ` and be inert. Abundance differences re-enter only through induction — **and only in the QSS path** (§5.12). |

### 3.5 Cytokine constants — **CONTAMINATION SECTION, READ IN FULL**

| Symbol | Value | Units | Tag | Source | Verdict |
|---|---|---|---|---|---|
| `CYTO_HIER` (`:22`) | IL6 1.0, IFN 0.36, TNF 0.31, IL2 0.18 | — (ratios) | **[UNSOURCED — TBD]** | code comment: *"project-standard CYTO_HIER, mosunetuzumab-anchored"*. **No PMID, no figure, no table cited.** The word "mosunetuzumab-anchored" points at the same drug whose in-code IL-6 anchor (570) the audit found to have no source. | **LIVE and consumed** (`:374`, `:477`; re-exported to `coupled_percell_pd.py:14`). The four ratios have **no verifiable provenance in the code**. |
| `CYTO_IL6_CLINICAL_ANCHOR_PGML` (`:23`) | **570.0** | pg/mL | **[FITTED to a FABRICATED anchor] — CONTAMINATED** | In-code comment claims *"mosunetuzumab peak IL-6 ~570 pg/mL (Hosseini 2020 Fig5A)"*. **The 2026-07-13 provenance audit established that 570 has NO SOURCE.** The only valid clinical IL-6 anchors are **mosunetuzumab = 152** and **teclistamab = 21**, both population means. The driver's own comment now says so explicitly: *"That 570 has since been shown to have NO SOURCE — it is a fabricated anchor"* (`run_tce_pd_reval.py:205–206`). | **ORPHANED-DEAD but still physically present.** Repo-wide grep (this task): the symbol has **zero consumers**. **DO NOT RE-WIRE. DO NOT CITE. It should be deleted from the module.** |
| `cytokine_to_pgml()` (`:32–35`) | function | — | **[DEAD CODE]** | — | **Zero callers anywhere in the live tree** (verified by grep this task). It is the function that used to convert `cyto_rate` → pg/mL using the fitted scale. The path it served was **deleted from the driver on 2026-07-13** (`run_tce_pd_reval.py:203–221`), which now **raises a hard error** rather than falling back to it. |
| `IL6_SCALE_kin` / `IL6_SCALE_KIN` (`pd_model_config.py:45,59,63`) | **0.05473008541177651** | pg/mL per raw unit | **[FITTED: BY CONSTRUCTION to the fabricated 570]** | `params/kinetic_calib.json` contains `"IL6_SCALE_kin": 0.05473008541177651`, `"il6_raw_peak": 10414.747130603793`, `"clin_il6": 570.0`. **I verified the arithmetic this task: 570 / 10414.747130603793 = 0.05473008541177651 — EXACT.** The "calibration" is literally `clinical_anchor / model_raw_peak`. | **The number is a one-point division by a number with no source.** It is loaded at import (`pd_model_config.py:75`) and assigned in the driver (`run_tce_pd_reval.py:41`) but — verified by grep this task — **never consumed** thereafter. **Dead, but loaded.** |
| `K_CYTO_DESENS` (`:29`) | 30.0 | day⁻¹ | **[FITTED: "ported from qsp_costim_window_v2 cytokine calib" (`:26`)]** | no independent source | LIVE (`:383`, `:488`). |
| `K_CYTO_RECOV` (`:30`) | 0.003 | day⁻¹ | **[FITTED: same]** | no independent source | LIVE. |
| `DESENS_ENG_REF` (`:31`) | 1.0e4 | raw eng-sum | **[FITTED: comment says "~ mosun 1mg peak" (`:31`)]** — i.e. normalised to a mosunetuzumab model output | no independent source | LIVE. It only sets *where on the dose axis* desensitisation kicks in; it is **not** an IL-6 concentration and does not carry the 570 contamination directly, but it is anchored to the same molecule's simulated peak. |
| `_MYELOID_COUNT_SCALE` / `myeloid_count_scale` (`:42–70`, `:144`) | per-organ; **1.0 if census absent** | — | **[DATA-GATED — census file]** | Loaded from `handoff/organ_myeloid_counts.json`; **present in this tree** (RUN-verified this task: 12 organs loaded, e.g. spleen = 290,206). | Scales sampled myeloid counts → physiological counts. Plasma IL-6 is **exactly linear** in it (`:41`), so a run at 1.0 is analytically re-scalable. The loader **announces on stderr** and warns loudly if the file is missing (`:66–68`) — a deliberate anti-silent-failure guard installed after a measured ~2.9e5× IL-6 under-prediction (`:44–46`). **This is the one place in the module where a missing input screams instead of defaulting quietly.** |

> **Also inherited, not owned:** IL-6 clearance (used by the plasma ODE downstream of the `il6_prod_pg_hr`
> T6 emits) is cited in the codebase to **PMID 31268236**. Per the 2026-07-13 audit, that paper (Chen 2019,
> *Clin Transl Sci*) is a **semi-mechanistic modelling** paper and reports **no measured** IL-6 clearance.
> **Tag it [FITTED], never [MEASURED]. Human IL-6 clearance appears to be UNMEASURED in the literature.**
> T6 does not own that constant, but T6's `il6_prod_pg_hr` is its only input, so any IL-6 concentration
> quoted from this pipeline inherits that gap.

### 3.6 The honest account of `k_death` — read this before defending the model

`k_death` is **set in exactly one place**: `params/kinetic_calib.json` → `"k_death": 1.0`, loaded by
`pd_model_config._load_calib()` (`:59`), forwarded by the driver (`run_tce_pd_reval.py:41,173`) into
`CoupledPerCellPD.run(k_death=…)` → `OrganPD.step(…, k_death)` (`:297`).

**What the code comments claim, and why you cannot trust them.** `wholebody_pd.py` contains **two stale
docstrings that contradict the live calibration file**:
- `:299` — *"k_death (1/day) = trimer->death rate constant, calibrated to the validated tumor 28.7%
  (calibrate_kdeath)"*
- `:361` — *"kkill is calibrated against the validated tumor 28.7% (see calibrate_kkill)"*

Neither `calibrate_kdeath` nor `calibrate_kkill` exists in this module or anywhere in the live tree
(grep, this task). More importantly, **the live calibration file explicitly repudiates that anchoring**:

> `params/kinetic_calib.json` → `k_death_provenance`: *"k_death=1.0 LOCKED. PRIMARY anchor
> (compartment-independent): engaged-CTL serial ceiling k_hit·koff/(k_hit+koff) = 11.6/day at KD_CD3 = 40nM;
> emergent engaged rate 2.7/day median inside Halle 2-16/day band, NO fitting. SECONDARY check:
> mosunetuzumab full PBPK-PD at k_death=1.0 gives IL6 609 (clin 570), organ B-cell depl 0.93, heme blast
> depl 0.53. … the prior depletion-fit 0.5 was an artifact of forcing heme clinical depletion through a
> solid-geometry+ratio engine. ONE shared value, all engagers."*

**My verdict, stated for a committee:**

1. The **serial-ceiling arithmetic checks out.** RUN-verified this task: `kon = 1e5 M⁻¹s⁻¹ = 8.64 nM⁻¹day⁻¹`;
   `koff_CD3 = kon·KD_CD3 = 345.6 day⁻¹`; `k_hit·koff/(k_hit+koff) = 12·345.6/(12+345.6) = **11.597 /day**`
   — the JSON's "11.6/day". So the *consistency argument* is real and reproducible.
2. **But that is a consistency check, not a measurement.** `k_death = 1.0` is not derived from any measured
   trimer→apoptosis rate constant. No such measurement is cited anywhere in the code. The value 1.0 is a
   **round number that was checked against clinical outputs and locked**.
3. **The clinical cross-check named in the JSON is against the contaminated 570 anchor** ("IL6 609 (clin
   570)"). That leg of the argument is **void**.
4. **There IS a documented depletion-fitting history.** `engine/calib_kdeath.py` — a **DEAD module, correctly
   excluded from the live import graph, and NOT documented here as part of the model** — exists and reads:
   *"k_death calibration vs epcoritamab depletion time-course (day7 0.30 / day14 0.90 / day28 0.94). ONE fresh
   process per k_death"*, computing an RMSE per `k_death` value (`calib_kdeath.py:1–2`, RMSE at `:28`, written
   out `:29–31`). It is cited **only** as historical provenance evidence for `k_death`, never as live
   mechanism. The JSON's own provenance string refers to *"the prior depletion-fit 0.5"* — i.e. **a depletion
   fit was in fact performed, produced 0.5, and was then rejected as an artifact and replaced by 1.0.**

**Therefore: `k_death` = [FITTED / CALIBRATED-AND-LOCKED].** It is a calibrated constant that was selected by
comparing model output to clinical depletion and IL-6, and it should **never** be presented as a measured or
literature parameter. The most defensible thing that can be said for it is: *a single shared value of 1.0 is
used for every engager, it is consistent with an independent serial-killing rate ceiling, and no per-drug or
per-compartment tuning is applied.* That is a genuinely strong position — **but it is a position about
parsimony, not about provenance.**

**And `k_hit` = 12/day is [UNVERIFIED CITATION].** The code asserts three times that it is "FIXED from
serial-killing literature, NOT fitted" and **never once names the literature**. The only trace of a source is
the surname "Halle" inside a JSON string, with no PMID, which I did not verify in this task. A reviewer will
ask for the citation, and the code cannot currently produce it.

---

## 4. WHAT IS EMERGENT vs IMPOSED

### Genuinely EMERGENT (computed from mechanism — these are the claims the model can defend)

| Property | Emerges from | Evidence |
|---|---|---|
| **Prozone / hook effect** | Receptor conservation in the Schropp QE solve (EQ-2): at high drug both arms saturate as separate binary complexes, starving the bridge. | **RUN-VERIFIED this task**: Cb peaks at **C = 7.58 nM (Cb = 0.1994 nM)** and falls ~1/C to 1.28e-3 nM at 1e4 nM. Nothing in the code contains a hook, a bell curve, or a descending limb — it falls out of the quadratic. |
| **Abundance-graded depletion (the counterscreen signal)** | Absolute copies → nM (EQ-1) with **no per-organ normalisation**. | RUN-verified: BCMA 11k copies → 0.257 nM (≲ KD_TAA) vs CD20 95k → 2.22 nM (≫ KD_TAA). Low-copy targets form a genuinely weaker bridge and deplete exposure-dependently instead of saturating identically. |
| **Spatial heterogeneity of kill** | `syn_TAA_mean` is a *neighbourhood* mean (EQ-6); kill is apportioned through the *actual* KD-tree contact graph (EQ-5/EQ-9). | A T cell in an antigen-poor pocket forms a weak trimer even at saturating drug. Organ-to-organ and intra-organ kill differences require no organ-specific parameter. |
| **Treg protection of Treg-rich tissue** | Local Treg census (EQ-7) → hyperbolic damping (EQ-8). | No global immunosuppression scalar. (But `TREG_K` itself is imposed and unsourced.) |
| **Tumour-conditionality of inducible costim arms** — 🔴 **QSS PATH ONLY; NOT ACTIVE IN THE CANONICAL ENGINE** | `p_eng` (this T cell's own trimer) drives receptor induction (EQ-11, `:333–334`) — but `costim_ind.step()` is called **only** inside `step()` (QSS), never inside `_step_kinetic`. | In the QSS path, a T cell that never engages never upregulates 4-1BB/OX40/ICOS/GITR. **In the canonical kinetic path the induction is never executed and every inducible arm runs at its resting density.** Grep-verified this task. **This is a code defect, not an emergent property — see §5.12.** |
| **Cis/trans coincidence** | `p_cis` from arm-span geometry × this T cell's own CD3 engagement (EQ-12b). | Coincident signalling is a geometric consequence, not a flag. |
| **Serial killing & dwell-time effects (canonical engine)** | `kinetic_synapse` engage/hit/detach cycle with `per_target_surv` fed back (EQ-16, `:468`). | A slow-koff T cell that cannot detach is throughput-limited — the model can *lose* efficacy by raising CD3 affinity. |
| **Cytokine self-limitation (spike then collapse)** | The `(1−Dcyto)` structural saturation in EQ-15. | Output does not track bound drug. |
| **Saturating IL-6 with per-drug differences** | Delegated to `myeloid_il6`: a **finite, spatially-distributed** myeloid pool contacted by *engaged* T cells. T6's contribution is `p_eng` (`:379`, `:485`). | No Emax, no EC50, no fitted IL-6 scale in the live path. |

### IMPOSED (handed in as a constant — this is where emergence stops)

| Imposed | Where | Why it matters |
|---|---|---|
| **`k_death` = 1.0** | the trimer→death potency | The **entire absolute kill scale**. Everything above is *shape*; this is the *magnitude*. It is calibrated (§3.6). |
| **`k_hit` = 12/day** | serial-kill ceiling | Sets the maximum throughput of an engaged synapse. Uncited (§3.6). |
| **`NM_PER_COPY` (71 pL synapse)** | copies → nM | Sets **where every receptor sits relative to its KD**. Move it and every trimer, every prozone position and every kill fraction moves with it. Derived from an unsourced compartmental scale (§3.2, §5.1). |
| **`TREG_K` = 0.25, `R_TREG` = 50 µm, `R_SYN` = 30 µm** | suppression & synapse geometry | The whole spatial structure. All unsourced. |
| **`CYTO_HIER` ratios** | cytokine stoichiometry | The model has **one** cytokine degree of freedom; the other three species are IL-6 × a constant. |
| **`kE_gain`=0.55, exhaustion 0.30, cytokine weights 0.45/0.32/0.18** | costim → PD coupling | These convert the signalling subsystem's abstract program levels into kill and cytokine. **They are the entire quantitative link between the costim screen and its readout, and they are fitted/unsourced.** |
| **`α` = 1 (no cooperativity)** | EQ-2 | An assumption presented as a default. |
| **`kon` = 1e5 M⁻¹s⁻¹ (both arms)** | EQ-16 | With `koff = kon·KD`, *all* the kinetics of the canonical kinetic engine derive from this one assumed number. |

**The honest one-line summary:** *the SHAPE of the dose-response (prozone), the SPATIAL structure of killing,
and the ABUNDANCE-ranking of targets are emergent and defensible. The MAGNITUDE of killing, the strength of
Treg suppression, the reaction volume, and the costim→PD coupling are imposed constants, most of them
uncited.*

---

## 5. KNOWN LIMITATIONS / OPEN QUESTIONS

### 5.1 The 71 pL "synapse" and the circular receptor-capacity pin — **the biggest attack surface**

`NM_PER_COPY = 6.0/257000` is presented (`:75–83`) as *"ONE physical constant: the synapse reaction volume."*
It is not measured. It is `Rcap_TAA` — a **receptor-capacity multiplier from a compartmental QSP model**
(`qsp_costim_window_v2.py:227`, itself uncited) — divided by a CEACAM5 copy number. Inverting it gives a
**71 pL** reaction volume, which is **two orders of magnitude larger than an entire T cell**. A reviewer will
immediately ask what physical object has a volume of 71 pL, and the honest answer is: *none — it is an
effective volume that absorbs the unmodelled synapse geometry and was pinned so that the previous model's
validated kill scale was reproduced.*

**Why it is still the right choice:** it is a **single, molecule-independent, organ-independent** constant, so
it cannot be used to tune one drug against another — which is exactly what a counterscreen requires. Its
consequence (BCMA sits below its KD, CD20 above) is *falsifiable* and *correct in direction*.

**What would fix it:** derive the synapse reaction volume geometrically (contact-area × cleft-height) as the
kinetic engine already does for its `ageff_nM` ceff, and check whether the two agree. If they disagree by 100×,
the model is carrying two incompatible volume conventions.

### 5.2 Three fitted constants anchored to a value that does not exist

`CYTO_IL6_CLINICAL_ANCHOR_PGML = 570.0` (`:23`), `IL6_SCALE_kin = 0.05473` (= 570/10414.75, verified exactly),
and `DESENS_ENG_REF = 1e4` ("~mosun 1mg peak") are all tied to a mosunetuzumab IL-6 peak of **570 pg/mL that
has no source**. The valid clinical anchors are mosunetuzumab **152** and teclistamab **21** (population
means). The first two are now **dead code** — the fitted IL-6 path was deleted from the driver and replaced by
a hard error (`run_tce_pd_reval.py:214–221`) — but **they are still in the module**, and `pd_model_config.py`
still *loads* `IL6_SCALE_KIN` at import. **Recommended action: delete `CYTO_IL6_CLINICAL_ANCHOR_PGML`,
`cytokine_to_pgml()` and `IL6_SCALE_KIN` outright.** A dead fitted constant in a live file is precisely how a
page number became a concentration.

`DESENS_ENG_REF` is **not** dead — it is live at `:382`/`:487` — and it is normalised to a mosunetuzumab model
peak. It does not carry a fabricated *concentration*, but its provenance chain ends in the same place.

### 5.3 The Treg suppression program damps killing but NOT cytokine

EQ-8 (`:351`) uses `1/(1 + TREG_K·n_treg·(1 + supp_extra))`. EQ-13 (`:364`) uses
`1/(1 + TREG_K·n_treg)` — **`supp_extra` is missing**. So a costim arm that drives the Treg suppression
program reduces killing but leaves cytokine untouched, which would make that arm look **artificially bad on
efficacy and artificially unchanged on tox**. This is almost certainly an oversight, and it biases the exact
comparison the counterscreen exists to make. **This should be verified against intent before any ranking is
published.**

### 5.4 `k_death` has different dimensions in the two engines

QSS: `H += dt·k_death·dkill` with `dkill` in **nM** ⇒ `k_death` is **nM⁻¹·day⁻¹**.
Kinetic: `dH = dt·k_death·serial_rate` with `serial_rate` in **targets·day⁻¹** ⇒ `k_death` is **dimensionless**.
The docstring at `:299` says "1/day", which is correct for **neither**. The same numeric value 1.0 is used in
both, so nothing is *numerically* wrong today — but the constant is **not transferable between engines**, and a
future re-calibration of one silently mis-scales the other.

### 5.5 The canonical (kinetic) engine LOSES the spatial resolution of the Treg and costim effects

In the QSS path, `g_eff` and the Treg damping are **per T cell** (`:345`, `:351`). In the **canonical** kinetic
path they are collapsed to **organ-mean scalars** before being applied (`gscale = mean(g_eff)`, `:463`;
`treg_damp` uses `mean(n_treg)`, `:464`). So in the engine that actually runs:
- a highly-activated CD8 next to a target confers **no more** kill on *that* target than a quiescent T cell across the organ;
- a target sitting inside a Treg cluster is protected by the **organ-average** Treg count, not its own.

The module's own comment concedes this (`:459`: *"applied to the per-T serial output via a global suppression
factor"*). **The spatial Treg biology — one of the model's headline features — is effectively mean-field in the
canonical engine.** Fixing it requires pushing `g_eff` and the per-T `treg_damp` **into** `KineticSynapse.step`
before the apportionment, rather than scaling the apportioned hazard afterwards.

### 5.6 Migration chemotaxis follows the GLOBAL centroid, not a local gradient

`move_immune` (`:411`) computes the bias vector toward `mean(x[tgtidx]), mean(y[tgtidx])` — the **global**
target centroid of the organ. The comment says *"vector to mean of currently-reachable targets (from existing
syn_nb); fall back to global centroid"*, but **only the fallback is implemented**. In a tissue with several
antigen-dense foci, every T cell walks toward the centre of mass of all of them — potentially an antigen-*empty*
point between foci. This will systematically distort late-time depletion in spatially heterogeneous organs.

### 5.7 The kinetic engine's kinetics are derived from an assumed kon

The canonical engine exists to resolve **dwell time** — and `koff = kon·KD` with `kon` **assumed** at
1e5 M⁻¹·s⁻¹ for both arms (`pd_model_config.py:33–36`; `wholebody_pd.py:234–235`; the literature-koff override
hook is `None` for every run). So the affinities are per-drug and real, but the **kinetics — the entire reason
the kinetic law was made canonical — are set by one uncited number applied identically to every construct.**
Two constructs with the same KD but different koff are, at present, **indistinguishable to the canonical
engine**. This is the single most important experimental input the model is missing.

### 5.8 Tregs are effectors

`is_treg ⊂ is_T` (`:156–159`), so Tregs enter `Tidx`, form trimers (EQ-2), and contribute to `kill_T`
(EQ-8) at baseline gain (they are excluded from the CD8 effector boost, `:345`, but not from killing).
A Treg-rich organ therefore gets **both** more suppression **and** more killers. Whether that is intended
(TCEs do redirect Tregs, and Treg-mediated killing is reported) or an artefact of the label mask needs an
explicit decision.

### 5.9 One cytokine degree of freedom

All four species in `CYTO_HIER` are the same curve × a constant (EQ-14). The model **cannot** represent an
IFN-γ/IL-6 kinetic dissociation, cannot represent a species-specific clearance, and cannot be validated
against IFN-γ or TNF time-courses independently. Only IL-6 has an independent, mechanistic path (via
`myeloid_il6`). **The other three species should be reported as *relative indices*, never as concentrations.**

### 5.10 Free drug is not depleted by binding within the PD step

EQ-2 conserves receptors but treats `C` as an inexhaustible free-ligand reservoir. At high target burden and
low dose (exactly the antigen-sink regime that matters clinically), the trimer will be **over-estimated**
unless the upstream transport/sink layer has already debited `C`. The seam between "who debits the drug" (PD or
transport) should be stated explicitly and tested with a conservation check.

### 5.11 Label-substring cell typing

`is_T`/`is_CD8`/`is_treg` are substring matches on free-text labels (`:156–159`). A relabelling upstream
(e.g. a cell type whose name happens to contain `cd4`) silently changes who kills. This class of bug —
substring matching on identity — has bitten this project before and should be replaced with exact
set-membership against a controlled vocabulary.

### 5.12 🔴 THE COSTIM RECEPTOR IS NEVER INDUCED IN THE CANONICAL ENGINE — a live code defect

**This is the most consequential finding of the 2026-07-13 adversarial re-verification, and it is a defect in
the CODE, not merely in the documentation.**

`CostimInduction` is constructed for every organ (`:220–224`) and `self.R_costim` is explicitly made a private
copy *because* "activation-induction WRITES into self.R_costim each step" (`:195–197`). But
`self.costim_ind.step(...)` is called in **exactly one place**: line **334**, inside `step()` — the **QSS**
path. Grep-verified this task: `costim_ind` occurs only at `:190`, `:222`, `:332`, `:334`.

`_step_kinetic` — **the canonical engine** (`pd_model_config.PD_ENGINE = "kinetic"`) — reads
`Rc_T = self.R_costim[Tidx]` at `:447` and **never writes it**. Therefore, in every canonical run:

- `R_costim` is **frozen at the resting scRNA-seq density** for the entire simulation;
- `occ` (`:448`) is a pure Langmuir occupancy × a **static** relative copy number;
- **4-1BB / OX40 / ICOS / GITR are read at their resting (≈0) densities** — which is precisely the failure the
  code comment at `:322–331` says the induction machinery was built to prevent: *"systematically UNDER-RATES
  exactly the arms that matter, while OVER-RATING the constitutive ones (CD28/CD2/CD27) → a
  confidently-wrong ranking ('CD2 beats 4-1BB')."*

**The bug the comment describes is still live, in the engine that actually runs.** The comment sits in the
QSS path, where the fix was applied; the fix was never ported to `_step_kinetic`.

**Consequences.** Any costim-arm ranking produced by the canonical engine is, at present, the **resting-density
ranking** — the exact artefact the induction subsystem exists to remove. Constitutive arms (CD28/CD2/CD27,
fold = 1.0) are unaffected; **inducible arms are systematically under-rated**. The clinical validation
molecules (`costim_arm=None`) are **not** affected — `self.sig is None` for them, so the whole block is
skipped.

**Fix (one line, in `_step_kinetic`, before `:447`):** replicate `:332–334`, driving induction with the
kinetic path's own engagement measure —
```python
if self.costim_ind is not None:
    p_eng_ind = np.clip(self.kin.B2/np.maximum(self.kin.RC, 1e-30), 0.0, 1.0)
    self.R_costim[Tidx] = self.costim_ind.step(dt, p_eng_ind)
```
**Until that lands, §4's "tumour-conditionality is emergent" claim is TRUE ONLY OF THE NON-CANONICAL QSS PATH
and MUST NOT be presented as a property of the model that produced any kinetic-engine result.**

### 5.13 Three different quantities are all called `p_eng`

`:333` uses `Cb/(Cb + R_A + 1e-12)` (a saturating fraction) to drive **induction**; `:379` uses
`clip(Cb/R_A, 0, 1)` (a *ratio*, clipped) to drive the **myeloid IL-6** emitter in the same QSS step; `:485`
uses `clip(B2/RC, 0, 1)` in the kinetic step. The first two are **different functions of the same two
numbers** and agree only in the limit `Cb ≪ R_A`. Nothing forces them to stay consistent, and a reader (or a
future edit) will assume "the engaged fraction" is one quantity. They should be unified, or renamed so the
three are visibly distinct.

### 5.14 What a reviewer will attack first, in order

1. *"Your canonical engine never induces the costim receptor."* → §5.12. **This one is a bug, not a
   limitation. Fix it before any costim ranking is shown to anyone.**
2. *"What is `k_death` and what did you fit it to?"* → §3.6. Answer honestly: **calibrated and locked**, one
   value for all engagers, consistent with an independent serial-killing ceiling; **not** measured.
3. *"Where does 12 hits/day come from?"* → §3.6. **The code cannot currently name the paper.** Get the PMID or
   re-tag the parameter.
4. *"A 71 pL synapse?"* → §5.1.
5. *"Your IL-6 was fitted to 570 pg/mL — where is that from?"* → §5.2. Answer: **it isn't from anywhere; that
   path is deleted; the live IL-6 is mechanistic myeloid.** Then be ready for: *"then why is 570 still in
   your source file?"*
6. *"Your model's headline is spatial Treg suppression — but your canonical engine averages it."* → §5.5.

---

## 6. QUICK REFERENCE — live call chain into T6

```
run_tce_pd_reval.py
  ├─ pd_model_config.PD_ENGINE = "kinetic"           (engine/pd_model_config.py:29)
  ├─ pd_model_config._load_calib()  → K_DEATH=1.0, K_HIT=12.0   (params/kinetic_calib.json)
  ├─ CoupledPerCellPD.attach_pd(KD_CD3_nM=40, KD_TAA_nM=<per-drug>, costim_arm=…, kin_params=KINETIC)
  │     └─ OrganPD(...)                               (engine/wholebody_pd.py:110)
  │           ├─ _build_neighborhoods()  → W, Wt_norm, syn_TAA_mean, n_treg    (:250–279)
  │           ├─ MyeloidIL6(...)                                                (:128–129)
  │           ├─ PerCellSignaling(...)  + CostimInduction(...)                  (:213, :222)
  │           └─ KineticSynapse(...)     [pd_kinetics=True → CANONICAL]         (:241–248)
  └─ CoupledPerCellPD.run(k_death=K_DEATH)
        └─ per step:  OrganPD.step(C_percell, dt, k_death)                      (:297)
              ├─ QSS      → EQ-2 … EQ-15                                         (:304–387)
              └─ KINETIC  → _step_kinetic → KineticSynapse.step + T6 modifiers   (:438–491)
```

**DEAD / NOT in the execution path (do not document, do not cite):** `cytokine_pbpk.py`, `il6_pbpk.py`,
`unified_binding.py`, `multiarm_kinetic.py`, `biexact_solver.py`, `rna_to_receptor.py`,
`convert_copies_ALL.py`, `calib_kdeath.py`. (`calib_kdeath.py` is referenced in §3.6 **only** as historical
evidence of a depletion-fitting exercise — its code is not part of the model.)

> In the **kinetic** branch, `CostimInduction` is **constructed (`:222`) but never stepped** — the arrow from
> `CostimInduction` into the per-step loop exists only on the QSS branch (`:334`). **§5.12.**

**DEAD INSIDE A LIVE FILE (delete on sight):** `CYTO_IL6_CLINICAL_ANCHOR_PGML` (`wholebody_pd.py:23`),
`cytokine_to_pgml()` (`:32–35`), `RCAP_CD3`/`RCAP_TAA` (`:73–74`, zero consumers),
`IL6_SCALE_KIN` (`pd_model_config.py:45,59,63`, loaded but never consumed).
All four re-confirmed zero-consumer by repo-wide grep during the 2026-07-13 adversarial re-verification.

---

## 7. CORRECTIONS APPLIED BY THE 2026-07-13 ADVERSARIAL RE-VERIFICATION

Every equation, line citation, parameter value and source claim in the first draft was re-checked against
`engine/wholebody_pd.py` and `engine/pd_model_config.py`. Findings:

**Citations: CLEAN.** Every source named in this doc was traced to a real string in the repo — Schropp 2019
(`wholebody_pd.py:86`), Hosseini 2020 Fig5A (`:23`, quoted as a *disputed* in-code claim, not endorsed),
Miller 2002/2004 (`:393`), "Halle 2-16/day band" (`params/kinetic_calib.json:20`), PMID 31268236
(`myeloid_il6.py:64`, `coupled_percell_pd.py:277`), mosunetuzumab 152 / teclistamab 21 + PMID 38831634
(`run_tce_pd_reval.py:74–77`). **No fabricated PMID, author or figure reference was found.** All arithmetic
re-ran exactly: `NM_PER_COPY` = 2.3346303e-5, `V_syn` = 7.1126e-11 L, serial ceiling = 11.597/day,
570/10414.747130603793 = 0.05473008541177651 (exact), the full prozone table, the myeloid census (12 organs,
spleen = 290,206).

**Substantive errors found and FIXED in this file:**

| # | Was | Is |
|---|---|---|
| 1 | EQ-16 claimed the kinetic path runs the costim block "identical" to QSS; §4 listed inducible-arm tumour-conditionality as EMERGENT. | **FALSE.** `costim_ind.step()` exists only at `:334` (QSS). The canonical engine never induces. Doc now carries a top-of-file danger callout, a corrected EQ-11/EQ-16, a demoted §4 row, and **§5.12** with the fix. |
| 2 | EQ-11 wrote the occupancy normaliser as `mean(R_cos over T cells)`. | **WRONG — that would make induction cancel out.** It is `_costim_anchor`, a **constant** fixed in `__init__` from the **resting** copies (`:200–201`). Corrected in EQ-11 and §3.4. |
| 3 | `AVO` tagged **[MEASURED: SI definition]** — the doc's only `[MEASURED]` tag. | Avogadro is an **exact SI definition**, not a measurement. Retagged **[DEFINED — SI exact]**. **No parameter in T6 is [MEASURED].** |
| 4 | EQ-2 prozone table labelled `C = 10 nM` "(PEAK)"; §4 said "peaks at C ≈ 10 nM". | It is the **grid max**. A 400-point scan puts the true peak at **C = 7.58 nM, Cb = 0.1994 nM**. Corrected in both places. |
| 5 | `KD_TAA_nM = 1.45` called "never used … dead". | Overstated. Dead **on the driver path only**; it remains the live fallback default of `OrganPD.__init__` and `attach_pd` for any other caller. |
| 6 | `summary()`/`summary_full()` both described as "O(1)". | `summary_full()` recomputes `exp(−H)` over every cell (`:499`). Corrected. |
| 7 | Line refs: Schropp `:87–88`; `_apply_cis_coincidence` `:289–295`; `speed_um_per_min`/`chemotax` `:390`; `calib_kdeath` RMSE `:29–32`; `pd_model_config` docstring `:3–5`. | Corrected to `:86–87`; `:282–295`; `:389`; `:28`; `:3–4`. |
| 8 | Three distinct quantities all called `p_eng` (`:333`, `:379`, `:485`) were documented as if one. | New **§5.13** tabulates all three and shows `Cb/(Cb+R_A)` ≠ `clip(Cb/R_A,0,1)`. |

**Not changed (checked, correct as written):** the Schropp algebra and every `a`/`b`/`d`/`R_Bf`/`R_Af`/`Cb`
line ref; EQ-3…EQ-10, EQ-12, EQ-13, EQ-15, EQ-17, EQ-18; the §5.3 `supp_extra` asymmetry (`:351` vs `:364`);
the §5.4 `k_death` dimensional inconsistency; the §5.5 mean-field collapse; the §5.6 global-centroid
chemotaxis; the zero-consumer status of `CYTO_IL6_CLINICAL_ANCHOR_PGML`, `cytokine_to_pgml`, `RCAP_CD3`,
`RCAP_TAA`, `IL6_SCALE_KIN`; the non-existence of `calibrate_kdeath`/`calibrate_kkill`; the `k_death`
provenance quotation (verbatim against `params/kinetic_calib.json`). No DEAD module is documented as live —
`calib_kdeath.py` is cited **only** as historical provenance evidence for `k_death` and is labelled DEAD at
every mention.
