---
title: "T5 — Kinetic multivalent-Rhoden immune synapse (per cell)"
subsystem: T5
model: costim_engager_counterscreen
source_file: engine/kinetic_synapse.py
source_md5: dc8872ae8d9e7a84ddf21a5b39a0cb60
source_lines: 242
date: 2026-07-13
generated_by: workflow-subagent T5
live_import_graph: "run_tce_pd_reval -> coupled_percell_pd -> wholebody_pd -> kinetic_synapse"
engine_status: "CANONICAL (pd_model_config.PD_ENGINE = 'kinetic', pd_model_config.py:29)"
---

up:: [[00_INVENTORY_AND_MAP]]
tags:: #atlas/generated #model-doc
dates:: 2026-07-13

> [!warning]+ PROVENANCE DISCIPLINE FOR THIS DOC
> Every number below was **read from `engine/kinetic_synapse.py` or from a file explicitly cited by
> file:line**, or was **computed by executing the live module** (all such results are marked **RUN-VERIFY**,
> with the exact configuration stated). Nothing is quoted from memory. Where the code's own comment asserts a
> literature provenance but names no PMID/DOI, the row is tagged **[UNVERIFIED CITATION]** and the comment's
> claim is reported verbatim rather than endorsed — **this doc does not guess at which paper a bare surname
> means.** Two in-code assertions were **tested and found NOT to hold** (§4.3 prozone; §5.2 "window optimum")
> — both are documented as measured, not asserted.
>
> **ADVERSARIAL RE-VERIFICATION, 2026-07-13** (source md5 `dc8872ae…` unchanged). Every equation, line
> citation, and parameter value was re-checked against the live code, and every RUN-VERIFY table was
> re-executed. **Nine defects were found and fixed in place.** The four that changed a stated number or claim:
> 1. **§5.9 / EQ-6 discriminant** — the previous "1 of 18,000 points, worst −8.97×10⁻⁶, essentially never
>    complex" was a **grid artifact**. The negative-discriminant region is a sharp *resonance* the old grid
>    stepped over. Re-measured: **2.93% of a resonance-resolved grid, worst −1.21×10⁻³**, and **cevostamab
>    (a molecule in the validation set) sits inside the gating condition.** Materially understated before.
> 2. **§4.3 / §4.2 tables** — low-`C` and low-`kint` rows were reported from an **unconverged 4.5-day probe**
>    but tabulated as steady state. Re-run to convergence; corrected. *Conclusions unchanged.*
> 3. **`k_death` / IL-6 570** — the doc previously stated the 570 pg/mL anchor "has NO SOURCE". **The code
>    does cite one** (`wholebody_pd.py:23`: *"Hosseini 2020 Fig5A"*). The correct statement is that the
>    in-repo audit found that attribution **does not survive checking** — which is a different, and citable,
>    claim. Now cited to `docs/PROVENANCE_AND_VALIDATION.md`.
> 4. **Two uncited names removed** — a co-attribution ("Vauquelin") and an uncited B-cell radius, neither of
>    which appears in the code. The code's own bare string *"Rhoden-2016"* is retained, tagged UNVERIFIED.

<!-- ATAGLANCE:START -->
> [!abstract]+ At a glance
> **T5 — the literal per-T-cell bond ODE that replaces the QSS/Schropp equilibrium trimer.**
> - **File:** `engine/kinetic_synapse.py` (242 lines; 3 module functions — `ageff_nM`, `cleft_feasibility`,
>   `_expm2x2_apply` — + 1 class `KineticSynapse` with 5 methods).
> - **States:** `B1` (armed: drug·CD3), `B2` (bridged trimer: CD3·drug·TAA), both nM per T cell. Per-target
>   `surv_j ∈ [0,1]` is owned by the caller.
> - **Integrator:** EXACT 2×2 matrix exponential (`_expm2x2_apply`, :74–102), coefficients frozen per step.
> - **Owned constants:** 8 (`NM_PER_COPY`, `R_CELL_UM`, `SPAN_BRIDGE_DEFAULT_NM`, `SPAN_CIS_DEFAULT_NM`,
>   `CLEFT_MIN_NM`, `CLEFT_MAX_NM`, `K_HIT_DEFAULT`, `AVO`) + 2 unnamed shape constants in the cleft ramp.
> - **Injected (owner-supplied):** 14 constructor args (`:109–112`) + `kdeg_CD3` via `set_turnover`. Live
>   values in §3.2. Two of the 14 (`span_cis_nm`, and `cis_avidity` at its live 0.0) are inert — see EQ-7.
> - **The honest headline:** the ODE genuinely delivers *avidity* (380× apparent-potency gain, §4.2) and
>   *serial-killing throughput* (43× spread across the 22 clinical engagers, §4.4). It does **NOT** deliver
>   a prozone/hook — that is **structurally absent** (§4.3), verified analytically and numerically. At
>   saturating drug the ODE collapses to `B2 = min(R_CD3, TAA_cap)` to within 0.11% (§4.5).
<!-- ATAGLANCE:END -->

# T5 — Kinetic Multivalent-Rhoden Immune Synapse (Per Cell)

> Siblings: `wholebody_pd.py` (the OWNER — builds the synapse graph, the QSS alternative kill law, the
> costim/cytokine layer), `costim_induction.py` (consumes this module's engaged fraction),
> `myeloid_il6.py` (consumes this module's engaged fraction), `multiarm_binding.py` (the same Rhoden
> geometry, format/valency layer), `kinetic_rhoden_percell.py` (the PK-side twin of the same geometry).

---

## 1. PURPOSE & DATAFLOW POSITION

### 1.1 What it does

T5 is the **PD reaction core**: for one organ, it advances a *literal, time-integrated bond ODE* for every
individual T cell in that organ's synapse graph. It answers one question per T cell per step:

> *Of this T cell's CD3 capacity, what fraction is right now committed to a **bridged trimer** with a live
> neighbouring target cell?*

That fraction, `p_eng = B2/RC` (:211), is the **single scalar that the rest of the PD model reads**. It is
not an intermediate — it *is* the model's definition of "this T cell is engaged," and it drives three
different downstream subsystems (§1.3).

The module deliberately **does not** impose quasi-steady-state (QSS). The header states the design intent
(`kinetic_synapse.py:4–7`): *"The QSS limit EMERGES automatically when koff is fast vs the PD timescale; a
slow (high-affinity) arm produces a genuine lag — no QSS assumption is imposed anywhere."* This claim is
correct and is the module's principal justification: it is verified in §4.1.

### 1.2 What feeds it

| Input | Source | Units | Line |
|---|---|---|---|
| `C_percell_T` — per-T local **free** drug | transport grid → `wholebody_pd.OrganPD.step` → `_step_kinetic` (`wholebody_pd.py:461`) | nM | :148, :153 |
| `R_CD3_nM` — per-T CD3 capacity | `wholebody_pd.py:168,237` = per-cell `CD3E` copies × `NM_PER_COPY` | nM | :113 |
| `R_TAA_nM_target` — per-target TAA capacity | `wholebody_pd.py:169,238` = per-cell antigen copies × `NM_PER_COPY` | nM | :118 |
| `Wt_norm` — antigen-weighted, row-normalized T×target apportionment | `wholebody_pd.py:269` (KDTree, `R_SYN_UM=30 µm`) | — | :116 |
| `W_incidence` — raw 0/1 T×target incidence | `wholebody_pd.py:263,245` | — | :124 |
| `syn_TAA_dens_copies` — per-T mean neighbour TAA density | `wholebody_pd.py:240` = `syn_TAA_mean / NM_PER_COPY` | copies/cell | :117 |
| `kon/koff` (both arms) | `pd_model_config.KINETIC` (:32–42), overridden per molecule by `run_tce_pd_reval.py:138–141` from `handoff/eng_params_normalized.json` | /nM/day, /day | :125–126 |
| `kint_bridge` | `pd_model_config.py:37` = 0.9/day | /day | :126 |
| `surv_target` — per-target survival | `wholebody_pd.py:468` (`exp(−kill_hazard)`), fed back each step | — | :148 |

**Unit conversion happens in the OWNER, not here** (`wholebody_pd.py:230–235`):
`kon [/nM/day] = kon [/M/s] × 86400 / 1e9`. At the generic default `kon = 1e5 /M/s` this is **8.64 /nM/day**
(RUN-VERIFY). When a `koff` is not supplied it is derived KD-consistently, `koff = kon × KD`
(`wholebody_pd.py:234–235`).

### 1.3 What it feeds — three consumers of the SAME engaged fraction

`p_eng = clip(B2 / max(RC, 1e-30), 0, 1)` (:211, and recomputed identically at `wholebody_pd.py:485`) is the
model's engagement currency. It fans out to:

1. **KILLING** — via the serial-kill hazard returned by `step()` (:234) → `wholebody_pd.py:461,467` →
   `kill_hazard` → `per_target_surv = exp(−hazard)` (`wholebody_pd.py:468`) → **fed back into the SAME
   module next step** as `surv_target`. This closed loop is what makes serial killing emergent (§4.4).
2. **COSTIM INDUCTION** — `costim_induction.CostimInduction.step(dt, p_eng)` (`costim_induction.py:73–81`):
   `da/dt = k_on·p_eng·(1−a) − k_off·a`; `R_costim(t) = R_rest·(1 + (FOLD−1)·a)`. A T cell that never
   bridges never upregulates 4-1BB/OX40/ICOS/GITR. **The tumour-conditionality of the inducible costim arms
   is therefore emergent from THIS module's B2.** (`wholebody_pd.py:485` supplies exactly `B2/RC`; the QSS
   path supplies a different surrogate, `Cb/(Cb+RA)` at `wholebody_pd.py:333` — see §5.8.)
3. **CYTOKINE / IL-6** — `myeloid_il6.MyeloidIL6.step(dt, x_T, y_T, p_eng)` (`wholebody_pd.py:486`): each
   myeloid agent integrates contact with *engaged* T cells. Separately, `engaged_dwell_rate()` (:236–239)
   returns `Σ_i p_eng_i` and is read at `wholebody_pd.py:471` as the organ cytokine drive.

**Where it sits in the life of the molecule:** the drug has already been dosed, distributed (PBPK), and
extravasated; `C_percell` is the *free* local concentration at this T cell. T5 is the last step before
biology: it converts local free drug + two receptor fields + a spatial neighbour graph into *bound trimer*,
and bound trimer into *dead targets and cytokine*.

### 1.4 Engine selection

T5 is the **CANONICAL** PD engine: `pd_model_config.py:18,29` (`PD_ENGINE = "kinetic"`). The alternative
`'qss'` path (`wholebody_pd.ternary_equilibrium`, the Schropp closed form, `wholebody_pd.py:85–105`) is
still present and reachable via `PD_ENGINE_OVERRIDE` (`run_tce_pd_reval.py:39`). **The two are not
equivalent** — most importantly on the prozone question (§4.3).

**Live PD step size:** `dt = 0.02 d` base × `pd_every = 3` (`run_tce_pd_reval.py:143,173`;
`coupled_percell_pd.py:271,273`) ⇒ **dt_PD = 0.06 d = 1.44 h**. The exact-exponential integrator (EQ-6)
exists precisely because this step is far longer than the fast bond timescales (armed-state lifetime
`1/koff_CD3 = 4.17 min` at the mosunetuzumab operating point — RUN-VERIFY).

---

## 2. GOVERNING EQUATIONS

**Notation.** `RC` = per-T CD3 capacity (nM); `B1` = armed drug·CD3 (nM); `B2` = bridged CD3·drug·TAA (nM);
`C` = local free drug (nM); `kon_*` (/nM/day); `koff_*`, `kint`, `k_hit`, `kdeg` (/day); `dt` (days);
`RTAA_j` = per-target TAA capacity (nM); `surv_j ∈ [0,1]`; `W` = raw T×target incidence; `Wt_norm` =
antigen-weighted row-normalized apportionment; `N_A = 6.02214076e23` (:34).

---

### EQ-1 — Rhoden geometric effective 2nd-arm concentration `c_eff` (kinetic_synapse.py:53–65)

```
r_Ab      = max(r_arm_nm, 1e-3)/1000                      [µm]        (:58)
Ag_bulk   = dens_copies · N/N_A · 1e9,   N = 1e9                       (:59–60)
SA_cell   = 4·π·r_cell²                                   [µm²]       (:61)
SA_Ab     = π·r_Ab²                                       [µm²]       (:62)
V_Ab      = (2/3)·π·r_Ab³                                 [µm³]       (:63)
Am_cell   = Ag_bulk/N = dens_copies/N_A · 1e9             [nmol/cell] (:64)
c_eff     = (Am_cell / SA_cell) · (SA_Ab / V_Ab) · 1e15   [nM]        (:65)
```

**Algebraic reduction** (`SA_Ab/V_Ab = π r² / ((2/3)π r³) = 3/(2 r_Ab)`), which makes the physics explicit:

```
c_eff  =  σ_2D / ( (2/3)·r_Ab )      where  σ_2D = (surface copies) / (4π r_cell²)
       =  [ number of partners inside the arm's reach disk ] / [ volume of the arm's reach hemisphere ]
       =  σ_2D · π r_Ab² / ( (2/3) π r_Ab³ )
```

- **(b) Biological meaning.** A drug molecule already anchored by one arm cannot sample the whole tissue —
  its free arm sweeps a **hemisphere of radius `r_arm`** above the membrane. The partners it can reach are
  the receptors lying in the projected disk of that hemisphere. Dividing the partner count by the swept
  volume gives the *local effective concentration* the tethered arm actually experiences. It is what turns a
  bivalent molecule's second binding event from a bimolecular into a pseudo-unimolecular reaction.
  (The only provenance the live code offers for this construction is the bare string *"Rhoden-2016"* at
  `multiarm_binding.py:21` — no DOI, no PMID. **[UNVERIFIED CITATION]**; no other author or paper is named
  anywhere in the live code, and none is asserted here.)
- **(c) Mechanistic rationale / alternative rejected.** The rejected alternative is to use the *bulk* drug
  concentration for the second arm (which would make the bridge rate `∝ C` and destroy avidity entirely),
  or to use a fitted "avidity factor." Here the avidity magnitude is *computed from geometry*: `c_eff ∝ 1/r_arm`,
  so a **longer arm DILUTES** the local concentration (bigger swept volume). The docstring states this
  explicitly (:57): *"Larger arm reach → larger explored shell → LOWER c_eff (dilution); shorter → higher."*
  This is the correct sign and is the physical hook by which the format/span sweep acts.
  `N = 1e9` (:59) **cancels exactly** between :60 and :64 (`Ag_bulk/N`) — it is a no-op scaffold, not a parameter.
- **(d) Units.** `dens_copies` [copies/cell]; `r_cell_um` [µm]; `r_arm_nm` [nm]; `c_eff` [nM].
- **CROSS-ENGINE IDENTITY.** The docstring (:55–56) asserts this is *"IDENTICAL convention to
  percell_binding.percell_ageff_nM."* `percell_binding.py` is in `_archived_2026-07-13/` and is **not live**.
  The *live* twins are `kinetic_rhoden_percell.geo_ageff_nM` (`:37–38`) and `multiarm_binding.geo_ageff_nM`
  (`:29–30`), both of which assert identity to this function. Identity of the *live* pair was not
  byte-verified in this task — **[UNVERIFIED — code comment claims it]**.
- **RUN-VERIFY (live module, this task):**
  | copies/cell | `ageff_nM(·, 8 µm, 12.5 nm)` | at 25 nm |
  |---|---|---|
  | CEACAM5 257,000 (the `NM_PER_COPY` anchor) | **63,675.7 nM** | 31,837.8 nM |
  | CD20 ~95,000 | 23,537.7 nM | 11,768.9 nM |
  | BCMA ~11,000 | 2,725.4 nM | 1,362.7 nM |
  Exactly `∝ 1/r_arm` as claimed. **Note the magnitude: 63.7 µM.** This is 10,600× the *same* 257,000
  copies expressed on the `NM_PER_COPY` basis (6.0 nM). The two bases coexist in this module and the
  consequence is §4.5.

---

### EQ-2 — Emergent synapse cleft + bridge feasibility (kinetic_synapse.py:67–72, :128–129)

```
cleft_nm  = clip( span_bridge_nm , CLEFT_MIN_NM=13.0 , CLEFT_MAX_NM=40.0 )        (:128)
g         = span_bridge_nm / max(cleft_nm, 1e-6)                                   (:71)
feas      = clip( (g − 0.6) / 0.4 , 0 , 1 )                                        (:72)
```

- **(b) Biological meaning.** Two apposed membranes are held at an intermembrane distance (the *cleft*). A
  cross-cell bridging arm must physically span it. The cleft is not a fixed number: the code lets it
  **relax toward the size of the bound complex** (:70, :127) — a bridging molecule with a long span pushes
  the membranes apart; a short one pulls them together — but bounded below by the TCR–pMHC dimension
  (13 nm) and above by the distance past which the bond is mechanically unfavourable (40 nm).
- **(c) Mechanistic rationale / alternative rejected.** The rejected alternative is a fixed cleft with a
  hard pass/fail gate. The `clip(span, 13, 40)` self-consistency means a construct is *never* penalised for
  being longer than the resting cleft (the cleft simply opens) — the penalty applies only to constructs
  **shorter than 13 nm**, which cannot open the TCR-scale gap. The `(g−0.6)/0.4` ramp then gives: `g < 0.6`
  ⇒ no bridge; `g ≥ 1.0` ⇒ full bridge; linear in between.
- **(d) Units.** `span_bridge_nm`, `cleft_nm` [nm]; `g`, `feas` [dimensionless].
- **⚠ The two ramp constants `0.6` and `0.4` (:72) are UNSOURCED shape parameters.** They are not named,
  not commented, and carry no citation. They set the entire span→feasibility transfer function. See §3.1.
- **RUN-VERIFY (live module):**
  | span (nm) | cleft (nm) | feas |
  |---|---|---|
  | **6.5** (`DART_Fc` format, `multiarm_binding.py:107`) | 13.0 | **0.0000 — hard zero, no bridge at all** |
  | **12.5** (the live default, `pd_model_config.py:38`) | 13.0 | **0.903846** |
  | 13.0 (`BiTE`/`IgG_1x1`, `multiarm_binding.py:103,105`) | 13.0 | 1.000000 |
  | 20 / 40 / 50 | 20 / 40 / 40 | 1.000000 |
  The live default span of 12.5 nm therefore runs at **90.4% feasibility, not 100%** — a fact not stated in
  the code. And a 6.5 nm span is a *hard zero*: this ramp makes short-span formats **unkillable**, not
  merely weaker. That is a strong structural claim resting on two unsourced constants.

---

### EQ-3 — Trans effective ALIVE TAA concentration (kinetic_synapse.py:131, :164–166)

```
c_eff,trans   = ageff_nM(dens, R_CELL_UM, span_bridge_nm) · feas        (nM, fixed at construction)  (:131)
alive_frac_i  = clip( Σ_j Wt_norm[i,j] · surv_j , 0 , 1 )               (per T cell, each step)      (:164–165)
Tfree_i       = c_eff,trans_i · alive_frac_i                            (nM)                         (:166)
```

- **(b) Biological meaning.** `c_eff,trans` is the local TAA concentration the armed CD3 arm samples across
  the cleft (EQ-1 × EQ-2). `alive_frac` is the antigen-weighted fraction of this T cell's synapse partners
  that are **still alive** — as the T cell's neighbours die, the bridgeable antigen in its synapse falls.
  `Tfree` is the effective concentration of *live, reachable* TAA.
- **(c) Mechanistic rationale.** This is the **first of two channels** by which killing feeds back on
  engagement (the second is the TAA cap, EQ-10). Because `Wt_norm` is row-normalized and antigen-weighted,
  `alive_frac` is a *survival-weighted mean* over the T cell's neighbours, not a raw count — a T cell
  surrounded by one high-antigen survivor is not treated as "empty."
- **(d) Units.** `c_eff,trans`, `Tfree` [nM]; `alive_frac` [dimensionless, 0–1].
- **RUN-VERIFY:** at the canonical config (257k-copy target, 8 µm cell, 12.5 nm span):
  `c_eff,trans = 63,675.7 × 0.903846 = ` **57,553.0 nM = 57.6 µM**. (Live-module value, exact.)

---

### EQ-4 — CD3 loading rate and bridge-forming rate (kinetic_synapse.py:167–168)

```
rate_on_i = kon_CD3 · C_i                    (/day)   — CD3 arm loading (pseudo-first-order in drug)  (:167)
kf_i      = kon_TAA · Tfree_i                (/day)   — 2nd-arm (bridge) forming rate                 (:168)
```

- **(b) Biological meaning.** `rate_on` is the rate at which a free CD3 receptor on this T cell captures a
  drug molecule from the local free pool. `kf` is the rate at which an *already-armed* drug closes its
  second arm onto a neighbouring TAA — pseudo-**uni**molecular, because the partner concentration is the
  *tethered* `c_eff`, not the bulk.
- **(c) Mechanistic rationale.** This asymmetry (`rate_on ∝ C` bimolecular; `kf ∝ c_eff` unimolecular) is
  exactly what avidity means, and it is the reason the trimer is so much more potent than either KD implies
  (§4.2).
- **(d) Units.** `kon_*` [/nM/day]; `C`, `Tfree` [nM]; `rate_on`, `kf` [/day].
- **RUN-VERIFY (canonical, mosunetuzumab-like):** `kon = 8.64 /nM/day`, `Tfree(alive) = 57,553 nM` ⇒
  **`kf = 497,258 /day`**. This is an enormous forward rate — see §4.5.

---

### EQ-5 — THE CORE: the two-state bond ODE (kinetic_synapse.py:169–182)

```
dB1/dt = rate_on·(RC − B1 − B2)  −  koff_CD3_eff·B1  −  kf·B1  +  koff_TAA·B2      (:170)
dB2/dt =                             kf·B1           −  koff_TAA·B2  −  kint·B2     (:171)
```

recast as the linear system `dX/dt = M·X + b`, `X = [B1, B2]ᵀ`:

```
m11 = −(rate_on + koff_CD3_eff + kf)                                               (:178)
m12 =  (koff_TAA − rate_on)                                                        (:179)
m21 =  kf                                                                          (:180)
m22 = −(koff_TAA + kint)                                                           (:181)
b   = [ rate_on · RC , 0 ]ᵀ                                                        (:182)
```

- **(b) Biological meaning — the trimer, arm by arm.** Three species exist physically:
  free CD3 (= `RC − B1 − B2`), the **armed** T cell (`B1` = drug held by CD3 only, dangling its TAA arm into
  the cleft), and the **bridged trimer** (`B2` = CD3·drug·TAA, the physical immune synapse bond). The system
  is written in only 2 states because free CD3 is not independent — it is `RC − B1 − B2`, which is why
  `rate_on` appears on BOTH the `m11` diagonal AND in `m12` (:176–177, :179): a molecule sitting in `B2` is
  *also* consuming a CD3 receptor, so it depletes the free pool that feeds `B1`.
- **(c) Mechanistic rationale — the KEY avidity construction (:172–175), and the alternative rejected.**
  Look at `m22`: the trimer decays at `(koff_TAA + kint)` — **`koff_CD3` is deliberately ABSENT**. The code
  states why (:172–175): *"from the BRIDGED trimer, drug can only leave the synapse by the TAA arm releasing
  (koff_TAA → back to armed B1) or by internalization (kint). The CD3 arm releasing does NOT dissolve the
  trimer (drug still held by TAA) — so the trimer lifetime is set by the SLOWER (TAA) arm, i.e. bivalent
  avidity. This is why B2 must NOT decay at koff_CD3."*

  This is **correct bivalent physics and it is the single most important line in the module.** The rejected
  alternative — letting `B2` decay at `koff_CD3 + koff_TAA` — would make the trimer a mere product of two
  independent monovalent equilibria and would destroy avidity. It is also *the* reason `k_hit` has anything
  to race against (EQ-13): the synapse persists on the TAA clock while the CD3 arm cycles.

  **But note the honest consequence:** the CD3 arm's *release from the trimer* is not represented at all.
  A CD3 release event from `B2` should return the complex to a "TAA-only bound" state, from which it could
  re-arm a *different* CD3 — that state does not exist here. §4.3 shows this omission is exactly what
  removes the prozone.
- **(d) Units.** `B1`,`B2`,`RC` [nM]; `rate_on`,`kf`,`koff_*`,`kint` [/day]; `m*` [/day]; `b` [nM/day].
- **RUN-VERIFY (canonical operating point):**
  - trimer mean lifetime `1/(koff_TAA + kint) = 1/(12.528 + 0.9) = 0.0745 d = ` **1.79 h**
  - armed-state lifetime `1/koff_CD3 = 1/345.6 d = ` **4.17 min**
  - the trimer therefore outlives the CD3 bond by **26×** — this ratio *is* the avidity.

---

### EQ-6 — Exact 2×2 matrix-exponential integrator (kinetic_synapse.py:74–102)

```
tr   = m11 + m22;   det = m11·m22 − m12·m21                                        (:78–79)
disc = max(tr² − 4·det, 0);   s = √disc                                             (:80–81)
λ1,2 = ½(tr ± s)                                                                    (:82)
e_k  = exp( clip(λ_k·dt, −50, +50) )                                                (:83)

Sylvester (distinct λ):   a1 = (e1 − e2)/(λ1 − λ2)                                  (:88)
                          a0 = (λ1·e2 − λ2·e1)/(λ1 − λ2)                            (:89)
Degenerate (|λ1−λ2|<1e-9): a1 = dt·e1 ;  a0 = e1 − λ1·dt·e1                          (:88–89)

expm(M·dt) = a0·I + a1·M                                                            (:90)
X(t+dt)    = expm(M·dt)·X(t)  +  M⁻¹·(expm(M·dt) − I)·b                             (:92–101)
M⁻¹        = (1/det)·[[m22, −m12], [−m21, m11]]                                     (:95, :99–101)
```

- **(b) Biological meaning.** None — this is pure numerics. It exists because the biology it integrates has
  a **stiffness ratio of ~10⁵**: `kf ≈ 5×10⁵/day` versus `kint = 0.9/day`, over a PD step of 0.06 d. An
  explicit Euler step would need `dt < 2/kf ≈ 4×10⁻⁶ d` (≈ 0.35 s) for stability — 15,000× smaller than the
  PD step actually used.
- **(c) Mechanistic rationale / alternative rejected.** The rejected alternative is exactly that: sub-stepping,
  or an implicit solver. The closed-form 2×2 exponential is **exact** for constant `M` and `b` over `dt`, is
  unconditionally stable for the stable `M` here (`tr < 0`, `det > 0`), costs O(1) per T cell, and is
  vectorized across all T cells at once. The `a0/a1` Sylvester decomposition is the standard 2×2 result and
  is **correct as written** (verified by hand this task against `f(M) = [f(λ1)(M−λ2 I) − f(λ2)(M−λ1 I)]/(λ1−λ2)`;
  the degenerate branch reduces correctly to `e^{λt}(I + (M−λI)t)`).
- **(d) Units.** `λ` [/day]; `dt` [days]; `a0` [dimensionless]; `a1` [days]; `X` [nM]; `b` [nM/day].
- **⚠ APPROXIMATION — the `max(·, 0)` on the discriminant (:80) SILENTLY DISCARDS COMPLEX EIGENVALUES.**
  `disc = (m11−m22)² + 4·m12·m21`, and `m12 = koff_TAA − rate_on` **goes negative whenever
  `rate_on > koff_TAA`** (i.e. at any appreciable drug level). When `disc < 0` the code clamps it to 0, which
  forces `λ1 = λ2 = tr/2` and routes the step through the **degenerate real branch** (`near` is then true at
  :85) — the oscillatory part of the true solution is silently dropped.
- **The exact condition (derived this task, not asserted).** Writing `u = rate_on − koff_TAA` and
  `d = koff_CD3_eff − kint`:
  ```
  disc  =  (u − kf)²  +  d²  +  2·d·(u + kf)
  ```
  so `disc < 0` requires **BOTH** (i) `d < 0`, i.e. **`koff_CD3_eff < kint`**, and (ii) `u ≈ kf`, i.e. drug
  tuned to the **resonance** `rate_on ≈ kf + koff_TAA`. At the live `kint = 0.9/day` and `kon = 8.64/nM/day`,
  condition (i) means **`KD_CD3 < 0.104 nM`**.
  - **RUN-VERIFY (this task).** On a grid that *resolves the resonance* (spans {13…40} nm × copies
    {10³…10⁶} × KD_CD3 {0.01…1000} nM × KD_TAA {0.01…100} nM × C on a log grid **refined around
    `C* = (kf + koff_TAA)/kon`**; 97,020 points): **negative discriminant at 2,844 points (2.93%)**, worst
    normalized excursion **`min(disc/tr²) = −1.21×10⁻³`**.
  - **⚠ A REAL MOLECULE SITS IN THE GATING REGIME.** Of the 22 clinical engagers, exactly one satisfies
    condition (i): **cevostamab** (`koff_CD3 = 0.285/day < kint = 0.9/day`; `KD_CD3 = 0.033 nM`,
    `rundir/handoff/eng_params_normalized.json`). At cevostamab's arm kinetics the worst excursion is
    `disc/tr² = −9.9×10⁻⁴`, reached at ~3×10² TAA copies/cell and **C ≈ 70 nM — a physiologically reachable
    plasma level**, not an extrapolation.
  - **Verdict:** the *magnitude* of the error is small (|Im λ|/|Re λ| ≲ 2%, so the discarded oscillation is a
    minor perturbation on a mode that is fully relaxed within one PD step anyway). But the clamp is **not**
    "essentially never active": it is gated by a sharp, checkable condition that one molecule in the
    validation set meets. It is an undocumented approximation, and any future high-CD3-affinity construct
    (`KD_CD3 < 0.1 nM`) enters this regime **by design**, which is exactly the corner a costim/affinity sweep
    will explore. It should be documented in the code and, if the sweep goes there, replaced with a proper
    complex-eigenvalue branch.

---

### EQ-7 — cis-avidity retention of the CD3 arm (kinetic_synapse.py:132–134)

```
koff_CD3_eff = koff_CD3 · ( 1 − clip(cis_avidity, 0.0, 0.95) )                     (:134)
```

- **(b) Biological meaning.** In a trispecific (CD3 × TAA × costim), the costim arm binds a receptor on the
  **same** T cell as the CD3 arm (a *cis* interaction). That second same-cell contact holds the drug on the
  T cell even when the CD3 bond breaks, so the *effective* CD3 detachment rate falls. This is the avidity
  retention that a costim×CD3 arm pair buys, and it is the hook the trispecific design sweep is built on.
- **(c) Mechanistic rationale / HONEST CAVEAT.** The intended mechanism is stated at :132–133:
  *"avidity factor = c_eff,cis/(c_eff,cis + KD_costim-ish); passed in as cis_avidity in [0,1)."* — i.e. the
  *intended* form is a mechanistic occupancy computed from the cis span geometry.
  **But that computation does not happen in this module.** `cis_avidity` arrives as a **bare
  dimensionless number** from the caller and is applied as a **direct multiplicative scalar on koff** (:134).
  `span_cis_nm` is accepted by the constructor (:111) and **never used anywhere in the file** —
  grep-verified: it appears at :40 (default), :111 (signature), and nowhere else. So the geometry→avidity
  link that the docstring promises (:23–26) is, in this module, **not wired**.
  This is a scalar-on-koff, which is precisely the class of construction a mechanistic model should avoid.
  It is **inert at the live operating point** (`cis_avidity = 0.0`, `pd_model_config.py:40` ⇒
  `koff_CD3_eff = koff_CD3`, byte-identical), so nothing published to date depends on it — but it must not
  be switched on in its current form without first computing the occupancy from `span_cis_nm`.
  The geometry needed to do so is already available in the live graph (`multiarm_binding.geo_ageff_nM`,
  `:29–30`, plus `multiarm_binding._cis_feasibility`, which `wholebody_pd.py:208` already imports).
- **(d) Units.** `cis_avidity` [dimensionless, 0–0.95]; `koff_CD3_eff` [/day].
- **Ceiling:** the `0.95` clamp caps the retention at a **20× koff reduction**. Unsourced.

---

### EQ-8 — Free-CD3 receptor turnover (kinetic_synapse.py:143–146, :157–161)

```
RC_free      = max( RC − (B1 + B2) , 0 )                                            (:158)
KSYN         = RC0 · kdeg_CD3                                                       (:159)
RC_free_new  = ( RC_free + dt·KSYN ) / ( 1 + dt·kdeg_CD3 )      [backward Euler]    (:160)
RC           = RC_free_new + (B1 + B2)                                              (:161)
```

- **(b) Biological meaning.** CD3 is synthesized and degraded. The *free* pool relaxes toward a set-point
  `RC0` (the initial capacity, :114) at rate `kdeg_CD3`, with zero-order synthesis `KSYN = RC0·kdeg` chosen
  so the drug-free steady state reproduces `RC0` exactly.
- **(c) Mechanistic rationale.** Backward Euler (unconditionally stable, positivity-preserving) is used
  because turnover is *slow* relative to the binding solve — the code calls this out as an operator split
  (:154). Bound CD3 (`B1+B2`) is explicitly **not** degraded here; its fate is `koff` release or `kint`
  internalization in the binding solve (:155–156).
- **(d) Units.** `kdeg_CD3` [/day]; `RC`, `RC0`, `KSYN·dt` [nM]; `dt` [days].
- **⚠ THIS ENTIRE BLOCK IS DEAD AT THE LIVE OPERATING POINT.** `wholebody_pd.py:247,433` calls
  `set_turnover(float(c.get('kdeg_CD3_perday', 0.0)))`, and **`kdeg_CD3_perday` is not a key in
  `pd_model_config.KINETIC`** (`pd_model_config.py:32–42` — grep-verified, the only `kdeg_CD3` occurrences
  in the live graph are inside this module plus those two `.get(..., 0.0)` calls). Therefore
  `self.kdeg_CD3 = 0.0` always, the `if self.kdeg_CD3 > 0.0` guard at :157 never fires, and **CD3 capacity
  is static in every published run.**
- **Behavioural note if it were switched on:** at steady state the recursion gives `RC = RC0 + (B1+B2)`
  (algebra, this task) — i.e. total capacity grows to accommodate bound drug, so drug binding would **not**
  deplete the free CD3 pool at steady state. That is a defensible modelling choice (synthesis replaces
  bound receptor) but it is not what "receptor turnover" usually means, and it should be stated before use.

---

### EQ-9 — CD3-side conservation clamp (kinetic_synapse.py:184–189)

```
tot = B1 + B2
if tot > RC:   scale = RC/max(tot,1e-30);   B1 ← B1·scale;   B2 ← B2·scale          (:185–188)
B1 ← max(B1, 0);   B2 ← max(B2, 0)                                                  (:189)
```

- **(b) Biological meaning.** A T cell cannot have more drug bound than it has CD3 receptors.
- **(c) Mechanistic rationale.** The exact exponential (EQ-6) *should* respect `B1+B2 ≤ RC` automatically for
  the true continuous system — the receptor-depletion terms in `M` enforce it. The clamp is a **numerical
  guard** against (i) the coefficient-freezing operator split, (ii) the `disc` clamp of EQ-6, and (iii) the
  `RC` update of EQ-8. Proportional rescaling (rather than truncating one state) preserves the `B1:B2` ratio
  and therefore does not bias the engaged fraction.
- **(d) Units.** all [nM]; `scale` [dimensionless].

---

### EQ-10 — TAA-side (two-sided) conservation cap — **THE ABUNDANCE GRADER** (kinetic_synapse.py:190–199)

```
taa_cap_i = Σ_j W[i,j] · RTAA_j · surv_j                     (nM, per T cell)       (:194)
if B2_i > taa_cap_i:
    excess_i = B2_i − taa_cap_i                                                     (:197)
    B2_i     = min(B2_i, taa_cap_i)                                                 (:198)
    B1_i     = B1_i + excess_i          # drug returns to ARMED (CD3 still holds it) (:199)
```

- **(b) Biological meaning.** A T cell cannot bridge more TAA than physically exists, alive, inside its
  synapse. When the trimer overshoots the reachable live antigen, the excess drug is **not destroyed** — it
  falls back to the armed state `B1`, because the CD3 arm still holds it. Total bound drug `B1+B2` is
  therefore exactly preserved by this operation (verified by inspection: `−excess` from `B2`, `+excess` to
  `B1`), so EQ-9's `B1+B2 ≤ RC` invariant survives it.
- **(c) Mechanistic rationale — and the honest admission the code itself makes.** The rationale is stated
  verbatim at :119–123: *"Without this cap the bridge is limited only by CD3 (RC) and the huge geometric ceff
  saturates every synapse → abundance signal lost (low-copy BCMA killed like high-copy CD20). WITH it,
  high-abundance → CD3-limited (full), low → TAA-limited (partial). This is the linear-in-R_TAA behaviour
  the Schropp QSS engine had."*

  **Read that carefully.** The cap exists because the *rate law* (EQ-4/EQ-5) has lost the ability to grade by
  antigen abundance — `kf = 497,258/day` swamps every other rate by 4–5 orders of magnitude, so mass action
  drives `B2 → RC` for essentially any target. The cap is a **hard `min()` clamp bolted on to restore a
  behaviour the kinetics no longer produce.** It works, and the resulting abundance grading is correct and
  necessary — but it is **IMPOSED, not emergent**, and §4.5 shows it dominates the model's antigen response.
  The root cause is the two-basis mismatch of §4.5, not a flaw in the clamp itself.
- **(d) Units.** `taa_cap`, `excess`, `B1`, `B2` [nM]; `W` [0/1]; `surv` [0–1].
- **Note on the basis:** `taa_cap` sums `RTAA_j` (which are on the **`NM_PER_COPY` basis**, 6.0 nM for a
  257k-copy cell) while `B2` was produced by a rate law whose partner concentration is on the **Rhoden
  `c_eff` basis** (57,553 nM for the same cell). These are different physical quantities (bulk capacity vs
  local tethered concentration), which is legitimate — but the cap is where the two bases collide.

---

### EQ-11 — Survivor-renormalized apportionment (kinetic_synapse.py:200–208)

```
Ws[i,j]   = Wt_norm[i,j] · surv_j                                                   (:205)
rs_i      = Σ_j Ws[i,j] ;   rs_i ← 1.0 where rs_i < 1e-12                           (:206)
Ws        ← diag(1/rs) · Ws                # rows renormalized onto SURVIVORS       (:207)
has_live_i = 1 if Σ_j Ws[i,j] > 1e-9 else 0                                         (:208)
```

- **(b) Biological meaning.** A T cell's killing output must be delivered to *live* neighbours. As its
  neighbours die, its lethal output **concentrates on the survivors** rather than being wasted on corpses.
  `has_live` zeroes the output entirely for a T cell whose whole synapse is dead.
- **(c) Mechanistic rationale / alternative rejected.** The rejected alternative — apportioning by the static
  `Wt_norm` — would spend a T cell's kill capacity on already-dead cells, artificially flattening the kill
  curve at high depletion. Renormalizing *per step* is what makes the last surviving target in a synapse
  receive the T cell's full attention. Together with EQ-3 (`alive_frac`) and EQ-10 (`taa_cap`, which also
  carries `surv_j`), this is the third of three places survival feeds back.
- **(d) Units.** dimensionless throughout; `Ws` rows sum to 1 over live neighbours (or 0 if none).

---

### EQ-12 — Engaged fraction (kinetic_synapse.py:211, :236–239)

```
p_eng_i        = clip( B2_i / max(RC_i, 1e-30) , 0 , 1 )                            (:211)
engaged_dwell_rate() = Σ_i clip(B2_i/max(RC_i,1e-30), 0, 1)                         (:239)
bridged_total()      = Σ_i B2_i                                                     (:242)
```

- **(b) Biological meaning.** The fraction of this T cell's CD3 that is committed to a bridged synapse —
  the model's operational definition of "engaged." `engaged_dwell_rate()` is the organ-level sum, read as
  the cytokine emission proxy (`wholebody_pd.py:471`).
- **(c) Mechanistic rationale.** Normalizing by `RC` rather than reporting raw `B2` makes the quantity a
  *per-cell state* (0–1) rather than an extensive amount, so it can drive per-cell processes (costim
  induction memory; myeloid contact activation) that saturate at "fully engaged."
- **(d) Units.** dimensionless (0–1); `engaged_dwell_rate` [T-cell-equivalents]; `bridged_total` [nM·cells].
- **`bridged_total()` (:241–242) has ZERO callers in the live import graph** (grep-verified) — a diagnostic
  accessor only.

---

### EQ-13 — Serial-killing throughput: the engage / hit / detach RACE (kinetic_synapse.py:209–228)

```
P(productive hit before detach)  =  k_hit / (k_hit + koff_CD3_eff)
cycling rate (detach & re-engage) ≈ koff_CD3_eff

serial_rate_i =  koff_CD3_eff · [ k_hit/(k_hit + koff_CD3_eff) ] · p_eng_i · has_live_i
              =  [ k_hit·koff_CD3_eff / (k_hit + koff_CD3_eff + 1e-30) ] · p_eng_i · has_live_i   (:228)
```
(targets killed per day per T cell — the **harmonic mean** of the two competing clocks.)

- **(b) Biological meaning.** Once bridged, **two clocks race.** One is the lethal hit: an engaged CTL needs
  sustained contact to polarize its granules and deliver a lethal dose (`k_hit`). The other is detachment
  (`koff_CD3_eff`). Only if the hit lands before detachment is the engagement *productive*. But the T cell
  must **also** detach in order to move on to the next target — so the throughput is the product of
  (cycling rate) × (probability each cycle is productive). Hence the harmonic form.
- **(c) Mechanistic rationale / limits (the code states them at :220–222).**
  - **Slow koff** (very high CD3 affinity): every engagement is productive (`P → 1`) but the T cell is
    **glued to one target and cannot cycle** ⇒ `serial_rate → koff` (low).
  - **Fast koff** (low CD3 affinity): the T cell cycles constantly but **detaches before hitting**
    (`P → 0`) ⇒ `serial_rate → k_hit` (saturating).
  This is the *serial-killing* rationale for why CD3 affinity should not simply be maximized, and it is the
  single mechanistic reason the model can distinguish a CD3-detuned engager from an affinity-matured one.
- **(d) Units.** `k_hit`, `koff_CD3_eff`, `serial_rate` [/day]; `p_eng`, `has_live` [dimensionless].
- **⚠ STRUCTURAL NOTE — this is an ALGEBRAIC race, not an ODE.** There is **no hit-state compartment**.
  `k_hit` is not a state variable and no "delivering a lethal hit" species is integrated. The race
  probability is a closed-form steady-state expression multiplied by the instantaneous `p_eng`. This is a
  *quasi-steady-state* treatment of the hit process — inside a module whose entire premise (:4–7) is *"no
  QSS assumption is imposed anywhere."* The bond kinetics are literal; the **hit kinetics are not**.
- **⚠ THE IN-CODE "WINDOW OPTIMUM" CLAIM IS FALSIFIED BY ITS OWN FORM.** Line 212 asserts the race is *"the
  mechanism that makes CD3 affinity a WINDOW optimum, not a monotone knob."* It is not: the harmonic
  `k_hit·koff/(k_hit+koff)` is **monotonically increasing and saturating in koff, with no interior maximum**
  (elementary; `d/dkoff = k_hit²/(k_hit+koff)² > 0` always). The code's own later text concedes this
  (:224–226: *"the efficacy endpoint is a broad plateau, not a sharp koff~k_hit optimum"*). See §5.2 for the
  full RUN-VERIFY. The *window* (kill ÷ cytokine dwell) may still be non-monotone, but that is a property of
  the ratio computed in `wholebody_pd`, **not of this equation.**
- **RUN-VERIFY (live module + the 22 clinical engagers' measured `koff_CD3` from
  `rundir/handoff/eng_params_normalized.json`, at `k_hit = 12/day`):**
  | molecule | measured `koff_CD3` (/s) | `koff_CD3` (/day) | serial ceiling (/day) | P(hit \| engaged) |
  |---|---|---|---|---|
  | cevostamab | 3.3e-6 | 0.285 | **0.279** | 0.977 |
  | talquetamab / forimtamig | 1.9e-4 | 16.4 | 6.93 | 0.422 |
  | glofitamab | 4.5e-4 | 38.9 | 9.17 | 0.236 |
  | tarlatamab | 1.49e-3 | 128.7 | 10.98 | 0.085 |
  | teclistamab | 3.30e-3 | 285.0 | 11.52 | 0.040 |
  | **mosunetuzumab** (the calibration anchor) | 4.0e-3 | **345.6** | **11.597** | 0.034 |
  | blinatumomab / solitomab | 2.6e-2 | 2246.4 | 11.94 | 0.005 |
  | REGN5459 (deliberately CD3-detuned) | 5.0e-2 | 4320.0 | 11.97 | 0.003 |

  **This is where the equation earns its keep.** The serial ceiling spans **0.279 → 11.97 /day, a 43× range**,
  driven purely by the measured CD3 off-rate. Cevostamab (KD_CD3 = 0.033 nM) is *throughput-crippled*: it
  binds CD3 so tightly it cannot cycle. That is a real, mechanistically-derived, clinically-relevant
  prediction that the QSS engine structurally cannot make. **Note also that the mosunetuzumab ceiling
  (11.597/day) is the "11.6/day" quoted in `handoff/kinetic_calib.json` as the PRIMARY `k_death` anchor** —
  reproduced exactly here.

---

### EQ-14 — Engaged-dwell accumulation (kinetic_synapse.py:229–231)

```
dwell_engaged_i  +=  dt · p_eng_i                        (T-cell-days)              (:231)
```

- **(b) Biological meaning.** Cumulative time this T cell has spent engaged — the intended
  *synapse-stability* cytokine driver (a slow-koff, over-stable synapse spends more time engaged per kill
  and therefore emits more cytokine per target killed: the over-stable-synapse CRS mechanism, :229–230).
- **(c) Mechanistic rationale.** This is the CRS axis of the design window: efficacy scales with *kills*,
  toxicity with *engaged time*. Separating them is the whole point.
- **(d) Units.** [T-cell-days].
- **⚠ `dwell_engaged` IS WRITE-ONLY — it has ZERO readers in the live import graph** (grep-verified across
  `engine/*.py`). The cytokine drive actually used is `engaged_dwell_rate()` (:236–239, EQ-12), which returns
  the **instantaneous** `Σ p_eng`, **not** the accumulated dwell. So the cumulative-dwell CRS mechanism
  described at :229–230 is **built but not wired**. (This is not necessarily wrong — `wholebody_pd` integrates
  the instantaneous rate over time itself at `wholebody_pd.py:478` — but the module's own dwell accumulator
  is dead code and should either be wired or removed.)

---

### EQ-15 — Kill-hazard apportionment (kinetic_synapse.py:232–234)

```
dkill_tgt_j  =  Σ_i Ws[i,j] · serial_rate_i                    (targets/day, per target)  (:233)
return          dt · k_death · dkill_tgt                       (hazard increment)         (:234)
```
Consumed at `wholebody_pd.py:461–468`: `kill_hazard[tgt] += dH·g_eff·treg_damp`;
`per_target_surv = exp(−kill_hazard)` → **fed straight back into this module's `surv_target` next step.**

- **(b) Biological meaning.** Each T cell's productive serial-kill output is distributed across the live
  targets in its synapse, in proportion to their (survival-weighted, antigen-weighted) share. Accumulated as
  a hazard, so target survival is `exp(−H)`.
- **(c) Mechanistic rationale.** The hazard formulation (rather than a direct fractional kill) makes the
  killed fraction depend on **both** drug level and **time**, so a low-exposure organ reaches an incomplete,
  drug-graded plateau rather than the same endpoint more slowly.
- **(d) Units.** `serial_rate` [/day]; `dt` [days]; `dkill_tgt` [/day]; hazard increment [dimensionless].
- **⚠ `k_death` CARRIES DIFFERENT DIMENSIONS IN THE TWO ENGINES.** Here `dkill_tgt` is already a **rate**
  [/day], so `dt·k_death·dkill_tgt` is dimensionless only if **`k_death` is dimensionless**. In the QSS path
  (`wholebody_pd.py:362`) the identical expression `dt·k_death·dkill` is applied to a `dkill` that is a
  **concentration** [nM] (`wholebody_pd.py:350–355`: `dkill` is built from `Cb`, the trimer concentration),
  so there `k_death` must have units [1/(nM·day)]. The same symbol, loaded from the same file
  (`handoff/kinetic_calib.json`), means two different things. It is **numerically silent** because
  `k_death = 1.0` in both — but it means the calibrated `k_death` is **not transferable between the two
  engines**, and a future re-fit on one path must not be quoted for the other.

---

## 3. PARAMETERS OWNED

### 3.1 Module-level constants (defined IN `kinetic_synapse.py`)

| Symbol | Value | Units | Line | Provenance tag | Source (verbatim from code, where any) | Mechanistic rationale / audit |
|---|---|---|---|---|---|---|
| `AVO` | 6.02214076e23 | /mol | :34 | **[EXACT — SI-defined constant]** (NOT `[MEASURED]`: it is a definition, not a measurement, and needs no citation) | none in code, none required | Avogadro constant. The value is exact by the SI definition of the mole; it is not a model parameter and carries no provenance risk. |
| `NM_PER_COPY` | `6.0/257000` = **2.3346303501945526e-5** | nM per receptor copy | :37 | **[FITTED: back-derived from the tumor kill calibration]** | code comment (:36–37): *"pinned by the validated tumor (CEACAM5 257,000 copies → Rcap_TAA 6.0 nM). SAME constant as wholebody_pd.NM_PER_COPY."* `wholebody_pd.py:75–83` elaborates: *"ONE physical constant: the synapse reaction volume … (equivalently a ~71 pL synapse …)"* | **This is the most load-bearing and most attackable constant in the module.** It is **not** an independently measured physical volume: it is `Rcap_TAA / (CEACAM5 copies)`, where `Rcap_TAA = 6.0` is a *calibration scale* imported from the prior tumor ABM (`wholebody_pd.py:74`). RUN-VERIFY of the implied volume: 1 copy → 2.3346e-14 M ⇒ `V = (1/N_A)/2.3346e-14 = ` **7.11e-11 L = 71.1 pL = 71,100 µm³** — a **sphere of radius 25.7 µm.** For comparison, a physical immune-synapse cleft (a ~3 µm-radius contact disc × the module's own 13 nm cleft, EQ-2) is ≈ **0.37 µm³** — the "synapse reaction volume" is **~2×10⁵× larger than a synapse.** It is a fitted concentration scale wearing a physical name. It is *self-consistent* (the same constant converts CD3 and TAA), and it correctly preserves absolute abundance ordering — which is what it is really for — but it must **not** be described to a committee as a measured synapse volume. |
| `R_CELL_UM` | 8.0 | µm | :38 | **[UNVERIFIED CITATION]** | code comment: *"target-cell radius (µm), Rhoden default"* — **no PMID, no DOI, no year.** | Sets `SA_cell = 4π(8)² = 804 µm²` in EQ-1, hence `c_eff ∝ 1/r_cell²`. The only corroborating string in the live code is *"Rhoden-2016"* (`multiarm_binding.py:21`) — also with no DOI. **No source for the 8.0 µm value exists in any live code file, and none is asserted here.** Sensitivity is quadratic and one-line-checkable: **halving `r_cell` quadruples `c_eff`.** Whether 8 µm is right for a B-lymphoblast target is an **open question for the committee — this doc does not supply an uncited "true" radius to compare it against** (see §5.11). |
| `SPAN_BRIDGE_DEFAULT_NM` | 12.5 | nm | :39 | **[ASSUMED]** | comment: *"trans CD3↔TAA arm span (nm); AF3/format override per construct"* — no source | Feeds EQ-1 (`c_eff ∝ 1/span`) and EQ-2 (feasibility). **The promised per-construct AF3/format override is NOT WIRED in the live runner:** `run_tce_pd_reval.py:137–141` overrides only `kon_*`/`koff_*` in `_kpm`; `span_bridge_nm` is never touched, so **every one of the 22 clinical engagers runs at 12.5 nm** regardless of format (BiTE vs IgG vs DART). Consequently `feas = 0.9038` for all of them and the geometry axis of the sweep is, at present, **unexercised**. |
| `SPAN_CIS_DEFAULT_NM` | 12.5 | nm | :40 | **[ASSUMED — and UNUSED]** | comment: *"cis CD3↔costim arm span (nm)"* | **Accepted by the constructor (:111) and never read anywhere in the file** (grep-verified). See EQ-7: the cis geometry→avidity computation the docstring promises is not implemented here. |
| `CLEFT_MIN_NM` | 13.0 | nm | :41 | **[ASSUMED: TCR–pMHC dimension]** | comment: *"immune-synapse cleft floor (~TCR-pMHC dimension)"* — no PMID | The ~13–15 nm TCR–pMHC span is a well-known structural figure, but **no citation is present in the code**. Load-bearing: it is what makes a 6.5 nm span a hard zero (EQ-2). |
| `CLEFT_MAX_NM` | 40.0 | nm | :42 | **[ASSUMED]** | comment: *"cleft ceiling before bond is mechanically unfavourable"* — no source | Sets the span above which no further cleft relaxation is allowed. Never binding at the live 12.5 nm span. |
| feasibility ramp `0.6` | 0.6 | — | :72 | **[UNSOURCED — TBD]** | *no comment at all* | Below `0.6×cleft` the bridge probability is **zero**. This single unnamed number decides whether short-span formats are killable. |
| feasibility ramp `0.4` | 0.4 | — | :72 | **[UNSOURCED — TBD]** | *no comment at all* | Width of the linear ramp from 0 → 1. Together with `0.6`, defines the entire span→feasibility transfer function. |
| `K_HIT_DEFAULT` | **12.0** | /day | :48 | **[UNSOURCED — TBD]** ⚠ | comment (:43–51), verbatim: *"Literature: a CTL delivers a lethal hit in ~tens of min once engaged; serial killing ~ few targets/day at saturating engagement. … K_HIT_DEFAULT = 12.0 /day (~1 lethal hit per 2 h engaged) — **FIXED from serial-killing literature, NOT fitted.**"* — **NO PMID, NO DOI, NO AUTHOR, NO YEAR is given anywhere in the file.** | **This is the exact pattern the provenance audit exists to catch: a constant that ASSERTS a literature provenance while naming no literature.** The claim "NOT fitted" cannot be checked against a source that is not cited. `handoff/kinetic_calib.json` (`k_death_provenance`) refers to a *"Halle 2-16/day band"* — a bare author surname; **no PMID/DOI/year appears in any code or params file**, so it too is **[UNVERIFIED CITATION]** (this doc does **not** guess at which paper is meant). `k_hit` sets the ceiling of every kill rate in the model (EQ-13). **It must be sourced with a PMID or explicitly relabelled as an assumption before this goes to committee.** |

**Numerical guards (code-internal, not parameters, but they shape behaviour):**
`1e-3` nm floor on `r_arm` (:58); `1e-6` nm floor on `cleft` (:71); `1e-9` eigenvalue-degeneracy threshold
(:85); `±50` exponent clip (:83); `1e-30` determinant guard (:99); `1e-30` `RC` guard (:187, :211, :239);
`1e-12` row-sum guard (:206); `1e-9` `has_live` threshold (:208); `1e-30` harmonic-denominator guard (:228);
`0.95` `cis_avidity` clamp (:134). All are [CODE-INTERNAL]. The `±50` exponent clip is a genuine safety net:
at `λ = −(rate_on + koff + kf) ≈ −5×10⁵/day` and `dt = 0.06 d`, `λ·dt ≈ −3×10⁴`, so **the clip fires on every
step** — correctly, since `exp(−3e4)` underflows to 0 anyway and the fast mode is fully relaxed within one PD
step. (This is the mathematical statement of "the QSS limit emerges automatically," EQ-6 / §4.1.)

### 3.2 Injected parameters (owned by the caller; live values traced to source)

| Symbol | Live value | Units | Set at | Provenance tag | Notes |
|---|---|---|---|---|---|
| `kon_CD3` | `1e5 /M/s` → **8.64 /nM/day** | /nM/day | `pd_model_config.py:33`; converted `wholebody_pd.py:231` | **[ASSUMED: "standard mAb assoc"]** — code comment, no PMID | Uniform 1e5 /M/s for **20 of 22** engagers in `eng_params_normalized.json` (RUN-VERIFY, this task: `Counter{1e5: 20, 1.157e5: 1, 9.954e5: 1}`); only **teclistamab** (1.1574e5) and **elranatamab** (9.9537e5) carry molecule-specific values. |
| `kon_TAA` | `1e5 /M/s` → 8.64 /nM/day (default) | /nM/day | `pd_model_config.py:34`; per-molecule override `run_tce_pd_reval.py:138` | **[CLASS-MEASURED / DERIVED — NOT molecule-measured]** ⚠ | RUN-VERIFY of the JSON: 15 of 22 sit at the generic 1e5; the 4 anti-CD20s at 4.3e5; teclistamab 1.2847e6; elranatamab 9.9537e5; catumaxomab 6.1e4. **Do not read the non-default values as measurements *of that molecule*.** The JSON's own `measured_vs_derived` string for mosunetuzumab says: *"CD20 kon=class-measured(rituximab SPR 4.3e5) … CD20 koff=DERIVED(kon*KD)"* — i.e. the on-rate is **transferred from rituximab**, and the off-rate is **computed**, not measured. The authority for any given molecule is that JSON string, which must be read before quoting a value as measured. |
| `koff_CD3` | `None` → `kon×KD_CD3` = **345.6 /day** at KD=40 nM | /day | `pd_model_config.py:35`; derived `wholebody_pd.py:234`; per-molecule override `run_tce_pd_reval.py:141` | **[DERIVED: koff = kon·KD]** unless measured | Range across the 22 engagers: **0.3 → 4320 /day** (cevostamab → REGN5459). This is the module's most discriminating input (EQ-13). |
| `koff_TAA` | `None` → `kon×KD_TAA` = **12.528 /day** at KD=1.45 nM | /day | `pd_model_config.py:36`; `wholebody_pd.py:235`; override `run_tce_pd_reval.py:139` | **[DERIVED: koff = kon·KD]** unless measured | Sets the trimer lifetime with `kint` (EQ-5, `m22`). |
| `kint_bridge` | **0.9** | /day | `pd_model_config.py:37` | **[UNSOURCED — TBD]** — comment says only *"trimer internalization (/day)"* | ⚠ **`eng_params_normalized.json` carries a per-molecule `kint_perday` (range 0.05–3.3 /day), and `run_tce_pd_reval.py:137–141` DOES NOT THREAD IT into `_kpm`.** Every molecule therefore runs the synapse at `kint = 0.9/day` regardless of its measured internalization. This matters more than it looks: §4.2 shows `kint` — not either KD — sets the engagement EC50. |
| `span_bridge_nm` | **12.5** | nm | `pd_model_config.py:38` | [ASSUMED] | Never overridden per construct (see §3.1). |
| `span_cis_nm` | 12.5 | nm | `pd_model_config.py:39` | [ASSUMED] | **Unused** (EQ-7). |
| `cis_avidity` | **0.0** | — | `pd_model_config.py:40` | [ASSUMED — inert at 0] | ⇒ `koff_CD3_eff = koff_CD3` exactly. The trispecific hook is present but not exercised in any validation run. |
| `k_hit_perday` | **12.0** | /day | `pd_model_config.py:41`, overwritten from `handoff/kinetic_calib.json` at `pd_model_config.py:59–60` | **[UNSOURCED — TBD]** (see §3.1) | Both sources carry 12.0, so the value is the same either way. Note the loader (`pd_model_config.py:55–56`) searches **`handoff/kinetic_calib.json`** relative to the module dir — the live copy is `rundir/handoff/kinetic_calib.json`. `params/kinetic_calib.json` is a **byte-identical duplicate** that the loader never reads; do not cite it as the live source. |
| `k_death` | **1.0** | dimensionless *in this path* (see EQ-15) | `handoff/kinetic_calib.json` (loaded at `pd_model_config.py:55–59`) | **[FITTED — with a contaminated secondary check]** | The JSON's own `k_death_provenance` string (quoted verbatim, `handoff/kinetic_calib.json`) claims a **PRIMARY** anchor (*"engaged-CTL serial ceiling k_hit*koff/(k_hit+koff)=11.6/day at KD_CD3=40nM; emergent engaged rate 2.7/day median inside Halle 2-16/day band, NO fitting"*) and a **SECONDARY** check (*"mosunetuzumab full PBPK-PD at k_death=1.0 gives IL6 609 (clin 570)"*). **On the 570 pg/mL secondary anchor, be precise about what the code actually says:** `wholebody_pd.py:23` *does* carry an attribution — `# mosunetuzumab peak IL-6 ~570 pg/mL (Hosseini 2020 Fig5A)`. It is therefore **not** true that the number is uncited in code. What *is* true, per the in-repo provenance audit (`docs/PROVENANCE_AND_VALIDATION.md` §1.1–1.3, 2026-07-13), is that **the attribution does not survive checking**: that audit records 570 as *"No source exists. The real reported population mean is **152**"*, and retains only two clinical IL-6 anchors — mosunetuzumab **152** pg/mL and teclistamab **21** pg/mL (the latter sourced to **PMID 38831634**). Both figures here are cited **to that audit document**, not asserted from memory. The primary anchor's *"Halle 2-16/day band"* names an author but **no PMID/DOI appears in any code or params file** → **[UNVERIFIED CITATION]**. **`k_death = 1.0` must be treated as [FITTED] until the primary anchor is sourced.** (The 11.6/day figure itself IS reproducible from this module — RUN-VERIFY, EQ-13 — so only the *biological band* it is compared against is unverified.) |
| `kdeg_CD3` | **0.0** | /day | `wholebody_pd.py:247` (`.get(..., 0.0)`, key absent) | [CODE-INTERNAL — dead] | EQ-8 never executes. |
| `R_CELL_UM` in `ageff_nM` call | 8.0 | µm | :131 (passes the module constant) | see §3.1 | Note: this is the **TARGET** cell radius, applied to the target's TAA surface density — correct usage. |

**NOT owned here (for the record):** `KD_CD3 = 40 nM` and `KD_TAA = 1.45 nM` defaults live at
`wholebody_pd.py:111`; the synapse reach `R_SYN_UM = 30 µm`, Treg constants `R_TREG_UM = 50`, `TREG_K = 0.25`
live at `wholebody_pd.py:36,71,72`. `CYTO_IL6_CLINICAL_ANCHOR_PGML = 570.0` (`wholebody_pd.py:23`) carries an
in-code attribution — *"(Hosseini 2020 Fig5A)"* — which the 2026-07-13 provenance audit
(`docs/PROVENANCE_AND_VALIDATION.md` §1.1) found **does not hold** ("No source exists"); it is *not* consumed
by this module in any case.

---

## 4. WHAT IS EMERGENT vs IMPOSED

This section is the reason the doc exists. The module's selling point is emergence; here is exactly where it
delivers and exactly where it stops. **Every claim below is either an analytic derivation performed in this
task or a RUN-VERIFY on the live module — with the configuration stated.**

**Standard RUN-VERIFY configuration used throughout §4** unless stated: 1 T cell (CD3 = 92,000 copies →
`RC = 2.1479 nM`), 4 target neighbours (CEACAM5 = 257,000 copies → `RTAA = 6.0 nM` each,
`taa_cap = 24.0 nM`), `kon = 8.64 /nM/day` both arms, `KD_CD3 = 40 nM` (`koff = 345.6/day`),
`KD_TAA = 1.45 nM` (`koff = 12.528/day`), `kint = 0.9/day`, `span = 12.5 nm` (`c_eff,trans = 57,553 nM`),
integrated to steady state.

---

### 4.1 ✅ EMERGENT — the QSS limit (the module's founding claim, and it holds)

The header claims (:4–7): *"The QSS limit EMERGES automatically when koff is fast vs the PD timescale; a slow
(high-affinity) arm produces a genuine lag — no QSS assumption is imposed anywhere."*

**This is TRUE, and it is a direct consequence of EQ-6.** The matrix exponential carries both eigenvalues.
The fast eigenvalue is `λ_fast ≈ −(rate_on + koff_CD3 + kf) ≈ −5×10⁵/day`; over `dt = 0.06 d`,
`exp(λ_fast·dt) = exp(−3×10⁴) → 0` — the fast mode relaxes fully within one step and the solution *becomes*
its quasi-steady value, with **no QSS assumption written anywhere.** The slow eigenvalue is
`λ_slow ≈ −(koff_CD3/ρ + kint) = −0.909/day` (τ ≈ 1.1 d) — a genuine, integrated lag that a QSS engine would
erase.

**RUN-VERIFY of the lag:** at `kint = 0` the slow eigenvalue collapses to `−koff_CD3/ρ = −0.0087/day`
(τ ≈ **115 days**). My steady-state probe (4.5-day integration) had *not converged* at `kint = 0` — which is
itself the proof: the module really is carrying a slow mode, and its relaxation time is set by
**internalization**, not by either affinity.

---

### 4.2 ✅ EMERGENT — bivalent avidity, quantified: a **380×** apparent-potency gain

The trimer's `m22 = −(koff_TAA + kint)` (EQ-5) — with `koff_CD3` deliberately absent — produces a genuinely
emergent avidity that is *nowhere written as a parameter*. Solving EQ-5 at steady state (derived this task):

```
ρ         ≡  kf / (koff_TAA + kint)                            [dimensionless bridging propensity]
B2_max    =  RC · ρ/(1+ρ)                                      [nM]
EC50_eng  =  ( koff_CD3/ρ + kint ) / ( kon_CD3 · (1 + 1/ρ) )   [nM]   →  ≈ (koff_CD3/ρ + kint)/kon_CD3  for ρ≫1

⇒   B2(C)  =  B2_max · C / ( C + EC50_eng )
```

**RUN-VERIFY (live module, integrated to convergence — the exact integrator makes `M` constant, so a long
step is exact; see the note below):**

| `kint` (/day) | ρ | `EC50_eng` predicted (nM) | measured `B2/B2max` at that C | slow-mode τ | verdict |
|---|---|---|---|---|---|
| 0.90 (**live**) | 37,031 | **0.105244** | **0.5000** | 1.1 d | formula **exact** |
| 2.00 | 34,228 | 0.232643 | 0.5000 | 0.50 d | exact |
| 10.0 | 22,073 | 1.159167 | 0.5000 | 0.10 d | exact |
| 0.10 | 39,377 | 0.012590 | 0.5000 | 9.2 d | exact |
| 0.00 | 39,692 | 0.001008 | 0.5000 | 115 d | exact |

> **Correction to a prior version of this doc.** The last two rows were previously reported as 0.312 and
> 0.038 ("not converged"). Those were artifacts of a **4.5-day** integration probe — at `kint ≤ 0.1` the slow
> mode has τ = 9.2–115 d, so 4.5 days is nowhere near steady state. Integrated to convergence the closed form
> is **exact at every `kint`**, as shown. The non-convergence is not a defect in the formula; it *is* the slow
> mode of §4.1, and it is real.

**The result:** at the live operating point the T cell reaches half-maximal engagement at
**`C = 0.1052 nM`, against a nominal `KD_CD3 = 40 nM`** — a **380-fold avidity enhancement**
(40 / 0.105244 = 380.1), computed from mechanism (geometry → `c_eff` → `ρ`), not fitted. **This is the model's
single best emergent result** and it is exactly the sort of thing an equilibrium engine cannot produce.

**And the sting in the tail:** because `ρ = 37,031 ≫ 1`, `EC50_eng ≈ kint/kon_CD3 = 0.9/8.64 = 0.104 nM`.
**The engagement potency of the synapse is set by INTERNALIZATION, not by either arm's affinity.** That is a
genuine, defensible mechanistic prediction — and it means the single unsourced `kint = 0.9/day` (§3.2) is the
most potency-determining number in the PD engine, while the *measured* per-molecule `kint` values sitting in
`eng_params_normalized.json` are **not being used** (§3.2). That is the highest-value fix in this subsystem.

---

### 4.3 ❌ **NOT EMERGENT — THE PROZONE / HOOK EFFECT IS STRUCTURALLY ABSENT**

**This must be stated plainly, because the model's own documentation asserts the opposite.** The QSS sibling
`wholebody_pd.ternary_equilibrium` genuinely has a hook — `wholebody_pd.py:87–91` explains it correctly:
*"free receptors R_A, R_B are solved from the QE quadratic BEFORE RC_AB, so the prozone/hook is emergent (at
high C both arms saturate as binary complexes → trimer collapses)."* And `pd_model_config.py:14–15` describes
the `'qss'` engine as *"Prozone emergent."*

**The `'kinetic'` engine — the CANONICAL one — does not have it.**

**(a) Analytic proof.** Solving EQ-5 at steady state (this task):
```
B2(rate_on)  =  rate_on · RC / [ rate_on·(1 + 1/ρ)  +  koff_CD3/ρ  +  kint ]
```
This is a **Michaelis (monotone-increasing, saturating) hyperbola in `rate_on`, hence in C.** It has **no
interior maximum and no descending limb** for any positive parameter values. As `C → ∞`,
`B2 → RC·ρ/(1+ρ)` — it saturates, it does not collapse.

**(b) RUN-VERIFY (live module, integrated to convergence, drug swept over 14 orders of magnitude):**

| C (nM) | 1e-4 | 1e-2 | 0.1 | 1 | 10 | 100 | 1e4 | 1e6 | **1e10** |
|---|---|---|---|---|---|---|---|---|---|
| `B2` (nM) | 0.002039 | 0.18637 | 1.04646 | 1.94328 | 2.12543 | 2.14554 | 2.14778 | 2.14780 | **2.14780** |
| `p_eng` | 0.0009 | 0.0868 | 0.4872 | 0.9048 | 0.9896 | 0.9989 | 1.0000 | 1.0000 | **1.0000** |

**Monotone increasing at every point. No hook. No down-turn. Not even at 10 µM (10¹⁰ nM).** The numerical
limit (2.147802) matches the analytic limit `RC·ρ/(1+ρ) = 2.147802` to 7 significant figures.

> **Correction to a prior version of this doc.** The three low-`C` entries were previously given as
> 0.00199 / 0.1829 / 1.0456 (and `p_eng` 0.085). Those came from a **4.5-day** integration that had not
> reached steady state, but were tabulated as if it had. The converged values are above. **The conclusion is
> unchanged and in fact strengthened** — the curve is monotone in both the converged and unconverged data.

**(c) WHY it is absent — the missing state.** The classical prozone requires that, at high drug, each CD3 is
occupied by *its own* drug **and each TAA is occupied by its own drug**, so no single molecule can bridge
both. That requires a **drug·TAA binary complex** that *sequesters TAA away from the bridge*. In this
two-state reduction that species **does not exist**:
- Free drug `C` enters the ODE **only** through `rate_on = kon_CD3·C` (:167). It never touches the TAA arm.
- The bridge-forming rate `kf = kon_TAA · c_eff,trans · alive_frac` (:168) is **completely independent of
  C**. It is a function of geometry and survival only.
- `RTAA_j` (:118) is set once from `wholebody_pd.py:169` and is **never reduced by drug occupancy** — not by
  this module, and not by the PK engine, which does model drug·antigen binding but never feeds the resulting
  occupancy back into `R_TAA`.

So the TAA arm is **never blocked by excess free drug**. The competition that *is* the hook is not in the
model.

**(d) What this means for the science.** The kinetic engine will predict **monotone-increasing efficacy with
dose, forever.** It cannot reproduce the bell-shaped dose–response that is a *defining, clinically observed*
feature of T-cell engagers (and the reason step-up dosing exists). Any dose-optimization, any "is there an
optimal dose" claim, any high-dose extrapolation made on the canonical `'kinetic'` engine is **structurally
incapable of finding a hook** — and its absence in a result is therefore **not evidence of its absence in
biology.** This is the first thing a reviewer will ask about, and the honest answer is: *the canonical engine
does not model it; the non-canonical QSS engine does.*

**(e) The fix (specified, not implemented).** Add a third state `B_T` = drug bound to TAA **only** (no CD3),
which (i) forms from free drug at `kon_TAA·C·(free TAA)`, (ii) **subtracts from the TAA available to `kf`**,
and (iii) can be re-captured by CD3 to form `B2`. That single state restores the competition and the hook
emerges automatically, without any imposed bell curve. It also fixes the CD3-release-from-trimer gap noted in
EQ-5(c).

---

### 4.4 ✅ EMERGENT — serial killing (the feedback loop is real and closed)

Serial killing is **not** parameterized. It emerges from a closed loop with **three** distinct channels by
which a target's death changes its killer's behaviour:

1. `alive_frac` (EQ-3, :164) shrinks `Tfree` → shrinks `kf` → shrinks `B2`;
2. `taa_cap` (EQ-10, :194) carries `surv_j` → the TAA ceiling falls as neighbours die;
3. `Ws` (EQ-11, :205–207) **renormalizes onto survivors** → the T cell's remaining output concentrates on
   the living rather than being wasted on the dead.

and the loop is closed by `wholebody_pd.py:468` feeding `exp(−kill_hazard)` straight back in as `surv_target`.

The **throughput** of that loop is set by EQ-13's race and, as §EQ-13's RUN-VERIFY table shows, spans a **43×
range across the 22 clinical engagers** driven by their *measured* CD3 off-rates. Cevostamab's 0.279/day
throughput ceiling — from `KD_CD3 = 0.033 nM` — is a real prediction of a real mechanism. **This is the second
of the module's two genuinely emergent, genuinely valuable results.**

*Caveat (from EQ-13):* the hit itself is a QSS race probability, not an integrated state. The *cycling* is
emergent; the *hit* is algebraic.

---

### 4.5 ⚠️ IMPOSED (and dominant) — the abundance response is a `min()` clamp, not mass action

The code comment at :119–123 admits the mechanism (quoted in full at EQ-10). Here is what it costs, measured.

Because `kf = 497,258/day` while `(koff_TAA + kint) = 13.43/day`, the bridging propensity is
**`ρ = 37,031`** — the second-arm reaction sits **4.6 orders of magnitude into saturation**. *Essentially
every armed CD3 is instantaneously bridged.* The ODE therefore degenerates:

```
B2  ≈  min(  RC · C/(C + 0.105 nM)  ,  taa_cap  )
```

**RUN-VERIFY (live module, saturating drug C = 100 nM, `RC = 2.1479 nM`):**

| TAA copies/cell | `taa_cap` (nM) | `B2` (nM) | `min(RC, taa_cap)` | `B2 / min` |
|---|---|---|---|---|
| 1,000 | 0.0934 | 0.0934 | 0.0934 | **1.0000** |
| 10,000 (BCMA-like) | 0.9339 | 0.9339 | 0.9339 | **1.0000** |
| 30,000 | 2.8016 | 2.1449 | 2.1479 | 0.9986 |
| 100,000 (CD20-like) | 9.3385 | 2.1454 | 2.1479 | 0.9989 |
| 257,000 (CEACAM5) | 24.0000 | 2.1455 | 2.1479 | 0.9989 |
| 1,000,000 | 93.3852 | 2.1456 | 2.1479 | 0.9989 |

**The "literal kinetic multivalent-Rhoden ODE" reproduces `min(R_CD3, TAA_cap)` to within 0.11% across three
orders of magnitude of antigen abundance.** The entire antigen-abundance response of the canonical PD engine
is carried by the hard clamp at lines 193–199 — **an imposed `min()`, not binding thermodynamics.**

**Corollary — the TAA arm's affinity is a near-null knob.** RUN-VERIFY (high-abundance target, C = 10 nM):

| `KD_TAA` (nM) | 0.01 | 0.1 | **1.45** | 10 | 100 | 1,000 | 10,000 |
|---|---|---|---|---|---|---|---|
| `B2` (nM) | 2.1257 | 2.1257 | 2.1254 | 2.1239 | 2.1076 | 1.9574 | 1.1430 |
| `p_eng` | 0.990 | 0.990 | 0.990 | 0.989 | 0.981 | 0.911 | 0.532 |

**A 10,000-fold change in TAA affinity (0.01 → 100 nM) moves the engaged fraction by 0.9%.** The TAA arm only
becomes rate-limiting when `KD_TAA` approaches `c_eff,trans = 57.6 µM` — i.e. outside any physical antibody
range. **For a tool whose stated purpose is an affinity + format design sweep, the TAA-affinity axis is,
at physiological antigen density, effectively inert.**

**Root cause (the two-basis mismatch).** The same 257,000 copies are represented **twice, incompatibly**:
- as a **capacity**: `257,000 × NM_PER_COPY = 6.0 nM` (the `taa_cap` basis, a 71 pL volume);
- as a **local tethered concentration**: `ageff_nM(257,000) = 63,676 nM` (the Rhoden hemisphere basis, a
  `(2/3)π(12.5 nm)³ = 4.1×10⁻⁶ µm³` volume).

The ratio is **10,600×**. Both are individually defensible physical quantities — but the model uses the huge
one for the **rate** and the small one for the **cap**, so the rate always wins and the cap always decides.
The `min()` clamp is not a bug; it is a **patch for a unit-basis mismatch**, and the code comment at :119–123
says so in as many words. A reviewer who spots this will ask why the bridge reaction is 4.6 logs from
equilibrium, and the answer must be this paragraph, not a hand-wave.

---

### 4.6 Summary ledger

| Quantity | Emergent? | Evidence |
|---|---|---|
| QSS limit / slow-arm lag | ✅ **EMERGENT** | §4.1 — from the matrix exponential's two eigenvalues; no QSS written |
| Bivalent avidity (380× potency gain) | ✅ **EMERGENT** | §4.2 — closed form derived + RUN-VERIFY exact to 4 dp |
| Trimer lifetime set by the slower arm | ✅ **EMERGENT** | EQ-5 `m22`; 1.79 h vs 4.17 min = 26× |
| Serial-killing cycle & its 43× koff spread | ✅ **EMERGENT** | §4.4 / EQ-13 — 3-channel closed feedback loop |
| Costim tumour-conditionality (via `p_eng`) | ✅ **EMERGENT** | §1.3 — `costim_induction.py:73–81` driven by B2/RC |
| Avidity magnitude (`c_eff`) | ✅ **EMERGENT from geometry** | EQ-1 — computed, not fitted |
| Cleft (relaxes to the bound complex) | ⚠️ **HALF** — the relaxation is emergent; the 13/40 nm bounds and the 0.6/0.4 ramp are imposed & unsourced | EQ-2 |
| **Prozone / hook** | ❌ **STRUCTURALLY ABSENT** | §4.3 — analytic proof + 14-decade sweep |
| **Antigen-abundance grading** | ❌ **IMPOSED** (`min()` clamp) | §4.5 — reproduces `min(RC, cap)` to 0.11% |
| **TAA-arm affinity response** | ❌ **~INERT** at physiological density | §4.5 — 0.9% over 4 logs of KD |
| The lethal hit (`k_hit`) | ❌ **IMPOSED** constant, **and QSS** (not a state) | §3.1, EQ-13 |
| `cis`-avidity (trispecific hook) | ❌ **IMPOSED** bare scalar on koff; geometry not wired; inert at 0.0 | EQ-7 |
| CD3 receptor turnover | ❌ **DEAD CODE** at the live config | EQ-8 |

---

## 5. KNOWN LIMITATIONS / OPEN QUESTIONS

*What a reviewer will attack, ranked by how much it hurts.*

### 5.1 The canonical engine cannot produce a hook (§4.3) — **the headline exposure**
A T-cell-engager PD model that is *structurally incapable* of a bell-shaped dose–response is a serious gap,
because the hook is (i) clinically real, (ii) the reason step-up dosing exists, and (iii) exactly the sort of
"hidden interaction" the module's own header (:9–12) says it was built to catch. The QSS engine that *does*
have it has been demoted to non-canonical. **Fix specified in §4.3(e): one extra state (`B_T` = drug·TAA
binary). Until then, no dose-optimum or high-dose claim may be made on this engine.**

### 5.2 The "window optimum" in the code comment is not in the equation
Line 212 asserts the engage/hit/detach race *"makes CD3 affinity a WINDOW optimum, not a monotone knob."*
**RUN-VERIFY (live module, C = 10 nM, high-abundance target):**

| `KD_CD3` (nM) | 0.1 | 1 | 10 | **40** | 100 | 1,000 |
|---|---|---|---|---|---|---|
| `koff_CD3` (/day) | 0.9 | 8.6 | 86.4 | 345.6 | 864 | 8,640 |
| `B2` (nM) | 2.1257 | 2.1257 | 2.1256 | 2.1254 | 2.1251 | 2.1200 |
| **`serial_rate` (/day)** | 0.798 | 4.971 | 10.427 | 11.476 | 11.710 | **11.828** |

`serial_rate` is **monotone increasing and saturating in `koff` — no interior optimum.** (`B2` itself is
essentially flat, so CD3 affinity acts *entirely* through the cycling term, not through the bound amount —
another consequence of §4.5's ρ ≫ 1.) The comment's own later text concedes this (:224–226), and further
claims an *"EMPIRICAL sweep (this build, k_hit=12/day): endpoint targets-killed peaks at koff ~ 346/day."*
**That peak is not reproducible from this module in isolation** (it would have to arise from the full
kill-feedback + PK loop), and **346/day is exactly the koff of mosunetuzumab, the calibration anchor** — a
coincidence that a reviewer will notice. The claim is marked **[UNVERIFIED IN-CODE CLAIM]** and should either
be reproduced with a committed sweep script or removed.

### 5.3 `K_HIT = 12/day` asserts a literature provenance and cites nothing (§3.1)
The comment says *"FIXED from serial-killing literature, NOT fitted"* — **with no PMID, DOI, author or year,
anywhere in the file.** `k_hit` is the ceiling on every kill rate in the model. Either source it or relabel it
`[ASSUMED]`. The related *"Halle 2-16/day band"* in `handoff/kinetic_calib.json` is likewise a bare surname
with no PMID anywhere in the repo.

### 5.4 `NM_PER_COPY` is a fitted scale wearing the name of a physical volume (§3.1)
`6.0/257,000` back-derives from a *calibration* constant (`Rcap_TAA = 6.0`), and the "71 pL synapse" it
implies is **~2×10⁵× the volume of a physical synaptic cleft** (RUN-VERIFY: 71,100 µm³ = a 25.7 µm-radius
sphere). It is internally consistent and it does the job it is actually there for (preserving absolute
abundance ordering across targets) — but it must not be presented as a measured synapse volume.

### 5.5 `k_death`'s clinical cross-check leans on a contaminated anchor
`handoff/kinetic_calib.json`'s `k_death_provenance` validates `k_death = 1.0` against *"IL6 609 (clin 570)"*.
The code **does** attribute that 570 — `wholebody_pd.py:23` reads *"(Hosseini 2020 Fig5A)"* — but the in-repo
provenance audit (`docs/PROVENANCE_AND_VALIDATION.md` §1.1–1.3, 2026-07-13) found the attribution does not
hold and records: *"No source exists. The real reported population mean is 152."* The two surviving clinical
IL-6 anchors in that audit are mosunetuzumab **152** pg/mL and teclistamab **21** pg/mL (**PMID 38831634**).
The *primary* anchor (the 11.6/day serial ceiling) **is** reproducible from this module (EQ-13 RUN-VERIFY) —
but the *"Halle 2-16/day"* biological band it is compared against carries **no PMID anywhere in the repo**.
`k_death = 1.0` should be tagged **[FITTED]**.

### 5.6 Three wiring gaps silently flatten the design sweep
- **`span_bridge_nm` is never overridden per construct** (`run_tce_pd_reval.py:137–141` threads only
  `kon`/`koff`), so all 22 engagers run at 12.5 nm and `feas = 0.9038`. **The format/geometry axis of the
  sweep is currently unexercised.** (Note this also means the ~2× `c_eff` difference between a 6.5 nm DART
  and a 13 nm BiTE — and the *hard zero* the DART would take from EQ-2 — never enters a validation run.)
- **`kint_bridge` is never overridden per construct**, despite `eng_params_normalized.json` carrying measured
  per-molecule `kint` (0.05 → 3.3 /day). Given §4.2 (`EC50_eng ≈ kint/kon_CD3`), this is the **single
  highest-impact unwired input in the subsystem** — it means every molecule is given the same synapse potency
  where the data says they differ by 66×.
- **`span_cis_nm` is accepted and never read** (EQ-7); the cis geometry→avidity computation the docstring
  promises is not implemented. `cis_avidity` is instead a bare scalar on `koff`, inert only because it is 0.0.

### 5.7 Dead code inside the module
`dwell_engaged` (:136, :139, :231) is **write-only** — the synapse-stability CRS mechanism it implements has
**zero readers** (§EQ-14). `self.surv` (:140) and `self._nTarget` (:141) are set and never used.
`bridged_total()` (:241) has no callers. EQ-8 (CD3 turnover) never executes. None of these are *wrong*; they
are unexercised surface area that a reviewer will read as live and ask about.

### 5.8 The two engines define `p_eng` differently
Kinetic: `p_eng = B2/RC` (:211, `wholebody_pd.py:485`). QSS: `p_eng = Cb/(Cb + RA)` (`wholebody_pd.py:333`)
for costim induction, and `Cb/RA` (`wholebody_pd.py:379`) for myeloid IL-6. These are **not the same
function** — `Cb/(Cb+RA)` saturates at 1 only as `Cb → ∞`, whereas `B2/RC` reaches 1 when the CD3 pool is
full. Any A/B comparison of `'kinetic'` vs `'qss'` on a costim or IL-6 endpoint is therefore comparing two
different definitions of engagement, not two kill laws.

### 5.9 The `max(disc, 0)` clamp is an undocumented approximation with a **named molecule inside its regime**
Complex eigenvalues are silently projected onto the degenerate real branch (:80). The condition is exact and
checkable (EQ-6): `disc < 0` ⟺ `koff_CD3_eff < kint` **and** `rate_on ≈ kf + koff_TAA`. At the live
`kint = 0.9/day` that first clause means `KD_CD3 < 0.104 nM` — and **cevostamab (`KD_CD3 = 0.033 nM`,
`koff_CD3 = 0.285/day`) is inside it.** **RUN-VERIFY (97,020-point resonance-resolved grid): 2,844 negative
discriminants (2.93%), worst `disc/tr² = −1.21×10⁻³`; at cevostamab's own arm kinetics, −9.9×10⁻⁴ at
C ≈ 70 nM.** The magnitude of the induced error is small (|Im λ|/|Re λ| ≲ 2% on a mode that relaxes fully
within one PD step), so this is **not** an emergency — but the earlier characterisation of the clamp as
"essentially never active" was an artifact of a grid too coarse to resolve the resonance, and any
sub-0.1 nM-CD3 construct the design sweep proposes lands in this regime **by construction**.

### 5.10 Operator splitting is first-order despite the "exact" integrator
EQ-6 is exact for **constant** `M` and `b` over `dt`. But `alive_frac` (:164), `taa_cap` (:194) and `RC`
(:161) are **frozen at the step's start** and updated between steps. Over `dt_PD = 1.44 h`, during which up
to ~2% of the trimer's own targets can die, the coupling between binding and killing is therefore
**first-order accurate, not exact.** No convergence study in `dt` exists. One should be run — it is cheap, and
it is the obvious question at a defence.

### 5.11 `R_CELL_UM = 8.0 µm` is uncited and quadratically load-bearing
`c_eff ∝ 1/r_cell²` (EQ-1), and **no source for the 8.0 µm value appears in any live code file** — the only
provenance string is *"Rhoden default"* (`:38`). *This doc deliberately does not offer a substitute radius:
supplying an uncited "true" B-cell radius would be the same provenance failure in the other direction.* What
can be stated without a source is the **sensitivity**: halving `r_cell` **quadruples** `c_eff` and hence `ρ`.
(In practice this is absorbed — ρ is already 37,031, so quadrupling it changes nothing; §4.5's degeneracy
makes the model insensitive to its own geometry. **That insensitivity is itself the finding**, and it is why
the uncited radius is not currently a numerical risk — only a provenance one.) **Action: source it or tag it
`[ASSUMED]`.**

### 5.12 Open questions for the committee
1. Will the `B_T` third state be added (§4.3(e)), restoring the hook to the canonical engine?
2. Can `k_hit` be sourced to a specific serial-killing measurement, or must it be relabelled `[ASSUMED]`?
3. Should `NM_PER_COPY` be re-derived from a *measured* synapse contact area × cleft, and the model
   re-calibrated on that basis — which would collapse `ρ` toward O(1) and let mass action (rather than a
   `min()` clamp) carry the abundance response?
4. Given §4.5, what *is* the model's claimed sensitivity to TAA affinity, and should the TAA arm be dropped
   from the design sweep until the basis mismatch is resolved?
5. Should `kint` (per molecule, already measured and sitting unused in `eng_params_normalized.json`) be
   threaded — given that §4.2 shows it, not affinity, sets synapse potency?
