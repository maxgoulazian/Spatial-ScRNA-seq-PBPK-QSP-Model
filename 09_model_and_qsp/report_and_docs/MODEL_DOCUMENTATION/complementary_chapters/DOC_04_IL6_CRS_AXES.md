# DOC_04 — IL-6 / CRS mechanism and the three-axis counter-screen

**Costim-engager counter-screen · QSP model documentation system**
Chapter 4 of the engine documentation. Companion to the binding/PK/PD chapters.

> **Scope — in-use final components only.** This chapter documents the live execution
> path only. The live path is the 12 engine files named in `SUBMISSION_MANIFEST.json →
> LIVE_FINAL_PATH.in_use_engine_files`, established from the runtime import trace of
> `run_tce_pd_reval.py`. For the IL-6/CRS mechanism the live files are
> **`myeloid_il6.py`** (the emitter + plasma ODE), **`coupled_percell_pd.py`** (plasma
> wiring), and **`wholebody_pd.py`** (per-organ instantiation + engagement driver). The
> retired fitted-scale files (`il6_pbpk.py`, `cytokine_pbpk.py`) are DEAD and are not
> documented here.
>
> **Ground-truth policy.** Every equation, constant, and parameter below was read from
> live source and is cited `file:line`. Where a number could not be found in source or
> the manifest, that is stated rather than filled from memory.

---

## 0. Three corrections to the incoming brief (live source wins)

The task brief and older docs carry three values that **do not match the live source**.
Per the ground-truth policy, the live runtime values are authoritative and are used
throughout this chapter. All three are flagged here up front.

| Quantity | Brief / older-doc value | **Live runtime value** | Where verified |
|---|---|---|---|
| Per-cell secretion `S_MAX_PG_PER_HR` | 0.0196 pg/hr/cell | **0.001331 pg/hr/cell** | `myeloid_il6.py:43` (`S_MAX_MOLEC_PER_S=10.6`) → `:106-111`; confirmed by importing the module |
| Myeloid↔engaged-T contact radius | 30 µm | **14.1 µm** | `myeloid_il6.py:101` (`R_CONTACT_UM=14.1`); imported + passed at `wholebody_pd.py:128-129` |
| Fitted `IL6_SCALE` path | "retired" | **deleted; now a hard `RuntimeError`** | `run_tce_pd_reval.py:203, 217-222` |

Detail:

1. **`S_MAX` = 0.001331 pg/hr/cell, not 0.0196.** The 0.0196 figure corresponds to the
   **retired** per-cell rate of 156 molec/s (the high-secretor tail). On 2026-07-13 the
   emitter rate was corrected down to the **measured population mean over secreting cells,
   10.6 molec/s** (`myeloid_il6.py:43-51`), which the unit-conversion `s_max_pg_per_hr()`
   turns into **0.001331 pg/hr/cell**. Importing the live module returns
   `S_MAX_PG_PER_HR = 0.001330689586870434`. **The source comment at `myeloid_il6.py:111`
   still reads "~0.0196 pg/hr/cell (verified this session)" — that comment is stale**; it
   predates the same file's own 156→10.6 correction two dozen lines above it. The code
   value (0.001331) is what runs.
2. **Contact radius = 14.1 µm, not 30 µm.** `MyeloidIL6.__init__` has a constructor
   *default* of `r_contact_um=30.0` (`myeloid_il6.py:124`), but the **only live call site
   overrides it** by importing and passing the module constant `R_CONTACT_UM = 14.1`
   (`wholebody_pd.py:128-129`). 14.1 µm is the physical two-sphere contact distance
   (macrophage radius 10.6 µm + T-cell radius 3.5 µm; `myeloid_il6.py:92-101`). The 30 µm
   value was the T-cell:target *synapse reach* `R_SYN_UM` (`wholebody_pd.py:36`), a
   reachability radius, not a contact distance; using it over-counted activation ~3×
   (`myeloid_il6.py:97-99`).
3. **The fitted `IL6_SCALE` path is deleted, not merely retired.** `run_tce_pd_reval.py`
   no longer falls back to it; an empty mechanistic IL-6 array raises
   `RuntimeError` (`run_tce_pd_reval.py:217-222`). See §A.1.

---

# PART A — Mechanistic myeloid IL-6 (the CRS magnitude engine)

`myeloid_il6.py` is the **final** IL-6 path. IL-6 is emitted per myeloid cell, spatially,
with **no fitted scale, no Emax, no EC50, no Hill term**. Saturation and the per-molecule
IL-6 spread emerge from the finite, spatially distributed myeloid pool plus a measured
first-order clearance. This section walks the mechanism, sources every parameter, and
documents why IL-6 is consumed by the counter-screen as a *between-molecule ratio* rather
than an absolute magnitude.

## A.1 Why the fitted `IL6_SCALE` was retired — and then deleted

The legacy engine computed `il6_pgml = engaged_dwell_rate * IL6_SCALE`, where `IL6_SCALE`
was one constant fit so that mosunetuzumab → 570 pg/mL (`myeloid_il6.py:3-6`). That was
broken three ways (`myeloid_il6.py:5-13`):

1. **Wrong emitter.** IL-6 is not made by T cells; it is myeloid-derived — "macrophages …
   the main overall source of IL-6" (Giavridis, *Nat Med* 2018, **PMID 29808005**) and
   "human monocytes were the major source of IL-1 and IL-6" (Norelli, *Nat Med* 2018,
   **PMID 29808007**). T cells are the *activating input*, not the emitter
   (`myeloid_il6.py:7-9`).
2. **No clearance.** The engine multiplied a production *rate* (pg/mL/hr) by a scalar and
   called the result a *concentration* (pg/mL). `IL6_SCALE` was silently standing in for
   the missing `1/kdeg` (`myeloid_il6.py:10-11`).
3. **Constructs washed out.** One global scale (or one global Emax) saturates every
   molecule to the same ceiling, erasing the between-construct differences a counter-screen
   exists to rank (`myeloid_il6.py:12-13`).

**Deletion, not just replacement.** As of 2026-07-13 the fitted legacy path is deleted from
the runner (`run_tce_pd_reval.py:203`). The old code contained a "safety" fallback —
`if il6.size==0: il6 = il6_legacy` — that on *any* failure of the mechanistic recorder
silently substituted the fitted 570-anchored number and emitted it under the same field
name (`run_tce_pd_reval.py:204-212`). That failure mode is exactly how a page number
(elranatamab "191", a dot-leader from an FDA Table of Figures) once circulated as a
clinical IL-6 value. The fallback is now a **hard error**: an empty mechanistic array
raises `RuntimeError` refusing to substitute a constant (`run_tce_pd_reval.py:217-222`).
The output field is stamped `il6_method="mechanistic_myeloid_percell"` as provenance
(`run_tce_pd_reval.py:235`).

> **Vestigial code note.** `wholebody_pd.py` still *defines* `CYTO_IL6_CLINICAL_ANCHOR_PGML
> = 570.0` (`:23`) and `cytokine_to_pgml()` (`:32-35`), but a full-tree grep shows **neither
> is called anywhere in the live path** — they are dead definitions left in the module, not
> live wiring. The `CYTO_HIER` cytokine-*rate* proxy (`wholebody_pd.py:22`) and its
> desensitization gate are a separate per-T cytokine tracker; the plasma IL-6 magnitude
> that reaches the QSP comes only from the myeloid path documented below.

## A.2 The mechanism (fully emergent, per cell)

Three coupled equations (`myeloid_il6.py:23-25`). Each myeloid agent *i* carries its own
activation state `a_i ∈ [0,1]`, driven by its own local contact with engaged T cells:

```
da_i/dt   = k_on · contact_i · (1 − a_i)  −  k_off · a_i          [structural saturation via (1−a_i)]
IL6_prod  = count_scale · Σ_i ( a_i · S_MAX · is_secretor_i )      [pg/hr]
dC/dt     = IL6_prod / V_ECF  −  kdeg · C                          [pg/mL/hr → pg/mL]
```

Implemented as:

- **Contact** (`myeloid_il6.py:176-181`): a `cKDTree` over T-cell coordinates; each myeloid
  cell sums the engaged-synapse level `p_eng` of the T cells within its own contact radius
  `r_contact` (= 14.1 µm live). `contact_i` is that local engaged-T load.
- **Activation ODE** (`myeloid_il6.py:184`): forward-Euler update
  `a += dt_hr·(K_ON·contact·(1−a) − K_OFF·a)`, then clipped to `[0,1]` (`:185`). The
  `(1−a)` term is the structural ceiling — no Emax/Hill/EC50.
- **Secretion** (`myeloid_il6.py:191`): `prod_pg_per_hr = Σ (a · S_MAX_PG_PER_HR · cs ·
  is_secretor)`. Only the intrinsic IL-6-secretor subset (`is_secretor`, §A.3) emits; each
  emitter secretes `a_i · S_MAX`, lifted by its own per-cell count scale `cs`.
- **Plasma ODE** (`PlasmaIL6.step`, `myeloid_il6.py:206-214`): exact exponential step of
  `dC/dt = R − kdeg·C` with `R = prod/V_ECF` — unconditionally stable, `C = Css + (C −
  Css)·exp(−kdeg·dt)`, `Css = R/kdeg` (`:213-214`). **This is the piece the engine never
  had**: it converts a production rate into a concentration (`myeloid_il6.py:200-202`).

**Emergent properties (nothing fitted)** (`myeloid_il6.py:27-34`):

- *Saturation* emerges from the finite myeloid pool + the per-cell `(1−a)` ceiling.
- The measured ~3.9% single-cell secretor fraction is imposed as a cell-intrinsic identity
  (§A.3), not fit.
- *Per-molecule differences emerge from anatomy*: a CD20 engager (mosunetuzumab) engages a
  large B-cell field in the spleen (16% myeloid, ~54,373 resident macrophages) → many
  myeloid activated → high IL-6; a BCMA engager (elranatamab) engages rare plasma cells
  whose mass sits in bone marrow (1.8% myeloid, ~3,000 cells) → few myeloid activated →
  lower IL-6 (`myeloid_il6.py:31-34`).

## A.3 Every parameter is literature-measured

All emitter/clearance constants live in `myeloid_il6.py:41-104`. None is fitted to model
output.

| Constant | Live value | Line | Source / meaning |
|---|---|---|---|
| `IL6_MW_DA` | 21 000 Da | `:41` | mature human IL-6 ~21 kDa |
| `S_MAX_MOLEC_PER_S` | **10.6 molec/s** | `:43-45` | mean per-cell secretion *over actively-secreting* monocytes, 12 h value (**PMID 20376398**) |
| `S_PEAK_MOLEC_PER_S` | 156 molec/s | `:52` | high-secretor tail (**PMID 37533643**); retained for a future heterogeneity refinement, **not** the emitter rate |
| `S_MAX_PG_PER_HR` | **0.001331 pg/hr/cell** | `:106-111` | `s_max_pg_per_hr(10.6)`; unit conversion `molec/s → pg/hr`. *(Source comment "~0.0196" is stale — see §0.)* |
| `SECRETOR_FRACTION` | 0.039 | `:53-61` | ~3.9% of stimulated monocytes are active IL-6 secretors even under maximal LPS stimulation → cell-intrinsic, not spatial (**PMID 37533643**) |
| `T_TO_MAX_MIN` | 150 min | `:62-63` | measured time to maximal secretion (**PMID 37533643**) → sets `k_on` |
| `K_ON_PER_HR` | **1.2 /hr** | `:112` | `3.0/(T_TO_MAX_MIN/60)`; reaches ~95% of max at the measured `t_to_max` |
| `K_OFF_PER_HR` | 0.10 /hr | `:102-104` | deactivation; slower than activation (IL-6 outlasts the stimulus), from the secretor-fraction decline (**PMID 37533643**) |
| `R_CONTACT_UM` | **14.1 µm** | `:92-101` | macrophage r 10.6 µm (**PMID 9400735**) + T-cell r 3.5 µm (**PMID 30571054**); CD40L–CD40 is membrane-bound → cells must touch |
| `KDEG_IL6_PER_HR` | 0.20 /hr | `:64` | IL-6 first-order elimination, 0.18–0.25 /hr, t½ ~2.8–3.8 h (**PMID 31268236**) |
| `V_PLASMA_ML` | 11 650 mL | `:65-91` | IL-6 ECF distribution volume = interstitial 8.55 L + plasma 3.10 L = 11.65 L; taken from the model's own PBPK volumes (§A.5) |

**Intrinsic secretor identity** (`myeloid_il6.py:141-146`): each myeloid cell is assigned
secretor/non-secretor **once**, deterministically (seeded RNG), with population fraction =
0.039. So a cell either is or is not an emitter; this is a measured constant restored, not a
knob (`:141-144`). Omitting it over-produced IL-6 ~23× (elranatamab 22,956 vs clinical ~191;
`:60`).

## A.4 `count_scale` and the blood-myeloid gate

The ABM is a spatial *subsample*; `count_scale` lifts each sampled myeloid agent to the
physiological cell count it represents. It is applied **per cell** (`self.cs`,
`myeloid_il6.py:140, 191`) because blood carries a per-lineage scale array — a population
mean would apply a lymphocyte-weighted factor to monocytes (`:137-139`).

**Tissue myeloid `count_scale`** is a drug-*independent* cellularity property, loaded from a
citation-gated organ census `handoff/organ_myeloid_counts.json`
(`wholebody_pd.py:42-70, 144`). Absent → every organ falls back to 1.0, and the module
emits a loud warning that IL-6 will be ~1e5× too low (`wholebody_pd.py:66-67`). Because
plasma IL-6 is **exactly linear** in this per-organ scale, a run at 1.0 can be rescaled
analytically post-hoc from the per-organ production trace with no re-run
(`wholebody_pd.py:40-41`; the unscaled per-organ trace is recorded at
`coupled_percell_pd.py:386`).

**Blood myeloid are GATED OFF** (`coupled_percell_pd.py:188`,
`self.blood_pd.myeloid.set_count_scale(0.0)`). The mechanistic justification
(`coupled_percell_pd.py:175-186`):

- IL-6 induction requires *sustained* CD40L–CD40 contact (Giavridis **PMID 29808005**). In
  flowing blood, leukocytes are ~50 µm apart and in motion → transient collisions, not
  synapses. Sustained contact occurs only where myeloid are **adherent or resident**.
- Per **PMID 3944542**, ~60% of blood monocytes are **marginating** (adherent to
  endothelium) and only ~40% freely circulating. The marginating/extravasated monocytes
  **are the myeloid already resident in the organ ABMs** (spleen 54,373 macrophages, marrow
  3,179). Counting blood myeloid too would **double-count** them.
- Worse, the blood ABM's coordinates are a synthetic 2D grid (~15 µm spacing) whose spatial
  contact is non-physical: it activates 98.6% of 2.0e9 monocytes and yields 61,874 pg/mL vs
  a severe-CRS ceiling of ~10–20k (`coupled_percell_pd.py:183-185`). Setting the blood scale
  to 0 removes a **geometry artifact**, not biology (`:185-186`).

## A.5 Plasma IL-6 ECF ODE — clearance and distribution volume

`PlasmaIL6` (`myeloid_il6.py:200-214`) integrates the measured first-order clearance:

- **`kdeg` = 0.20 /hr** (**PMID 31268236**, Chen 2019 fitted) → t½ = ln2/0.20 ≈ **3.5 h**.
- **`V_ECF` = 11.65 L** = interstitium 8.55 L + plasma 3.10 L (`myeloid_il6.py:65-67`). The
  denominator is the **ECF distribution space**, not plasma alone, because IL-6 is a 21 kDa
  cytokine that (a) extravasates more freely than a 150 kDa IgG, (b) has no FcRn recycling
  to retain it in circulation, and (c) is **produced in the interstitium** by tissue
  macrophages — it does not even start in plasma (`myeloid_il6.py:68-76`). The clinical assay
  reads the plasma concentration of a molecule distributed across that whole space, so the
  ODE denominator must be the distribution volume.
- **Implied clearance `CL` = kdeg · V_ECF = 0.20 /hr · 11.65 L · 24 = 55.9 L/day.** This
  clearance is **KEPT, not fitted** to peaks; the manifest notes the model's derived-CL
  (~76 L/day) is an *output cross-check*, not an input (`IL6_VALIDATION.clearance`).

**Documented open uncertainty (stated, not hidden)** (`myeloid_il6.py:82-91`): the correct
denominator for a tissue-produced 21 kDa cytokine measured in plasma is genuinely
unresolved. Plasma volume (3 L) is the conservative, transparent floor; a true Vd is
plausibly larger (interstitial distribution), which would *lower* predicted plasma IL-6. An
earlier revision briefly set 6.4 L citing an "IL-6 Vss" — that value is a ~150 kDa mAb
volume (almost certainly tocilizumab, anti-IL-6R), mis-attributed to the 21 kDa cytokine,
and was reverted (`:82-86`). This is a **named model uncertainty, not a knob to tune to the
clinical value** (`:87-91`). The mechanistically-complete version — routing IL-6 through the
antibody PBPK transport with its own reflection coefficient — is scoped but not yet done;
the lumped ECF compartment is a well-mixed approximation of it (`:77-81`).

## A.6 Why IL-6 is scored as a RATIO between molecules

**The absolute magnitude runs high; the mechanism-relevant signal is the between-molecule
ordering.** The counter-screen's job is to rank constructs, so IL-6 is consumed as a ratio
against clinical anchors that share a statistic and a step-up structure — never as an
absolute pg/mL value to be hit.

**Verified digitized clinical anchors** (`SUBMISSION_MANIFEST.json →
IL6_VALIDATION.anchors_VERIFIED_digitized`; anchor doc
`model/IL6_ANCHORS_VERIFIED_2026-07-13.md`):

| Molecule | Statistic | Peak IL-6 (pg/mL) | Regimen | Source |
|---|---|---|---|---|
| Mosunetuzumab | **mean** | 127.4 | 1/2/60 mg C1 step-up | Chen 2023 Fig 2 (n=212) |
| Glofitamab | **median** | 30.2 | 2.5/10/30 step-up | Djebli 2023 Table 2 |
| Talquetamab | **median** | 19.8 | 0.4 mg/kg QW | Willemin 2024 Fig 1a |
| Blinatumomab | median | 640 / 370 | — | Hosseini 2020 |

**Cleanest like-for-like test** (`IL6_VALIDATION.cleanest_test`): the **median pair
glofitamab 30.2 / talquetamab 19.8 = 1.53×** — same statistic (median) *and* same step-up
structure, so a model that reproduces this ratio is reproducing the biology, not the assay
convention. (Mosunetuzumab is a *mean*, not comparable like-for-like to the medians; the
blina/glofit median ratio is ~12–21× depending on which blina median.)

**Quarantined fabricated values — documented as NON-data** (`IL6_VALIDATION.QUARANTINED_
FABRICATED`). These are **not** clinical measurements and must never be used as anchors:

| Value | What it actually is |
|---|---|
| **570** | no source (was the wired mosun anchor for the retired `IL6_SCALE`) |
| **340** | a MagnetisMM figure-of-figures priming peak / page-number confusion |
| **230** | no source |
| **191** | a **page number** — an FDA BLA 761345 dot-leader |
| 366.88 | no source |

Teclistamab has **no digitized IL-6 curve** in the panel → dropped from the like-for-like
test; loose "21 mean / 288 individual" MajesTEC-1 values are not independently verified here
and are not cited as sourced anchors (`IL6_VALIDATION.teclistamab_note`).

## A.7 How the myeloid IL-6 is wired into plasma (`coupled_percell_pd.py`)

Per organ, each `OrganPD` instantiates its own `MyeloidIL6` emitter with the live 14.1 µm
contact radius (`wholebody_pd.py:128-129`) and, each PD step, drives it with the T cells'
engaged-synapse fraction `p_eng`:

- **Kinetic (canonical) path** (`wholebody_pd.py:485-486`): `p_eng = B2/RC` — each T cell's
  own engaged-synapse fraction from the exact `KineticSynapse` state.
- **QSS path** (`wholebody_pd.py:378-380`): `p_eng = Cb/(Cb+RA)` — engaged fraction from the
  Schropp equilibrium trimer.

The organ then reports `il6_prod_pg_hr` (`wholebody_pd.py:130, 380/486`). `CoupledPerCellPD`
sums every compartment's production and integrates the plasma ODE **once per PD step**
(`coupled_percell_pd.py:285-297`):

```python
_il6_prod = 0.0
for o in self.organs:                                            # coupled_percell_pd.py:286-288
    _cs = float(getattr(self.pd[o], "myeloid_count_scale", 1.0) or 1.0)   # tissue cellularity
    _il6_prod += _cs * float(getattr(self.pd[o], "il6_prod_pg_hr", 0.0))
# + blood_pd (scalar myeloid_count_scale) and heme_pd, same pattern     # :289-295
self._il6_prod_pg_hr = _il6_prod                                 # :296
self._plasma_il6.step(dt*pd_every, _il6_prod)                    # :297  -> PlasmaIL6 ECF ODE
```

`PlasmaIL6` is instantiated once at the top of the run loop
(`coupled_percell_pd.py:219-220`), and the physical plasma concentration is recorded every
step (`:380`, `il6_plasma_rec.append(self._plasma_il6.C)`), alongside the unscaled per-organ
production trace for post-hoc census rescaling (`:386`). The multiplier `_cs` here **must**
be the myeloid cell-count scale, **not** `graphs[o].count_scale` (which is antigen-derived
and differs per drug, 0.18×–5.34×) — using the antigen scale would scale the same monocytes
differently for mosun vs elran, a per-drug artifact that would corrupt the counter-screen
(`coupled_percell_pd.py:279-284`).

---

# PART B — The three-axis counter-screen and the 6-axis liability veto

The IL-6/CRS engine of Part A supplies the *magnitude* of the cytokine-release liability.
Part B is the **decision logic**: how each candidate costim receptor is scored on three
axes, and how a hard liability veto gates the panel **upstream of the QSP model** so that
QSP never re-decides the nomination.

## B.1 The three axes

Each axis answers one question about agonizing a candidate costim receptor. The scores are
in `COSTIM_FINAL_3AXIS_SCORE_v7.csv` (artifact `2490744f-…`); the scored panel is the
core-11 costim receptors.

**Axis 1 — Effector benefit (does agonism help the CD8 killers?).** Read from the
Schmidt et al. 2022 CRISPRa CD8 IFN-γ screen (GSE174255), column `E_schmidt_z`
(`A1_effector_axis_final.csv ← cd8_effector_scores.csv`, PIN v5f280455). This is the
**gain-of-function** anchor: its CRISPRa arm directly identifies costimulatory TNFRSF
members whose *overexpression* promotes IFN-γ — the agonism-direction evidence a CD4 CRISPRi
knockdown screen structurally cannot supply. A receptor with a low/negative `E_schmidt_z`
is effector-gated-out: it cannot supply signal-2 no matter how clean its toxicity axes.

**Axis 2 — Suppression liability (does its CD4 wiring feed the IL-10 / Treg-suppressive
program?).** Read from the CD4 Perturb-seq differential-expression matrix, column
`SUPP_agon` with BH-FDR `SUPP_q_BH` and a direction-concordance flag `SUPP_call`
(`A2_nomination_Stim48hr_3axis.csv`, within-donor Cliff's-delta agonism = −1×KD, Stouffer,
BH-FDR over core-11, 15-gene SUPP set). A Treg-fraction Mantel-Haenszel odds ratio across
4 donors (`tregfrac_OR_MH`) is carried as a corroborating readout.

**Axis 3 — CRS liability (does its CD4 wiring drive the storm cytokines?).** The storm set
is **TNF, IL-2, IFN-γ** — the CD4 helper-cell cytokines that seed the myeloid IL-6/IL-1
amplification of Part A. Read from the same CD4 DE matrix, column `CRS_agon` /
`CRS_q_BH` / `CRS_call`. The model's cytokine hierarchy `CYTO_HIER = {IL6:1.0, IFN:0.36,
TNF:0.31, IL2:0.18}` (`wholebody_pd.py:22`) encodes the same storm species as a per-T
production-rate proxy.

**Scope discipline (honest framing).** CD4 is the **counter-screen, not the effector axis**.
CD8 does the killing; the CD4 screen is used only for what it uniquely resolves — the
suppressive and cytokine-release programs. The effector axis comes from the separate CD8
screen. No claim is made that CD8 killing is measured in CD4 data.

## B.2 The 6-axis liability veto — gates UPSTREAM of QSP

The nomination rule (`FINAL_NOMINATION_v7.md`; manifest `NOMINATION.rule`):

> **The 6-axis liability veto is upstream and drive-independent; effector NEVER offsets a
> liability. QSP does not re-decide the nomination.**

An arm with **any** liability-up axis is eliminated. The six axes are:

1. **CRS** — storm cytokines up (TNF/IL-2/IFN-γ), Part A / Axis 3.
2. **SUPP** — IL-10 / Treg-suppressive program up, Axis 2.
3. **HELP-erosion** — erosion of the beneficial CD4-help / Tfh program.
4. **PROLIF** — proliferation-program up (uncontrolled expansion).
5. **EXH** — exhaustion program driven (liability) vs negated (favorable).
6. **DD_SUPP** — a data-driven suppression axis that **adds no independent gate call** —
   every arm it would flag is already SUPP-gated; both co-leads pass it (4-1BB q=0.79, CD27
   q=0.42) (`FINAL_NOMINATION_v7.md`).

The gate is materialized as the single canonical column `gate_status`, computed from the
per-axis veto calls and assertion-guarded (`COSTIM_FINAL_3AXIS_SCORE_v7.csv` header). **The
veto runs on the DE/GRN wiring alone; it does not consult QSP window or effector z.** This
is the architectural point: liability is decided *before* the QSP model sees the arm, so a
large effector benefit can never buy back a toxic arm.

## B.3 Result — CLEAN survivors and the CD28 proof-case

Verified directly from `COSTIM_FINAL_3AXIS_SCORE_v7.csv` (values below read in-kernel):

| Receptor | Effector `E_schmidt_z` | CRS call | SUPP call | Other flags | `gate_status` |
|---|---:|---|---|---|---|
| **4-1BB** (TNFRSF9) | **3.74** | ns | ns | HELP up (ns), EXH flat | **CLEAN** |
| **CD27** | **4.28** | ns | **down·conc** (favorable) | EXH-negating (favorable) | **CLEAN** |
| CD28 | 12.11 | **up·conc** | **up** | PROLIF up·conc | GATED[CRS,SUPP,PROLIF] |
| ICOS | −0.39 | down | down·conc | HELP down, PROLIF up | GATED[HELP,PROLIF] |
| DNAM1 (CD226) | 0.64 | ns | up·conc | EXH-driving | GATED[SUPP,EXH] |
| OX40 (TNFRSF4) | 2.07 | ns | up | EXH-driving | GATED[SUPP,EXH] |
| GITR (TNFRSF18) | 0.09 | ns | up | — | GATED[SUPP] |
| HVEM (TNFRSF14) | 0.02 | ns | up | EXH-driving | GATED[SUPP,EXH] |
| DR3 (TNFRSF25) | 1.58 | ns | up | EXH-driving | GATED[SUPP,EXH] |
| CD30 (TNFRSF8) | 3.22 | ns | down·conc | HELP down, PROLIF up·conc | GATED[HELP,PROLIF] |
| CD40 | 2.65 | ns | down | HELP down·conc, PROLIF up·conc; APC-side (not cis-costim) | GATED[HELP,PROLIF] |

**CLEAN set = {4-1BB, CD27}** — the only two arms passing all six axes and being
surface-bindable cis-costim T-cell receptors. The nomination is derived four independent
ways (effector-first genome-wide scan; CD8-vs-CD4 differential; network selectivity;
genome-scale GRN) and holds under both DE and GRN drive
(`FINAL_NOMINATION_v7.md`; manifest `NOMINATION.derived_n_independent_ways = 4`).

**CD28 is the proof that the counter-screen — not effector — drives the nomination.** CD28
has the **highest effector z in the panel (12.11)**, yet it is **gated anyway** on
[CRS, SUPP, PROLIF]. If effector could offset liability, CD28 would win; it does not. This
is the clearest demonstration that the veto is upstream and effector-independent
(`FINAL_NOMINATION_v7.md`; manifest `NOMINATION.gated.CD28`).

## B.4 QSP reads the survivors — it does not re-rank the panel

After the gate, the frozen QSP produces a therapeutic-window number per arm
(`A31b_QSP_rerun_GRN_vs_DE.csv`, artifact `88c264d7-…`). For the co-leads
(`SUBMISSION_MANIFEST.json → QSP_WINDOW_RESULTS.headline`):

- **4-1BB:** `qsp_window` +1.57 (DE) / +1.57 (GRN); TI 62.4 → 744.4; **cap = none** — the
  widest, CRS-coldest arm.
- **CD27:** `qsp_window` −2.37 (DE) → **+1.27 (GRN)**; the network (GRN) view reads CD27 as
  *safer* than per-gene DE; CRS-capped at monotherapy dose under both.

**Post-gate reading is mandatory** (`QSP_WINDOW_RESULTS.post_gate_warning`): raw GRN
`window_rank` puts **CD30 (0.746)** and **CD28 (0.724)** *above* the co-leads — but both are
hard-gated. Reading the window without the veto gate resurrects the "CD28 looks clean-ish"
artifact. **Gate first, then rank the CLEAN survivors.** The QSP window quantifies *how much*
window a surviving arm buys; it never resurrects a gated one.

## B.5 Residual liabilities carried forward (honest)

The CD4 counter-screen cannot see every liability (`FINAL_NOMINATION_v7.md`):

- **4-1BB:** urelumab-class agonists drive CD8/NK bystander proliferation and
  hepatotoxicity via **liver-myeloid 4-1BB**. The CD4 screen cannot see NK/myeloid — this is
  4-1BB's real residual liability, expected to **emerge from the PBPK-QSP liver
  compartment**, not a CD4-screen miss.
- **CD27:** clean on CD4, free of hepatotox baggage, broader healthy-tissue spatial coverage
  (9/11 organs vs 4-1BB 7/11) — the more spatially-verifiable co-lead.

---

## Appendix — the static-`R_costim` limitation (MUST disclose)

`R_costim` is set **once at initialization from resting copy numbers** and read unchanged at
every step (`SUBMISSION_MANIFEST.json → ARCHITECTURE_CONSTANTS.static_R_costim_LIMITATION`).
The engine captures costim *conditionality* through binding geometry (the cis-gate — costim
only co-engages a T cell already forming a CD3 trimer) but **not activation-induced receptor
upregulation**.

**Consequence, stated plainly:** a resting-copy ranking **under-rates 4-1BB, OX40, and
ICOS** — receptors that are low on resting T cells but strongly induced on activation — and
yields a **spurious "CD2 wins" ordering** (CD2 is constitutively high at rest). This is a
known limitation of the static-`R_costim` version, disclosed here as required. It does **not**
change the CLEAN nomination {4-1BB, CD27}, which is derived from the DE/GRN wiring and the
liability veto — not from resting-copy magnitude — but it does mean the resting-copy *drive
magnitude* for the induced receptors is a floor, not their activated value.

**Provisional items deliberately omitted.** Per the truth policy, claims that are
in-chat-only and not artifact-backed are excluded here: notably the
**OX40/GITR net-negative-kill in Treg-rich settings** claim, which is structurally
impossible under static `R_costim` and is labeled PROVISIONAL / omitted
(`SUBMISSION_MANIFEST.json → PROVISIONAL_DO_NOT_ASSERT`). Activation-induced `R_costim` is a
future upgrade, canonical only if wired with literature-sourced induction kinetics and a
version tag; until then the static version is canonical.

---

## Source ledger (files read for this chapter)

Live engine source at
`/media/balthasar-lab/RAID4/costim_engager_counterscreen/model/engine/`:

- `myeloid_il6.py` (216 lines) — emitter constants, `MyeloidIL6`, `PlasmaIL6` ECF ODE.
- `coupled_percell_pd.py` (395 lines) — blood-myeloid gate (`:188`), plasma IL-6 summation
  loop (`:285-297`), `PlasmaIL6` instantiation (`:219-220`), recorders (`:380, 386`).
- `wholebody_pd.py` (503 lines) — per-organ `MyeloidIL6` instantiation with 14.1 µm contact
  (`:128-129`), organ myeloid census loader (`:42-70, 144`), `p_eng` drivers (`:380, 486`),
  vestigial `CYTO_IL6_CLINICAL_ANCHOR_PGML`/`cytokine_to_pgml` (`:23, 32-35`).
- `run_tce_pd_reval.py` (harness) — IL6_SCALE deletion + hard-error fallback (`:203,
  217-222, 235`).
- `pd_model_config.py` — canonical engine selection (kinetic); k_hit=12, k_death=1.

Artifacts: `SUBMISSION_MANIFEST.json` (`70c9f973-…`), `COSTIM_FINAL_3AXIS_SCORE_v7.csv`
(`2490744f-…`), `FINAL_NOMINATION_v7.md` (`750853c9-…`), `A31b_QSP_rerun_GRN_vs_DE.csv`
(`88c264d7-…`), anchor doc `model/IL6_ANCHORS_VERIFIED_2026-07-13.md`.

All numeric values in this chapter were verified against live source or the manifest; none
was filled from memory. Where a source value is internally inconsistent (the stale
`S_MAX_PG_PER_HR` comment) or unresolved (the IL-6 distribution volume), that is flagged in
place rather than smoothed over.
