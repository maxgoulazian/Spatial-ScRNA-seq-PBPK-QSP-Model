---
title: "T8 — Mechanistic CRS / IL-6 (myeloid-derived, contact-gated, per-cell)"
subsystem: T8
model: costim_engager_counterscreen
source_of_truth: engine/myeloid_il6.py (216 lines, read in full this task)
consumers: engine/wholebody_pd.py, engine/coupled_percell_pd.py, engine/run_tce_pd_reval.py
date: 2026-07-13
generated_by: workflow-subagent T8
adversarially_verified: 2026-07-13 (second pass — 10 errors found and fixed; see APPENDIX B)
anchor_status: DISPUTED — see danger box item 3 before quoting any absolute pg/mL
---

# T8 — Mechanistic CRS / IL-6

> **Scope.** This document covers **`engine/myeloid_il6.py` only**, plus the exact lines in its live consumers
> (`wholebody_pd.py`, `coupled_percell_pd.py`, `run_tce_pd_reval.py`) that read from or write to it. The modules
> `il6_pbpk.py` and `cytokine_pbpk.py` are **NOT in the live import graph** and are therefore **not documented,
> not cited, and not treated as part of the model**, regardless of what their contents suggest.

> [!danger]+ READ THIS BEFORE ANY NUMBER IN THIS SUBSYSTEM IS QUOTED
> Three things in this subsystem are commonly mis-stated. All three are established below from live code:
> 1. **`kdeg = 0.20/hr` is `[FITTED]`, not measured.** Its in-code citation (PMID 31268236) is to a
>    *semi-mechanistic PK/PD **modelling*** paper (verified this task, Europe PMC: Chen, Kamperschroer, Wong,
>    Xuan, *Clin Transl Sci* 2019;12:600–608, "A Modeling Framework to Characterize Cytokine Release upon
>    T-Cell-Engaging Bispecific Antibody Treatment"; MeSH major topic **"Models, Biological"**). It reports **no
>    measured IL-6 clearance**. **Human IL-6 clearance appears to be UNMEASURED in the primary literature.**
>    This single constant, together with `V = 11.65 L`, sets the *entire absolute scale* of the IL-6 arm.
> 2. **The fitted `IL6_SCALE` path is DELETED**, not merely disabled (`run_tce_pd_reval.py:219–235`). An empty
>    mechanistic IL-6 array is now a **hard `RuntimeError`**, not a silent fallback to the fitted constant.
> 3. **Clinical IL-6 anchors: DISPUTED INSIDE THIS REPO — do not quote any single value as settled.**
>    The values **570, 340, 230, 366.88 have no source**; **191 is a page number** (a dot-leader from the Table
>    of Figures of FDA BLA 761345); **elranatamab has no clinical IL-6 value in existence.** That much is
>    consistent across every in-repo source. **But the *surviving* anchors are NOT settled**, and two in-repo
>    documents dated the same day disagree:
>    - The **live code** (`run_tce_pd_reval.py:62,74–77`, `il6_obs=152.0`) and `docs/PROVENANCE_AND_VALIDATION.md:46,84`
>      carry **mosunetuzumab = 152** and **teclistamab = 21** as "population means". **`152` carries NO citation
>      anywhere in the repository** — not a PMID, not a figure, not a table. Only **21** is verbatim-sourced
>      (PMID 38831634, MajesTEC-1).
>    - **`model/IL6_ANCHORS_VERIFIED_2026-07-13.md`** — the *digitized-database* anchor set (`params/mab_tce_pkpd.sqlite`)
>      — gives the mosunetuzumab **MEAN as 127.4 pg/mL** (Chen 2023 Fig2, n=212, 1/2/60 mg C1 step-up), **DROPS
>      teclistamab entirely** ("NO IL-6 curve in the digitized db"), and states in terms that the
>      **"21 mean / 288 individual" values "are NOT independently verified here … do not cite them as sourced
>      anchors."**
>
>    **⚠ THIS DOCUMENT PREVIOUSLY ASSERTED "the only valid ones are mosunetuzumab 152 and teclistamab 21" and
>    called 288 "real". That was an overclaim** — it propagated an uncited number (152) as sourced, and elevated
>    to "valid" two values that the repo's own verified-anchor file explicitly forbids citing. **Until the
>    conflict is adjudicated, quote the anchor with its provenance attached and its rival value beside it.**
>    Every absolute-scale statement in §4.3 inherits this uncertainty. Several **stale comments still inside the
>    live source** additionally repeat 570 and 191 (`myeloid_il6.py:6, 34`; `wholebody_pd.py:23`). **They must
>    not be propagated.** See §5.1.

---

## 1. PURPOSE & DATAFLOW POSITION

### 1.1 What this subsystem does

T8 answers one question: **given where a T-cell engager forms synapses in a spatially-resolved, per-cell whole-body
model, what plasma IL-6 concentration results?** It is the CRS (cytokine release syndrome) arm of the counterscreen,
and it is deliberately built so that CRS magnitude is a **consequence of anatomy** (where the drug's target antigen
lives, and how much myeloid tissue sits next to it) rather than a fitted per-molecule constant.

The design premise, stated in the module docstring (`myeloid_il6.py:1–37`) and **confirmed against the primary
literature this task**, is that **IL-6 in CRS is myeloid-derived, not T-cell-derived**:

- **Giavridis et al., *Nat Med* 2018;24:731–738 (PMID 29808005)** — VERIFIED. Abstract, verbatim: CRS severity is
  "mediated not by CAR T cell-derived cytokines, but by IL-6, IL-1 and nitric oxide (NO) produced by **recipient
  macrophages**". **⚠ Caveat that belongs in the body, not just the appendix: this is a MURINE CRS model — the
  "recipient macrophages" are MOUSE macrophages.** It establishes the myeloid-origin *mechanism*; it does not
  supply a human quantity.
- **Norelli et al., *Nat Med* 2018;24:739–748 (PMID 29808007)** — VERIFIED. Abstract, verbatim: "**Human monocytes
  were the major source of IL-1 and IL-6** during CRS." (Humanised-mouse model, but the monocytes are human — this
  is the stronger of the two for the human claim.)

T cells are therefore the **activating input**, not the emitter. The subsystem models a population of individual
myeloid agents, each of which:
1. senses **its own local contact** with engaged T cells (contact-gated — see EQ-3),
2. integrates that into **its own activation state** `a_i ∈ [0,1]` (EQ-4),
3. secretes IL-6 **only if it is one of the intrinsic ~3.9% secretor cells** (EQ-2, EQ-7),
4. contributes to a single well-mixed **plasma/ECF IL-6 compartment** with first-order elimination (EQ-9).

### 1.2 Position in the life of the molecule

```
                        [T1–T7: PK/transport + per-cell binding + synapse formation]
                                             │
                                             ▼
                        per-T-cell ENGAGED-SYNAPSE level  p_eng ∈ [0,1]
                        ├─ kinetic path : p_eng = B2 / RC          (wholebody_pd.py:486)
                        └─ QSS path     : p_eng = Cb / RA          (wholebody_pd.py:379)
                                             │
                                             ▼
  ┌──────────────────────────── T8 (this subsystem) ────────────────────────────┐
  │  MyeloidIL6.step()  — ONE INSTANCE PER ORGAN  (constructed wholebody_pd.py:129) │
  │                                                                             │
  │   myeloid agent coords (x,y) + cell-type labels  ─┐                         │
  │   engaged-T coords (T_x,T_y) + p_eng ────────────┤                         │
  │                                                   ▼                         │
  │        contact_i  = Σ p_eng over engaged T within R_CONTACT_UM  (EQ-2/EQ-3) │
  │        a_i        ← activation ODE, structural saturation (1−a_i)   (EQ-4)  │
  │        prod_organ = Σ_i a_i · S_MAX · cs_i · is_secretor_i     [pg/hr] (EQ-7)│
  │                                                   │                         │
  │   (returned to OrganPD.il6_prod_pg_hr — wholebody_pd.py:380 / :486)         │
  └───────────────────────────────────┬─────────────────────────────────────────┘
                                      ▼
       CoupledPerCellPD sums organs × myeloid_count_scale  (coupled_percell_pd.py:285–296)
       BLOOD myeloid are gated to zero:  set_count_scale(0.0)  (coupled_percell_pd.py:188)
                                      ▼
  ┌──────────────────────────── T8 (this subsystem) ────────────────────────────┐
  │  PlasmaIL6.step()  — ONE GLOBAL INSTANCE  (constructed coupled_percell_pd.py:220) │
  │        dC/dt = prod_total / V  −  kdeg · C      analytic step         (EQ-9) │
  └───────────────────────────────────┬─────────────────────────────────────────┘
                                      ▼
       il6_plasma_pgml[t]  →  il6_peak  →  CRS ranking / therapeutic-window readout
       (recorded coupled_percell_pd.py:380–386; consumed run_tce_pd_reval.py:215–235)
```

**Inputs consumed (none of them owned here):** myeloid agent coordinates `x,y` and cell-type labels `labs` from the
organ agent tables; engaged-T coordinates and `p_eng` from the synapse subsystem; the per-organ myeloid census
`myeloid_count_scale` (a **tissue** property, drug-independent, loaded from `handoff/organ_myeloid_counts.json`,
`wholebody_pd.py:42–70,144`).

**Outputs produced:** `il6_prod_pg_hr` per organ (pg/hr, **unscaled** by the myeloid census inside this module),
and the single systemic `PlasmaIL6.C` (pg/mL) — the CRS observable the counterscreen ranks on.

### 1.3 Why the fitted path was deleted

`run_tce_pd_reval.py:199–214` documents the removal in the code itself. The legacy engine computed
`il6_pgml = sys_cyto_rate['IL6'] × IL6_SCALE`, where `IL6_SCALE` was one constant fitted so that mosunetuzumab
returned 570 pg/mL. That is broken structurally (a production **rate** multiplied by a scalar and called a
**concentration** — the scalar was silently standing in for the missing `1/kdeg`), and broken empirically (570 has
**no source**). Worse, the old code fell back to the fitted value on any failure of the mechanistic recorder and
emitted it **under the same field name**. It is now a hard error:

```python
il6 = np.array(r.get('il6_plasma_pgml') or [])
if il6.size == 0:
    raise RuntimeError(... "Refusing to fall back to the retired fitted IL6_SCALE path" ...)
                                                             # run_tce_pd_reval.py:215–220
```
`il6_legacy_pgml` / `il6_legacy_peak` / `IL6_SCALE` are **deleted from the output record**
(`run_tce_pd_reval.py:232–233`), and the record carries `il6_method="mechanistic_myeloid_percell"`
(`run_tce_pd_reval.py:235`).

> **Residual:** the symbol `IL6_SCALE` still *exists* elsewhere (`pd_model_config.py:45,59,63`;
> `run_tce_pd_reval.py:34,41`) and is still loaded from a calibration JSON. It is **not** used by the IL-6 output
> path any more. It is a live symbol on a dead IL-6 path — a trip hazard, not a contamination of the current
> number. See §5.10.

---

## 2. GOVERNING EQUATIONS

Notation: `i` indexes **myeloid agents**; `j` indexes **T-cell agents**; `dt_d` = PD step in **days**;
`dt_h = 24·dt_d` = step in **hours**. All rates in this module are **per hour**.

---

### EQ-1 — Myeloid identity (substring token match) (`myeloid_il6.py:122, 128–130`)

```
MYELOID_TOKENS = ("macrophage", "monocyte", "myeloid", "kupffer", "microglia")

is_mye_i  =  OR over tokens t of  [ t ∈ lowercase(label_i) ]        (substring containment)
midx      =  { i : is_mye_i }
```

- **Biological meaning.** Selects, from the organ's real agent table, the cells that are competent to be the IL-6
  source. The module docstring (`:118–121`) records the actual labels observed in the tables: `'macrophage'`,
  `'classical monocyte'`, `'monocyte'`, `'intermediate monocyte'`, `'Myeloid'`, `'myeloid dendritic cell'`.
- **Mechanistic rationale.** The emitter population is not invented — it is **the real cells already present in the
  spatial agent tables** at their real coordinates. This is what makes the spatial contact gate (EQ-2) meaningful:
  the myeloid cells are *where the tissue actually puts them*.
- **Rejected alternative.** A bulk "myeloid density" scalar per organ (which is what the well-mixed QSP does, via
  `PB.myeloid` in `qsp_costim_window_v2.py`). Rejected because a density scalar cannot express *proximity to an
  engaged T cell*, and proximity is the entire mechanism (EQ-3).
- **Units.** Dimensionless boolean mask.
- **⚠ This is a substring match, not exact set-membership.** Any label containing `"myeloid"` — e.g. a myeloid
  *progenitor* in bone marrow, or `'myeloid dendritic cell'` (which the docstring confirms IS matched) — is
  admitted and then secretes at the **monocyte** rate. See §5.6.

---

### EQ-2 — Intrinsic IL-6-secretor identity (Bernoulli, drawn once per cell) (`myeloid_il6.py:145–146`)

```
seed_organ   = |hash(("il6_secretor", n_myeloid))| mod 2^32
is_secretor_i ~ Bernoulli( SECRETOR_FRACTION = 0.039 )      drawn ONCE at construction
```

- **Biological meaning.** **Most monocytes are simply not IL-6 secretors, even when maximally stimulated.** Only
  ~3.9% of them are. A cell either *is* or *is not* an emitter — this is a fixed cell property, assigned once, not
  a bulk multiplier applied to the population mean.
- **Mechanistic rationale — and an explicitly refuted earlier hypothesis.** The code (`:54–61`) records that this
  fraction was **previously assumed to EMERGE spatially** ("only the myeloid near engaged T cells activate"). **The
  source refutes that**: the 3.9% was measured by droplet microfluidics on **LPS-stimulated** monocytes — *every*
  cell was maximally stimulated, with **no spatial constraint** — and still only 3.9% secreted IL-6. It is
  therefore **cell-intrinsic heterogeneity**, NOT a spatial consequence. This is the single most important honesty
  correction in the subsystem: **the model's headline "emergence" claim does not extend to the secretor fraction.**
  It is an imposed constant (see §4).
- **Effect of omitting it** (recorded in-code, `:60`): over-produced IL-6 by ~23×.
- **Units.** Dimensionless boolean per cell; `SECRETOR_FRACTION` dimensionless.
- **⚠ The seed is NOT reproducible across processes — VERIFIED THIS TASK.** `hash()` on a tuple containing a `str`
  is salted by `PYTHONHASHSEED`, which Python randomises per interpreter by default. Three runs of
  `abs(hash(("il6_secretor", 5000))) % 2**32` in this session returned **1966803464, 2249461615, 1998730173**
  (with `PYTHONHASHSEED=0` fixed, they collapse to a single value). The in-code claim "deterministic seed …
  Deterministic per organ so runs are reproducible" (`:143–144`) is **false as written**. See §5.5.

---

### EQ-3 — Cell–cell CONTACT distance (`myeloid_il6.py:92–101`)

```
R_CONTACT_UM = r_macrophage + r_Tcell = 10.6 + 3.5 = 14.1 µm
```

- **Biological meaning.** **CD40L–CD40 is a membrane-bound ligand–receptor pair: the two cells must physically
  touch.** Two spheres touch when their centre-to-centre distance ≤ the sum of their radii. Myeloid activation is
  therefore gated on *contact*, not on being "nearby".
- **🔴 THE CD40L–CD40 ATTRIBUTION ITSELF HAS NO CITATION ANYWHERE IN `myeloid_il6.py` — VERIFIED THIS TASK.** The
  axis is named in the constructor docstring (`:125–127`) and in the module header (`:19–20`), and it is the *sole*
  mechanistic justification for choosing a **contact** radius over a **reach** radius — i.e. for the single most
  load-bearing geometric constant in the subsystem (§5.2, §5.3). It is asserted, not sourced. The nearby citation
  (PMID 29808005) is invoked for *proximity-dependence*, and the specific line the code quotes from it — *"IL-6
  induction and myeloid activation require proximity of CAR T cells and myeloid cells"* — **is not in that paper's
  abstract** (it would be in the body; unverified here). **A reviewer can ask "why CD40L?" and the code has no
  answer.** If the contact-gate mechanism is instead IL-1/inflammasome-, NO-, or soluble-mediator-driven, the
  radius is not `r_mac + r_T` at all and EQ-3 is the wrong functional form, not merely the wrong number.
- **Mechanistic rationale, and the alternative that was REJECTED.** This constant was previously set to
  `R_SYN = 30 µm` — the model's **T-cell : TARGET-cell synapse REACH**. A reach radius is a *reachability* radius
  (how far a motile T cell can find and engage a target), **not** a contact distance. Applying a reach radius to a
  contact-dependent interaction over-counts activation. The in-code correction note (`:97–99`) records the measured
  consequence on real spleen data: **producing fraction 1.92% at 30 µm → 0.60% at 13 µm** (a ~3.2× over-count).
- **Units.** µm. `r_macrophage = 10.6 µm` (= 21.2/2), `r_Tcell = 3.5 µm` (= 7/2).
- **⚠ TWO SOURCING PROBLEMS, both verified this task (see §3 rows and §5.3):**
  - **PMID 9400735** is **Krombach et al. 1997, *Environ Health Perspect* 105(Suppl 5):1261–1263, "Cell size of
    alveolar macrophages: an interspecies comparison."** The value **21.2 ± 0.3 µm is VERBATIM in the abstract** —
    but it is the diameter of a human **ALVEOLAR** macrophage (bronchoalveolar lavage, healthy non-smoking
    volunteers, Coulter volumetry). The code labels it "tissue macrophage". **The organs that dominate IL-6 in this
    model are spleen and bone marrow, not lung.** Alveolar macrophages are among the *largest* macrophages; using
    their radius inflates `R_CONTACT` and hence contact counts.
  - **PMID 30571054** is **Sauls, McCausland & Taylor, "Histology, T-Cell Lymphocyte", StatPearls (NCBI
    Bookshelf)** — a **tertiary study guide**, not a primary measurement. **Its abstract contains no cell
    diameter at all.** The 7 µm value is a commonplace textbook figure and is biologically uncontroversial, but the
    citation does not support it.
  - The 0.60% producing fraction was measured **at 13 µm**, not at the live **14.1 µm** (`:99`).

---

### EQ-4 — Per-cell local engaged-T contact load (`myeloid_il6.py:174–181`)

```
eng_j     = clip( T_engaged_j , 0 , 1 )                                        (:174)
tree      = cKDTree( {(T_x_j , T_y_j)} )                                       (:176, :160–162)
nb_i      = tree.query_ball_point( (x_i , y_i) , r = R_CONTACT_UM )            (:179)
contact_i = Σ_{j ∈ nb_i}  eng_j                                                (:180–181)
```

- **Biological meaning.** Each myeloid cell counts, within its own contact shell, the **engaged** T cells touching
  it — weighting each by *how engaged* it is (`p_eng ∈ [0,1]`, the fraction of that T cell's CD3 that is in a
  drug-bridged synapse). A T cell that is present but not engaged contributes nothing: it is not displaying CD40L.
  `contact_i` is therefore the **local engaged-T load on THIS myeloid cell** — an "engaged-T-equivalents" count.
- **Mechanistic rationale.** This is the entire per-molecule discrimination mechanism. A drug whose target antigen
  sits in a myeloid-rich compartment (spleen) puts engaged T cells inside many macrophages' contact shells; a drug
  whose target sits in a myeloid-poor compartment does not. **Nothing about the drug enters this equation except
  through `p_eng` and the T cells' positions.**
- **Rejected alternative.** An organ-level "engaged T-cell fraction × myeloid density" product. Rejected because it
  destroys the spatial covariance between where T cells engage and where myeloid cells sit — which is precisely
  the signal the counterscreen must resolve.
- **Units.** `contact_i` is dimensionless (a sum of engagement fractions, so ∈ [0, N_T_in_shell]); `R_CONTACT_UM`
  in µm; coordinates in µm.
- **⚠ `contact_i` is a SUM, not a mean or a fraction.** A myeloid cell touched by 20 fully-engaged T cells gets
  `contact = 20`, which drives EQ-5's activation **20× faster** than the rate at which `K_ON_PER_HR` was
  calibrated (which implicitly assumes `contact = 1`). Because activation saturates (EQ-5), this mostly does not
  change the *steady state* — but it does mean the 150-min time-to-max is only recovered for a cell with exactly
  one fully-engaged neighbour. See §5.4.
- **⚠ Coordinates are 2D and SUBSAMPLED.** The touch criterion `r_mac + r_T` is a **3D physical** criterion applied
  to a **2D projected, subsampled** coordinate frame. The number of T neighbours inside a 14.1 µm disc is therefore
  a property of the ABM's sampling density, not of physiology. This is the deepest structural criticism of the
  subsystem — see §5.2.

---

### EQ-5 — Per-cell myeloid activation ODE (`myeloid_il6.py:183–185`; k_on at `:112`)

```
da_i/dt  =  K_ON · contact_i · (1 − a_i)   −   K_OFF · a_i              [/hr]     (:184, docstring :23)

explicit-Euler realisation:
    a_i  ←  a_i + dt_h · ( K_ON·contact_i·(1 − a_i) − K_OFF·a_i )                 (:184)
    a_i  ←  clip( a_i , 0 , 1 )                                                    (:185)

K_ON  = 3.0 / (T_TO_MAX_MIN / 60)  =  3.0 / 2.5  =  1.2 /hr                        (:112)
K_OFF = 0.10 /hr                                                                   (:102)
```

- **Biological meaning.** A myeloid cell that is touched by engaged T cells switches on its IL-6 program (rate
  `K_ON·contact`), and relaxes back to rest when the stimulus goes away (rate `K_OFF`). `a_i` is the fraction of
  that cell's maximal IL-6 program that is running.
- **Mechanistic rationale — this is where SATURATION comes from.** The `(1 − a_i)` factor is a **structural**
  ceiling: a cell cannot be more than 100% activated. Combined with a **finite myeloid pool**, this produces a
  saturating dose-response **with no Emax, no Hill coefficient, and no EC50 anywhere in the model**
  (`:28`). The rejected alternative — an `Emax·C/(EC50+C)` cytokine term — was rejected because *a single global
  Emax makes every molecule saturate to the same ceiling, so construct differences wash out* (`:12–13`), which is
  fatal for a counterscreen whose entire job is to rank constructs.
- **Derivation of `K_ON` (verified this task).** The activation limb alone (`contact=1`, ignoring `K_OFF`) rises as
  `1 − exp(−K_ON·t)`. At `t = T_TO_MAX = 150 min = 2.5 h`: `1 − exp(−1.2 × 2.5) = 0.9502` — i.e. **95% of maximum at
  the measured time-to-max**, exactly as the code comment claims (`:112`). The `3.0` is the "three time-constants
  ≈ 95%" convention. **`K_ON` is therefore [DERIVED] from `T_TO_MAX_MIN`, not independently sourced.**
- **True steady state (computed this task).** Including `K_OFF`, at `contact = 1`:
  `a_ss = K_ON/(K_ON + K_OFF) = 1.2/1.3 = 0.923`. At `contact = 10`, `a_ss = 12/12.1 = 0.992`. **Activation is
  effectively a switch:** any myeloid cell with ≥1 engaged T-cell neighbour goes to ~0.92–1.0.
- **Units.** `a_i` dimensionless ∈[0,1]; `K_ON` has units **hr⁻¹ per unit engaged-T-contact**; `K_OFF` hr⁻¹;
  `dt_h` hr.
- **⚠ `K_OFF = 0.10/hr` is the weakest-sourced rate in the module.** The comment says it is "set from the measured
  decline of the secretor fraction with prolonged stimulation (PMID 37533643)". **No decline rate is derivable
  from that paper's abstract, and no derivation is shown in code.** Tag: **[ASSUMED]**, not [MEASURED]. It is,
  however, only weakly load-bearing: because activation saturates, `K_OFF` sets the *decay tail* after
  disengagement, not the peak. See §5.7.
- **⚠ Numerics.** Explicit Euler, then clipped. With `contact` large and `dt_h` not small, `K_ON·contact·dt_h` can
  exceed 1 and the un-clipped update would overshoot; the `clip` (`:185`) makes this **stable but rate-inexact** —
  the cell simply snaps to `a=1` in one step. Because the physiological answer in that regime is also ≈1, this is
  benign for the *peak*, but the model does **not** faithfully resolve the *rise time* of a heavily-contacted cell.
  Contrast EQ-9, which is solved **analytically** and is unconditionally exact.

---

### EQ-6 — Per-cell maximal secretion rate: molecules/s → pg/hr (`myeloid_il6.py:106–109, 111`)

```
S_MAX_PG_PER_HR  =  S_MAX_MOLEC_PER_S · (IL6_MW_DA / N_A) · 1e12 · 3600

                 =  10.6 · (21000 / 6.02214076e23) · 1e12 · 3600
                 =  1.3307e-3  pg/hr/cell            (recomputed this task)
```

- **Biological meaning.** Pure unit conversion: molecules of a 21 kDa protein leaving one cell per second →
  picograms per hour. `IL6_MW_DA = 21000` (`:41`), `N_A = 6.02214076e23` (`:42`).
- **Mechanistic rationale.** The per-cell secretion rate is the **only** quantity that converts "how many myeloid
  cells are activated" (which the ABM computes) into "how much protein enters the plasma" (which the ODE needs).
  There is no other scale factor anywhere in the chain — this is what replaced `IL6_SCALE`.
- **Units.** `S_MAX_MOLEC_PER_S` molecules·s⁻¹·cell⁻¹; `IL6_MW_DA` g/mol; `N_A` mol⁻¹; result pg·hr⁻¹·cell⁻¹.
- **⚠ THE IN-CODE COMMENT AT `:111` IS STALE — VERIFIED BY EXECUTION.** It reads
  `S_MAX_PG_PER_HR = s_max_pg_per_hr()  # ~0.0196 pg/hr/cell (verified this session)`. Importing the live module
  and printing the value returns **`0.0013306895868704338 pg/hr/cell`**. The value **0.0196 pg/hr/cell is what you
  get at `S_MAX = 156 molec/s`** (I recomputed: `s_max_pg_per_hr(156) = 0.019584`) — i.e. the comment is a
  left-over from the **retired peak-tail value**. **The live number is 14.7× smaller than the comment says.** The
  same stale `0.0196 pg/hr/cell` also appears in the consumer, `wholebody_pd.py:124`. The **code is correct; the
  comments are wrong.** See §5.1.
- **The 156 → 10.6 correction (in-code, `:46–51`).** `S_MAX` was previously 156 molec/s — the **peak /
  high-secretor-tail** rate. That is the rate of the top few percent of cells; the *same* paper measures only
  ~3.9% of cells as active secretors at peak frequency. Applying a top-tail rate to *every* activated cell
  over-states secretion by **156/10.6 = 14.7×**. The defensible per-activated-cell rate is the measured **mean over
  secreting cells (10.6)**. `S_PEAK_MOLEC_PER_S = 156.0` is **retained but unused** (`:52`) as the seed of a future
  per-cell heterogeneity refinement — explicitly **"NOT a scale knob"**.

---

### EQ-7 — Organ IL-6 production (`myeloid_il6.py:191`)

```
prod_organ  =  Σ_i  a_i · S_MAX_PG_PER_HR · cs_i · is_secretor_i           [pg/hr]
```

- **Biological meaning.** Every **activated** (`a_i`) **secretor-competent** (`is_secretor_i`) myeloid cell emits at
  `a_i` × its maximal rate. `cs_i` lifts each sampled agent to the number of physiological cells it stands for.
- **Mechanistic rationale for `cs_i` being PER-CELL, not a population mean** (`:137–139, :148–157`): blood carries a
  **per-lineage** count-scale array. Applying a population-mean scale would apply a **lymphocyte-weighted** factor
  to **monocytes** — the wrong cell type. `set_count_scale()` (`:148–157`) accepts a scalar, a myeloid-length array,
  or a full agent-table-length array (from which it selects `self.midx` — i.e. *this cell type's own* scales).
- **Composition note.** The effective producing fraction of the myeloid pool is
  `frac_activated × 0.039` — the secretor gate multiplies **on top of** the spatial activation gate. Both must fire.
- **Units.** `a_i` dimensionless; `S_MAX_PG_PER_HR` pg·hr⁻¹·cell⁻¹; `cs_i` physiological-cells per agent;
  `is_secretor_i` ∈ {0,1} → `prod_organ` in pg/hr.
- **⚠ Latent double-count hazard.** `MyeloidIL6.__init__` takes `count_scale` (default `1.0`, `:124`) **and** the
  consumer multiplies by `pd[o].myeloid_count_scale` again (`coupled_percell_pd.py:287–288`). At the live operating
  point this is safe — `OrganPD` constructs `MyeloidIL6` **without** passing `count_scale` (`wholebody_pd.py:129`),
  so `cs = 1.0` and the census is applied exactly once, downstream. But **passing `count_scale` at construction
  would silently square the census**. See §5.8.

---

### EQ-8 — Null-input decay branch (`myeloid_il6.py:171–173`)

```
if  n_myeloid == 0  or  len(T_x) == 0:
        a_i  ←  a_i · exp( −K_OFF · dt_h )        (exact decay)
        return 0.0                                 (production reported as ZERO)
```

- **Biological meaning.** With no myeloid cells, or no T cells at all in the organ, there is nothing to activate;
  existing activation relaxes.
- **⚠ Inconsistency (honest note).** The activation state is decayed **exactly** here (`exp(−K_OFF·dt_h)`) but by
  **explicit Euler** in the main branch (EQ-5). More importantly, this branch **returns production `0.0` even if
  `a_i > 0`** — i.e. if an organ's T-cell population were ever to empty while myeloid cells were still activated,
  their secretion would be **abruptly zeroed** rather than decaying with `a_i`. In the live runs organ T-cell
  populations do not empty, so this branch is not exercised on the production path; it is documented because a
  reviewer will find it.

---

### EQ-9 — Plasma / ECF IL-6 mass balance (`myeloid_il6.py:200–215`)

```
ODE                   :  dC/dt  =  R  −  kdeg · C ,        R = prod_total / V     [pg/mL/hr]

analytic step (:212–214, exact for piecewise-constant R over dt):
    C_ss  =  R / kdeg
    C(t+dt)  =  C_ss  +  ( C(t) − C_ss ) · exp( −kdeg · dt_h )

degenerate branch (:210–211):  if kdeg ≤ 0 :  C ← C + R · dt_h
```

- **Biological meaning.** **This is the piece the engine never had.** It converts a *production rate* (pg/hr) into a
  *concentration* (pg/mL). Myeloid cells in the interstitium secrete IL-6 into an extracellular fluid space of
  volume `V`; IL-6 is eliminated from that space first-order.
- **Mechanistic rationale for the ANALYTIC solve.** `dC/dt = R − kC` with piecewise-constant `R` has an exact
  closed-form solution, so the step is **unconditionally stable and exact** regardless of `dt` (`:207`). An explicit
  Euler step on a stiff first-order decay would impose a `dt < 1/kdeg` constraint on the whole whole-body PD loop.
  This is the correct engineering choice and is *not* an approximation.
- **Mechanistic rationale for `V = 11,650 mL = the ECF distribution space, NOT plasma** (`:65–91`). Three physical
  arguments, all in-code:
  1. **IL-6 is 21 kDa** vs a 150 kDa IgG → it extravasates **more** freely, not less;
  2. **it has no FcRn recycling** to retain it in the circulation;
  3. **it is PRODUCED IN THE INTERSTITIUM** (by tissue macrophages) — it does not even *start* in plasma.

  The clinical assay reads the **plasma** concentration of a molecule distributed across that **whole** space, so the
  ODE denominator must be the **distribution volume**. Using plasma volume (3.1 L) over-states IL-6 by
  **11.65/3.1 = 3.76×** (the code's own factor, `:76`).
- **`V` is taken from the MODEL'S OWN PBPK volumes, not from an external literature number** (`:66–67`):
  `sum(PB.Vis) + V_PLASMA` from `qsp_costim_window_v2._PBPKArrays`. **Recomputed this task by importing the live
  module: `PB.Vis.sum() = 8.5508 L`, `V_PLASMA = 3.1 L`, total = **11.6508 L** — the hard-coded `11650.0 mL`
  matches to 0.007%.** ✅ **[DERIVED — VERIFIED against the live PBPK arrays]**
- **Units.** `C` pg/mL; `R` pg·mL⁻¹·hr⁻¹; `prod_total` pg/hr; `V` mL; `kdeg` hr⁻¹; `dt_h` hr.
- **🔴 `kdeg = 0.20/hr` IS THE SINGLE LARGEST ASSUMPTION IN THE ENTIRE IL-6 ARM.** See EQ-9a.

---

### EQ-9a — The IL-6 clearance constant (the load-bearing assumption) (`myeloid_il6.py:64`)

```
KDEG_IL6_PER_HR = 0.20        # in-code comment: "IL-6 first-order elimination;
                              #  0.18-0.25 /hr (PMID 31268236) -> t1/2 ~2.8-3.8 h"
```

**Internal arithmetic (recomputed this task — the code's own comment is self-consistent):**
- `t½ = ln2 / 0.20 = 3.465 h` ✓ (the band: `ln2/0.25 = 2.77 h`, `ln2/0.18 = 3.85 h` → "2.8–3.8 h" ✓)
- **Implied systemic clearance: `CL = kdeg · V = 0.20 /hr × 11.65 L = 2.33 L/hr = 55.9 L/day`.**

**PROVENANCE — the citation does NOT support the number. VERIFIED THIS TASK against Europe PMC:**

| Claim | In-code | Reality |
|---|---|---|
| PMID 31268236 | cited as the source of "0.18–0.25 /hr" | **Chen X, Kamperschroer C, Wong G, Xuan D. *Clin Transl Sci* 2019;12(6):600–608.** "A Modeling Framework to Characterize Cytokine Release upon T-Cell-Engaging Bispecific Antibody Treatment: **Methodology and Opportunities**." Abstract: *"A 'fit-for-purpose' **semimechanistic** pharmacokinetic/pharmacodynamic **model** was developed…"*. MeSH **major topic: "Models, Biological"**. It **reports no measured IL-6 clearance.** |

**Therefore: `kdeg` is `[FITTED]` — a model parameter borrowed from another model, not a measurement.** Tagging it
`[MEASURED]` (as the `PlasmaIL6` docstring at `:201` still does — *"Systemic plasma IL-6 with **MEASURED**
first-order clearance"*) is **incorrect and must be corrected in the source.**

**Human IL-6 clearance appears to be UNMEASURED in the primary literature.** No radiolabelled-IL-6 or exogenous-
IL-6 infusion disposition study is cited anywhere in this module, and none was found in the provenance audit. This
is a genuine, stated gap — not a gap that is being papered over.

**Why this matters quantitatively.** Plasma IL-6 at steady state is
```
C_ss  =  prod_total / (V · kdeg)  =  prod_total / (11650 mL × 0.20 /hr)  =  prod_total / 2330
```
so **`C_ss` is exactly inversely proportional to `kdeg`** and exactly inversely proportional to `V`.
`V·kdeg = CL` is a **single lumped constant** (55.9 L/day) that the model's *entire absolute* IL-6 magnitude hangs
on, and **neither factor is measured**: `V` is [DERIVED] from the model's own PBPK volumes (defensible, but it is
the model's assumption about IL-6's distribution, not a measured IL-6 Vd), and `kdeg` is [FITTED] to another
model. A 2× error in `CL` is a 2× error in every predicted IL-6 concentration and could invert a construct ranking
against a clinical anchor. **The per-molecule ORDERING, however, is invariant to `CL`** — see §4.3.

**Sourcing history to NOT repeat (in-code, `:82–86`):** `V` was briefly set to 6.4 L citing an "IL-6 Vss". That
value (central 3.5 L + peripheral 2.9 L) is a textbook **~150 kDa monoclonal antibody** volume — almost certainly
**tocilizumab (anti-IL-6R)**, not IL-6 (21 kDa). Mis-attributing a mAb's Vd to the cytokine is a sourcing error,
and *it happened to shrink the residual 2.1× — exactly when to be most suspicious of a number.* It was reverted.

---

### EQ-10 — Whole-body production summation + the blood gate (consumer: `coupled_percell_pd.py:285–297`)

```
prod_total  =  Σ_organs  myeloid_count_scale[o] · il6_prod_pg_hr[o]                 (:286–288)
            +  myeloid_count_scale[blood] · il6_prod_pg_hr[blood]                   (:290–292)
            +  myeloid_count_scale[heme]  · il6_prod_pg_hr[heme]                    (:294–295)

PlasmaIL6.step( dt·pd_every , prod_total )                                          (:297)

BLOOD GATE:   blood_pd.myeloid.set_count_scale(0.0)     →  blood contributes EXACTLY 0   (:188)
```

- **Biological meaning of the blood gate.** IL-6 induction requires **sustained** CD40L–CD40 contact. In *flowing*
  blood, leukocytes are far apart and in motion — contacts are transient collisions, not synapses. Sustained
  contact occurs where myeloid cells are **adherent or resident**. The in-code argument (`:175–186`) is that ~60% of
  blood monocytes are **marginating** (adherent to endothelium) and those cells **are already counted** as the
  resident myeloid in the organ ABMs — so adding the flowing-blood pool would **double-count** them.
- **Why `myeloid_count_scale` is a TISSUE property and NOT `graphs[o].count_scale`** (`wholebody_pd.py:131–144`;
  `coupled_percell_pd.py:279–284`): `graphs[o].count_scale` is **antigen-derived** and was **measured to be
  drug-dependent (0.18×–5.34× between MS4A1 and TNFRSF17 in the same organ)**. Using it would scale the *identical
  monocyte population* differently for mosunetuzumab vs elranatamab — **a per-drug artifact that would corrupt the
  counterscreen**. The myeloid census must be drug-independent.
- **🔴 THE BLOOD GATE'S STATED NUMERICAL JUSTIFICATION IS STALE — VERIFIED THIS TASK.** The comment
  (`coupled_percell_pd.py:183–185`) justifies the gate with: *"it activates 98.6% of 2.0e9 monocytes and yields
  **61,874 pg/mL** (measured this session) vs a severe-CRS ceiling of ~10–20k."* **That number belongs to the OLD
  parameterisation.** Reconstructing it with the retired parameters
  (`S=156 molec/s → 0.0196 pg/hr/cell`, `V = 3.1 L`, **no secretor gate**):
  `2.0e9 × 0.986 × 0.0196 / (3100 × 0.20) = 62,340 pg/mL` — **reproduces the quoted 61,874 to within 1%**, confirming
  the reconstruction. Recomputing the **same blood configuration at the LIVE parameters**
  (`S_MAX = 1.3307e-3 pg/hr/cell`, secretor gate 0.039, `V = 11,650 mL`):
  `2.0e9 × 0.986 × 0.039 × 1.3307e-3 / (11650 × 0.20) = **43.9 pg/mL**`.
  **The artifact the gate was justified by is ~1,400× smaller than stated and is no longer an artifact-scale
  number.** The *mechanistic* argument for the gate (margination double-count; the blood ABM's synthetic 2D grid
  is non-physical) is **independent and may well still stand** — but the **quantitative** argument in the comment
  is dead and must not be quoted. See §5.9.

---

### EQ-11 — Diagnostic readouts (`myeloid_il6.py:194–197`)

```
n_myeloid       =  |midx|
frac_activated  =  mean_i [ a_i > 0.05 ]        (threshold 0.05 is a REPORTING threshold only)
mean_activation =  mean_i [ a_i ]
```
- Pure readouts; **do not feed back** into any state. The `0.05` threshold is a reporting convention and appears
  nowhere in the dynamics. Note `frac_activated` is computed over **all** myeloid, **not** gated on
  `is_secretor` — so it is **not** the producing fraction (the producing fraction is `≈ frac_activated × 0.039`).

---

## 3. PARAMETERS OWNED

Every module-level constant in `myeloid_il6.py`. **Provenance tags are per rule 4 of the doc brief. Every PMID
below was resolved against Europe PMC in this task; where the paper exists but does not (from its abstract)
support the specific number, that is stated explicitly.**

### 3.1 Secretion & cell biology

| Symbol | Value | Units | Tag | Source (verified this task) | Mechanistic rationale / why this value |
|---|---|---|---|---|---|
| `IL6_MW_DA` (`:41`) | 21000 | Da | **[TEXTBOOK — NO SOURCE IN CODE]** (*not* `[MEASURED]`: an earlier version of this table tagged it so, but the code carries **no PMID and no reference** for it) | mature human IL-6 ≈21 kDa. **No citation anywhere in the module.** | Sets molecules→mass conversion (EQ-6). Uncontroversial and a reviewer will not attack it — but by the tagging rules of this repo, an uncited number is not `[MEASURED]`. |
| `NA` (`:42`) | 6.02214076e23 | mol⁻¹ | **[DEFINED: SI exact]** (*not* `[MEASURED]` — Avogadro has been an **exactly defined** constant since the 2019 SI redefinition; it is no longer a measured quantity) | SI 2019 defining constant. | — |
| `S_MAX_MOLEC_PER_S` (`:43`) | **10.6** | molec·s⁻¹·cell⁻¹ | **[MEASURED: PMID 20376398]** — *paper VERIFIED; the specific value NOT confirmable from the abstract* → **[UNVERIFIED CITATION at the value level]** | **Han Q, Bradshaw EM, Nilsson B, Hafler DA, Love JC. *Lab Chip* 2010;10:1391–1400** — "Multidimensional analysis of the frequencies and rates of cytokine secretion from single cells by quantitative microengraving." VERIFIED: correct paper, correct modality (single-cell **rates** in molec/s), MeSH major topic **Monocytes/immunology**, IL-6 LOD 0.5–4 molec/s (so 10.6 is comfortably above LOD). The in-code claim of "6.5 ± 3.2 at 3 h → **10.6 ± 7.1 at 12 h**" is **not in the abstract** (it would be in a figure/table) — I could not confirm it. | **The MEAN over ACTIVELY-SECRETING cells**, deliberately *not* the peak. See EQ-6: the previous 156 was the **high-secretor tail**, 14.7× too high, and applying a top-tail rate to every activated cell double-counts the same heterogeneity that `SECRETOR_FRACTION` already encodes. |
| `S_PEAK_MOLEC_PER_S` (`:52`) | 156.0 | molec·s⁻¹·cell⁻¹ | **[MEASURED: PMID 37533643]** — value not confirmable from abstract → **[UNVERIFIED CITATION at value level]** | Portmann/Eyer 2023 (below). | **DEFINED BUT NEVER READ** (grep-confirmed: no consumer). Retained as the seed of a future per-cell secretion-rate **heterogeneity** distribution. The code explicitly states it is **"NOT a scale knob"** (`:51`). |
| `SECRETOR_FRACTION` (`:53`) | **0.039** | dimensionless | **[MEASURED: PMID 37533643]** — *paper VERIFIED; the specific 3.9% NOT confirmable from the abstract* → **[UNVERIFIED CITATION at value level]** | **Portmann K, Linder A, Oelgarth N, Eyer K. *Cell Rep Methods* 2023;3(7):100502** — "Single-cell deep phenotyping of cytokine release unmasks stimulation-specific biological signatures and distinct secretion dynamics." VERIFIED: correct paper, correct modality (**droplet microfluidics**, single-cell cytokine secretion, CRA context) — exactly what the code describes. | **Cell-INTRINSIC**, not spatially emergent. Measured on **LPS-stimulated** monocytes where *every* cell was maximally stimulated with **no** spatial constraint — and still only 3.9% secreted. This **refutes** the model's earlier hypothesis that the fraction would emerge from spatial proximity (see EQ-2 and §4.2). Omitting it over-produced IL-6 ~23× (`:60`). |
| `T_TO_MAX_MIN` (`:62`) | 150.0 | min | **[MEASURED: PMID 37533643]** — value not confirmable from abstract → **[UNVERIFIED CITATION at value level]** | Portmann/Eyer 2023 ("maximal IL-6 rate reached ~150 min", per in-code quote). | Sets `K_ON` (EQ-5). Only quantity that fixes the **timescale** of myeloid activation. |
| `K_ON_PER_HR` (`:112`) | **1.2** | hr⁻¹ per unit contact | **[DERIVED: from `T_TO_MAX_MIN`]** | `3.0/(150/60)`. **Verified this task: `1 − exp(−1.2×2.5) = 0.9502` → 95% of max at 150 min ✓.** | The `3.0` encodes the "≈3 time-constants → 95%" convention. **Not independently sourced** — it inherits `T_TO_MAX`'s provenance and adds a convention. |
| `K_OFF_PER_HR` (`:102`) | 0.10 | hr⁻¹ | **🔴 [ASSUMED: "IL-6 secretion outlasts the stimulus"]** — the in-code comment attributes it to "the measured decline of the secretor fraction with prolonged stimulation (PMID 37533643)" but **no decline rate is shown, derived, or confirmable** → **[UNVERIFIED CITATION]** | — | Chosen **slower than activation** (0.10 vs 1.2 → 12× slower) so the IL-6 program outlasts the stimulus. Sets `a_ss = K_ON/(K_ON+K_OFF) = 0.923` at unit contact, and the post-disengagement decay tail. Weakly load-bearing on the **peak** (activation saturates), materially load-bearing on the **duration**. |
| `S_MAX_PG_PER_HR` (`:111`) | **1.3307e-3** | pg·hr⁻¹·cell⁻¹ | **[DERIVED: EQ-6 from `S_MAX_MOLEC_PER_S`]** | Recomputed by importing the live module this task. | **⚠ The in-code comment on this very line says "~0.0196 pg/hr/cell (verified this session)" — that is the value at the RETIRED 156 molec/s, and is 14.7× too high. The code is right; the comment is stale.** Same stale value at `wholebody_pd.py:124`. |

### 3.2 Contact geometry

| Symbol | Value | Units | Tag | Source (verified this task) | Mechanistic rationale |
|---|---|---|---|---|---|
| `R_CONTACT_UM` (`:101`) | **14.1** | µm | **[DERIVED: r_mac + r_T]** | = 10.6 + 3.5 | **Contact**, not reach. CD40L–CD40 is membrane-bound → the cells must touch (EQ-3). Replaced `R_SYN = 30 µm`, which is a **reachability** radius and over-counted macrophage activation ~3× (producing fraction 1.92% → 0.60%, `:99`). |
| `r_macrophage` (`:94`, in comment) | 10.6 (= 21.2/2) | µm | **🔴 [MEASURED: PMID 9400735] — but for the WRONG MACROPHAGE COMPARTMENT** | **Krombach F, Münzing S, Allmeling AM, Gerlach JT, Behr J, Dörger M. *Environ Health Perspect* 1997;105(Suppl 5):1261–1263** — "**Cell size of alveolar macrophages**: an interspecies comparison." **The 21.2 ± 0.3 µm value IS VERBATIM in the abstract** (human, n=10, bronchoalveolar lavage, healthy non-smokers, Coulter volumetry). | The value is real and correctly transcribed — but it is a human **ALVEOLAR** macrophage. The code labels it "tissue macrophage". **The organs that dominate IL-6 here are spleen and bone marrow.** Alveolar macrophages are among the largest; using their radius **inflates** `R_CONTACT` and therefore contact counts. See §5.3. |
| `r_Tcell` (`:95`, in comment) | 3.5 (= 7/2) | µm | **🔴 [UNVERIFIED CITATION — TERTIARY SOURCE]** | **PMID 30571054 = Sauls RS, McCausland C, Taylor BN, "Histology, T-Cell Lymphocyte", StatPearls (NCBI Bookshelf NBK535433)** — a **study guide**, not a primary measurement. **Its abstract contains no cell diameter.** | The ~7 µm (range 6–10) figure is a standard, uncontroversial textbook value; but the cited source does not establish it. Should be re-cited to a primary morphometry paper. |
| `MyeloidIL6.__init__(r_contact_um=…)` default (`:124`) | **30.0** | µm | **🔴 [STALE / HAZARDOUS DEFAULT]** | — | **The constructor default is still the RETIRED 30 µm synapse reach**, and its docstring (`:125–127`) still says *"Default = the model's existing synapse reach R_SYN so contact is defined identically to the T:target synapse"*. **The live call site passes `r_contact_um=R_CONTACT_UM` explicitly (`wholebody_pd.py:129`), so live behaviour is correct** — but any new caller that omits the argument silently reverts to the 3×-over-counting radius. See §5.3. |

### 3.3 Systemic disposition

| Symbol | Value | Units | Tag | Source | Mechanistic rationale |
|---|---|---|---|---|---|
| `KDEG_IL6_PER_HR` (`:64`) | **0.20** | hr⁻¹ | **🔴 [FITTED: borrowed from another model — NOT measured]** | In-code cites **PMID 31268236 = Chen X, Kamperschroer C, Wong G, Xuan D, *Clin Transl Sci* 2019;12:600–608** — VERIFIED to be a **semi-mechanistic PK/PD MODELLING** paper (MeSH major topic "Models, Biological") that reports **no measured IL-6 clearance**. **Human IL-6 clearance appears to be UNMEASURED in the literature.** | `t½ = ln2/0.20 = 3.47 h`; implied **`CL = kdeg·V = 55.9 L/day`**. **THE SINGLE LARGEST ASSUMPTION IN THE IL-6 ARM.** The `PlasmaIL6` docstring at `:201` calls this clearance "**MEASURED**" — that word is **wrong** and should be removed from the source. |
| `V_PLASMA_ML` (`:65`) | **11650.0** | mL | **[DERIVED: model's own PBPK volumes]** — **VERIFIED this task** | `sum(qsp_costim_window_v2.PB.Vis) + V_PLASMA` = **8.5508 + 3.1 = 11.6508 L** (recomputed by importing the live module; hard-coded 11650.0 mL matches to 0.007%). | **The ECF DISTRIBUTION SPACE, not plasma.** IL-6 is 21 kDa (extravasates freely), has **no FcRn** to retain it, and is **produced in the interstitium**. The clinical assay samples plasma but the molecule is distributed across the whole space → the ODE denominator must be the **distribution volume**. Using plasma (3.1 L) over-states IL-6 by **3.76×**. The name `V_PLASMA_ML` is now a **misnomer** — it is `V_ECF`. |
| — (`:82–91`, comment only) | 6.4 L | L | **[RETIRED — SOURCING ERROR, do not restore]** | — | A briefly-used "IL-6 Vss" (3.5 central + 2.9 peripheral) that is in fact a **~150 kDa mAb** volume — almost certainly **tocilizumab (anti-IL-6R)**, not IL-6. It conveniently shrank the residual 2.1×. Reverted. Kept here as a provenance warning. |

### 3.4 Cell-type selection

| Symbol | Value | Tag | Rationale |
|---|---|---|---|
| `MyeloidIL6.MYELOID_TOKENS` (`:122`) | `("macrophage","monocyte","myeloid","kupffer","microglia")` | **[ASSUMED: substring token list]** | Observed labels in the real agent tables (docstring `:118–121`): `'macrophage'`, `'classical monocyte'`, `'monocyte'`, `'intermediate monocyte'`, `'Myeloid'`, `'myeloid dendritic cell'`. **Substring matching, not exact set-membership** — see §5.6. |
| `MyeloidIL6.__init__(count_scale=…)` default (`:124`) | 1.0 | **[CODE-INTERNAL]** | Identity. Live call site does not override it (`wholebody_pd.py:129`), so the census is applied exactly once, downstream. See §5.8. |

### 3.5 Parameters this subsystem does NOT own but is fully hostage to

| Symbol | Owner | Why it matters here |
|---|---|---|
| `myeloid_count_scale[organ]` | `wholebody_pd.py:42–70, 144` (`handoff/organ_myeloid_counts.json`) | **Plasma IL-6 is EXACTLY LINEAR in this scale** (`wholebody_pd.py:40–41,142–143`). If the census file is missing, every organ falls back to **1.0** (sampled counts) and IL-6 collapses. **⚠ The code states the size of that collapse TWICE, INCONSISTENTLY:** the docstring at `wholebody_pd.py:46` says **"~290,000× too low (measured: 2.38 pg/hr instead of ~5e5)"** while the runtime warning at `wholebody_pd.py:66–68` says **"~1e5× TOO LOW"**. These differ by ~3×; the `:46` figure is the one with a measured instance attached. **Neither has been re-derived here — do not quote either as established.** This is a **DATA GATE**, not a model parameter — and because the dependence is exactly linear, a run made at 1.0 can be **rescaled analytically** from the recorded per-organ `il6_prod` trace (`coupled_percell_pd.py:382–386`) without re-running. |
| `p_eng` (per-T engaged fraction) | `wholebody_pd.py:379` (QSS: `Cb/RA`) / `:486` (kinetic: `B2/RC`) | The **sole** drug-dependent input to T8. Everything the counterscreen discriminates on enters through `p_eng` and the T cells' positions. |

---

## 4. WHAT IS EMERGENT vs IMPOSED

The model's selling point is emergence. Here is precisely where it stops.

### 4.1 Genuinely EMERGENT (computed from mechanism, not handed in)

| Property | How it emerges | Where |
|---|---|---|
| **Saturation of the IL-6 dose-response** | From the **finite myeloid pool** × the per-cell `(1 − a_i)` structural ceiling. **There is no Emax, no Hill coefficient and no EC50 anywhere in this module** — grep-confirmed. | EQ-5, `:184` |
| **Per-molecule (per-construct) IL-6 differences** | From **WHERE each drug engages**. Drug enters only via `p_eng` and the T cells' coordinates; a drug whose target antigen is expressed in a myeloid-rich compartment places engaged T cells inside many macrophages' 14.1 µm contact shells. **The IL-6 difference between two constructs is an ANATOMY calculation, not a parameter.** | EQ-4, `:176–181` |
| **The spatial heterogeneity of activation** | Each myeloid cell integrates **its own** local contact — nothing is averaged. `frac_activated` is an output, never an input. | EQ-4/EQ-5 |
| **The rate→concentration conversion** | A real ODE with a real volume and a real elimination rate — the thing `IL6_SCALE` was silently standing in for. | EQ-9 |
| **The producing fraction** | `frac_activated × 0.039` — the *activated* part is emergent; the *0.039* part is not (below). | EQ-7 |

### 4.2 IMPOSED (constants handed to the subsystem — be honest about these)

| Imposed quantity | Value | Why it is NOT emergent |
|---|---|---|
| **The 3.9% secretor fraction** | 0.039 | **This was ONCE claimed to be emergent and the claim was RETRACTED IN THE SOURCE ITSELF** (`:54–61`). The measurement was made under **maximal LPS stimulation with no spatial constraint** — every cell was stimulated and still only 3.9% secreted. It is **cell-intrinsic** heterogeneity. **The module docstring at `:29–31` still claims it emerges. That claim is FALSE and is contradicted 25 lines later.** |
| **Per-cell maximal secretion rate** | 10.6 molec/s | An experimental constant. Nothing in the model derives a secretion rate from mechanism. |
| **Activation on/off rates** | `K_ON` 1.2/hr (derived from a measured 150 min), `K_OFF` 0.10/hr (assumed) | The activation *kinetics* are imposed; only the *driving input* (`contact_i`) is emergent. |
| **Contact radius** | 14.1 µm | A geometric constant from two measured cell diameters. Correctly *derived* (r+r), but the two radii are imposed — and one of them is from the wrong macrophage compartment (§5.3). |
| **IL-6 clearance and distribution volume** | `kdeg` 0.20/hr [FITTED], `V` 11.65 L [DERIVED from the model's own PBPK] | **Their product `CL = 55.9 L/day` sets the ENTIRE absolute scale.** Nothing about IL-6 disposition is emergent: there is **no** per-organ interstitial→lymph→plasma routing of IL-6 (see §5.11); it is one well-mixed compartment. |
| **Blood myeloid contribute zero** | `set_count_scale(0.0)` | A hard mechanistic gate, imposed by hand (`coupled_percell_pd.py:188`). |
| **Which cells are myeloid** | substring token list | Imposed. |

### 4.3 The one thing that is invariant to the worst assumption

Because
`C_ss = prod_total / (V·kdeg)` and `prod_total = Σ_organs (census × Σ_i a_i·S_MAX·cs_i·is_secretor_i)`,
**`V`, `kdeg` and `S_MAX` are EXACT global multiplicative constants**: they factor cleanly out of `prod_total` and
are identical for every construct. Therefore:

> **The per-construct RANKING of IL-6 (which is what a counterscreen delivers) is EXACTLY INVARIANT to `kdeg`,
> `V` and `S_MAX`. Only the ABSOLUTE pg/mL depends on them.**

**⚠ `SECRETOR_FRACTION` is NOT in that list, and an earlier version of this document wrongly included it.**
`is_secretor_i` is a **per-cell Bernoulli mask ∈ {0,1}** (EQ-2), not a scalar multiplier:
`Σ_i a_i·S_MAX·cs_i·is_secretor_i ≠ SECRETOR_FRACTION · Σ_i a_i·S_MAX·cs_i` — the two are equal only **in
expectation**. Two consequences the invariance argument must not paper over:
1. The mask is **spatially correlated with activation** by construction (it multiplies `a_i` cell-by-cell), so it
   does not factor out.
2. The mask's seed depends on `hash()`, which **varies per process** (§5.5). Each construct is scored from a
   **separate** `run_tce_pd_reval` invocation → **a different realised mask**. Per §5.5 the resulting per-organ
   noise is 2.2% at `n=50,000` but **9.1% at n=3,000 and 29% at n=300** — comparable to, or larger than, the
   inter-construct differences the counterscreen is trying to resolve in small organs. **The ranking is invariant
   to `SECRETOR_FRACTION` only in expectation, and only once the seed is made stable.** Fix the seed (§5.5) and
   this becomes an exact statement; until then it is a statistical one.

With that correction, the honest defence of the subsystem is: *the mechanistically load-bearing, drug-discriminating
part of T8 is the spatial contact calculation (EQ-4), which contains no fitted constants; the unsourced disposition
constants (`kdeg`, `V`) live entirely in a global scale factor which cancels exactly in any comparison.* The
corollary is equally honest: **any comparison to an absolute clinical anchor is a test of the global scale factor —
i.e. of `CL = 55.9 L/day` and of the myeloid census — as much as it is a test of the mechanism.**

**Derived scale check (arithmetic, not a verdict — and the anchor itself is DISPUTED, see the danger box).**
The implied activated-myeloid requirement is **exactly linear in whichever anchor you adopt**, so it is stated for
both rival values rather than for one:
```
generic:  activated_myeloid = C_ss × V × kdeg / (S_MAX_PG_PER_HR × SECRETOR_FRACTION)
                            = C_ss × 11650 mL × 0.20/hr / (1.3307e-3 × 0.039)
                            = C_ss × 4.49e7  cells per (pg/mL)

anchor = 152 pg/mL   (run_tce_pd_reval.py:77 / PROVENANCE_AND_VALIDATION.md — UNCITED)
    prod_needed  = 152 × 11650 × 0.20        = 3.542e5 pg/hr
    secretors at a=1 = 3.542e5 / 1.3307e-3   = 2.66e8 cells
    implied ACTIVATED myeloid (÷ 0.039)      = 6.8e9 cells

anchor = 127.4 pg/mL (IL6_ANCHORS_VERIFIED_2026-07-13.md — digitized db, Chen 2023 Fig2, n=212)
    prod_needed  = 127.4 × 11650 × 0.20      = 2.968e5 pg/hr
    secretors at a=1 = 2.968e5 / 1.3307e-3   = 2.23e8 cells
    implied ACTIVATED myeloid (÷ 0.039)      = 5.7e9 cells
```
**The model must therefore activate ~6 × 10⁹ myeloid cells (5.7–6.8 × 10⁹ across the two candidate anchors) to
reproduce the mosunetuzumab observation.** Whether that is physiologically attainable **cannot be decided from
anything in this module** — it requires a sourced whole-body macrophage + monocyte census, which the model does not
currently have (`organ_myeloid_counts.json` is the open data gate, §3.5). **This arithmetic is stated, not
resolved.** It is the single most useful falsification test available for the IL-6 arm, and it should be run the
moment the census lands — **against an adjudicated anchor, not against whichever value is nearest to hand.**

---

## 5. KNOWN LIMITATIONS / OPEN QUESTIONS

*Written adversarially: this is what a reviewer will attack, stated before they get to.*

### 5.1 🔴 Stale comments inside the live source contradict the live values (fix in code)

**Six** separate stale statements, all verified this task:

| Location | Says | Truth |
|---|---|---|
| `myeloid_il6.py:111` | `# ~0.0196 pg/hr/cell (verified this session)` | Live value is **1.3307e-3 pg/hr/cell** (printed from the live module). 0.0196 is the value at the **retired 156 molec/s** — **14.7× too high**. |
| `wholebody_pd.py:124` | `(0.0196 pg/hr/cell, PMID 37533643)` | Same stale value. Also mis-attributes the rate to PMID 37533643 when the live rate comes from **PMID 20376398**. |
| `myeloid_il6.py:29–31` | "THE ~3.9% secretor fraction … **EMERGES** as 'only the myeloid cells that happen to sit near engaged T cells are activated' — spatial heterogeneity IS the secretion heterogeneity" | **Explicitly refuted at `:54–61` of the same file.** It is cell-intrinsic, not spatial. |
| `myeloid_il6.py:36` | "**EVERY PARAMETER IS LITERATURE-MEASURED (none fitted, none tuned to model output)**" | **False.** `kdeg` is [FITTED] (§EQ-9a). `K_OFF` is [ASSUMED]. `K_ON` is [DERIVED]. |
| `myeloid_il6.py:201` | `PlasmaIL6`: "with **MEASURED** first-order clearance" | **False.** See EQ-9a. |
| `myeloid_il6.py:88–91` | "**Plasma volume (3 L) is the conservative, transparent choice**" | Contradicts the live `V = 11650 mL`. This paragraph is a leftover from the pre-correction state. |

**Contaminated clinical anchors still printed in the live source and MUST NOT be propagated:**
- `myeloid_il6.py:6` — "IL6_SCALE was ONE constant fit so mosunetuzumab -> **570** pg/mL" (historical, correctly framed as the *broken* legacy — acceptable as history).
- `myeloid_il6.py:34` — "**The clinical 570-vs-191 split falls out of ANATOMY.**" **🔴 BOTH numbers are contaminated: 570 has NO SOURCE; 191 IS A PAGE NUMBER (a dot-leader from the Table of Figures of FDA BLA 761345). Elranatamab has NO clinical IL-6 value in existence.** This sentence describes a validation target that does not exist and must be deleted.
- `wholebody_pd.py:23` — `CYTO_IL6_CLINICAL_ANCHOR_PGML = 570.0  # mosunetuzumab peak IL-6 ~570 pg/mL (Hosseini 2020 Fig5A)`. **570 has no source.** This constant is on the **legacy `cytokine_to_pgml` path**, not the mechanistic IL-6 path, but it is a live symbol carrying a fabricated anchor.
- **What the SURVIVING anchors are is DISPUTED inside this repo (see the danger box, item 3).** 570, 340, 230 and 366.88 have no source, and 191 is a page number — that is agreed everywhere. But:
  - `run_tce_pd_reval.py:62,74–77` + `docs/PROVENANCE_AND_VALIDATION.md:46,84` + `reference_unified_binding/score_il6_t24.py:37–39` carry **mosunetuzumab 152 / teclistamab 21** as the "two surviving population means". **`152` has NO citation anywhere in the repository.** Only `21` is verbatim-sourced (PMID 38831634).
  - `IL6_ANCHORS_VERIFIED_2026-07-13.md` — the **digitized-database** anchor set — gives mosunetuzumab MEAN = **127.4** (Chen 2023 Fig2, n=212), **drops teclistamab** (no IL-6 curve in the db), and says the "21 / 288" pair is **"NOT independently verified … do not cite them as sourced anchors."** It also supplies anchors this document did not previously mention at all: glofitamab 30.2 (median), talquetamab 19.8 / 7.9 (median), blinatumomab 640 / 370 (median).
  - **288** is described as the *highest individual patient* Cmax in `PROVENANCE_AND_VALIDATION.md:47` and `run_tce_pd_reval.py:68–70` (an order statistic that scales with cohort N → **not comparable across trials**), but `IL6_ANCHORS_VERIFIED` disputes that it is independently verified at all. **Do not assert it is "real" without resolving that.**
  - **Action:** adjudicate the two anchor files before any absolute-scale claim is made. Do not quote a single number as "the valid anchor".

### 5.2 🔴 The contact gate is evaluated on a subsampled 2D coordinate frame

`R_CONTACT = r_mac + r_T` is a **3D physical touch** criterion. It is applied to **2D projected coordinates of a
SUBSAMPLED cell population**. Two consequences:

1. **The number of engaged-T neighbours within 14.1 µm is a property of the ABM's sampling density, not of
   physiology.** If the ABM samples 1 in *k* T cells, contact counts fall ~*k*-fold. **⚠ This claim is NOT
   demonstrated anywhere in the source, and an earlier version of this section wrongly offered the radius
   sensitivity as evidence for it.** The only sensitivity datum in the code — producing fraction **1.92% at
   30 µm → 0.60% at 13 µm** (`:99`), a ratio of 3.2 against a naive areal scaling of `(30/13)² = 5.3`
   (sub-quadratic because activation saturates and cells cluster) — measures the response to **the contact
   RADIUS at fixed sampling**, which is a *different* independent variable. It establishes that **the producing
   fraction is a strong, roughly-areal function of a geometric radius applied to a synthetic coordinate frame**;
   it says **nothing** about the response to sampling density. **The sampling-density sensitivity has never been
   run** (see §5.14 Q5) — it is a structural concern, correctly raised, but currently *unquantified*.
2. **Because activation saturates (`a_ss ≈ 0.92` at contact=1), the model is far more sensitive to the *fraction of
   myeloid with ≥1 engaged neighbour* than to *how many* neighbours they have.** That fraction is precisely the
   quantity most contaminated by sampling density and by 2D projection.

**This is the deepest structural criticism of the subsystem.** The mitigation the code relies on — that the *ranking*
is preserved because the sampling geometry is identical across constructs — is sound (§4.3), but it does mean
the **absolute** producing fraction is not a physiological prediction.

**Also note:** the 0.60% figure was measured **at 13 µm**, not at the live **14.1 µm** (`:99`). No producing fraction
at the live radius is recorded in the source.

### 5.3 🔴 The contact radius is anchored on the WRONG macrophage, and its constructor default is still 30 µm

- **`r_macrophage = 10.6 µm` comes from ALVEOLAR macrophages** (PMID 9400735, Krombach 1997 — **verified: the
  21.2 ± 0.3 µm value IS in the abstract, and IS an alveolar macrophage from bronchoalveolar lavage**). The
  IL-6-dominant organs in this model are **spleen and bone marrow**. Alveolar macrophages are among the largest
  macrophage populations; if splenic red-pulp macrophages are smaller, `R_CONTACT` is inflated and (per §5.2, ~areal
  scaling) the producing fraction is inflated roughly as the square of the error.
- **`r_Tcell = 3.5 µm` is cited to a StatPearls study guide** (PMID 30571054) whose abstract contains no diameter.
  The value is uncontroversial; the citation is not admissible.
- **The `MyeloidIL6.__init__` default is STILL `r_contact_um=30.0`** (`:124`) with a docstring that still argues *for*
  the 30 µm reach (`:125–127`). Live behaviour is correct because `wholebody_pd.py:129` passes `R_CONTACT_UM`
  explicitly — but **the retired, 3×-over-counting value is one omitted keyword argument away from returning
  silently.** This should be changed to `r_contact_um=R_CONTACT_UM` (or made a required argument).

### 5.4 The activation input `contact_i` is a SUM whose units are implicitly assumed to be 1

`K_ON` was derived from a 150-min time-to-max measured on **maximally stimulated** monocytes. In the model, the
driving term is `K_ON · contact_i`, where `contact_i` is a **sum** of engagement fractions over neighbours. A cell
with 20 fully-engaged neighbours therefore activates **20× faster** than the measured maximal rate. The
calibration implicitly assumes `contact = 1` corresponds to "maximal stimulation". **There is no measurement
linking the number of contacting engaged T cells to the myeloid activation rate**, and no saturation on the *input*
side (only on the *output* side, via `(1−a)`). Because activation saturates, this distorts the **rise time**, not the
plateau — but a reviewer will ask why an input that can reach 20 is multiplied by a rate calibrated at 1.

### 5.5 🔴 The "deterministic" secretor draw is NOT reproducible across processes — VERIFIED

```python
_rng = np.random.default_rng(abs(hash(("il6_secretor", int(self.n_myeloid)))) % (2**32))   # :145
```
Python's `hash()` of a `str` (and hence of a tuple containing one) is **salted by `PYTHONHASHSEED`, which is
randomised per interpreter by default**. Three fresh interpreters in this session returned three different seeds
(**1966803464 / 2249461615 / 1998730173**); with `PYTHONHASHSEED=0` they collapse to one value. **The in-code claim
"deterministic seed … Deterministic per organ so runs are reproducible" (`:143–144`) is false as written.**

**Magnitude of the resulting run-to-run noise.** The realised secretor count is Binomial(`n_myeloid`, 0.039), so the
relative SD of an organ's IL-6 production is `sqrt((1−p)/(p·n)) = sqrt(0.961/(0.039·n))`:
- `n = 50,000` → **2.2%** (negligible)
- `n = 3,000` → **9.1%**
- `n = 300` → **29%**
- `n = 50` → **70%** (and a non-trivial chance of **zero** secretors → that organ emits **nothing**)

So for large sampled myeloid populations this is cosmetic, but **for small organs it is a material, unreported
stochastic term in a model that is otherwise deterministic.** Fix: `np.random.default_rng(<literal int seed>)`, or
seed on a stable hash (`hashlib`), or seed on the organ name with a fixed salt.

**Second-order issue:** the seed depends **only on `n_myeloid`** — so two different organs that happen to have the
same sampled myeloid count receive the **identical** secretor mask.

### 5.6 Myeloid identity is a SUBSTRING match, not exact set-membership

`MYELOID_TOKENS` are matched by `t in s` on the lowercased label (`:129`). Consequences:
- `'myeloid dendritic cell'` **is** matched (the docstring confirms it) and then secretes IL-6 at the **monocyte**
  rate with the **monocyte** secretor fraction. mDCs are not monocytes.
- Any bone-marrow label containing `"myeloid"` (e.g. a myeloid **progenitor**) would also be admitted. I **could not
  confirm** which labels exist in the marrow agent table from this module alone — **this is flagged as a risk, not
  asserted as a fact.**
- Substring matching on cell-type identity is a known fabrication-class hazard elsewhere in this lab's models
  (substring gene-symbol matching swapped gene identity). **The fix is exact set-membership against an explicit,
  audited label list.**

### 5.7 `K_OFF` is an assumption wearing a citation

The comment attributes 0.10/hr to "the measured decline of the secretor fraction with prolonged stimulation
(PMID 37533643)". **No such rate is derivable from that paper's abstract and no derivation appears in the code.**
It is tagged **[ASSUMED]** here. It is not load-bearing on the peak (activation saturates), but it *is* load-bearing
on the **duration** of the IL-6 excursion — and duration is what distinguishes a Grade 1 from a Grade 3 CRS
clinically. **If any CRS-grade readout is built on IL-6 AUC or time-above-threshold rather than peak, `K_OFF`
becomes a first-order assumption and must be sourced.**

### 5.8 Latent census double-count

`MyeloidIL6.__init__(count_scale=1.0)` (`:124`) multiplies production by `cs_i` **inside** the module (`:191`), and
`CoupledPerCellPD` multiplies by `myeloid_count_scale` **again** downstream (`coupled_percell_pd.py:287–288`). At the
live operating point this is **safe** (`wholebody_pd.py:129` does not pass `count_scale`, so `cs = 1.0`), but the two
scales occupy the same conceptual slot and **passing the census at construction would silently square it.** The two
should be unified or one should be removed.

### 5.9 🔴 The blood-myeloid gate's numerical justification is stale by ~1,400×

The gate (`coupled_percell_pd.py:188`, `set_count_scale(0.0)`) is justified in-comment by: *"it activates 98.6% of
2.0e9 monocytes and yields **61,874 pg/mL** … vs a severe-CRS ceiling of ~10–20k"* (`:183–185`). **Reconstructed this
task:** that number is reproducible **only under the retired parameterisation** (`0.0196 pg/hr/cell`, `V = 3.1 L`,
**no** secretor gate → **62,340 pg/mL**, within 1% of the quoted figure). At the **live** parameters the identical
blood configuration yields **43.9 pg/mL**.

**The "geometry artifact" the gate was justified by does not exist at the current parameters.** The *mechanistic*
arguments for the gate — (i) marginating monocytes are already counted as tissue-resident myeloid, so including
blood double-counts them; (ii) the blood ABM's coordinates are a **synthetic 2D grid** (mean spacing ~15 µm) whose
"contacts" are non-physical — are **independent of the arithmetic and may well still hold.** But the comment's
quantitative case must be re-derived or withdrawn, and **the gate should be re-justified on mechanism alone.**

### 5.10 `IL6_SCALE` is dead on the IL-6 path but is still a live symbol

The fitted path is deleted from the IL-6 output (`run_tce_pd_reval.py:215–235`) and an empty mechanistic array is a
**hard error** — this is the right design. But the fitted scale is still **loaded from a calibration JSON and exposed
as a module global**: `run_tce_pd_reval.py:34` (`IL6_SCALE=CAL['il6_eng_scale_pgml_per_raw']`) and `:41`
(`IL6_SCALE=PDC.IL6_SCALE_KIN`), fed by **`pd_model_config.py:45,59,63`** — where the symbol is spelled
**`IL6_SCALE_KIN`**, not `IL6_SCALE` (an earlier version of this section conflated the two names), with a documented
fallback literal of **0.05473** (`:63`). Separately, `wholebody_pd.py:23,32–35` still carries
`CYTO_IL6_CLINICAL_ANCHOR_PGML = 570.0` and a `cytokine_to_pgml()` converter that multiplies by it.
**Nothing on the live IL-6 path reads any of them** — but a fabricated anchor that is still a live, importable symbol
is a loaded gun. They should be deleted, not merely bypassed.

> Note: the *dead* module `calib_kdeath.py` (`:9,30`) also still reads `H.IL6_SCALE`. It is on the not-in-use
> register and is therefore **not documented here** — but it is a reason the symbol cannot simply be deleted without
> retiring that file too.

### 5.11 IL-6 is NOT routed through the PBPK the antibody uses (stated in-code, `:77–81`)

The model has, in principle, the machinery to route IL-6 the way it routes the antibody — **per-organ interstitial
production → extravasation/lymph → plasma, with IL-6's own reflection coefficient**. That machinery is **currently
BYPASSED**. What T8 implements is a **single, lumped, well-mixed ECF compartment** — an approximation of it.

Consequences a reviewer will raise:
- **No production-site → sampling-site gradient.** A cytokine made in the spleen interstitium and measured in
  plasma is treated as instantaneously well-mixed across 11.65 L. Real cytokine kinetics have a **lymphatic transit
  delay** (hours) and organ-level interstitial concentrations far above plasma.
- **No distinction between "IL-6 in the interstitium" (which drives local biology) and "IL-6 in plasma" (which the
  assay reads).** The model reports one number for both.
- **Elimination is applied to the whole ECF pool**, whereas physiologically IL-6 is cleared largely by
  receptor-mediated and hepatic/renal routes that act on specific compartments.
- The lumped treatment is **defensible** for a **peak-plasma** readout with a ~3.5 h half-life and a multi-day
  simulation (the compartment equilibrates fast relative to the CRS timecourse), and is honestly labelled as an
  approximation in the source. It is **not** defensible for any tissue-level IL-6 claim.

### 5.12 Two different single-cell assays are multiplied together

`prod ∝ SECRETOR_FRACTION × S_MAX`, i.e. **(fraction of cells that are secretors, from PMID 37533643 — droplet
microfluidics, LPS-stimulated) × (mean rate over secreting cells, from PMID 20376398 — microengraving, IL-6 LOD
0.5–4 molec/s)**. **This product is only valid if the two assays define "a secreting cell" identically.** They almost
certainly do not: different platforms, different detection limits, different stimuli, different incubation windows.
If PMID 20376398's "secreting cells" are a **broader** set than PMID 37533643's 3.9%, then the mean rate 10.6 is a
mean over a broader (and therefore lower-rate) population and the product **under**-counts; if narrower, it
**over**-counts. **The entire absolute magnitude of the IL-6 arm rides on this cross-assay composition, and it is
not validated.** The clean fix is to take **both** the secretor fraction **and** the per-secretor rate from **one**
assay — PMID 37533643 measures both (it is the source of the 3.9% *and* of the 156 molec/s tail), so a
self-consistent (fraction, mean-rate) pair from that single paper would remove this risk entirely.

### 5.13 Explicit-Euler activation vs analytic plasma step

EQ-9 (plasma) is solved **exactly**; EQ-5 (activation) is **explicit Euler + clip**. With large `contact` and a
non-small PD step, the activation update overshoots and is caught by the clip — stable, but the **rise time** of a
heavily-contacted myeloid cell is not resolved. Given the CRS timescale (days) vs the activation timescale (2.5 h),
and given that the plateau is what sets peak IL-6, this is likely benign — **but it has not been demonstrated.** A
`dt`-halving convergence check on `il6_peak` would settle it and is not in the source.

### 5.14 Open questions (the honest list)

0. **🔴 WHICH CLINICAL ANCHOR IS CORRECT — 152 or 127.4?** Two in-repo documents dated the same day disagree, the
   live code uses the **uncited** one (152), and the repo's own **digitized-database** file gives 127.4 and forbids
   citing the teclistamab pair at all. **Nothing downstream of "is the absolute scale right?" can be settled until
   this is adjudicated.** This is now the top open question in the subsystem, ahead of the census. (Danger box
   item 3; §5.1.)
1. **What IS human IL-6 clearance?** No measurement was found. Until one exists, `CL = 55.9 L/day` is an assumption
   and every absolute pg/mL inherits it.
2. **Does the whole-body myeloid census support the ~6 × 10⁹ activated cells that the mosunetuzumab anchor
   requires** (§4.3; 5.7e9 at anchor 127.4, 6.8e9 at anchor 152)? This is the falsification test, and it is blocked
   on **both** `organ_myeloid_counts.json` **and** Q0.
2b. **What is the mechanistic basis of the contact gate?** The CD40L–CD40 axis — the sole justification for using a
   *contact* radius rather than a *reach* radius — is **uncited anywhere in the module** (EQ-3). If the gate is
   actually IL-1/NO/soluble-mediator-driven, EQ-3 has the wrong functional form.
3. **What is the diameter of a SPLENIC red-pulp macrophage and a MARROW macrophage?** `R_CONTACT` should be built
   from the macrophages that actually dominate the signal, not from alveolar macrophages.
4. **Are the 3.9% secretors a fixed subset, or a stochastic per-episode state?** The model assumes a **fixed cell
   identity**. If secretion is instead a stochastic burst state that any monocyte can enter, the *fraction* is right
   but the *spatial correlation* with engaged T cells is wrong.
5. **Does the producing fraction converge as ABM sampling density increases?** Untested (§5.2).
6. **Should `K_OFF` be sourced?** Required if any readout is AUC- or duration-based rather than peak-based (§5.7).
7. **Is `frac_activated`'s 0.05 reporting threshold (`:196`) ever mistaken for a producing fraction?** It is not
   secretor-gated; the true producing fraction is `frac_activated × 0.039`.

---

## APPENDIX — VERIFICATION LOG (this task)

```
SOURCE READ IN FULL:  engine/myeloid_il6.py  (216 lines, top to bottom)
CONSUMERS READ:       engine/wholebody_pd.py (:1-150, :365-390, :480-495)
                      engine/coupled_percell_pd.py (:160-230, :270-300, :380-391)
                      engine/run_tce_pd_reval.py (:195-240)
                      engine/qsp_costim_window_v2.py (:120-165, :362)
DEAD MODULES NOT DOCUMENTED (per brief): cytokine_pbpk.py, il6_pbpk.py, unified_binding.py,
    multiarm_kinetic.py, biexact_solver.py, rna_to_receptor.py, convert_copies_ALL.py, calib_kdeath.py

EXECUTED (conda claude-skills python, live modules imported):
  qsp_costim_window_v2.PB.Vis.sum()          = 8.5508 L
  qsp_costim_window_v2.V_PLASMA              = 3.1 L
  sum                                        = 11.6508 L   vs V_PLASMA_ML=11650.0 mL   -> MATCH (0.007%)
  myeloid_il6.S_MAX_PG_PER_HR                = 1.3306895868704338e-3 pg/hr/cell
  myeloid_il6.s_max_pg_per_hr(156.0)         = 1.9583733542621478e-2 pg/hr/cell   [== the STALE :111 comment]
  ratio 156/10.6                             = 14.72x
  myeloid_il6.K_ON_PER_HR                    = 1.2 /hr
  1 - exp(-1.2*2.5)                          = 0.9502            [confirms the ":112" 95%-at-150min claim]
  a_ss at contact=1 = 1.2/(1.2+0.10)         = 0.9231
  ln2/0.20                                   = 3.465 h           [confirms "t1/2 ~3.5 h"]
  CL = 0.20/hr * 11.65 L                     = 2.33 L/hr = 55.92 L/day
  C_ss per unit prod = 1/(11650*0.20)        = 4.2918e-4 pg/mL per (pg/hr)

  PYTHONHASHSEED determinism test (3 fresh interpreters):
    abs(hash(("il6_secretor",5000)))%2**32   = 1966803464 / 2249461615 / 1998730173   -> NOT deterministic
    with PYTHONHASHSEED=0                    = 733582674 / 733582674                  -> deterministic

  Blood-gate arithmetic reconstruction:
    RETIRED params (0.0196 pg/hr/cell, V=3.1L, no secretor gate):
       2.0e9 * 0.986 * 0.0196 / (3100*0.20)  = 62,340 pg/mL   [comment says 61,874 -> reproduces to 1%]
    LIVE params (1.3307e-3 pg/hr/cell, 0.039 secretor gate, V=11650 mL):
       2.0e9 * 0.986 * 0.039 * 1.3307e-3 / (11650*0.20) = 43.9 pg/mL

  mosunetuzumab-152 requirement:
    prod_needed = 152 * 11650 * 0.20         = 3.542e5 pg/hr
    secretors at a=1 = 3.542e5 / 1.3307e-3   = 2.66e8 cells
    implied activated myeloid (/0.039)       = 6.82e9 cells

PMIDs RESOLVED against Europe PMC (bc_get_europepmc_articles, EXT_ID query):
  29808005  VERIFIED  Giavridis T et al. Nat Med 2018;24:731-738. "CAR T cell-induced cytokine release
                      syndrome is mediated by macrophages and abated by IL-1 blockade."
                      Abstract verbatim: severity "mediated not by CAR T cell-derived cytokines, but by
                      IL-6, IL-1 and nitric oxide (NO) produced by recipient macrophages".
                      NOTE: the code's quoted line "IL-6 induction and myeloid activation require proximity
                      of CAR T cells and myeloid cells" is NOT in the abstract -> quote UNVERIFIED (would be
                      in the body). The CD40L-CD40 attribution has NO citation anywhere in myeloid_il6.py.
                      NOTE: murine CRS model (recipient MOUSE macrophages).
  29808007  VERIFIED  Norelli M et al. Nat Med 2018;24:739-748. "Monocyte-derived IL-1 and IL-6 are
                      differentially required for cytokine-release syndrome and neurotoxicity due to CAR
                      T cells." Abstract VERBATIM: "Human monocytes were the major source of IL-1 and IL-6
                      during CRS."  [exact match to the in-code quote]
  20376398  VERIFIED (paper) / UNVERIFIED (value)
                      Han Q, Bradshaw EM, Nilsson B, Hafler DA, Love JC. Lab Chip 2010;10:1391-1400.
                      "Multidimensional analysis of the frequencies and rates of cytokine secretion from
                      single cells by quantitative microengraving." Correct modality (molec/s per cell);
                      MeSH major topic Monocytes/immunology; IL-6 LOD 0.5-4 molec/s. The specific
                      "6.5+/-3.2 @3h -> 10.6+/-7.1 @12h" is NOT in the abstract.
  37533643  VERIFIED (paper) / UNVERIFIED (value)
                      Portmann K, Linder A, Oelgarth N, Eyer K. Cell Rep Methods 2023;3(7):100502.
                      "Single-cell deep phenotyping of cytokine release unmasks stimulation-specific
                      biological signatures and distinct secretion dynamics." Droplet microfluidics,
                      single-cell cytokine secretion -> exactly the modality the code describes. The
                      specific 3.9% / 150 min / 156 molec/s values are NOT in the abstract.
  31268236  VERIFIED -- AND IT REFUTES THE IN-CODE TAG.
                      Chen X, Kamperschroer C, Wong G, Xuan D. Clin Transl Sci 2019;12(6):600-608.
                      "A Modeling Framework to Characterize Cytokine Release upon T-Cell-Engaging
                      Bispecific Antibody Treatment: Methodology and Opportunities."
                      Abstract: "A 'fit-for-purpose' SEMIMECHANISTIC pharmacokinetic/pharmacodynamic MODEL
                      was developed..."  MeSH MAJOR TOPIC: "Models, Biological".
                      -> reports NO measured IL-6 clearance. kdeg = 0.20/hr is [FITTED], NOT [MEASURED].
  9400735   VERIFIED -- BUT WRONG MACROPHAGE COMPARTMENT.
                      Krombach F, Munzing S, Allmeling AM, Gerlach JT, Behr J, Dorger M.
                      Environ Health Perspect 1997;105(Suppl 5):1261-1263.
                      "Cell size of ALVEOLAR macrophages: an interspecies comparison."
                      Abstract VERBATIM: human AM diameter "21.2 +/- 0.3 microns (n = 10)", from
                      bronchoalveolar lavage of healthy non-smoking volunteers, Coulter volumetry.
                      -> value real, correctly transcribed, but it is an ALVEOLAR macrophage, not the
                         splenic/marrow macrophage that dominates IL-6 in this model.
  30571054  VERIFIED -- TERTIARY SOURCE, DOES NOT SUPPORT THE VALUE.
                      Sauls RS, McCausland C, Taylor BN. "Histology, T-Cell Lymphocyte." StatPearls
                      (NCBI Bookshelf NBK535433). A STUDY GUIDE. Its abstract contains NO cell diameter.
                      -> the 7 um T-cell diameter is uncontroversial but this citation does not establish it.

NOT VERIFIABLE FROM THIS MODULE:
  - the CD40L-CD40 mechanistic attribution for contact-gating: NO citation anywhere in myeloid_il6.py.
  - IL6_MW_DA = 21000: no PMID in code (uncontroversial).
  - K_OFF_PER_HR = 0.10: no derivation shown; the attributed source does not (from its abstract) support it.
```

---

## APPENDIX B — ADVERSARIAL RE-VERIFICATION (second pass, 2026-07-13)

An independent adversarial pass re-checked **every** equation, line citation, parameter value, provenance tag and
PMID in this document against the live source. Recorded here in full, including what it found wrong.

**RE-EXECUTED AND REPRODUCED (all confirm the doc):**
```
myeloid_il6.S_MAX_PG_PER_HR        = 0.001330689586870434   -> doc "1.3307e-3"      MATCH
myeloid_il6.s_max_pg_per_hr(156.0) = 0.019583733542621478   -> doc "0.0196"         MATCH (the stale :111 comment)
156/10.6                           = 14.717                 -> doc "14.7x"          MATCH
K_ON_PER_HR                        = 1.2                    -> doc                  MATCH
1-exp(-1.2*2.5)                    = 0.950213               -> doc "0.9502"         MATCH
a_ss(contact=1)=1.2/1.3            = 0.923077 ; (contact=10)=0.991736  -> doc "0.923 / 0.992"  MATCH
ln2/0.20                           = 3.4657 h               -> doc "3.465 h"        MATCH
CL = 0.20 * 11.65 * 24             = 55.92 L/day            -> doc "55.9 L/day"     MATCH
11.65/3.1                          = 3.7581                 -> doc "3.76x"          MATCH
qsp_costim_window_v2.PB.Vis.sum()  = 8.5508 L               -> doc "8.5508 L"       MATCH
qsp_costim_window_v2.V_PLASMA      = 3.1 L                  -> doc "3.1 L"          MATCH
  sum = 11.6508 L vs hard-coded 11650.0 mL                  -> doc "0.007%"         MATCH
blood @ RETIRED params             = 62,340 pg/mL           -> doc "62,340"         MATCH
blood @ LIVE params                = 43.92 pg/mL            -> doc "43.9"           MATCH
binomial relSD n=50k/3k/300/50     = 2.2 / 9.1 / 28.7 / 70.2 %  -> doc "2.2/9.1/29/70%"  MATCH
PYTHONHASHSEED=0 seed              = 733582674              -> doc "733582674"      MATCH
  (3 unseeded interpreters -> 3 different seeds: NON-DETERMINISM CONFIRMED INDEPENDENTLY)
grep Emax|EC50|Hill in myeloid_il6.py -> comments only, never in the dynamics       doc claim UPHELD
grep S_PEAK_MOLEC_PER_S repo-wide     -> defined, zero consumers                    doc claim UPHELD
```

**ALL 7 PMIDs RE-RESOLVED against Europe PMC. ZERO FABRICATED CITATIONS.** Every PMID in this document appears in
`myeloid_il6.py`, and every bibliographic expansion (authors, journal, volume, pages, year, MeSH, abstract quotes)
was independently confirmed correct — including the two load-bearing adversarial findings (31268236 is MeSH-major
"Models, Biological" and reports no measured clearance; 9400735's 21.2 ± 0.3 µm is verbatim in the abstract and IS
alveolar). The Giavridis and Norelli abstract quotes are **verbatim exact**.

**ERRORS FOUND AND FIXED IN THIS DOCUMENT:**
```
1. [CRITICAL] Clinical anchors. The doc asserted "the only valid anchors are mosunetuzumab 152 and
   teclistamab 21" and called 288 "real". FALSE/OVERCLAIMED:
     * 152 has NO citation anywhere in the repo (not in code, not in PROVENANCE_AND_VALIDATION.md).
     * model/IL6_ANCHORS_VERIFIED_2026-07-13.md -- the digitized-db anchor set, which this doc never
       mentioned -- gives mosunetuzumab MEAN = 127.4 (Chen 2023 Fig2, n=212), DROPS teclistamab, and says
       the 21/288 pair is "NOT independently verified ... do not cite them as sourced anchors."
   -> The doc committed, inside its own anti-fabrication danger box, the exact error class it warns about.
   FIXED: danger box item 3 rewritten; §5.1 anchor bullet rewritten; §4.3 arithmetic restated for BOTH
   candidate anchors (5.7e9 vs 6.8e9 activated myeloid); new top open question §5.14 Q0.
2. [MAJOR] §4.3 claimed SECRETOR_FRACTION is a "GLOBAL multiplicative constant" and the ranking is
   "INVARIANT" to it. FALSE: is_secretor_i is a per-cell Bernoulli MASK, not a scalar; it does not factor
   out of the sum, it correlates with a_i cell-by-cell, and its seed varies per PROCESS (§5.5) so separate
   per-construct runs draw DIFFERENT masks. §4.3 directly contradicted §5.5. FIXED (V/kdeg/S_MAX remain
   exact; SECRETOR_FRACTION is invariant only in expectation, and only once the seed is stabilised).
3. [MAJOR] §5.2 offered the RADIUS sensitivity datum (1.92%@30um -> 0.60%@13um) as evidence for the
   SAMPLING-DENSITY claim. Different independent variable; the datum does not support the claim. The
   sampling-density sensitivity has never been run. FIXED.
4. [MODERATE] §3.5 attributed "IL-6 ~1e5x too low ... 2.38 pg/hr instead of ~5e5" to wholebody_pd.py:66-68.
   The 2.38 pg/hr instance and the figure "~290,000x" are at :46; :66-68 says "~1e5x". The CODE ITSELF
   carries two inconsistent magnitudes (~2.9e5x vs ~1e5x) and the doc reported only one. FIXED (both cited,
   inconsistency flagged, neither quoted as established).
5. [MODERATE] EQ-3 presented CD40L-CD40 as established mechanism. It is UNCITED anywhere in myeloid_il6.py
   (the doc's own Appendix A said so, but the body did not) -- and it is the sole justification for a CONTACT
   radius over a REACH radius, i.e. for the most load-bearing constant in the subsystem. FIXED (flagged
   inline; new open question §5.14 Q2b).
6. [MINOR] Mistagged provenance: NA = 6.02214076e23 tagged "[MEASURED: SI definition]" -- Avogadro has been
   an EXACTLY DEFINED constant since the 2019 SI redefinition, not a measured one -> [DEFINED: SI exact].
   IL6_MW_DA tagged "[MEASURED]" while the doc itself admitted "No PMID in code" -> [TEXTBOOK - NO SOURCE
   IN CODE]. FIXED.
7. [MINOR] Line-citation errors. §1.3 cited the hard-error block as run_tce_pd_reval.py:216-221; actual is
   215-220 (il6=np.array at :215). The removal narrative is :199-214, not :199-223. §3.5 cited
   coupled_percell_pd.py:383-389 for the per-organ recorder; actual :382-386. FIXED.
8. [MINOR] §5.1 said "Four separate stale statements" above a table containing SIX. FIXED.
9. [MINOR] §5.10 said "the symbol IL6_SCALE ... (pd_model_config.py:45,59,63)". The symbol there is
   IL6_SCALE_KIN. FIXED (both names given; fallback literal 0.05473 at :63 noted; the dead calib_kdeath.py
   dependency noted as the reason the symbol cannot just be deleted).
10. [MINOR] §1.1 presented Giavridis as the human premise without noting it is a MURINE model (the caveat
    was in Appendix A only). FIXED.

DEAD-MODULE CHECK: PASS. No content from cytokine_pbpk, il6_pbpk, unified_binding, multiarm_kinetic,
   biexact_solver, rna_to_receptor, convert_copies_ALL or calib_kdeath is documented. calib_kdeath is
   mentioned in §5.10 ONLY as a dead-symbol dependency and is explicitly labelled not-documented.
"VALIDATED" CHECK: PASS. The doc claims no validation it does not have; the one green check (V = 11.65 L
   against the live PBPK arrays) reproduces exactly.
```
