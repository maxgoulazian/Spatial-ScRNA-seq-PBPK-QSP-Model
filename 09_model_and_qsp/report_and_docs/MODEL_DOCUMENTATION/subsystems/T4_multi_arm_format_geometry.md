---
title: "T4 — Multi-arm format geometry (cis/trans co-engagement, reach kernel, linker span, cleft feasibility)"
subsystem: T4
model: costim_engager_counterscreen
primary_source: engine/multiarm_binding.py (128 lines, read in full)
cross_referenced_live_source: engine/kinetic_synapse.py, engine/wholebody_pd.py, engine/pd_model_config.py
date: 2026-07-13
generated_by: workflow-subagent T4
adversarially_verified: 2026-07-13 — every equation, file:line, and parameter value re-checked against live code; 7 defects found and fixed (see AUDIT note below)
provenance_rule: "every number below was read out of, or computed by executing, the live code in this task. NO number may be introduced from memory, from another agent's brief, or from an un-cited 'audit' — including numbers offered as CORRECTIONS to code values."
---

> [!danger]+ AUDIT 2026-07-13 — defects found and corrected in this doc
> 1. **FABRICATED-PROVENANCE NUMBERS (most serious).** §3.6 previously asserted *"Per the provenance audit of 2026-07-13, 570 has no source; the only valid IL-6 clinical anchors are mosunetuzumab **152** and teclistamab **21** (both population means)."* **Neither 152 nor 21 exists anywhere in this repository** (grepped: engine, docs, params, SHARED_LOG, CHANGELOG). They were stated as established fact, attributed to a named audit, in a doc that owns no IL-6 parameter — a textbook provenance fabrication, and a violation of this doc's own frontmatter rule. **Both numbers removed.** A number offered as a *correction* is still a number and still needs a source.
> 2. **Wrong constant in a RUN-VERIFY** (EQ-1): the closed-form check used `2.490811`; the true value is `2.49080860076077`. With the correct constant the closed form matches the code to **1 ULP**, not "5 significant figures" as claimed. Corrected.
> 3. **Off-by-one file:line citations** (×3): the live cis gate is `wholebody_pd.py:295`, not `:296`; range `282–295`, not `282–296`; §4.1 pointed at `:288` (a docstring line) instead of the computation at `:294–295`. Corrected.
> 4. **Internal contradiction**: EQ-9(c) listed five formats then called them "the four formats". Corrected to five.
> 5. **Mis-tagged provenance**: `AVO` was tagged **[MEASURED]**. Avogadro's number is *exact by SI definition* since 2019 — a defined constant, not a measurement. Retagged **[DEFINED — SI exact]**. No [MEASURED] tag now survives in T4, which is the correct state: **T4 owns no measured parameter.**
> 6. **Under-verified citation**: the doc declared "Rhoden 2016" unverifiable. It is in fact cited in full in the live codebase (`qsp_costim_window_v2.py:704`: *J Biol Chem 291:11337-47*). Upgraded to [SOURCED-IN-REPO, UNSOURCED-IN-FILE] — while keeping the honest caveats that (a) it was not checked against the paper, and (b) it is attached to the *binding law*, not to the 12.5 nm reach or the 8.0 µm radius, which remain uncited.
> 7. **Rounding**: EQ-9 table `s=8` fraction 0.063 → 0.0625 (exact).
>
> **Checked and found CORRECT (no change):** all EQ-1…EQ-10 equations against `multiarm_binding.py`; every other file:line; every parameter default (`KD_CD3=3.0`, `KD_TAA=0.3`, `KD_cos=1.0`, `rec_*`, `r_cell_T=4.0`, `r_cell_tum=8.0`, spans); the entire EQ-9 span table and its 13.0 nm argmax; the EQ-5 bivalent table; all `p_cis` and feasibility values; the EQ-3a/EQ-3b discrepancy (0.5 vs 0.0 at 6.5 nm) — **real, reproduced**; the claim that no live module calls `multiarm_bound`/`bound_for_format`/`FORMATS` — **verified by grep** (`wholebody_pd.py:208` imports only `_cis_feasibility`). **No dead module** (`cytokine_pbpk`, `il6_pbpk`, `unified_binding`, `multiarm_kinetic`, `biexact_solver`, `rna_to_receptor`, `convert_copies_ALL`, `calib_kdeath`) is documented here.

up:: [[00_INDEX]]
tags:: #atlas/generated #model-doc

> [!warning]+ Read this box before reading the doc — it defines what is LIVE
> `multiarm_binding.py` is **two things in one file**, and only one of them executes in the canonical run:
> 1. **`_cis_feasibility`** (`multiarm_binding.py:39–46`) — **LIVE**. It is the *only* symbol any module in the declared live import graph imports from this file (`wholebody_pd.py:208`, inside a `try:`). It sets `OrganPD.p_cis`, which gates costim occupancy on same-cell CD3 engagement (`wholebody_pd.py:282–296`).
> 2. **`geo_ageff_nM` / `_bridge_feasibility` / `bound_arm` / `multiarm_bound` / `bound_for_format` / `FORMATS`** — a self-contained **format-driven binding calculator**. Verified by grep over the live engine: **no live module calls any of them.** It is the design-sweep / counterscreen entry point, not part of the clinical-validation execution path.
>
> **But the physics of item 2 is NOT dead.** `geo_ageff_nM` (`multiarm_binding.py:29–37`) is *numerically identical* to `kinetic_synapse.ageff_nM` (`kinetic_synapse.py:53–65`) — verified this task: both return `61226.6206637434…` nM at `(rec=257000, r_cell=8 µm, span=13 nm)`, agreeing to **1 ULP** (relative difference `2.4e-16`; the two files reduce the same algebra in a different operation order, so they are float-equal but not always bit-equal). The 1/span dilution law documented here **is** the law that runs live inside `KineticSynapse`, at `kinetic_synapse.py:131`, where it is multiplied by the cleft-feasibility gate. **That product is where the analytic span optimum bites the live model** (§2, EQ-9). This doc therefore documents the law in its home module (`multiarm_binding.py`) and follows it into its live consumer.

---

# T4 — Multi-Arm Format Geometry

## 1. PURPOSE & DATAFLOW POSITION

### 1.1 What the subsystem does

T4 answers one question: **given the physical architecture of a T-cell engager — how many arms, and how far apart they sit — what fraction of the binding it could do does it actually do?**

Everything upstream of T4 (PK, per-cell drug concentration, receptor copies) is *chemistry and abundance*. Everything downstream (kill hazard, serial killing, cytokine, costim program) is *signalling*. T4 is the layer in between where **format becomes physics**. It converts three architectural distances into three dimensionless modifiers:

| Construct distance (nm) | Physical question it answers | Modifier it produces |
|---|---|---|
| `span_bridge_nm` | Can the T-arm and the tumor-arm hold open the immune-synapse cleft between two apposed membranes? | cross-cell bridge feasibility → `Cb_kill`, `Cb_costTAA` |
| `span_coeng_T_nm` | Can *one* molecule reach a CD3 epitope **and** a costim epitope on the **same** T-cell surface? | `p_cis` (cis/trans coincidence) → costim gating |
| `span_coeng_tumor_nm` | Can *one* molecule reach two TAA epitopes on the **same** tumor cell? | bivalent avidity via `c_eff` |

The same three spans also set the **effective local concentration** (`c_eff`, "ageff") that a *second, already-tethered* arm experiences — and that is where the counter-intuitive result of this subsystem lives (§2, EQ-1 and EQ-9).

### 1.2 Position in the life of the molecule

```
   PK / transport            T4 (this subsystem)                  PD / signalling
  ─────────────────   ───────────────────────────────────   ────────────────────────────
   C_free (nM) per cell ─┐
   R_CD3, R_TAA, R_cos  ─┼─► [1/span dilution law: geo_ageff_nM]
   (copies/cell)         │        └─► c_eff (2nd-arm local conc, nM)
   r_cell (µm)          ─┘                │
                                          ├─► × cleft_feasibility(span_bridge)   ──► c_eff,trans
                                          │      (kinetic_synapse.py:131)             │
                                          │                                            ▼
                                          │                            kf = kon_TAA · c_eff,trans
                                          │                            → bridged trimer B2 → KILL, CRS
                                          │
                                          ├─► bivalent avidity quadratic (bound_arm) ──► Cb3/CbT/CbC
                                          │
                                          └─► p_cis = _cis_feasibility(span_coeng_T)  ──► LIVE:
                                                 wholebody_pd._apply_cis_coincidence
                                                 occ_eff = occ·[(1−p_cis) + p_cis·f_cd3]
                                                 → costim program (effector / supp / exh)
```

**Feeds T4:** per-cell free drug `Cfree` (nM), per-cell receptor copies (`rec_CD3`, `rec_TAA`, `rec_cos`), cell radii (`r_cell_T`, `r_cell_tum`, µm), arm affinities (`KD_*`, nM), and the three construct spans (nm) — supplied per format from `FORMATS` (`multiarm_binding.py:101–120`) or, per the module docstring (`multiarm_binding.py:5`), by an **AF3-derived override** per construct.

**T4 feeds:** (a) `Cb_kill` (CD3·drug·TAA bridged trimer → kill + signal-1 cytokine), (b) `Cb_costim` / `Cb_costTAA` (costim program, TAA-anchored signal-2), (c) `Cb_cis` (coincident signal 1+2 on the same T cell), and — in the live path — (d) the scalar `p_cis` consumed by `wholebody_pd.OrganPD`.

### 1.3 Live wiring, exactly

| Live line | What happens | Consequence at the canonical operating point |
|---|---|---|
| `wholebody_pd.py:208` | `from multiarm_binding import _cis_feasibility as _cisf` | the single live import from this module |
| `wholebody_pd.py:209` | `self.p_cis = _cisf(span_coeng_T_nm) if (n_costim>0 and span_coeng_T_nm is not None) else 0.0` | `span_coeng_T_nm` defaults to `None` (`wholebody_pd.py:116`) and **is not set by `pd_model_config.KINETIC`** (`pd_model_config.py:32–42` contains `span_bridge_nm`, `span_cis_nm` — **no** `span_coeng_T_nm`). ⇒ **live `p_cis = 0.0`**: costim drive is cell-autonomous (trans/legacy) unless a construct sweep explicitly supplies the T-side co-engagement span. |
| `wholebody_pd.py:282–295` | `_apply_cis_coincidence`: `occ_eff = occ·[(1−p_cis) + p_cis·f_cd3]` (`f_cd3` at `:294`, the gate itself at `:295`) | with `p_cis=0` this returns `occ` unchanged (byte-identical, early-return at `:291`). The cis machinery is **built and wired but default-inert**. |
| `kinetic_synapse.py:128–131` | `cleft = clip(span_bridge, 13, 40)`; `feas = cleft_feasibility(span_bridge, cleft)`; **`ceff_trans = ageff_nM(dens, 8 µm, span_bridge) · feas`** | this is the **live** realisation of T4's central product. `span_bridge = 12.5 nm` (`pd_model_config.py:38`) ⇒ `feas = 0.9038`, `c_eff,trans = 0.9038 × ageff(12.5 nm)`. |

> **Honest statement of scope:** the `FORMATS` library and `multiarm_bound` are *code-complete and unit-consistent but unexercised by the live clinical run*. They are the machinery a format counterscreen would call. Nothing in this doc should be read as "the validated clinical model predicts X for a BiTE" — it should be read as "the format layer, as written, predicts X, and here is the line it is written on."

---

## 2. GOVERNING EQUATIONS

Notation: `s` = span (nm); `rec` = receptor copies/cell; `r_c` = cell radius (µm); `C` = free drug (nM); `KD` = arm equilibrium dissociation constant (nM); `N_A = 6.02214076e23` (`multiarm_binding.py:24`). All concentrations nM unless stated.

---

### EQ-1 — The 1/span dilution law: geometric effective 2nd-arm concentration (`multiarm_binding.py:29–37`)

```python
SA_cell = 4π·r_c²                    # µm²   (whole-cell surface)
r_Ab    = max(s, 1e-3)·1e-3          # µm    (span converted nm→µm; 1e-3 nm floor)
SA_Ab   = π·r_Ab²                    # µm²   (footprint DISK the tethered arm can sweep)
V_Ab    = (2/3)·π·r_Ab³              # µm³   (HEMISPHERE the tethered arm explores)
Am_cell = rec·1e9/N_A                # nmol/cell   (receptor copies → nanomoles)
c_eff   = (Am_cell/SA_cell)·SA_Ab / V_Ab · 1e15     # nM
```

**Closed form** (algebraically reduced from the code; `SA_Ab/V_Ab = 3/(2·r_Ab)` exactly):

$$
c_{\text{eff}}(s)\;[\text{nM}] \;=\; \underbrace{\frac{3\times10^{24}}{2\,N_A}}_{=\,2.4908086}\;\cdot\;\frac{rec}{4\pi r_c^2 \cdot s_{\mu m}}
\;\;\propto\;\; \frac{\text{surface density}}{s}
$$

**RUN-VERIFY (executed this task, conda `claude-skills`):** the constant is `3e24/(2·AVO) = 2.49080860076077` (computed from `multiarm_binding.AVO`). The closed form reproduces the code to **1 ULP, not merely to a few significant figures**: `2.49080860076077 × 257000 / (4π·8²·0.0125) = 63675.68549029321` vs code `geo_ageff_nM(257000, 8.0, 12.5) = 63675.68549029322` nM. And the product `c_eff·s` is **invariant**: at `rec=1e5, r_c=8 µm`, `c_eff·s = 309706.64149` for every `s ∈ {5, 6.5, 7.8, 10, 12.5, 13, 20, 30, 40, 60}` nm — i.e. **exactly 1/s, no approximation.**

- **(b) BIOLOGICAL MEANING.** One arm of the engager is already bound. The other arm is now not free in bulk solution — it is *tethered*, and it samples only the volume its linker lets it reach: a hemisphere of radius `s` sitting on the membrane. Inside that hemisphere it "sees" the receptors that lie under the hemisphere's footprint. `c_eff` is the concentration those receptors represent **from the point of view of the tethered arm**: (moles of receptor under the footprint) ÷ (litres of explored volume). It is the quantity that replaces bulk `C` in the second binding event. At the tumor anchor (`rec = 257,000` CEACAM5 copies, `r_c = 8 µm`, `s = 12.5 nm`) it is **63,676 nM ≈ 64 µM** — four to five orders of magnitude above any achievable plasma concentration. *That* is avidity.

- **(c) MECHANISTIC RATIONALE — and why longer is worse.** The footprint area grows as `s²`; the explored volume grows as `s³`. Concentration = (things ∝ s²) / (volume ∝ s³) = **∝ 1/s**. A longer linker does bring more receptors into range — but it dilutes itself faster than it gains them, because it must search a volume that grows one power of `s` faster than the membrane it is searching. **A flexible arm is its own worst enemy: every extra nanometre of reach is paid for in three dimensions and repaid in only two.**
  *Alternative rejected:* a 3-D "swept volume ∝ number of accessible partners" formulation (which would make `c_eff` span-independent) is wrong here because the partners are **confined to a 2-D membrane**, not distributed through the swept volume. The code makes exactly this choice explicit — the numerator is a *2-D areal* quantity (`SA_Ab`), the denominator a *3-D volume* (`V_Ab`).
  *Alternative also rejected:* worm-like-chain / Gaussian-chain end-to-end distributions for the linker. The code uses a **uniform hemisphere** — a hard, unweighted reach envelope. See §5(1).

- **(d) UNITS.** `s`, `r_Ab` [nm→µm]; `r_c` [µm]; `SA_cell`, `SA_Ab` [µm²]; `V_Ab` [µm³]; `rec` [copies/cell]; `Am_cell` [nmol/cell]; `c_eff` [nM]. The `1e15` is µm³→L; the `1e9` is mol→nmol. The intermediate `_cc = 1e9` (`multiarm_binding.py:32`) multiplies at line 35 and divides at line 36 — **it cancels exactly** and is a no-op carried over from the sibling implementation.

- **Provenance of the construction:** the module docstring (`multiarm_binding.py:21, 30`) attributes it to **"the SAME Rhoden-2016 geometric ageff used by percell_binding / kinetic_synapse."** No PMID/DOI appears anywhere in **this file**. It does, however, appear **elsewhere in the live codebase** — `engine/qsp_costim_window_v2.py:704` carries the full reference *"Rhoden JJ, Dyas GL, Wroblewski VJ (2016) J Biol Chem 291:11337-47"*, and `:231` cites a preprint recipe DOI `10.1101/2022.09.12.507653`. → tag: **[SOURCED-IN-REPO, UNSOURCED-IN-FILE: Rhoden 2016 JBC 291:11337-47, per `qsp_costim_window_v2.py:704`]**. Two caveats a reviewer should hold: (i) the identifier is *transcribed from the codebase*, **not** independently checked against the paper in this task — I did not read Rhoden 2016; (ii) what Rhoden supplies is the *bivalent-binding QSS with a geometric effective concentration*; that the specific **uniform-hemisphere** reduction used here (`SA_Ab/V_Ab = 3/2r`) is Rhoden's own construction is **asserted by the docstring, not demonstrated**. The *geometry* is closed-form and self-evident from the code; the *attribution* is repo-internal hearsay until someone opens the paper.

---

### EQ-2 — Cis co-engagement feasibility (Gaussian span-match) (`multiarm_binding.py:39–46`)

```python
if span_coeng_nm is None: return 0.0
d = (span_coeng_nm − gap_match_nm) / max(tol_nm, 1e-6)
p_cis = exp(−0.5·d²)
```

$$p_{\text{cis}}(s)=\exp\!\left[-\tfrac{1}{2}\left(\frac{s - 12.5}{8.0}\right)^{2}\right] \in (0,1]$$

- **(b) BIOLOGICAL MEANING.** Can **one** engager molecule simultaneously grip CD3 and a costim receptor on the **same T cell**? Only if the distance between its two T-side binders matches the actual gap between the two epitopes as they present on the membrane. The gap depends on **epitope height**: a membrane-proximal costim epitope sits next to CD3 and a compact construct bridges them; a **tall** epitope (the module cites CRD1 of 4-1BB) sits far above the membrane, the required geometry no longer matches a compact span, and the molecule cannot co-grip. `p_cis → 0` means the costim arm can only work in **trans** (engaging a costim receptor on a *different* cell, or acting cell-autonomously); `p_cis → 1` means the construct delivers **coincident signal 1 + signal 2 to one T cell** — the entire design intent of a CD3×TAA×costim trispecific.
- **(c) MECHANISTIC RATIONALE.** A Gaussian in span, not a hard cutoff, because both the linker and the membrane-proximal receptor stalks are flexible: a mismatch of a few nm is absorbed by conformational slack, a mismatch of many nm is not. `tol_nm = 8.0` sets how much slack. Rejected alternative: a step function (would make `p_cis` discontinuous in the design variable and destroy any gradient a format sweep could follow). **This is the model's cis/trans switch, and it is a *design* switch, not a fitted one.**
- **(d) UNITS.** `s`, `gap_match_nm`, `tol_nm` [nm]; `p_cis` dimensionless [0,1].
- **RUN-VERIFY (executed):** `p_cis(12.5) = 1.0` (exact match to gap) · `p_cis(8) = 0.8537` · `p_cis(20) = 0.6444` · `p_cis(28.5) = 0.1353` · `p_cis(60) = 2.21e-8`.
- **⚠ UNITS DISCREPANCY — flagged, not resolved.** The docstring (`multiarm_binding.py:42`) motivates the trans case as *"a deliberate epitope-HEIGHT mismatch (e.g. tall CRD1 4-1BB, **~60A**)"* — 60 **ångström = 6.0 nm**. But the format that implements it, `tetravalent_Cterm_trans` (`multiarm_binding.py:118`), sets `span_coeng_T_nm = **60.0**`, which the function reads as **60 nm**. Executed: `p_cis(60.0 nm) = 2.2e-8` (fully trans) whereas **`p_cis(6.0 nm) = 0.719`** (substantially *cis*). The "TRANS" designation of that format therefore rests on a 10× unit slip between the comment and the value, **or** on an unstated (and un-derived) mapping from epitope height to required co-engagement span. Either way the 60.0 value is **[ASSUMED — chosen to force p_cis→0; the ångström figure in the comment does not derive it]**.

---

### EQ-3 — Cross-cell bridge feasibility — **two inconsistent versions live in the codebase**

**EQ-3a — this module's version** (`multiarm_binding.py:48–53`):
```python
cleft = clip(span_bridge_nm, cleft_min=13.0, cleft_max=40.0)
g     = span_bridge_nm / max(cleft, 1e-6)
f     = clip(g, 0.0, 1.0)                       # ⇒ f = min(span/13, 1)
```

**EQ-3b — the version that actually runs** (`kinetic_synapse.py:67–72`, called at `kinetic_synapse.py:129`):
```python
g = span_bridge_nm / max(cleft_nm, 1e-6)
return clip((g − 0.6)/0.4, 0.0, 1.0)            # "<0.6x cannot bridge, >=1.0x full"
```

$$
f_{3a}(s)=\min\!\left(\frac{s}{13},\,1\right)
\qquad\qquad
f_{3b}(s)=\mathrm{clip}\!\left(\frac{s/\mathrm{clip}(s,13,40) - 0.6}{0.4},\,0,\,1\right)
=\begin{cases}0, & s \le 7.8\\[2pt] \dfrac{s-7.8}{5.2}, & 7.8<s<13\\[4pt] 1, & s\ge 13\end{cases}
$$

- **(b) BIOLOGICAL MEANING.** Two apposed membranes — T cell and tumor cell — are held a certain distance apart (the **immune-synapse cleft**, floored at the TCR–pMHC dimension, `CLEFT_MIN_NM = 13.0 nm`, `kinetic_synapse.py:41`). A drug that bridges them must physically span that gap. Too short → the two arms cannot simultaneously touch two membranes that cannot come closer together → **no trimer, no kill**. Long enough → the bridge forms; being *even longer* buys nothing, because reaching further across a gap you already crossed is not a benefit. The cleft is **emergent** in the sense that it *relaxes toward the bound-complex size*: `cleft = clip(span, 13, 40)` (`kinetic_synapse.py:128`) — the membranes settle at whatever separation the drug imposes, within a physical window.
- **(c) MECHANISTIC RATIONALE.** Feasibility is a **saturating, bounded gate**: it is a probability, it cannot exceed 1. This boundedness is the whole reason the optimum in EQ-9 exists. The `(g−0.6)/0.4` rescale in EQ-3b encodes a *soft* lower bound: below 0.6× the cleft, bridging is declared impossible; between 0.6× and 1.0×, partial (membrane ruffling / receptor tilt can close part of the gap). EQ-3a omits that rescale entirely.
- **(d) UNITS.** `s`, `cleft_*` [nm]; `f` dimensionless [0,1].
- **🚩 DISCREPANCY — MATERIAL, must be flagged to a reviewer.** `multiarm_binding.py:50` explicitly claims *"Same convention as kinetic_synapse.cleft_feasibility (span/cleft ramps to 1)"*. **It is not the same convention.** Executed this task at the `DART_Fc` span (`span_bridge = 6.5 nm`, `multiarm_binding.py:107`):
  - `multiarm_binding._bridge_feasibility(6.5) = **0.5**` → bridge forms at half efficiency (`KD_TAA` inflated 0.3 → 0.6 nM, EQ-6).
  - `kinetic_synapse.cleft_feasibility(6.5, cleft=13.0) = **0.0**` → **no bridge at all → zero killing.**
  Two live modules, one shipped format, **opposite predictions** (a partially active drug vs a completely inert one). The live clinical engine uses EQ-3b, so the *live* verdict on any construct with `span_bridge < 7.8 nm` is **total loss of efficacy**. See §5(4).

---

### EQ-4 — Monovalent arm binding (`multiarm_binding.py:60`)

```python
if n_arm == 1: return R_local·C/(KD + C + 1e-12)
```
$$B_{\text{mono}} = R_{\text{local}}\cdot\frac{C}{K_D + C}$$

- **(b) BIOLOGICAL MEANING.** Textbook 1:1 Langmuir occupancy: a single-arm binder on a cell with local receptor capacity `R_local` (nM, on the synapse reaction-volume basis) at free drug `C`. **No geometry enters.** A monovalent arm has nothing to co-engage, so no span, no reach, no `c_eff`.
- **(c) MECHANISTIC RATIONALE.** Equilibrium (QSS) rather than kinetic, because the arm-level binding equilibrates fast against the PD timestep; the *kinetic* treatment is reserved for the cross-cell bridge, where dwell time is the entire point (`kinetic_synapse.py`, T-subsystem for synapse kinetics). The `1e-12` is a divide-by-zero floor, not physics.
- **(d) UNITS.** `R_local`, `C`, `KD` [nM]; result [nM].
- **`n_arm ≤ 0` → returns 0** (`multiarm_binding.py:59`): the arm is **absent** from the construct. Valency 0 is how the format library says "this molecule has no costim arm" — it is not a zero-affinity arm, it is *no arm*.

---

### EQ-5 — Bivalent arm binding: the avidity quadratic (`multiarm_binding.py:61–66`)

```python
span  = coeng_span_nm if coeng_span_nm is not None else DEFAULT_ARM_REACH_NM   # 12.5 nm
ageff = geo_ageff_nM(rec_pc, r_cell_um, span)                                  # EQ-1
a = 2·ageff·C/KD² ;  b = 1 + 2·C/KD
x = (−b + sqrt(b² + 4a)) / (2a)          # positive root; fallback 1/b when a→0
x = clip(x, 0, 1)
return R_local·( 2·x·C/KD  +  ageff·x²·C/KD² )
```

The root solves the **receptor conservation identity** (verified algebraically and numerically this task):

$$
\underbrace{x}_{\text{free}} \;+\; \underbrace{\tfrac{2xC}{K_D}}_{\text{singly attached}} \;+\; \underbrace{2\cdot\tfrac{a_{\text{eff}}x^{2}C}{K_D^{2}}}_{\text{doubly attached (2 receptors each)}} \;=\; 1
\qquad\Longleftrightarrow\qquad a x^{2} + b x - 1 = 0
$$

and the **returned quantity is bound *drug*** (complexes), not bound receptor:

$$
B_{\text{bival}} = R_{\text{local}}\left(\underbrace{\frac{2xC}{K_D}}_{\text{1-armed complexes}} + \underbrace{\frac{a_{\text{eff}}\,x^{2}C}{K_D^{2}}}_{\text{2-armed complexes}}\right)
$$

- **(b) BIOLOGICAL MEANING.** A two-armed binder on one cell surface. The **factor 2** on the singly-bound term is the two-fold statistical degeneracy of *which* arm binds first. The second term is the avidity term: having landed, the free arm faces not bulk `C` but the tethered concentration `a_eff` (EQ-1) — hence `C/KD × a_eff/KD`. `x` is the free-receptor fraction, depressed by the fact that each doubly-engaged drug consumes **two** receptors (the factor 2 in the conservation identity).
- **(c) MECHANISTIC RATIONALE.** Solving the quadratic rather than assuming `x≈1` is what makes avidity **self-limiting**: at high `a_eff` the receptors get consumed two-at-a-time, `x` collapses, and the avidity term saturates. A linear ("infinite receptor") avidity model would let `a_eff/KD ≈ 10⁵` run away. **RUN-VERIFY (executed, `KD=0.3 nM`, `rec=1e5`, `r_c=8 µm`, `R_local=1`):**

  | `C` (nM) | mono `B` | bival `B` @ s=6.5 nm | bival `B` @ s=60 nm | avidity gain (×mono) 6.5 → 60 nm | receptor occupancy 6.5 → 60 nm |
  |---|---|---|---|---|---|
  | 0.003 | 0.00990 | 0.4914 | 0.4743 | **49.63× → 47.90×** | 0.9824 → 0.9476 |
  | 0.03  | 0.0909  | 0.4978 | 0.4933 | 5.475× → 5.426× | 0.9944 → 0.9831 |
  | 3.0   | 0.9091  | 0.5053 | 0.5159 | 0.556× → 0.567× | 0.9994 → 0.9983 |

  Two things to read off this table. **(i)** Sub-saturating (`C ≪ KD`), bivalency buys ~**50×** — and the gain **falls** as the span lengthens, exactly as EQ-1 demands. **(ii)** At saturating `C`, bivalent *bound-drug* is **lower** than monovalent (0.505 vs 0.909) and *rises* with span. That is **not a bug**: at saturation the receptor pool is the limiting resource, and a doubly-bound drug occupies two receptors, so fewer drug molecules fit. Receptor **occupancy** (right column) still behaves correctly — higher at short span. **A reader comparing `Cb` values across valencies must know they are counting molecules, not sites.**
- **(d) UNITS.** `a_eff`, `C`, `KD`, `R_local` [nM]; `a`, `b`, `x` dimensionless; result [nM].
- **Numerical guards** (`multiarm_binding.py:64`): `a > 1e-30` branch, `2a + 1e-300` denominator floor, `max(b, 1e-30)` in the monovalent fallback. All [CODE-INTERNAL].

---

### EQ-6 — Cross-cell bridge occupancy (`multiarm_binding.py:82–83`)

```python
f_bridge = _bridge_feasibility(span_bridge_nm)                 # EQ-3a
occ_b    = Cfree / (KD_TAA/max(f_bridge, 1e-3) + Cfree + 1e-12)
```
$$occ_b = \frac{C}{\dfrac{K_{D,\text{TAA}}}{f_{\text{bridge}}} + C}\qquad\Longleftrightarrow\qquad K_D^{\text{app}} = \frac{K_{D,\text{TAA}}}{f_{\text{bridge}}}$$

- **(b) BIOLOGICAL MEANING.** Geometric infeasibility is expressed as an **apparent affinity penalty on the cross-cell bond**: a construct whose arms struggle to span the cleft binds the *trans* partner as though its `KD` were worse by `1/f_bridge`. The intra-cell (mono) bindings `Cb3`, `CbT` cross no cleft and are untouched — the penalty is applied *only* where a membrane–membrane gap must be crossed.
- **(c) MECHANISTIC RATIONALE.** An effective-KD widening (rather than a multiplicative haircut on the bound species) keeps the penalty **concentration-dependent and saturable**: a geometrically-poor construct can still be driven to full bridge occupancy by enough drug, which is the correct physics for a *strained but possible* bond. Rejected alternative: `Cb_kill × f_bridge`, which would impose a hard ceiling no dose could overcome.
- **(d) UNITS.** `C`, `KD_TAA` [nM]; `f_bridge`, `occ_b` dimensionless.
- **⚠ FAIL-OPEN FLOOR.** `max(f_bridge, 1e-3)` caps the penalty at **1000×** `KD`. An infeasible bridge (`f_bridge → 0`) therefore still yields `occ_b > 0` — the geometry gate can never fully close in this module. (`kinetic_synapse`'s EQ-3b gate *can* — it multiplies `c_eff` by a feasibility that reaches exactly 0.) Tag: **[ASSUMED — numerical floor with behavioural consequences]**.
- **⚠ ASYMMETRY.** `occ_b` is built from `KD_TAA` and is applied to **both** `Cb_kill` (CD3×TAA) **and** `Cb_costTAA` (costim×TAA). `KD_CD3` and `KD_cos` never enter the bridge occupancy. Defensible (the TAA arm anchors the bridge on the tumor side) but it means **CD3 affinity has no effect on bridge feasibility in this module** — a design axis the counterscreen presumably wants to sweep.

---

### EQ-7 — Bridged species: the min() bottleneck (`multiarm_binding.py:84–85`)

```python
Cb_kill    = min(Cb3, CbT)·occ_b   if (n_CD3>0 and n_TAA>0) else 0
Cb_costTAA = min(CbC, CbT)·occ_b   if (n_cos>0 and n_TAA>0) else 0
```
$$C_{b,\text{kill}} = \min\!\left(C_{b,3},\,C_{b,T}\right)\cdot occ_b$$

- **(b) BIOLOGICAL MEANING.** A cross-cell trimer needs *both* ends. The number of bridges cannot exceed the scarcer of the two singly-bound pools — a low-copy TAA (BCMA) is bridge-limiting no matter how much CD3 is around, and vice-versa. This is the **abundance grading** that makes low-antigen tumors harder to kill.
- **(c) MECHANISTIC RATIONALE.** `min()` is a **hard bottleneck heuristic**, not a mass-action ternary equilibrium. The rigorous alternative (a Schropp/Douglass ternary closed form, which the live QSS path *does* use — `wholebody_pd.py:100–104`, Eq. 26–33) solves the coupled ternary explicitly. `min()` is chosen here for a format *screen*, where the ranking, not the absolute magnitude, is the deliverable.
- **⚠ DOUBLE-COUNTING OF `C`.** `Cb3` and `CbT` already contain `C/(KD+C)` saturation; multiplying by `occ_b`, which is *another* `C/(KD+C)` factor, makes `Cb_kill` scale as roughly `C²` at low dose instead of `C`. **This is a known structural approximation of the screen layer** and is one reason it is not the clinical-validation engine. Tag **[ASSUMED — screening closure]**.
- **(d) UNITS.** all [nM] × dimensionless → [nM].

---

### EQ-8 — Cis-coincident species (`multiarm_binding.py:87–91`)

```python
if n_CD3>0 and n_cos>0:
    p_cis  = _cis_feasibility(span_coeng_T_nm, gap_match_nm=12.5, tol_nm=8.0)
    Cb_cis = min(Cb3, CbC)·p_cis
else:
    Cb_cis = 0 ; p_cis = 0.0
```
$$C_{b,\text{cis}} = \min\!\left(C_{b,3},\,C_{b,C}\right)\cdot p_{\text{cis}}$$

- **(b) BIOLOGICAL MEANING.** The subset of engager molecules that are gripping **CD3 and costim on the same T cell simultaneously** — i.e. delivering signal 1 and signal 2 *coincidently, to the cell that is doing the killing*. This is the mechanistic definition of a "cis-costim" trispecific and the reason such a molecule should be safer than a free costim agonist: costim is only delivered where CD3 is engaged.
- **(c) MECHANISTIC RATIONALE.** Same `min()` bottleneck as EQ-7, now scaled by the *geometric probability* that one molecule can reach both. Both arms are on the T-cell side, so the relevant span is `span_coeng_T_nm` — **not** `span_bridge_nm`. In the live engine this same idea is implemented as an **occupancy gate** rather than a species: `occ_eff = occ·[(1−p_cis) + p_cis·f_cd3]` (`wholebody_pd.py:295`), so `p_cis=0` reduces to cell-autonomous costim and `p_cis=1` makes costim strictly conditional on that T cell's own CD3 engagement.
- **(d) UNITS.** [nM] × dimensionless → [nM]; `p_cis` dimensionless.

---

### EQ-9 — **THE CENTRAL RESULT: bounded gate × unbounded dilution ⇒ analytic optimum at span = cleft**

The live composite (`kinetic_synapse.py:128–131`):
```python
self.cleft_nm   = clip(span_bridge_nm, CLEFT_MIN_NM=13.0, CLEFT_MAX_NM=40.0)
self.feas       = cleft_feasibility(span_bridge_nm, self.cleft_nm)      # EQ-3b
self.ceff_trans = ageff_nM(dens, R_CELL_UM=8.0, span_bridge_nm) * self.feas
```

$$
\boxed{\;c_{\text{eff,trans}}(s)\;=\;\underbrace{\frac{A}{s}}_{\text{EQ-1: dilution, unbounded below}}\;\times\;\underbrace{f_{3b}(s)}_{\text{EQ-3b: feasibility, bounded above by 1}}\;}
\qquad A \equiv 2.4908086\cdot\frac{rec\cdot 10^{3}}{4\pi r_c^{2}}\;[\text{nM}\cdot\text{nm}]
$$

**Piecewise solution (exact, no numerics needed):**

$$
c_{\text{eff,trans}}(s)=
\begin{cases}
0 & s \le 7.8\ \text{nm} \quad (=0.6\times \text{cleft — cannot bridge})\\[6pt]
\dfrac{A}{s}\cdot\dfrac{s-7.8}{5.2} \;=\; \dfrac{A}{5.2}\left(1-\dfrac{7.8}{s}\right) & 7.8 < s < 13\ \text{nm} \quad \textbf{strictly increasing}\\[10pt]
\dfrac{A}{s} & s \ge 13\ \text{nm} \quad \textbf{strictly decreasing}
\end{cases}
$$

Both branches meet continuously at `s = 13` with value `A/13`. The function rises on `(7.8, 13)` and falls as `1/s` on `[13, ∞)`.

$$\Longrightarrow\qquad \boxed{\;s^{*} = \text{CLEFT\_MIN\_NM} = 13.0\ \text{nm}\;}\qquad c_{\text{eff,trans}}^{\max} = \frac{A}{13}$$

**RUN-VERIFY (executed this task; `rec = 257,000`, `r_c = 8 µm`; grid `s ∈ [1, 80]` nm at 1e-3 nm resolution):**

| span `s` (nm) | `ageff(s)` (nM) | `f_3b(s)` | `c_eff,trans` (nM) | fraction of max |
|---|---|---|---|---|
| 5.0 | 159,189.2 | 0.0000 | **0** | 0 |
| 7.8 | 102,044.4 | 0.0000 | **0** | 0 |
| 8.0 | 99,493.3 | 0.0385 | 3,826.7 | 0.0625 |
| 10.0 | 79,594.6 | 0.4231 | 33,674.6 | 0.550 |
| **12.5** *(LIVE default)* | 63,675.7 | **0.9038** | **57,553.0** | **0.940** |
| **13.0** *(= CLEFT_MIN)* | 61,226.6 | **1.0000** | **61,226.6** | **1.000 ← MAXIMUM**|
| 14.0 | 56,853.3 | 1.0000 | 56,853.3 | 0.929 |
| 20.0 | 39,797.3 | 1.0000 | 39,797.3 | 0.650 |
| 40.0 | 19,898.7 | 1.0000 | 19,898.7 | 0.325 |
| 60.0 | 13,265.8 | 1.0000 | 13,265.8 | 0.217 |
| 80.0 | 9,949.3 | 1.0000 | 9,949.3 | 0.162 |

Numerical argmax over the 79,001-point grid: **`s* = 13.000 nm`, `max = 61,226.62 nM`** — the analytic prediction, to the grid resolution. Above the optimum the decay is *exactly* `13/s` (check: `20 → 0.650`, `40 → 0.325`, `60 → 0.2167`, `80 → 0.1625` — all equal `13/s` to 4 decimals).

#### (b) BIOLOGICAL MEANING — what is physically happening

A drug that must bridge two cells has to do **two things at once**, and they pull in opposite directions:

1. **Reach across the cleft.** The membranes cannot come closer than ~13 nm (TCR–pMHC dimension). An arm pair that spans less than ~7.8 nm simply cannot touch both membranes → the trimer never forms → **no killing at all**. Getting longer helps — *until you can reach*.
2. **Find the partner once you're anchored.** The moment the first arm binds, the second arm stops being a solute and becomes a **tethered searcher**. Its search volume is `(2/3)πs³`, its target is a 2-D sheet of receptors covering `πs²`. Getting longer **hurts** — the partner it needs is diluted as `1/s`.

Requirement (1) is **satisfiable and then irrelevant** — once you can cross the gap, crossing it "more" is meaningless; feasibility is a probability capped at 1.0. Requirement (2) is **never satisfied** — every additional nanometre costs concentration, forever, with no floor.

**A benefit that saturates, multiplied by a cost that does not, is maximised exactly at the point where the benefit saturates.** That point is the cleft: **13.0 nm**.

#### (c) MECHANISTIC RATIONALE — why this is a real prediction and not an artefact

- **The result is structural, not tuned.** It follows from two independently-motivated pieces (a `1/s` Rhoden-type tethered-concentration law; a `[0,1]` steric-feasibility gate). Neither was written to produce an optimum. The optimum is **emergent from their product**, and it lands on a parameter (`CLEFT_MIN_NM`) that was set from an entirely different consideration (the TCR–pMHC dimension). Nothing was fitted to make `s* = 13`.
- **It inverts a common design intuition.** The engineering reflex — "add a longer, more flexible linker so the arms can find their targets" — is **wrong in this model, and wrong for a physical reason**: flexibility buys reach in three dimensions but the target lives in two. Every nm of extra linker triples the haystack and only doubles the needles.
- **It is consistent with the observed clinical format landscape.** Compact bridging formats (BiTE ~55 kDa, DART, tandem-scFv) and Fab–Fab spans cluster around ~13 nm, and the module's own `FORMATS` library independently assigns `span_bridge_nm = 13.0` to BiTE, IgG_1x1, IgG_2TAA_1CD3 and both trispecifics (`multiarm_binding.py:103–118`) — i.e. **the five formats the model treats as functional all sit exactly at the analytic optimum**, and the one that does not (DART_Fc, 6.5 nm) is the one the live gate would declare non-bridging (§5(4)). *Caveat:* this is consistency, **not validation** — the spans were assigned by architecture estimate, not measured, so this is not an independent test.
- **The clinically-canonical run sits 6% off the peak, on the steep side.** `span_bridge = 12.5 nm` (`pd_model_config.py:38`) gives `c_eff,trans = 0.940 × max` (executed). Moving the span 12.5 → 13.0 nm (+4% in length) **raises the second-arm concentration by +6.4%**; moving it 12.5 → 10.0 nm **drops it by 41%**; moving it to 7.8 nm **abolishes killing entirely**. The live model is therefore sitting on the *rising* limb, ~1 nm from a cliff — the single most sensitive geometric parameter in the engine.
- **Rejected alternative formulation:** making the feasibility gate *unbounded* (e.g. `f ∝ s`, "longer reaches better") would remove the optimum and make bridging span-independent (`(A/s)·(s/13) = A/13` = constant). The model deliberately does not do this, because a bond that already spans the gap gains nothing from further length — and the code says so in words at `kinetic_synapse.py:69–70` ("Full reach (span>=cleft) -> 1").

#### (d) UNITS
`s`, `CLEFT_MIN_NM`, `CLEFT_MAX_NM` [nm]; `A` [nM·nm]; `c_eff,trans` [nM]; `f_3b` dimensionless. `c_eff,trans` then enters the live bridge-formation rate as `kf = kon_TAA · (c_eff,trans · alive_frac)` (`kinetic_synapse.py:166–168`).

---

### EQ-10 — Format dispatch (`multiarm_binding.py:122–128`)

```python
if fmt not in FORMATS: raise KeyError(...)
g = {k:v for k,v in FORMATS[fmt].items() if k != 'mw_kda'}
g.update(kd_rec_over)
return multiarm_bound(Cfree, R_CD3, R_TAA, R_cos, **g)
```
- **Statement.** Named format → (valencies `n_*`, three spans) → `multiarm_bound`. Caller-supplied `KD_*`/`rec_*` override the defaults (target-specific).
- **`mw_kda` is stripped before the call (`:126`) — it is metadata and has NO effect on any binding computation in this module.** (Molecular weight matters for PK/extravasation, which lives elsewhere.) Documented explicitly so nobody assumes format MW is influencing the binding result.
- **Tag:** [CODE-INTERNAL] dispatch.

---

## 3. PARAMETERS OWNED

Provenance tags per the model's convention: **[MEASURED: source]**, **[DERIVED: from what]**, **[FITTED: fit to what]**, **[ASSUMED: rationale]**, **[UNSOURCED — TBD]**, plus **[UNVERIFIED CITATION]** where an in-code comment names a source I could not confirm from the code. Two tags added this audit: **[DEFINED — SI exact]** (a defined physical constant, which is *not* a measurement) and **[SOURCED-IN-REPO, UNSOURCED-IN-FILE]** (the identifier is absent here but present verbatim elsewhere in the live codebase, cited to file:line).

> **No [MEASURED] tag survives in this subsystem.** Every number T4 owns is [ASSUMED], [UNSOURCED], [CODE-INTERNAL], or [DEFINED]. That is the honest headline of this section: **T4 has no measured parameters.**

### 3.1 Module-level constants (`multiarm_binding.py`)

| Symbol | Value | Units | file:line | Provenance | Source / rationale |
|---|---|---|---|---|---|
| `AVO` | 6.02214076e23 | mol⁻¹ | :24 | **[DEFINED — SI exact]** | Avogadro constant. **Not [MEASURED]**: since the 2019 SI redefinition this is an *exact defined* value, not the result of a measurement with an uncertainty. Tagged separately so the provenance vocabulary is not diluted. |
| `DEFAULT_ARM_REACH_NM` | 12.5 | nm | :27 | **[UNSOURCED — TBD]** | Comment says only *"arm reach when a span is unset (compact within-module reach, nm)"*. **No citation in this file.** It equals `kinetic_synapse.SPAN_BRIDGE_DEFAULT_NM` (=12.5, `kinetic_synapse.py:39`) and `pd_model_config.span_bridge_nm` (=12.5, `pd_model_config.py:38`), so the three are consistent — but **consistency is not sourcing.** I grepped the repo: **no file attaches any citation to the number 12.5** (the Rhoden reference at `qsp_costim_window_v2.py:704` is attached to the *binding law*, not to this reach value). So 12.5 nm is an **uncited default replicated in three places**, which is a single unsourced assumption wearing three hats. |
| `_cc` | 1e9 | — | :32 | **[CODE-INTERNAL]** | Scratch unit factor; multiplied at :35 and divided at :36 — **cancels exactly**, no numerical effect. |

### 3.2 `_cis_feasibility` parameters (`multiarm_binding.py:39`)

| Symbol | Value | Units | file:line | Provenance | Source / rationale |
|---|---|---|---|---|---|
| `gap_match_nm` | 12.5 | nm | :39 (default), :88 (call site) | **[ASSUMED — architectural]** | The inter-epitope gap at which one molecule can co-grip CD3 + costim on one T cell. No citation in code. Same 12.5 nm as `DEFAULT_ARM_REACH_NM` — i.e. the model assumes the compact-construct reach *is* the height-matched gap. Circular but self-consistent. |
| `tol_nm` | 8.0 | nm | :39 (default), :88 | **[ASSUMED — flexibility slack]** | Gaussian width. Sets how much span/gap mismatch a flexible construct absorbs. No citation. **Load-bearing**: it alone decides whether a given epitope-height mismatch reads as cis or trans (`p_cis(20 nm)=0.64` at `tol=8` vs `0.044` at `tol=3`). |

### 3.3 `_bridge_feasibility` parameters (`multiarm_binding.py:48`) — **and their live twins**

| Symbol | Value | Units | file:line | Provenance | Source / rationale |
|---|---|---|---|---|---|
| `cleft_min` | 13.0 | nm | `multiarm_binding.py:48` | **[ASSUMED — structural dimension]** | Comment (`:49`): *"the T-arm<->tumor-arm span must hold the ~13-40nm immune synapse cleft"*. **No citation.** Its live twin `CLEFT_MIN_NM = 13.0` (`kinetic_synapse.py:41`) is commented *"immune-synapse cleft floor (~TCR-pMHC dimension)"* — again **no PMID**. The TCR–pMHC intermembrane dimension being ~13–15 nm is a standard structural-immunology figure, but **the code cites nobody**, so per the provenance rule this is **[ASSUMED]**, not [MEASURED]. **This is the single most consequential unsourced number in T4: it *is* the analytic optimum** (EQ-9). |
| `cleft_max` | 40.0 | nm | `multiarm_binding.py:48` | **[ASSUMED]** / **partly INERT** | Twin: `CLEFT_MAX_NM = 40.0` (`kinetic_synapse.py:42`), commented *"cleft ceiling before bond is mechanically unfavourable"*. **It imposes no mechanical penalty.** It only clips the cleft used in the ratio `g = s/cleft`; since feasibility saturates at 1 for any `s ≥ cleft`, a 60 nm or 80 nm span still returns `f = 1.0` (executed: `f_3b(60) = f_3b(80) = 1.0`). Long constructs are penalised **only** by `1/s` dilution — **the "mechanically unfavourable" claim in the comment is not implemented.** |
| `0.6 / 0.4` rescale | — | — | `kinetic_synapse.py:72` | **[ASSUMED — soft steric cutoff]** | *"<0.6x cannot bridge, >=1.0x full"*. No citation. Sets the hard-zero cutoff at **0.6 × 13.0 = 7.8 nm**. **Absent from `multiarm_binding._bridge_feasibility` — see EQ-3 discrepancy.** |
| `1e-3` feasibility floor | 1e-3 | — | `multiarm_binding.py:83` | **[ASSUMED — numerical, behavioural]** | Caps the apparent-KD penalty at 1000×; prevents the geometry gate from ever fully closing in this module. |

### 3.4 `multiarm_bound` signature defaults (`multiarm_binding.py:68–73`)

These are **function defaults**, overridable by every caller. They are **not** the clinical operating point.

| Symbol | Value | Units | file:line | Provenance | Source / rationale |
|---|---|---|---|---|---|
| `n_CD3` | 1 | valency | :69 | [CODE-INTERNAL default] | 1+1 baseline. |
| `n_TAA` | 1 | valency | :69 | [CODE-INTERNAL default] | |
| `n_cos` | 0 | valency | :69 | [CODE-INTERNAL default] | costim arm **absent** by default → `Cb_costim = Cb_cis = 0`. |
| `KD_CD3` | 3.0 | nM | :70 | **[UNSOURCED — TBD]** | **⚠ Inconsistent with the live engine**, which uses `KD_CD3_nM = 40.0` (`wholebody_pd.py:111`). A 13× difference. Whichever is right, they disagree. |
| `KD_TAA` | 0.3 | nM | :70 | **[UNSOURCED — TBD]** | **⚠ Live engine uses `KD_TAA_nM = 1.45`** (`wholebody_pd.py:111`). ~5× difference. |
| `KD_cos` | 1.0 | nM | :70 | **[UNSOURCED — TBD]** | Matches `wholebody_pd.KD_costim_nM = 1.0` (`wholebody_pd.py:112`) — also unsourced there. |
| `rec_CD3` | 1e4 | copies/cell | :71 | **[ASSUMED — placeholder]** | Overridden per-cell by real IHC/scRNA-derived copies in any live use. |
| `rec_TAA` | 1e5 | copies/cell | :71 | **[ASSUMED — placeholder]** | cf. the validated tumor anchor 257,000 CEACAM5 copies (`kinetic_synapse.py:36–37`). |
| `rec_cos` | 1e4 | copies/cell | :71 | **[ASSUMED — placeholder]** | |
| `r_cell_T` | 4.0 | µm | :72 | **[ASSUMED]** | T-cell radius. No citation here. (`wholebody_pd.py:126` independently uses `r_Tcell = 3.5 µm` with PMIDs 9400735 / 30571054 for the *myeloid contact* radius — **a different number in a different module**; I did not verify either against those PMIDs.) |
| `r_cell_tum` | 8.0 | µm | :72 | **[ASSUMED]** | Target-cell radius; matches `kinetic_synapse.R_CELL_UM = 8.0` (`:38`), commented *"Rhoden default"*. The Rhoden reference exists in-repo (`qsp_costim_window_v2.py:704`) but **is not attached to this radius** — no file states that 8.0 µm came from it. → **[ASSUMED — attribution asserted, not sourced]**. |
| `span_bridge_nm` | 12.5 | nm | :73 | **[UNSOURCED — TBD]** | Same value as the live `pd_model_config.span_bridge_nm` (`:38`). |
| `span_coeng_T_nm` | `None` | nm | :73 | [CODE-INTERNAL] | `None` ⇒ that cell-side has **no co-engagement** ⇒ `p_cis = 0` (`:44`) and bivalent arms fall back to `DEFAULT_ARM_REACH_NM`. |
| `span_coeng_tumor_nm` | `None` | nm | :73 | [CODE-INTERNAL] | as above. |

### 3.5 `FORMATS` library (`multiarm_binding.py:101–120`)

The header comment (`:97–98`) states the provenance in full: *"Spans are literature/architecture geometry estimates (Fab-Fab ~13nm; tandem-scFv flexible ~13nm; C-term fusion end ~5-6nm; Fc top-to-bottom ~10-15nm); AF3-derived spans OVERRIDE per construct."* **No PMID, DOI or structure ID appears.** Every span below is therefore **[ASSUMED — architecture estimate; AF3 override intended but not present in this file]**. `p_cis` and `f_bridge` columns were **computed by executing the code this task**.

| Format | `n_CD3` | `n_TAA` | `n_cos` | `span_bridge` (nm) | `span_coeng_T` (nm) | `span_coeng_tumor` (nm) | `mw_kda` | ⇒ `p_cis` | ⇒ `f_bridge` (EQ-3a) | ⇒ live `f_3b` (EQ-3b) | line |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `BiTE` | 1 | 1 | 0 | **13.0** | None | None | 55 | 0.0 | 1.0 | **1.0** | :103 |
| `IgG_1x1` | 1 | 1 | 0 | **13.0** | None | None | 150 | 0.0 | 1.0 | **1.0** | :105 |
| `DART_Fc` | 1 | 1 | 0 | **6.5** | None | None | 105 | 0.0 | **0.5** | **0.0 ⚠** | :107 |
| `IgG_2TAA_1CD3` | 1 | **2** | 0 | **13.0** | None | 13.0 | 195 | 0.0 | 1.0 | **1.0** | :110 |
| `tetravalent_Cterm_cis` | 1 | **2** | 1 | **13.0** | **12.5** | 13.0 | 200 | **1.000** | 1.0 | **1.0** | :114 |
| `tetravalent_Cterm_trans` | 1 | **2** | 1 | **13.0** | **60.0** | 13.0 | 200 | **2.21e-8** | 1.0 | **1.0** | :117 |

Notes on this table:
- **All five "working" formats sit at `span_bridge = 13.0 nm` = the analytic optimum** (EQ-9). That is a *design assignment*, not a fit — but it means the format library cannot, as written, demonstrate the optimum, because it never samples off it (except DART_Fc, which samples the infeasible side).
- **`DART_Fc` is the only format where EQ-3a and EQ-3b disagree — and they disagree maximally** (0.5 vs 0.0: "half-active" vs "completely inert"). Real DART molecules are clinically active bridgers; a live engine that returns zero kill for them is making a falsifiable — and probably false — prediction. See §5(4).
- **`mw_kda` is never used in any binding computation** (stripped at `:126`).

### 3.6 Parameters T4 *reads* but does not own

| Symbol | Value | Owner | Why T4 cares |
|---|---|---|---|
| `CLEFT_MIN_NM` | 13.0 nm | `kinetic_synapse.py:41` | **is the optimum** (EQ-9) |
| `CLEFT_MAX_NM` | 40.0 nm | `kinetic_synapse.py:42` | clips the emergent cleft; **imposes no penalty** |
| `R_CELL_UM` | 8.0 µm | `kinetic_synapse.py:38` | sets `A` in EQ-9 (`c_eff ∝ 1/r_c²`) |
| `NM_PER_COPY` | 6.0/257000 = 2.335e-5 nM/copy | `kinetic_synapse.py:37` | *"pinned by the validated tumor (CEACAM5 257,000 copies -> Rcap_TAA 6.0 nM)"* — **[DERIVED: from a stated 257,000-copy CEACAM5 anchor; the anchor itself carries no citation in code → UNVERIFIED]** |
| `span_bridge_nm` | 12.5 nm | `pd_model_config.py:38` | **the live operating point** — 0.940 × optimum |
| `span_cis_nm` | 12.5 nm | `pd_model_config.py:39` | feeds `KineticSynapse` cis avidity (`cis_avidity = 0.0` live → inert) |

> **UNSOURCED-NUMBER NOTICE (scope: not T4's, do not propagate as sourced).** `pd_model_config.py:22` and `:44` describe the calibration anchor as *"mosunetuzumab IL-6 **570** pg/mL"*. A grep of the repository this task finds **570 carries no citation anywhere** — no PMID, DOI, or reference; the only other appearance is `model/CHANGELOG_2026-07-13.md:30`, which merely *asserts* "clinical 570". So: **570 is [UNSOURCED — TBD] as far as the code can show.** That is the entire verifiable claim. **T4 owns no IL-6 parameter and this doc uses none**, so no IL-6 value — 570 or any replacement — is stated here. Anyone needing the correct anchor must source it from the clinical literature and record the identifier in `pd_model_config.py`; this doc deliberately proposes no substitute value, because it has no code-verifiable basis for one.

---

## 4. WHAT IS EMERGENT vs IMPOSED

### 4.1 Genuinely EMERGENT (computed from mechanism, not handed in)

| Quantity | Emerges from | Why it counts as emergence |
|---|---|---|
| **`c_eff(s) ∝ 1/s`** | receptor surface density + hemispherical reach volume (EQ-1) | Not a fitted curve. The exponent `−1` is *forced* by the 2-D-target/3-D-search mismatch. Nobody chose it; the geometry did. Verified exactly (`c_eff·s` invariant to machine precision). |
| **The span optimum `s* = 13.0 nm`** | product of EQ-1 (unbounded `1/s` decay) and EQ-3b (bounded `[0,1]` gate) | **The headline result.** It is not a parameter, not a fit, and nobody wrote "the best span is 13". It falls out of two independent mechanisms whose only shared symbol is `s`. It lands on `CLEFT_MIN_NM` — a constant introduced for an unrelated reason (membrane separation). |
| **"Longer linkers reduce bridging"** | EQ-9 above the optimum: `c_eff,trans = A/s` exactly | A counter-intuitive, falsifiable prediction (executed: `13/s` decay — a 60 nm construct retains **21.7%** of the optimal second-arm concentration; an 80 nm construct **16.2%**). |
| **Avidity self-limitation** | the conservation quadratic (EQ-5) | The `x` root, and therefore the saturation of avidity as receptors are consumed two-at-a-time, is *solved for*, not assumed. Executed: bivalent gain collapses from 49.6× at `C = 0.01·KD` to 0.56× at `C = 10·KD`. |
| **Abundance grading of bridging** | `min(Cb3, CbT)` (EQ-7) + per-cell receptor copies | A low-copy target becomes bridge-limiting without any per-target tuning constant. |
| **cis vs trans behaviour of a trispecific** | `p_cis` from span geometry × real per-cell CD3 engagement (`wholebody_pd.py:294–295`; the design intent is stated in the docstring at `:288`) | The costim/CD3 coincidence is a *geometric probability × a computed binding state*, not a phenomenological "coincidence factor". |

### 4.2 IMPOSED (constants handed to the subsystem)

| Quantity | Imposed value | Honest status |
|---|---|---|
| `CLEFT_MIN_NM = 13.0` | the location of the optimum | **[ASSUMED — no citation in code].** The *existence* of an optimum is emergent; its *location* is imposed by this one number. Move it to 15 nm and the optimum moves to 15 nm. **The emergence claim is about the shape of the curve, not the position of its peak.** |
| The 0.6× hard-zero cutoff | 7.8 nm | **[ASSUMED].** Sets which formats are declared dead. |
| The hemisphere reach envelope | uniform, radius `s` | **[ASSUMED].** A real linker samples a chain-statistics distribution, not a uniform ball. This is where the `1/s` law is most vulnerable (§5(1)). |
| `gap_match = 12.5`, `tol = 8.0` | the cis/trans switch | **[ASSUMED × 2].** `p_cis` is *computed*, but from two invented constants. |
| All `KD_*`, `rec_*`, `r_cell_*` defaults | see §3.4 | **[UNSOURCED / ASSUMED placeholders]**, two of which (`KD_CD3`, `KD_TAA`) **contradict the live engine's values**. |
| `FORMATS` spans | 6.5 / 13.0 / 12.5 / 60.0 nm | **[ASSUMED — architecture estimates].** The docstring promises *"AF3-derived spans OVERRIDE per construct"* — **no AF3-derived span is present in this file.** The AF3 hook is a hook, not a value. |
| `mw_kda` | 55–200 kDa | **[ASSUMED]**, and **unused** (stripped at `:126`). |

### 4.3 The honest summary sentence

> **T4 derives a non-trivial, counter-intuitive shape (bridging strength peaks at a finite span and decays as 1/span beyond it) from mechanism, and then pins that shape's peak to a single uncited constant.** The *physics* is emergent; the *number 13* is an assumption. A reviewer is entitled to accept the first and demand a citation for the second — and they would be right to.

---

## 5. KNOWN LIMITATIONS / OPEN QUESTIONS

1. **The reach envelope is a uniform hemisphere, not a polymer.** `V_Ab = (2/3)πs³` treats the tethered arm as being *equally likely to be anywhere* within radius `s`. A real flexible linker (or a Fab on a hinge) has an end-to-end distribution — worm-like-chain or Gaussian-chain — that concentrates probability at intermediate extensions and *penalises full extension entropically*. Consequences: (a) the true dilution exponent is not exactly `−1`; (b) a real construct pays an **additional** entropic cost for a long linker (it must stretch), which the model **does not charge** — so the model, if anything, **under-states** the penalty for long spans. The direction of the headline conclusion is therefore robust to this limitation, but its magnitude is not. **This is the most likely first attack from a biophysics reviewer.** *(The sibling ADC model in this lab implements a WLC/Shaw composite reach kernel; T4 does not.)*
2. **`CLEFT_MIN_NM = 13.0` carries no citation anywhere in the code.** The whole optimum rides on it. It should be sourced to a structural/immunology measurement of the T:target intermembrane distance (a TCR–pMHC-set cleft), and the sensitivity of every format ranking to it should be reported. **Until then, "the optimal span is 13 nm" must be stated as "the optimal span is the cleft dimension, which the model takes to be 13 nm."**
3. **The `0.6×` hard-zero cutoff (7.8 nm) is invented.** It determines, by itself, that any construct with `span_bridge < 7.8 nm` has **exactly zero** efficacy. That is a very strong claim resting on an uncited numerator.
4. **🚩 Two live modules implement the cleft gate differently, and it changes a shipped format's verdict.** `multiarm_binding._bridge_feasibility` (linear ramp, no 0.6 offset) vs `kinetic_synapse.cleft_feasibility` (0.6/0.4 rescale) — despite `multiarm_binding.py:50` asserting they are "the same convention." At the `DART_Fc` span of 6.5 nm they return **0.5 and 0.0**. Since DART-format engagers are clinically active bridgers, the live gate's verdict (zero kill) is probably **wrong**, which means either (a) the 6.5 nm span assigned to DART_Fc is wrong, or (b) the 0.6× cutoff is too aggressive, or (c) the cleft floor is too high. **This must be reconciled before any format counterscreen is run.** It does not currently corrupt the clinical validation run, because that run uses `span_bridge = 12.5 nm` and never touches `FORMATS`.
5. **`CLEFT_MAX_NM = 40 nm` is documented as a mechanical-unfavourability ceiling but implements no penalty.** Executed: `f_3b(80 nm) = 1.0`. A very long construct is penalised only by dilution. If long spans are additionally destabilising (membrane tension, bond-force geometry, poor mechanical coupling to TCR triggering — all plausible), the model **under-penalises** them.
6. **In the live engine the bridge may be geometrically saturated, blunting the equilibrium consequences of the span.** At the live point, `c_eff,trans / KD_TAA = 57,553 / 1.45 ≈ **3.97e4**` (executed). At ~40,000× KD, the *equilibrium* bridge occupancy is saturated, so the `1/s` law cannot move the equilibrium magnitude much — it moves the bridge-**formation rate** `kf = kon_TAA · c_eff,trans` (`kinetic_synapse.py:168`), which matters only because the synapse is treated kinetically (dwell time, serial killing, internalisation `kint` competing). **So: the span optimum is a *kinetic* prediction in the live model, not an equilibrium one.** A reviewer will ask for the sensitivity analysis (`kill_frac` and `IL-6` vs `span_bridge` over 8–40 nm) that turns this analytic claim into a model output. **That run has not been done and is the obvious next experiment.**
7. **`multiarm_bound` double-counts drug concentration** (EQ-7): `Cb3` and `CbT` are already saturating in `C`, and `occ_b` multiplies in *another* `C/(KD+C)`, so `Cb_kill ~ C²` at low dose. Acceptable for a *ranking* screen; **not** acceptable as a dose-response predictor. The live clinical path avoids this by using either the Schropp ternary QSS (`wholebody_pd.py:100–104`) or the kinetic bond ODE.
8. **`min()` is not mass action.** The bridged-species bottleneck is a heuristic; the true ternary complex has a prozone (hook effect) that `min(·)·occ_b` reproduces only qualitatively.
9. **The trispecific "trans" format's span (60.0 nm) is 10× the epitope height its own comment cites (~60 Å = 6 nm).** Executed: `p_cis(6.0 nm) = 0.719` — i.e. if the comment's number were used as the span, the "trans" construct would be **substantially cis**. Either the comment or the value is wrong, or an unstated height→span mapping is intended. **The entire cis-vs-trans contrast in the format library depends on which.**
10. **The bivalent return value is bound *drug*, not bound *receptor*.** At saturating dose, a bivalent arm returns a *smaller* number than a monovalent one (0.505 vs 0.909, executed) because each drug consumes two receptors. Anyone comparing `Cb` across valencies without reading EQ-5 will conclude that bivalency *hurts*. It does not; the metric changes. This is a **documentation-level trap**, and it is now documented.
11. **`p_cis` is live-imported but returns 0.0 in the canonical run** (`span_coeng_T_nm` is never set by `pd_model_config.KINETIC`). The entire cis-coincidence mechanism is built, wired, and **default-inert**. It is honest to say the model *can* represent cis costim delivery; it is **not** honest to say the validated clinical runs *exercise* it. They do not.
12. **No AF3-derived span exists in the code**, despite the module docstring promising per-construct AF3 overrides. Every span in `FORMATS` is an architecture estimate. Until real predicted (or crystallographic/SAXS) spans land, the format layer's *inputs* are the weakest link in an otherwise mechanistic chain — and no amount of rigor in EQ-9 can compensate for a wrong `s`.

---

## Appendix — reproduction commands (executed this task)

Environment: `/usr/local/Caskroom/miniconda/base/envs/claude-skills/bin/python`, cwd = `/Volumes/RAID4/costim_engager_counterscreen/model/engine`.

```python
import numpy as np, sys; sys.path.insert(0,'.')
from multiarm_binding import geo_ageff_nM, _cis_feasibility, _bridge_feasibility, bound_arm, FORMATS, AVO
from kinetic_synapse import cleft_feasibility, ageff_nM, CLEFT_MIN_NM, CLEFT_MAX_NM

# EQ-1 closed form vs code: agrees to 1 ULP (NOT "5 sig figs" — that claim was a typo'd constant)
K = 3e24/(2*AVO)                                     # -> 2.49080860076077
assert abs(K*257000/(4*np.pi*8.0**2*0.0125) - geo_ageff_nM(257000,8.0,12.5)) < 1e-9   # 63675.68549029321

# 1/span law is exact: c_eff * s is invariant (309706.641489753 at rec=1e5, r_c=8 um)
assert all(abs(geo_ageff_nM(1e5,8.0,s)*s - 309706.641489753) < 1e-6 for s in [5,6.5,10,13,20,30,40,60])

# geo_ageff_nM  ==  kinetic_synapse.ageff_nM   (same physics, two files; equal to 1 ULP)
a, b = geo_ageff_nM(257000,8.0,13.0), ageff_nM(257000,8.0,13.0)
assert abs(a-b)/b < 1e-15        # 61226.6206637434... ; rel_diff = 2.4e-16 (op-order, not physics)

# the composite optimum
P = lambda s: geo_ageff_nM(257000,8.0,s)*cleft_feasibility(s, float(np.clip(s,CLEFT_MIN_NM,CLEFT_MAX_NM)))
g = np.linspace(1,80,79001); v = np.array([P(s) for s in g])
print(g[v.argmax()], v.max())        # -> 13.0 nm, 61226.62 nM   == CLEFT_MIN_NM
print(P(12.5)/P(13.0))               # -> 0.9400   (live operating point = 94.0% of max)

# the EQ-3 discrepancy
print(_bridge_feasibility(6.5), cleft_feasibility(6.5, 13.0))   # -> 0.5   vs   0.0
```
