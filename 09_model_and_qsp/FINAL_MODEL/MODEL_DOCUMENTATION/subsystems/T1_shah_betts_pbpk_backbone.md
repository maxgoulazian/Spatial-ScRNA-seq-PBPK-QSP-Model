---
title: "T1 — Shah & Betts platform PBPK backbone (organ physiology, 2-pore transport, FcRn/renal clearance)"
subsystem: T1 (canonical; the ORIGIN layer — everything else is built on this)
model: costim_engager_counterscreen
source_file: engine/qsp_costim_window_v2.py
date: 2026-07-13
generated_by: workflow-subagent T1
scope: "LAYER 1 (lines 62-295, 342-363, 890-934) + the PK readouts (1262-1282). The PD/cytokine/window layers are NOT in scope (see T2/T3/T4/T5)."
verification: "All arrays and derived scalars in this doc were PRINTED FROM THE LIVE MODULE, and all validation numbers were MEASURED BY RUNNING IT (conda claude-skills, scipy LSODA + exact eigen-decomposition). Nothing here is quoted from a code comment without an independent check; where a comment could not be reproduced, that is stated."
adversarial_review: "2026-07-13 — re-verified line-by-line against the live module by an independent adversarial pass. 6 defects found and FIXED in place: (1) the claim that the Tabula Sapiens provenance CSV does not exist was FALSE (it exists and all 30 densities reproduce from it); (2) the Shah-Betts sigma_V ranges quoted in §3.1 were FABRICATED numbers (real: sigma_1 0.883-0.987, sigma_2 0.311-0.837); (3) BiTE t50 was stated as 1.80 h, measured 1.78 h; (4) four off-by-one line citations; (5) a dangling 'EQ-17' cross-reference to an equation that does not exist; (6) an internal 5-vs-6 dead-parameter count contradiction. See §7."
---

up:: [[00_INDEX]]
tags:: #model-doc #pbpk #shah-betts

> [!abstract]+ At a glance
> **T1 is the foundation layer of the entire model — the Shah & Betts (2012) platform PBPK reproduction that every other subsystem was built on top of.**
> - **Owns:** 15-organ physiology (V, Q, L, σ_V, V_v, V_is, portal), 2 immune-density arrays, the plasma/lymph pools, the 2-pore convective transport law, and the format-dependent (FcRn + renal) clearance.
> - **Equations:** EQ-1 … EQ-16 (13 in code, 3 derived here) · **Parameters:** 30 documented, of which **6 entries (9 symbols) are confirmed DEAD** (defined, never read — §3.6)
> - **Provenance:** every number tagged `[MEASURED]` / `[DERIVED]` / `[FITTED]` / `[ASSUMED]` / `[UNSOURCED — TBD]`.
> - **Reach:** these arrays are the single source of physiological truth for BOTH execution paths — the abstract window-scoring path AND the per-cell clinical-validation path (`run_tce_pd_reval.py:8` instantiates `q._PBPKArrays()` directly).
> - **Headline honest findings:** σ_L = 0.85 is FITTED and is **4.25× the platform value it claims lineage from** (Shah-Betts σ_L = 0.2); the "terminal t½ = 16.1 d exactly" claim is **dose-conditional**, not a model property; `(_PLASMA_CO, _LYMPH_RATIO, k_dist)` are **provably non-identifiable** in this module (only their product matters); 6 dead parameter entries; the σ_V array is the layer's most influential biodistribution input and is **unsourced in code**.
> - **Retracted finding (adversarial pass 2026-07-13):** an earlier revision of this doc claimed the Tabula Sapiens provenance CSV "does not exist in the repo". **That was false.** It exists at `deliverables/06_qsp_science-qsp/v2_enhanced/qsp_tabula_sapiens_densities.csv` and **all 30 immune densities reproduce from it exactly** (§3.4). The immune-density provenance is SOUND.

---

## 0. BUILD HISTORY — why this layer exists and what it is a reproduction OF

This document carries the build-history opening for the whole model, because this layer *is* the origin.

The project did not begin as a costim counter-screen. It began as a **reproduction of the Shah & Betts (2012) platform PBPK model** — the canonical whole-body monoclonal-antibody distribution model from the Balthasar-adjacent PK literature — and everything else in the codebase (trimer synapse, cytokine cascade, Treg suppression, therapeutic-window scoring) was subsequently grafted onto that backbone. The module's own header states the lineage explicitly:

> `# LAYER 1 — FULL-BODY PBPK (Shah & Betts 2012 platform; 2-pore + FcRn)` (`engine/qsp_costim_window_v2.py:63`)

**The reference being reproduced.** The PMID and DOI below **are present in the code** (`:893`); they were re-checked against PubMed/the publisher record in the 2026-07-13 adversarial pass and are correct:

- Shah DK, Betts AM. *Towards a platform PBPK model to characterize the plasma and tissue disposition of monoclonal antibodies in preclinical species and human.* **J Pharmacokinet Pharmacodyn 2012; 39(1):67–86. PMID 22143261. DOI 10.1007/s10928-011-9232-2.** ✅ **CITATION VERIFIED** — and, importantly, **the citation is in the code**, so it is not a doc-side invention.

Properties of the real Shah-Betts platform (checked against the published record in the 2026-07-13 pass; **note these are literature claims, NOT things the repo can prove**):
1. **15 tissues + a carcass + a tumor compartment.** ✅ confirmed.
2. **Lymph flow set at 0.2% of plasma flow** — i.e. `L = Q/500`. ✅ confirmed; a *genuine* platform relation, and the code reproduces it exactly (`_LYMPH_RATIO = 1.0/500.0`, line 101). *(An earlier revision presented this as a verbatim quotation from the paper; it is a paraphrase of the platform relation, and is stated as such here.)*
3. **Vascular reflection coefficients adjusted per tissue according to vascular pore size**, in a genuine TWO-pore split: tight-pore σ₁ ≈ **0.883–0.987**, leaky-pore σ₂ ≈ **0.311–0.837**. ⚠️ reproduced only *in spirit* (a single per-organ `sigV`, no small/large pore split) — see §3.1.
4. **Lymphatic capillary reflection coefficient σ_L = 0.2**, assigned (not fitted) in the platform. ❌ **NOT reproduced** — the code uses `sigL = 0.85` (see §3.3 and §5.1; this is the single largest deviation).
5. Each tissue subdivided into **five** spaces (plasma, blood-cell, endothelial, interstitial, cellular), with **FcRn binding described by explicit on/off rates**, recycled bound mAb, and first-order elimination of unbound. ❌ **NOT reproduced** — the code collapses this to **two** spaces (vascular, interstitial) and replaces the endosomal FcRn ODE with a quasi-steady-state *recycled-fraction* scalar. The code is candid about this (lines 69–72):
   > *"FcRn salvage is modelled as an endosomal recycled fraction (quasi-steady-state reduction of catabolism) rather than an explicit on/off endosome ODE, for numerical robustness while preserving the platform's long IgG half-life."*

So: **T1 is a *reduced* Shah-Betts, not a faithful one.** It keeps the platform's *topology* (organ interstitia fed by lymph-flow-limited 2-pore convection, returning through a lymph pool) and its signature `L = Q/500` relation, and discards the platform's 5-space endosomal machinery. A reviewer must be told this in one sentence; §5 does so.

**What was subsequently built ON this backbone** (each a separate doc): the emergent bivalent-binding TMDD sink (§2 EQ-13 here, mechanism doc elsewhere), the format-aware PD driver, the resolved CD8/CD4-conv/Treg/myeloid cell layer, the cytokine cascade, and the window-scoring metrics. T1 supplies all of them with their compartment volumes and their drug concentrations. **Nothing downstream can be more correct than this layer.**

---

## 1. PURPOSE & DATAFLOW POSITION

### What T1 does
T1 converts an **administered dose** (IV bolus or continuous infusion) into a **per-organ interstitial drug concentration time-course**, `Cis_i(t)`, plus a plasma concentration `C_pl(t)`. That is its entire job. Every downstream subsystem consumes one or both of those.

### Where it sits in the life of the molecule
```
   IV bolus (A_pl(0)=dose)  /  infusion (mg/day)
              |
              v
      [ PLASMA POOL  A_pl ]  <-- V_PLASMA = 3.1 L
        |            ^
        | J_extrav   | k_lymph_return * A_ly       <-- EQ-6, EQ-11
        | (2-pore    |
        |  convection)
        v            |
  [ 15 ORGAN INTERSTITIA  A_is_i ]  <-- V_is,i     <-- the site where TARGETS and CELLS live
        |            |
        | J_return   | tmdd_sink (internalisation of bound complex)
        v            v
   [ LYMPH POOL A_ly ] ---> back to plasma          (elimination)
```
- **Organ VASCULAR spaces are NOT integrated.** They equilibrate ~10⁵× faster than the antibody timescale (`Q_i/Vv_i` up to ~1e5/day vs pharmacology over days), so they are treated as quasi-steady-state and **lumped into the plasma pool** (rationale at lines 301–308, 890–902). This is the model's stiffness-removal move and it is physically defensible for perfusion-limited vascular spaces.
- Elimination is **plasma-side** (catabolic + renal, EQ-8/9/10) and **interstitium-side** (TMDD internalisation, EQ-13).

### What feeds T1
Nothing upstream — T1 is the head of the graph. Its inputs are the dose, the frozen physiology, and (optionally) a construct-format CSV that overrides `mw_kda`, `has_fc`, valency and `reach_gate` (`apply_construct_format`, line 1573).

### What T1 feeds — and this is the load-bearing claim of this document
T1's arrays are the **single source of physiological truth for BOTH execution paths in the live model**, which I verified by tracing the import graph rather than assuming it:

| Consumer | How it gets T1 | Evidence |
|---|---|---|
| **Abstract window-scoring path** (the costim nomination) | in-module: `rhs()` integrates the PK block directly | `qsp_costim_window_v2.py:890–934` |
| **Per-cell clinical-validation path** (`run_tce_pd_reval` → `coupled_percell_pd` → `coupled_percell_pk`) | **imports and instantiates T1's array class**, then hands Q, L, σ_V, V_is, V_v to the per-cell PK engine | `run_tce_pd_reval.py:6` `import qsp_costim_window_v2 as q` → `:8` `pb=q._PBPKArrays()` → `:51` `Q,L,sigV,Vis,Vv=arr(pb.Q),arr(pb.L),arr(pb.sigV),arr(pb.Vis),arr(pb.Vv)` → `:142–143` passed into `CoupledPerCellPD(...)` |

So a physiology error here propagates into *everything*, including the clinical PK/PD validation that is supposed to check the model. That is the strongest possible argument for auditing this layer hardest — and it is why this doc runs the numbers rather than trusting the comments.

### Live-code scope boundary
In scope: lines **62–295** (physiology + PBPK dataclass), **342–363** (state layout, plasma/lymph volumes), **890–934** (the PK block of the RHS), **1224–1282** (dose entry + PK readouts), **1542–1595** (the format/spatial loaders that mutate T1 state). Everything from line 936 (`L2 engagement`) onward is PD and belongs to other subsystem docs.

---

## 2. GOVERNING EQUATIONS

Notation: `i` indexes organ (0..14). Amounts `A` in model-amount units (mg when a real `dose_mg` is supplied — the linear-PK rescale + µg/mL and nM conversion live at `:1272–1279`, not in any equation below); volumes in L; flows in L/day; rate constants in /day; concentrations in amount/L.

---

### EQ-1 — Organ blood-flow assignment with loop-conserving normalisation (`qsp_costim_window_v2.py:138–140`)
```
mask_i   = (name_i != "lung")
Qf_i     = Qf_i / Σ_{j∈mask} Qf_j            for i ∈ mask      (line 139)
Q_i      = Qf_i · CO_plasma                  for i ∈ mask      (line 140)
Q_lung   = CO_plasma                                            (line 140)
```
- **Biological meaning.** The lung sits in **series** with the whole circulation (all cardiac output traverses it), while the other 14 organs sit in **parallel** and share that same output between them. Renormalising the parallel fractions to sum to 1 guarantees the circulatory loop conserves flow: Σ(parallel Q) ≡ Q_lung.
- **Mechanistic rationale / alternative rejected.** The raw `Qfrac` column in `_PBPK_TISSUES` does *not* sum to 1 (the literal fractions sum to ≈0.945 excluding lung), so using them unnormalised would silently leak ~5.5% of cardiac output. The normalisation is the fix. **RUN-VERIFIED:** Σ Q over the 14 parallel organs = **5000.000000000001 L/day**, Q_lung = **5000.0 L/day** — the loop closes exactly.
- **Units.** `Qf` dimensionless; `CO_plasma` L/day; `Q_i` L/day.
- ⚠️ **CRITICAL — `Q` is never used again inside this module.** Its *only* consumer here is EQ-2. Verified exhaustively: `grep "PB.Q"` over the live file returns **zero** hits outside the class constructor. The blood-flow topology the header advertises at lines 72–73 — *"Portal organs drain to liver; liver + all organs drain to a venous pool -> lung (series) -> arterial pool"* — **is not implemented in the live RHS.** See §5.2. (`Q` *does* acquire an independent role downstream in the per-cell path, `wholebody_percell.py:168`.)

---

### EQ-2 — Lymph flow (the signature Shah-Betts relation) (`:141`, constant at `:101`)
```
L_i = Q_i · (1/500)
```
- **Biological meaning.** Lymph is the return route for the protein-rich fluid that convects out of the vasculature. In the antibody-transport literature the lymphatic flow of a tissue is ~500-fold smaller than its plasma flow; because antibody extravasation is **convection-dominated and lymph-flow-limited** (not blood-flow-limited, unlike a small molecule), `L` — not `Q` — is the transport-relevant flow for an IgG.
- **Mechanistic rationale.** This single line is what makes the model an *antibody* PBPK rather than a small-molecule one. Rejecting the alternative (perfusion-limited, `Q`-driven distribution) is the correct call for a 150 kDa protein that cannot cross a continuous endothelium by diffusion.
- **Provenance.** ✅ **[MEASURED: Shah & Betts 2012, PMID 22143261 — PMID/DOI present in the code at `:893`]**. The platform sets lymph flow at 0.2% of plasma flow (= Q/500); confirmed against the published record in the 2026-07-13 adversarial pass. This is a real platform value, faithfully reproduced. *(Not a verbatim quotation — an earlier revision of this doc wrongly presented it as one.)*
- **Units.** L/day. **RUN-VERIFIED:** L_lung = 10.0, L_muscle = 1.7989, L_tumor = 0.2116 L/day.

---

### EQ-3 — Vascular and interstitial sub-volumes (`:142–143`)
```
Vv_i  = fV_i  · V_i        (organ vascular sub-volume)
Vis_i = fIS_i · V_i        (organ interstitial sub-volume)
```
- **Biological meaning.** An organ's total volume splits into a vascular fraction (plasma inside its capillaries) and an interstitial fraction (the extravascular, extracellular space where a therapeutic antibody actually meets its target). `Vis` is the **pharmacologically relevant volume** — it is the denominator that turns an extravasated amount into the concentration that drives binding, TMDD and PD.
- **Mechanistic rationale.** Antibodies are essentially excluded from the intracellular space, so the model deliberately carries no cellular sub-volume: the drug's accessible distribution space is vascular + interstitial only.
- **Units.** L. **RUN-VERIFIED:** Σ Vv = **2.0559 L**, Σ Vis = **8.5508 L**, Σ V = **62.43 L**.
- ⚠️ **`Vv` is not used by this module's own RHS.** Its only in-module reader is the `V_cen` convenience property (`:295`), which the RHS never calls (the RHS uses `V_PLASMA` instead, `:903`). `Vv` **is** live downstream — it is exported to the per-cell PK engine (`run_tce_pd_reval.py:51`). So: dead-here, live-downstream. Do not delete it.

---

### EQ-4 — Plasma concentration (`:903–904`)
```
V_pl  = V_PLASMA = 3.1 L
C_pl  = A_pl / V_pl
```
- **Biological meaning.** The measurable clinical quantity. The comment at lines 898–902 documents a real modelling decision worth preserving: the organ vascular spaces hold drug *at plasma concentration and are part of the same circulating plasma*, so they are **not** added as a separate central volume. Doing so had put V_c at 5.16 L and depressed C₀ by ~40% below the clinical value. The 3.1 L physiological plasma volume is what 2-compartment clinical fits actually recover (~3.0–3.5 L).
- **Mechanistic rationale / alternative rejected.** The rejected alternative (V_c = V_PLASMA + ΣVv) double-counts: the vascular sub-volumes are *inside* the plasma pool conceptually, not appended to it. This is correct and is one of the better-reasoned choices in the layer.
- **Units.** A_pl [amount]; V_pl [L]; C_pl [amount/L].

---

### EQ-5 — Interstitial concentration (`:905`)
```
Cis_i = A_is_i / Vis_i
```
- **Biological meaning.** Converts the extravasated amount in organ *i* into the local free-drug concentration that every downstream mechanism reads (target binding, TMDD, trimer formation, liver tox).
- **Units.** amount/L.

---

### EQ-6 — Two-pore convective extravasation, plasma → interstitium (`:908`)
```
J_extrav,i = k_dist · L_i · (1 − σ_V,i) · C_pl
```
- **Biological meaning.** Antibody leaves the blood by **convection with the filtered fluid**, not by diffusion. The fluid flux out of the capillary is the lymph flow `L_i`; the fraction of *antibody* that rides along with it is `(1 − σ_V,i)`, where σ_V is the **vascular reflection coefficient** — the fraction of protein "reflected" (held back) by the endothelial pore. σ_V → 1 is a tight endothelium (brain, 0.99: almost nothing crosses); σ_V → 0.75 is a leaky one (tumour: the EPR-like leak that lets an antibody reach a solid tumour at all).
- **Mechanistic rationale.** This is the *entire* reason a whole-body antibody model needs organ-specific σ_V: the biodistribution differences between brain, muscle and tumour are almost entirely a reflection-coefficient story, not a blood-flow story. The rejected alternative — permeability-surface-area diffusion (`PS·(C_pl − Cis)`) — is the small-molecule law and is wrong for a protein whose transport is convection-dominated.
- **Units.** `k_dist` dimensionless; `L_i` L/day; `σ_V` dimensionless; `C_pl` amount/L → `J_extrav` amount/day.
- **Note the structure:** `k_dist` multiplies `L`, always. It never appears alone. This is the root of the identifiability problem in §5.3.

---

### EQ-7 — Lymphatic return, interstitium → lymph pool (`:909`)
```
J_return,i = k_dist · L_i · (1 − σ_L) · Cis_i
```
- **Biological meaning.** Interstitial fluid drains into the initial lymphatics and is carried back to the blood. `σ_L` is the **lymphatic reflection coefficient** — how much antibody the lymphatic wall holds back on the way out. Lymphatic capillaries are far more permeable than blood capillaries, so physiologically σ_L should be *small* (protein returns freely).
- **Mechanistic rationale.** Together EQ-6/EQ-7 are the "2-pore" convective couple: in at (1−σ_V), out at (1−σ_L). Their **ratio sets the steady-state interstitial:plasma partitioning** (EQ-14), and their absolute scale sets the distribution *rate*. That clean separation of "extent" from "rate" is the design intent of having both σ_L and k_dist.
- ⚠️ **σ_L is a single scalar shared by all 15 organs** (it is a `PBPK` dataclass field, not a per-organ column), whereas σ_V is per-organ. Shah-Betts likewise assigns one σ_L to all tissues, so the *structure* matches — but the *value* does not (§3.3, §5.1).
- **Units.** as EQ-6.

---

### EQ-8 — Total plasma-side clearance rate constant (`:930`)
```
k_cat = CLup · (1 − fFcRn_eff) + k_renal(mw)
```
- **Biological meaning.** Two mechanistically distinct, **format-dependent** elimination routes, summed:
  1. **Endosomal catabolism of the non-salvaged fraction.** Vascular endothelium pinocytoses plasma at rate `CLup`. Of what is taken up, a fraction `fFcRn_eff` is bound by FcRn in the acidified endosome and **recycled back out intact**; the remainder `(1 − fFcRn_eff)` is trafficked to the lysosome and destroyed. This is *the* reason an IgG lives for weeks.
  2. **Size-gated renal filtration + proteolysis** (EQ-10).
- **Mechanistic rationale / what was rejected.** The platform (Shah-Betts) writes FcRn as an explicit endosomal compartment with on/off rates and a recycling flux. This code replaces that with a **quasi-steady-state recycled fraction**: if endosomal binding equilibrates fast relative to pharmacology, the net effect of the whole endosome is simply to scale down the catabolic rate by `(1 − fFcRn)`. That is a legitimate QSS reduction and it buys a large numerical-stiffness win. **What it costs:** FcRn saturation is now impossible to represent (a high-dose IgG cannot out-compete endogenous IgG for FcRn in this model), and `Kdeg` becomes meaningless (§3.3 — it is dead).
- **Units.** /day. **RUN-VERIFIED (default IgG):** `CLup·(1−fFcRn)` = **0.03503/day**, `k_renal` = **0.004545/day**, `k_cat` = **0.039575/day**.
- 🚩 **STALE COMMENT.** Line 188 asserts *"CLup*(1-fFcRn) = 0.0385/day"*. The live value is **0.03503/day**. The comment is arithmetically inconsistent with the live `CLup = 0.3503` and `fFcRn = 0.90` (0.3503 × 0.10 = 0.03503). It appears to be a leftover from an earlier `CLup ≈ 0.385`. Harmless to the model (the comment is not executed) but it means the comment block cannot be trusted as documentation — which is precisely why this doc re-derives everything.

---

### EQ-9 — Effective FcRn salvage (the Fc switch) (`:212–214`)
```
fFcRn_eff = fFcRn        if has_fc
          = 0            otherwise
```
- **Biological meaning.** FcRn salvage requires an intact Fc. A tandem-scFv BiTE has none, so **every** pinocytosed molecule is catabolised. This one boolean is half of the explanation for why T-cell-engager plasma half-lives span two orders of magnitude across formats.
- **Mechanistic rationale.** Encoding the format as a *structural switch* rather than as a free per-molecule clearance knob is the right call: it forces the half-life difference to *emerge* from the molecule's architecture instead of being fitted per molecule. This is the strongest "emergence" claim T1 can legitimately make (§4).
- **Units.** dimensionless.

---

### EQ-10 — Size-gated renal clearance (glomerular sieving Hill) (`:216–218`)
```
k_renal(mw) = k_renal_max / ( 1 + (mw / mw50_renal)^hill_renal )
```
- **Biological meaning.** The glomerulus sieves by size with a sharp cutoff near serum albumin (~66–69 kDa). Below it, a protein is filtered into urine and proteolysed in the tubule; above it, essentially not at all. A ~54 kDa BiTE is *below* the cutoff and is renally cleared; a 147 kDa IgG is far above it.
- **Mechanistic rationale.** This is the **second half** of the no-Fc story, and it matters: losing FcRn alone takes an IgG's half-life only from ~weeks to ~days. Getting all the way down to blinatumomab's ~2 h *requires* renal filtration of the small scaffold. Writing it as a Hill in MW (rather than a per-molecule clearance) again keeps it structural.
- **Units.** `k_renal_max` /day; `mw`, `mw50_renal` kDa; `hill_renal` dimensionless.
- **RUN-VERIFIED:** IgG (146.9 kDa) → **0.004545/day**; BiTE (54 kDa) → **8.0097/day** — a 1762× swing across the cutoff. The Hill exponent 10 makes the switch appropriately sharp.
- 🚩 **The comment's "~0 for a 150 kDa IgG" is an overstatement.** 0.004545/day is not zero: it is **11.5%** of the IgG's total `k_cat` (0.004545/0.039575), and it shortens the IgG's elimination half-life from ln2/0.03503 = 19.8 d to ln2/0.039575 = 17.5 d. Small, but not the no-op the comment implies, and it is silently absorbed into the `CLup` calibration.

---

### EQ-11 — Plasma pool ODE (`:931–932`)
```
infusion = inf_rate   if t ≤ inf_duration   else 0                    (:931)
dA_pl/dt = infusion − k_cat·A_pl − Σ_i J_extrav,i + k_lymph_return·A_ly
```
- **Biological meaning.** Mass balance on the central pool: drug in (infusion; a bolus enters instead as the initial condition `A_pl(0)=dose`, `:1228`), drug destroyed (catabolism + renal), drug convected out to all 15 interstitia, drug returned from the lymph pool.
- **Mechanistic rationale.** Elimination is **first-order on the plasma amount** — i.e. non-saturable. All *nonlinearity* in the PK (TMDD) lives on the interstitial side (EQ-13). That is a deliberate and defensible split: catabolic/renal clearance really is linear over the clinical range, while target-mediated clearance really is saturable.
- **Units.** amount/day.

---

### EQ-12 — Interstitial and lymph ODEs (`:933–934`)
```
dA_is,i/dt = J_extrav,i − J_return,i − tmdd_sink_i
dA_ly/dt   = Σ_i J_return,i − k_lymph_return · A_ly
```
- **Biological meaning.** Each organ's interstitium fills by convection, empties by lymphatic drainage, and loses drug to target internalisation. All 15 lymphatic drains collect into **one systemic lymph pool** (the thoracic duct, abstracted), which empties into plasma at `k_lymph_return`.
- **Mechanistic rationale.** The lymph pool is what makes the return **delayed** rather than instantaneous, and that delay is a real second distribution phase. Rejected alternative: draining each interstitium straight back to plasma (no pool) — that would remove a genuine kinetic phase and collapse the model to a 2-compartment system.
- ⚠️ **The lymph pool is tracked as an AMOUNT with no volume.** `V_LYMPH = 2.6` L is declared at `:363` and **never read anywhere** (verified: `grep V_LYMPH` → the definition line only). It is a dead parameter (§3.3).
- **Units.** amount/day; `k_lymph_return` /day.

---

### EQ-13 — Target-mediated (TMDD) sink — the PK↔binding coupling (`:916–922`)
```
R_CD3,i = Rcap_CD3 · tcell_i                                   (:916)
R_TAA,i = Rcap_TAA · taa_i          (taa = 1 in tumour, else 0) (:917)
(Cb_mono, Cb_bridge) = bivalent_binding(Cis, R_CD3, R_TAA, p)  (:918)
tmdd_sink_i = (kint_mono·Cb_mono,i + kint_cplx·Cb_bridge,i) · Vis_i   (:922)
```
- **Biological meaning.** Drug that is bound to a receptor gets **internalised with the receptor** and destroyed. Because the receptor pool is *finite*, this elimination route **saturates** — so fractional clearance *falls* as dose rises. That is target-mediated drug disposition, and here it is not imposed: it emerges from the binding solution.
- **Mechanistic rationale.** The sign is the whole point. A phenomenological nonlinear-clearance term could be *fitted* to give TMDD; deriving it from a saturable receptor-occupancy solution makes the nonlinearity a *consequence* of finite antigen rather than a fitted curve. This is the layer's best emergence claim after the Fc switch.
- **Scope note.** The binding solver (`bivalent_binding`, `:689–761`) and its parameters (KD_CD3, KD_TAA, Rhoden avidity geometry) are **owned by the binding subsystem doc, not T1.** T1 owns only the *mass-balance coupling* — that the sink is subtracted from `A_is` and scaled by `Vis_i`.
- ⚠️ **The receptor field is STATIC.** `R_CD3,i = Rcap_CD3 · tcell_i` is a **constant** built from the frozen density array — it is **not** connected to the dynamic CD8/CD4/Treg state variables the model integrates. So T-cell expansion (which the PD layer simulates, sometimes 5–10×) does **not** increase the CD3 TMDD sink, and tumour killing does not shrink the TAA sink. The PK is therefore **one-way decoupled** from the PD. See §4 and §5.5.
- **Units.** `tmdd_sink` amount/day; `kint_*` /day; `Cb_*` amount/L.

---

### EQ-14 — [DERIVED, this task] Steady-state interstitial:plasma partitioning
Not a line of code — it is the *closed-form consequence* of EQ-6 = EQ-7, and it is the relation the `sigL` calibration comment (`:170–175`) invokes. I derive it here because the comment's claim needed checking.

Setting `J_extrav,i = J_return,i` at distribution equilibrium:
```
k_dist·L_i·(1−σ_V,i)·C_pl = k_dist·L_i·(1−σ_L)·Cis_i
  ⇒   Cis_i / C_pl = (1 − σ_V,i) / (1 − σ_L)          [k_dist and L_i CANCEL]
```
and the whole-body apparent volume ratio (adding the lymph pool, whose SS content is `Σ J_return / k_lymph_return`):
```
Vss/Vc = 1 + [ Σ_i Vis_i·(1−σ_V,i)/(1−σ_L)  +  k_dist·Σ_i L_i·(1−σ_V,i)/k_lymph_return ] / V_pl
              \_______ interstitial term, k_dist-FREE _______/   \___ lymph term, k_dist-LINEAR ___/
```
- **RUN-VERIFIED at defaults:** interstitial equivalent volume = **3.2182 L**; lymph equivalent volume = **0.1616 L**; ⇒ **Vss/Vc = 1 + (3.2182 + 0.1616)/3.1 = 2.0903.**
- ✅ **The `sigL = 0.85` calibration target is met:** the code comment claims "Vss/Vc ≈ 2.1"; the exact analytic value is **2.09**. The fit is real and it works.
- 🚩 **But the comment's k_dist-invariance claim is only half true.** Line 174 states the SS ratio is *"independent of k_dist, so this knob sets DISTRIBUTION EXTENT (Vss) without touching the rate."* The **per-organ** ratio `(1−σ_V)/(1−σ_L)` is indeed k_dist-free ✅ — but the **whole-body Vss/Vc is not**, because the lymph-pool term scales linearly with k_dist. **MEASURED:** Vss/Vc = **2.056 / 2.090 / 2.195** at k_dist = 1 / 3 / 9. A ~7% drift across a 9× k_dist change — small, but the clean "extent knob vs rate knob" orthogonality the design leans on is **approximate, not exact.**
- **Per-organ equilibrium partitioning at the live σ_L = 0.85** (RUN-VERIFIED) — note the values >1:

| organ | σ_V | Cis/C_pl |
|---|---|---|
| brain | 0.99 | 0.067 |
| lung, heart, muscle, skin, adipose, bone | 0.95 | 0.333 |
| kidney, stomach, small_int, large_int, pancreas | 0.90 | 0.667 |
| liver | 0.85 | **1.000** |
| spleen | 0.80 | **1.333** |
| **tumor** | 0.75 | **1.667** |

  A tumour interstitial concentration **1.67× plasma at equilibrium** is not a physiological antibody biodistribution. This is a direct consequence of σ_L = 0.85 (§5.1).

---

### EQ-15 — [DERIVED, this task] Terminal half-life is an EIGENVALUE, not `ln2/k_cat`
The TMDD-free PK block is a linear 17-state system (plasma + 15 interstitia + lymph). Its half-life structure is therefore the eigen-spectrum of that matrix, **not** `ln2/k_cat` — a distinction the code comments repeatedly blur.

**RUN-VERIFIED (exact eigen-decomposition, cross-checked against the LSODA solution):**

| | IgG (146.9 kDa, has_fc) | BiTE (54 kDa, no Fc) |
|---|---|---|
| `k_cat` | 0.03958 /day | 8.35996 /day |
| `ln2/k_cat` (elimination rate constant) | 17.51 d | **1.99 h** |
| **True terminal t½ (slowest eigenvalue)** | **38.11 d** | 5.13 d |
| AUC fraction in that terminal phase | 95.8% | **0.34%** |
| **Time for plasma to fall to 50% of C₀** | 3.69 d | **1.78 h** |

*(Both columns are the TMDD-FREE linear block. With TMDD on, the same two runs give BiTE 1.78 h — unchanged, because TMDD is an interstitial sink and does not touch the plasma ODE directly — and IgG 2.64 d at dose = 1, where TMDD is far from saturated. RE-MEASURED 2026-07-13: an earlier revision of this doc reported the BiTE value as 1.80 h; the reproduced value is **1.7830 h** (linear block) / **1.7829 h** (full model).)*

Reading this table correctly matters, and it *reverses* a naive verdict:
- **BiTE — the calibration is SOUND.** The `k_renal_max = 8.70` fit targets `ln2/k_cat = 1.99 h`, and the model's **clinically observable** half-life (dominant phase; time to 50% of C₀) is **1.78 h** against the BLINCYTO label's **2.11 h (SD 1.42)** — ✅ label value independently confirmed in the 2026-07-13 pass; the code comment only says "~2.1 h", so the 2.11/1.42 precision is a doc-side literature lookup, not a code value. ✅ Well within the clinical SD. The 5.13-day "terminal" eigenvalue is a **numerical artifact carrying 0.34% of AUC** — a low-amplitude redistribution tail from interstitial drug trickling back through the lymph pool, since a no-Fc construct has no interstitial elimination of its own. It is invisible in practice but it is there.
- **IgG — the calibration claim does NOT hold as stated.** See EQ-16 and §5.4.

---

### EQ-16 — [DERIVED, this task] The IgG terminal half-life is DOSE-DEPENDENT
Because TMDD (EQ-13) is saturable, the IgG's apparent terminal half-life is **not a fixed model property** — it slides with dose.

**MEASURED** (backbone arm, log-linear fit over days 60–120):

| dose (abstract units) | apparent terminal t½ | tumour:plasma AUC₀₋₅₆ |
|---|---|---|
| 0.1 | 9.11 d | 0.65 |
| 1.0 | 9.11 d | 0.73 |
| 10 | 9.12 d | 1.18 |
| 100 | 9.42 d | 1.58 |
| 1 000 | 19.27 d | 1.66 |
| 10 000 | 35.24 d | 1.66 |
| 100 000 (TMDD-saturated limit) | 37.81 d | 1.66 |

- TMDD-dominated (low dose) → **9.1 d**. TMDD-saturated (high dose) → **38.1 d**, converging exactly on the linear eigenvalue of EQ-15. ✅ (The two independent methods agree, which validates both.)
- 🚩 **The `CLup` comment (`:184–187`) claims the backbone terminal t½ "matches the mosunetuzumab clinical anchor 16.1 d exactly (FDA LUNSUMIO label)."** The **16.1 d anchor itself is REAL** — ✅ independently verified against the FDA LUNSUMIO label (terminal t½ 16.1 d, population-PK estimate at steady state). But the model **only produces 16.1 d at one particular dose**: solving for it gives **dose ≈ 695** abstract units. The scored dose grid is `logspace(-0.5, 3.0)` = 0.316 → 1000 (`:1346`), so 695 sits *inside* the grid — but at the grid's low end the model's IgG half-life is 9.1 d, not 16.1 d. **"Matches exactly" is true only at an unstated dose.** It is not a property of the calibration. See §5.4.
- 🚩 **The same comment's "tumor/plasma AUC ratio 0.234" could NOT be reproduced at any dose** (measured range 0.65 → 1.66, RE-CONFIRMED 2026-07-13). Marking **[UNREPRODUCED]**. *Inference (flagged as such, not fact):* with the *platform's* σ_L = 0.2 the model's equilibrium tumour ratio would be (1−0.75)/(1−0.2) = **0.3125** — far closer to the stale comment than the live 1.667, which is consistent with the 0.234 figure predating the σ_L = 0.85 re-fit. I cannot confirm that from the code and do not assert it. *(A previous revision also asserted 0.234 was "suspiciously close to the tumour antibody-biodistribution-coefficient (~0.24) from the Shah/Betts ABC work". That numeric attribution was **not verified against the ABC paper** and has been removed rather than propagated.)*

---

## 3. PARAMETERS OWNED

### 3.1 Frozen organ physiology — `_PBPK_TISSUES` (`:82–99`)
15 organs × 6 columns. **Values below are read directly from the code table; `Q`, `L`, `Vv`, `Vis` are the RUN-VERIFIED derived values printed from the live module.**

| organ | V (L) | Qfrac | **Q (L/d)** | **L (L/d)** | σ_V | fV | **Vv (L)** | fIS | **Vis (L)** | portal |
|---|---|---|---|---|---|---|---|---|---|---|
| lung | 0.50 | 1.000 | 5000.00 | 10.0000 | 0.95 | 0.105 | 0.0525 | 0.188 | 0.0940 | F |
| heart | 0.33 | 0.040 | 211.64 | 0.4233 | 0.95 | 0.157 | 0.0518 | 0.320 | 0.1056 | F |
| kidney | 0.31 | 0.190 | 1005.29 | 2.0106 | 0.90 | 0.105 | 0.0325 | 0.200 | 0.0620 | F |
| brain | 1.45 | 0.120 | 634.92 | 1.2698 | 0.99 | 0.037 | 0.0536 | 0.150 | 0.2175 | F |
| muscle | 30.0 | 0.170 | 899.47 | 1.7989 | 0.95 | 0.026 | 0.7800 | 0.120 | 3.6000 | F |
| skin | 3.40 | 0.050 | 264.55 | 0.5291 | 0.95 | 0.038 | 0.1292 | 0.302 | 1.0268 | F |
| adipose | 13.0 | 0.050 | 264.55 | 0.5291 | 0.95 | 0.020 | 0.2600 | 0.135 | 1.7550 | F |
| bone | 10.0 | 0.050 | 264.55 | 0.5291 | 0.95 | 0.041 | 0.4100 | 0.100 | 1.0000 | F |
| stomach | 0.15 | 0.010 | 52.91 | 0.1058 | 0.90 | 0.038 | 0.0057 | 0.170 | 0.0255 | **T** |
| small_int | 0.65 | 0.100 | 529.10 | 1.0582 | 0.90 | 0.038 | 0.0247 | 0.200 | 0.1300 | **T** |
| large_int | 0.37 | 0.040 | 211.64 | 0.4233 | 0.90 | 0.038 | 0.0141 | 0.200 | 0.0740 | **T** |
| pancreas | 0.18 | 0.010 | 52.91 | 0.1058 | 0.90 | 0.038 | 0.0068 | 0.180 | 0.0324 | **T** |
| spleen | 0.19 | 0.030 | 158.73 | 0.3175 | 0.80 | 0.110 | 0.0209 | 0.200 | 0.0380 | **T** |
| liver | 1.80 | 0.065 | 343.92 | 0.6878 | 0.85 | 0.115 | 0.2070 | 0.200 | 0.3600 | F |
| **tumor** | 0.10 | 0.020 | 105.82 | 0.2116 | **0.75** | 0.070 | 0.0070 | 0.300 | 0.0300 | F |

**Provenance for the whole table: `[ASSUMED: "reference-human (ICRP/Brown 71 kg)" per the header comment at line 75]` — with the following honest caveats:**
- The code names its source only as *"reference-human (ICRP/Brown 71 kg)"* (`:75`). There is **no PMID, no DOI, and no per-value citation in the code.** "Brown" is presumably the standard physiological-parameter compilation used across the PBPK field, but **I will not invent the citation.** → **[UNVERIFIED CITATION — the source is named but not cited; individual values are not traceable from the code.]**
- **σ_V (all 15 values): [ASSUMED — UNSOURCED IN CODE].** No citation accompanies them. They are *ordered* physiologically sensibly (brain 0.99 tightest → tumour 0.75 leakiest), which is the right qualitative structure and matches Shah-Betts' "adjusted according to vascular pore size." But the specific numbers have no in-code source. For comparison, the platform's two-pore coefficients are **σ₁ ≈ 0.883–0.987 (tight pore)** and **σ₂ ≈ 0.311–0.837 (leaky pore)** *(literature values, re-checked 2026-07-13; NOT in the code)*: the code's single 0.75–0.99 band straddles the two and does **not** reproduce the platform's small+large pore split. **[UNSOURCED — TBD]**
  > 🚩 **RETRACTED NUMBERS (adversarial pass 2026-07-13).** An earlier revision of this doc stated these platform coefficients as "σ₁ ≈ 0.69–0.999 (mean 0.908) and σ₂ ≈ 0.258–0.841 (mean 0.579)". **Those four numbers appear nowhere in the code, nowhere in the source, and were fabricated.** They are replaced above with the published ranges. No mean is quoted, because no mean was verified.
- **`portal` column: DEAD.** Every value is read into `PB.portal` / `PB.i_portal` (`:136`, `:147`) and then **never consumed anywhere** (verified: zero hits). It encodes the splanchnic→liver drainage topology that the header describes but the RHS does not implement (§5.2).

### 3.2 Global scalars

| symbol | value | units | tag | source | rationale |
|---|---|---|---|---|---|
| `_PLASMA_CO` (`:100`) | 5000.0 | L/day | **[ASSUMED — UNSOURCED IN CODE]** | comment: *"cardiac output × plasma fraction"*; no citation | A 71-kg reference human has CO ≈ 5 L/min ≈ 7200 L/day; × plasma fraction (1−Hct ≈ 0.55) ⇒ ≈ **3960 L/day**. The code's **5000** is ~26% higher and is a round number. **It does not matter** — see §5.3: it is exactly degenerate with `k_dist`. |
| `_LYMPH_RATIO` (`:101`) | 1/500 | — | ✅ **[MEASURED: Shah & Betts 2012, PMID 22143261]** | verified verbatim in the source | The signature platform relation. Faithfully reproduced. |
| `V_PLASMA` (`:362`) | 3.1 | L | **[MEASURED/standard physiology]** — reference-human plasma volume | no in-code citation, but 3.1 L is the standard value and the comment's defence (2-cpt clinical fits recover 3.0–3.5 L) is sound | Chosen deliberately over `V_PLASMA + ΣVv = 5.16 L`, which depressed C₀ ~40% (EQ-4). Good call. |
| `V_LYMPH` (`:363`) | 2.6 | L | **💀 DEAD — never read** | — | Declared and never referenced. The lymph pool is amount-only. |

### 3.3 `PBPK` dataclass — the tunable system parameters (`:166–295`)

| symbol | value | units | tag | source / rationale |
|---|---|---|---|---|
| `sigL` (`:170`) | **0.85** | — | 🔴 **[FITTED: to Vss/Vc ≈ 2.1]** | The code is **honest** that this is calibrated ("CALIBRATED so the steady-state interstitial:plasma amount ratio gives Vss/Vc ≈ 2.1"). ✅ The fit works (analytic Vss/Vc = **2.0903**, EQ-14). ❌ **But Shah & Betts assign σ_L = 0.2** (VERIFIED against the source). The code's value is **4.25× the platform's** and makes the leaky-organ interstitium *concentrate* drug above plasma (tumour 1.67× — EQ-14). This is the single largest deviation in T1. See §5.1. Named comparators in the comment ("pembrolizumab 2.17, trastuzumab 2.7, mosunetuzumab 2.1") are **[UNVERIFIED CITATION]** — I did not confirm them. |
| `k_dist` (`:176`) | **3.0** | — | 🔴 **[FITTED: to the α-phase, "calibrated vs pembrolizumab day-1/3/7 fall"]** | A pure rate multiplier on BOTH EQ-6 and EQ-7. Intent: set distribution *rate* without touching *extent*. **That orthogonality is approximate, not exact** (EQ-14: Vss/Vc drifts 2.056→2.195 over k_dist 1→9). The pembrolizumab day-1/3/7 target is **[UNVERIFIED CITATION]**. |
| `fFcRn` (`:181`) | 0.90 | — | **[ASSUMED — UNSOURCED IN CODE]** | "FcRn salvage fraction of pinocytosed mAb". No citation. Physiologically plausible (FcRn salvage is highly efficient) but **not sourced**, and it is fully confounded with `CLup` — only the *product* `CLup·(1−fFcRn)` enters EQ-8, so `fFcRn` and `CLup` are **not separately identifiable** from plasma PK either. |
| `CLup` (`:182`) | 0.3503 | /day | 🔴 **[FITTED: to mosunetuzumab terminal t½ = 16.1 d]** | Anchor ✅ **VERIFIED REAL** (FDA LUNSUMIO label, terminal t½ 16.1 d, pop-PK at steady state). 🚩 **But the claim "matches 16.1 d exactly" is DOSE-CONDITIONAL** — the model gives 9.1 d at low dose, 38.1 d TMDD-saturated, and hits 16.1 d only at dose ≈ 695 (EQ-16). The suspiciously precise 4-significant-figure value (0.3503) is a fitted artefact, not a measurement. |
| `Kdeg` (`:189`) | 26.0 | /day | **💀 DEAD — never read** (grep: hits only at `:76` comment and `:189` definition) | Comment calls it "catabolism of UNsalvaged endosomal mAb (fast)" and the header (`:76`) lists it as one of "the 4 platform system parameters." **It is not in any equation.** It is a vestige of the explicit-endosome formulation that EQ-8's QSS reduction replaced — once catabolism is written as `CLup·(1−fFcRn)`, the endosomal degradation rate is assumed infinitely fast and drops out. Honest statement: **the header's "4 platform parameters" is really 3.** |
| `k_lymph_return` (`:190`) | 24.0 | /day | **[ASSUMED — UNSOURCED IN CODE]** | "fast turnover" — i.e. the lymph pool empties with a ~1 h time constant. No citation. Sets the depth of the second distribution phase. Also appears in the Vss/Vc lymph term (EQ-14). |
| `mw_kda` (`:198`) | 146.9 | kDa | **[MEASURED — standard IgG mass]** | Overridable per construct format (`apply_construct_format`, `:1586`). Drives EQ-10 and the µg/mL↔nM readout conversion (`:1279` — a readout, not one of the equations above; an earlier revision of this doc cross-referenced a nonexistent "EQ-17"). |
| `has_fc` (`:199`) | True | bool | **[STRUCTURAL — construct property]** | The Fc switch (EQ-9). Not a fitted knob — an architectural fact about the molecule. |
| `k_renal_max` (`:205`) | 8.70 | /day | 🟡 **[FITTED: to blinatumomab t½ ≈ 2.1 h]** | Anchor ✅ **VERIFIED REAL** (FDA BLINCYTO label: mean t½ **2.11 h, SD 1.42**; MW ~54 kDa — label value confirmed 2026-07-13; the *code* comment says only "~2.1 h"). ✅ **And this fit actually reproduces:** the model's BiTE falls to 50% of C₀ in **1.78 h** (EQ-15). This is the best-validated fitted parameter in T1. |
| `mw50_renal` (`:209`) | 69.0 | kDa | **[MEASURED — albumin ~66–69 kDa glomerular sieving midpoint]** | Physiologically the right anchor; no in-code citation but the value is textbook. |
| `hill_renal` (`:210`) | 10.0 | — | **[ASSUMED: "sharp glomerular threshold"]** | A steepness knob. n=10 is very sharp; no source. It is not sensitive at the two formats actually used (both are far from the midpoint), so it is effectively a shape assumption with little leverage. |

**TMDD-coupling parameters read by the PK mass balance** (owned by the binding subsystem doc; listed for completeness because EQ-13 is a PK sink):

| symbol | value | units | tag | note |
|---|---|---|---|---|
| `Rcap_CD3` (`:226`) | 2.0 | conc | **[ASSUMED — UNSOURCED]** | "CD3 receptor capacity scale (× tissue T-cell density)". Abstract units. |
| `Rcap_TAA` (`:227`) | 6.0 | conc | **[ASSUMED — UNSOURCED]** | tumour-only TAA capacity. |
| `kint_cplx` (`:287`) | 0.90 | /day | **[ASSUMED — UNSOURCED]** | internalisation of the bridged complex. |
| `kint_mono` (`:288`) | 0.25 | /day | **[ASSUMED — UNSOURCED]** | internalisation of singly-bound drug. Bridged > mono is mechanistically sensible (crosslinking accelerates internalisation) but the 3.6× ratio has no source. |

### 3.4 Immune-density arrays — `_TCELL_DENSITY` (`:111–115`), `_MYELOID_DENSITY` (`:119–123`)
Two 15-element arrays, spleen-normalised to 1.0. They set (a) where the CD3-arm TMDD sink is strongest (EQ-13) and (b) the liver-tox substrate (Kupffer density).

Values (read from code): T-cell — lung 0.309, heart 0.117, kidney 0.292, brain 0.020, muscle 0.183, skin 1.510, adipose 0.767, bone 0.860, stomach 1.028, small_int 1.240, large_int 2.359, pancreas 0.182, spleen 1.000, liver 0.429, tumor 0.400. Myeloid — lung 0.924, heart 0.306, kidney 0.106, brain 0.050, muscle 0.178, skin 0.572, adipose 0.564, bone 1.287, stomach 0.305, small_int 0.035, large_int 0.017, pancreas 0.394, spleen 1.000, liver **0.839** (Kupffer), tumor 0.300.

- **Tag: ✅ [MEASURED: Tabula Sapiens Consortium, Science 2022, DOI 10.1126/science.abl4896 — DOI present in the code at `:110`; provenance CSV present and REPRODUCED].** The code states these were computed from the Tabula Sapiens single-cell atlas (~1.14M cells, CZ CELLxGENE Census stable 2025-11-08) and normalised to spleen, and points to a provenance file: *"See qsp_tabula_sapiens_densities.csv for provenance"* (`:109`).
- 🚩 **RETRACTION (adversarial pass 2026-07-13) — the previous revision of this section was WRONG.** It asserted, in bold, that *"THAT FILE DOES NOT EXIST"*, that a repo-wide `find` returned nothing, and that a grep for the literal values (`2.359`, `0.839`) found "only coincidental hits". **All three of those statements are false.** The file exists:
  ```
  deliverables/06_qsp_science-qsp/v2_enhanced/qsp_tabula_sapiens_densities.csv
  ```
  Its columns are `organ, ts_tissue, f_Tcell, f_CD8, f_CD4, f_Treg, f_Myeloid, n, tcell_density, myeloid_density` — i.e. exactly the per-tissue immune-cell-**fraction** table the code claims, with the Tabula Sapiens tissue name it was mapped from and the cell count `n` behind each fraction (e.g. spleen n=70,448; lung n=65,847; large_int n=10,387).
- ✅ **ALL 30 DENSITIES REPRODUCE FROM IT, EXACTLY** (re-derived this pass): `density_organ = f_organ / f_spleen`, per array. Spot checks — large_int T-cell 0.6058/0.2568 = **2.359** ✅; liver myeloid 0.3076/0.3666 = **0.839** ✅; lung T-cell 0.0794/0.2568 = **0.309** ✅; large_int myeloid 0.0061/0.3666 = **0.017** ✅. 15/15 T-cell and 15/15 myeloid values match the code arrays to the printed precision. **The immune-density provenance is SOUND and needs no regeneration.**
- **brain (0.020 / 0.050) and tumor (0.400 / 0.300) are explicitly NOT from Tabula Sapiens** — and the CSV agrees: both rows carry `ts_tissue = (none)`, `n = 0`, and blank fractions, with the density column carrying the prior directly. The code says brain is "immune-privileged; absent in TS" and tumor uses a "TIL prior from project brief." → **[ASSUMED: literature/project prior]** for those 4 values. This is the one honest gap here, and the CSV itself marks it.
- `taa` (`:154`): 1.0 in tumour, 0 elsewhere. **[STRUCTURAL — by construction]**, the definition of a tumour-associated antigen in this abstraction.

### 3.5 Spatial hook
`spatial_exposure` (`:161`): 15×1.0 by default. **[STRUCTURAL — identity by default]**, RUN-VERIFIED all-ones. Loaded from `spatial_exposure_perorgan.csv` by `apply_spatial_exposure` (`:1542`), so it can only *reduce* exposure. With the default all-ones the model is exactly well-mixed. Clean interface; nothing to audit until the CSV lands.
- ⚠️ **Interval nit (code vs comment).** The code comment and the loader docstring both say the value is clipped to **(0,1]**, but the implementation is `min(max(v, 0.0), 1.0)` (`:1554`) → the admissible range is **[0,1]**: a CSV value of `0` is accepted and would zero an organ's engagement-driving concentration entirely. Not a bug today (no CSV), but the contract as written is not the contract as coded.

### 3.6 💀 Confirmed DEAD parameters — 6 entries / 9 symbols (defined, never read — verified by exhaustive grep, re-run 2026-07-13)

| symbol | line | why it's dead |
|---|---|---|
| `Kdeg` | 189 | superseded by the QSS FcRn reduction in EQ-8 |
| `V_LYMPH` | 363 | lymph pool is amount-only; no volume needed |
| `avidity` | 228 | explicitly marked "DEPRECATED heuristic cap… retained for back-compat"; replaced by the geometric Rhoden avidity |
| `hill_CD3` / `hill_TAA` | 255–256 | "avidity Hill exponent per extra arm (set at build)" — never read |
| `n_costim` | 250 | settable via `apply_construct_format` (`:1592`) but has **no consumer**; reserved for a future format library |
| `portal` / `i_portal` / `i_parallel` | 136, 147–148 | the blood-flow topology described in the header but not implemented |

None of these affect results. They are documented because a reader tracing the header's "4 platform parameters" or its portal-drainage description **will otherwise conclude the model does something it does not do.**

---

## 4. WHAT IS EMERGENT vs IMPOSED

The model's selling point is emergence. Here is precisely where T1 earns that and precisely where it stops.

### ✅ Genuinely EMERGENT in T1
1. **The two-order-of-magnitude half-life span across formats.** This is T1's strongest and most defensible claim. An IgG's ~weeks and a BiTE's ~2 h come out of **two structural facts** — `has_fc` (EQ-9) and `mw_kda` (EQ-10) — not from a per-molecule clearance knob. Change the molecule's architecture and the PK follows. **RUN-VERIFIED (TMDD-free linear block, so the contrast is pure format):** flipping only `has_fc=False` and `mw_kda=54` moves the time-to-50%-of-C₀ from **3.69 d → 1.78 h** with no other change. That is real mechanistic emergence.
2. **TMDD nonlinearity, with the correct sign.** Fractional clearance *falls* as dose rises because the receptor pool saturates (EQ-13). Not imposed as a nonlinear clearance curve — it is a consequence of finite antigen. **RUN-VERIFIED:** apparent terminal t½ 9.1 d → 38.1 d across a 10⁶ dose range.
3. **Organ-differential biodistribution.** The 6–25× spread in interstitial exposure (brain 0.067 → tumour 1.667, EQ-14) emerges from per-organ σ_V and volumes; no organ has a fitted "penetration factor."
4. **The multi-phasic plasma curve.** α/β phases are eigenvalues of the transport topology (EQ-15), not fitted exponentials.

### ⛔ IMPOSED (handed to T1 as constants)
1. **All physiology** — V, Qfrac, σ_V, fV, fIS. Frozen at import, never adapted.
2. **σ_L = 0.85 and k_dist = 3.0** — **both FITTED**, to Vss and to the α-phase respectively. These are the two knobs that set the entire distribution behaviour, and neither is a measurement. Everything "emergent" about biodistribution is emergent *conditional on these two fitted numbers*.
3. **CLup = 0.3503** — FITTED to a clinical half-life. So the IgG half-life is **not** emergent; it is fitted. Only the *format-dependence* of half-life (the IgG↔BiTE contrast) is emergent, because the BiTE inherits `CLup` unchanged and gets its short life from structure.
4. **The receptor field is static.** `R_CD3 = Rcap_CD3 · tcell` is a **frozen array** (EQ-13). It does not respond to the T-cell expansion the PD layer simulates, nor does the TAA sink shrink as the tumour is killed. **The PK is one-way decoupled from the PD.** So "TMDD emerges" is true *within a fixed receptor field*; the receptor field itself is imposed and static. This is the sharpest limit on the emergence claim in this layer.
5. **The immune densities** (Tabula Sapiens) are data-grounded *inputs*, not emergent. Their provenance is intact and reproducible (§3.4) — but they are still an imposed, static field, and the 4 brain/tumor values are priors, not atlas data.

### Honest one-line summary
> T1 **imposes** the physiology and **fits** three transport/clearance constants (σ_L, k_dist, CLup); from those it **emerges** the format-dependence of PK, saturable TMDD, organ-differential exposure, and the multi-exponential curve shape. The absolute IgG half-life is fitted; the *contrast* between formats is earned.

---

## 5. KNOWN LIMITATIONS / OPEN QUESTIONS — what a reviewer will attack

### 5.1 🔴 σ_L = 0.85 is 4.25× the platform value, and it makes leaky organs concentrate antibody
**This is the most serious scientific issue in T1 and the first thing a PK reviewer will find.**

Shah & Betts assign a **lymphatic reflection coefficient σ_L = 0.2** (VERIFIED against the source). The code uses **0.85**, and is candid that it is *calibrated to make Vss/Vc ≈ 2.1*. The fit succeeds numerically (EQ-14: 2.0903 ✅) — but consider what it does mechanistically. The equilibrium partitioning is `(1−σ_V)/(1−σ_L)`:

| | with platform σ_L = 0.2 | with code's σ_L = 0.85 |
|---|---|---|
| tumour (σ_V=0.75) | 0.25/0.80 = **0.3125** | 0.25/0.15 = **1.667** |
| spleen (σ_V=0.80) | **0.250** | **1.333** |
| liver (σ_V=0.85) | **0.1875** | **1.000** |
| muscle (σ_V=0.95) | **0.0625** | **0.333** |

Setting σ_L ≈ σ_V (0.85 vs 0.85 in liver!) means the lymphatic wall reflects antibody **as strongly as the blood capillary does** — which is not what a lymphatic capillary is. The consequence is that drug that convects into a leaky interstitium **cannot easily get back out**, so it accumulates to *above-plasma* concentrations. A tumour interstitium at **1.67× plasma at equilibrium** is not an antibody biodistribution any experimentalist would recognise: measured tumour antibody-biodistribution coefficients are well below 1 (order 0.1–0.3) — *this range is a literature recollection, **NOT** verified against a source in this pass and **not** present in the code; treat it as a direction-of-error statement, not a citable number.* The direction of the error is what matters and it is unambiguous: σ_L = 0.85 makes leaky organs **concentrate** antibody above plasma, which is the wrong sign of physiology.

**Why it was done is obvious and defensible-in-isolation:** with σ_L = 0.2, the interstitial equivalent volume would be `Σ Vis_i(1−σ_V,i)/0.8` = 0.4827/0.8 = 0.60 L, giving Vss/Vc ≈ 1.2 — far below the class-typical ~2.1. So the model **cannot simultaneously reproduce a physiological σ_L and a class-typical Vss** with only these two knobs. σ_L was sacrificed to save Vss.

**What that tells us:** the real reason Vss/Vc ≈ 2.1 for a mAb is that the *available* interstitial volume for an antibody is large (the ~12–15 L of accessible extracellular fluid), whereas this model's Σ Vis is only **8.55 L** and its (1−σ_V)-weighted version is only **0.48 L**. The structure is under-supplied with distribution volume, and σ_L was pushed to 0.85 to manufacture the missing Vss. **The honest fix is structural (more accessible interstitial volume / a proper 2-pore small+large pore split), not a further σ_L tweak.**

**Impact on conclusions:** the model **over-predicts interstitial (and therefore tumour and liver) exposure relative to plasma by roughly 5×** (denominator 0.15 vs 0.80). Because the costim window score is a *comparative* ranking across arms that all share this physiology, the **nomination is likely robust** to it — but any **absolute** tumour- or liver-exposure statement from this model should be treated as unreliable.

### 5.2 🔴 The advertised blood-flow topology is not implemented
The header (lines 72–73) states: *"Portal organs drain to liver; liver + all organs drain to a venous pool → lung (series) → arterial pool."* **None of this is in the live RHS.** `PB.Q`, `PB.portal`, `PB.i_portal` and `PB.i_parallel` are computed and **never read** (verified). The live model is a **star topology**: one plasma pool feeding 15 independent interstitia in parallel, all draining to one lymph pool.

Consequences: (a) there is **no hepatic first-pass** and no splanchnic→liver sequence, so the liver sees plasma drug directly, not gut-effluent drug — relevant because the liver compartment is exactly where the 4-1BB hepatotox readout lives; (b) the lung's series position (which for an IV bolus means the whole dose transits the lung first) is absent. For an IgG distributing over days these are probably second-order — but the docstring **claims a structure the code does not have**, and that is a documentation defect that must be fixed before this is shown to a committee.

### 5.3 🔴 `_PLASMA_CO`, `_LYMPH_RATIO` and `k_dist` are provably non-identifiable
Within this module, `Q` is used **only** to build `L` (EQ-2), and `L` **only ever appears** as the product `k_dist · L` (EQ-6/7). Therefore the PK depends on the triple only through the single combination `k_dist × _PLASMA_CO × _LYMPH_RATIO`.

**PROVEN BY EXECUTION** (not argued): holding the product constant while varying the factors gives a **bit-identical** trajectory —

| CO | k_dist | product | max relative deviation from reference |
|---|---|---|---|
| 5000 | 3.0 | 15000 | 0.000e+00 (reference) |
| 10000 | 1.5 | 15000 | **0.000e+00** |
| 2500 | 6.0 | 15000 | **0.000e+00** |
| 5000 | 6.0 | 30000 | 2.495e-01 (control — the product *did* change) |

So the ~26%-high `_PLASMA_CO = 5000` (§3.2) is **harmless but also meaningless**: any error in it is exactly absorbed by the fitted `k_dist`. The apparent physiological grounding of the flow scale is **illusory** — one fitted number (`k_dist`) is doing all the work. *(Caveat, verified: this degeneracy is scoped to THIS module. In the per-cell path `Q` acquires an independent role in a vascular-mixing QSS, `wholebody_percell.py:168`, which breaks the degeneracy there.)*

### 5.4 🟠 The "16.1 d exactly" calibration claim is dose-conditional
The mosunetuzumab 16.1 d anchor is real (FDA label ✅) but the model reproduces it **only at dose ≈ 695** abstract units (EQ-16). Across the scored dose grid the IgG half-life ranges **9.1 → 19.3 d**, and in the TMDD-free limit it is **38.1 d** — i.e. the model's *intrinsic* (linear) IgG half-life is more than **2× longer** than the anchor, and it is TMDD, not `CLup`, that drags the apparent value down into the clinical range at the fitted dose. The header's separate claim of "IgG t1/2 ~ 2-3 wk" (`:77`) is likewise met only over part of the dose range.

**What to do:** state the dose at which the calibration holds, or re-anchor `CLup` on the TMDD-saturated (linear) half-life so the claim is dose-independent. As written, the comment overstates what was achieved.

### 5.5 🟠 The PK is one-way decoupled from the PD (static receptor field)
`R_CD3 = Rcap_CD3 · PB.tcell` is a **constant array**. The model simultaneously simulates a T-cell population that expands substantially under engagement — and that expansion **does not feed back** into the CD3 TMDD sink. Likewise the TAA sink does not shrink as the tumour is killed. This is internally inconsistent: the model knows the tumour is disappearing (`dTumor`) but the drug keeps being cleared by a tumour-sized antigen sink. For a TCE — where target-cell depletion is precisely the therapeutic effect and is well known to *change the PK over cycles* — this is a real structural gap. It is the single most valuable next improvement to this layer.

### 5.6 🟡 FcRn cannot saturate
The QSS recycled-fraction reduction (EQ-8) makes `fFcRn` a constant. Real FcRn is a **finite, saturable** receptor pool — which is why high-dose IgG (and IVIG) shortens IgG half-life, and why anti-FcRn agents work at all. This model structurally cannot express that. Acceptable for the doses simulated; must be stated.

### 5.7 🟡 Provenance gaps that must be closed before publication
- ~~`qsp_tabula_sapiens_densities.csv` does not exist~~ — **RETRACTED, this was false.** The CSV exists (`deliverables/06_qsp_science-qsp/v2_enhanced/`), and all 30 densities were re-derived from it this pass (§3.4). **No action needed.** *(Residual, minor: the code's `:109` pointer is a bare filename with no path, which is what made it look missing. Worth making the path explicit in the comment.)*
- **σ_V (15 values): no source in code** — with the Tabula Sapiens issue retracted, this is now **the single biggest open provenance gap in T1**: the most influential biodistribution parameters in the layer are unsourced, and they do not reproduce the platform's two-pore split (tight σ₁ 0.883–0.987 / leaky σ₂ 0.311–0.837).
- `fFcRn`, `k_lymph_return`, `hill_renal`, `Rcap_*`, `kint_*` — no sources.
- The named PK comparators in the `sigL` comment ("pembrolizumab 2.17, trastuzumab 2.7, mosunetuzumab 2.1") and the "pembrolizumab day-1/3/7 fall" target for `k_dist` are **[UNVERIFIED CITATION]** — I did not confirm them.
- 🚩 **Contamination note (do not propagate).** The comment at `:988–989` of this file asserts *"clinical 717 > mosun 570"* for IL-6. Per the provenance audit of 2026-07-13, **570 has NO SOURCE** (the only valid IL-6 anchors are mosunetuzumab **152** and teclistamab **21**, both population means). That comment is contaminated. It is in the PD/cytokine block, not the PBPK block, so it does not touch T1's mechanism — but it sits in this file and must not be quoted. See the IL-6 subsystem doc.

### 5.8 🟢 What is solid and should be defended
- `L = Q/500` — a faithful, verified reproduction of the platform's signature relation.
- The Fc/MW format switch (EQ-9/EQ-10) — genuinely mechanistic, and **the BiTE half-life validates** (1.78 h modelled vs 2.11 h ± 1.42 clinical ✅).
- The Tabula Sapiens immune densities — provenance CSV present, all 30 values reproduced from it (§3.4).
- The QSS-vascular lumping (EQ-4) and the `V_c = 3.1 L` decision — correct, well-reasoned, and defensible.
- Emergent, correctly-signed TMDD (EQ-13).

---

## 6. REPRODUCIBILITY

Every number in §2–§5 tagged **RUN-VERIFIED / MEASURED / PROVEN** was produced by executing the live module (`/usr/local/Caskroom/miniconda/base/envs/claude-skills/bin/python`, scipy LSODA `rtol=1e-9, atol=1e-12`, plus an exact `numpy.linalg.eig` decomposition of the linear PK block cross-checked against the ODE solution). Probe scripts:

- `/Volumes/T7_SSD/claude-tmp/2026-07-13/t1_pbpk_probe.py` — prints the frozen arrays + derived scalars (§3.1)
- `/Volumes/T7_SSD/claude-tmp/2026-07-13/t1_pk_validate2.py` — dose-dependence of t½, k_dist cancellation (EQ-16)
- `/Volumes/T7_SSD/claude-tmp/2026-07-13/t1_eigen.py` — eigen-phase decomposition, IgG + BiTE (EQ-15)
- `/Volumes/T7_SSD/claude-tmp/2026-07-13/t1_vss.py` — exact analytic Vss/Vc; dose solving t½=16.1 d (EQ-14, §5.4)
- `/Volumes/T7_SSD/claude-tmp/2026-07-13/t1_ident.py` — the identifiability proof (§5.3)

**Independent adversarial re-run (2026-07-13), which found the defects in §7:**
- `/Volumes/T7_SSD/claude-tmp/2026-07-13/t1_adversarial_verify.py` — re-prints every array, every derived scalar, k_cat/k_renal for both formats, the analytic Vss/Vc at k_dist ∈ {1,3,9}, and an independent eigen-decomposition + LSODA t50 for IgG and BiTE
- `/Volumes/T7_SSD/claude-tmp/2026-07-13/t1_adv_dose.py` — re-runs the full nonlinear dose sweep (EQ-16 table) and re-solves dose(t½ = 16.1 d)
- `/Volumes/T7_SSD/claude-tmp/2026-07-13/t1_adv_ident.py` — re-runs the CO/k_dist degeneracy control **and re-derives all 30 immune densities from `qsp_tabula_sapiens_densities.csv`**

**Everything in §2–§5 reproduced exactly** — the frozen physiology table, Σ Q = 5000.000000000001, Σ Vv = 2.0559 L, Σ Vis = 8.5508 L, k_cat(IgG) = 0.039575/day, k_renal(BiTE) = 8.0097/day, Vss/Vc = 2.0903, the eigen t½ (38.11 d / 5.13 d), the AUC-fraction split (95.8% / 0.34%), the whole EQ-16 dose table (9.11 → 37.81 d; AUC 0.65 → 1.66), dose(16.1 d) = 695.09, and the identifiability control (0.000e+00 / 0.000e+00 / 2.495e-01) — **except the BiTE t50, which is 1.78 h, not 1.80 h.**

**Literature status (be precise about what is and is not in the code):**
| claim | in the code? | checked against source? |
|---|---|---|
| Shah & Betts 2012, PMID 22143261, DOI 10.1007/s10928-011-9232-2 | ✅ yes, `:893` | ✅ confirmed |
| Platform: 15 tissues + carcass + tumor; 5 sub-spaces; explicit FcRn on/off | ✗ no | ✅ confirmed |
| Platform: lymph = 0.2% of plasma flow (`L = Q/500`) | ✅ yes (as `_LYMPH_RATIO`, `:101`) | ✅ confirmed |
| Platform: σ_L = 0.2 | ✗ no | ✅ confirmed |
| Platform: σ₁ 0.883–0.987 / σ₂ 0.311–0.837 | ✗ no | ✅ confirmed (and this **corrects** fabricated ranges in the prior revision) |
| Tabula Sapiens, Science 2022, DOI 10.1126/science.abl4896 | ✅ yes, `:110` | ✅ paper real; **and the 30 derived values reproduce from the in-repo provenance CSV** |
| FDA LUNSUMIO: mosunetuzumab terminal t½ 16.1 d | ✅ yes, `:184–185` | ✅ confirmed |
| FDA BLINCYTO: blinatumomab t½ **2.11 h (SD 1.42)**, MW ~54 kDa | ⚠️ partial — the code says only "~2.1 h" (`:207`) | ✅ confirmed (2.11 / 1.42 is a doc-side lookup, not a code value) |
| "pembrolizumab 2.17, trastuzumab 2.7, mosunetuzumab 2.1" Vss/Vc; "pembrolizumab day-1/3/7 fall" | ✅ yes, in the `sigL`/`k_dist` comments (`:172–180`) | ❌ **NOT checked — [UNVERIFIED CITATION]** |
| tumour antibody-biodistribution coefficient ≈ 0.1–0.3 | ✗ no | ❌ **NOT checked — do not cite from this doc** |

---

## 7. ADVERSARIAL AUDIT LOG — defects found in this doc and fixed (2026-07-13)

This doc was itself audited against the live code. Six defects were found **in the documentation** (none in the model) and are corrected above. They are recorded here rather than silently overwritten, because a doc that has been wrong once must show its scar tissue.

| # | severity | defect | status |
|---|---|---|---|
| 1 | 🔴 **CRITICAL** | §3.4/§4/§5.7 + the At-a-glance asserted that `qsp_tabula_sapiens_densities.csv` **"DOES NOT EXIST"** anywhere in the repo, that a repo-wide `find` returned nothing, and that a grep for the literal values found only coincidental hits. **All false.** The file is at `deliverables/06_qsp_science-qsp/v2_enhanced/`, and all 30 densities re-derive from it exactly (`f_organ / f_spleen`). The doc invented a provenance defect that does not exist and demanded work ("regenerate and commit this CSV") that is not needed. | **RETRACTED + corrected** |
| 2 | 🔴 **CRITICAL** | §3.1 attributed to Shah & Betts the vascular coefficients "σ₁ ≈ 0.69–0.999 (mean 0.908), σ₂ ≈ 0.258–0.841 (mean 0.579)". Those four numbers are in neither the code nor the source — **fabricated**. Published: σ₁ 0.883–0.987, σ₂ 0.311–0.837 (no mean verified, so none is quoted). | **corrected** |
| 3 | 🟠 | BiTE time-to-50%-of-C₀ was stated as **1.80 h** in four places. Re-measured: **1.7830 h** (linear block) / **1.7829 h** (full model). | **corrected to 1.78 h** |
| 4 | 🟡 | Four off-by-one line citations into the module header: FcRn passage `68–72`→**69–72**; ICRP/Brown `:76`→**`:75`**; "4 platform system parameters" `:77`→**`:76`**; "IgG t1/2 ~2-3 wk" `:78`→**`:77`**. | **corrected** |
| 5 | 🟡 | §3.3 (`mw_kda`) and the §2 notation preamble cross-referenced **"EQ-17"**, which does not exist (the doc has EQ-1…EQ-16). The µg/mL↔nM conversion is a readout at `:1279`. | **corrected** |
| 6 | 🟡 | Internal contradiction: the At-a-glance said "**6** are confirmed DEAD" in one bullet and "**5** parameters are dead" in the next. §3.6 lists **6 entries / 9 symbols**. | **corrected** |
| — | ℹ️ | Also newly noted (not a doc defect — a code/comment mismatch): `apply_spatial_exposure` clips to **[0,1]**, not "(0,1]" as both the code comment and this doc claimed (§3.5). | **logged** |

**What survived the audit unchanged:** every equation-to-line-number mapping in §2 (EQ-1…EQ-13 all verified at their cited lines), every parameter value in §3.2/§3.3 (all match the code exactly), every dead-parameter claim in §3.6 (all re-grepped, all still dead), the σ_L = 0.85 / Vss = 2.09 analysis (§5.1), the non-identifiability proof (§5.3), the dose-conditional 16.1 d finding (§5.4), the static-receptor-field finding (§5.5), and the cross-file reach claims (`run_tce_pd_reval.py:6/8/51/142–143`, `wholebody_percell.py:168` — all confirmed at those exact lines). **The model itself was found to be exactly as this doc describes it. The errors were all doc-side.**
