# MASTER MODEL DOCUMENTATION
## TCE Spatially-Resolved Counterscreen — authoritative, provenance-tagged, multi-level reference

> **Read `00_START_HERE.md` first.** It carries the eight invariants; this document assumes them.
>
> **Levels.** L0 scope · L1 parameters · L2 equations · L3 modules · L4 subsystems · L5 integration &
> dataflow · L6 provenance ledger · L7 known-unknowns · L8 validation status · L9 update protocol.

---

## L0. SCOPE & PROVENANCE

**What the model is.** A whole-body, per-cell, mechanistic simulation of a T-cell engager: a Shah-Betts PBPK
backbone carrying antibody to organs whose interstitia are populated by **individual cells** with **scRNA-seq-
derived receptor copy numbers**, bound by **literal kinetic multivalent chemistry**, producing **emergent**
killing and **emergent, contact-gated, myeloid-derived** cytokine release.

**What it is for.** A **counterscreen**: rank construct designs (TAA × costim arm × affinity × valency ×
linker length × architecture) by therapeutic window — *including constructs that have never been built.*

**What it is NOT.**
- It is **not** validated on absolute IL-6 (see L8, and `PROVENANCE_AND_VALIDATION.md` §2.1 — absolute IL-6 is
  clearance-limited and clearance is fitted, so absolute concentration **cannot** validate the mechanism).
- It does **not** give each cytokine its own PBPK compartment set. `cytokine_pbpk.py` exists and **is not
  wired in** (`NOT_IN_USE_REGISTER.md` §2a).
- Its **kill-rate constants are fitted**, not measured (L6).

**Code pin.** The live path is 12 of 20 modules in `engine/`. Re-derive it — do not trust this list —
per `NOT_IN_USE_REGISTER.md` §3.

---

## L1. PARAMETERS

Every parameter carries exactly one provenance tag: `[MEASURED]` · `[DERIVED]` · `[FITTED]` · `[ASSUMED]` ·
`[UNSOURCED — TBD]`. Definitions in `PROVENANCE_AND_VALIDATION.md` §0.

The **complete** parameter tables live in the subsystem docs (`subsystems/T1…T9`, §3 of each). The
**load-bearing** ones — the parameters whose value most changes the answer — are:

| parameter | value | tag | why it matters |
|---|---|---|---|
| `k_death`, `k_hit` | see `pd_model_config.py` | **`[FITTED]`** — calibrated to the epcoritamab depletion time-course (`calib_kdeath.py`, offline) | The efficacy scale. Everything downstream of the kill law inherits this tag. |
| `KDEG_IL6_PER_HR` | `0.20 /hr` | **`[FITTED]`** — its in-code citation (PMID 31268236) is to a **modeling** paper | With `V`, sets the **entire absolute IL-6 scale**. Human IL-6 clearance appears **unmeasured**. |
| `V_PLASMA_ML` | `11 650 mL` (ECF) | `[ASSUMED — physical]` — IL-6 is 21 kDa, made in the interstitium, has no FcRn | Together with `k_deg` gives implied `CL = 55.9 L/day`. |
| `R_CONTACT_UM` | `14.1 µm` | `[DERIVED]` — r_macrophage 10.6 µm + r_Tcell 3.5 µm | CD40L–CD40 is membrane-bound → cells must **touch**. (30 µm — a *synapse reach* — was wrong.) |
| `SECRETOR_FRACTION` | `0.039` | `[MEASURED]` — PMID 37533643 | Cell-intrinsic; measured under maximal stimulation, so **not** something the spatial model re-derives. |
| `S_MAX` | `10.6 molec/s` | `[MEASURED]` — mean of an actively-secreting monocyte | An earlier `156` was a **peak-tail**, applied to every cell → **14.7× too high**. |
| `COSTIM_INDUCTION[*].fold` | — | **`[UNSOURCED — TBD]`** | **BLOCKS the costim screen.** Code **fails closed** rather than run inducible arms at resting density. |
| `RCAP_CD3`, `RCAP_TAA` | `2.0`, `6.0` | see `T6` | Receptor-capacity scaling onto the validated tumour basis. |
| Shah-Betts `Q, L, Vis, Vv, σ_V, σ_L` | see `T1` | `[MEASURED]` (Shah & Betts) | The physiology the model was built by **reproducing**. |

**Clinical IL-6 anchors — the only two that survived provenance audit:**

| molecule | value | statistic | source |
|---|---|---|---|
| teclistamab | **21** pg/mL | population **MEAN** peak Cmax | PMID 38831634 |
| mosunetuzumab | **152** pg/mL | population **MEAN** peak Cmax | — |

→ **target ratio 7.24×**. **Elranatamab has no clinical IL-6 value in existence.** Everything else previously
used (`570`, `340`, `230`, `366.88`, `191`, `288`, `18.2`) was fabricated, unsourced, or a **different
statistic** — including one that was **a page number**. See `PROVENANCE_AND_VALIDATION.md` §1.

---

## L2. EQUATIONS

Full derivations, biological meaning, and mechanistic rationale for **every** equation are in the subsystem
docs. The equations that carry the model's scientific claims:

**Vascular quasi-steady state** (`T2`) — forced, not chosen: `Q/Vv ~ 10⁵/day` vs slow tissue dynamics makes
explicit integration catastrophically stiff.

    C_vasc = (Q·C_plasma + PS·C_is) / (Q + PS)

**Rhoden bivalent binding** (`T3`) — first and second bonds solved explicitly, per cell. A **backward-Euler**
inner solver is used because the explicit one, with *measured* `kon`, pinned the substep count at 5651 and made
a single step take 7.83 s (BE: 0.026 s, and flat in `kon`).

**Effective 2nd-arm concentration & the linker optimum** (`T4`):

    C_eff ∝ 1/span                      (dilution — a longer linker searches a bigger shell)
    feasibility(span) = 0 below ~0.6·cleft ; → 1 at the cleft
    ⇒ product maximised at span = cleft ≈ 13 nm

> **Prediction: past the cleft distance, LONGER LINKERS REDUCE BRIDGING.** Reach was never the constraint;
> *concentration* was.

**Ternary equilibrium & emergent prozone** (`T5`, `T6`): the trimer term is quadratic in free arms, so at high
drug — when both arms saturate *separately* — bridging collapses **on its own**. **There is no prozone
equation.** This is what licenses extrapolation to unbuilt formats; an imposed hook would need re-fitting for
each one.

**Activation-induced costim** (`T7`):

    da_i/dt = k_on·p_eng_i·(1 − a_i) − k_off·a_i        per-T-cell activation memory
    R_i(t)  = R_rest,i · (1 + (FOLD − 1)·a_i)

`p_eng_i` is **that cell's own** engaged fraction → tumour-conditionality is **emergent**. `FOLD = 1`
reproduces the old static behaviour exactly (constitutive arms are unaffected).

**Mechanistic CRS IL-6** (`T8`): per-macrophage, contact-gated (`R_contact = 14.1 µm`), intrinsic-secretor-
gated, then

    dC/dt = production/V − k_deg·C          (solved analytically, not by explicit Euler)

---

## L3. MODULES (live only)

`run_tce_pd_reval` · `qsp_costim_window_v2` · `pd_model_config` · `coupled_percell_pd` · `coupled_percell_pk` ·
`wholebody_percell` · `wholebody_pd` · `kinetic_rhoden_percell` · `kinetic_synapse` · `multiarm_binding` ·
`myeloid_il6` · `costim_induction`

**8 modules in `engine/` do NOT run.** → `NOT_IN_USE_REGISTER.md`

---

## L4. SUBSYSTEMS

| | subsystem | code |
|---|---|---|
| `T1` | Shah-Betts PBPK backbone | `qsp_costim_window_v2` |
| `T2` | Whole-body per-cell PK | `coupled_percell_pk`, `wholebody_percell` |
| `T3` | Rhoden bivalent binding core | `kinetic_rhoden_percell` |
| `T4` | Multi-arm format geometry | `multiarm_binding` |
| `T5` | Kinetic immune synapse | `kinetic_synapse` |
| `T6` | Per-cell PD — killing | `wholebody_pd`, `pd_model_config` |
| `T7` | Costim signaling & induction | `costim_induction` |
| `T8` | Mechanistic CRS IL-6 | `myeloid_il6` |
| `T9` | Integration & driver | `coupled_percell_pd`, `run_tce_pd_reval` |

---

## L5. INTEGRATION & DATAFLOW — end to end

    dose (IV bolus | SC: F_sc, ka)
      → plasma
      → per organ: vascular pool  [QUASI-STEADY]
      → extravasation (σ_V)  →  interstitium        ← where the cells actually are
      → per-cell drug exposure
      → first bond (kon/koff, measured where available)
      → second bond: avidity, format geometry, linker span            [T4]
      → trimer  →  IMMUNE SYNAPSE  (prozone EMERGES)                  [T5]
      ├─→ kill hazard  →  target depletion  (Treg-damped)             [T6]
      ├─→ costim occupancy → per-cell programs → activation-induced R  [T7]
      └─→ engaged T cell TOUCHES a macrophage (≤14.1 µm)              [T8]
             → that macrophage secretes (if an intrinsic secretor)
             → Σ over organs, lifted by the MYELOID census
             → plasma IL-6 ODE (ECF volume, first-order clearance)
             → **the clinical readout**

**The single most important structural fact:** efficacy and toxicity emerge from **the same geometry**. A T
cell engaged beside a tumour cell kills it; the same T cell engaged beside a macrophage causes CRS.
**Therapeutic window is one mechanism read two ways** — which is exactly what makes a counterscreen possible.

**Two traps in the driver** (`T9`):
- **Myeloid `count_scale` is a TISSUE property** (drug-independent). The target-cell `count_scale` is
  **drug-dependent** (it is `total_antigen/sampled_antigen`, and varies 0.18×–5.34× between antigens in the
  same organ). Using the latter for myeloid would scale the *same monocytes* differently for different drugs —
  **a per-drug artifact that would corrupt the counterscreen outright.**
- **Step-up ladders.** `TSIM_DAYS=7` gives mosunetuzumab **only its 1 mg priming dose** (60 mg is on day 15)
  while giving teclistamab its **full 120 mg** (day 7). **Use `TSIM_DAYS=24`.**

---

## L6. PROVENANCE LEDGER

See `PROVENANCE_AND_VALIDATION.md` — the authoritative record. Summary of what is **not** measured:

- **`k_death` / `k_hit`** — `[FITTED]` to epcoritamab depletion. The efficacy scale.
- **`KDEG_IL6_PER_HR`** — `[FITTED]`; its citation is to a **modeling** paper. The IL-6 scale.
- **`COSTIM_INDUCTION.fold`** — `[UNSOURCED]`. **Blocks the costim screen.**
- **Liver Kupffer census** — could not be sourced → **liver excluded** from the myeloid census → **biases IL-6
  downward.** Stated, not hidden.

---

## L7. KNOWN-UNKNOWNS

1. **Human IL-6 clearance is unmeasured.** The single largest assumption in the model.
2. **Costim induction folds are unsourced.** The costim screen cannot produce a valid ranking until they are.
3. **Liver Kupffer count** — unsourced; liver excluded; IL-6 biased low.
4. **Several CD3 affinities are class-estimates**, not molecule-specific SPR. Tagged in `T6`.
5. **Does the model reproduce step-up protection?** Clinically, **CRS is worst at the FIRST dose** — that is
   why step-up dosing exists. Whether this model reproduces that (rather than IL-6 simply climbing with dose)
   is **an open question and a first-order test of its validity.**

---

## L8. VALIDATION STATUS — stated honestly

| arm | status |
|---|---|
| **PK** (clinical plasma, 22 molecules) | Independent of IL-6. See `FIGURES_AND_VALIDATION.md`. |
| **IL-6 / CRS** | **NOT VALIDATED.** The decisive test — model vs clinical **7.24×** — is defined and running but has **not returned a trustworthy result.** Every prior IL-6 "validation" was scored against contaminated anchors and is **void**. |
| **Costim** | **BLOCKED** on unsourced induction folds. Code fails closed. |

**Runs that must NOT be trusted** (they look good and are wrong):
- **Pre-census runs** (`tce_pd_rd_*`) — `count_scale` was 1.0, IL-6 ~290,000× too low. Teclistamab ≈ 25 vs a
  clinical 21 there is **a bug that looks like a triumph.**
- **Any `TSIM_DAYS=7` run** — dose-stage mismatch (above).

---

## L9. UPDATE PROTOCOL

1. Any mechanistic change **must** update the matching `subsystems/TX_*.md` **in the same change.** A drifted
   doc is worse than none — it is a confident lie.
2. Re-derive the live/dead module split (`NOT_IN_USE_REGISTER.md` §3) **before** editing docs.
3. When a parameter's provenance changes, update `PROVENANCE_AND_VALIDATION.md` — **not just the code comment.**
   The code comment is where the last fabrication hid.
4. **Never** add a fallback that substitutes a constant for a failed calculation. Crash instead.
