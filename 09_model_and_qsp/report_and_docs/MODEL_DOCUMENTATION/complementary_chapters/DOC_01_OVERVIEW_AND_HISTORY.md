# The QSP Costim-Engager Counter-Screen Model — Overview & Build History

**Document 1 of the model documentation system.** This chapter is the spine: it
introduces the model, traces the *life of a molecule* through it end-to-end, and
narrates how the model was built from a reproduction of the Shah & Betts (2012)
platform PBPK up to the current whole-body single-cell PBPK-PD system. The three
companion chapters give the deep mechanistic/mathematical treatment of each arm:

- **DOC_02_PK_DISTRIBUTION.md** — plasma → tissue → per-cell interstitium → recirculation, and the per-cell Rhoden TMDD sink.
- **DOC_03_BINDING_PD.md** — synapse formation, serial killing, effector/suppression/exhaustion programs.
- **DOC_04_IL6_CRS_AXES.md** — mechanistic myeloid IL-6/CRS and the three-axis counter-screen.

> **Scope.** This documentation covers only the **in-use final execution path** — the
> 12 engine files that `run_tce_pd_reval.py` actually loads at runtime (verified by an
> import trace, recorded in `SUBMISSION_MANIFEST.json → LIVE_FINAL_PATH`). Retired and
> experimental files (`unified_binding.py`, `biexact_solver.py`, `multiarm_kinetic.py`,
> `il6_pbpk.py`, `cytokine_pbpk.py`, `convert_copies_ALL.py`, `rna_to_receptor.py`,
> `calib_kdeath.py`) are **not** part of the live model and are excluded.

---

## 1. What the model is, in one paragraph

The model is a whole-body physiologically-based pharmacokinetic / pharmacodynamic
(PBPK-PD) system for T-cell engagers (TCEs) and costimulatory bispecifics, in which
the classical *well-mixed interstitial compartment* of a Shah-Betts PBPK is **replaced,
per organ, by a spatial single-cell agent layer**: every cell from a Xenium/scRNA-seq
overlay is an individual agent carrying its own receptor copy numbers, running its own
Rhoden binding kinetics, at its own grid location. Plasma drug feeds each tissue through
its blood endothelium (BEC), diffuses and binds cell-by-cell across the interstitial
grid, and drains back through the lymphatic endothelium (LEC) into a systemic lymph pool
that recirculates to plasma. Target-mediated drug disposition (TMDD), cytotoxic synapse
formation, cell killing, and cytokine (IL-6/CRS) release are therefore all **emergent
from per-cell binding on real tissue geometry**, not imposed by compartment-level rate
constants. The same binding engine runs in the PK and the PD arms, so the two are
mechanistically identical by construction.

## 2. The life of a molecule (end-to-end trace)

1. **Administration.** A dose is delivered as an IV bolus, a timed IV infusion, or a
   subcutaneous depot (bioavailability `F_sc`, absorption rate `ka`). Dosing is
   mass-exact — the integral of the infusion equals `mg/mw·10³` regardless of step size.
   *(DOC_02 §Dosing.)*
2. **Plasma & whole-body distribution.** The molecule distributes across 15 tissues
   (+tumor), each a vascular + interstitial pair connected by blood flow `Q` and lymph
   flow `L = Q/500` (Shah-Betts). It leaves the blood by **2-pore convective
   extravasation** governed by the vascular reflection coefficient `σ_V`. *(DOC_02 §PBPK.)*
3. **Entry into tissue (BEC).** In a spatially-resolved organ, drug enters the
   interstitial grid at blood-endothelial cells and spreads by diffusion. *(DOC_02 §Transport.)*
4. **Per-cell binding (the sink).** At every cell, the drug engages CD3 (T side) and/or
   TAA (target side) and/or costim, by **Rhoden bivalent kinetic binding** with
   receptors-as-states, free-receptor turnover (`KSYN = R_tot·k_deg`), and internalization
   (`k_int`). The bound-complex internalization is the physical origin of TMDD. *(DOC_02 §TMDD.)*
5. **Synapse & killing.** Where a T cell and a target cell are within synapse range
   (30 µm incidence graph `W`), the CD3·drug·TAA trimer forms (`B2`), gated by arm-span
   geometry (bridge/cis spans, cleft feasibility 13–40 nm) and cis/trans coincidence
   (`p_cis`). An engaged T cell delivers lethal hits at a serial rate capped by
   `k_hit = 12/day`; killing removes TAA and frees the synapse for re-bridging. *(DOC_03.)*
6. **Effector programs.** Costim occupancy `occ = (C_d/(C_d+K_D))·(R_c/anchor)` drives
   signaling programs; the effector program multiplies kill by `g_eff = exp(0.55·eff_p)`
   on CD8 T cells, exhaustion attenuates it, and Treg suppression damps the neighborhood.
   *(DOC_03 §PD.)*
7. **Cytokine / CRS.** Engaged myeloid cells emit IL-6 per-cell (Giavridis/Norelli
   biology); IL-6 accumulates in an extracellular-fluid (ECF) plasma pool with
   first-order clearance. CRS magnitude is emergent. *(DOC_04.)*
8. **Drainage & recirculation (LEC).** Unbound drug drains through lymphatic endothelium,
   distance-graded, into the systemic lymph pool and returns to plasma; catabolism
   `k_cat = CL_up·(1−fFcRn)` sets the terminal half-life, with FcRn salvage as a
   quasi-steady recycled fraction. *(DOC_02 §Recirculation.)*

## 3. Build history — from Shah & Betts to now

The model was built by reproducing an established platform and then progressively
replacing its lumped compartments with mechanistic, data-grounded single-cell layers.

**Origin — Shah & Betts (2012) platform PBPK.** The whole-body skeleton is the
Shah & Betts 2-pore + FcRn platform: 15 reference-human tissues (ICRP/Brown 71 kg),
each vascular + interstitial, connected by blood flow `Q` and lymph `L = Q/500`, with
antibody extravasation by 2-pore convection and FcRn salvage giving IgG t½ ≈ 2–3 weeks.
(`qsp_costim_window_v2.py`, LAYER 1.) Tissue T-cell and myeloid densities were later
grounded in the Tabula Sapiens single-cell atlas (Science 2022) rather than left as priors.

**Interstitium → per-cell spatial layer.** The well-mixed interstitial compartment was
replaced, per organ, by a spatial agent grid built from Xenium morphology + scRNA-seq
overlays (CRC tumor first, then 12 organs + additional solid tumors + heme + blood).
Each barcode became an agent with its own receptor copies (HPA-IHC-anchored Glassman
conversion) and its own Rhoden binder. Transport was wired blood-in (BEC) / lymph-out
(LEC) with an implicit-Euler diffusion solver on a 100 µm grid.

**PK/PD binding unification.** The PK TMDD sink and the PD kill synapse were unified onto
one kinetic binding engine (kon/koff, not QSS), with free-receptor turnover
(`KSYN = R_tot·k_deg`) everywhere and multi-arm valency (CD3/costim/TAA each 0–2, up to
tetravalent). A backward-Euler receptors-as-states solver made the bivalent path stable
and fast. This is the "binding is identical in PK and PD" invariant.

**Validation.** PK was validated against digitized clinical concentration-time curves
(teclistamab SC, AFE 1.29×; elranatamab matched PK+IL-6, AFE 2.06×). PD kill timescale
was anchored to a shared `k_death = 1.0` with `k_hit = 12/day` fixed from serial-killing
literature. IL-6 emission was moved from a fitted scale to a **mechanistic per-cell
myeloid model** and scored as a between-molecule ratio against digitized clinical anchors.

**Counter-screen nomination.** On top of the mechanistic model sits the three-axis
costim nomination: effector benefit (Schmidt CRISPRa CD8 IFN-γ), suppression liability
(CD4 IL-10/Treg program), and CRS liability (TNF/IL-2/IFN-γ). A 6-axis liability veto
gates candidates *upstream* of the QSP window; the co-leads are **4-1BB (TNFRSF9)** and
**CD27**.

*(A dated, blow-by-blow record of every fix in the final build day is in
`CHANGELOG_2026-07-13.md`; the frozen headline numbers and their source artifacts are in
`SUBMISSION_MANIFEST.json`.)*

## 4. Honest limitations (carried into every chapter)

- **Static R_costim.** Costim receptor copies are set once from *resting* numbers and read
  unchanged each step. The model captures tumor-conditionality through binding geometry
  (the cis gate) but **not** through activation-induced receptor upregulation. A resting-copy
  ranking therefore under-rates 4-1BB/OX40/ICOS and can yield a spurious "CD2 wins" ordering
  on the effector axis — which is exactly why the liability veto, not the effector drive,
  decides the nomination.
- **CRISPRi is loss-of-function; an engager arm is gain-of-function.** No screen can directly
  prove agonism helps; the screen supplies a validated state change and the QSP model
  translates it into a predicted window.
- **CD4 is the counter-screen, not the effector axis.** CD8 does the killing; the CD4 screen
  is used only for what it uniquely resolves — the suppressive and cytokine-release programs.
- **Provisional results are excluded or labeled.** The OX40/GITR net-negative-kill-in-Treg-rich
  figure is provisional (structurally impossible under static R_costim) and is omitted from
  the validated deliverables.
