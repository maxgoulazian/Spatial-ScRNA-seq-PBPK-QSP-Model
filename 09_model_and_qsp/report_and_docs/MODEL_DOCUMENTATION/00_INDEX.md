# 00 — INDEX

**TCE Spatially-Resolved Counterscreen — complete documentation set**

---

## Read in this order

| # | document | what it gives you |
|---|---|---|
| 1 | **[`00_START_HERE.md`](00_START_HERE.md)** | Orientation, the **8 non-negotiable invariants**, glossary. *Every invariant was earned by a real failure.* |
| 2 | **[`NOT_IN_USE_REGISTER.md`](NOT_IN_USE_REGISTER.md)** | **8 of the 20 modules in `engine/` do not run.** Read before you document, cite, or trust any code. |
| 3 | **[`PROVENANCE_AND_VALIDATION.md`](PROVENANCE_AND_VALIDATION.md)** | What is measured, what is fitted, **what was fabricated**. Read before you trust any number. |
| 4 | **[`MASTER_MODEL_DOCUMENTATION.md`](MASTER_MODEL_DOCUMENTATION.md)** | L0–L9: scope · parameters · equations · modules · subsystems · dataflow · provenance · known-unknowns · validation status · update protocol. |
| 5 | **[`TCE_LIFECYCLE_WALKTHROUGH.md`](TCE_LIFECYCLE_WALKTHROUGH.md)** | **The life of the molecule.** Syringe → plasma → extravasation → first bond → avidity → synapse → kill → myeloid contact → IL-6 → blood test. Opens with the Shah-Betts reproduction and the build history. **Start here if you want to *understand* the model rather than operate it.** |

---

## The subsystems (deep docs)

Each carries: purpose & dataflow position · **every governing equation** (with `file.py:line`, biological
meaning, mechanistic rationale, units) · **every parameter** (value, units, provenance tag, source, rationale)
· what is **emergent vs imposed** · known limitations.

| | subsystem | live code | doc |
|---|---|---|---|
| **T1** | Shah-Betts PBPK backbone | `qsp_costim_window_v2.py` | [`subsystems/T1_shah_betts_pbpk_backbone.md`](subsystems/T1_shah_betts_pbpk_backbone.md) |
| **T2** | Whole-body per-cell PK | `coupled_percell_pk.py`, `wholebody_percell.py` | [`subsystems/T2_whole_body_per_cell_pk.md`](subsystems/T2_whole_body_per_cell_pk.md) |
| T2b_receptor_overlay_imputation.md — **LOAD-BEARING**: scVI Tabula-Sapiens overlay that imputes every receptor onto every cell (the step that makes the model runnable) + RNA→copies conversion + retention gap.
| **T3** | Rhoden bivalent binding core | `kinetic_rhoden_percell.py` | [`subsystems/T3_rhoden_bivalent_binding_core.md`](subsystems/T3_rhoden_bivalent_binding_core.md) |
| **T4** | Multi-arm format geometry | `multiarm_binding.py` | [`subsystems/T4_multi_arm_format_geometry.md`](subsystems/T4_multi_arm_format_geometry.md) |
| **T5** | Kinetic immune synapse | `kinetic_synapse.py` | [`subsystems/T5_kinetic_immune_synapse.md`](subsystems/T5_kinetic_immune_synapse.md) |
| **T6** | Per-cell PD — killing | `wholebody_pd.py`, `pd_model_config.py` | [`subsystems/T6_per_cell_pd_killing.md`](subsystems/T6_per_cell_pd_killing.md) |
| **T7** | Costim signaling & **activation-induced** receptor density | `costim_induction.py` | [`subsystems/T7_costim_signaling_activation_induced_receptor_density.md`](subsystems/T7_costim_signaling_activation_induced_receptor_density.md) |
| **T8** | Mechanistic CRS IL-6 | `myeloid_il6.py` | [`subsystems/T8_mechanistic_crs_il_6.md`](subsystems/T8_mechanistic_crs_il_6.md) |
| **T9** | Integration & driver | `coupled_percell_pd.py`, `run_tce_pd_reval.py` | [`subsystems/T9_integration_driver.md`](subsystems/T9_integration_driver.md) |

---

## The five things a reader most needs to know

1. **The model's selling point is emergence.** Prozone, CRS saturation, tumour-conditional costim, and the
   per-molecule differences are **computed**, not fitted. That is what licenses ranking constructs that have
   never been built. → `T5`, `T7`, `T8`

2. **Efficacy and toxicity come out of the same geometry.** A T cell engaged beside a tumour cell kills it;
   the same T cell engaged beside a macrophage causes CRS. **Therapeutic window is one mechanism read two
   ways.** → `TCE_LIFECYCLE_WALKTHROUGH.md` §10

3. **The linker prediction is counter-intuitive and falsifiable.** Effective 2nd-arm concentration falls as
   `1/span`; cleft feasibility is zero below ~0.6× the cleft. Their product peaks at **span = cleft ≈ 13 nm** —
   so **past the cleft distance, longer linkers make bridging *worse*.** → `T4`

4. **Two constants carry the model's scales, and both are FITTED.** `k_death` (efficacy, calibrated to
   epcoritamab depletion) and `KDEG_IL6_PER_HR` (the entire absolute IL-6 scale, whose in-code citation is to
   a *modeling* paper). **Human IL-6 clearance appears to be unmeasured.** → `PROVENANCE_AND_VALIDATION.md` §2

5. **The IL-6 arm is not yet validated, and the costim screen is blocked.** Stated plainly in
   `MASTER_MODEL_DOCUMENTATION.md` §L8 rather than buried. The decisive IL-6 test is a single number
   (**7.24×**); the costim blocker is four unsourced fold-upregulations, and the code **fails closed** rather
   than guess them.

---

## Regenerating / maintaining

- **Live/dead module split:** `NOT_IN_USE_REGISTER.md` §3 (a script — do not trust a stale list).
- **Doc drift:** any mechanistic change **must** update its `subsystems/TX_*.md` in the **same** change.
  A drifted doc is worse than no doc: it is a confident lie.
- **Provenance changes** go in `PROVENANCE_AND_VALIDATION.md`, **not** just a code comment. *The code comment
  is where the last fabrication hid.*
