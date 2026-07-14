# MODEL DOCUMENTATION — costim-engager whole-body QSP PBPK-PD (unified, validated)

Two complementary documentation systems, both validated against the LIVE engine source
(cross-checked: both independently identify the same NON-RUNTIME modules and the SAME live binding
path = kinetic_rhoden_percell + multiarm_binding + kinetic_synapse). IN-USE FINAL components only.
NOTE: NOT_IN_USE_REGISTER separates 5 DEAD modules (cytokine_pbpk, il6_pbpk, unified_binding,
multiarm_kinetic, biexact_solver — never wired, do not cite) from 3 OFFLINE-but-LEGITIMATE build-time
modules (rna_to_receptor, convert_copies_ALL, calib_kdeath — real, correct, in-use provenance/build code,
NOT runtime and NOT dead). Do not conflate the two.

## PRIMARY SYSTEM — subsystem deep-docs (start here)
Authoritative, file:line-cited, adversarially reviewed (each T-doc carries its own defect log).

  00_START_HERE.md              — 8 invariants, each earned by a real failure. Read first.
  MASTER_MODEL_DOCUMENTATION.md — L0–L9 architecture map.
  TCE_LIFECYCLE_WALKTHROUGH.md  — life of the molecule: Shah-Betts reproduction -> production.
  NOT_IN_USE_REGISTER.md        — 5 DEAD modules (never wired, do not cite) + 3 OFFLINE build-time
                                  modules (rna_to_receptor/convert_copies_ALL/calib_kdeath — real,
                                  correct, legitimately in-use at build time, just not runtime).
  PARAMETER_AUDIT_2026-07-13.md — every parameter: FITTED vs MEASURED vs DERIVED, with source.
  PROVENANCE_AND_VALIDATION.md  — validation record + provenance.
  RESULTS_SUMMARY.md            — headline results.
  subsystems/
    T1_shah_betts_pbpk_backbone.md        — 2-pore + FcRn PBPK backbone (the foundation layer)
    T2_whole_body_per_cell_pk.md          — per-cell interstitial transport (replaces well-mixed IS)
    T2b_receptor_overlay_imputation.md — **LOAD-BEARING**: scVI Tabula-Sapiens overlay that imputes every receptor onto every cell (the step that makes the model runnable) + RNA→copies conversion + retention gap.
    T3_rhoden_bivalent_binding_core.md    — Rhoden kinetic bivalent binding core
    T4_multi_arm_format_geometry.md       — 0/1/2 valency per arm, cis/trans spans, tetravalent max
    T5_kinetic_immune_synapse.md          — B1/B2 synapse, K_HIT, cleft geometry
    T6_per_cell_pd_killing.md             — per-cell kill law
    T7_costim_signaling_activation_induced_receptor_density.md — costim signaling layer
    T8_mechanistic_crs_il_6.md            — myeloid IL-6 / CRS mechanism
    T9_integration_driver.md              — the integration driver tying it together

## COMPLEMENTARY CHAPTERS — narrative "life of the molecule" + report
  complementary_chapters/DOC_01..04       — narrative walkthrough (overview/PK/binding-PD/IL6-axes)
  complementary_chapters/QSP_REPORT_SECTION.md — publication Results+Methods

## EVIDENCE / PROVENANCE
  evidence/QSP_VERIFICATION_LEDGER.csv    — one row per load-bearing claim, PASS/FLAG + source
  evidence/IL6_VALIDATION_RESULT.json     — IL-6 gate: FAIL on 3 verified anchors (talq inverts)
  evidence/IL6_ANCHORS_VERIFIED_2026-07-13.md — sourced clinical IL-6 anchors
  evidence/CHANGELOG_2026-07-13.md        — dated build history

## HONEST STATUS (do not overclaim)
- IL-6 validation FAILS on the 3 manifest-verified anchors (mosun/glofit/talq); talquetamab inverts
  (GPRC5D on plasma/keratin, not myeloid-contacting -> compartment gap). teclistamab 21 is NOT a verified
  anchor and any tecli-based test is PROVISIONAL. Nomination is independent of this layer (liability veto upstream).
- cytokine_pbpk.py / il6_pbpk.py are DEAD — the model does NOT give each cytokine its own PBPK compartment.
- Two IL-6 pre-registration files (endpoint + prediction) were deliberately EXCLUDED from this folder.
- Some parameters are FITTED (kdeg IL-6, sigma_L) — see PARAMETER_AUDIT; not presented as measured.

## INTEGRATION NOTE
- complementary_chapters/DOC_03 has a CITATION CAVEAT banner: its narrative is correct but some
  wholebody_pd.py line numbers drifted (its reviewer flagged p_eng cited at :485 = wrong line).
  The subsystems/T5,T6 docs are authoritative for binding+kill file:line citations (snapshot-verified).
