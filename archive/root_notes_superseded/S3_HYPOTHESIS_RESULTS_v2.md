# S3 — Hypothesis Test Results v2 (scoped re-run)
_Module S3' (v2). Generated 2026-07-10T02:27:03Z. Scoped refresh: 12 verdicts CARRIED VERBATIM from v1 (costim/A-track inputs unchanged), H07 re-ingested (revised leads), H15 computed (new). Every statistic computed/read from a named source this run; no fabrication._

## Summary
- **14 testable-now hypotheses** (H02–H15). H01 is bench-only (`testable_now=no`), carried, NOT tested.
- **12 supported · 2 inconclusive · 0 refuted.** The two inconclusive (H08, H10) are v1-carried, limited by data coverage / power, not by contradicting evidence.
- **Change tags:** 12 carried verbatim, 1 revised (H07), 1 new (H15).
- **Multiplicity:** the BH-FDR discovery family is the k=6 affirmative-direction p-value tests (H02,H03,H05,H06,H09,H11) — all CARRIED from v1, all q<0.05 (family max q=1.01e-05, H02). H07-v2 and H15 are threshold/categorical (no primary p) and are correctly OUTSIDE the BH family, so the family and its correction are unchanged from v1.
- **v2 provenance:** H07 re-ingests B3-v2 healthy breadth (leads now CEACAM6/CEACAM5; v1 co-lead LY6E dropped → H15). H15 ingests B8 per-cell-type expression (the S3 'second data wave' for the TAA-breadth question) — verdict on the replicated off-tumor flag, not a recomputation.

## Verdict table
| H | Category | Verdict | tag | Primary statistic | value | p | BH q | Evidence (one line) |
|---|---|---|---|---|---|---|---|---|
| H02 | mechanistic | ✅ supported | carried | CD27 SUPP agonism within FOXP3-negative cells, Stouffer p (deepen, 4-donor) | -0.0998 | 1.01e-05 | 1.01e-05 | CD27 SUPP−negation in FOXP3⁻ cells (−0.0998), Treg fraction unchanged (OR 1.014) → spares Treg |
| H03 | mechanistic | ✅ supported | carried | OX40 SUPP agonism Stouffer p (deepen 4-donor Stim48hr) | 0.1386 | 1.67e-09 | 3.34e-09 | OX40 CRS ns but SUPP up (+0.139) → cytokine-only screen misses it |
| H04 | mechanistic | ✅ supported | carried | 4-1BB mean NF-kB signalosome cosine (A4); shared_signed CRS-axis rank | -0.0728 | — | — | 4-1BB NF-κB cosine −0.073, all q ns → orthogonal, not CD28 phenocopy |
| H05 | mechanistic | ✅ supported | carried | CD30 HELP agonism Stouffer p at Stim48hr (deepen 4-donor) | -0.1163 | 8.22e-07 | 1.23e-06 | CD30 HELP flat early → sig-toxic at 48h (−0.116); late help-erosion disqualifies |
| H06 | safety_resistance | ✅ supported | carried | CD28 dual-liability: max(CRS_p, SUPP_p) at Stim48hr (both must be significant) | 0.1518 | 2.01e-06 | 2.41e-06 | CD28 CRS up AND SUPP up → dual liability (TGN1412 TRUE-NEG) |
| H07 | safety_resistance | ✅ supported | revised | CEACAM6 vital_max vs v2 TAA-panel median (Tabula Sapiens healthy breadth, B3-v2) | 0.7563 | — | — | CEACAM6 vital_max 0.756 > v2 panel median 0.395 (gating needed); CEACAM5 0.099 + 0 copies/cell → escape-resistant. LY6E dropped→H15 |
| H08 | safety_resistance | ⚠️ inconclusive | carried | Shifrut costim-arm coverage (n SUPP-up arms measurable) + A11 contradiction count | 0 | — | — | only 4-1BB of panel in Shifrut library → direct-killing claim untestable; finalists not contradicted |
| H09 | patient_stratification | ✅ supported | carried | 4-1BB (CD137 ADT + TNFRSF9) enrichment in dysfunctional-high CD8, Cliff's delta / Stouffer p | 0.4356 | 1.62e-304 | 9.74e-304 | 4-1BB (CD137 ADT + TNFRSF9) enriched in dysfunctional vs naive CD8 (d=+0.44) both samples |
| H10 | patient_stratification | ⚠️ inconclusive | carried | baseline suppressive-tone heterogeneity (Treg fold-range) + tone->induction Spearman (n=4, underpowered) | 4.719 | — | — | Treg tone 4.7× range (premise ✓) but n=4 tone→induction underpowered (p≥0.2) |
| H11 | mechanistic | ✅ supported | carried | DNAM1 within-donor SUPP agonism Stouffer p (deepen 4-donor Stim48hr) | 0.1478 | 1.19e-19 | 3.57e-19 | DNAM1 surface Treg-spared but functional SUPP up (+0.148) → surface≠function |
| H12 | mechanistic | ✅ supported | carried | Pearson r(effector, CRS agonism) across panel | 0.2953 | 0.327 | — | E~CRS r=0.30 ns, E~SUPP r=−0.19 ns → effector-tox decoupled (null-confirmation) |
| H13 | de_novo | ✅ supported | carried | n reproducible among ranked de-novo surface leads | 0 | — | — | 6 de-novo leads clean but 0 cross-donor reproducible → exploratory only |
| H14 | de_novo | ✅ supported | carried | IL4R presence among CD8 effector hits (Schmidt panel) | 0 | — | — | IL4R top clean_score but absent from CD8 hits + CRS-driven → score≠mechanism |
| H15 | safety_resistance | ✅ supported | new | n replicated off-tumor cell types (B8, both cohorts): CXCL16, LY6E | 6 | — | — | CXCL16 replicated myeloid flag (~55% macrophages; DUAL w/ tumor-epi) + LY6E 6 off-tumor types incl ~50% T cells → mis-targeting confirmed; CEACAM6/5 clean |

_BH q shown only for the k=6 discovery family (H02,H03,H05,H06,H09,H11), all carried from v1. 'raw'/'—' = null-confirmation or threshold/categorical test (H04,H07,H08,H10,H12,H13,H14,H15)._

## What changed vs v1
- **H07 (revised):** verdict unchanged (**supported** — gating necessary), but the leads updated to CEACAM6/CEACAM5; v1 co-lead LY6E is demoted (TME off-tumor expression) and its off-tumor risk is now tested explicitly by H15. Threshold recomputed on the v2 panel (median vital_max 0.395); call robust to threshold.
- **H15 (new):** the demoted leads' mis-targeting prediction, decidable now on B8 per-cell-type data. **Supported** — CXCL16 replicated myeloid_macrophage flag (dual with tumor epithelium), LY6E 6 replicated off-tumor cell types including T cells; CEACAM6/CEACAM5 clean (contrast control).
- **Everything else carried verbatim** — the costim/A-track results (incl. H02 CD27-Treg, H05 CD30 help-erosion) did not change.

## Provenance
- Re-run/new results: `S3_hypothesis_tests_v2/H07_v2_result.csv`, `H15_result.csv` (+ headers). Carried results: `S3_hypothesis_tests/H{02..14}_result.csv` (v1, unchanged).
- Master v2 table: `S3_hypothesis_tests_v2/S3_verdict_master_v2.csv` (change_tag column).
- Scripts: `S3_hypothesis_tests_v2/S3_build_rerun_v2.py`, `S3_render_v2.py`. Env `tools/sc-analysis-venv`. Seed=0.