# S3 — FINAL SYNTHESIS v3 (widen-the-slate re-frame; the C1-v3 ingestion target)
<!-- module S3 (v3) | S3v3_write_synthesis.py | generated 2026-07-10T14:41:12Z | env tools/sc-analysis-venv | RE-FRAME of S3_SYNTHESIS_FINAL_v2.md: carries the widened top-15 surface-QC'd TAA slate x the 3 effector-hit non-T4 costim slate. NO scRNA re-analysis, NO hypothesis re-test. Every value READ from a verified module CSV (B6v2/B5/B8/B3v2/B7/A8). Carried v2 hypothesis verdicts embedded BYTE-IDENTICAL (src S3_HYPOTHESIS_RESULTS_v2.md sha256=c76a4e6477ef601e15544c6bb1132bc62d2d06bfd819bdda11176c8136c08d48). Reads N_TAA/N_COSTIM as parameters. NAMES NO WINNER. -->

**What this is.** The v2 synthesis, RE-FRAMED to carry the full analysis rather than a few named leads. The v2 file distilled the TAA arm to a short lead list; the genome-wide analysis scored all 33,694 genes and the ranking holds many more surface-valid candidates. This v3 re-frame carries the **top-15 surface-valid TAA finalists** (mechanically = the first 15 rows by `rank_safety_v2` that pass surface-QC) x the **3 effector-hit non-T4 costim arms** (tier-ordered), each with its honest flags, straight through to the QSP hand-off. It re-frames the synthesis FRAMING only; it re-runs no analysis and re-tests no hypothesis. **This synthesis names no winner** — the slate emerges from rank + flags; QSP is the discriminator and Track D makes the final pick.

---

## TAA SLATE — top-15 surface-valid finalists (rank_safety_v2 order)

The candidate ordering is the verified v2 genome-wide per-cell-type re-rank (`B6_taa_ranked_genomewide_v2.csv`, column `rank_safety_v2`, rank 1 = best) — **not re-derived here**. Surface-QC then removes predicted-surface false-positives whose UniProt topology is internal-compartment / cytoplasmic-catalytic / secreted-only (see the QC audit + method note below). The first 15 PASS rows, in rank order, are the slate. Rank + flags speak; no row is editorialized as a 'lead' or a 'drop'.

- **Slate size:** 15 carried (target N_TAA=15; shortfall=0). 26 of 31 ranked rows pass surface-QC.
- **Surface-QC anchors (must exclude where in-window):** EBP=EXCLUDE, PIGT=EXCLUDE, RNF149=EXCLUDE — EBP (Q15125, ER/nuclear membrane, rank 4) and PIGT (Q969N2, ER membrane, rank 11) are excluded from the slate they would otherwise have entered; RNF149 (Q8NC42, cytoplasmic-facing RING E3 ligase, rank 29) also excluded.

| # | gene | rank_safety_v2 | rank_balanced_v2 | z_restriction (B8) | off_tumor_flag (B8) | healthy vital_max (B3v2) | top_location (B3v2) | density copies/cell (B7, QSP-input) | label_source (B5) | surface-QC |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | **CEACAM6** | 1 | 1 | 2.3725 | none | 0.7563 | Skin:Epithelium | TBD | GPI (UniProt) | PASS |
| 2 | **CEACAM5** | 2 | 2 | 2.3929 | none | 0.0992 | Large_Intestine:Epithelium | 65000-257000 | GPI (UniProt) | PASS |
| 3 | **SLC52A2** | 3 | 5 | 0.3724 | fibroblast_stromal | 0.1452 | Skin:Endothelium | TBD | pos. trainingset | PASS |
| 4 | **ITGB4** | 5 | 3 | 1.3943 | endothelial | 0.2074 | Heart:Neural | TBD | pos. trainingset | PASS |
| 5 | **TSPAN6** | 6 | 11 | 1.0549 | none | 0.1586 | Ovary:Epithelium | TBD | machine learning | PASS |
| 6 | **CXCL16** | 7 | 7 | 0.222 | myeloid_macrophage | 0.6054 | Testis:Germline | TBD | pos. trainingset | PASS |
| 7 | **SLC2A1** | 8 | 4 | 1.4883 | none | 0.0688 | Bladder:Epithelium | TBD | pos. trainingset | PASS |
| 8 | **SDC1** | 9 | 9 | 1.0232 | fibroblast_stromal | 0.5878 | Bladder:Epithelium | TBD | machine learning | PASS |
| 9 | **CD320** | 10 | 13 | 0.2242 | endothelial | 0.3271 | Pancreas:Endothelium | TBD | machine learning | PASS |
| 10 | **DPEP1** | 12 | 12 | 1.9735 | none | 0.3055 | Small_Intestine:Epithelium | TBD | GPI (UniProt) | PASS |
| 11 | **CD24** | 13 | 10 | 1.5383 | B_plasma | 1.387 | Prostate:Epithelium | TBD | GPI (UniProt) | PASS |
| 12 | **LY6E** | 14 | 8 | -0.1991 | B_plasma;T_cell;endothelial;fibroblast_stromal;mast;myeloid_macrophage | 1.2012 | Ovary:Epithelium | TBD | GPI (UniProt) | PASS |
| 13 | **CD46** | 15 | 14 | 0.3642 | endothelial;fibroblast_stromal | 1.0242 | Skin:Epithelium | TBD | machine learning | PASS |
| 14 | **SERINC3** | 16 | 19 | -0.129 | endothelial;fibroblast_stromal;mast;myeloid_macrophage | 0.8573 | Pancreas:Endothelium | TBD | machine learning | PASS |
| 15 | **ATRAID** | 17 | 18 | -0.1347 | B_plasma;T_cell;endothelial;fibroblast_stromal;mast;myeloid_macrophage | 0.5756 | Stomach:Stromal | TBD | pos. trainingset | PASS |

_Evidence-field provenance (per carried TAA): `rank_safety_v2`/`rank_balanced_v2` ← B6_taa_ranked_genomewide_v2.csv; `z_restriction`(=`z_tumor_cell_restriction_combined`)/`off_tumor_flag`(=`off_tumor_celltype_flag_replicated`) ← B8_taa_restriction_scores.csv; healthy breadth (`vital_max`,`max_healthy_mean`,`n_expr_compartments`,`top_location`) ← B3v2_tabula_sapiens_healthy_breadth.csv; `density_copies_cell` ← B7_antigen_densities.csv (**QSP INPUT parameter only, never a selector**; TBD where unsourced — only 16 targets carry a curated absolute density); surface flags (`label_source`,`er_secretory_flag`,`uniprot_accession`) ← B5_surface_flags.csv + the S3v3 surface-QC audit. Full per-row table: `S3v3_taa_slate.csv`. QC decisions: `S3v3_surface_qc_audit.csv`._

**Honest flags shown, not editorialized.** Several carried rows have replicated off-tumor expression flags (e.g. LY6E and ATRAID flag across 6 non-tumor compartments incl. T cells; CXCL16 flags myeloid/macrophage). These are carried in rank order WITH their flags visible — the flags inform QSP; they are not used here to drop a row.

### Surface-QC method (agent-level UniProt lookup; script did no network I/O)

Per the orchestration header, the UniProt subcellular-location lookup is an **agent-level pre-fetch**, not an in-script REST call (the box is headless with no established outbound HTTPS; an in-script call would silently no-op and default every predicted-surface candidate to PASS — the exact failure that pulls EBP/PIGT back in). The Claude lane pre-fetched UniProtKB subcellular-location + membrane-topology for every `machine learning` candidate accession via the connected **genes-ontologies** MCP server (`get_uniprot_entries`, fields mode — the in-environment equivalent of the spec's `bc_get_uniprot_protein_info`, hitting the identical UniProtKB source) and wrote a STATIC table (`S3v3_gene_localization_table.csv`); the build script reads only that table.

**Decision rule (cell-surface topology, not lumen exposure):** a `pos. trainingset` / `GPI (UniProt)` label is experimentally grounded → PASS. A `machine learning` label is confirmed against UniProt: **KEEP** if it carries any 'Cell membrane'/'Cell surface' annotation, or an Extracellular topological domain plus a transmembrane segment (a single/multi-pass PM protein whose ectodomain is extracellular — even with a co-annotated shed/secreted isoform). **EXCLUDE** only when it has no plasma-membrane/cell-surface localization AND is sole-internal (ER/Golgi/mito/nuclear/endosomal), secreted-only with no membrane form, or cytoplasmic-catalytic. No annotation → HOLD (not a silent PASS).

## COSTIM SLATE — 3 effector-hit non-T4 arms (data_driven_tier order)

From `COSTIM_FINAL_3AXIS_SCORE.csv` (A8): arms that are effector hits (`E_hit==True`) AND not `data_driven_tier=='T4_liability'`, ranked T1_favorable_robust > T2_favorable > T3_neutral_or_insufficient. Each carries its data-driven caveat read straight from A8 (no fabricated caveat; TBD where none recorded). **No arm is pre-picked** — QSP adjudicates the trade-offs across the slate.

| tier | receptor | alias | E_schmidt_z | data-driven caveats (from A8 HELP/SUPP/tregfrac fields) |
|---|---|---|---|---|
| T1_favorable_robust | **CD27** | CD27 | 4.278 | SUPP-lowering: SUPP_call=down*conc, SUPP_agon=-0.0819, q=6.259755763277066e-05 (A8/A2); Treg-fraction (4-donor MH, neutral flag for QSP): OR=1.0139, q=0.9328, concordant=no |
| T2_favorable | **TNFRSF8** | CD30 | 3.224 | HELP-eroding: HELP_call=down, HELP_agon=-0.1163, q=2.260430284395217e-06 (A8/A2); SUPP-lowering: SUPP_call=down*conc, SUPP_agon=-0.2015, q=4.8958431071715954e-17 (A8/A2); Treg-fraction (4-donor MH, neutral flag for QSP): OR=1.808, q=0.0, concordant=no |
| T3_neutral_or_insufficient | **TNFRSF9** | 4-1BB | 3.741 | Treg-fraction (4-donor MH, neutral flag for QSP): OR=0.8344, q=0.30835, concordant=no |

_Slate size: 3 carried (target N_COSTIM=3; shortfall=0). Per-arm axes/full fields: `S3v3_costim_slate.csv`._

**T4 / effector-gated eliminations (per A8, reasons stated):**

| receptor | alias | E_hit | tier | elimination reason |
|---|---|---|---|---|
| CD28 | CD28 | True | T4_liability | T4_liability -- feeds CRS + SUPP (TGN1412-class superagonist risk) |
| ICOS | ICOS | False | T4_liability | effector-gated-out -- below E_hit threshold, cannot supply signal-2 |
| CD226 | DNAM1 | False | T4_liability | effector-gated-out -- below E_hit threshold (DNAM1) |
| TNFRSF4 | OX40 | True | T4_liability | T4_liability -- OX40: dominant private module is Treg-suppression (feeds SUPP) |
| TNFRSF18 | GITR | False | T4_liability | effector-gated-out -- below E_hit threshold (GITR) |
| TNFRSF14 | HVEM | False | T4_liability | effector-gated-out -- below E_hit threshold (HVEM) |
| TNFRSF25 | DR3 | False | T4_liability | effector-gated-out -- below E_hit threshold (DR3) |
| CD40 | CD40 | True | T4_liability | T4_liability -- APC-side receptor, not a cis-costim on the engaged T cell (mechanism-gated) |

_CD28 / OX40 (TNFRSF4) feed the CRS / SUPP liability axes (TGN1412-class); CD40 is APC-side, not a cis-costim on the engaged T cell; ICOS / DNAM1(CD226) / GITR(TNFRSF18) / HVEM(TNFRSF14) / DR3(TNFRSF25) are below the effector hit threshold and cannot supply signal-2. Full elimination table: `S3v3_costim_eliminations.csv`._

## HAND-OFF TO C1-v3 / QSP

C1-v3 ingests this section as authoritative: the **top-15 surface-valid TAA finalists** (rank-ordered, flag-annotated) x the **3 effector-hit non-T4 costim arms** (tier-ordered, elimination-reasoned). Density is a QSP **input parameter**, not a selector; TBD where unsourced. Literature caveats (soluble-antigen sink on CEA-family members, glycoform-blind epitopes, no-HPA-antibody targets, all-agonism-is-computational, tumor-conditional gating requirement) travel as QSP-input context, never as demotions. **QSP is the discriminator across the slate; Track D makes the definitive post-QSP pick. This synthesis names no winner.**

---

## APPENDIX — carried-forward verified v2 hypothesis verdicts (NOT re-tested)

These verdicts are carried **byte-identical** from the verified v2 hypothesis run — they are NOT recomputed or re-tested in v3 (the reuse rule forbids re-running the retired tests). The block below is the verbatim content of `S3_HYPOTHESIS_RESULTS_v2.md` (sha256=c76a4e6477ef601e15544c6bb1132bc62d2d06bfd819bdda11176c8136c08d48), reproduced between sentinels for the C3 byte-diff check.

<!-- BEGIN CARRIED_V2_VERDICTS (byte-identical to S3_HYPOTHESIS_RESULTS_v2.md) -->
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
<!-- END CARRIED_V2_VERDICTS -->

_Verdict master (unchanged): S3_hypothesis_tests_v2/S3_verdict_master_v2.csv (12 supported / 2 inconclusive / 0 refuted; BH family k=6 all q<0.05). Conclusions: S1_CONCLUSIONS_v2.md. Evidence matrix: S1_taa_evidence_matrix_v2.csv._
