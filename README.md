# Costim Counter-Screen — FINAL SUBMISSION

**This is the single authoritative directory. Everything final lives here.**

**A toxicity counter-screen for costimulatory T-cell engagers: nominating a costim arm that amplifies CD8
effector function without feeding the CD4 suppressive and cytokine-release programs.**

Final nomination: liability-CLEAN co-leads **4-1BB (TNFRSF9) + CD27**.
This directory is a clean, organized view of all final deliverables (copies — originals remain in place elsewhere
in the project directory; nothing was deleted). The PBPK/QSP model is left in its own location — see 08_model_pointer/.

## Start here
- **01_report/COSTIM_COUNTERSCREEN_RESEARCH_REPORT.md** — the full research report (8 sections, including the mechanistic QSP translation)

## Layout
```
01_report/               Master research report, 7 section sources, nomination/audit/framing docs,
                         QSP report section + verification ledger
02_analysis_results/     Canonical result tables, organized by axis:
   nomination/           COSTIM_FINAL_3AXIS_SCORE_v7.csv (the 3-axis score), finalists, QSP input matrix
   toxicity_axes/        A25 CRS master-regulator scan
   effector_expression/  A29 network selectivity, A30 nomination recheck, A34 CD8-vs-CD4, A21 breadth/prior
   discovery_modules/    A22 redirector co-expression, A24 kinetic onset
   taa_selection/        TAA_finalists_6.csv
   verification/         S3 verdict master, C3 verification, V7 verification ledger
   model_output_validation/  20-molecule full-model output analysis (2026-07-14): cross-molecule
                         depletion (FIG_depletion_cross_molecule), spatial penetration->kill metrics,
                         IL-6 provenance correction. See MODEL_OUTPUT_ANALYSIS_2026-07-14.md.
   model_pkpd/           Model PK/PD cross-check: per-arm costim screen (signal-2 kill rescue + storm cytokines) + analysis memo
03_reproduction/         Reproduction notebooks: data-analysis (36/36 checks passed, *_executed.ipynb carries stored outputs),
                         binder-design (23/23 checks passed, *_executed.ipynb carries stored outputs), QSP repro notebook, requirements.txt
04_binder_design/        AF3 campaign results (funnel, scored, panels), 39 finalist structures (CIF), sequences
05_grn/                  ALL GRN data: hero full edges, QSP operators, GOF effector nets, Tier-2 A26, source, tooling
06_figures_and_media/    Key figures (GRN rings, discovery, binder, effector, model PD) + the 4-arm diffusion animation
07_datasets_reference/   DATASETS.md manifest, full data bundle zip, canonical artifact map
08_model_pointer/        Pointer to the PBPK/QSP model (NOT moved — self-contained, would break if relocated)
09_model_and_qsp/        PBPK/QSP model lane's complete bundle: runnable model (FINAL_MODEL/), QSP research
                         report + reproduction notebook, model documentation, 311 figures, generators (869M)
10_full_binder_campaign/ COMPLETE RFantibody→ProteinMPNN→AF3 campaign I/O + model weights (4.8G)
11_animation_source/     COMPLETE 4-arm diffusion timelapse source: frames, structures, trajectories, videos (1.2G)
archive/                 Superseded root notes, scratch scripts, old deliverables-layout tree (copies; nothing lost)
```

## The nomination at a glance
Of 11 costimulatory arms, **4-1BB** and **CD27** are the only two that clear the six-axis liability veto.
CD28 (top effector, z=12.11) is gated on CRS+SUPP+PROLIF. The GRN/QSP layers additionally carry CD2 as a 12th
network source (pan-lineage, CRS-gated) — see report §4.

## Model & figures (09 + 06)
- **09_model_and_qsp/FINAL_MODEL/** — the runnable whole-body single-cell PBPK/QSP model:
  `python run.py <engager>`, 13 engine modules, 31 handoff inputs, 18 agent grids, OPERATION.md, 54 model-documentation files.
- **Simulation-output figures (all from full-model runs):** 262 spatial 5-panel overlays
  (cell-type / drug / bound-receptor / receptor / kill), 20 per-organ kill, 11 PK-vs-clinical,
  2 cell-type overlays (PD/IL-6 overlays being finalized separately). Plus 13 TIFF-underlay spatial, 66 antigen maps, GRN rings, 2 diffusion animations (06).
- **Model-output validation (2026-07-14):** the 20-molecule sweep shows the heme-vs-solid depletion split
  (0.77 vs 0.18) is *emergent from geometry* under one shared kill law; spatial kill is 10-26x enriched in
  high-drug regions (drug-penetration-limited); IL-6 magnitude is read only from the validated `fin_` set.
  Full detail: **02_analysis_results/model_output_validation/**.

## Reproducibility
Every analysis number is re-derived by 03_reproduction/costim_counterscreen_reproduction.ipynb (36/36 inline
checks pass end-to-end from deposited artifacts, including the GRN runs). The binder campaign is reproduced by 03_reproduction/binder_design_reproduction.ipynb; the executed copy
(binder_design_reproduction_executed.ipynb) carries the stored '*** ALL CHECKS PASSED *** (23/23)' output, 0 errors.

## What is NOT here (intentionally)
- The live/authoring model tree (32 GB) — see 08_model_pointer/ (left in place; would break if moved).
  A self-contained RUNNABLE copy of the model + its final outputs is included under 09_model_and_qsp/.
- Raw datasets (1.7 TB in ../data/) — manifest + bundle in 07_datasets_reference/

Everything else — reports, analysis, reproduction notebooks, the FULL binder campaign, and the FULL
animation source — is inside this directory. archive/ holds superseded standalone layouts (GRN_DATA_COMPLETE,
FINAL_DELIVERABLES) whose contents are already live under 05_grn/ and the curated sections; nothing was deleted.
