# RESUME STATE — 2026-07-13 evening (account-switch safe)

## WHERE THINGS ARE
- Project root: /media/balthasar-lab/RAID4/costim_engager_counterscreen
- Xeon (ssh:balthasar-lab), 80 cores. CUDA_VISIBLE_DEVICES="". sudo=Balthasar99.
- totalvi-venv: tools/totalvi-venv/bin/python

## DELIVERABLE 1: SPATIAL_VIEWER_BUNDLE/ (569MB+, project root) — FOR THE INTERACTIVE DASHBOARD
The user's dashboard vision: model structure diagram -> click compartment -> spatial tissue/receptor
render on TIFF -> filter by gene-expr / cell-type (granular|general) -> select validation molecule ->
PK/PD overlay + that tissue's drug-distribution & cell-kill maps -> (later) simulated variants.
Contents (ALL FULL COUNT, no subsample):
- per_cell_csv/         21 tissues, 3,399,519 cells (cell_id,tissue,x_um,y_um,cell_type,<receptors>)
- tissue_images/        9 organ Xenium-morphology PNGs (base layers)
- model_matrix/         receptor_by_celltype_by_tissue.csv (2576 rows) = linked heatmap
- validation_data/pk_pd_json/     6 molecules: Cplasma_ugml, t_pk, kill_final/traj, il6, cyto
- validation_data/spatial_fields/ 6 mol x 11 organs: x,y,C_nM,bound_nM,kill,is_target,labs (.npz float32)
- validation_data/validation_molecule_list.json
- pk_datasets/          mab_tce_pkpd.sqlite (76 drugs/183 curves/936 pts) + CSV dumps + all specs
- model_structure_topology.json, receptors_manifest.json, cell_type_palette.json,
  IMAGE_COORDINATE_TRANSFORM.json, per_cell_manifest.json, README.md
DASHBOARD-READY MOLECULES NOW (real kill+distribution, NOT waiting on big run):
  mosunetuzumab, teclistamab, glofitamab, talquetamab, elranatamab, epcoritamab
  (each: PK C(t) + PD + 11-organ per-cell drug conc/bound/kill at day1 & day7)
Solid tumors (lung/ovary/prostate/skin) = raw expression NOT copies (conversion not run) — flagged.

## DELIVERABLE 2: reproducibility/finalization package (PBPK_QSP_LANE_PROMPT.md task)
Frozen: SUBMISSION_MANIFEST.json v3 = 70c9f973-dd3c-41c0-a302-5ffc4bf8e207
Done: requirements.txt, env_snapshot.txt, .gitignore, GITHUB_LAYOUT.md, README.md,
  QSP_PARAMETERIZATION_TABLE.csv (550f6044), QSP_VERIFICATION_LEDGER.csv v2 (08da5c74/2a1df74c, 55 rows),
  Fig_counterscreen_window.png (c9fa7648), FINAL_FIGURE_CATALOG.csv, DOC_01_OVERVIEW_AND_HISTORY.md.
IN FLIGHT (4 OPERON children, may still be running on resume): QSP_REPORT_SECTION.md,
  DOC_02_PK_DISTRIBUTION.md, DOC_03_BINDING_PD.md, DOC_04_IL6_CRS_AXES.md.
  Frames: 092e7650 / 4054d80d / 728bacbc / c8f1c8e1. Collect on resume; integrate DOC_01-04 + report.
NOT STARTED: qsp_reproduction.ipynb, End-to-end verify + underway note.

## IL-6 VALIDATION — DONE, verdict FAIL on verified anchors (model/IL6_VALIDATION_RESULT.json)
Emergent mech-ECF peaks: mosun 723.7, tecli 429.0, glofit 650.7, talq 772.4.
VERIFIED digitized anchors ONLY (per manifest v3): mosun 127.4 mean, glofit 30.2 median, talq 19.8 median.
  teclistamab EXCLUDED — 21 mean is a loose MajesTEC-1 value NOT a digitized curve; manifest says do not cite.
  ANY tecli-based ratio (glofit/tecli, mosun/tecli) is PROVISIONAL, not a validated result.
Rank model: tecli<glofit<mosun<talq ; clinical(verified 3): talq<glofit<mosun.
VERDICT: FAIL on verified anchors. Cleanest same-stat pair glofit/talq (both medians) WRONG DIRECTION
  (0.84x vs 1.53x); mosun/glofit dir OK but magnitude-compressed (1.11x vs 4.22x). talquetamab INVERTS
  (GPRC5D on plasma/keratin not myeloid-contacting; myeloid IL-6 term doesn't condition on target compartment).
  NOT rescued. Nomination INDEPENDENT of this layer (liability veto decides upstream).

## OTHER-LANE HANDOFF (project root): HANDOFF_DATA_MAP.md + OTHER_LANE_RUN_PROMPT.md (66 lines, sqlite backbone appended)
Massive PBPK validation runpod. 22 TCE set run-ready (20/22); 55 PK-breadth set harness GAP (build_um_any not in model/).

## NOMINATION (unaffected by IL-6): co-leads 4-1BB(TNFRSF9)+CD27; CD28 top-effector GATED. FINAL_NOMINATION_v7.md 750853c9.


## 2026-07-13 LATE — TUMOR Rtot COMPLETE + REPRODUCTION NOTEBOOK BUILT+VERIFIED
- ALL 6 solid-tumor receptor pools written on real HPA-anchored basis (nM):
  CEACAM5 6.4587(crc), ERBB2 1.9188(crc), FOLH1 3.5518(prostate), EPCAM 486.5278(ovary),
  DLL3 4.3481(lung/SCLC), PMEL 9.4729(skin/melanoma). Basis verified: kidney EPCAM reproduces 247.637.
- HPA anchors pulled by 3 sub-agents (DLL3 ENSG00000090932, EPCAM ENSG00000119888,
  PMEL ENSG00000185664 Melanocytes nCPM=3051.8 orchestrator-verified); merged into noihc_ncpm_anchor.json.
- Cell-Ontology->model-lineage crosswalk applied to lung/skin/ovary/prostate builds (celltype_ontology_to_lineage.json).
- EPCAM re-extracted into ovary ABM via exact CP10K-log1p (CEACAM5/ERBB2 reproduce corr=1.0000).
- KNOWN GAP (documented T2b §5): EPCAM/DLL3/PMEL ORGAN pools not overlaid (ref AnnData not retained ->
  scVI re-train needed). TUMOR pools complete. Heme panel + CRC pair (CEACAM5/ERBB2) unaffected.
- NEW DOC: MODEL_DOCUMENTATION_FINAL/subsystems/T2b_receptor_overlay_imputation.md (157 lines) — the
  LOAD-BEARING overlay/imputation method (scVI topk_spread K=15) + RNA->copies + retention gap. Registered in indices.
- NEW DELIVERABLE: qsp_reproduction.ipynb (21 cells) — EXECUTED end-to-end, 0 errors, all assertions pass,
  figure displays. Sections: env/inventory, receptor overlay method, Rtot, GRN nomination, QSP window,
  PK/PD, IL-6 gate (verdict FAIL reported as-is), figures. Builder: build_reproduction_notebook.py.
  Run: CUDA_VISIBLE_DEVICES='' QSP_ROOT=<root> jupyter execute qsp_reproduction.ipynb
