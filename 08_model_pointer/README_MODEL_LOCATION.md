# PBPK/QSP Model — Location Pointer (NOT moved)

The PBPK/QSP model is DELIBERATELY LEFT IN PLACE and NOT copied into REPO/, because it is a
live, self-contained pipeline with internal relative paths and run-staging that would break if moved.

## Where it lives (unchanged)
- `../model/`                    — the live model (engine, config, run staging) — 24 GB
- `../model_spec/`               — model specification + figure renderers (grn_ring.py etc.)
- `../model_snapshots/`          — pre-unification snapshots
- `../MODEL_DOCUMENTATION_FINAL/`— the current model documentation (MASTER_MODEL_DOCUMENTATION.md, START_HERE.md)
- `../MODEL_DOCUMENTATION/`      — earlier model doc chapters (DOC_01..05)
- `../spatial_submodel/`, `../tumor_builds/`, `../SPATIAL_VIEWER_BUNDLE/` — spatial/tumor model layers
- `../qsp_reproduction.ipynb`    — QSP reproduction notebook (also copied to REPO/03_reproduction/ as reference)

## Verified
The live model reads NOTHING outside `model/` by absolute or ../ path (checked: it loads via its own
config). Archiving/reorganizing the analysis dirs does not affect it.

## Integration contract (how the model consumes the analysis outputs)
- magnitude  <- REPO/05_grn/ drive tables (cd4_grn_qsp_drive_full.csv, per_arm_drive_magnitude_uncertainty)
- shape      <- REPO/05_grn/02_qsp_operators/grn_operator_shared_backbone.parquet on the 3 CITE-seq baselines
- timescale  <- hero Rest/8hr/48hr kinetics
- veto       <- applied UPSTREAM (nomination = CLEAN {4-1BB, CD27}); model is downstream window characterization
See REPO/01_report/PBPK_QSP_LANE_PROMPT.md for the full handoff.
