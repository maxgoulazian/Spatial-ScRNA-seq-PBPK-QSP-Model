# Final Deliverables — Whole-Body Per-Cell Spatial PBPK–PD/QSP Model for Costim T-Cell Engagers
Generated 2026-07-13. All figures produced read-only from the live model results
(`model/rundir/handoff/tce_pd_*.json` and `tce_spatial_*.pkl`).

## Contents

### `FINAL_MODEL/`  — the runnable model
Self-contained copy of the 13 live engine modules + all input data. Runs from inside the folder:
```bash
cd FINAL_MODEL && python run.py teclistamab      # -> handoff/tce_pd_teclistamab.json
```
See `FINAL_MODEL/OPERATION.md` for the full operating manual.

### `figures/`
| Subfolder | Count | What |
|---|---|---|
| `spatial_resolved/` | 273 | Per-cell 5-panel overlays: cell types \| interstitial drug(nM) \| bound antibody(nM) \| receptor(copies) \| cell kill(%). One per molecule × organ (11 heme organs + tumor for 9 solid engagers). |
| `per_organ_kill/` | 20 | Per-molecule 3-panel: per-organ kill trajectory \| final kill bar (target-cell-weighted) \| per-organ IL-6. |
| `PK_overlays/` | 11 | Model plasma curve vs clinical PK (elranatamab = full digitized curve; others = label-extracted anchor points, marked in legend). |
| `PD_overlays/` | 5 | Model vs clinical: IL-6 (mosun/talq/tebentafusp) + B-cell depletion (epcoritamab/blinatumomab). |
| `celltype_overlays/` | 2 | Dense per-cell cell-type maps (bone marrow, adipose). |

### `report_and_docs/`
- `QSP_RESEARCH_REPORT.md` — publication-grade writeup (7 sections: architecture, receptor layer, binding engine, synapse/window/nomination, CRS/IL-6 gate, validation/reproducibility).
- `qsp_reproduction.ipynb` — reproduces every analysis number from committed artifacts (executes end-to-end, 0 errors).
- `MODEL_DOCUMENTATION/` — full method: subsystem docs T1–T9 + T2b, every equation at file:line.

### `figure_generators/`
Reusable scripts (read-only on the live handoff dir) — rerun against the latest sweep results anytime:
- `gen_spatial_resolved.py` + `run_all_spatial.py` + `render_one_mol.py` — spatial 5-panel overlays
- `gen_perorgan_figs.py` — per-organ kill/cytokine
- `gen_pk_all.py`, `gen_pkpd_overlays_all.py` — PK + PD clinical overlays

## Key results (honest state)
- **Nomination:** 4-1BB (TNFRSF9) + CD27 co-leads (CLEAN); CD28 gated out on CRS + suppression + proliferation.
- **Depletion** tracks target biology: heme BCMA/CD20 engagers 0.82–0.94; solid CEACAM5/DLL3/HER2 lower (ECM-throttled, on-target/off-tumor characterized).
- **IL-6 CRS axis:** mechanistic per-cell myeloid model; FAILS a validation gate (talquetamab inverts — GPRC5D on plasma/keratin, not myeloid-contacting compartments). Reported straight; does not affect the GRN/DE nomination.
- **Costim receptor density** is static resting (4-1BB induction fold unsourced → engine refuses to guess); reported windows are conservative floors for inducible arms.

## Full spatial set
273 figures also bundled at `../spatial_resolved_all_2026-07-13.tar.gz` (777 MB).
