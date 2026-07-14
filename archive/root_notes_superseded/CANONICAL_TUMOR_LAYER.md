# CANONICAL TUMOR LAYER — costim engager counter-screen
**User-approved as THE tumor layer for all lanes (2026-07-08).**

The spatially-explicit tumor substrate is built on the **10x Genomics public Xenium** dataset
*FFPE Human Colorectal Cancer* (Human Immuno-Oncology Profiling Panel + custom add-on, Xenium
Onboard Analysis 2.0.0, **CC BY 4.0**). Section P1, 307,762 cells, ~6.6 × 7.2 mm, 0.2125 µm/px.
Cell types = scANVI 10-class label-transfer from the Pelka CRC atlas.

## Use these artifacts
| role | filename | version_id |
|---|---|---|
| **per-cell data package (the layer)** | xen_p1_cells_enriched.npz | `3288057b-2ca1-42f8-8087-5db1662591d0` |
| rasterized grid fields (50 µm, for PDE/QSP) | spatial_grid_fields.npz | `d999df20-bdfe-4875-94c2-05089bba5869` |
| Fig 1 — tumor architecture (broad + fine) | xenium_fig1_tumor_architecture.png | `7f4f9d3b-e5ef-444a-af66-86a7aeba11e2` |
| Fig 2 — T-cells + vessels/lymphatics | xenium_fig2_tcell_vasculature.png | `82591676-bdf5-496e-83ea-6ceaadfe4ba8` |
| captions / provenance / caveats | xenium_overlay_captions.md | `b172030b-7f6e-4726-afa3-fca6eeb79097` |

## Per-cell package schema (`xen_p1_cells_enriched.npz`, numpy, allow_pickle=True)
- `x`, `y` — float cell centroids (µm)
- `lab` — uint8 class index into `labels`
- `prob` — uint8 = scANVI max-probability × 255
- `labels` = [B_plasma, CD4conv_T, CD8_T, Endothelial, Epithelial, Fibroblast_stromal, Myeloid, Other, Treg, otherT_NK]
- `endo_subtype` — uint8: 0 = non-endothelial, 1 = blood endothelium (BEC, **drug IN**), 2 = lymphatic endothelium (LEC, **drug OUT**)

## Class counts
Epithelial 135,861 · Other 71,510 · Fibroblast_stromal 38,985 · CD4conv_T 15,191 · Myeloid 13,592 ·
Endothelial 11,202 · B_plasma 8,108 · CD8_T 7,172 · Treg 6,000 · otherT_NK 141. BEC 10,068 · LEC 1,134.

## Key spatial finding
CD8 killers concentrate in a centre-left immune zone and are largely **excluded from the tumour
glands** — this exclusion is what the spatial transport/kill model quantifies.

## Caveats (load-bearing — keep in any downstream use)
1. scANVI over-calls rare classes; **Treg** calls are lower-confidence (mean max-prob 0.855; ~70% at
   prob ≥ 0.8). For QSP proportions use the **Pelka reference split, NOT raw scANVI Treg counts**.
2. **LEC = Endothelial ∩ LYVE1⁺ only** (1,134 cells; PROX1 too noisy). The lymphatic-drain field is
   the weakest data element.

## Lineage note
`qsp_pk_multimolecule_validation.png` = artifact `50ed0a88` (latest v `6867a444`).
`qsp_spatial_overlay.png` (`6c555e10`) is the OLD 3-zone schematic, **not** the tumor layer.
