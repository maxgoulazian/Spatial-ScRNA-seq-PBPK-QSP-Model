# Costim-Arm Sweep — Solid CRC Antigens

Whole-body per-cell PBPK–PD sweep of CD3×TAA(×costim) constructs across 11 costim arms +
no-costim control on three solid CRC antigens (EGFR, CEACAM5, CEACAM6), CD3-off (`bi`) and
CD3-on (`tri`) formats. Tests the costim nomination on an on-target/off-tumor toxicity axis
the CD4 counter-screen structurally cannot see.

## Headline
1. **Efficacy CD3-gated & arm-invariant** — 9.1× solid / 32× heme gate; no-costim = armed pack.
2. **On-target/off-tumor TI** — EGFR 1.29, CEACAM6 0.59, CEACAM5 0.54 (CEACAMs kill normal
   colon harder than tumor; off-tumor lands in large intestine — clinically correct).
3. **CD28 dominant off-tumor amplifier** — +40–57% off-tumor kill, zero on-tumor benefit;
   ~5× next arm (CD2); inducible arms (4-1BB/OX40/GITR/ICOS) at baseline. Independently
   confirms the CD28 veto (GRN/DE z=12.11) on a separate whole-body axis; reproduces the
   TGN1412 signature from expression alone.

## Limitation
IL-6/CRS term is identically zero for all 36 solid constructs (compartment-coverage gap —
solid ABM does not populate the myeloid pool the IL-6 engine reads). Solid ranking rests on
kill-based off-tumor TI only.

## Files
- `FIG_costim_offtumor_window.png` — 3-panel: TI ladder · off-tumor-by-arm · per-arm Δoff by receptor class
- `FIG_EGFR_sweep.png`, `FIG_costim_efficacy_invariance.png`, `FIG_costim_sweep.png`
- `COSTIM_OFFTOX_ANALYSIS.csv` — per-construct on/off-tumor kill + TI (3 TAAs × 12 arms)
- `EGFR_SWEEP_ANALYSIS.csv`, `COSTIM_SWEEP_ANALYSIS.csv`
- `COSTIM_SWEEP_SOLID_RESULTS.md` — full results writeup
- Spatial overlays: `../../06_figures_and_media/costim_sweep_spatial_overlays/` (432 PNGs, 3 TAA × 12 arms × 12 organs)
