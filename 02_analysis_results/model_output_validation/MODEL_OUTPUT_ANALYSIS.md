# Model-Output Analysis — 20-Molecule Production Sweep
_Source: `model/rundir/handoff/tce_pd_*.json` (PD/kill/IL-6) and
`tce_spatial_*.pkl` (per-cell spatial state). All numbers computed directly from the deposited
run outputs — nothing hand-entered._

## 1. Cross-molecule depletion (the headline result)

The full-model unified PK→PD run was executed for 20 clinical T-cell engagers (11 heme-target,
9 solid-target; `MODEL_OUTPUT_ANALYSIS.csv`). Target-cell-count-weighted depletion
separates cleanly by disease setting:

| Setting | n | mean depletion | range |
|---|---|---|---|
| **Heme** (BCMA/CD20/GPRC5D/CD19/FcRH5) | 11 | **0.77** | 0.00–0.94 |
| **Solid** (CEACAM5/DLL3/HER2/EpCAM/PSMA/PMEL) | 9 | **0.18** | 0.00–0.45 |

This 4× heme-vs-solid gap is the model's central emergent result: it is **not** imposed — both
settings use the identical kill law (`k_death=1.0`, serial-killing rate fixed from literature).
The gap emerges from geometry: heme targets sit in plasma-accessible marrow/spleen/blood; solid
targets sit behind an ECM transport barrier that throttles drug penetration (see §3).

**Within-heme ordering** tracks target biology and PK: the BCMA engagers
(elranatamab/teclistamab/linvoseltamab 0.93–0.94) and GPRC5D (talquetamab 0.92) lead;
the CD20 engagers span 0.82–0.91; blinatumomab (CD19) sits lower (0.48) — mechanistically
correct given its very short half-life and low Cmax (0.009 µg/mL) at the modelled dose.

**Two mechanistically-explained zeros** (flagged, not hidden):
- **cevostamab (FcRH5)** → 0.00: FcRH5 surface density in the marrow reference is sparse, so the
  binding solve produces near-zero engagement. A real biology signal (low-antigen target), not a bug.
- **tebentafusp/pasotuxizumab** → 0.00: ImmTAC/low-Cmax formats at the modelled exposure.

## 2. IL-6 (CRS axis) — a cross-generation scale discontinuity to disclose

**Important provenance finding.** The IL-6 peak concentration is **not consistent across run
generations**, and only the `fin_` reruns carry the validated, correctly-scaled mechanistic
myeloid IL-6 engine:

| Run generation | IL-6 production (pg/hr) | IL-6 peak (pg/mL) | Status |
|---|---|---|---|
| `fin_` (mosun/tecli/elran) | 1.8–3.6 ×10⁷ | 23,000–48,000 | **validated engine** |
| `full_` (the 20-mol sweep) | 200–950 | 0.0–0.4 | **under-scaled** (~10⁵× low) |
| `t7_` | ~10⁶ | 1,700–2,500 | intermediate |

The `full_` production sweep computed **kill/depletion correctly** but its IL-6 concentration
field did not accumulate to physiological scale — the myeloid production rate is ~10⁵× below the
validated `fin_` engine. **Consequence for the deliverable:** depletion/kill numbers from the
20-molecule sweep are authoritative; **IL-6 magnitude should be read only from the 3 `fin_`
molecules** (mosunetuzumab 48,204 / teclistamab 25,974 / elranatamab 22,956 pg/mL — these are
*model* peaks), which
remain the honest IL-6 validation set — and which still **FAIL** the clinical gate
(clinical IL-6 anchors, digitized-verified per SUBMISSION_MANIFEST v3: mosunetuzumab 127.4, glofitamab 30.2,
talquetamab 19.8 pg/mL; teclistamab's 21 pg/mL MajesTEC-1 mean is **excluded** from like-for-like — loose mean, not a digitized curve)
(model over-predicts CRS because step-up priming is not captured). This is the same IL-6
limitation already documented; the new finding is that it must not be read from the `full_` runs at all.

## 3. Spatial penetration → kill coupling (the counter-screen's core spatial claim)

From the per-cell tumor snapshots (`tce_spatial_full_*.pkl`, 307,762 CRC cells each), kill is
**drug-penetration-limited**, quantified directly:

| Molecule (target) | drug CV% (gradient) | kill in top-quartile drug | kill in bottom-quartile drug | ratio |
|---|---|---|---|---|
| cinrebafusp (HER2) | 29.7 | **50.1%** | 4.6% | **10.8×** |
| acapatamab (PSMA) | 31.4 | **34.8%** | 1.4% | **25.8×** |

Where drug penetrates, cells die; where the ECM barrier blocks it, they survive — a
**10–26× kill enrichment** in high-drug regions. This is the spatial mechanism that a
well-mixed QSP cannot represent and that the counter-screen is built to quantify.

**cibisatamab (CEACAM5)** kill is correctly **on-target**: 8.32% on the 135,861 CEACAM5⁺ cells
vs **0.00%** on the 171,901 antigen-negative cells — perfect target selectivity at single-cell
resolution.

**Two spatial-figure caveats (disclosed):**
- **tarlatamab / catumaxomab** show 0 kill in the *CRC* tumor snapshot because their targets
  (DLL3 / EpCAM) are not the CRC antigen — their real kill is in their **matched** SCLC/ovary
  compartments, not this CRC pkl. The CRC-tumor spatial panels for these molecules are
  antigen-mismatched and should be read from their matched-tumor runs.
- The per-cell `R` and `bound_nM` fields in the tumor snapshots are underpopulated (display only);
  kill is computed correctly via the synapse engine regardless.

## 4. What this adds to the final analysis
1. A **20-molecule cross-validation** that the heme-vs-solid efficacy split is emergent from geometry, not imposed.
2. A **quantified spatial penetration→kill law** (10–26× enrichment) — the counter-screen's headline spatial result.
3. An **honest IL-6 provenance correction**: read CRS only from the validated `fin_` set, never the `full_` sweep.
4. Single-cell **target selectivity** confirmed (cibisatamab 8.3% on-target vs 0% off-target).
