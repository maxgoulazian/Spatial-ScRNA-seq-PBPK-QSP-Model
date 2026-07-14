# Model PK/PD/spatial analysis — costim-arm screen

*Analysis of the PBPK/QSP model's own per-arm output (85 `tce_pd_*.json` runs in the live model
rundir), read read-only. This is the mechanistic model's independent read on the nomination.
Source: `09_model_and_qsp/FINAL_MODEL` engine; per-arm results table `model_costim_arm_pd_screen.csv`.*

## What was run
The model executed a costimulatory-arm screen on a **CD3×CD19** backbone (a B-cell-malignancy
surrogate where the model has validated clinical PK anchors): a **bispecific** baseline (signal 1
only) and, for each costim receptor, a **trispecific** CD3×CD19×`<costim>` construct adding signal 2.
Each run emits plasma PK, the four storm cytokines (IL-6/IFN/TNF/IL-2, systemic cumulative), per-organ
and target-organ kill, and heme depletion.

## Finding 1 — the model confirms the signal-2 rationale
Adding **any** costim arm raises target-organ kill from **0.086 → ~0.975** (≈11×) and heme depletion
from 0.68 → 0.999. This is the core mechanistic premise of the whole project — a CD3-only engager
under-kills without costimulation — reproduced by an independent PBPK/QSP model that never saw the
Perturb-seq screen.

## Finding 2 — at matched dose the arms are near-identical *in this model*, and that is expected
Across the 7 costim arms the model resolves, target-organ kill spans only **0.9749–0.9758** and total
storm-cytokine burden spans only **~4%** (4907–5106 a.u.). The arms are not separable on the model's
systemic-cytokine output. This is not a null result against the nomination — it is the expected
consequence of what the model does and does not encode:

- The model uses **static resting receptor density** and a systemic myeloid IL-6 term. It does **not**
  encode the CD4-subset-specific suppressive/Treg wiring or the per-receptor differential regulation
  of the storm cytokines — which is precisely the axis the CD4 Perturb-seq counter-screen resolves and
  the model cannot.
- Therefore the arm-discriminating signal lives **upstream**, in the screen-derived liability veto
  (4-1BB and CD27 the only CLEAN arms), not in this model's systemic PD. The model corroborates
  *efficacy* and *PK*; the screen supplies *arm selection*. The two are complementary, not redundant.

## Finding 3 — PK is validated; the IL-6/CRS axis is a characterized honest negative
- **PK validated independently:** route-matched digitized-curve overlays give **teclistamab SC AFE
  1.29×** and **elranatamab (matched PK+IL-6) AFE 2.06×**; the BiTE renal-clearance term predicts a
  fall to 50% of C₀ in **1.78 h vs the BLINCYTO label's 2.11 h**.
- **IL-6 CRS axis FAILS its validation gate** (`IL6_VALIDATION_RESULT.json`): among the three
  manifest-verified digitized anchors, the cleanest same-statistic pair (glofitamab/talquetamab
  medians) is **wrong direction** (model 0.84× vs clinical 1.53×), driven by **talquetamab inversion**
  — GPRC5D sits on plasma cells / keratinized tissue, not myeloid-contacting compartments, and the
  myeloid IL-6 term does not condition on which compartment killing occurs in. This is reported
  straight, not buried.

## Bottom line for the nomination
The model **strengthens** the case on the two axes it can speak to (signal-2 efficacy; validated PK)
and is **appropriately silent** on arm discrimination — the differences among CLEAN and gated arms are
resolved by the counter-screen, which the model structurally cannot see. Because the nomination is
decided upstream of the QSP IL-6 output, the CRS-axis validation failure **does not affect** the
4-1BB + CD27 nomination. Treat the model's per-arm CRS numbers as **non-discriminating** and cite the
screen-side liability axes for arm selection.

## Files
- `model_costim_arm_pd_screen.csv` — per-arm PD scalars (kill, depletion, storm cytokines, Cmax)
- `../06_figures_and_media/figures/model_pd_costim_screen.png` — signal-2 rescue + per-arm efficacy-vs-CRS
