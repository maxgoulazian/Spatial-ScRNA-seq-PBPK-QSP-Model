# Model off-tumor / therapeutic-index sweep — solid-tumor antigens

*Analysis of the PBPK/QSP model's whole-body antigen×arm sweep on three solid-tumor antigens
(EGFR, CEACAM5, CEACAM6) × 12 costim arms + no-costim baseline. Source table:
`COSTIM_OFFTOX_ANALYSIS.csv` (artifact ceb9250c), read read-only from the model lane.
All numbers below re-derived directly from that table.*

## What was run
A CD3×TAA engager was simulated on each solid antigen, with each costim arm added as a third
specificity, scoring **on-tumor kill** vs **off-tumor (normal-tissue) kill** across the whole-body
compartment set. The therapeutic index **TI = on-tumor / off-tumor kill**. This is the model's
independent whole-body read on the *antigen× arm* window — the complement to the CD3×CD19 heme
screen (`model_costim_arm_pd_screen.csv`).

## Finding 1 — efficacy is arm-invariant (the counter-screen premise, whole-body)
On-tumor kill varies only **1.3–2.0%** across all 12 arms for every antigen. Adding a costim arm does
**not** change how much tumor is killed — efficacy is CD3-gated. Arm choice is therefore a pure
**toxicity-minimization** problem, which is the entire rationale for using a CD4 toxicity screen to
pick the arm. The sweep demonstrates this directly.

## Finding 2 — the costim arm's dominant effect is to widen the off-tumor gap, and CD28 is the worst offender
Increase in cumulative off-tumor kill vs no-costim baseline:

| antigen | CD28 | 4-1BB | CD27 |
|---|---|---|---|
| EGFR | **+56.5%** | +8.2% | +6.8% |
| CEACAM6 | **+48.5%** | +7.2% | +5.5% |
| CEACAM5 | **+39.8%** | +6.1% | +4.4% |

CD28 amplifies off-tumor killing by **40–57%** with no on-tumor benefit; the CLEAN co-leads (4-1BB,
CD27) stay near baseline (**4–8%**). This is a second, whole-body-geometry axis of CD28 elimination,
independent of the GRN/DE effector-liability axis (where CD28 is z≈12 but GATED[CRS,SUPP,PROLIF]).
Two different methods, same verdict.

## Finding 3 — the model discriminates antigens on the window (antigen×arm co-selection)
Therapeutic index (no-costim baseline): **EGFR 1.29** (on-tumor > off-tumor, favorable) >
**CEACAM6 0.59** > **CEACAM5 0.54** (both < 1: off-tumor ≥ on-tumor). The **antigen sets the window**;
the arm modulates within it. CD28 is the only arm that materially degrades it — it pushes EGFR from
1.29 to **0.89, across the TI=1 line** — while co-leads preserve it (EGFR ~1.22–1.24). The framework
is thus not only an arm-selector but an **antigen×arm co-selection tool**: it ranks which TAA
tolerates which costim arm.

## Interpretation and honest scope
- This is a first-principles result from **expression + geometry with no fitting**. That CEACAM5/6
  come out TI < 1 is **consistent with** the on-target/off-tumor gut toxicity that has burdened real
  clinical CEACAM CRC programs, and that CD28 is the worst off-tumor amplifier is **consistent with
  the direction** of the TGN1412 CD28-superagonist experience. These are directional/rank
  corroborations of known liabilities — the model was **not fit** to TGN1412 or to any CEACAM
  toxicity readout, and no quantitative clinical toxicity value is being matched.
- The window here is driven by **static resting expression geometry**; it does not encode the
  CD4-subset suppressive/Treg wiring (that remains the screen's job). The two axes agree on CD28 and
  on the co-leads, which is the point.
- Scope: this table covers the three **solid** antigens. Heme (CD19) window is characterized in the
  separate CD3×CD19 screen.

## Bottom line
The sweep **reinforces** the nomination on an independent whole-body axis: efficacy is arm-invariant
(so arm choice = tox minimization), CD28 is re-eliminated (worst off-tumor amplifier), and 4-1BB/CD27
preserve the window. It also adds a genuinely new capability — **antigen×arm window ranking**.

## Files
- `COSTIM_OFFTOX_ANALYSIS.csv` — 3 antigens × 12 arms (on/off kill, TI, depletion)
- `../06_figures_and_media/figures/model_offtumor_TI_screen.png` — CD28 off-tumor amplification + TI ladder
