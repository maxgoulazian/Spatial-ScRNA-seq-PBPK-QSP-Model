# IL-6 clinical anchors — VERIFIED from digitized database (2026-07-13)

Source: model/params/mab_tce_pkpd.sqlite (curves + timeseries tables, digitized from primary literature with PMIDs/figure refs).
This supersedes ALL loose IL-6 numbers carried in engine comments (570, 340, 230, 191, 366.88 — all unsourced/fabricated/mixed-statistic).

## The ONE rule: score like-for-like (same summary_stat, comparable dose structure). NEVER mix mean vs median vs individual-peak.

## VERIFIED anchor curves (analyte='IL-6', digitized peak of the summary curve):

| molecule       | statistic | dose regimen             | clinical peak (pg/mL) | source (in db)                         |
|----------------|-----------|--------------------------|-----------------------|----------------------------------------|
| mosunetuzumab  | MEAN      | 1/2/60 mg C1 step-up     | 127.4  (day15, post-60mg) | Chen 2023 Fig2, n=212 (study B11)  |
| glofitamab     | MEDIAN    | 2.5/10/30 mg step-up     | 30.2   (day1)         | Djebli 2023 Table 2                    |
| talquetamab    | MEDIAN    | 0.4 mg/kg QW             | 19.8   (day5-6)       | Willemin 2024 Fig1a                    |
| talquetamab    | MEDIAN    | 0.8 mg/kg Q2W            | 7.9                   | Willemin 2024 Fig2a                    |
| blinatumomab   | MEDIAN    | 5/15/60 ug/m2/d step-up  | 640                   | Hosseini 2020                          |
| blinatumomab   | MEDIAN    | 60 ug/m2/d continuous    | 370                   | Hosseini 2020                          |

INDIVIDUAL-PATIENT-PEAK curves (talquetamab 213/3390, catumaxomab 37k-60k) are ORDER STATISTICS that scale with
cohort N — usable only as a worst-case envelope, NOT as a central validation target.

teclistamab: NO IL-6 curve in the digitized db, so teclistamab is DROPPED from the like-for-like IL-6 test.
(Loose "21 mean / 288 individual" MajesTEC-1 values appear in a prior lane's web-search notes but are NOT
independently verified here and NOT in this digitized panel — do not cite them as sourced anchors.)

## Like-for-like test pairs AVAILABLE:
- MEDIAN pair (cleanest, same statistic + same step-up structure): glofitamab 30.2 vs talquetamab 19.8 -> ratio 1.53x
- MEAN: mosunetuzumab 127.4 (only clean digitized mean with full dosing)

## Scale constant (KEEP, do not fit to these peaks):
kdeg_IL6 = 0.20/hr (PMID 31268236, Chen 2019 fitted) x ECF Vd 11.65 L (interstitium 8.55 + plasma 3.10)
-> CL_IL6 = 55.9 L/day. IL-6 is 21 kDa, made in interstitium, no FcRn -> ECF is the correct distribution volume
(committed on physiology BEFORE seeing results; plasma-only 3.10 L reported as labelled sensitivity only).
Model's derived-CL cross-check (~76 L/day from peaks) is an OUTPUT, not an input — reporting it is a result, using it is circular.
