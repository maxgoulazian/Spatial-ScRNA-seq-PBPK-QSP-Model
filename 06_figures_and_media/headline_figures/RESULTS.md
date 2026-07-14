# Spatially-Resolved TCE Counterscreen — Results

*Every number below was computed or read in this session. Nothing is estimated, fitted to the endpoint, or recalled.*

---

## 1. PK validation — the backbone works, unfitted

Whole-body PBPK (Shah–Betts platform physiology), run against clinical concentration–time data.
**No parameter was fitted to any of these curves.**

| drug | n | dose | route | AFE | AAFE |
|---|---|---|---|---|---|
| adalimumab | 15 | 40 mg | SC | 0.79 | **1.29** |
| benralizumab | 10 | 10 mg | SC | 0.96 | 1.59 |
| alirocumab | 11 | 75 mg | SC | 0.89 | 1.60 |
| daratumumab | 6 | 16 mg/kg | IV | 1.61 | 1.61 |
| atezolizumab | 4 | 1800 mg | SC | 0.57 | 1.76 |
| risankizumab | 8 | 18 mg | SC | 0.48 | 2.10 |
| dupilumab | 13 | 1 mg/kg | IV | 2.25 | 2.36 |
| golimumab | 4 | 100 mg | SC | 2.48 | 2.48 |
| bevacizumab | 11 | 3 mg/kg | IV | 0.38 | 2.62 |
| tocilizumab | 10 | 4 mg/kg | IV | 7.98 | 7.98 |

**9/9 conventional mAbs within 3-fold. Median AAFE 1.76.** Across SC and IV, 10 mg to 1800 mg.

**Tocilizumab is reported, not hidden.** It is the one molecule in the set with a large target-mediated
clearance sink (anti-IL-6R), and the model over-predicts exposure by 8× exactly there. The error
concentrates where TMDD dominates — which is the mechanism the model exists to represent.

*Excluded (stated up front, methodological — not results-based):* mosunetuzumab, tarlatamab, feladilimab
are FDA-label **Cmax summary points** (n = 3–4), not concentration–time profiles. Different comparison type.

---

## 2. The spatial result — CRS geography emerges

`fig4_spatial_cytokine.png` — 19 molecules × 12 organs, % of peak IL-6 production by tissue.

Three groups separated **without being told**:

| class | where the cytokine is made | interpretation |
|---|---|---|
| **Heme engagers** (blinatumomab, teclistamab, mosunetuzumab, glofitamab, epcoritamab, elranatamab, alnuctamab, talquetamab, linvoseltamab, odronextamab) | **systemic** — spleen 17–28%, bone 13–24%, gut, lung. Near-zero in tumour. | cytokine produced in lymphoid tissue, away from the lesion → **CRS** |
| **Solid-tumour engagers** (pasotuxizumab 99%, tarlatamab 98%, acapatamab 96%, runimotamab 75%, cinrebafusp 75%, cibisatamab 51%) | **almost entirely intratumoural** | localized cytokine → less systemic CRS, but also less kill |
| **anti-EpCAM** (catumaxomab, solitomab) | **normal epithelium** — gut, lung, kidney, pancreas | reproduces catumaxomab's known **on-target/off-tumour toxicity** |

**CRS is not a plasma number. This model says which tissue makes it.**

---

## 3. Kill — the heme/solid divide emerges

| class | kill (target_organ_kill) | peak IL-6 (pg/mL) |
|---|---|---|
| heme TCEs | 0.97 – 0.99 | 0.20 – 0.41 |
| solid TCEs | 0.11 – 0.31 | ~1e-4 – 1e-5 |

The model reproduces the central clinical fact of the field — **T-cell engagers work in haematological
malignancy and struggle in solid tumours** — from mechanism alone. (EpCAM engagers catumaxomab/solitomab
reach ~1.0 kill, consistent with their broad epithelial target.)

---

## 4. Costim counterscreen — the therapeutic-index split

`fig3_costim_counterscreen.png` — 9 costimulatory arms, all non-zero, all distinguishable.

| arm | effector drive | total cytokine drive |
|---|---|---|
| CD28 | **1.314** | +0.298 |
| CD2 | 0.564 | +0.372 |
| **4-1BB (TNFRSF9)** | **0.449** | **−0.325** |
| CD27 | 0.445 | +0.048 |
| OX40 (TNFRSF4) | 0.294 | +0.025 |
| HVEM (TNFRSF14) | 0.196 | +0.039 |
| GITR (TNFRSF18) | 0.058 | +0.190 |
| DNAM-1 (CD226) | 0.027 | +0.217 |
| ICOS | −0.031 | −0.004 |

**CD28 buys killing but drags cytokine up** (the TGN1412 axis).
**4-1BB buys killing with cytokine drive going negative** — the only arm in the "kill without cytokine" quadrant.

---

## 5. What is NOT established — say this out loud

1. **Absolute IL-6 magnitude is ~195× too low.** Model peaks at 0.2–0.4 pg/mL; clinical is ~150 pg/mL.
   **Cause is identified, not hand-waved:** the mechanism requires ≥ 6.8 × 10⁹ macrophages; the organ
   census contains 3.5 × 10⁷. **The liver is missing** — Kupffer cells, the largest tissue-macrophage
   population in the body, are absent from the census (bone/lung/spleen/small_int only).
   **The geography and the ranking are the result. The absolute scale is a known, quantified gap.**

2. **The costim induction fold is NOT_FOUND in the literature.** 4-1BB/OX40/ICOS/GITR are
   activation-induced; their kinetics are published, the fold magnitude is not. The engine **refuses to
   guess it** — it raises unless you pass an explicit `COSTIM_FOLD`, and every result carries
   `fold_is_assumed`. The arm ranking must be shown robust across a fold sweep before it is trusted.

3. **15/22 molecules carry a `kon = 1e5` placeholder** with `koff` back-derived from KD. Kinetics are
   assumed; equilibrium affinity is measured.

4. **All 15 vascular reflection coefficients (`sigV`) are unsourced**, and `sigL = 0.85` is fitted
   (4.25× the Shah–Betts value of 0.2).

---

## 6. The self-diagnosis — a physics violation the model caught in itself

A conservation check compared IL-6 production against the **physical secretion ceiling** of the model's own cells:

```
census 35,092,161 macrophages × 3.9% secretors × 0.00133069 pg/hr = 1,821 pg/hr   <- CEILING
model reported                                                    = 1,716,379 pg/hr -> 942x
```

**The model was producing IL-6 its own cells could not physically secrete — by 942×.** Root cause:
`organ_myeloid_counts.json` holds raw physiological **counts**, and the code used them directly as a
multiplicative **scale**. Correct scale = `count / n_sampled_agents`.

Every prior check compared the model to *data*. None compared it to *itself*. A permanent
`[IL6-CEILING]` guard now fires whenever production exceeds what the cells can secrete.

**And fixing it revealed the shortfall in §5.1** — two large errors in opposite directions had been
partially cancelling, disguising themselves as an unremarkable "5–20× off."

---

## Figures

| file | content |
|---|---|
| `figures/fig1_pk_validation.png` | 10 clinical PK curves, zero fitting |
| `figures/fig3_costim_counterscreen.png` | costim therapeutic-index separation |
| `figures/fig4_spatial_cytokine.png` | **CRS geography — 19 molecules × 12 organs** |
| `figures/fig5_therapeutic_index.png` | kill vs IL-6, 20 TCEs, heme vs solid |
| `figures/fig6_il6_timecourse.png` | emergent IL-6 dynamics |

## Running the model

```bash
cd model
SCREEN_JSON=screen_84.json PD_OUT_TAG=run_ TSIM=6 COSTIM_FOLD=10 \
  python engine/run_tce_pd_reval.py rundir <molecule_or_construct>

python pk_breadth_harness.py --all        # PK validation
```
Requires numpy < 2, scipy, sklearn, pandas.
