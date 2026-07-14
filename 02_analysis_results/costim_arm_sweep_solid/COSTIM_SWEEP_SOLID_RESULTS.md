# Solid-Tumor Costim-Arm Sweep — Results

**Substrate.** A stage-1 flat screen of CD3×TAA(×costim) constructs was run through the full
per-cell PBPK–PD engine on three solid CRC-relevant antigens — **EGFR, CEACAM5, CEACAM6** —
across all 11 costim arms plus a no-costim control, in both CD3-off bispecific (`bi`) and
CD3-on trispecific (`tri`) formats (23 constructs per antigen; 69 solid + 24 heme CD19 for
reference). Each construct is a full whole-body run; the readouts below are count-weighted
fractional cell kill per organ.

## 1. Efficacy is CD3-gated and arm-invariant
Killing is switched on by the CD3 arm, not the costim arm. With CD3 affinity dropped to
1e9 nM (`bi`), whole-body depletion collapses to the baseline; turning CD3 on (`tri`,
KD 40 nM) restores it — a **9.1× gate on EGFR** (solid, penetration-limited) versus **32×
on heme CD19**. Within the CD3-on set, the **no-costim control lands inside the armed pack**
(EGFR tri: armed 0.098–0.115, no-costim 0.099; 1.7% total spread). The costim arm does not
change how much the engager kills. This is the counter-screen premise proven in the model's
own hands: arm choice is a pure toxicity-minimization problem, which is exactly why a CD4
toxicity screen — not an efficacy screen — is the right instrument to choose it.

## 2. On-target/off-tumor toxicity is real and antigen-dependent (TI)
Because EGFR and the CEACAMs are expressed on normal epithelium, each construct carries an
on-target/off-tumor liability. Scoring a therapeutic index (TI = on-tumor kill ÷ worst
off-tumor organ, no-costim baseline):

| TAA | on-tumor | worst off-tumor | TI |
|---|---|---|---|
| EGFR | 0.167 | 0.129 (colon) | **1.29** |
| CEACAM6 | 0.105 | 0.178 (colon) | **0.59** |
| CEACAM5 | 0.105 | 0.194 (colon) | **0.54** |

For both CEACAMs the model kills normal tissue **harder than the tumor** (TI < 1). The
dominant off-tumor organ for all three antigens is **large intestine** — the model places
the toxicity in the tissue where EGFR (cetuximab-associated GI toxicity) and CEACAM CRC
programs actually manifest clinically. This falls out of expression + spatial geometry with
**no fitting**.

## 3. The costim arm widens the off-tumor window; CD28 is the dominant amplifier
Tumor kill is flat across arms (penetration-limited), so any arm effect shows up off-tumor.
Ranking the mean change in summed off-tumor kill (7 normal organs) vs no-costim, averaged
across the three antigens:

- **CD28: +0.093** — a **+40–57% increase** in off-tumor kill (EGFR +56.5%, CEACAM6 +48.6%,
  CEACAM5 +39.7%), with **zero on-tumor benefit**. A singular outlier, ~5× the next arm.
- **CD2: +0.018** — a distant second.
- 4-1BB (+0.014), CD27 (+0.011), CD40 (+0.007) — small positive.
- The activation-induced TNFRSF arms (OX40, GITR, ICOS, DR3, HVEM) and CD226 sit **at or
  below** the no-costim baseline (CD226 net-negative, −0.006).

**Interpretation, and its limit.** The two largest amplifiers (CD28, CD2) are constitutively
expressed on resting T cells, so they can add costimulation to bystander T cells wherever the
drug distributes, including normal gut; the activation-induced receptors are near-absent until
a T cell is engaged and therefore fire mainly at the tumor synapse. The grouped
constitutive-vs-induced contrast is 12–16×, but this is **dominated by CD28 alone** —
CD226 (constitutive) is the *lowest* arm, so the effect is not a clean receptor-class law. The
defensible statement is: **CD28 is a singular off-tumor amplifier, and the activation-induced
arms (including the nomination co-lead 4-1BB) spare the off-tumor window.** Notably, the
model's static resting-copy R_costim under-rates inducible arms *on-tumor* (a disclosed
limitation), but **off-tumor the resting density is physiologically correct** — bystander T
cells genuinely do not express 4-1BB/OX40 — so the off-tumor sparing of inducible arms is
precisely the regime where the model is right.

This independently confirms **gating CD28** — already gated on the GRN/DE axis (z ≈ 12) — now
on a completely separate whole-body toxicity axis, and it quantitatively reproduces the
CD28-superagonist off-tumor amplification signature (TGN1412 mechanism) from expression alone.

## Limitation — CRS axis is uninformative for the solid sweep
The myeloid IL-6 term returns **identically zero for all 36 solid-tumor constructs**
(EGFR/CEACAM5/CEACAM6), versus non-zero (0.29–0.44 pg/mL, under-scaled) for heme CD19. The
IL-6 engine reads the myeloid compartment of the heme/blood pools, which the solid-tumor ABM
does not populate in these runs. **The solid-TAA arm ranking therefore rests entirely on the
kill-based off-tumor TI (§2–3), not on CRS**; the zero should not be read as "no CRS
liability." This is a compartment-coverage gap, not a claim.

## Files
- `FIG_costim_offtumor_window.png` — 3-panel: TI ladder, off-tumor-by-arm, per-arm Δoff ranked by receptor class
- `COSTIM_OFFTOX_ANALYSIS.csv` — per-construct on/off-tumor kill + TI, all 3 solid TAAs × 12 arms
- `EGFR_SWEEP_ANALYSIS.csv`, `COSTIM_SWEEP_ANALYSIS.csv` — per-arm depletion/IL-6 tables
- Spatial overlays: `spatial_resolved/spatialres_cs_cs_tri_{EGFR,CEACAM5,CEACAM6}_*` (432 PNGs, 12 organs each)
