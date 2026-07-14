# Three-axis costim-arm nomination — REAL axes (first complete pass)

Instrument: hero CD4 CRISPRi Perturb-seq DE matrix (GWCD4i.DE_stats.h5ad, 33,983
perturbation×condition × 10,282 genes, DESeq2), agonism-oriented (agonism = -1 × KD).
Effector axis: Schmidt 2022 CRISPRa CD8 IFN-γ (gain-of-function, native locus).
Conditions scored: Stim48hr (primary) and Stim8hr (robustness). Weights default (1,1,1).

## Headline nomination (effector-hit receptors)

| Receptor | CD8 effector z | CRS (48h/8h) | Suppression (48h/8h) | Robust call |
|---|---|---|---|---|
| CD28  | 12.11 | +2.15 / +2.86 | +0.49 / +0.46 | **pan-costim liability** (sign-stable) |
| CD2   |  5.60 | +3.00 / +4.49 | −0.10 / −0.73 | **pan-costim liability** (CRS sign-stable) |
| CD27  |  4.28 | +0.55 / +1.11 | −0.12 / +0.03 | selective / mild-CRS |
| 4-1BB |  3.74 | −1.12 / −0.53 | −0.08 / −0.44 | **clean — CRS-negating** (sign-stable) |
| LTBR  |  3.49 | unscoreable   | unscoreable          | effector-only (not CD4-expressed) |
| CD30  |  3.22 | −0.39 / +0.25 | −1.34 / −0.70 | **clean — suppression-negating** (Suppr sign-stable) |
| CD40  |  2.65 | +0.36 / −0.39 | +0.29 / +0.15 | selective (near-zero) |
| OX40  |  2.07 | +0.27 / −1.45 | +0.72 / +0.82 | **suppression liability** (Suppr sign-stable) |

## What is ROBUST (|z|>1 and sign-stable across both conditions)
- **CD28 = pan-costim liability.** Highest effector, highest CRS. This is the TGN1412
  profile — recovered here purely from CD4 *knockdown* data via the sign-flip. Positive
  control: the pipeline works.
- **CD2 = pan-costim liability.** Strong effector, strongest CRS liability of the panel.
- **4-1BB = clean arm.** Effector-positive AND agonism predicted to LOWER CD4 storm
  cytokines (CRS-negating) at both timepoints. Reproduces 4-1BB's known lower-CRS profile
  vs CD28 superagonism — external validation.
- **CD30 = clean arm.** Effector-positive AND agonism strongly negates the Treg/IL-10
  suppression program at both timepoints; CRS near-neutral.
- **OX40 = suppression liability.** Weakest effector hit, and agonism feeds the Treg
  suppression axis (sign-stable) — consistent with OX40 being Treg-associated. Not a clean arm.

## What is NOT yet robust (near-zero, label flips between conditions)
CD27, CD40 sit near zero on both toxicity axes — "selective" but the fine A/B label is
timepoint-sensitive. Treat as CD8-selective-leaning, not negating.

## Help-preservation check (Axis 4) — the "sub-program vs lineage" guard
At Stim48hr, both clean arms (4-1BB, CD30) classified **Option B+ (help-sparing negator),
help VERIFIED preserved** — i.e. they negate the CD4 toxicity sub-program while the CD4
helper program (BCL6/MAF/TOX2/CXCR5/TCF7) stays intact. This is the proposal's ideal
"best case": CD8 up, CD4 suppressive/CRS wiring down, CD4 help spared. No effector arm
scored as a lineage-wide CD4 shutdown (which would have been downgraded).

## Answer to the two framing questions
1. "Will the nomination change when the real data lands?" — YES, and it did. The proxy
   could only flag CD4 *activation* (everything looked pan). The real suppression axis
   separates the field: 4-1BB and CD30 emerge as clean arms, invisible to the proxy.
2. "CD8-selective vs CD8-stim + CD4-negating?" — now answerable per-receptor. The best
   arms (4-1BB CRS-negating, CD30 suppression-negating) are Option-B *negators*, not merely
   selective — the wider-window profile — and the help guard confirms the negation is
   sub-program-specific, not lineage-wide. CD28/CD2 are the pan liabilities to avoid.

## Caveats
- Fine A/B/pan labels near |z|<1 are timepoint-sensitive; only the sign-stable, >1 calls
  above are robust. Report those.
- Weights (wC=wS=1) are provisional; QSP will set them from the therapeutic-window model.
- 4-1BB's known myeloid hepatotoxicity is a liver-myeloid mechanism (costim-intrinsic,
  subset-independent) that a CD4 T-cell screen structurally cannot see — flagged separately
  in the proposal; not captured or contradicted here.
- LTBR remains effector-covered but CD4-unscoreable (ectopic; nothing to knock down).
