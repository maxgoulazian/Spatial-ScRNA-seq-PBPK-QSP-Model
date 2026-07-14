## 6. Tumor-antigen (TAA) selection for the redirector arm

The costimulatory nomination (Sections 1–5) fixes *which* second signal the molecule
should deliver — the liability-clean 4-1BB (TNFRSF9) + CD27 co-leads — but a costim-armed
CD3 engager is only as safe as the coordinate that tells it *where* to fire. In the
tumor-conditional format this report adopts, the tumor-associated antigen (TAA) arm sets
that coordinate twice over: it is the binding event that redirects the T cell onto the
malignant cell, and, in a coincidence-gated geometry, it is the localizer that restricts
delivery of signal 2 to the tumor bed. On-target/off-tumor risk for the whole construct
is therefore inherited from the TAA arm's healthy-tissue expression, which makes TAA
selection a safety decision, not only an efficacy one. This section nominates the TAA
slate; it does not design binders against it (see scope note, §6.5).

### 6.1 Antigen atlas and balanced-rank scoring

Candidate antigens were scored on a human colorectal-cancer (CRC) single-cell reference —
Lee et al. 2020 (*Nat Genet*, doi:10.1038/s41588-020-0636-z), comprising two independent
patient cohorts, SMC (GSE132465) and KUL3 (GSE144735). CRC is the disease anchor because
it is where CD4⁺ effector help and direct CD4 cytotoxicity are documented, and where the
CEA-directed engager class already has clinical footing (§6.3). Each antigen carries a
composite **balanced-rank score** (`score_balanced_v2`) integrating three verifiable axes
from the atlas: tumor-cell **coverage** (fraction of malignant cells expressing),
**selectivity** (`sel_mean`, malignant- vs non-malignant-cell contrast within the atlas),
and **tumor-restriction** (`z_restriction`, a z-score of tumor expression against a healthy
tissue atlas). Surface accessibility (`is_surface`) and reproducibility across both cohorts
(`replicate`) were applied as gates. The six retained finalists occupy balanced ranks 1, 2,
3, 8, 11 and 12 — a curated slate rather than a raw top-6, with the three lower-ranked
members retained for specific and stated reasons (§6.4).

### 6.2 The finalist slate

All six finalists are surface-expressed (`is_surface = True`) and five of six reproduce
across the SMC and KUL3 cohorts (TSPAN6 is the single non-replicating member, retained on
safety grounds). By clinical maturity the slate splits 3 clinical / 2 preclinical / 1
untargeted.

| Gene | Balanced rank | Restriction (z) | Tumor coverage | Selectivity | Top healthy location | Clinical maturity |
|------|:---:|:---:|:---:|:---:|------|------|
| *CEACAM6* | 1 | +2.37 | 0.63 | 0.43 | Skin: epithelium | Clinical — dual CEACAM5/6 ADC (EBC-129, FDA Fast Track 2025) |
| *CEACAM5* | 2 | +2.39 | 0.56 | 0.14 | Large intestine: epithelium | Clinical — CEA×CD3 TCE (cibisatamab) + ADCs (M9140, tusamitamab rav., SGN-CEACAM5C) |
| *ITGB4* | 3 | +1.39 | 0.50 | 0.26 | Heart: neural | Preclinical — anti-CD3/ITGB4 bispecific-armed T cells (mouse) |
| *LY6E* | 8 | −0.20 | 0.57 | 0.59 | Ovary: epithelium | Clinical — ADC (DLYE5953A, Ph1); **no TCE** |
| *TSPAN6* | 11 | +1.05 | 0.23 | 0.10 | Ovary: epithelium | Untargeted (novel) |
| *DPEP1* | 12 | +1.97 | 0.20 | 0.37 | Small intestine: epithelium | Early/preclinical — CRC-associated, GPI-anchored |

*Values read from `TAA_finalists_6.csv`. Restriction z is against a healthy-tissue atlas;
coverage and selectivity are within-CRC-atlas metrics (Lee et al. 2020).*

Bindability, from the extracellular-domain architecture: CEACAM5 and CEACAM6 are rated
**EASY** (immunoglobulin-domain ECDs); ITGB4 (large integrin β4 ECD), LY6E (small,
GPI-anchored) and DPEP1 (GPI-anchored dipeptidase, accessible ECD) are **MODERATE**; TSPAN6
is **HARD** (a four-transmembrane tetraspanin presenting only a small EC2 loop).

The two-axis landscape (Figure 6.1) separates the slate cleanly. The CEACAM pair sits at
the top-right — high restriction *and* high coverage — while DPEP1 trades coverage for the
third-highest restriction, ITGB4 is intermediate on both, TSPAN6 is a low-coverage/moderate-
restriction outlier, and LY6E is the informative failure of the restriction axis: highest
selectivity within the tumor atlas (0.59) yet a *negative* healthy-atlas restriction z
(−0.20), meaning it is not tumor-restricted against normal tissue despite discriminating
malignant from non-malignant cells inside the tumor.

![Figure 6.1. Balanced-rank CRC redirector finalists positioned by tumor-restriction (z vs. healthy atlas) against tumor-cell coverage, colored by clinical maturity. The CEACAM5/CEACAM6 anchor pair occupies the high-restriction/high-coverage quadrant; LY6E falls left of the semantic-zero restriction line despite high within-tumor selectivity.]({{artifact:1ac9574d-44d5-4187-88cb-bf30353ce928}})

### 6.3 CRC anchor logic: CEACAM5/CEACAM6

CEACAM5 (carcinoembryonic antigen, CEA) and CEACAM6 are the anchor pair because they lead
every axis simultaneously: the two highest tumor-restriction z-scores in the slate
(CEACAM5 +2.39, CEACAM6 +2.37), high coverage (0.56 and 0.63), and reproducibility across
both cohorts. They also carry the most clinical de-risking, which is what separates the two
in role.

**CEACAM5 is the nominated clinical anchor** despite ranking second on the composite score,
because it is the class's positive control. Cibisatamab (CEA-TCB, RG7802) is a
clinical-stage CEA×CD3 T-cell engager (*Nat Commun* 2024, doi:10.1038/s41467-024-48479-8),
and CEACAM5 additionally has multiple ADC programs (M9140; tusamitamab ravtansine;
SGN-CEACAM5C). Critically, the exact architecture this report proposes — a CEA-directed
CD3 engager paired with tumor-localized 4-1BB costimulation — has now been tested in the
clinic: cibisatamab plus FAP-4-1BBL (a fibroblast-activation-protein-targeted 4-1BB ligand)
in microsatellite-stable metastatic CRC, phase 1b, NCT04826003 (*Nat Med* 2026,
doi:10.1038/s41591-026-04380-z). Across 52 patients the combination was tolerable
(dose-limiting toxicities in 2/52, 3.8%; CRS in 30/52, 57.7%, with grade ≥3 in only 2/52),
produced confirmed partial responses in 7/52 (13.5%), and drove intratumoral CD8⁺ and
CD8⁺Ki67⁺ infiltration — direct external validation that redirector + tumor-localized
costim is buildable and active, and that 4-1BB (a co-lead here) is a rational costim arm
for it.

That same trial also validates the safety concern this section raises. Its dominant
target-related toxicities were gastrointestinal — colitis in 7/52 (13.5%), including
immune-mediated enterocolitis and one fatal CMV colitis — which is precisely concordant
with CEACAM5's top healthy-tissue location in our atlas (large-intestine epithelium). The
mapping between a TAA's healthy expression site and the observed on-target/off-tumor organ
toxicity is not incidental; it is the mechanism, and it is why healthy-location is a
first-class column in the finalist table.

**CEACAM6 is the top balanced-rank finalist** (score 9.33 vs. CEACAM5's 7.30), driven by
the highest coverage in the slate (0.63). It is clinically anchored through the dual
CEACAM5/6 ADC EBC-129 (FDA Fast Track, 2025) but has no CD3-engager program of its own,
placing it between anchor and whitespace: a de-risked antigen with an open T-cell-engager
lane.

### 6.4 Whitespace versus crowded targets

The slate is deliberately balanced between de-risked anchors and whitespace. **CEACAM5 is
the crowded target** — its value is precisely that it is crowded, giving a validated
positive control against which a new construct can be benchmarked, but it offers little
differentiation. The remaining finalists open whitespace of different kinds:

- **DPEP1** — the strongest whitespace nomination. A GPI-anchored dipeptidase with the
  third-highest tumor-restriction in the slate (+1.97), reproducing across both cohorts,
  CRC-associated (a marker of low- to high-grade intraepithelial-neoplasia transition and
  adverse prognosis; *Br J Cancer* 2013, doi:10.1038/bjc.2013.363) yet only early/preclinical
  as a target and with no T-cell-engager program. Its low coverage (0.20) is the cost of its
  high restriction.
- **ITGB4** — preclinical only (mouse anti-CD3/ITGB4 bispecific-armed T cells), moderate on
  both restriction and coverage; whitespace for a clinical engager but with a healthy-tissue
  flag (heart) that warrants attention.
- **TSPAN6** — untargeted/novel and the safest on paper (lowest healthy-tissue RNA in the
  slate, non-replicating so retained conservatively), but **HARD** to drug: a four-TM
  tetraspanin whose only accessible epitope is the small EC2 loop.
- **LY6E** — an instructive case: an ADC target (DLYE5953A, Ph1) with **no TCE**, so nominally
  whitespace, and the highest within-tumor selectivity in the slate — but its negative
  healthy-atlas restriction z means it is broadly expressed in normal tissue, a poor fit for
  a costim-amplified engager where off-tumor engagement is the central liability.

The through-line: for a costim-armed CD3 engager, whitespace is only worth pursuing where
tumor-restriction against *healthy* tissue is genuine (DPEP1, CEACAM6), not merely where the
antigen discriminates malignant from non-malignant cells inside the tumor (LY6E). CEACAM5
supplies the crowded, de-risked benchmark; DPEP1 and CEACAM6 supply the differentiated,
restriction-genuine whitespace.

### 6.5 Scope note: binder design is deferred for all six

In `TAA_finalists_6.csv` the `rfdiffusion_campaign` field reads *"deferred (stretch goal —
not run now)"* for **all six finalists, CEACAM5 included**. This section nominates and
prioritizes antigens; it does not deliver binders against them. The de novo AF3/RFdiffusion
binder campaign reported elsewhere (Section 7) is a separate track and is not to be read as
a binder result for this TAA slate. Where CEACAM5 is described here as a "positive control,"
that refers to its clinical precedent (cibisatamab) as a benchmark antigen, not to any
in-house binder having been generated for it in this work.
