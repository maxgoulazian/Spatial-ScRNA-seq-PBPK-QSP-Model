# 1. Introduction, hypothesis, and the three-axis nomination

## 1.1 The clinical problem: costimulation amplifies whatever CD4 does

CD3-based T-cell engagers (TCEs) recruit polyclonal T cells to a tumor antigen by delivering T-cell-receptor "signal 1" in the absence of a costimulatory "signal 2." T cells activated through signal 1 alone activate transiently, fail to sustain proliferation, and drift toward exhaustion, which bounds the depth and durability of the redirected response. The field's response is to add a costimulatory arm — most commonly CD28 or 4-1BB (TNFRSF9) — and tumor-targeted costimulatory bispecifics are now advanced in preclinical and clinical development.

The unsolved problem is toxicity, and its origin is structural rather than incidental. Because the engaging arm is anti-CD3, an antigen present on all T cells, the molecule cannot avoid engaging CD4⁺ T cells; a costimulatory arm then amplifies whatever program those CD4 cells run. Three mechanisms translate that engagement into harm, and they do not share a single point of control:

1. **Cytokine-release syndrome (CRS).** CD4 helper cells are disproportionate cytokine producers, so broad CD4 activation seeds an IFN-γ/TNF/IL-2 cascade that myeloid cells amplify through IL-6 and IL-1. In CAR-T models, monocyte-derived IL-1 and IL-6 are the required amplifiers of CRS (Norelli et al. 2018, *Nat. Med.*, doi:10.1038/s41591-018-0036-4; Giavridis et al. 2018, *Nat. Med.*, doi:10.1038/s41591-018-0041-7, PMID 29808005). This axis scales with the magnitude and systemic breadth of activation, not with the identity of the costimulatory receptor alone.
2. **Regulatory-T-cell (Treg) expansion.** Tregs are CD4⁺ and express CD28 and 4-1BB at high levels, so a costimulatory arm risks expanding the immunosuppressive compartment — the one origin of toxicity that is genuinely a CD4-subset problem and therefore addressable by receptor choice.
3. **Costimulation-intrinsic toxicity, independent of subset.** The CD28 superagonist TGN1412 caused near-fatal CRS in healthy volunteers (Suntharalingam et al. 2006, *NEJM*, doi:10.1056/NEJMoa063842; PMID 16908486), and 4-1BB agonists such as urelumab cause dose-limiting hepatotoxicity via 4-1BB on liver-resident myeloid cells (Segal et al. 2017, *Clin. Cancer Res.*, doi:10.1158/1078-0432.CCR-16-1272; PMID 27756788).

The central premise of this work follows from that decomposition. **The enemy is not the CD4 lineage but a specific CD4 sub-program — the Treg-suppressive plus high-cytokine wiring.** CD4 help is beneficial: it licenses dendritic cells, sustains CD8 cytotoxic differentiation and memory, and CD4 cells can act as direct effectors in some tumors (Borst et al. 2018, *Nat. Rev. Immunol.*, doi:10.1038/s41577-018-0044-0; PMID 30057419). A design that deletes or blocks CD4 help would forfeit that benefit; the goal is to spare the help-and-effector wiring while starving the suppressive-and-cytokine wiring. Of the three toxicity origins, only Treg expansion is soluble at the receptor layer; CRS breadth is magnitude-driven and myeloid-amplified, and costimulation-intrinsic toxicity is subset-independent. Receptor selection therefore addresses roughly one-third of the toxicity problem directly and de-risks a second third (CRS wiring) — the residual (magnitude control, coincidence-gated delivery) is left to the delivery format and to the downstream quantitative model, and is out of scope here.

## 1.2 Why the CD4 Perturb-seq is the right instrument for the counter-screen

The apparent paradox — using a CD4 screen to nominate a costimulatory arm whose therapeutic payload is CD8 killing — resolves once the screen is read as a **toxicity counter-screen** rather than an effector screen. The suppressive and cytokine-release programs that gate costim-engager toxicity are precisely the programs a CD4⁺ compartment expresses and resolves. CD4 is the one compartment in which effector/help wiring can be separated from suppressive/cytokine-release wiring at single-cell, genome-scale resolution.

The anchor dataset is the Marson/Pritchard genome-scale CD4⁺ CRISPRi Perturb-seq screen (Zhu et al. 2025, bioRxiv, doi:10.64898/2025.12.23.696273), which perturbs 11,526 regulators, reads out 10,282 genes, and resolves three activation timepoints. It scores each regulator's effect on the IL-10/Treg-suppressive program and across the storm-cytokine set (TNF, IL-2, IFN-γ), making it a direct map of the two receptor-soluble liability programs. Two orthogonal screens complete the instrument: the Marson-lab CD8 CRISPRa/CRISPRi Perturb-seq (Schmidt et al. 2022, *Science*, doi:10.1126/science.abj4008), whose gain-of-function CRISPRa arm sorted on IFN-γ in CD8⁺ cells and directly identified costimulatory TNF-receptor-superfamily members that raise IFN-γ when overexpressed; and the SLICE genome-wide CRISPR screen with a cancer-cell-killing readout (Shifrut et al. 2018, *Cell*, doi:10.1016/j.cell.2018.10.024).

An honest scope note is load-bearing here. CRISPRi is loss-of-function, whereas an engager costimulatory arm is gain-of-function, so no knockdown screen can *prove* that agonizing a receptor helps. The screens supply a validated, directional state change per receptor; the direction of benefit is anchored on the Schmidt CRISPRa gain-of-function data, and translation of the validated state change into a predicted efficacy and therapeutic window is reserved for the quantitative-systems-pharmacology (QSP) layer downstream of this section.

## 1.3 The three-axis scoring logic and the liability veto gate

Each candidate costimulatory receptor is scored on three axes:

- **Effector benefit** — does agonizing it enhance CD8 effector function? Read from the Schmidt CD8 CRISPRa anchor as a signed effector z-score (`E_schmidt_z`), the gain-of-function evidence a CD4 CRISPRi knockdown screen structurally cannot supply.
- **Suppression liability (SUPP)** — does the receptor's CD4 wiring feed the IL-10/Treg-suppressive program? Read from the CD4 Perturb-seq differential-expression matrix (15-gene Treg/IL-10 signature).
- **CRS liability (CRS)** — does its CD4 wiring drive the storm cytokines (TNF, IL-2, IFN-γ)? Read from the same matrix.

The nomination rule is a **liability veto, not a trade-off.** A receptor is scored on a six-axis liability gate — CRS, SUPP, HELP-erosion, PROLIF (cell-cycle), and EXH (exhaustion), plus a data-driven suppression axis (DD_SUPP) that adds no independent gate call — and **any** liability-up axis eliminates the arm. Effector benefit never offsets a liability: a high effector score cannot rescue a receptor that feeds suppression, CRS, help-erosion, uncontrolled proliferation, or exhaustion. This ordering is deliberate. Reading a composite "therapeutic-window" rank *without* applying the veto first reorders hard-gated arms above the clean survivors and resurrects a "CD28-looks-acceptable" artifact; the gate is therefore applied first, and only the clean survivors are ranked on effector benefit.

## 1.4 The nomination

Applying the veto to the scored panel of 11 costimulatory receptors leaves exactly two clean arms — **4-1BB (TNFRSF9) and CD27** — and they are the co-leads. The table below is the full effector-ranked panel with the canonical gate status for each arm (values from `COSTIM_FINAL_3AXIS_SCORE_v7.csv`).

| Rank | Arm | Gene | Effector z (`E_schmidt_z`) | Gate status |
|---:|---|---|---:|---|
| 1 | CD28 | CD28 | 12.11 | **GATED** [CRS, SUPP, PROLIF] |
| 2 | **CD27** | CD27 | 4.28 | **CLEAN** |
| 3 | **4-1BB** | TNFRSF9 | 3.74 | **CLEAN** |
| 4 | CD30 | TNFRSF8 | 3.22 | GATED [HELP, PROLIF] |
| 5 | CD40 | CD40 | 2.65 | GATED [HELP, PROLIF] |
| 6 | OX40 | TNFRSF4 | 2.07 | GATED [SUPP, EXH] |
| 7 | DR3 | TNFRSF25 | 1.58 | GATED [SUPP, EXH] |
| 8 | DNAM1 | CD226 | 0.64 | GATED [SUPP, EXH] |
| 9 | GITR | TNFRSF18 | 0.09 | GATED [SUPP] |
| 10 | HVEM | TNFRSF14 | 0.02 | GATED [SUPP, EXH] |
| 11 | ICOS | ICOS | −0.39 | GATED [HELP, PROLIF] |

The single most informative row is **CD28**. It carries the top effector score in the entire panel (`E_schmidt_z` = 12.11) and is nonetheless gated out on three liabilities (CRS, SUPP, PROLIF). This is the clearest demonstration that the counter-screen — not the effector axis — drives the nomination: the arm that would win any effector-only ranking is the arm the toxicity gate most decisively rejects, consistent with the TGN1412 experience with CD28 superagonism.

The two clean co-leads are effector-beneficial and, on complementary measures, actively CD8-biased rather than merely liability-neutral:

- **CD27** is the more effector-differentiated on a per-gene basis (CD8-minus-CD4 effector differential +3.40) and leads on receptor **expression** selectivity (CD8−CD4-conventional = +0.46 vs +0.04 for 4-1BB and −0.09 for CD28; `A34_expression_selectivity_CD8vsCD4.csv`). Its CD4 wiring is suppression-*negating* (SUPP down, BH-q = 6.3×10⁻⁵), i.e. it reads as anti-suppressive in the CD4 compartment, and it is help-spared.
- **4-1BB** is fully clean on CRS/SUPP/HELP (all non-significant) and leads on effector-**network** selectivity: its CD8 effector drive is 57.4 versus 3.0 in CD4 (net +54.4, ≈19×; `A29_network_effector_selectivity.csv`). Both co-leads also pass the data-driven suppression axis (DD_SUPP BH-q: 4-1BB 0.79, CD27 0.42), which raises no independent gate call.

The two arms are CD8-biased by *different* measures — CD27 by where it is expressed, 4-1BB by how its network is wired — so the positive criterion is stronger than cleanliness alone: the costimulatory boost is directed onto the killers rather than merely kept off the brakes. Each carries a residual liability that is out of the CD4 screen's field of view and is handed to the downstream model: 4-1BB agonism drives CD8/NK bystander proliferation and liver-myeloid hepatotoxicity (the urelumab class), and CD27 is CRS-capped at monotherapy-equivalent dose in the network read. The nomination that these axes hand forward is therefore a **liability-clean co-lead pair, 4-1BB + CD27**, to be dose- and window-characterized by the QSP layer and to be pursued at agonist-compatible epitopes by the binder-design track.
