# 3. Effector benefit and expression selectivity

The counter-screen (Sections 1–2) establishes which candidate costimulatory arms are *quiet* on the CD4⁺ suppressive (Treg/IL-10) and cytokine-release (CRS) programs. That is a necessary but not sufficient condition for a nomination: an arm that is liability-clean but does nothing for CD8⁺ effector function is not worth building. This section supplies the positive axis — evidence that agonizing the receptor **enhances CD8⁺ effector output** — and confirms that the two liability-clean co-leads, **4‑1BB (TNFRSF9)** and **CD27**, are wired and expressed *selectively* toward the CD8⁺ effector compartment rather than the CD4⁺ subsets a CD3-engaging arm would co-recruit.

## 3.1 Why a gain-of-function CD8 anchor is required

The hero CD4⁺ Perturb-seq map (Zhu et al., 2025) is a CRISPR **interference** (loss-of-function) screen. Knocking a gene down reports whether it is *required* for a program under the screen conditions; it cannot report whether *over*-expressing or agonizing the receptor would *push* the effector program. A therapeutic costim arm is a gain-of-function agonist, so the direction the engager will drive is structurally invisible to a knockdown screen. The agonism-direction evidence must come from a matched gain-of-function readout.

We anchor that axis on **Schmidt et al., 2022** (*Science*, DOI [10.1126/science.abj4008](https://doi.org/10.1126/science.abj4008)), the Marson-lab genome-scale CRISPR-activation (CRISPRa) and CRISPRi platform in primary human T cells, which sorted CD8⁺ cells on IFN‑γ production. The anchor is uniquely fit for purpose because its CRISPRa arm made exactly the gain-of-function observation a CD4 CRISPRi screen cannot: CRISPRa "selectively detected a set of tumor necrosis factor superfamily receptors … including 4‑1BB, CD27, CD40, and OX40" that "were not individually required for signaling … but could promote IFN‑γ when overexpressed." That is the agonism-direction signal — a costim receptor whose effect appears on overexpression and is silent on knockdown — read directly from primary human CD8⁺ T cells.

## 3.2 Per-arm effector benefit

We score each of the 11 candidate arms by its Schmidt CD8 CRISPRa IFN‑γ activation z (Figure 3a, x-axis). The ordering is:

| Arm | Effector z (Schmidt CD8 CRISPRa) | Status |
|---|---|---|
| CD28 | **12.11** | GATED [CRS, SUPP, PROLIF] |
| CD27 | **4.28** | **CLEAN — co-lead** |
| 4‑1BB (TNFRSF9) | **3.74** | **CLEAN — co-lead** |
| CD30 | 3.22 | gated |
| CD40 | 2.65 | gated |
| OX40 | 2.07 | gated |
| DR3 | 1.58 | gated |
| DNAM1 | 0.64 | gated |
| GITR | 0.09 | gated |
| HVEM | 0.02 | gated |
| ICOS | −0.39 | gated |

CD28 is the strongest effector hit by a wide margin (z = 12.11), but it is vetoed on every liability axis (Section 2) — the same biology that makes CD28 a potent activator makes it the CRS and Treg-expansion hazard that the CD28 superagonist **TGN1412** realized clinically, causing near-fatal cytokine-release syndrome in six healthy volunteers (Suntharalingam et al., 2006, *NEJM*, DOI [10.1056/NEJMoa063842](https://doi.org/10.1056/NEJMoa063842)). Among the liability-clean arms, **CD27 and 4‑1BB are the two highest effector scores** (z = 4.28 and 3.74), sitting immediately below CD28 and above every other candidate. The nomination therefore recovers near-CD28-level effector benefit while discarding the CD28 liability profile.

![Figure 3. Effector benefit and CD8-selective wiring of the liability-clean co-leads. (a) Effector benefit (Schmidt 2022 CD8 CRISPRa IFN-γ activation z, x-axis) versus expression selectivity (CD8 minus CD4-conventional single-cell expression, RTCC CITE-seq GSE292621; y-axis). Liability-clean co-leads 4-1BB and CD27 in focal color; gated arms (liability veto, Section 2) in grey. CD28 is the top effector score but is expression-nonselective (below zero) and gated. (b) Network effector-hub selectivity (A29): signed GRNBoost2 importance onto the CD8 cytotoxic-effector target set (GZMB/PRF1/IFNG/…) in a CD8-gain-of-function network minus the CD4-gain-of-function network. 4-1BB is the dominant CD8-selective effector hub (+54.4); CD27 is near-neutral on this axis (+0.49) and carries its selectivity through expression instead.]({{artifact:art_30f87d11-5b2b-40a3-9592-a7402a2041e4}})

## 3.3 Network effector-hub selectivity

Effector z scores the receptor's *own* IFN‑γ effect. To ask whether the receptor sits *upstream of the killing program* and does so *selectively in CD8⁺ cells*, we scored network topology (A29): the summed signed regulatory importance (GRNBoost2) flowing from each receptor onto the CD8 cytotoxic-effector target set (GZMB, PRF1, IFNG and related genes) in a CD8-gain-of-function gene-regulatory network minus the same quantity in a CD4-gain-of-function network. The signed difference, *net effector selectivity*, is positive when the receptor drives the cytotoxic program more strongly in CD8⁺ than in CD4⁺ cells (Figure 3b).

**4‑1BB is the dominant CD8-selective effector hub: net effector selectivity = +54.4**, the largest of any arm and far above the next candidate (ICOS, +8.3). Its raw drive onto the cytotoxic-effector targets is **~19-fold higher in CD8⁺ than in CD4⁺ cells (57.4 vs 3.04)**, and its overall effector-hub connectivity is likewise CD8-weighted (net hub selectivity +155.0; CD8 hub strength 173.1 vs CD4 18.0). This is network-topology evidence, independent of the per-cell expression and per-gene effector axes, that 4‑1BB agonism lands on the CD8 killing program specifically.

CD27 behaves differently and instructively: its net effector selectivity is **+0.49 — essentially neutral** — and its overall hub connectivity is actually CD4-weighted (net hub selectivity −59.0). CD27 is therefore *not* a network-level CD8-selective effector hub; its subset selectivity is carried by the expression axis (Section 3.5), not the network. The two co-leads are thus selective by **complementary and non-redundant mechanisms** — 4‑1BB through effector-network topology, CD27 through compartment-restricted expression — which makes the joint nomination robust to the failure of any single axis. (CD27's CD4-side network connectivity was screened against the suppressive/CRS programs in Section 2 and cleared; it is not re-litigated here.)

## 3.4 Concordance 2×2 and nomination robustness

We cross-checked the anchor against a second, independent gain-of-function readout and against the loss-of-function frame (A30). Each arm was placed in a 2×2 on an effector axis (CD8/IFN‑γ) and an activation axis (CD4/IL‑2), scoring a CRISPRa gain-of-function z (GOF) and a CRISPRi loss-of-function "agonism-frame" z (LOF):

- **4‑1BB: GOF z = +3.75, LOF-agonism z = −0.01 → "LOF-uninformative."** This is the counter-screen's structural limitation made quantitative: 4‑1BB is a clear gain-of-function effector hit whose knockdown carries no signal — precisely the "not required, but promotes IFN‑γ when overexpressed" behavior Schmidt reported. It is the concrete proof that the CD4 CRISPRi map cannot, by construction, supply 4‑1BB's agonism-direction evidence, and that the Schmidt GOF anchor is required.
- **CD27: GOF z = +3.75, LOF-agonism z = +2.28 → "concordant."** CD27 shows effector signal in both frames.
- Both co-leads are **quiet on the CD4/IL‑2 activation axis** (activation GOF z = 0.30 for 4‑1BB, 0.24 for CD27), in sharp contrast to CD28 (activation GOF z = 10.90). The co-leads raise CD8 effector output without raising the CD4 activation axis that seeds CRS.

A 20,000-draw Monte-Carlo re-check propagated the effector-z uncertainty and stress-tested the co-leads against the gated but effector-competitive challengers (DR3, HVEM, CD40). Under effector confidence alone, the co-leads' joint probability of occupying the top two effector slots is **0.528**; **weighting by CD8 expression breadth raises it to 0.836**, and CD27's probability of the single top slot is **0.74**. The main challenger, CD40, is both gated and expression-suppressed (CD8 breadth 9.5%, expression selectivity −0.064), so it churns under effector CI alone but collapses once expression is accounted for. CD27 also carries the tighter effector estimate (z SD = 0.74 vs 2.38 for 4‑1BB), consistent with its higher expression and breadth.

## 3.5 Expression selectivity: CD8 versus CD4

Subset selectivity of a costim arm is only real if the receptor is preferentially expressed or wired on CD8⁺ effectors rather than on the CD4⁺/Treg compartments a CD3 arm co-engages. We confirmed the expression differential in an orthogonal single-cell surface-and-transcript dataset — RTCC CITE-seq (GSE292621, 25,348 cells; gating CD8 = 5,412, CD4-conventional = 3,156, Treg = 2,262). A34 reports the CD8-minus-CD4 expression selectivity (Figure 3a, y-axis):

| Arm | CD8 − CD4conv | CD8 − Treg | CD8 − mean(tox compartments) | CD8 breadth (% positive) |
|---|---|---|---|---|
| **CD27** | **+0.464** | +0.400 | +0.432 | **84.4%** |
| **4‑1BB** | **+0.038** | +0.104 | +0.071 | 57.5% |
| CD28 | **−0.089** | +0.120 | +0.015 | 86.3% |

**CD27 is strongly CD8-selective on expression (+0.464 over CD4-conventional, +0.400 over Treg)** and broadly expressed on CD8⁺ cells (84.4% positive), the mirror image of its neutral network-effector profile — it is the expression-selective co-lead. **4‑1BB is near-parity on bulk expression (+0.038)** but, as Section 3.3 showed, is the network-effector-selective co-lead; its per-compartment levels (CD8 1.04, CD4conv 1.00, Treg 0.93) show a modest CD8/Treg gradient consistent with its inducible-costim biology. **CD28 is expression-nonselective, trending CD4-biased (−0.089)** — it is expressed across all three compartments at high breadth (86% CD8, 94% CD4conv, 94% Treg), which is the expression correlate of its Treg-expansion and CRS liability. This is the expression-layer confirmation that CD28's problem is not merely functional but distributional: it cannot be engaged CD8-selectively.

The two co-leads are thus selective on *different* axes — CD27 on expression, 4‑1BB on effector-network topology — with neither showing the pan-compartment distribution that makes CD28 unsuitable. This is the subset-selectivity requirement the whole strategy rests on: a CD3-engaging molecule co-recruits CD4⁺ cells, so a costim arm de-risks toxicity only if its benefit is wired or expressed toward CD8⁺ effectors. Both nominated arms satisfy that condition, by complementary means.

## 3.6 Summary

The gain-of-function CD8 anchor supplies the agonism-direction evidence a CD4 CRISPRi counter-screen structurally cannot: among liability-clean arms, **CD27 (effector z = 4.28) and 4‑1BB (z = 3.74) are the two strongest CD8 effector hits**, recovering near-CD28-level benefit without the CD28 liability profile. Their subset selectivity is confirmed on two independent axes and is non-redundant — **4‑1BB is the dominant CD8-selective effector-network hub (+54.4, ~19× CD8-vs-CD4 drive), CD27 is the CD8-selective on expression (+0.46 over CD4)** — and both are quiet on the CD4/IL‑2 activation axis. Nomination robustness holds under Monte-Carlo effector-uncertainty and rises when expression breadth is weighted in (joint top-two probability 0.528 → 0.836). The effector-plus-selectivity axis is consistent with the residual clinical liability of pan-costim agonism — for example the dose-dependent hepatotoxicity of the 4‑1BB agonist urelumab, where significant transaminitis was strongly associated with doses ≥1 mg/kg (Segal et al., 2017, *Clin. Cancer Res.*, DOI [10.1158/1078-0432.CCR-16-1272](https://doi.org/10.1158/1078-0432.CCR-16-1272)) — which motivates the tumor-conditional delivery format that layers orthogonally on top of receptor choice (Section 5).

---

*Data sources (verified against artifacts): A29 network effector selectivity; A30 2×2 concordance classification and nomination re-check; A34 expression selectivity (RTCC CITE-seq GSE292621); A21 per-compartment expression breadth. Effector anchor: Schmidt et al., 2022 (Science, 10.1126/science.abj4008). Effector-z estimates are gain-of-function CRISPRa activation z; the Schmidt CD8 IFN‑γ anchor is the primary axis (Section 3.2), cross-checked in A30 against an independent CRISPRa GOF / CRISPRi LOF 2×2 (Section 3.4).*
