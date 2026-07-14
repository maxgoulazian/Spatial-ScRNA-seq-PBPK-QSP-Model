# Project framing — receptor **and** format (v2)

**One-line thesis.** Costim-engager toxicity is *two* problems — a **wiring** problem (which costim program a receptor feeds) and a **where/when** problem (whether signal 2 fires only where signal 1 does). A CD4⁺ single-cell counter-screen solves the first; a per-cell spatial QSP model solves the second. The "wrong subset" screen turns out to be the right instrument for the wiring half.

---

## Why reframe

The original framing is **accurate and should be kept at its core**, but it is a *single-lever* framing (pick a clean receptor) describing a *two-lever* result (clean receptor **+** coincidence-gated format). The strongest things we found — active CD8-selectivity and the cis/trans format axis — have no home in the single-lever version.

### Original opening (for reference)
> CD3 engagers deliver signal 1 but no signal 2 — so the field adds a costimulatory arm. The unsolved problem is toxicity. Because the engaging arm is anti-CD3 — present on every T cell — the molecule cannot avoid CD4⁺, and costimulation amplifies whatever CD4⁺ does. The enemy is a CD4 sub-program — the Treg-suppressive, high-cytokine wiring — not the CD4 lineage. CD4 help is beneficial.

### What is load-bearing and stays (verbatim)
- **"The enemy is a sub-program, not a lineage."** The crown jewel — counterintuitive, correct, differentiating. The nomination separates effector wiring from suppressive/CRS wiring *within* CD4, and it held (4-1BB + CD27 clean on all six liability axes, effector-beneficial, survived the unbiased genome-wide start).
- **"CD4 help is beneficial."** Separates us from the naive Treg-depletion crowd.
- **Toxicity as the unsolved problem** (TGN1412, urelumab hepatotox).

## What the single-lever framing undersells

1. **Our positive criterion is stronger than "clean."** The sharpest results are *active CD8-selectivity*, not just liability-cleanliness:
   - **CD27** leads on receptor **expression** selectivity (CD8−CD4 ≈ +0.46; 4-1BB ≈ +0.04; CD28 negative).
   - **4-1BB** leads on effector-**network** selectivity (~19× CD8-vs-CD4 drive; net +54).
   - → The co-leads are CD8-biased by **complementary** measures (CD27 by *where expressed*, 4-1BB by *how wired*). "The boost lands on the killers" beats "the boost doesn't hit the brakes."

2. **Toxicity is multi-origin, and only one origin is receptor-soluble.** Of the three mechanisms, only **Treg expansion** is truly a CD4-subset problem. **CRS breadth** is magnitude-driven / myeloid-amplified; **costim-intrinsic toxicity** (TGN1412, 4-1BB hepatotox) is subset-independent. So "the enemy is a CD4 sub-program" is honestly ~⅓ of the toxicity. The other ⅔ is a *where/when/how-much* problem — magnitude control, coincidence-gating, tumor-conditional delivery — which is exactly the **cis/trans QSP axis**. In the single-lever framing that work is a homeless bolt-on.

---

## The reframe: two-part nomination

**Part 1 — the receptor.** Sub-program-clean **and** CD8-selective → **4-1BB / CD27**. Scored by the CD4 Perturb-seq counter-screen (6 liability axes + effector benefit + expression/network selectivity).

**Part 2 — the format.** Deliver it coincidence-gated so signal 2 fires only where signal 1 does → **cis geometry / tumor-conditional**. Scored by the per-cell spatial QSP model (cis/trans span sweep, Treg-axis differential).

**Unifying line.** *Costim-engager toxicity is a wiring problem and a where/when problem; single-cell data solves the first, the QSP model solves the second, and the "wrong subset" screen turns out to be the right instrument for the wiring half.*

### Three evidence pillars
| Pillar | Claim | Primary instrument |
|---|---|---|
| **Wiring (receptor)** | 4-1BB/CD27 are sub-program-clean + CD8-selective | CD4 Perturb-seq counter-screen; Schmidt CD8 effector anchor |
| **Where/when (format)** | Cis (height-matched) gating narrows the therapeutic window vs trans | Per-cell spatial QSP, cis/trans span sweep |
| **Buildability** | Both co-leads are structurally druggable at agonist-compatible epitopes | RFdiffusion/RFantibody de novo binders + clustering-compatibility check |

## Positioning notes
- **Supporting, not headline:** the redirector×costim co-expression axis (CD3 vs alternative redirectors) is a *third* selectivity layer — keep as a supporting result, it's the least mature of the three.
- **Instrument-forward variant** ("the wrong-subset screen is the right instrument") is elegant but abstract — use as a mid-talk reveal, not the cold open. Lead with the therapeutic payoff.
- **Honest scope carried forward:** CRISPRi is loss-of-function while an engager arm is gain-of-function (no screen proves agonism helps — QSP translates the validated state change); CD4 is the counter-screen, not the effector axis; cis epitope is hypothetical (no designed binder yet), CRD1 is the one with an actual design.

*This framing is tighter and more novel than the original, not looser: it is honest about multi-origin toxicity, promotes the QSP cis/trans work to a co-equal pillar, and keeps the counterintuitive core intact.*