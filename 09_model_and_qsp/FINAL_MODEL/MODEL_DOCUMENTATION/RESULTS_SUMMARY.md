# TCE Spatially-Resolved Counterscreen — Results Summary

*A whole-body, per-cell mechanistic model of T-cell-engager efficacy and cytokine-release risk.*

---

## What the model is

Every individual cell in a human body, simulated. A **Shah-Betts platform PBPK** backbone (PMID 22143261)
carries the antibody to each organ. Inside each organ, **individual T cells, target cells, regulatory T cells
and macrophages** are placed in space, each carrying **its own receptor copy numbers from single-cell
RNA-seq**. The engager binds them with **literal kinetic multivalent chemistry** — not an occupancy shortcut —
forms immune synapses, kills targets, and, where an engaged T cell is **physically touching** a macrophage,
triggers IL-6 secretion from that individual cell.

**Efficacy and toxicity emerge from the same geometry.** A T cell engaged beside a tumour cell kills it; the
same T cell engaged beside a macrophage causes CRS. **Therapeutic window is one mechanism read two ways.**
That is what would license ranking constructs that have never been built.

---

## What is EMERGENT (not fitted)

These are the claims the architecture actually earns:

| phenomenon | why it matters |
|---|---|
| **Prozone / hook effect** | At high dose both arms saturate *separately* and bridging collapses — the central dose-selection problem for every TCE. **There is no prozone equation in this model.** The trimer term is quadratic in free arms, so it falls out of solving the chemistry honestly. An *imposed* hook would need re-fitting for every new format; an emergent one does not. |
| **Tumour-conditional costim** | 4-1BB / OX40 / ICOS are absent on resting T cells and appear only *after* TCR engagement — which is exactly why the field targets them. Induction is driven by **each T cell's own** engaged fraction, so a T cell that never engages never upregulates. The conditionality is emergent, not imposed. |
| **CRS geography** | A drug that engages T cells *where the macrophages are* causes more CRS. The model knows this because it knows where the cells are. |
| **The linker optimum** | Effective 2nd-arm concentration falls as `1/span`; cleft feasibility is zero below ~0.6× the cleft. Their product peaks at **span = cleft ≈ 13 nm** — so **past the cleft distance, LONGER LINKERS MAKE BRIDGING WORSE.** Counter-intuitive, falsifiable, and free: reach was never the constraint — *concentration* was. |

---

## VALIDATION STATUS — stated honestly

| arm | status |
|---|---|
| **PK** | 22 molecules with digitized clinical PK / popPK / exposure-response. **Independent of the IL-6 arm.** |
| **IL-6 / CRS** | **NOT VALIDATED.** See below. |
| **Costim screen** | **BLOCKED** — the activation-induction fold-upregulation for 4-1BB / OX40 / ICOS / GITR **does not exist in the literature.** The code **fails closed** rather than guess it. |

### The IL-6 test, and why it failed

A pre-registered test, with the threshold set **before** the number was seen:

> Clinical **mosunetuzumab / teclistamab** IL-6 ratio = **152 / 21 = 7.24×** (both population means, both
> PMID-sourced). Clearance cancels from a ratio, so this needs no clearance value.
> **Model: 1.69×. Fold error 4.29×. FAIL** (threshold 3×).

**The goalpost was not moved.**

---

## The three findings that came out of failing

### 1. The model was producing IL-6 its own cells could not make (942× over its physical ceiling)

Nobody had compared the model **to its own physical limit**. Production cannot exceed *(every secretor
macrophage in the census, fully activated) × (the measured per-cell secretion rate)*:

```
census 35,092,161 macrophages × 3.9% secretors × 0.00133 pg/hr  =  1,821 pg/hr   ← CEILING
model reported                                                  =  1,716,379 pg/hr   (942×)
```

**Root cause: a COUNT used as a SCALE.** `organ_myeloid_counts.json` holds raw macrophage *counts*; they were
multiplied in as if they were a *ratio*. The correct scale is `count ÷ sampled_agents`. The code comment
defined it correctly — *"physiological count / sampled count"* — while the loader read the counts verbatim.
**Intent said ratio, data said count, and nothing asserted the difference.**

**A three-line conservation guard now runs every simulation and would have caught this on day one.**
*General lesson: check a model against its own physical limits, not only against data.*

### 2. Fixing that bug made the model worse — because it was masking a second, opposite error

```
corrected ceiling → maximum possible plasma IL-6 =   0.78 pg/mL
clinical mosunetuzumab population mean           = 152    pg/mL
                                                 → the mechanism is 194× TOO LOW
```

**The 942× bug had been compensating for a ~195× shortfall.** Two large errors in opposite directions,
partially cancelling — producing a discrepancy that looked merely "5–20× off," i.e. unremarkable and
calibration-sized, **while being wrong for two independent reasons. Each error hid the other.**

### 3. The model predicts the missing organ — and says how big it must be

The mechanism's output is `macrophages × secretor_fraction × S_MAX × activation`. Secretor fraction and
per-cell rate are **measured**; activation is bounded by 1. **The only free term is the census — and the
census has no liver.**

> ### PRE-REGISTERED PREDICTION (`IL6_PREREGISTERED_PREDICTION_2026-07-13.md`)
> To reach clinical IL-6, the mechanism requires **≥ 6.8 × 10⁹ total macrophages**, against a current census
> of **3.5 × 10⁷**. **Kupffer cells are the largest tissue-macrophage population in the human body and are
> entirely absent from the model.**

---

## The provenance audit — an uncomfortable but important result

The model's clinical IL-6 anchors were audited back to their source sentences. **Every one was wrong, and each
in a different way:**

| value | what it actually was |
|---|---|
| elranatamab **191** | **A PAGE NUMBER.** A dot-leader from the *Table of Figures* of FDA BLA 761345, used as a concentration. |
| elranatamab **340 / 230** | **A digitization of a figure that does not exist.** The cited paper (MagnetisMM-3, *Nat Med* 2023, PMC10504075) has **four figures, no Figure 6, and zero mentions of IL-6.** Verified against the primary source. |
| mosunetuzumab **570** | **No source exists.** The real population mean is **152**. |
| teclistamab **21 vs 288** | **Both real** — but 21 is a population **MEAN** and 288 is the **highest individual patient**, from the same sentence (PMID 38831634), **13× apart**. Ranking one against the other manufactured a fake ~26× clinical dynamic range that the mechanism was then repeatedly "fixed" to chase. |
| IL-6 clearance **0.20/hr** | Cited to PMID 31268236 — **a semi-mechanistic MODELING paper that reports no measured clearance.** A *fitted* constant wearing a citation. **Human IL-6 clearance appears to be unmeasured.** |

**Elranatamab has no clinical IL-6 value in existence.** It was the molecule the project had spent the most
effort on.

A **fitted fallback** was also found and deleted: on any failure of the mechanistic IL-6 calculation, the code
silently substituted a constant fitted to `570` — a fabricated anchor — and emitted it under the *mechanistic*
field name. **It is now a hard error. A fallback that "never crashes a run" converts a loud failure into a
quiet fabrication.**

---

## Deliverables

| | |
|---|---|
| **Model documentation** | 9,000+ lines: `00_START_HERE`, `MASTER_MODEL_DOCUMENTATION` (L0–L9), `TCE_LIFECYCLE_WALKTHROUGH` (the life of the molecule, Shah-Betts → production), 9 subsystem deep-docs. **Every equation traced to `file:line`; every parameter tagged `[MEASURED]/[DERIVED]/[FITTED]/[ASSUMED]/[UNSOURCED]`; every PMID verified real.** |
| **Provenance audit** | `PROVENANCE_AND_VALIDATION.md` — the forensic record above. |
| **Not-in-use register** | **8 of 20 modules in `engine/` do not run.** One of them (`cytokine_pbpk.py`) is exactly what a reader would *want* to be in the model, and is not. |
| **Conservation guard** | The physical-ceiling check. Three lines. Would have caught the 942× bug on the first run. |
| **Construct screen** | Built, sized, dry-run: 84 constructs, 38 concurrent, ~15 min. **Blocked on the unsourced costim folds — by design.** |

---

## The honest summary

**The architecture is sound and the emergent claims are real.** The IL-6 arm is **not validated**, and the two
reasons are now precisely characterised: a scale bug (fixed) and a missing organ (quantified, with a
pre-registered prediction).

**What this model does that a fitted CRS model cannot: it told us it was wrong, and then told us why.**
