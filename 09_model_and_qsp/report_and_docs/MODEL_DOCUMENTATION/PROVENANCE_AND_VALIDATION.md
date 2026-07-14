# PROVENANCE & VALIDATION — what is measured, what is fitted, what is fabricated

> **Read this before you trust any number this model produces.**
>
> This document exists because on 2026-07-13 a full provenance audit of the IL-6 arm found that **a page
> number from an FDA table of contents had been circulating as a clinical IL-6 concentration**, that a
> constant labelled "literature-sourced" was somebody else's *fitted* model parameter, and that a fitted
> scale factor could silently substitute itself for a mechanistic calculation and report itself under the
> mechanistic field name. None of these were caught by review. All of them were caught by *tracing every
> number back to the sentence it came from.*
>
> That exercise is now permanent policy, and this file is its output.

---

## 0. The provenance vocabulary

Every parameter in this model carries exactly one of these tags. There is no "probably fine" tag.

| tag | meaning |
|---|---|
| `[MEASURED]` | A number reported from an experiment, with a PMID/DOI/label and a **verbatim quote** containing the number. |
| `[DERIVED]` | Computed from `[MEASURED]` quantities by a stated formula. The formula is given. |
| `[FITTED]` | Adjusted so the model reproduces an observation. **Not** evidence — the observation it was fitted to cannot then be used to validate it. |
| `[ASSUMED]` | Chosen on physical/mechanistic grounds with no direct measurement. The rationale is stated. |
| `[UNSOURCED — TBD]` | We do not know where this came from. **It is a liability until resolved.** |

**A `[FITTED]` parameter with a PMID next to it is the most dangerous object in a model.** It looks like
evidence and behaves like a free parameter. Section 2 is a case study.

---

## 1. THE CLINICAL IL-6 ANCHORS — a forensic record

The model is validated against clinical peak serum IL-6 for approved T-cell engagers. Before 2026-07-13 the
repository carried the following values. **Every one of them was wrong, and they were wrong in different
ways.** They are recorded here — rather than quietly deleted — because the *failure modes* generalise.

### 1.1 The values that were in the code, and what they actually were

| molecule | value carried | what it **actually** is |
|---|---|---|
| elranatamab | **191** | **A PAGE NUMBER.** A dot-leader page reference from the *Table of Figures* of FDA BLA 761345. It is not a concentration. It was being used as one. |
| elranatamab | **340**, **230** | **FABRICATED — verified twice, directly against the primary source.** `params_myeloma.json` claims *"IL-6 digitized (Nat Med 2023 Fig6): priming-only peak ~340 pg/mL, priming+premed peak ~230 pg/mL."* The paper (Lesokhin et al., MagnetisMM-3, *Nat Med* 2023 — **PMC10504075**) has **Figures 1–4 plus Extended Data 1–2: there is no Figure 6**; **IL-6 does not appear anywhere in the paper**; and it **reports no cytokine concentrations in any units.** These are digitizations of a figure that does not exist, from a paper containing no IL-6 data. |
| elranatamab | **366.88** | **No source exists.** (Surfaced from a web search during the audit and could not be traced to any primary document.) |
| mosunetuzumab | **570** | **No source exists.** The real reported population mean is **152**. |
| teclistamab | **288** | **Real** — but it is the *highest individual patient* Cmax. See §1.3. |
| teclistamab | **21.5** | **Real** (the value is 21) — the *population mean*. See §1.2. |
| talquetamab | **18.2** | **Real** — but a *median*, not a mean. Not comparable to the mean pair. |

### 1.2 The statistic-mixing failure (the one that did the most damage)

From MajesTEC-1 (**PMID 38831634**), in a single sentence:

> *"The **mean** IL-6 peak concentration (Cmax) was **21 pg/mL**, while the **highest IL-6 Cmax observed
> among patients** receiving the recommended phase II dose was **288 pg/mL**."*

**21 and 288 are both true, differ by 13×, and describe different things.** The repository held teclistamab's
**mean** (21.5) alongside what were taken to be *peaks* for the other molecules. Ranking a population mean
against an individual peak manufactured an apparent **~26× clinical dynamic range** that no mechanism could
reproduce — and the mechanism was then repeatedly "corrected" to chase it. Secretor fraction, contact radius,
distribution volume, and clearance were all adjusted, at various points, against a target that was an artifact
of mixing two statistics.

**Rule, now enforced in code:** every anchor in a comparison set must be the **same statistic** at a
**comparable dose stage**. The validation harness *refuses to score* a mixed-statistic set rather than
silently producing a ratio.

### 1.3 Why "highest individual patient" cannot be a validation target

The maximum of a sample is an **order statistic**: its expectation grows with the sample size N. A 165-patient
trial will report a higher "highest patient" than a 40-patient trial *from sampling alone*, with no difference
in the underlying biology. Cross-trial comparison of maxima therefore measures cohort size as much as drug
effect. Such values are usable as a **worst-case envelope**; they are **not** a central-tendency target.

### 1.4 What survives

**Only two clinical IL-6 anchors survived the audit** — same statistic, same kind of measurement, both
verbatim-sourced:

| molecule | value | statistic | source |
|---|---|---|---|
| teclistamab | **21** pg/mL | population MEAN IL-6 peak Cmax | PMID 38831634 (verbatim quote above) |
| mosunetuzumab | **152** pg/mL | population MEAN IL-6 peak Cmax | — |

> **CLINICAL RATIO mosunetuzumab / teclistamab = 7.24×.** This is the single quantitative target the IL-6
> arm is validated against.

**ELRANATAMAB HAS NO CLINICAL IL-6 VALUE IN EXISTENCE** — not in any paper, FDA review, or label that the
audit could locate. It **cannot** be used to validate the IL-6 arm, and has been removed from the anchor set
(`il6_obs=None`). This is the molecule the project had spent the most effort on.

---

## 2. THE IL-6 CLEARANCE CONSTANT — a fitted parameter wearing a citation

`engine/myeloid_il6.py`:

```python
KDEG_IL6_PER_HR = 0.20    # cited in-code to PMID 31268236
V_PLASMA_ML     = 11650.0 # ECF: interstitial 8.55 L + plasma 3.10 L
                          # -> implied CL = 0.20/hr x 11.65 L = 55.9 L/day
```

**PMID 31268236 is Chen et al. 2019, *Clinical and Translational Science*, "A Modeling Framework to
Characterize Cytokine Release upon T-Cell-Engaging Bispecific Antibody Treatment."** It is a *semi-mechanistic
PK/PD **modeling** paper*. **It reports no measured IL-6 clearance.** The 0.20/hr is a value *fitted inside
somebody else's model*, which this model then adopted as though it were an experimental measurement.

Searching for a real one:

* Human **IV** IL-6 phase-I trials report toxicity and dosing but **no PK parameters**.
* Human **SC** IL-6 gives t½ ≈ 4.2–5 h — but SC kinetics are **absorption-limited** (flip-flop), so this is
  *not* an elimination half-life.
* The rat study (**PMID 3263918**) gives biphasic half-lives (3 min / ~55 min) and shows clearance is ~80%
  **hepatic, receptor-mediated** (*"about 80% of the ¹²⁵I-rhIL-6 that had disappeared from the circulation
  was found in the liver"*) — but reports **no volumetric clearance**.

> ### **Human IL-6 systemic clearance appears to be UNMEASURED.**
> Consequently `KDEG_IL6_PER_HR` is tagged **`[FITTED]`**, never `[MEASURED]`, and it is **the single largest
> assumption in the IL-6 arm.**

### 2.1 The identifiability consequence — why this matters more than it looks

Plasma concentration is set by
`C_peak = production_peak / CL`.

With **one** molecule this is degenerate: a production that is 10× too high and a clearance that is 10× too
high are **indistinguishable**. Choosing CL to make the model land on a clinical value **is fitting** — it
merely relocates the free parameter from the emitter to the eliminator, where it is harder to see.

**Therefore:**
1. **Absolute** plasma IL-6 cannot validate the mechanism on its own.
2. **Ratios between molecules CAN** — CL is a property of *IL-6*, not of the drug, so it **cancels**:
   `C_i / C_j = prod_i / prod_j`. This is why §1.4's ratio is the stated target.
3. If several molecules independently imply the *same* CL, that agreement is a genuine test with real degrees
   of freedom. **Agreement is the validation; reproducing the peaks is only the fit.**

---

## 3. THE SILENT-FALLBACK DEFECT (deleted 2026-07-13)

`run_tce_pd_reval.py` previously contained:

```python
il6        = np.array(r.get('il6_plasma_pgml') or [])          # mechanistic
il6_legacy = np.array(r['sys_cyto_rate']['IL6']) * IL6_SCALE   # retired FITTED path
if il6.size == 0: il6 = il6_legacy        # "safety: never crash a run"
```

`IL6_SCALE` was **one constant fitted so that mosunetuzumab → 570 pg/mL** — i.e. fitted to a value that
§1.1 shows **has no source at all.**

The fallback meant that if the mechanistic recorder ever failed, the run would **silently emit the fitted
constant under the mechanistic field name** (`il6_peak`), with no warning, no flag, and no difference in the
output schema. A downstream reader — human or script — could not tell.

**Status: DELETED.** An empty mechanistic array is now a **hard `RuntimeError`**. A missing IL-6 calculation
is a bug to fix, not a number to substitute.

> **The general lesson, which is the reason this section exists:** a fallback that "never crashes a run" is a
> fallback that **converts a loud failure into a quiet fabrication.** Crashing is a feature.

---

## 4. VALIDATION STATUS — stated honestly

| arm | status |
|---|---|
| **PK** (22 molecules vs clinical plasma) | Independent of the IL-6 arm. See `FIGURES_AND_VALIDATION.md`. |
| **IL-6 / CRS** | **NOT YET VALIDATED.** The decisive test — model vs clinical **7.24×** ratio (§1.4) — is defined and runnable but had not produced a trustworthy number as of this writing. Every earlier IL-6 "validation" in this project was scored against contaminated anchors and is **void**. |
| **Costim arms** | **BLOCKED.** Activation-induction fold-upregulation for 4-1BB / OX40 / ICOS / GITR is **unsourced**; the code **fails closed** (raises) rather than running them at resting density. See `subsystems/T7_*`. |

### 4.1 Runs that must NOT be trusted (they will look good and they are wrong)

* **Any run predating the myeloid census** (`organ_myeloid_counts.json`). Before it loaded, `count_scale`
  defaulted to `1.0` and IL-6 was **~290,000× too low**. Files named `tce_pd_rd_*` show teclistamab ≈ 25 pg/mL
  against a clinical mean of 21 — an apparently excellent result **produced by a bug.** If these are found
  later they will look like the best data in the project.
* **Any run at `TSIM_DAYS=7`.** The step-up regimens mean a 7-day window gives **mosunetuzumab only its 1 mg
  priming dose** (its 60 mg dose is on **day 15**) while giving **teclistamab its full 120 mg** (day 7).
  Molecules are then compared at different rungs of their own dose ladders. Use `TSIM_DAYS=24`.

---

## 5. STANDING RULES (violating these is how today happened)

1. **Trace every number to the sentence it came from.** Not to a comment. Not to a variable name. Not to a
   prior version of this file. To a *sentence in a source document containing that number.*
2. **Never mix statistics.** A mean is not a median is not a maximum. State the statistic beside every value.
3. **A citation next to a number does not make it measured.** Open the paper. Chen 2019 is a modeling paper.
4. **No silent fallbacks.** Fail loudly.
5. **A fit cannot validate itself.** If a constant was tuned to an observation, that observation is spent.
6. **Suspicion scales with convenience.** A weak source plus a flattering result is when to look hardest —
   both `570`→`618.6` and the pre-census `21`→`25` looked like triumphs and were artifacts.
