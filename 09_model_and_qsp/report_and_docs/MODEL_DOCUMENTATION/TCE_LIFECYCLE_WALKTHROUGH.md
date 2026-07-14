# The Life of a T-Cell Engager — a mechanistic walkthrough

> This document follows **one molecule** of a T-cell engager from the syringe to the cytokine it causes,
> and explains at every stage **what is happening physically, why the model represents it that way, and what
> was rejected.** It is written to be read straight through.
>
> **Every stage names its subsystem doc (`T1`–`T9`) and the live code file.** Nothing here describes code
> that does not run — see `NOT_IN_USE_REGISTER.md`.
>
> **Honesty note, up front:** where a quantity is *fitted* or *assumed* rather than measured, this document
> says so in the text, not in a footnote. The IL-6 arm in particular rests on one large assumption
> (§9.4), and the reader is entitled to know that before the narrative seduces them.

---

## 0. Where this model came from — the build history

The model was not designed top-down. It was grown, and the order matters because each layer is only
trustworthy to the extent the one beneath it was validated first.

**Stage 1 — Reproduce Shah & Betts.**
The foundation is the platform PBPK model of Shah & Betts: a human body as a set of organs, each with a
blood flow `Q`, a lymph flow `L`, a vascular volume `Vv`, an interstitial volume `Vis`, and vascular/lymphatic
reflection coefficients `σ_V`, `σ_L` that determine how readily a 150 kDa antibody leaks out of a capillary.
**Nothing was invented at this stage.** The goal was to reproduce a published, validated antibody
biodistribution model exactly, so that everything built on top inherits a physiology that someone else already
defended. → **`T1`**, `engine/qsp_costim_window_v2.py`

**Stage 2 — Make the tissue cellular.**
Shah-Betts treats an organ's interstitium as a well-mixed box. A T-cell engager does not act on a box: it
acts on the *encounter between two specific cells*. So each organ's interstitium was populated with
**individual cells** — T cells, target cells, regulatory T cells, monocytes/macrophages — each placed in space,
each carrying its **own receptor copy numbers derived from single-cell RNA-seq**. → **`T2`**, **`T3`**

**Stage 3 — Bind them properly.**
An engager is bivalent by construction: one arm to CD3, one to the tumour antigen. Its behaviour is dominated
by **avidity, geometry and prozone**, none of which survive an occupancy shortcut. So the binding is solved as
**literal Rhoden kinetics** — explicit association/dissociation of first and second bonds, per cell, every step.
→ **`T3`**, **`T4`**, **`T5`**

**Stage 4 — Kill, and count the cost.**
Productive trimers drive a kill hazard; Tregs damp it; engaged T cells emit cytokine. → **`T6`**, **`T7`**

**Stage 5 — Make CRS mechanistic.**
The original cytokine readout was `IL6 = engaged_dwell_rate × IL6_SCALE`, where `IL6_SCALE` was **one constant
fitted to a single molecule**. That is not a model of cytokine release; it is a units conversion with a
tuning knob. It was replaced by per-cell, contact-gated myeloid emitters. → **`T8`**
*(That fitted path has since been deleted outright — `PROVENANCE_AND_VALIDATION.md` §3.)*

**Stage 6 — Counterscreen.**
With efficacy and toxicity both emerging from the same spatial mechanism, the model can be asked the question
it exists for: **for a construct nobody has built, what is the therapeutic window?**

---

## 1. The dose — arrival in plasma

A dose enters either **intravenously** (instantly into the plasma compartment) or **subcutaneously**, in which
case it appears through a first-order absorption process with bioavailability `F_sc` and rate `ka`. Route is
taken from each molecule's clinical label — mosunetuzumab and glofitamab are IV; the BCMA and GPRC5D engagers
are SC.

**Why this matters more than it looks:** SC absorption *slows and flattens* the plasma peak. Since cytokine
release is driven by how fast T cells become engaged, **route is a CRS lever, not a bookkeeping detail.**
The model gets this for free by representing absorption rather than assuming a plasma profile.

**A trap that has bitten this project:** the clinical regimens are **step-up ladders**, e.g. mosunetuzumab
`1 mg (d1) → 2 mg (d8) → 60 mg (d15) → 60 mg (d22)`. A simulation shorter than 15 days gives mosunetuzumab
**only its 1 mg priming dose** while giving teclistamab (whose 120 mg lands on d7) its **full** dose. Comparing
the two then compares different rungs of different ladders. **Use `TSIM_DAYS=24`.** → `T9`

---

## 2. Distribution — how the antibody reaches a cell

Each organ has a **vascular** and an **interstitial** compartment.

The vascular pool is treated as **quasi-steady**. This is not laziness — it is forced. Blood flow `Q` exceeds
the extravasation permeability-surface product `PS` by orders of magnitude, so `Q/Vv ~ 10⁵ /day`. Integrating
that explicitly alongside the slow tissue dynamics would be **catastrophically stiff**: you would either take
absurdly small timesteps or watch the solution explode. Solving the vascular pool at its algebraic steady state
each step is both faster **and** more accurate. → **`T2`**, `engine/wholebody_percell.py`

Antibody crosses into the interstitium by **convective extravasation**, resisted by the vascular reflection
coefficient `σ_V` (≈0.95 for a 150 kDa IgG — capillaries hold most of it back), and leaves via **lymph**. What
remains in plasma is protected from catabolism by **FcRn recycling**, which is why an IgG lives for weeks.

**The consequence that drives everything downstream:** only a *small fraction* of the dose ever reaches the
interstitium where the cells actually are. The concentration a T cell sees is **not** the plasma concentration.

---

## 3. The cells it finds — receptor copies from real tissue

Every cell in the model carries its own **surface receptor copy number**, converted from single-cell RNA-seq
rather than assumed. A T cell has *its* CD3 count; a target cell has *its* TAA count; and — critically — those
counts are **distributions, not means.** Two T cells in the same organ do not see the same drug in the same way.

This is the reason the model can rank constructs it has never seen: a construct's performance depends on the
*shape* of the receptor distribution it must work against, and only a per-cell model has that shape.
→ **`T3`** (data provenance: `rna_to_receptor.py`, **offline** — see `NOT_IN_USE_REGISTER.md` §2b)

---

## 4. The first bond

One arm of the engager finds a receptor. This is ordinary bimolecular kinetics — `kon`, `koff`, a `KD` — and
where the affinity is *measured by SPR* it is used as measured. **Where it is not, the model says so.** (Several
CD3 affinities are class-estimates rather than molecule-specific measurements; each is tagged in `T6`.)

---

## 5. The second bond — where a bispecific stops being two antibodies

Once the first arm is bound, the molecule is **tethered**. The second arm is no longer searching a 3D volume —
it is searching a **small tethered shell**, at an *effective concentration* far above the bulk. This is
**avidity**, and it is the entire point of the format.

The model computes that effective concentration from **geometry**: the arm's reach, the antigen's height above
the membrane, and the **span** of the linker between the arms. Two effects fight each other:

- **Dilution:** effective 2nd-arm concentration falls as **1/span** — a longer linker searches a bigger shell,
  so it is *less* concentrated where the partner actually is.
- **Feasibility:** if the span is *shorter* than the intermembrane cleft, the second arm **physically cannot
  reach** across. Below ~0.6× the cleft distance, feasibility is zero.

Their product has a **maximum where span equals the cleft distance (≈13 nm)**.

> **This is a real, falsifiable prediction, and it is counter-intuitive: past the cleft distance, LONGER
> LINKERS MAKE BRIDGING WORSE.** The instinct that "more flexibility = better reach" is wrong, because reach
> was never the binding constraint — *concentration* was. → **`T4`**

---

## 6. The synapse — and why prozone is not an equation

A T cell bound to an engager bound to a target cell is a **trimer**, and killing scales with trimers.

At high drug concentration, something famous happens: **both arms saturate separately.** Every CD3 is occupied
by an engager whose other arm is dangling; every TAA likewise. There is nobody left to bridge *to*. Efficacy
**falls** with increasing dose. This is the **prozone (hook) effect**, and it is the central dose-selection
problem for every TCE ever developed.

**The model does not contain a prozone equation.** It falls out of solving the binding chemistry honestly:
the trimer term is quadratic in free arms, so once arms are saturated, it collapses on its own. Prozone is
**emergent**. That is precisely what makes the model trustworthy for constructs that don't exist yet — an
imposed hook function would have to be re-fitted for every new format; an emergent one does not.
→ **`T5`**, `engine/kinetic_synapse.py`

---

## 7. The kill

Each engaged T cell accumulates a **kill hazard** proportional to its own trimer count, damped by local
regulatory T cells. Target cells die stochastically.

**Be honest about this layer:** the kill-rate constants (`k_death`, `k_hit`) are **`[FITTED]`** — calibrated
against the epcoritamab B-cell depletion time-course (`calib_kdeath.py`, offline). They are the model's main
efficacy tuning parameters, and **anything downstream of the kill law inherits that tag.** → **`T6`**

---

## 8. Costimulation — and the bug that was the whole point

A costim arm (4-1BB, OX40, ICOS, CD28, CD2, CD27, GITR) is added to make the T cell fight *harder* and *longer*.
Each T cell integrates its own costim receptor occupancy into signalling programs — **effector** (raises kill),
**exhaustion** (decays it), **suppression** (on Tregs).

**Here is the thing that makes costim arms interesting, and that the model originally got exactly backwards.**

4-1BB, OX40, ICOS and GITR are **activation-induced**. They are essentially **absent on resting T cells** and
appear only *after* TCR engagement. **That is the entire reason the field targets them** — the costim fires only
on T cells that have *already engaged tumour*, making it **tumour-conditional**. It is a safety feature written
into the biology.

The model read costim receptor density **once, from resting tissue, and never updated it.** So it saw 4-1BB at
≈0 and concluded 4-1BB was a poor arm — **penalising exactly the property that makes it a good one.** A screen
run on that model would have confidently ranked CD2 above 4-1BB, and any immunologist would have rejected the
result on sight.

It is now:

    da_i/dt = k_on · p_eng_i · (1 − a_i)  −  k_off · a_i        (per-T-cell activation memory)
    R_i(t)  = R_rest,i · (1 + (FOLD − 1) · a_i)                 (induced receptor density)

where `p_eng_i` is **that cell's own** engaged-synapse fraction. A T cell that never engages never upregulates.
**The tumour-conditionality is emergent, not imposed.**

> ⚠️ **OPEN BLOCKER, stated plainly:** the `FOLD` values are **not yet sourced.** The code **fails closed** —
> inducible arms **raise an error** rather than run at resting density, because running them at resting density
> is the very bug this fixes. **The costim screen cannot produce a valid ranking until those folds are
> sourced.** → **`T7`**
>
> *(Separately: `costim_arm` was previously never passed into the PD engine at all, so every costim construct
> silently ran as a plain TCE. Fixed. Any costim result predating that fix is void.)*

---

## 9. Cytokine release — the most scrutinised part of this model

### 9.1 Why the old approach was indefensible

`IL6 = engaged_dwell_rate × IL6_SCALE`, with `IL6_SCALE` fitted to one molecule. It cannot generalise, it
cannot rank a novel construct, and it silently absorbs any upstream error into itself. **Deleted.**

### 9.2 What the literature actually says

CRS IL-6 is **not made by T cells.** It is **myeloid-derived** — monocytes and macrophages (Giavridis, **PMID
29808005**: *"macrophages … the main overall source of IL-6"*; Norelli, **PMID 29808007**: *"Human monocytes
were the major source of IL-1 and IL-6"*). And induction is **contact-dependent**: it requires CD40L on the
activated T cell to touch CD40 on the myeloid cell.

### 9.3 The mechanism as built

**CD40L–CD40 is membrane-bound. The cells must physically TOUCH.** So the interaction radius is a *contact*
distance — the sum of the two cell radii:

    R_contact = r_macrophage (10.6 µm) + r_Tcell (3.5 µm) = 14.1 µm

*(An earlier version used 30 µm. That is the T-cell:target **synapse reach** — a different quantity entirely,
and far too permissive for a membrane-bound ligand pair.)*

Each macrophage agent then integrates **its own** local contacts with **engaged** T cells, and secretes if it
is one of the ~**3.9%** of monocytes that are intrinsic IL-6 secretors (PMID 37533643 — measured under maximal
stimulation, so this is a **cell-intrinsic** property, *not* something the spatial model should re-derive).

**Everything interesting emerges from this.** Saturation, per-molecule differences, the effect of tumour burden
and of myeloid geography — none are fitted. A drug that engages T cells *where the macrophages are* (bone
marrow: ~33.6 million; spleen: ~290 thousand) causes more CRS than one that engages them elsewhere, and the
model knows this **because it knows where the cells are.**

### 9.4 From secretion to a blood test — and the one big assumption

Secretion is a **rate** (pg/hr). A clinician measures a **concentration** (pg/mL). Converting requires a volume
and a clearance:

    dC/dt = production / V − k_deg · C

- **`V = 11.65 L`** — the **ECF** (interstitium 8.55 L + plasma 3.10 L). IL-6 is **21 kDa**, is **made in the
  interstitium**, and has **no FcRn** to retain it in circulation. It is *not* confined to plasma. `[ASSUMED —
  physical]`
- **`k_deg = 0.20 /hr`** — **`[FITTED]`.**

> ### **Read this before quoting any absolute IL-6 number.**
> `k_deg` is cited in the code to PMID 31268236. **That paper is a modeling paper. It reports no measured IL-6
> clearance.** A search for a real one found none: human IV trials report no PK; the SC half-life is
> absorption-limited; the rat study gives half-lives but no clearance. **Human IL-6 systemic clearance appears
> to be unmeasured.**
>
> Because `C_peak = production / CL`, a production that is 10× too high and a clearance that is 10× too high are
> **indistinguishable**. **Absolute plasma IL-6 therefore cannot validate this mechanism.**
>
> **Ratios can.** Clearance is a property of IL-6, not of the drug, so it **cancels** from any between-molecule
> comparison. The IL-6 arm is validated — or refuted — on **one number**:
>
> **clinical mosunetuzumab / teclistamab = 152 / 21 = 7.24×** *(both population means, both sourced)*
>
> Every other clinical IL-6 "anchor" this project has used was contaminated — including one that was **a page
> number from an FDA table of contents.** See `PROVENANCE_AND_VALIDATION.md` §1 before you trust any of them.

---

## 10. What the model is for

Efficacy and toxicity now emerge from **the same spatial mechanism**: a T cell engaged next to a tumour cell
kills it; the same T cell engaged next to a macrophage causes CRS. **Therapeutic window is not two models bolted
together — it is one geometry, read two ways.**

That is what licenses the counterscreen. For a construct that has never been built — a new target, a new costim
arm, a different affinity, a different linker length — the model can ask: *where do its T cells engage, and who
is standing next to them?*

**Provided** the mechanism has been shown to reproduce the ratios it claims to predict. **That test is defined
in §9.4 and, as of this writing, has not yet returned a trustworthy result.** Until it does, this document
describes a well-built machine of unproven accuracy — and says so.
