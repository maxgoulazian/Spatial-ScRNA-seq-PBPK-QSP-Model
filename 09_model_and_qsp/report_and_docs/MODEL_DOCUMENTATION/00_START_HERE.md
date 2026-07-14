# 00 — START HERE

**The TCE Spatially-Resolved Counterscreen Model**
*A whole-body, per-cell, mechanistic model of T-cell-engager efficacy and cytokine-release risk.*

---

## What this model is, in one paragraph

It simulates what happens to **every individual cell in a human body** when a T-cell engager is dosed. A
whole-body PBPK backbone (reproduced from Shah & Betts) carries the antibody to each organ; inside each organ,
**individual T cells, target cells, regulatory T cells and myeloid cells** are placed in space with receptor
copy numbers taken from single-cell RNA-seq. The engager binds them with **literal kinetic multivalent
chemistry** (not an occupancy shortcut), forms immune synapses, kills targets, and — where an engaged T cell
happens to be **physically touching** a monocyte — triggers IL-6 secretion from that individual myeloid cell.
Cytokine release is therefore **not a fitted function of drug concentration**: it is an emergent consequence of
where the cells are and who is touching whom. The purpose is a **counterscreen** — to rank candidate construct
designs (target × costimulatory arm × affinity × valency × linker length) by therapeutic window, *including
designs that have never been built.*

---

## If you are a new instance (Atlas / Codex / agent), read in this order

1. **`NOT_IN_USE_REGISTER.md`** — **first, always.** 8 of the 20 modules in `engine/` **do not run.** One of
   them (`cytokine_pbpk.py`) is exactly what you would *want* to be in the model and is not. Do not document,
   cite, or trust code before checking it against the live import graph.
2. **`PROVENANCE_AND_VALIDATION.md`** — what is measured, what is fitted, what was **fabricated**. This model
   has a documented history of a fitted constant carrying a PMID, and of **a page number circulating as a
   clinical IL-6 concentration**. Read it before you believe a number.
3. **`MASTER_MODEL_DOCUMENTATION.md`** — parameters, equations, modules, integration, dataflow.
4. **`TCE_LIFECYCLE_WALKTHROUGH.md`** — the life of the molecule, dose to cytokine, told as a narrative.
5. **`subsystems/T1…T9`** — the deep docs. One per live subsystem.

---

## The eight non-negotiable invariants

Violating any of these has, at some point, produced a confidently-wrong result in this project.

1. **Trace every number to the sentence it came from.** Not a code comment. Not a variable name. Not a prior
   doc. A *sentence in a source document containing that number.* (This rule exists because `191` — used as a
   clinical IL-6 concentration — turned out to be a **page number** from an FDA table of contents.)

2. **Never mix statistics.** A mean is not a median is not a maximum. A population mean ranked against an
   individual patient's peak produces a ratio that is **pure artifact**. State the statistic beside every
   clinical value, and refuse to compare across statistics.

3. **A citation next to a number does not make it measured.** Open the paper. The IL-6 clearance constant is
   cited to PMID 31268236 — a *modeling* paper that reports no such measurement.

4. **No silent fallbacks. Crashing is a feature.** A fallback that "never crashes a run" converts a loud
   failure into a **quiet fabrication.** The IL-6 path once silently substituted a fitted constant under the
   mechanistic field name. It is now a hard error.

5. **A fit cannot validate itself.** If a constant was tuned to reproduce an observation, that observation is
   spent and can no longer serve as evidence. Absolute IL-6 is clearance-limited and clearance is fitted;
   therefore **ratios between molecules**, from which clearance cancels, are the only honest IL-6 test.

6. **Suspicion scales with convenience.** A weak source plus a flattering result is precisely when to look
   hardest. Two of the best-looking results in this project's history (`570 → 618.6` and a pre-census
   `21 → 25`) were artifacts.

7. **Verify by name, never by exit status.** Parallel runs have repeatedly been launched twice for the same
   molecule, or silently OOM-killed. Confirm what is running with
   `ps aux | grep '[r]un_tce_pd_reval' | awk '{print $NF}' | sort | uniq -c`.

8. **Run the gate; do not reason about it.** In a coupled multi-cell model, "this flag obviously cannot
   affect that" has been wrong more often than right. Measure it.

---

## The subsystems (one line each → its deep doc)

| | subsystem | what it does |
|---|---|---|
| **T1** | Shah-Betts PBPK backbone | organ blood/lymph flows, volumes, reflection coefficients — the physiology the model was *built by reproducing* |
| **T2** | Whole-body per-cell PK | vascular (quasi-steady) → interstitium → lymph → plasma; FcRn; SC absorption |
| **T3** | Rhoden bivalent binding core | literal kinetic bivalent binding, per cell; backward-Euler solver |
| **T4** | Multi-arm format geometry | cis/trans co-engagement, reach kernel, **linker span optimum** |
| **T5** | Kinetic immune synapse | engage / hit / detach; trimer; **emergent prozone** |
| **T6** | Per-cell PD — killing | ternary equilibrium, kill hazard, Treg damping |
| **T7** | Costim signaling & induction | per-cell costim programs; **activation-induced** receptor density |
| **T8** | Mechanistic CRS IL-6 | contact-gated per-cell myeloid emitters → plasma ODE |
| **T9** | Integration & driver | how a run is assembled and stepped end-to-end |

---

## Glossary

| term | meaning |
|---|---|
| **TCE** | T-cell engager — a bispecific that bridges CD3 on a T cell to a tumour antigen (TAA) on a target cell. |
| **trimer** | The productive ternary complex T-cell·CD3 — engager — TAA·target. Killing scales with it. |
| **prozone / hook** | At high drug, both arms saturate *separately* and bridging **falls**. In this model it is **emergent** from the binding chemistry, not imposed by an equation. |
| **step-up dosing** | Small priming doses before the full dose, to blunt CRS. Clinically, **CRS is worst at the FIRST dose** — this is the single most important fact a TCE model must reproduce. |
| **CRS** | Cytokine release syndrome. Here: myeloid-derived, **contact-gated** IL-6. |
| **counterscreen** | Ranking construct designs by efficacy *vs* toxicity — the model's actual purpose. |
| **count_scale** | Sampled cells → physiological cell count. **Myeloid count_scale is a tissue property (drug-independent); target-cell count_scale is drug-dependent.** Confusing the two corrupts the counterscreen. |
| **ECF** | Extracellular fluid = interstitium (8.55 L) + plasma (3.10 L) = 11.65 L. IL-6's distribution space — it is made in tissue, is 21 kDa, and has no FcRn to hold it in blood. |

---

## Maintenance

- Any mechanistic change to the engine **must** update the corresponding `subsystems/TX_*.md` in the same
  change. A doc that has drifted from the code is worse than no doc — it is a confident lie.
- Re-derive the live/dead module split from the import graph (`NOT_IN_USE_REGISTER.md` §3) **before** editing
  docs. Modules get orphaned silently.
- When a parameter's provenance changes, update `PROVENANCE_AND_VALIDATION.md`, not just the code comment.
  The code comment is where the last fabrication hid.
