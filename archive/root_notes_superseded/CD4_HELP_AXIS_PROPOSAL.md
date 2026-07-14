# Proposed Axis 4: CD4 Help Preservation
(costim engager counter-screen — closing the "help vs. lineage-shutdown" gap)

## Why this axis is needed
The three-axis contract scores each costim receptor on:
  - Axis 1 CD8 effector benefit (want UP)
  - Axis 2 CRS liability        (want DOWN)
  - Axis 3 suppression liability (want DOWN)

A receptor whose AGONISM lowers the CD4 Treg/IL-10 program and storm cytokines
(negative Axis 2 & 3) scores as the ideal "Option B (CD4-negating)" arm — the best
therapeutic window. BUT axes 2 & 3 measure only the HARMFUL CD4 programs. They are
blind to CD4 HELP, which the project treats as PROTECTIVE:
  "CD4 help is beneficial; it drives CD8 durability and persistence, and CD4 cells
   can be direct effectors."  (project proposal)

Consequence: a receptor whose agonism broadly shuts CD4 down would look like a perfect
Option B (storm + suppression both fall) while actually destroying the CD4 help that
gives the CD8 boost its durability. Axis 4 distinguishes:
  - TRUE sub-program negation : agonism damps Treg/IL-10/storm, help PRESERVED   -> keep
  - lineage-wide CD4 shutdown : agonism damps everything incl. help              -> downgrade

This is exactly the proposal's own thesis — "the enemy is not the CD4 lineage but a
specific CD4 sub-program" — made into a scoreable guard.

## Gene sets (overlap-audited against axes 1-3; see integrate script)
Help mediators notoriously overlap storm/effector. Audit result:

  help_tfh_differentiation  [CLEAN — no overlap with any axis] -- PRIMARY guard backbone
    BCL6, MAF, TOX2, CXCR5, TCF7
    = the CD4 helper/follicular-helper differentiation state. Cleanly separable from
      storm and effector, so a reliable readout of "is the helper program intact."

  help_cd8_durability       [use, but FLAGGED] -- most CD8-relevant, least clean
    CD40LG  (overlaps: is itself a PANEL target -> self-exclude when scoring CD40LG)
    IL21    (overlaps: already in Axis-2 CRS-extended -> shared mediator; read in the
             help direction ONLY when the clean Tfh backbone agrees)
    = canonical CD4->CD8 help: CD40LG licenses DCs for better CD8 priming; IL21
      sustains CD8 persistence and counters exhaustion.

  help_shared_polyfunctional [REPORTED, NOT used in guard] -- inseparable
    IFNG, IL2, TNF
    = help/polyfunctionality markers that overlap other axes: IFNG and TNF are both
      effector (Axis 1) AND storm (Axis 2); IL2 is storm (Axis 2). All three excluded
      from the guard — using them would be circular; reported for transparency only.

  direction: HIGHER = more help preserved (want a costim arm to NOT tear down CD4 help).
  agonism-oriented, same convention as axes 2/3 (hero CRISPRi: agonism = -1 * KD effect).

## How it enters the nomination (guard, not a 4th maximization term)
Help is a PRESERVATION CONSTRAINT, not a target to maximize, so it does NOT go into the
linear window score W = zE - wC*zC - wS*zS (that would double-count with the guard, and
CD8-durability from help belongs in the QSP efficacy term as a persistence multiplier,
not the screen-level window). Instead H (agonism-oriented mean help z, primary = Tfh
backbone) refines the classification of effector-hit receptors:

  H_thr = tau (default 1.0), evaluated primarily on help_tfh_differentiation
  Among effector hits with a negating signal (C < -tau OR S < -tau):
    Option B+ (help-sparing negator)   : H >= 0        [ideal — dampens tox, spares/boosts help]
    Option B  (negator, help preserved): -tau <= H < 0 [acceptable — mild help dip]
    lineage-wide CD4 shutdown (DOWNGRADE): H < -tau     [tears down help too — NOT a clean win]
  Among selective (|C|,|S| <= tau):
    Option A  (CD8-selective)          : H >= -tau      [neutral on tox, help intact]
    help-eroding selective (CAUTION)   : H < -tau       [tox-neutral but quietly erodes help]

Every non-downgraded Option-B call previously stamped help_preservation=UNVERIFIED can now
be RESOLVED to VERIFIED / DOWNGRADED once Axis 4 is computed from the hero DE matrix.

## Data source
Same instrument as axes 2 & 3: the hero CD4 CRISPRi DE matrix (GWCD4i.DE_stats.h5ad).
All help genes are CD4-expressed and in the DE var set (to be confirmed on load; BCL6,
MAF, TCF7, CXCR5, CD40LG, IL21 are expressed in stimulated CD4 T cells). No new download.
NB: LTBR remains effector-only (not CD4-expressed) — Axis 4, like axes 2/3, cannot score it.
