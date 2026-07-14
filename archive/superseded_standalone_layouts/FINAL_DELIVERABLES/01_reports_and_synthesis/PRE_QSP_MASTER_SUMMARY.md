# PRE-QSP MASTER SUMMARY — complete validated state upstream of the QSP
_generated 2026-07-12T03:37:07Z · single source of truth for everything feeding the QSP · v7-consistent_

This digests the entire analysis pipeline (A-track axes, B-track TAA, S/C verification)
as it stands the moment before QSP dose-characterization. Nomination is FINAL and validated.

================================================================================
## 1. THE NOMINATION (final)
================================================================================
**Costim co-leads: 4-1BB (TNFRSF9) + CD27.** Derived FOUR independent ways, all concordant:
  (i)   3-axis hand-curated score (effector benefit vs SUPP + CRS liability)
  (ii)  genome-wide effector-first scan (panel-free, all 18,875 genes) — A17
  (iii) network effector selectivity (signed-edge sum onto cytotoxic targets) — A29 (4-1BB +54.4)
  (iv)  genome-scale hero-CD4 GRN re-derivation through the QSP — A31b (CLEAN set stable)

Full construct (trispecific):
  - Redirector / Signal-1 : CD3 (VH/VL, OKT3 epitope)
  - Costim (novel core)   : 4-1BB + CD27
  - TAA (tumor arm)       : CEACAM5 lead + 5 finalists (below)

================================================================================
## 2. THE 6-AXIS LIABILITY GATE (drive-independent veto)
================================================================================
An arm with ANY liability-up axis is eliminated; effector NEVER offsets a liability.
Axes: CRS / SUPP / HELP-erosion / PROLIF / EXH  (+DD_SUPP = no independent gate call).

  CLEAN  : 4-1BB, CD27
  GATED  : CD28[CRS,SUPP,PROLIF]  ICOS[HELP,PROLIF]  DNAM1[SUPP,EXH]  OX40[SUPP,EXH]
           GITR[SUPP]  HVEM[SUPP,EXH]  DR3[SUPP,EXH]  CD30[HELP,PROLIF]  CD40[HELP,PROLIF]  CD2[CRS]

  Signature result: CD28 has the TOP effector score (z=11.94) and is STILL gated out (CRS,SUPP,PROLIF)
  — the clearest demonstration that the counter-screen, not effector, drives the nomination.

================================================================================
## 3. AXIS INVENTORY (all current, all incorporate hero dataset)
================================================================================
EFFECTOR (benefit):
  A0eff  genome-wide effector pool (Schmidt CRISPRa CD8 IFN-γ, 18,875 genes)
  A0diff CD8-vs-CD4 effector differential — 4-1BB +3.81, CD27 +3.40 (top 2 of 12)
  A17    genome-wide effector-first novel-target scan (no novel arm beats co-leads on mechanism)
  A29    network effector selectivity (4-1BB +54.4; CD8 57.4 vs CD4 3.0, ~19x)
LIABILITY (6 counter-screen axes, hero CD4 Perturb-seq, genome-scale):
  CRS (TNF/IL2/IFNG) · SUPP (15-gene Treg/IL-10) · HELP · PROLIF (cell-cycle) · EXH (TOX/PDCD1/…) · DD_SUPP
  A24    kinetic CRS/SUPP onset, genome-wide (Stim8hr vs Stim48hr)
  A25    CRS master-regulator (storm-suppressor) genome-wide scan
EXPRESSION / SELECTIVITY:
  A21    per-compartment expression breadth (CD8 4-1BB 57.5% / CD27 84.4%)
  A34    CD8-vs-CD4 expression selectivity (CD27 +0.43, 4-1BB +0.07)
  A22/A32 redirector×costim single-cell co-expression
NETWORK (GRN):
  hero CD4 genome-scale GRN (1.65M TF→target edges) + A31b QSP DE-vs-GRN re-run
  Tier-2 costim GRN → honest scope note (A26) + genome-scale hero GRN delivered

================================================================================
## 4. A31b — GRN vs DE QSP RE-RUN (the network confirmation)
================================================================================
Nomination HOLDS under BOTH drives. CLEAN = {4-1BB, CD27} either way.
  4-1BB : window +1.57 (DE) / +1.57 (GRN) · TI 62.4 → 744.4 · cap none   [widest, CRS-coldest]
  CD27  : window -2.37 (DE) → +1.27 (GRN) · CRS-capped both              [network reads CD27 safer]
  POST-GATE reading mandatory: raw GRN window puts CD30 (0.75) & CD28 (0.72) above co-leads,
  but both are hard-gated — reading raw resurrects the "CD28 looks clean-ish" artifact.

================================================================================
## 5. TAA FINALISTS (B-track, Lee et al. 2020 CRC: SMC GSE132465 + KUL3 GSE144735)
================================================================================
  CEACAM6   balanced# 1  z_restr +2.37  repl=True   CLINICAL: dual CEACAM5/6 ADC EBC-129 (FDA Fast Track 2025)
  CEACAM5   balanced# 2  z_restr +2.39  repl=True   CLINICAL: TCE (cibisatamab) + ADCs (M9140, tusa rav, SGN-CEACAM5C)
  ITGB4     balanced# 3  z_restr +1.39  repl=True   PRECLINICAL: anti-CD3/ITGB4 bispecific-armed T cells (mouse)
  TSPAN6    balanced#11  z_restr +1.05  repl=False  UNTARGETED (novel)
  LY6E      balanced# 8  z_restr -0.20  repl=True   CLINICAL: ADC DLYE5953A Ph1 (no TCE — whitespace)
  DPEP1     balanced#12  z_restr +1.97  repl=True   EARLY/PRECLINICAL (CRC-associated; GPI-anchored)

  Lead pair CEACAM5+CEACAM6 (highest restriction, replicating, clinically precedented).
  3/6 clinically targeted (CEACAM5 TCE+ADC, CEACAM6 ADC, LY6E ADC); ITGB4/DPEP1 preclinical; TSPAN6 novel.
  Binder design DEFERRED (stretch goal) except CEACAM5 (queued as de novo positive control).

================================================================================
## 6. VERIFICATION STATUS
================================================================================
  V7_VERIFICATION_LEDGER.csv : 30 claims, ALL PASS (incl. A31b, DD_SUPP no-independent-gate, Lee-cohort fix)
  S3_verdict_master_v7.csv   : 25 hypotheses — 22 supported + 1 reconciled / 2 inconclusive / 0 refuted
  C3_v7_master_verification  : 10 primary-source re-derivations, all PASS
  Two auditor findings this session fixed: (1) nomination rule = liability gate (NOT gate∩positive-window);
  (2) ledger count = 30 (current), earlier 16 was a prior-session version.

================================================================================
## 7. WHAT REMAINS (downstream of this summary)
================================================================================
  - QSP TAA-density/selectivity term fold-in (CPU, model input — not an analysis gap)
  - Binder campaigns running (4-1BB→CD27→CD3→CEA5→ProteinMPNN); AF3/Boltz-2 iPTM screen (cloud)
  - Deck / consolidated construct figure (30% rubric — least built)
  - Spatial PBPK-QSP (9-organ, awaiting user's additional Xenium tissues → 11+Tumor)
