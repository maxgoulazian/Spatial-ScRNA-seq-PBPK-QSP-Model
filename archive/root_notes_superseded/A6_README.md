# A6 — Method validation vs clinical ground-truth  (CONFIRMATORY)

Generated: 2026-07-09T01:59:41Z · lane latent-final · env `tools/sc-analysis-venv` (pandas 3.0.3, py 3.12.3)
Build: `A6_build_method_validation.py` · Verify: `A6_verify.py` (5/5 PASS) · Wave-2 (depends on A1 + A2)

## WHAT
Cross every costim-panel receptor's **computational method-call** (our A1 effector + A2 CD4 tox wiring)
against its **independent clinical ground-truth** (approved / clinical-agonist / failed-trial / preclinical /
novel), and emit a data-driven `validation_outcome`. 14 receptors scored. This is the external-validity check
on the whole pipeline: if our method is calling costim biology correctly, it must — with no clinical input to
the scoring — independently recover the known winners and losers. It does:

| receptor | clinical ground-truth | our method_call | validation_outcome |
|---|---|---|---|
| 4-1BB (TNFRSF9) | validated_good (approved-class) | FAVORABLE (effector+, no CD4 liability) | **TRUE-POS** |
| CD27 | validated_clinical (varlilumab) | CLEAN LEAD (effector+, suppression-negating) | CONCORDANT |
| CD30 (TNFRSF8) | novel_whitespace (ADC-approved, costim unexplored) | CLEAN LEAD (suppression-negating) | **NOVEL / white-space** |
| CD40 | clinical_apc_side | APC-side / not cis-costim | CONCORDANT (APC-side) |
| CD28 | cautionary_bad (TGN1412) | PAN-LIABILITY (effector+ but CRS-driving) | **TRUE-NEG** |
| OX40 (TNFRSF4) | clinical_underperformer | Treg/SUPPRESSION liability | **EXPLAINS UNDERPERFORMER** |
| GITR (TNFRSF18) | clinical_underperformer | Treg/SUPPRESSION liability | **EXPLAINS UNDERPERFORMER** |
| ICOS | clinical_underperformer | HELP-ERODING / no effector benefit | **EXPLAINS UNDERPERFORMER** |
| DR3 (TNFRSF25) | preclinical_treg | Treg/SUPPRESSION liability | CONCORDANT (preclinical Treg) |
| DNAM-1 (CD226), HVEM, CD2, CD40L, LTBR | not_clean_arm | (various) | CONCORDANT (not a clean arm) |

## WHY
Nomination-decision served: this is the **method-validation** confirmatory anchor (spec L191). It licenses the
three-axis nomination by showing the axes independently reproduce clinical reality — the spec-named recoveries
(CD28=TGN1412 bad, 4-1BB approved-class good, OX40/GITR/ICOS Treg-associated underperformers) all land, and the
one deviation is informative not error: **CD30 is a favorable method-call with no clinical costim-agonist
program = genuine white-space**, the project's lead novel nomination. Evidence tier: **CONFIRMATORY** — the
method_call is derived from confirmatory A1/A2 scores; the clinical layer is external ground-truth, not modeled.

## DATA USED
- `deliverables/07_final_prequsp/A1_effector_axis_final.csv` — A1 effector axis (Schmidt CRISPRa CD8 IFN-γ z, Legut ORF, Shifrut SLICE; per-receptor effector_hit + z). Supplies effector benefit.
- `analysis/rerun_2026-07/04_nomination/A2_nomination_Stim48hr_3axis.csv` — A2 tox nomination, **core-11**, Stim48hr, 4-donor within-donor Stouffer (CRS/SUPP/HELP calls + BH-q on the core-11 family). Supplies CD4 tox wiring.
- `analysis/rerun_2026-07/02_scored/deepen_2donor_D1D2D3D4_Stim48hr.csv` — A2 scored raw (19 receptors); used only for the two **exploratory** rows (CD2, CD40LG) not in the core-11 BH family.
- `deliverables/03_toxicity_axes_science-tox/costim_clinical_status.csv` — clinical ground-truth, each row naming the drug/trial/status (Utomilumab/urelumab, varlilumab CDX-1127, TGN1412, MEDI0562, TRX518, feladilimab GSK3359609, brentuximab-vedotin, …).

## METHOD
- **method_polarity** = documented function of A1 `schmidt_E_hit` + A2 `{CRS,SUPP,HELP}_call` + `mechanism_gate`.
  Precedence: no-CD4-tox → APC-gate → CRS-up-liability → SUPP/Treg-up-liability → help-erosion(+no effector) →
  low-effector → favorable(clean if suppression-negating). **Honesty rule:** every label reflects the axis that
  actually drives it in our data — ICOS is `help_erosion` (its A2 SUPP is *down*conc*, i.e. it LOWERS
  suppression) + negative effector z, NOT a SUPP-up Treg liability; OX40/GITR/DR3 are the true SUPP-up liabilities.
- **Axis-call rule** = `assemble_nominate.py:137-143` verbatim (QSIG=0.05; sign→up/down, +*conc if concordant).
  core-11 q consumed from A2 (BH family = core-11, reproduced this run: CD28 CRS q=1.1045e-05 @ n=11).
  CD2/CD40LG use their own Stouffer p (flagged: NOT in the core-11 BH family; exploratory).
- **clinical_polarity** mapped from `costim_dev_status` + sourced notes; **validation_outcome** = documented
  truth-table crossing method_polarity × clinical_polarity. No hardcoded per-receptor verdicts.
- **LTBR** = not CD4-perturbed in the hero screen → tox axes NaN + reason (`effector_only_no_CD4_tox`), never imputed.

## SOURCES
Clinical labels are **sourced, not asserted** (VERIFY gate G2 = all 14 pass). Per-row source in `clin_source`.
Spot-check citations (VERIFY gate G3). NOTE: web_search this session returned result **titles + URLs only** (no abstract text); only facts inferable from those titles/URLs are asserted below — trial statistics, NCT ids, and clinical timelines seen in an earlier draft were fabricated-from-memory and were removed in a 2026-07-09 audit:
- **CD28 / TGN1412** — anti-CD28 agonist mAb; cytokine storm in a phase-1 trial (Suntharalingam et al., *NEJM* 2006, **PMID 16908486**) — sourced from the web_search result title+URL this session. (Trial-level timeline/volunteer-count and the specific CD28-bispecific NCT programs were NOT in the retrieved text and are not asserted.)
- **OX40 / MEDI0562 (tavolimab)** — humanized OX40 agonist mAb studied in advanced solid tumors; safety/tolerability in *Clin Cancer Res* 2022;28(17):3709 (AACR) — sourced from the web_search result title+URL this session. (irAE %, ORR, NCT ids, GSK3174998, and the intratumoral-Treg readout were NOT in the retrieved text and are not asserted.) The clinical-underperformer call rests on the sourced 'CLINICAL (agonist), efficacy modest; Treg-associated' status in the tox-lane clinical-status file.
- **4-1BB** — utomilumab & urelumab agonists clinical; tumor-targeted 4-1BB costim bispecifics well into clinic (e.g. cinrebafusp/PRS-343); the canonical costim-engager arm. Source: tox-lane `costim_clinical_status.csv` note (no web_search was run for 4-1BB).
Dataset sources for the method side: Zhu 2025 hero CD4 Perturb-seq (via A2), Schmidt 2022 GSE174255 / Legut 2022 GSE193736 / Shifrut 2018 GSE119450 (via A1).

## OUTPUTS
- `A6_method_validation_vs_clinical.csv` — 14-row table (symbol, alias, clinical_dev_status, clinical_ground_truth, method_call, method_polarity, effector_z/_hit, CRS/SUPP/HELP calls, mechanism_gate, evidence, tox_scope, validation_outcome, in_core11_nomination, clin_source) + provenance header. → `.prov.json` sidecar.
- `A6_verify_verdict.json` — 5-gate independent verify (G1 recoveries, G2 all-sourced, G3 CD28+OX40 spot-check, G4 method-reproduces, G5 no-label-contradicts-axis) — **overall_pass = True**.
- `A6_build_method_validation.py`, `A6_verify.py` — the exact generating + verifying scripts.

## PROVENANCE / DATE
2026-07-09T01:59:41Z · lane latent-final. Numbers moved vs the stale hand-built `analysis/latent/METHOD_VALIDATION_vs_clinical.csv`
(which had **no generator script** — the anti-pattern this module retires): (1) that file's ICOS row was
"(not scored)"; here ICOS is scored and correctly reads help-erosion + no-effector. (2) CD30/CD27 method_calls
now trace to this-run A1/A2 numbers, not prose. (3) LTBR/CD40L/HVEM/DNAM1 given explicit `not_clean_arm`
concordance instead of asserted labels. Every value computed/read this run; (2026-07-09 fabrication audit: NCT ids + trial-stat details struck from script/CSV/README/ledger; scores unaffected.) CD28 CRS q + OX40 SUPP reproduced
from A2 to full precision by the verify.
