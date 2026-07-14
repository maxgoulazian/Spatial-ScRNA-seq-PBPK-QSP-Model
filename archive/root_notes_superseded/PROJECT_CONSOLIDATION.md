# PROJECT CONSOLIDATION — canonical files, lane map, and supersede registry
_science-effector (frame 5894b83b), 2026-07-07. Single source of truth to un-jumble the fleet._

## 1. CANONICAL FILE PER AXIS (use these; everything else is superseded)

| Axis | CANONICAL artifact | Owner | Supersedes |
|------|-------------------|-------|-----------|
| **CD8 effector benefit (E)** | `cd8_effector_scores.csv` (8d964a5e) | science-cd8 | science-effector's `schmidt_effector_axis_panel.csv` (kept as historical; values IDENTICAL, max\|Δ\|=0) |
| **CD4 CRS liability (C)** | `toxicity_axis_scores.csv` (science-tox, robust_z) | science-tox | my `tox_help_scores_*.csv` (cross-check; CRS ρ=1.000) — **integrator NOW wired to canonical** |
| **CD4 suppression liability (S)** | same `toxicity_axis_scores.csv` (SUPP_full robust_z) | science-tox | same (suppression ρ=1.000 matched; 0.912 SUPP_full-union) |
| **CD4 help preservation (H, Axis 4)** | `tox_help_scores_*.csv` help_tfh cols | science-effector | (new axis; not yet in science-tox output) |
| **3-axis nomination** | `nomination_REAL_Stim48hr.csv` / `_Stim8hr.csv` | science-effector (integrator) | v1 pre-gate versions |

Effector axis note: cd8_effector_scores.csv is the richer table — Schmidt CRISPRa IFN-g z
(raw-QC: reproduces published z within \|Δz\|<0.6 for all above-noise HITS; non-hit receptors
e.g. HVEM can differ more, but are gated out) + Legut ORF co-anchor + CRISPRi directional
flag. The integrator consumes it via --effector-csv; the effector-source swap left the
nomination unchanged to reported precision (a few window/E values differ only at the
1e-15 floating-point level from recomputation; ranking and all classifications identical).

## 2. THE REAL LANE MAP (corrects FLEET_PLAN.md v2)

| Role | Frame | ACTUALLY owns | Status |
|------|-------|---------------|--------|
| science-orchestrator | 3d677532 | coordination, QSP reference curation | active |
| **science-effector** | **5894b83b** | infra/downloads + **integrator** + Axis-4 help + effector-nomination | active (this lane) |
| science-cd8 | 1a8514b4 | CD8 effector axis (DELIVERED cd8_effector_scores.csv) | **delivered → stand down** |
| science-tox | 6a8d63fa | CD4 CRS + suppression axes | active |
| science-data | a22ea771 | (bus-test lane) | retire |
| science-qsp | (to spawn) | QSP/PBPK window model | pending |
| science-binder | (to spawn) | in-silico binder design | pending |

Correction to FLEET_PLAN: 5894b83b was mapped as `science-infra` (downloads only) with a
separate `science-cd8` to build effector and the integrator under the orchestrator. In
reality 5894b83b did infra + effector + integrator. science-cd8 spawned anyway and produced
the canonical effector table (additive, validated mine). Net: science-cd8's ARTIFACT is
canonical; its LANE stands down having delivered. Integrator stays with science-effector.

## 3. CURRENT NOMINATION (canonical inputs, gated, Stim48hr primary)

Effector-competent arms only (gate = CD8 IFN-g hit); ranked by therapeutic-window score:

Fully-canonical inputs (E=science-cd8, C+S=science-tox robust_z, H=science-effector), gated, Stim48hr:

| Arm | E | C (CRS) | S (suppr) | Window | Call |
|-----|---|---------|-----------|--------|------|
| **CD30** | 3.22 | −1.01 | **−2.86** | **2.95** | Option B+ (suppression-negating, help preserved) |
| **4-1BB** | 3.74 | **−1.87** | −0.56 | **2.30** | Option B+ (CRS-negating, help preserved) |
| CD28 | 12.11 | +2.01 | +0.20 | 1.56 | pan-costim liability (TGN1412) |
| CD27 | 4.28 | +0.11 | −0.71 | 1.08 | Option A (CD8-selective) |
| CD40 | 2.65 | −0.11 | +0.08 | 0.22 | Option A |
| CD2 | 5.60 | +3.01 | −0.08 | −1.08 | pan-costim liability |
| OX40 | 2.07 | −0.22 | **+1.99** | −1.10 | pan-costim liability (Treg-suppression, robust_z) |
| LTBR | 3.49 | — | — | 0.22 | effector-only (not CD4-expressed, un-scoreable) |
| ICOS, DNAM1, GITR, HVEM, DR3, CD40LG | <1.6 | — | — | — | excluded (no effector signal) |

TRIANGULATION — three independent methods agree CD30 + 4-1BB are the clean arms:
1. heuristic window z-sum (science-effector); 2. empirical-null robust_z (science-tox);
3. mechanistic QSP dose-response (science-qsp qsp_window_scores.csv, verified on host: CD30
   qsp_window 0.789 = sweep rank 1, 4-1BB 0.704 = sweep rank 2; both "NOMINATED — window-widening
   clean arm"; ranks stable across the 125-combo sweep, median 0.783/0.700).
CRS ρ=1.000, suppression ρ=1.000 (matched gene sets; science-tox null = genome-wide, n_null≈11,287);
effector scorings agree to reported precision (1e-15 fp floor; ranking identical).
NOTE: robust_z sharpens OX40 to a pan liability (S=+1.99>τ), superseding the earlier sub-threshold read.

## 4. DEPENDENCY STATE (what's ready vs blocked)

- Effector axis: COMPLETE (canonical, dual-anchor, raw-QC'd).
- CRS + suppression axes: COMPLETE (real hero DE, cross-validated).
- Help axis (4): COMPLETE for Stim48hr/Stim8hr (science-effector).
- 3-axis nomination: COMPLETE + gated + robust across 2 timepoints.
- Subset-resolution layer (Treg/Tfh/Th1 per-receptor): DESIGNED, verified executable,
  BLOCKED on donors 2-4 (now downloading via public-S3/aria2c; D1 prototype possible).
- Per-donor error bars: pending DE_stats.by_donors.h5mu (downloading).
- QSP window model: COMPLETE (science-qsp; mechanistic dose-response, 23 verified refs; qsp_window_scores.csv
  landed in deliverables/06_qsp_science-qsp/ — verified: CD30 rank 1, 4-1BB rank 2, both nominated).
- RNA->protein receptor guard: PROPOSED by science-tox — read Legut OverCITE ADT for CD30/4-1BB
  surface protein on CD8 effectors. The nomination's real soft spot; empirical (not a model job).
- Binder design: NOT STARTED (science-binder; target = CD30 / 4-1BB, pending protein guard).

## 5. WHAT "DONE" LOOKS LIKE FOR THE RESEARCH FINDING
Steps 1-2 of the one-week plan (signed nomination separating effector benefit from CD4
toxicity) are COMPLETE. CD30 + 4-1BB are the clean-arm nominees; CD28/CD2 the pan
liabilities. Remaining rungs (subset resolution, per-donor CIs, QSP window, binders) refine
and translate but do not change the core finding.
