# VALIDATED NOW vs STILL UNDERWAY — honest submission-state note
2026-07-13 evening. Written to separate what is artifact-backed and reproducible right now
from what is in progress, so nothing is overclaimed at submission.

## VALIDATED NOW (artifact-backed, engine-faithful, reproducible)
- **Nomination** — 4-1BB(TNFRSF9) + CD27 co-leads; CD28 top-effector but GATED by the 6-axis
  liability veto. Decided UPSTREAM of the QSP/IL-6 layer, so it does NOT depend on the IL-6 result.
  (FINAL_NOMINATION_v7.md; A31b QSP window results.)
- **Model documentation system** — T1-T9 subsystem deep-docs + master + walkthrough + registers,
  every equation cited at file:line, adversarially reviewed, cross-validated against live engine
  (independent traces agree on the SAME 8 dead modules + SAME live binding path).
- **Live-path definition** — 12 IN-USE engine files; 8 DEAD modules named and excluded
  (unified_binding.py is a trap; cytokine_pbpk.py/il6_pbpk.py are NOT wired).
- **PK anchors** — teclistamab SC AFE 1.29x; elranatamab AFE 2.06x (route-matched, digitized-curve overlays).
- **Spatial viewer bundle** — 21 tissues, 3.4M cells (full count), 14 tissue images, 6-molecule
  validation layer with per-cell drug-distribution + kill fields.
- **Verification ledger** — one row per load-bearing claim, PASS/FLAG with source (NEW-21 corrected
  to engine-faithful-only, no un-reproducible external citation).

## VALIDATED WITH A KNOWN LIMITATION (report, do not hide)
- **IL-6 / CRS gate: FAIL on the 3 manifest-verified anchors** (mosun 127.4 / glofit 30.2 / talq 19.8).
  The cleanest same-statistic pair (glofit/talq, both medians) is WRONG DIRECTION; mosun/glofit is
  right-direction but magnitude-compressed. talquetamab INVERTS (model highest, clinical lowest):
  GPRC5D sits on plasma cells / keratinized tissue, not the myeloid-contacting compartments that drive
  CRS, and the myeloid IL-6 term does not condition on WHICH target compartment killing occurs in.
  This is a real mechanism gap, reported as-is (NOT rescued by denominator-shopping or dropping talq).
  teclistamab is NOT a verified anchor; any tecli-based ratio is PROVISIONAL.
- **Static R_costim** — resting-copy ranking under-rates inducible arms (4-1BB/OX40/ICOS) and yields a
  spurious 'CD2 wins' ordering. Disclosed everywhere; the costim signaling layer is WIRED BUT DORMANT
  in the 6 clinical CD3xTAA re-validation engagers.
- **Some parameters are FITTED, not measured** — IL-6 kdeg (0.20/hr), sigma_L (0.85, 4.25x the
  Shah-Betts platform value). Flagged in PARAMETER_AUDIT; never presented as measured.

## STILL UNDERWAY (NOT done at this submission state)
- **qsp_reproduction.ipynb** — the judge-runnable notebook reproducing every number from scratch
  (overlays, spatial-organ generation, voxel->receptor conversion, GRN, all figures). NOT built yet.
  This is the top outstanding deliverable. All INPUTS exist (engine, params, sqlite, bundles); the
  notebook that stitches them into a one-click reproduction is the remaining work.
- **Massive PBPK validation runpod** — the 22-TCE + 55-PK-breadth batch (other lane). Prompt +
  data map ready (OTHER_LANE_RUN_PROMPT.md); the 55-set harness gap (build_um_any not in model/) is
  the flagged risk.
- **Simulated construct variants** — the dashboard's future extension (schema-compatible with the
  validation layer); awaits the Modal construct sweep.
- **OX40/GITR net-negative-kill figure** — PROVISIONAL only (in-chat, structurally impossible under
  dormant signaling); omitted from artifacts.
