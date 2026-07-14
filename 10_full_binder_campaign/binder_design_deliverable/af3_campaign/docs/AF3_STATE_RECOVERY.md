# AF3 BINDER CAMPAIGN — STATE RECOVERY INDEX
_Written 2026-07-13 09:21 before account swap. Workspace is swept on swap; artifact store is durable. This doc indexes every durable artifact + open threads._

## PROJECT
CD4 Perturb-seq costim counter-screen. FINAL NOMINATION: co-leads **4-1BB (TNFRSF9) + CD27** (liability-CLEAN). Deadline Jul 13 9pm.
Workstation (A2000, NOT sandbox); env `sc-analysis`; ChimeraX /usr/bin/chimerax; RFantibody at /media/balthasar-lab/RAID1/RFantibody.

## AF3 DE NOVO BINDER CAMPAIGN — COMPLETE
Funnel: 3,964 designs (1-seed screen) -> 150 (3-seed refold) -> 39 (10-seed final, incl 15 panel expansion).
Ran on 3 RunPod A100 pods (auk2ec4ihqc8my, f61l8x1ym7qicu, l3d8msotdlky8d) — ALL DATA HARVESTED, pods safe to terminate.

### Screen completion: 3,964 / 3,974 (99.75%)
10 folds dropped on pod2 (never scored, not in funnel): 7x 4-1BB (r0053,r0095,r0368,r0557,r0683,r0788,r0809 backbones 96,96,68,28,78,86,62), 3x CD3 (r0032,r0200,r0326 backbones 3,80,87). Decision pending: leave as 99.75% or refold 10 (needs a GPU).

### DURABLE ARTIFACTS (artifact_id -> latest version_id)
RAW STRUCTURE HARVEST (checkpoints, .tgz — untar to recover .cif models):
- Screen 3964 full models: pod1 1795c55a-3819-44f4-b105-ddbd631c2be5 / pod2 de074ea8-b70e-49a2-b1c1-c23b13362fd3 / pod3 1c980ac9-d1a5-4e8e-ab45-702449f8257c (results_final.tgz, ~235MB each, 1325/1315/1324 model.cif+summary)
- Tenseed 39 full (per-seed): pod1 14de1e36-ddc6-4937-9439-18595002e0aa / pod2 76f840b6-f7bb-416a-807b-db6688585baf / pod3 aa2b34f4-06dd-46a9-9915-b15cd7affc4e
- Refold 150 (3-seed): artifact 9e8cbb8d-c1b7-43a7-9c32-cd217be5eaff (refold_all_results.tgz)
- Tenseed 24-orig: artifact 4908f65b-42a4-4641-b4f2-ee81f31018bd (tenseed_all_results.tgz)

SCORED TABLES (readable CSV):
- af3_master_scored.csv (3964 screen, 1-seed): version 93c7c8ac-bde7-4227-9b58-af859a7ee819
- refold_scored.csv (150, 3-seed): artifact bca29901-7029-49c2-9bfe-b887966c0571 version c57afcc8-6e3f-4878-b2bd-8d8743f6abe9
- tenseed_final_scored.csv (24 orig, 10-seed): artifact 28b9a066-84bf-424c-824b-768b0c7509d0 version 9f522fac-2197-4f16-a6e0-1e096032b405
- scrmsd_final.csv (fold scRMSD, CORRECTED v2): artifact e3940c50-e200-4fba-81d9-02f814219427 version 2a9fb3c7-6d05-4d31-a932-769bbb1640f2
- finalist_combined_verdict.csv (CORRECTED v2): artifact 1f47bb79-a871-414e-96b9-f8969b015fd2 version 03829d25-502f-4a4e-8297-0b7d71c112c4

DOCS + FIGURES:
- AF3_FINAL_FINALISTS.md: artifact 7bb3d5d1-f84d-4e2f-a5ce-f813d871d265 version e335a9b4-bce6-405d-b954-c0669b85e0e6 (seed-stability claim CORRECTED)
- AF3_SCREEN_RESULTS.md (v2 MIN): version 2e9de794-46b1-4619-b38b-82c3558b2bab
- fold_overlay_41bb_r0867.png (binder-on-binder scRMSD overlay): artifact 8c8fb2c8-6d1b-40a4-b495-628e01bd3b73 version c7d5d4fb-4929-436e-80ee-39d854fe002d
- Top finalist .cif (render in 3D): 41bb_r0867 27817fb8-b678-4811-b197-2332f798948e / cd27_r0011 51f0769e-4d52-4f23-ad63-5ae57889ec8d / cd3_r0478 78c46746-7621-4add-aab3-c2fb74cfb2d4 / cea5_r0075 ad4c7f5d-9bc8-45c6-b68e-63d532da7898

REPRODUCIBILITY (this save): af3_campaign_reproducibility.tgz (all runner scripts, pod maps, parse code, panel defs, handoffs), af3_design_backbone_map.txt (3974-row design->pod->backbone map, needed for RMSD), finalist_panels.csv (top-12/target panel defs).

## FINALIST VERDICT (validated, corrected)
Metrics: ipSAE(MIN-of-directions, primary interface), iPTM/pTM/pDockQ/LIS, 10-seed consistency (rank_mean/std), fold-scRMSD (binder-on-binder Kabsch, <2A=pass, ChimeraX-validated).
Per-target combined-pass (fold_scRMSD<2A AND ipSAE>0.35): CEA5 6/6, CD3 3/6, CD27 2/6, 4-1BB 1/6.
CO-LEAD PICKS (pass both gates):
- 4-1BB: 41bb_r0867_41bb_62 (VHH) ipSAE 0.492, fold 0.86A, pDockQ 0.43, 10-seed mean 0.80 std 0.023. Backup 41bb_r0025 (fold 0.92A, most seed-stable in group std 0.008).
- CD27: cd27_r0011_cd27_24 (VHH) ipSAE 0.553, fold 0.60A, 10-seed mean 0.79. Backup cd27_r0846 (fold 0.44A, ipSAE 0.451, passes both).
- CD3 (redirector): cd3_r0478_cd3_65 (VH/VL) ipSAE 0.491, fold 1.00A.
- CEA5 (TAA): cea5_r0075_cea5_10 (VH/VL) ipSAE 0.587, fold 1.03A, 10-seed mean 0.88 std 0.012 (cleaner than r0420 which has ipSAE 0.65 but std 0.114).

## RMSD METHOD (settled, don't re-litigate)
fold-scRMSD = superpose designed binder chain H onto AF3 binder chain H (BINDER ONLY, target excluded), one Kabsch, RMSD=sqrt(mean(sum((P-Q)^2))). This is fold self-consistency (does MPNN seq fold into RFdiff backbone). 41bb_r0867=0.86A verified vs ChimeraX matchmaker 0.887A.
dock_rmsd (secondary) = target displacement after binder-align = epitope/angle shift; large values EXPECTED (RFdiffusion samples diverse docking modes, AF3 is better docking predictor per RFantibody/Bennett Nature 2026). NOT a failure criterion.
BUG FIXED: earlier batch used sqrt((diff^2)).sum().mean()^0.5 (double-sqrt, wrong) — corrected to single sqrt. No verdict changed.

## OPEN THREADS
1. TERMINATE 3 PODS (all harvested). IDs above.
2. Decide: leave screen 99.75% or refold 10 stragglers.
3. Finalize per-target TESTABLE PANELS (user wants finalist GROUPS not single leads): panel defs in finalist_panels.csv (top-12/target), 39 have 10-seed. Honest graded sizes: 4-1BB thin (only 23 real binder candidates in whole screen), CEA5 deep.
4. Optional: CD27 co-lead overlay figure (like 41bb one); CDR-highlighted overlays for deck.
5. 4-1BB CRD1 agonist-epitope compatibility check on r0867 (vs urelumab non-blocking site) — not yet run this session.
