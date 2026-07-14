# HANDOFF DATA MAP — PBPK/QSP validation + exploratory run

All paths absolute under `/media/balthasar-lab/RAID4/costim_engager_counterscreen/`.
Written 2026-07-13 for the lane launching the massive validation + exploratory run.

## 1. Per-molecule kinetic parameters (22 molecules) — THE param source
`model/rundir/handoff/eng_params_normalized.json`
Per molecule: kon_CD3_perM_pers, koff_CD3_pers, kon_TAA_perM_pers, koff_TAA_pers,
KD_CD3_nM, KD_TAA_nM, kint_perday, kdeg_perday, ksyn_frac_perday,
n_arm_CD3, n_arm_TAA, mw_kda, fFcRn, route, regimen_mg_day, measured_vs_derived, sources.

## 2. Digitized clinical PK/PD — the validation ground truth
`model/params/mab_tce_pkpd.sqlite`  (114 PK curves, 69 PD curves, 936 timeseries pts, 76 drugs)
- 11 molecules with route-matched digitized PK C-t curve: blinatumomab, catumaxomab,
  elranatamab, epcoritamab, glofitamab, linvoseltamab, mosunetuzumab, talquetamab,
  tarlatamab, tebentafusp, teclistamab.
- 11 digitized IL-6 curves (the PD/CRS anchors).

## 3. The runner
`model/engine/run_tce_pd_reval.py`
Invoke:  python run_tce_pd_reval.py $KWS <molecule>   with  KWS=model/rundir
Writes handoff/tce_pd_<molecule>.json (parallel-safe, own file/molecule).
Dose/route/il6_obs from the ENG dict at top of runner; kinetics merged from eng_params_normalized.json.
Env: CUDA_VISIBLE_DEVICES="" (RTX A2000 driver too old), TSIM_DAYS=<n>, PD_OUT_TAG=<tag>.
THREADING: cap threads per process (OMP_NUM_THREADS/MKL_NUM_THREADS/OPENBLAS_NUM_THREADS)
  so N concurrent runs do not over-subscribe cores — a single uncapped run spawns ~137 threads.

## 4. Regimen schedules
`model/engine/regimen_schedules_final.json`  (mirrored in model/rundir/handoff/)

## 5. Per-organ single-cell ABM (spatial substrate)
`model/rundir/handoff/agents/*.npz`  (18 files: 12 organs + 5 tumors + heme + blood)

## READINESS TABLE (verified 2026-07-13)
RUN-READY (16): acapatamab, alnuctamab, catumaxomab, cevostamab, cibisatamab,
  cinrebafusp_alfa, elranatamab, epcoritamab, glofitamab, linvoseltamab, mosunetuzumab,
  odronextamab, runimotamab, talquetamab, tarlatamab, teclistamab

Fc-LESS — fFcRn=0.0 is CORRECT, do NOT "fix" (4): blinatumomab (54 kDa BiTE),
  solitomab (55 kDa BiTE), pasotuxizumab (55 kDa BiTE), tebentafusp (77 kDa ImmTAC).

NOT run-ready — HOLD or drop (2):
  - REGN5459  : regimen_mg_day = None; TAA KD class estimate; CD3 KD=500 placeholder.
  - forimtamig: regimen_mg_day = None (real GPRC5D KD, no dose schedule wired).

## Not re-verified this session
Per-molecule PK-overlay dose/route alignment for all 11 PK-curve molecules NOT re-run here
(only params + provenance confirmed present). Clean anchors already validated:
teclistamab (SC, AFE 1.29x), elranatamab (matched PK+IL-6, AFE 2.06x).
