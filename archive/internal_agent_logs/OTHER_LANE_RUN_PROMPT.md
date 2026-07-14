# START THE MASSIVE PBPK VALIDATION + EXPLORATORY RUN

All paths under /media/balthasar-lab/RAID4/costim_engager_counterscreen/.
Full data map + readiness table: HANDOFF_DATA_MAP.md (read it first).

## TWO SETS — DON'T CONFLATE THEM

1) TCE PK+PD set (22 molecules) — full per-cell kinetics, this is the PD-capable set.
   PARAMS:  model/rundir/handoff/eng_params_normalized.json
   RUNNER:  model/engine/run_tce_pd_reval.py
   INVOKE:  cd model/rundir && python engine/run_tce_pd_reval.py $PWD <molecule>
            (KWS = model/rundir; writes handoff/tce_pd_<molecule>.json, own file per molecule)
   NOTE: the runner's ENG dict is the authoritative dose/route/il6_obs source; kinetics are
         merged in from eng_params_normalized.json at cfg=ENG[name].
   READY: 16 fully ready + 4 Fc-less (blinatumomab/solitomab/pasotuxizumab/tebentafusp:
          fFcRn=0.0 is CORRECT, do NOT change). HOLD/drop 2: REGN5459, forimtamig (no regimen).

2) PK-breadth set (55 molecules, incl. non-TCE reference mAbs for PK validation only).
   SPECS:    model/params/all_specs_comprehensive.json  (55 mol: MW, fFcRn, F_sc, ka, route, dose)
             model/params/MOLPROPS.json  (58 mol: MW + fFcRn)
             model/params/pk_specs_expansion.json  (19 mol expansion)
   CLINICAL: model/params/mab_tce_pkpd.sqlite  (51 distinct drugs with digitized PK C-t curves,
             15 with PD curves) — the validation ground truth to overlay against.
   ⚠ HARNESS GAP: run_tce_pd_reval.py does NOT run the 55-set (its ENG dict is TCE-only).
     The expanded-PK harness (build_um_any path) is NOT in model/ — it was in the earlier
     workspace. You must either locate/port it, or confirm your pod infrastructure
     (reference_unified_binding/pod_shard_runner.py, xeon_screen_launcher.py) carries it.
     Do NOT assume the TCE runner will accept non-TCE mAb names — it will KeyError on ENG.

## HOW PK VALIDATION WORKS (architecture reminder)
PK and PD are ONE unified simulation: plasma C(t) drives PD in the same loop. "PK validation"
= run the molecule, overlay model plasma C(t) vs the digitized clinical curve from the sqlite,
report AFE. Clean anchors already validated: teclistamab SC (AFE 1.29x), elranatamab (2.06x).

## RUN HYGIENE (learned this session — do not skip)
- CUDA_VISIBLE_DEVICES=""   (RTX A2000 driver too old; CUDA init elsewhere makes ACC=cpu insufficient)
- CAP THREADS PER PROCESS: OMP_NUM_THREADS / MKL_NUM_THREADS / OPENBLAS_NUM_THREADS.
  A single uncapped run spawns ~137 threads; N concurrent uncapped runs THRASH the box
  (4 concurrent TSIM=24 runs took >1h wall vs ~23min solo). For a pod fan-out, one run per
  container with threads capped to that container's core count.
- TSIM_DAYS=<n> caps sim length; for IL-6/CRS peak use 24 (PK is TRUNCATED at 24d, don't quote steady-state).
- Each molecule = its own process, its own output JSON (parallel-safe). Save EVERYTHING (PK+PD+spatial).
- Launch via setsid; kill by EXACT pid only (never pattern — ps|grep self-matches the launcher).

## SUGGESTED FIRST WAVE
20 TCE runs (16 ready + 4 Fc-less) through run_tce_pd_reval.py, one process per molecule,
threads capped, TSIM=24. In parallel, stand up / confirm the 55-set PK harness against the
sqlite curves. Hold REGN5459 + forimtamig until regimens are entered.

## THE VALIDATION BACKBONE — model/params/mab_tce_pkpd.sqlite
A structured pharmacology parameter database built by reading clinical PK/PD literature +
FDA labels and extracting the numbers a QSP model is validated against. This IS the ground
truth you overlay model outputs against.

  Drugs curated                              76
  Literature sources                         93  (48 DOI, 43 PMID, 45 FDA labels)
  Clinical populations (dose/disease/N)      91
  PK/PD curves                               183 (114 PK + 69 PD)
    - digitized from published figures        70 curves
  Digitized time-concentration points        936 (673 flagged digitized)
  SPR / binding-kinetics records (kon,koff,KD) 14, verbatim from Biacore tables

Tables: sources, drugs, populations, curves, timeseries, kinetics.
Curves carry readout_class (PK|PD) + qc status; timeseries rows carry the digitized points.
Use kinetics (kon/koff/KD) to cross-check the per-molecule params in eng_params_normalized.json;
use curves+timeseries as the overlay target for every model PK/PD run (report AFE per curve).
