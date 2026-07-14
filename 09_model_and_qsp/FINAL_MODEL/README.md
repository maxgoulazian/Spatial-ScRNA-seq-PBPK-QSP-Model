# FINAL_MODEL — Whole-Body Per-Cell Spatial PBPK–PD/QSP Engine for Costim T-Cell Engagers

A self-contained, runnable copy of **only the components that execute in a model run**. No dead code, no
accumulated run outputs. Built as a clean duplicate so it never interferes with active production runs.

## Run it
```bash
python run.py teclistamab            # -> handoff/tce_pd_teclistamab.json
python run.py mosunetuzumab --tsim 24
```
See **OPERATION.md** for the full operating manual (options, engager list, folder layout, fidelity notes).

## What's here
| Path | Contents |
|---|---|
| `run.py` | Entry point — sets paths/env, calls the harness |
| `OPERATION.md` | Operating manual: exactly how to run and interpret it |
| `engine/` | The 13 live modules (import-traced; nothing dead) |
| `handoff/` | All input data (Rtot, params, anchors, regimens) + `agents/` (per-cell Xenium ABMs) |
| `params/` | Target routing, MW/fFcRn, clinical PK/PD database (`mab_tce_pkpd.sqlite`) |
| `MODEL_DOCUMENTATION/` | Full method: subsystem docs T1–T9 + T2b, every equation at file:line |
| `qsp_reproduction.ipynb` | Reproduces every analysis number from committed artifacts |
| `QSP_RESEARCH_REPORT.md` | Publication-grade writeup of the model |

## The model in one paragraph
The Shah & Betts 2-pore/FcRn PBPK backbone carries plasma and lymph, but each tissue's interstitial
compartment is replaced by a spatial agent-based model in which every real Xenium cell runs its own Rhoden
multi-arm kinetic binding ODE (true kon/koff, receptor turnover). TMDD, cytotoxic killing, and mechanistic
myeloid IL-6 all emerge from per-cell binding on the tissue grid. 15 PBPK tissues + a target-matched tumor
compartment, ~1.9M agents, PK and PD in one coupled pass. Every receptor is present on every cell via an scVI
Tabula-Sapiens overlay. Nominated arms: **4-1BB (TNFRSF9) and CD27**; CD28 gated out on CRS + suppression +
proliferation.

## Provenance & honesty
Load-bearing code comments (PMIDs, calibration rationale, the anti-fiction refusals like the unsourced
4-1BB induction fold) are **retained deliberately** — they are the model's scientific provenance, not clutter.
Known limitations (IL-6 CRS validation-gate failure, static costim receptor density, organ EPCAM/DLL3/PMEL
imputation gap) are documented in OPERATION.md §6 and MODEL_DOCUMENTATION/.
