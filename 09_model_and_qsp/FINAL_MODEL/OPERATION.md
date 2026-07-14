# FINAL_MODEL — Operating Manual

A self-contained, whole-body, per-cell **spatial PBPK–PD/QSP** engine for costimulatory CD3×TAA T-cell
engagers. This folder contains **only the components that execute in a model run** — no dead code, no
accumulated run outputs. It runs from inside this folder with a single command.

---

## 1. Quick start

```bash
# from inside FINAL_MODEL/, using the project's venv python:
python run.py teclistamab
```

Output is written to `handoff/tce_pd_teclistamab.json` — a single file containing the PK curve, the PD
depletion trajectory, systemic IL-6, and per-organ kill for that engager.

### Options
```
python run.py <engager> [--tsim DAYS] [--dt DT] [--graph-cache DIR|0]
  --tsim   simulated days (default: engine default; use 24 for a full step-up cycle)
  --dt     timestep in days (default 0.01; 0.04 is ~4x faster and unconditionally stable)
  --graph-cache  organ-graph cache dir, or 0 to rebuild each run (default 0)
```

### Environment
- `run.py` sets `CUDA_VISIBLE_DEVICES=""` automatically (the engine is CPU-only; the local GPU driver is too
  old for CUDA init). No GPU is required.
- A full run builds the spatial transport graph for every gridded organ on first use (minutes). Set
  `--graph-cache <dir>` to cache and reuse them across runs.

---

## 2. What runs when you call `run.py`

`run.py` → `engine/run_tce_pd_reval.py` (the harness) drives `CoupledPerCellPD`, which couples:

1. **PBPK transport backbone** (Shah & Betts 2-pore + FcRn) — plasma and lymph, 15 tissues + a
   target-matched tumor compartment.
2. **Per-cell spatial interstitium** — each tissue's interstitial compartment is replaced by an agent-based
   grid of **real Xenium cells**; every cell runs its own **Rhoden multi-arm kinetic binding ODE**
   (true kon/koff, receptor turnover ksyn=kdeg·Rtot).
3. **Emergent pharmacology** — TMDD, cytotoxic killing (2×2 kinetic synapse), and mechanistic per-cell
   **myeloid IL-6** all emerge from per-cell binding on the grid. PK and PD are computed in **one coupled
   pass**: the molecule's own plasma exposure drives its own per-cell pharmacology.

**Binding is identical in PK and PD** — the same engine solves the sink and the synapse, so efficacy and
toxicity are mechanistically consistent by construction.

---

## 3. Runnable engagers

The normalized parameter set (`handoff/eng_params_normalized.json`) carries:

REGN5459, acapatamab, alnuctamab, blinatumomab, catumaxomab, cevostamab, cibisatamab, cinrebafusp_alfa, elranatamab, epcoritamab, forimtamig, glofitamab, linvoseltamab, mosunetuzumab, odronextamab, pasotuxizumab, runimotamab, solitomab, talquetamab, tarlatamab, tebentafusp, teclistamab

Each entry supplies measured/derived kon, koff, KD (TAA and CD3), kint, kdeg, MW, fFcRn, route, and regimen.
Molecules with `regimen=None` (e.g. REGN5459, forimtamig) are held pending a sourced regimen.

---

## 4. Folder layout

```
FINAL_MODEL/
├── run.py                    # entry point (sets paths + env, calls the harness)
├── OPERATION.md              # this file
├── engine/                   # the 13 live modules (nothing dead)
│   ├── run_tce_pd_reval.py   #   harness / entry
│   ├── coupled_percell_pd.py #   PD coupler (top of the import tree)
│   ├── coupled_percell_pk.py #   per-cell PK
│   ├── wholebody_pd.py       #   whole-body PD, kill, IL-6 wiring
│   ├── wholebody_percell.py  #   spatial transport (grid Laplacian, BEC/LEC)
│   ├── kinetic_rhoden_percell.py  # per-cell Rhoden kinetic binding
│   ├── kinetic_synapse.py    #   2×2 matrix-exp synapse killing
│   ├── multiarm_binding.py   #   multi-arm valency/span geometry
│   ├── myeloid_il6.py        #   mechanistic per-cell myeloid IL-6
│   ├── costim_induction.py   #   activation-induced costim (wired, default-OFF)
│   ├── signaling_dynamics.py #   GRN-driven signaling shape
│   ├── pd_model_config.py    #   SINGLE SOURCE OF TRUTH for engine + calibration
│   └── qsp_costim_window_v2.py # PBPK arrays + therapeutic-window / nomination
├── handoff/                  # ALL input data the engine reads
│   ├── *.json                #   Rtot, eng params, kinetic calib, anchors, regimens, ...
│   ├── agents/               #   per-organ ABM npz (real Xenium cells) + tumor + heme
│   ├── signaling_surrogate/  #   GRN program kinetics + per-arm drive
│   └── tumor_abm_cells.npz   #   tumor ABM (harness reads at handoff/tumor_abm_cells.npz)
└── params/                   # target routing, MW/fFcRn props, clinical PK/PD sqlite
```

---

## 5. Reproduction & validation

- **Reproduction notebook:** `../qsp_reproduction.ipynb` reproduces every analysis number from committed
  artifacts (executed end-to-end, 0 errors).
- **Full method + every equation:** `../MODEL_DOCUMENTATION_FINAL/` (subsystem docs T1–T9 + T2b, each with
  equations at file:line).
- **Clinical validation data:** `params/mab_tce_pkpd.sqlite` (76 drugs, 183 curves, 936 timepoints).

---

## 6. Notes on fidelity (read before interpreting a run)
- **IL-6 CRS axis** is mechanistic but **fails a validation gate** on verified clinical anchors (talquetamab
  inverts — GPRC5D sits on plasma cells/keratinized tissue, not the myeloid-contacting compartments that seed
  CRS). Treat IL-6 output as mechanistic-but-unvalidated. See `../MODEL_DOCUMENTATION_FINAL/`.
- **4-1BB (TNFRSF9)** is activation-induced; its induction fold is not sourced, so `costim_induction.py`
  refuses to guess (raises under `strict=True`). Costim runs therefore use **static resting** receptor
  density — a conservative floor for inducible arms.
- **Organ-level EPCAM/DLL3/PMEL** receptor imputation is a known gap (scVI reference not retained); tumor
  pools for those targets are complete.
