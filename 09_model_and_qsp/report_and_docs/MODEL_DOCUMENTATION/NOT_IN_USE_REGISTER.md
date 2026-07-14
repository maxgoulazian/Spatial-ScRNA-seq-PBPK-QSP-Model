# NOT-IN-USE REGISTER — code that exists but does NOT run

> **Purpose:** this repository contains modules that look like part of the model and are not. Some were
> superseded, some were built and never wired in. **A reader who assumes "it's in `engine/`, so it runs" will
> document, cite, and trust code that never executes.** This file is the authoritative list.
>
> **Regenerate it** (do not hand-edit the verdicts) by tracing the import graph from the entry point —
> see §3.

---

## 1. THE LIVE EXECUTION PATH (12 modules — these, and only these, run)

Traced from the entry point `engine/run_tce_pd_reval.py`:

```
run_tce_pd_reval
├── qsp_costim_window_v2          Shah-Betts PBPK backbone (organ physiology arrays)
├── pd_model_config               which PD kill law + its parameters
└── coupled_percell_pd            THE INTEGRATOR — assembles and steps a run
    ├── coupled_percell_pk        whole-body PK (per-organ vascular QSS + extravasation + lymph)
    │   └── wholebody_percell     per-cell whole-body engine
    │       └── kinetic_rhoden_percell    Rhoden bivalent binding (BE solver)
    ├── kinetic_rhoden_percell
    ├── myeloid_il6               mechanistic CRS IL-6 (emitters + plasma ODE)
    ├── wholebody_percell
    └── wholebody_pd              per-cell PD (killing + cytokine + Treg suppression)
        ├── costim_induction      activation-induced costim receptor density
        ├── kinetic_synapse       kinetic multivalent immune synapse
        ├── multiarm_binding      format geometry (cis/trans, reach, linker span)
        └── myeloid_il6
```

Documented in `subsystems/T1_*` … `subsystems/T9_*`.

---

## 2. NOT IN THE EXECUTION PATH

### 2a. DEAD — superseded or never wired. **Do not document. Do not cite. Do not trust.**

| module | why it is dead |
|---|---|
| **`cytokine_pbpk.py`** | A full per-organ PBPK mirror for IL-6/IFN/TNF/IL-2 (interstitium + vascular QSS + lymph + plasma, per organ). **Built 2026-07-13 and never wired into the engine.** The IL-6 the model actually computes comes from the far simpler `PlasmaIL6` ODE inside `myeloid_il6.py` (production ÷ ECF volume, first-order clearance). **If you have been told the model gives every cytokine its own PBPK compartment set — it does not.** This is the single most misleading file in the repository, because it is exactly what someone would *want* to be true. Wiring it is an open work item, not a completed one. |
| **`il6_pbpk.py`** | IL-6-specific ancestor of `cytokine_pbpk.py`. Superseded the same day it was written. Dead twice over. |
| **`unified_binding.py`** | Despite the name and the docstring ("ONE Rhoden-kinetic binding solve for every cell, every location"), it is **not reachable from the entry point.** The live binding path is `kinetic_rhoden_percell` + `multiarm_binding` + `kinetic_synapse`. The name is a trap. |
| **`multiarm_kinetic.py`** | Not reachable. The live format-geometry module is `multiarm_binding.py`. |
| **`biexact_solver.py`** | An inner solver for `multiarm_kinetic.step_mul` — dead because its only caller is dead. |

### 2b. OFFLINE — real, correct, and *not runtime*. They produce inputs; they are **provenance, not mechanism**.

These are legitimately part of the pipeline and must be documented **as data provenance**, but they do not
execute during a simulation. Document them in the data/provenance appendix, never as subsystems.

| module | what it produces | consumed as |
|---|---|---|
| `rna_to_receptor.py` | per-cell surface-receptor copy numbers from scRNA-seq (PI method, absolute scaling) | pre-computed receptor-copy files |
| `convert_copies_ALL.py` | unified per-cell receptor-copies conversion, no lineage restriction | ditto |
| `calib_kdeath.py` | the `k_death` / `k_hit` calibration against the epcoritamab depletion time-course (day 7 = 0.30, day 14 = 0.90, day 28 = …) | constants in `pd_model_config.py` |

> ⚠️ **`calib_kdeath.py` is where `k_death` comes from — which means `k_death` is `[FITTED]`, not measured.**
> Anything downstream of the kill law inherits that tag. It also still imports the retired `IL6_SCALE`; do not
> resurrect that path (see `PROVENANCE_AND_VALIDATION.md` §3).

---

## 3. HOW TO REGENERATE THIS FILE (do not trust a stale copy)

```bash
cd engine
python - <<'PY'
import re, os, collections
ENTRY = "run_tce_pd_reval.py"
local = {f[:-3] for f in os.listdir('.') if f.endswith('.py')}
def imports(f):
    src = open(f).read()
    return {m.group(1) or m.group(2)
            for m in re.finditer(r'^\s*(?:from\s+(\w+)\s+import|import\s+(\w+))', src, re.M)
            if (m.group(1) or m.group(2)) in local}
seen, stack = set(), [ENTRY[:-3]]
while stack:
    cur = stack.pop()
    if cur in seen: continue
    seen.add(cur)
    stack += [d for d in imports(cur + ".py") if d not in seen]
print("LIVE:", *sorted(seen), sep="\n  ✓ ")
print("DEAD:", *sorted(local - seen), sep="\n  ✗ ")
PY
```

**Run this before writing or believing any documentation.** A module's presence in `engine/`, the confidence
of its docstring, and the plausibility of its name are all worthless as evidence that it runs.
