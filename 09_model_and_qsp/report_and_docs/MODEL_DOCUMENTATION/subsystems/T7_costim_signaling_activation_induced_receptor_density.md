---
title: "T7 — Costim signaling & activation-induced receptor density"
subsystem: T7
model: costim_engager_counterscreen (whole-body per-cell PBPK-PD)
date: 2026-07-13
generated_by: workflow-subagent T7
live_source_files:
  - engine/costim_induction.py   (86 lines, mtime 2026-07-13 13:54)
  - engine/wholebody_pd.py       (503 lines, mtime 2026-07-13 14:42)
status: "NEW 2026-07-13. Mechanism BUILT; NOT RUNNABLE (see §0 — FOUR BLOCKERS). Fold values NOT SOURCED."
revised: "2026-07-13 — adversarial verification pass. 4 defects fixed in-place (EQ-0 arm/slot labels were
  swapped; Euler bound 1.44 h -> 1.475 h; 'engine default dt' misattributed; R_costim mistagged [MEASURED]).
  2 new code defects found and added: BLOCKER-4 (4 of 11 screened arms absent from COSTIM_INDUCTION, fail
  OPEN at fold=1.0) and the QSS-only costim->IL-6 gain leak (EQ-5). All citations re-checked against the
  live tree: ZERO fabricated citations. See §6.1."
---

# T7 — Costim Signaling & Activation-Induced Receptor Density

> Sibling docs: the TCE kill core (Schropp ternary trimer, `ternary_equilibrium`, `NM_PER_COPY`, Treg
> damping), the myeloid IL-6 subsystem (`myeloid_il6.py`), the kinetic-synapse engine
> (`kinetic_synapse.py`), and the multi-arm format geometry layer (`multiarm_binding.py`). This doc owns
> ONLY the costim arm: per-cell costim occupancy → signaling programs → kill/cytokine modulation, and the
> activation-induction of the costim receptor itself.

---

## 0. STATUS — FOUR BLOCKERS (read this before anything else)

This subsystem is **code-complete but NOT EXECUTABLE at this HEAD.** All four blockers below were
established by *running the code in this task*, not by reading it. A reviewer will find them in minutes; we
state them first.

### BLOCKER-1 (FATAL) — `signaling_dynamics` does not exist

`engine/wholebody_pd.py:192` executes `import signaling_dynamics as _sigmod` **unguarded**, inside the
`if costim_arm is not None and R_costim_percell is not None:` branch (wholebody_pd.py:191).

**There is no `signaling_dynamics.py` anywhere in the repository.** Verified by `find` over the whole
project tree (zero hits) and by Python import resolution with both `engine/` and `rundir/handoff/` on
`sys.path` (`ModuleNotFoundError`).

RUN-VERIFIED this task (constructing a real `OrganPD` with 3 T cells / 2 targets / 1 macrophage):

| Construction | Result |
|---|---|
| `OrganPD(...)` **no costim arm** | **OK** — `kill_frac=0.0028471`, `sig=None`, `costim_ind=None` |
| `OrganPD(..., costim_arm='TNFRSF9', R_costim_percell=Rc)` | **`ModuleNotFoundError: No module named 'signaling_dynamics'`** |
| `OrganPD(..., costim_arm='CD28',    R_costim_percell=Rc)` | **`ModuleNotFoundError: No module named 'signaling_dynamics'`** |

Consequence: **every** costim construct — inducible *and* constitutive — raises at `OrganPD.__init__`.
The entire contents of §2 EQ-6…EQ-16 and the whole of `costim_induction.py` are **unreachable in the
current tree**, because `CostimInduction` is instantiated at wholebody_pd.py:220–222, which is *after* the
fatal import at :192. The plain-TCE (no-costim) path is unaffected and runs — so the clinical-validation
molecules are fine, and the counterscreen is the thing that is down.

`PerCellSignaling(n, arm, drive_source="lane")` (wholebody_pd.py:213) and the "lane drive table" its
docstring refers to (wholebody_pd.py:185) are therefore **SOURCE NOT FOUND IN CODE — TBD**. Every
program magnitude this doc describes (`effector`, `suppression`, `exhaustion`, `cyto_IFNG`, `cyto_TNF`,
`cyto_IL2`) is produced by that missing module. **We cannot document the program dynamics, their
timescales, or their provenance, because the code that computes them is not in the tree.** The doc below
documents what `wholebody_pd.py` *does with* those program outputs, which is all that is auditable.

### BLOCKER-2 (BY DESIGN — the honest one) — the FOLD values are NOT SOURCED

`costim_induction.py:39–42` sets `fold=None, source="PENDING"` for all four inducible arms. With the
default `strict=True` (costim_induction.py:48, and the call site wholebody_pd.py:222 passes `strict=True`
explicitly), the constructor **raises `ValueError`** rather than silently running them at the resting
density. RUN-VERIFIED this task:

| Arm | Gene | `fold` | `kind` | `source` | `strict=True` construction |
|---|---|---|---|---|---|
| 4-1BB | `TNFRSF9` | `None` | inducible | `"PENDING"` | **ValueError** |
| OX40 | `TNFRSF4` | `None` | inducible | `"PENDING"` | **ValueError** |
| ICOS | `ICOS` | `None` | inducible | `"PENDING"` | **ValueError** |
| GITR | `TNFRSF18` | `None` | inducible | `"PENDING"` | **ValueError** |

This is fail-closed **for the four arms that are in the table**: for those, the module refuses to guess. But
it means **those inducible arms — the entire scientific point of the counterscreen — cannot currently be run
at all.** This is an **OPEN BLOCKER, not a completed feature.** The module's own header says so
(costim_induction.py:31–32: *"PENDING SOURCING (agent running) — every entry MUST carry a PMID before this
is used for a claim"*).

⚠ **The fail-closed guard is NOT complete — see BLOCKER-4.** It protects only arms present in
`COSTIM_INDUCTION`. Four of the eleven screened arms are absent from the table entirely, and for those the
module **does** guess (fold=1.0), `strict=True` notwithstanding.

`t_peak_h` and `t_decay_h` are `None` for **every** arm including the constitutive ones — so even a sourced
`fold` would still run on the **fallback kinetics** `k_on=1.0/hr`, `k_off=0.05/hr` (costim_induction.py:70–71),
which are `[ASSUMED]` round numbers with no source.

### BLOCKER-3 (SILENT) — induction is inert on the kinetic-synapse path

`self.costim_ind.step(...)` is called at **exactly one** site: `wholebody_pd.py:334`, inside `step()` — the
**QSS** path. Grep of every `costim_ind` reference in the module:

```
190:  self.costim_ind=None          # declaration
221:  from costim_induction import CostimInduction as _CI
222:  self.costim_ind = _CI(costim_arm, len(self.Tidx), Rc[self.Tidx], strict=True)
224:  self.costim_ind = None        # ImportError fallback
332:  if self.costim_ind is not None:
334:  self.R_costim[Tidx] = self.costim_ind.step(dt, p_eng)     <-- ONLY step() call
```

`_step_kinetic()` (wholebody_pd.py:438–491) — the path taken when `pd_kinetics=True` — reads
`Rc_T=self.R_costim[Tidx]` at :447 **without ever calling `costim_ind.step()`**. So under
`PD_KINETICS=True`, `R_costim` stays frozen at the **resting** density forever: *the exact bug this module
was written to fix silently returns.* The QSS path is fixed; the kinetic path is not.

### BLOCKER-4 (SILENT, FOUND IN ADVERSARIAL REVIEW) — 4 of the 11 screened arms are NOT in the induction table

The screen runs **eleven** costim arms. The runner names them explicitly
(run_tce_pd_reval.py:150–151): *"11 costim arms have per-cell density in all 11 organs (verified 2026-07-13):
TNFRSF9 (4-1BB), TNFRSF4 (OX40), TNFRSF18 (GITR), CD28, ICOS, CD27, CD40, TNFRSF25, CD2, CD226 (DNAM-1),
TNFRSF14 (HVEM)."*

`COSTIM_INDUCTION` (costim_induction.py:33–43) contains **seven**. The missing four — **CD40, TNFRSF25 (DR3),
CD226 (DNAM-1), TNFRSF14 (HVEM)** — hit the `p is None` branch (costim_induction.py:53–57), which sets
`fold=1.0, kind="unknown", source="NOT IN TABLE"` and writes a **stderr warning only**. There is **no
`ValueError`, and `strict=True` does not apply to this branch at all** — the `strict` check lives inside the
`else:` (costim_induction.py:60–67), i.e. it is reachable *only* for arms already in the table.

RUN-VERIFIED this task (`strict=True` for every call):

| Arm (screened) | In `COSTIM_INDUCTION`? | Resulting `fold` | `kind` | `source` | Behaviour |
|---|---|---|---|---|---|
| CD40 | **no** | **1.0** | `unknown` | `NOT IN TABLE` | **runs, silently constitutive** (stderr only) |
| TNFRSF25 (DR3) | **no** | **1.0** | `unknown` | `NOT IN TABLE` | **runs, silently constitutive** |
| CD226 (DNAM-1) | **no** | **1.0** | `unknown` | `NOT IN TABLE` | **runs, silently constitutive** |
| TNFRSF14 (HVEM) | **no** | **1.0** | `unknown` | `NOT IN TABLE` | **runs, silently constitutive** |

**This directly contradicts the module's stated contract.** Its own header (costim_induction.py:23–24) says
*"Any arm not in the table falls back to FOLD=1 (constitutive / unknown) and is FLAGGED, never guessed"* —
but **fold=1.0 IS a guess**: it is the assertion *"this receptor is not activation-induced"*, made with no
data, for four receptors whose induction status was never assessed. **TNFRSF25 (DR3) and CD226 are not
obviously constitutive.** This is the *same class of error* T7 was built to eliminate (BLOCKER-2's
"under-rate the inducible arms"), merely relocated from the four table arms to the four non-table arms —
and unlike BLOCKER-2 it **fails OPEN**: the run proceeds and produces a ranking.

**Fix:** every arm the runner can screen must be an explicit row in `COSTIM_INDUCTION`; an unknown arm must
`raise` under `strict=True`, not warn. Until then, **any ranking that includes CD40 / TNFRSF25 / CD226 /
TNFRSF14 is asserting they are constitutive without evidence.**

---

## 1. PURPOSE & DATAFLOW POSITION

### What this subsystem does

T7 is the **costim arm** of the counterscreen. A costim engager is a multi-specific whose third arm
agonizes a T-cell costimulatory receptor (4-1BB, OX40, ICOS, GITR, CD28, CD2, CD27, …) on top of the
CD3×TAA bridge. T7 answers two questions the plain TCE model cannot:

1. **(a) Signaling** — given a T cell's *own* costim-receptor occupancy at its *own* local drug
   concentration, what does that do to its killing and its cytokine output? Occupancy drives per-cell
   **programs** (`effector`, `suppression`, `exhaustion`, and three cytokine programs), which multiply
   the trimer-derived kill and the cytokine emission. This **replaces the legacy scalar `costim_boost`**
   (a single `1+boost` multiplier, wholebody_pd.py:319–320).

2. **(b) Activation-induction of the receptor itself** — 4-1BB (TNFRSF9), OX40 (TNFRSF4), ICOS and GITR
   are **essentially absent on resting T cells** and appear **only after TCR engagement**. That
   conditionality **is why the field targets them**: the costim fires only on T cells that have already
   engaged a target, i.e. it is **tumour-conditional**.

### The bug this fixes (the reason T7 exists)

The model read costim receptor density **once**, from **resting** healthy-tissue scRNA-seq, and **never
updated it** (`self.R_costim` set in `__init__`, read thereafter). For the inducible arms this reads them
at ~zero forever. Two compounding errors (costim_induction.py:9–11, wholebody_pd.py:323–327):

- it **systematically UNDER-RATES** exactly the arms that matter (4-1BB, OX40, ICOS, GITR), and
- it **OVER-RATES** the constitutive ones (CD28, CD2, CD27), which *are* correctly described by a resting
  density.

The screen would have produced a **confidently-wrong ranking** — the code comment names it: *"CD2 beats
4-1BB"* — *"that any immunologist would reject on sight."*

### The second bug this fixes — `costim_arm` was never plumbed through

Independently of the biology: `attach_pd()` **has always accepted** `costim_arm`
(`coupled_percell_pd.py:17–18`), but **the runner never passed it**. `run_tce_pd_reval.py:146–149` states
it verbatim:

> *"COSTIM ARM: attach_pd has always supported it (costim_arm -> per-cell '<ARM>_copies' + the per-cell
> signaling engine), but the runner never PASSED it — so **every costim construct silently ran as a plain
> TCE**. Thread it from cfg so the costim screen is actually possible."*

It is now threaded at `run_tce_pd_reval.py:152–154` (`costim_arm=cfg.get('costim_arm')`). `cfg['costim_arm']=None`
(the clinical-validation molecules) reproduces the previous plain-TCE behaviour **exactly** — the change is
strictly **additive**. Note the grim implication: **before this fix, every costim result in the project's
history was a plain-TCE result wearing a costim label.** Note also that this fix is what *exposes*
BLOCKER-1: passing `costim_arm` for the first time is what first executes wholebody_pd.py:192.

### Position in the life of the molecule

```
   PK (coupled_percell_pk)  ──►  per-cell local free drug  C_i (nM)
                                        │
                    ┌───────────────────┼────────────────────┐
                    ▼                   ▼                    ▼
            CD3 arm (R_CD3)      TAA arm (R_TAA)      COSTIM arm (R_costim)   ◄── T7 OWNS THIS
                    └────────┬──────────┘                    │
                             ▼                               │
              Schropp ternary trimer Cb  (wholebody_pd.py:85–105)
                             │                               │
              ┌──────────────┴───────────┐                   │
              ▼                          ▼                   ▼
      p_eng = Cb/(Cb+R_CD3)        (kill species)      occ = C/(C+KD_cos)·(R_cos/anchor)
        (wholebody_pd.py:333)                            (wholebody_pd.py:336)
              │                                              │
              │  ┌───────────────────────────────────────────┤ cis/trans coincidence gate (:337)
              │  │                                           ▼
              │  │                          PerCellSignaling.step(occ, dt·24)  ◄── MODULE MISSING
              │  │                                           │
              │  │                       ┌───────────────────┼──────────────────┐
              │  │                       ▼                   ▼                  ▼
              │  │                 effector p_E        suppression p_S    exhaustion p_X
              │  │                  cyto_IFNG/TNF/IL2       │                  │
              ▼  ▼                       │                  │                  │
    ┌─────────────────────────┐          ▼                  ▼                  ▼
    │  CostimInduction.step   │    g_eff = e^{0.55 p_E}·e^{-0.30 p_X}   supp_extra
    │  da/dt = k_on·p_eng·(1-a)│         │                              │
    │        - k_off·a         │         └──────────┬───────────────────┘
    │  R(t)=R_rest(1+(F-1)a)   │                    ▼
    └───────────┬─────────────┘        kill_T = Cb·g_eff / (1 + 0.25·n_treg·(1+supp_extra))
                │  FEEDBACK LOOP                    │
                └──► R_costim (:334) ──► occ ──►    ▼
                                              hazard → kill_frac ;  cytokine × cyto_sig_gain
```

**The load-bearing structural fact:** induction is driven by **that T cell's OWN engaged fraction**
(`p_eng`), not by a global signal, not by a tumour flag, not by an organ label. A T cell that never engages
never upregulates. **The tumour-conditionality is therefore EMERGENT, not imposed** (costim_induction.py:17–18,
wholebody_pd.py:329–330). This is the single most important design claim in this subsystem, and it is the
one thing a reviewer should be told to check.

**Feeds T7:** per-cell local drug `C_i` (from the transport grid); per-cell costim copies
`R_costim` from the agent table column `'<ARM>_copies'` (`coupled_percell_pd.py:26–27, 38–43`); the per-cell
trimer `Cb` (kill core); the T-cell / CD8 / Treg masks (wholebody_pd.py:156–159).

**T7 feeds:** the per-cell kill hazard (wholebody_pd.py:350–362), the organ cytokine rate
(wholebody_pd.py:366–374), and — through the induced `R_costim` — **itself**, on the next step.

---

## 2. GOVERNING EQUATIONS

Notation: `i` = T-cell index; `a_i ∈ [0,1]` = activation memory; `p_eng,i ∈ [0,1)` = engaged fraction;
`R_rest,i` = resting costim copies/cell; `FOLD` = activation fold-upregulation (dimensionless);
`C_i` = local free drug (nM); `Cb_i` = ternary trimer (nM); `dt` = PD step (days); `dt_h = 24·dt` (hours).

---

### EQ-0 — INPUT: Schropp ternary trimer (wholebody_pd.py:98–104)

Not owned by T7 (it is the kill core), but reproduced because **it is the driver of both `p_eng` and the
cis-coincidence gate**, so T7 cannot be read without it.

```
a    = (1 + C/KD2)·C / (α·KD1·KD2)                      (Eq.28, :98)
b    = C·(R_A − R_B)/(α·KD1·KD2) + (1 + C/KD2)(1 + C/KD1)   (Eq.29, :99)
d    = −R_B·(1 + C/KD1)                                 (Eq.30, :100)
R_Bf = (−b + sqrt(max(b²−4ad, 0))) / (2a)               (Eq.27, :102)  free TAA-side receptor
R_Af = R_A / (1 + C/KD1 + R_Bf·C/(α·KD1·KD2))           (Eq.26, :103)  free CD3-side receptor
Cb   = C·R_Af·R_Bf / (α·KD1·KD2)                        (Eq.33, :104)  bridged trimer
```

**Arm→slot mapping (read this before the labels above — it is easy to get backwards).** The call site is
`ternary_equilibrium(Cd, RA, RB, self.KD_CD3, self.KD_TAA)` (wholebody_pd.py:313) with `RA=self.R_CD3[Tidx]`
(:308) and `RB=self.syn_TAA_mean` (:309). Therefore **`R_A` = CD3 with `KD1` = `KD_CD3`**, and **`R_B` = TAA
with `KD2` = `KD_TAA`**. So `R_Af` is the **free CD3** and `R_Bf` is the **free TAA** — the quadratic is
solved on the **TAA** side, not the CD3 side.

- **Biological meaning:** one drug molecule bridges a CD3 on the T cell to a TAA on a neighbouring target.
  `Cb` is the concentration of that bridged trimer in the ~71 pL synapse volume.
- **Mechanistic rationale:** free receptors are solved from the QE quadratic **before** the trimer, so the
  **prozone/hook effect is emergent** — at high drug both arms saturate as binary complexes and the trimer
  collapses. The code names the rejected alternative (:90–91): the reduced-linear form
  `C·R_A·R_B/(KD1·KD2)`, which has **no free-receptor depletion and therefore no prozone**, and
  "overstated kill/cytokine at high drug."
- **Units:** `C`, `KD1`, `KD2`, `R_A`, `R_B`, `Cb` all nM; `α` dimensionless (default `1.0`, :85).

---

### EQ-1 — Per-T-cell activation memory (costim_induction.py:79)

```
da_i/dt = k_on · p_eng,i · (1 − a_i)  −  k_off · a_i
```
Implemented as **explicit forward Euler in hours** (costim_induction.py:75, :79–80):
```
dt_h  = 24 · dt_days                                          (:75)
a_i  += dt_h · ( k_on · e_i · (1 − a_i) − k_off · a_i )       (:79)
a_i   = clip(a_i, 0, 1)                                       (:80)
e_i   = clip(p_eng,i, 0, 1)                                   (:76)
```

- **Biological meaning:** `a_i` is a **per-T-cell memory of TCR engagement** — the integrated
  transcriptional/translational state that licenses activation-induced costim receptor to appear on the
  surface. It **builds** while the cell is engaged (`p_eng > 0`) and **decays** when it disengages. This is
  a lumped one-state surrogate for {TCR signal → NF-κB/NFAT → transcription → translation → surface
  trafficking}.
- **Mechanistic rationale:** the `(1 − a_i)` factor makes induction **saturating** — a cell cannot
  upregulate past its maximal induced state, which is what makes `FOLD` a *bounded* fold-change rather
  than an unbounded ramp. The `−k_off·a_i` term gives **reversibility**: 4-1BB/OX40 surface expression is
  transient and wanes after antigen withdrawal, so a monotone latch would be wrong. The rejected
  alternative is explicit in the module header (costim_induction.py:4–5): a **static resting density read
  once and never updated**, which under-rates the inducible arms by the whole induction factor.
- **Units:** `a` dimensionless [0,1]; `k_on`, `k_off` hr⁻¹; `p_eng` dimensionless [0,1); `dt_h` hr.
- **Steady state (DERIVED, RUN-VERIFIED):** `a* = k_on·e / (k_on·e + k_off)`. At full engagement
  (`e=1`) with the fallback kinetics: `a* = 1.0/(1.0+0.05) = 0.952381`. Confirmed by integration this task
  (a → 0.9524 by t≈3 h and stationary thereafter).
- **NUMERICAL CAVEAT (DERIVED, RUN-VERIFIED):** forward Euler is **monotone only while
  `dt_h·(k_on·e + k_off) ≤ 1`**. The step actually used by the screen is `dt = 0.02 d → dt_h = 0.48 h`
  (**passed by the runner at run_tce_pd_reval.py:143**, *not* an engine default — the `CoupledPerCellPK`
  class default is `dt=0.01`, coupled_percell_pk.py:65, and `PD_DT` in the environment overrides both,
  coupled_percell_pk.py:74). At `dt_h=0.48` with the fallback `k_on=1.0`, the factor is `0.48·1.05 = 0.504`
  — safe. But `k_on = 3/t_peak_h` (EQ-3), so the monotonicity bound is
  `t_peak_h ≥ 3 / (1/0.48 − 0.05) = 1.475 h`: a **sourced `t_peak_h ≲ 1.48 h` would overshoot at the
  screen's step.** Demonstrated this task: forcing `dt_h=1.0` drives `a` to exactly 1.0 on the first step
  (saved only by the `clip` at :80) before relaxing back to 0.9524. If a fast `t_peak` is sourced, either
  shrink `dt` or move to an implicit/exponential update.

---

### EQ-2 — Activation-induced receptor density (costim_induction.py:81)

```
R_i(t) = R_rest,i · ( 1 + (FOLD − 1) · a_i )
```

- **Biological meaning:** the cell's **surface costim receptor copy number** rises from its resting value
  toward `FOLD × R_rest` as its activation memory fills. This is the equation that makes 4-1BB "appear"
  on an engaged T cell.
- **Mechanistic rationale:** written as a **linear interpolation in `a`** between `R_rest` (a=0) and
  `FOLD·R_rest` (a=1), so that **`FOLD = 1 ⟹ R(t) ≡ R_rest` identically, for any `a`, any history**. That
  algebraic identity is what makes the whole feature **strictly additive**: constitutive arms (CD28/CD2/CD27)
  and the no-costim clinical-validation molecules are **byte-identical** to the pre-fix model
  (costim_induction.py:20–22). RUN-VERIFIED this task: for CD28, CD2 and CD27 with `p_eng=[1.0, 0.5, 0.0]`,
  `R(t) = [5000, 5000, 5000] ≡ R_rest` exactly (`np.allclose → True`).
  A multiplicative form `R = R_rest·FOLD^a` would *also* satisfy the fold=1 identity, but would make the
  *rate* of appearance exponential in `a` rather than proportional to it; the linear form was chosen (no
  in-code rationale given — **SOURCE NOT FOUND IN CODE — TBD** as to why).
- **Preservation of per-cell heterogeneity:** `R_rest,i` is a **per-cell vector** (from the agent table),
  so the induced density inherits the receptor-copy spread across the T-cell population — the induction
  multiplies, it does not homogenize. (That spread is **scRNA-seq-derived, not measured protein**; see the
  `R_costim` row in §3.3 — do not call it "measured".)
- **Units:** `R` copies/cell; `FOLD` dimensionless; `a` dimensionless.
- **⚠ `FOLD` IS NOT SOURCED for any inducible arm (BLOCKER-2).** At this HEAD the four inducible arms
  carry `fold=None` and **raise**.

---

### EQ-3 — Induction on-rate from time-to-peak (costim_induction.py:70)

```
k_on = 3.0 / t_peak_h        if t_peak_h is given
     = 1.0 /hr               otherwise   (FALLBACK — currently ALWAYS taken)
```

- **Biological meaning:** `t_peak_h` is the measured time from TCR engagement to peak surface expression
  of the costim receptor.
- **Mechanistic rationale:** for a first-order approach to steady state, the state reaches
  `1 − e^{−k_on·t}` of its plateau; `k_on·t_peak = 3` gives `1 − e^{−3} = 0.9502`, i.e. **~95% of max at
  `t_peak`** — exactly what the code comment claims (costim_induction.py:68). The mapping is therefore
  **internally consistent [DERIVED]**, and would be a legitimate way to convert a literature
  time-to-peak into a rate constant *if a time-to-peak were sourced*.
- **Units:** `k_on` hr⁻¹; `t_peak_h` hr.
- **⚠ `t_peak_h = None` for EVERY arm in the table** (costim_induction.py:35–42). The `3.0/t_peak_h` branch
  is **dead code at this HEAD**; every arm takes the `1.0 /hr` fallback, which is an
  **[ASSUMED]** round number with no source.

---

### EQ-4 — Induction off-rate from decay time (costim_induction.py:71)

```
k_off = 1.0 / t_decay_h      if t_decay_h is given
      = 0.05 /hr             otherwise   (FALLBACK — currently ALWAYS taken)
```

- **Biological meaning:** the rate at which a disengaged T cell loses its induced costim receptor
  (internalization + shedding + transcriptional shut-off). Sets **how long the tumour-conditional window
  stays open** after a T cell leaves a synapse.
- **Mechanistic rationale:** `k_off = 1/t_decay` is the standard e-folding-time → rate-constant identity
  (`t_decay` = the 1/e time, not the half-life; if a literature value is a *half-life*, the correct
  conversion is `k_off = ln2/t_½`, and plugging a half-life straight into this slot would **over-estimate
  the persistence by 1/ln2 = 1.44×**). Flagged here because the code does not say which convention
  `t_decay_h` expects — **SOURCE NOT FOUND IN CODE — TBD**.
- **Units:** `k_off` hr⁻¹; `t_decay_h` hr.
- **⚠ The `0.05 /hr` fallback implies a 20 h decay time**, i.e. an induced receptor that persists ~20 h
  after disengagement. This is **[ASSUMED]** — no source. It is also **load-bearing**: with `k_on=1.0` it
  sets the plateau `a* = 0.952` (EQ-1), so it directly scales the induced density.

---

### EQ-5 — The engagement drive `p_eng` (wholebody_pd.py:333)

```
p_eng,i = Cb_i / ( Cb_i + R_CD3,i + 1e-12 )
```

- **Biological meaning:** the fraction of that T cell's CD3 that is committed to a **bridged, target-
  engaged** synapse. This is the model's proxy for "has this T cell seen TCR signal?", and it is the
  **sole driver** of EQ-1.
- **Mechanistic rationale:** this is the equation that makes conditionality **emergent**. It is per-cell
  and it is spatial: `Cb_i` is built from *this* T cell's local drug and *its own neighbours'* antigen
  (via `syn_TAA_mean`, wholebody_pd.py:309). A T cell in a healthy organ with no antigen-bearing neighbours
  has `Cb ≈ 0 ⟹ p_eng ≈ 0 ⟹ a → 0 ⟹ R ≡ R_rest`: **it never upregulates, so the costim arm never fires on
  it.** No tumour flag is consulted anywhere.
- **Units:** dimensionless [0,1). `Cb`, `R_CD3` both nM (already converted from copies by `NM_PER_COPY`).
- **⚠ INCONSISTENCY (3-way) — the same name `p_eng` denotes three different formulas:**

  | Site | Formula | Consumer |
  |---|---|---|
  | wholebody_pd.py:333 | `Cb/(Cb + R_CD3)` | **costim induction (EQ-1)** |
  | wholebody_pd.py:379 | `clip(Cb / R_CD3, 0, 1)` | myeloid IL-6 |
  | wholebody_pd.py:485 | `clip(B2 / RC, 0, 1)` | myeloid IL-6 (kinetic path) |

  The first is a **saturating Hill form** (`→ 0.5` when `Cb = R_CD3`); the second is a **ratio clipped at
  1** (`→ 1.0` when `Cb = R_CD3`). They are **not the same function** and differ by 2× at `Cb=R_CD3`.
  Worse, `costim_induction.py:17` **documents the drive as `B2/RC`** — the *third* form — while the QSS
  path actually feeds it the *first*. The code and its own docstring disagree. **No rationale in code for
  the divergence — SOURCE NOT FOUND IN CODE — TBD.** A reviewer will ask which one is intended.

- **⚠ SECOND-ORDER CONSEQUENCE (found in adversarial review): the costim effector gain leaks into IL-6 in
  the QSS path only.** `Cb` is **reassigned in place** at :350 (`Cb = Cb*g_eff`, EQ-13). The myeloid IL-6
  drive at :379 then reads that **effector-gained** `Cb` — so in the QSS path, costim agonism raises myeloid
  IL-6 **twice**: once through `cyto_sig_gain` (EQ-15, intended) and once, silently, through an inflated
  `p_eng` handed to `myeloid.step()` (:380). The kinetic path's IL-6 drive (:485, `B2/RC`) is **not**
  `g_eff`-scaled, so it has no such term. **The two paths therefore disagree on costim's effect on IL-6**,
  on top of the `supp_extra` disagreement in EQ-13. Whether the QSS double-count is intended is
  **SOURCE NOT FOUND IN CODE — TBD**; the comment at :377–378 describes `p_eng` as the plain engaged
  fraction, with no mention that `Cb` has already been multiplied by the costim gain.

---

### EQ-6 — Per-cell costim occupancy (wholebody_pd.py:336)

```
occ_i = [ C_i / (C_i + KD_costim) ] · [ R_costim,i / max(anchor, 1e-9) ]
anchor = mean( R_costim,j  over T cells j with R_costim,j > 0 )     (:200–201, computed ONCE at __init__)
```

- **Biological meaning:** the costim arm's engagement on cell `i` = (Langmuir occupancy of the costim
  paratope at the local free drug) × (how many costim receptors *this* cell has, relative to the
  population). A cell with twice the population-mean receptor density gets twice the drive at the same drug.
- **Mechanistic rationale:** factorizes the drive into a **drug term** (saturable, KD-set) and a
  **receptor-abundance term** (per-cell, from scRNA-seq-derived copies — *derived*, not measured; §3.3).
  The receptor term is what
  **couples EQ-2 back into the signaling**: as `R_costim` induces upward, `occ` rises — the tumour-
  conditional amplification.
- **CRITICAL — the anchor is computed ONCE, from RESTING densities** (`__init__`, :200–201) and is
  **never recomputed**. This is *deliberate and necessary*: it is what lets `occ` exceed 1 as the receptor
  induces. If the anchor were recomputed each step it would track the induced mean and **cancel the entire
  induction effect**. This is a subtle, correct, and completely undocumented-in-code design choice.
- **Units:** `C`, `KD_costim` nM; `R_costim`, `anchor` copies/cell; `occ` dimensionless.

---

### EQ-7 — Occupancy clip (wholebody_pd.py:338)

```
occ_i ← clip(occ_i, 0.0, 5.0)
```

- **Mechanistic rationale:** a numerical guard on a quantity that is a ratio to a population mean and can
  in principle be large for an outlier cell.
- **⚠ THIS CLIP SILENTLY CAPS THE INDUCTION (interaction with EQ-2/EQ-6 — flagged as a design hazard).**
  Because the anchor is the **resting** mean (EQ-6), the induced occupancy scales roughly as
  `(1 + (FOLD−1)·a)` times the resting-normalized occupancy. **The clip at 5.0 therefore truncates any
  `FOLD` above ~5×** for a cell of average resting density at saturating drug. Reported fold-upregulations
  for 4-1BB/OX40 are, as far as this subsystem is concerned, **UNSOURCED** — but **if the sourcing agent
  returns a fold of, say, 10–50×, this clip will quietly eat most of it**, and the model will once again
  under-rate the inducible arms — the very failure mode T7 exists to prevent. **This must be revisited at
  the same time the FOLD values land.** The value 5.0 is **[ASSUMED]** — no source, no comment.
- **Units:** dimensionless.

---

### EQ-8 — CIS/TRANS coincidence gating (wholebody_pd.py:294–295; `p_cis` from multiarm_binding.py:45–46)

```
f_cd3,i  = e_i / ( e_i + max(ref, 1e-12) ),     ref = median( e_j : e_j > 0 )      (:293–294)
occ_i    ← occ_i · [ (1 − p_cis) + p_cis · f_cd3,i ]                               (:295)

p_cis    = exp( −0.5 · ( (span_coeng_T_nm − 12.5) / 8.0 )² )    (multiarm_binding.py:45–46)
         = 0.0   if span_coeng_T_nm is None                     (multiarm_binding.py:44; the DEFAULT)
```
where `e` is the per-T CD3 engagement: the QSS trimer `Cb` (:337) or the kinetic bridged species `B2` (:449).

- **Biological meaning:** can **one molecule** simultaneously engage CD3 and the costim receptor on the
  **same** T-cell surface (**cis**), or must the costim arm act on a **different** cell / independently
  (**trans**)? A cis design makes costim agonism **coincident** with TCR engagement — costim signal only
  where there is already a synapse. A deliberate epitope-**height mismatch** (the code names tall 4-1BB
  CRD1, ~60 Å) forces the co-engagement span away from the ~12.5 nm inter-epitope gap, the arms cannot
  co-reach, and the design becomes trans.
- **Mechanistic rationale:** `p_cis` is a **Gaussian tolerance in span** about a height-matched gap —
  geometry, not a phenomenological constant. The gate is a **convex blend**: `p_cis=0` (trans / no costim /
  the default) returns `occ` **unchanged and byte-identical** (early return at :291), so this is additive;
  `p_cis=1` (perfectly height-matched cis) gates costim **entirely** on the same cell's CD3 engagement.
  Coincidence is thus **emergent from span geometry × real per-cell CD3 binding** (:287–288).
- **⚠ SELF-NORMALIZING REFERENCE:** `ref` is the **median engagement over currently-engaged T cells**,
  recomputed **every step**. So `f_cd3` is a *relative* rank-like quantity, not an absolute occupancy: a
  cell at the population median always gets `f_cd3 = 0.5`, **regardless of the absolute drug level**. This
  makes the cis gate **dose-insensitive in its normalization** — at 1000× the dose, the median cell still
  reads 0.5. Whether that is intended is **SOURCE NOT FOUND IN CODE — TBD**; the comment calls it
  "self-scaling" (:293) without justifying it. A reviewer will attack this.
- **Units:** `span` nm; `gap_match=12.5` nm; `tol=8.0` nm; `p_cis`, `f_cd3` dimensionless.
- **Default is OFF:** `span_coeng_T_nm=None` and `n_costim=0` (wholebody_pd.py:116, :209) ⟹ `p_cis=0` ⟹
  costim drive is **cell-autonomous (trans/legacy)**.

---

### EQ-9 — Signaling integration (wholebody_pd.py:339)

```
PerCellSignaling.step( occ, dt·24.0 )        # signaling kinetics are per-HOUR; dt is days
programs → effector p_E , suppression p_S , exhaustion p_X , cyto_IFNG , cyto_TNF , cyto_IL2   (:340–341, :368)
```

- **⚠ THE MODULE THAT COMPUTES THIS DOES NOT EXIST (BLOCKER-1).** `PerCellSignaling` is constructed at
  :213 with `drive_source="lane"`. Its state equations, its "lane drive table" magnitudes, its
  Rest/8 hr/48 hr "hero kinetics" timescales (all referenced in the comment at :183–187), and its
  per-cell heterogeneity model are **all in `signaling_dynamics.py`, which is not in the tree**.
- **Everything about the programs is therefore SOURCE NOT FOUND IN CODE — TBD.** We can document only the
  *consumption* of `p_E`, `p_S`, `p_X` (EQ-10…EQ-13), which is what the rest of this section does. We
  **cannot** state their units, their bounds, their sign convention, or their provenance. The equations
  below imply `p_E`, `p_X`, `p_S` are dimensionless and roughly O(1), and that positive = more agonism —
  but **that is inferred from the consuming code, not read from the producing code**, and is flagged as
  such.
- The `dt·24.0` conversion is the one thing that *is* verifiable: PD steps in **days**, signaling in
  **hours**.

---

### EQ-10 — Effector program → kill gain (wholebody_pd.py:343–345)

```
kE_gain = 0.55                                              (:343, "locked calib")
g_eff,i = exp( kE_gain · p_E,i )    for CD8/effector T cells only  (:345)
        = 1.0                        for all other T cells (incl. Tregs)  (:344)
```

- **Biological meaning:** costim agonism (4-1BB/OX40/CD28) makes an **effector** T cell a **better killer** —
  more granzyme/perforin, faster serial killing, more sustained synapses.
- **Mechanistic rationale:** **exponential**, so the gain is (i) **strictly positive** for any real `p_E`,
  and (ii) **multiplicative on the trimer** (EQ-13), i.e. costim *scales* the existing antigen-driven kill
  rather than adding a costim-independent killing term. This is the right structure: **costim without a
  TCR signal must not kill.** With `Cb=0`, `g_eff` multiplies zero. Restricting the gain to `_is_cd8_T`
  (:345, mask built at :215) prevents costim from spuriously making Tregs cytotoxic.
- **Units:** `g_eff` dimensionless multiplier; `kE_gain` dimensionless (per unit `p_E`); `p_E` dimensionless.
- **`kE_gain = 0.55` is [FITTED: target not named in code].** The comment says only *"locked calib"*
  (:343). No fit target, no dataset, no residual is stated anywhere in either live file. **A "locked calib"
  with no named target is exactly the fitted-constant-masquerading-as-a-value pattern this project has a
  history of.** It must be treated as **FITTED, provenance unknown**, until the fit is produced.

---

### EQ-11 — Exhaustion program → kill attenuation (wholebody_pd.py:347)

```
g_eff,i ← g_eff,i · exp( −0.30 · max(p_X,i , 0) )
```

- **Biological meaning:** **agonism that drives exhaustion loses durable killing.** This is the *tox/efficacy
  trade* that makes the counterscreen a counterscreen: an arm can win on acute effector gain and still
  lose on exhaustion. Chronic/excessive costim (notably 4-1BB and CD28 superagonism) drives terminal
  differentiation and functional exhaustion.
- **Mechanistic rationale:** exponential **decay**, multiplied onto the effector gain, so the two programs
  **compete on the same axis** — net gain `exp(0.55·p_E − 0.30·p_X)`. An arm with high effector *and* high
  exhaustion drive can net out at ≈1 (no benefit), which is precisely the discriminating behaviour the
  screen needs. `max(p_X, 0)` one-sidedly clamps: a *negative* exhaustion program cannot become a kill
  *bonus*.
- **Units:** dimensionless throughout.
- **The coefficient `0.30` is [UNSOURCED — TBD]** — a bare literal at :347 (and again at :455), with **no
  name, no comment, no source**. It sets the entire effector-vs-exhaustion exchange rate (0.55 vs 0.30 ⟹
  exhaustion is weighted 0.55× as strongly as effector gain). This single unsourced number decides whether
  4-1BB's known exhaustion liability is enough to beat its effector benefit — i.e. **it can flip the
  ranking the screen exists to produce.**

---

### EQ-12 — Suppression program on Tregs (wholebody_pd.py:349)

```
supp_extra = mean( max( p_S,j , 0 ) )   over T cells j that are Tregs        (:349)
           = 0.0                        if the organ has no Tregs
```

- **Biological meaning:** costim receptors are **not T-effector-exclusive** — CD28, GITR, OX40 and 4-1BB
  are all expressed on Tregs, and agonizing them can **expand/activate Tregs** and thereby *suppress* the
  antitumour response. This is a genuine, clinically-observed failure mode (the GITR/OX40 Treg paradox),
  and it is the mechanism by which a costim arm can be **net-negative**.
- **Mechanistic rationale:** the Treg suppression program is read from **only the Treg subset**
  (`_is_treg_T`, mask at :214) and **averaged into a single scalar per organ**, which then **amplifies the
  existing spatial Treg damping** (EQ-13).
- **⚠ MEAN-FIELD COLLAPSE:** the per-cell suppression program is reduced to **one scalar for the whole
  organ**, destroying the spatial structure that the rest of the module works so hard to preserve. The
  existing Treg term (`n_treg`) *is* spatial (a per-T-cell count within `R_TREG_UM=50 µm`,
  :275–277), but `supp_extra` is **global**: a Treg on the far side of the organ suppresses a synapse it
  can never reach. This is an **imposed, non-emergent shortcut** in an otherwise per-cell subsystem.
- **Units:** dimensionless.

---

### EQ-13 — Effector-gained trimer and suppressed kill (wholebody_pd.py:350–351)

```
Cb_i    ← Cb_i · g_eff,i                                                  (:350)
kill_T,i = Cb_i / ( 1 + TREG_K · n_treg,i · (1 + supp_extra) )            (:351)   TREG_K = 0.25 (:72)
```

- **Biological meaning:** the killing capacity of T cell `i` = its engaged trimer, **amplified** by its
  costim effector program, **damped** by the Tregs in its 50 µm neighbourhood, and damped **further** if
  the costim arm is also driving those Tregs' suppression program.
- **Mechanistic rationale:** the costim suppression enters as a **multiplier on the per-Treg suppression
  constant**, `TREG_K·n_treg·(1+supp_extra)`, not as an additive kill penalty. This is the right structure:
  **an organ with no Tregs (`n_treg=0`) is immune to costim-driven suppression** no matter how high
  `supp_extra` is — you cannot be suppressed by Tregs that are not there. That property is emergent from
  the product form.
- **Units:** `Cb`, `kill_T` nM; `n_treg` dimensionless count; `TREG_K` per-Treg (dimensionless);
  `supp_extra` dimensionless.
- **⚠ Note the asymmetry:** `Cb` is **reassigned** at :350, so the cytokine `eng` computed at :364 uses
  the **effector-gained** `Cb` — but **without** the `(1+supp_extra)` term that the kill gets
  (:364 vs :351). So in the QSS path **costim-driven Treg suppression damps KILLING but NOT CYTOKINE.**
  In the kinetic path it damps **both** (`treg_damp` at :464 includes `supp_extra` and is applied to
  `eng` at :471). **The two paths are mechanistically inconsistent.** No rationale in code —
  **SOURCE NOT FOUND IN CODE — TBD.** Biologically, Treg suppression of cytokine release is well-described,
  so the QSS path is the suspect one.

---

### EQ-14 — Hazard accumulation and kill fraction (wholebody_pd.py:354, :362, :385)

```
dkill_target = Wt_normᵀ · kill_T                                  (:354)   antigen-weighted apportionment
kill_hazard += dt · k_death · dkill                               (:362)
kill_frac    = mean( 1 − exp(−kill_hazard) )   over target cells  (:385)
```

- **Biological meaning:** each T cell's killing is **apportioned** among the target cells in its 30 µm
  synapse reach (`R_SYN_UM`, :36), weighted by each neighbour's antigen density and row-normalized
  (`Wt_norm`, :267–269). Cumulative hazard → survival `exp(−H)` → killed fraction.
- **Mechanistic rationale:** hazard is integrated as a **drug-graded RATE** (`dt·k_death·dkill`), not a
  saturating fixed hazard, so **low exposure ⟹ low rate ⟹ incomplete kill within the dosing window**
  (:356–361). This is what lets the model distinguish exposure-limited from receptor-limited arms.
- **T7's entry point** is solely through `kill_T` (EQ-13): costim modulates kill **only** by scaling the
  trimer and the Treg damping. It introduces **no** costim-specific death pathway.
- **Units:** `kill_hazard` dimensionless; `k_death` day⁻¹; `dt` day; `kill_frac` ∈ [0,1).

---

### EQ-15 — Costim → cytokine gain (wholebody_pd.py:368–374)

```
cs            = mean( 0.45·p_IFNG + 0.32·p_TNF + 0.18·p_IL2 )        (:370)
cyto_sig_gain = max( 0.2 , 1.0 + cs )                                (:371)
resp          = max( 1 − Dcyto , 0 )                                 (:373)
cyto_rate[k]  = CYTO_HIER[k] · eng · resp · cyto_sig_gain            (:374)
eng           = Σ_i  Cb_i / (1 + TREG_K · n_treg,i)                  (:364)
```

- **Biological meaning:** costim agonism **raises the cytokine storm** (this is the CRS liability of the
  costim arm — the tox half of the counterscreen), while a "cold" arm **lowers** it. `cyto_sig_gain > 1`
  ⟹ more CRS than the plain TCE at the same engagement.
- **Mechanistic rationale:** a **bounded** linear read-out of the three cytokine programs. The floor at
  0.2 (:371) prevents a strongly-negative program from driving cytokine to zero (a cold costim arm
  suppresses, but cannot abolish, TCE-driven CRS — correct: the CD3 arm alone still causes CRS).
- **⚠ The weights `0.45 / 0.32 / 0.18` are [UNSOURCED — TBD]** — bare literals at :370 (and again at
  :475), no comment, no source. Note they are **NOT** the `CYTO_HIER` weights (`IFN=0.36, TNF=0.31,
  IL2=0.18`, :22) — they are a *different, similar-looking* set applied at a *different* stage. Two
  unexplained near-duplicate weight vectors in one file is a provenance smell a reviewer will pull on.
- **Units:** all dimensionless; `cyto_rate` in raw engagement-sum units (converted to pg/mL only via
  `cytokine_to_pgml`, :32–35).

---

### EQ-16 — LEGACY scalar costim path (wholebody_pd.py:319–320)

```
g_eff = 1.0 + costim_boost        ONLY IF (self.sig is None) AND (costim_boost > 0)      (:319–320)
```

- **Mechanistic rationale:** the **pre-T7** costim model: a single organ-wide scalar multiplier, with no
  per-cell receptor density, no occupancy, no programs, no induction, and no conditionality. Retained
  purely for backward compatibility and **mutually exclusive** with the per-cell path (`self.sig is None`
  guard). `costim_boost = 0.0` is the default (:111), so this is **inert** unless explicitly set.
- **This is the thing T7 replaces**, and it is documented here only so the contrast is on the record.

---

### EQ-17 — Receptor write-back and the aliasing guard (wholebody_pd.py:198, :334)

```
self.R_costim   = Rc.copy()                                (:198)   ← the .copy() is LOAD-BEARING
self.R_costim[Tidx] = self.costim_ind.step(dt, p_eng)      (:334)   ← induced density written back
```

- **Mechanistic rationale:** `CostimInduction.step()` returns the induced density **over the T-cell subset
  only**, which is written back into the organ's `R_costim` at the T-cell indices. It is read on the *next*
  step by EQ-6 — closing the induction↔occupancy feedback loop.
- **The `.copy()` at :198 is a real bug-fix, not boilerplate** (:195–197): `np.asarray` does **not** copy
  an existing float array. Without the copy, **every organ would share and mutate the caller's array** and
  leak its induced densities into the other organs — a cross-organ contamination that would be almost
  impossible to detect in output. Worth keeping this comment.

---

## 3. PARAMETERS OWNED

### 3.1 `costim_induction.py` — the induction table (COSTIM_INDUCTION, :33–43)

| Symbol | Value | Units | Provenance tag | Source (as claimed in code) | Mechanistic rationale |
|---|---|---|---|---|---|
| `CD28.fold` | `1.0` | — | **[ASSUMED: constitutive]** | `source="TBD"` (:35) | CD28 is constitutive on resting T cells ⟹ the resting scRNA-seq density IS correct ⟹ fold=1 makes EQ-2 an identity. The *claim* "constitutive" is biologically standard but **carries no citation in code**. |
| `CD2.fold` | `1.0` | — | **[ASSUMED: constitutive]** | `source="TBD"` (:36) | as above |
| `CD27.fold` | `1.0` | — | **[ASSUMED: constitutive]** | `source="TBD"` (:37) | as above |
| `TNFRSF9.fold` (4-1BB) | **`None`** | — | **[UNSOURCED — TBD]** | `source="PENDING"` (:39) | **HARD BLOCKER.** Absent on resting T cells; induced post-TCR. Value required. |
| `TNFRSF4.fold` (OX40) | **`None`** | — | **[UNSOURCED — TBD]** | `source="PENDING"` (:40) | **HARD BLOCKER.** |
| `ICOS.fold` | **`None`** | — | **[UNSOURCED — TBD]** | `source="PENDING"` (:41) | **HARD BLOCKER.** |
| `TNFRSF18.fold` (GITR) | **`None`** | — | **[UNSOURCED — TBD]** | `source="PENDING"` (:42) | **HARD BLOCKER.** |
| `t_peak_h` (**all 7 arms**) | **`None`** | hr | **[UNSOURCED — TBD]** | — (:35–42) | Would set `k_on = 3/t_peak` (EQ-3). Never populated ⟹ EQ-3's sourced branch is **dead code**. |
| `t_decay_h` (**all 7 arms**) | **`None`** | hr | **[UNSOURCED — TBD]** | — (:35–42) | Would set `k_off = 1/t_decay` (EQ-4). Never populated ⟹ **dead code**. |
| **CD40, TNFRSF25, CD226, TNFRSF14** — **screened but ABSENT from the table** | **`1.0` (implicit)** | — | **[ASSUMED — NOT EVEN DECLARED]** | `source="NOT IN TABLE"`, set at runtime (:54) | **BLOCKER-4.** These four arms *are* in the 11-arm screen (run_tce_pd_reval.py:150–151) but have no table row. They take the `p is None` branch (:53–57) ⟹ `fold=1.0`, `kind="unknown"`, **stderr warning, NO raise, `strict` not consulted**. The model silently asserts they are constitutive. RUN-VERIFIED. |

**Every `source` field in the live table is `"TBD"` or `"PENDING"`. There is not one PMID in this module.**
The header (`:23–24`) *claims* "fold-upregulation + induction kinetics per arm are LITERATURE values" —
**that claim is false at this HEAD** and is contradicted 8 lines later by its own `*** PENDING SOURCING ***`
banner (:31). Do not propagate the header claim.

### 3.2 `costim_induction.py` — kinetics and internals

| Symbol | Value | Units | Provenance tag | Source | Mechanistic rationale |
|---|---|---|---|---|---|
| `k_on` (fallback) | `1.0` | hr⁻¹ | **[ASSUMED]** | none (:70) | Taken by **every arm** at this HEAD. Round number. Implies ~3 h to plateau. |
| `k_off` (fallback) | `0.05` | hr⁻¹ | **[ASSUMED]** | none (:71) | Taken by **every arm**. Implies a 20 h decay time. **Load-bearing**: with `k_on=1.0` it sets `a* = 0.952`. |
| `k_on` (sourced form) | `3.0 / t_peak_h` | hr⁻¹ | **[DERIVED: 1−e⁻³ = 95% of plateau at t_peak]** | internal (:68, :70) | Self-consistent conversion; **currently unreachable** (no `t_peak_h`). |
| `k_off` (sourced form) | `1.0 / t_decay_h` | hr⁻¹ | **[DERIVED: e-folding time → rate]** | internal (:71) | **Ambiguity:** does not state whether `t_decay_h` is a 1/e time or a half-life (see EQ-4). |
| `a(0)` | `0` (all cells) | — | **[ASSUMED]** | (:51) | T cells start naive/unengaged at dose. Reasonable for a first-dose sim; **wrong for a re-dose** into an already-activated repertoire (see §5). |
| `strict` | `True` | bool | **[CODE-INTERNAL]** | (:48; call site :222) | Fail-closed **only for arms present in the table**: the check sits inside the `else:` branch (:60–67), so an arm that is *not in the table at all* never reaches it and runs at fold=1.0 regardless (**BLOCKER-4**, :53–57). For table arms the behaviour is correct — **do not "fix" it by flipping it to False.** |
| `frac_induced` threshold | `0.1` | — | **[ASSUMED]** | (:86) | Reporting only (`summary()`); does not affect dynamics. |

### 3.3 `wholebody_pd.py` — costim-arm parameters

| Symbol | Value | Units | Provenance tag | Source | Mechanistic rationale |
|---|---|---|---|---|---|
| `KD_costim_nM` | `1.0` | nM | **[UNSOURCED — TBD]** | none (:112; default at run_tce_pd_reval.py:153) | Sets the drug-occupancy half-point of EQ-6. A **round default for all 11 arms**, i.e. the model currently assumes **every costim binder has identical affinity** — which erases a first-order design variable of the screen. |
| `kE_gain` | `0.55` | — | **[FITTED: target NOT NAMED in code]** | comment says only *"locked calib"* (:343, :453) | Effector→kill sensitivity. **No fit target, dataset, or residual is stated anywhere.** Treat as fitted-with-unknown-provenance. |
| exhaustion coefficient | `0.30` | — | **[UNSOURCED — TBD]** | **bare literal, no comment** (:347, :455) | Sets the effector-vs-exhaustion exchange rate. **Can flip the screen's ranking.** |
| cytokine program weights | `0.45 / 0.32 / 0.18` (IFNG/TNF/IL2) | — | **[UNSOURCED — TBD]** | **bare literals, no comment** (:370, :475) | Not the same as `CYTO_HIER`. Two near-duplicate weight vectors, one file. |
| `cyto_sig_gain` floor | `0.2` | — | **[ASSUMED]** | (:371, :475) | A cold arm suppresses but cannot abolish TCE-driven CRS. |
| occupancy clip max | `5.0` | — | **[ASSUMED]** | **bare literal, no comment** (:338, :450) | **Silently caps FOLD at ~5× (EQ-7). Revisit when FOLD lands.** |
| `costim_boost` | `0.0` | — | **[CODE-INTERNAL: legacy, inert]** | (:111) | Pre-T7 scalar path. Mutually exclusive with the per-cell path. |
| `n_costim` | `0` | count | **[CODE-INTERNAL default]** | (:116, :205) | ⟹ `p_cis = 0` ⟹ trans/legacy. |
| `span_coeng_T_nm` | `None` | nm | **[CODE-INTERNAL default]** | (:116, :206) | ⟹ `p_cis = 0`. Cis geometry **OFF by default**. |
| `gap_match_nm` | `12.5` | nm | **[UNVERIFIED CITATION]** | multiarm_binding.py:39 — no PMID given | Height-matched inter-epitope gap for cis feasibility. Dependency module, not owned here. |
| `tol_nm` | `8.0` | nm | **[ASSUMED]** | multiarm_binding.py:39 | Gaussian span tolerance. |
| `_costim_anchor` | `mean(R_costim[T] > 0)` | copies/cell | **[DERIVED: from per-cell data, at __init__]** | (:200–201) | **Resting** population mean. Frozen by design (EQ-6). |
| `R_costim` | per-cell, from `'<ARM>_copies'` | copies/cell | **[DERIVED — external RNA→protein pipeline; UNVERIFIED]** — **NOT [MEASURED]** | (coupled_percell_pd.py:26–27, :40) | These are **scRNA-seq transcript counts converted to protein copies**, not measured surface densities. The conversion factor, its calibration, and its per-gene validity are **outside both live files and were not audited here**. The only in-code provenance is the bare phrase "HPA/Glassman" (wholebody_pd.py:182) — **[UNVERIFIED CITATION]**, no PMID/DOI anywhere in either live file. Anyone tempted to write "measured receptor copies" in a figure legend: **you cannot support that from this code.** For an activation-induced receptor this is doubly fraught — resting transcript for TNFRSF9/TNFRSF4 sits at the detection floor, so `R_rest` is a **noise-dominated estimate of a near-zero quantity** (see §5 Q4). |

### 3.4 Co-resident constants in `wholebody_pd.py` that T7 multiplies into (NOT owned by T7)

Listed because T7's outputs pass through them; each is owned by the kill-core / cytokine doc.

| Symbol | Value | Units | Provenance tag | Note |
|---|---|---|---|---|
| `TREG_K` | `0.25` | — | **[FITTED: "validated tumor constant"]** (:72) | Amplified by `supp_extra` in EQ-13. |
| `R_SYN_UM` | `30.0` | µm | **[ASSUMED: "validated tumor value"]** (:36) | Synapse reach. |
| `R_TREG_UM` | `50.0` | µm | **[ASSUMED: "validated tumor value"]** (:71) | Treg suppression neighbourhood. |
| `NM_PER_COPY` | `6.0/257000 = 2.3346e-5` | nM/copy | **[FITTED: pinned by the tumor kill anchor]** (:83) | Comment (:76–82) is explicit that it is *"pinned by the validated tumor kill anchor: Rcap_TAA=6.0 nM at CEACAM5 257,000 copies/cell"* ⟹ a **fit**, presented as "ONE physical constant" (a ~71 pL synapse volume). The synapse-volume framing is a *reinterpretation* of a fitted ratio, not an independent measurement. |
| `KD_CD3_nM` | `40.0` | nM | **[UNVERIFIED CITATION]** (:111) | Default; overridden per-molecule from `eng_params`. |
| `KD_TAA_nM` | `1.45` | nM | **[UNVERIFIED CITATION]** (:111) | Default; overridden per-molecule. |
| `CYTO_HIER` | `IL6 1.0 / IFN 0.36 / TNF 0.31 / IL2 0.18` | — | **[UNSOURCED — TBD]** (:22) | Comment claims *"project-standard CYTO_HIER, mosunetuzumab-anchored"* — **no PMID. [UNVERIFIED CITATION].** |
| `K_CYTO_DESENS` | `30.0` | day⁻¹ | **[FITTED]** (:29) | Cytokine tachyphylaxis; comment says "ported from qsp_costim_window_v2 cytokine calib" ⟹ fitted. |
| `K_CYTO_RECOV` | `0.003` | day⁻¹ | **[FITTED]** (:30) | as above |
| `DESENS_ENG_REF` | `1.0e4` | raw eng-sum | **[FITTED/ASSUMED]** (:31) | "~mosun 1mg peak" — a normalizer chosen to make the build O(1). |
| `RCAP_CD3` / `RCAP_TAA` | `2.0` / `6.0` | nM | **[DEAD CODE]** (:73–74) | **Defined but never used anywhere** (grep-verified: only their definitions). The comment at :310–312 refers to `self.R_CD3_cap` / `self.R_TAA_cap`, **attributes that do not exist**. Stale comment + dead constants. |

### 3.5 ⚠ KNOWN-CONTAMINATED VALUE PRESENT IN THIS FILE

| Symbol | Value | Provenance tag | Action |
|---|---|---|---|
| `CYTO_IL6_CLINICAL_ANCHOR_PGML` | `570.0` (wholebody_pd.py:23) | **[UNSOURCED — TBD — KNOWN CONTAMINATION]** | The in-code comment claims *"mosunetuzumab peak IL-6 ~570 pg/mL (Hosseini 2020 Fig5A)"*. **Per the provenance audit of 2026-07-13, the value 570 has NO SOURCE.** The only valid IL-6 clinical anchors are **mosunetuzumab 152** and **teclistamab 21**, both **population means**. **DO NOT propagate 570.** Mitigating fact established this task: the constant is **never read** — grep across `engine/` and `rundir/handoff/` finds only its definition; `cytokine_to_pgml` (:32–35) takes `il6_eng_scale` as an *argument*. It is **dead but dangerous**: it will be picked up by the next person who needs an anchor. **Recommend deletion.** |

Also in this file: the comment at **wholebody_pd.py:484** asserts organ IL-6 is summed *"against the
**measured** IL-6 clearance"* and that *"**NOTHING here is fitted**"*. Per the same audit, the IL-6
clearance used in this project (0.20/hr, cited to PMID 31268236 / Chen 2019, *Clin Transl Sci*) is from a
**semi-mechanistic MODELING paper that reports no measured IL-6 clearance**; human IL-6 clearance appears
to be **unmeasured in the literature**. The clearance parameter itself lives in another module (out of
scope here), but **the claim made in this file is false** and should be corrected in place.

---

## 4. WHAT IS EMERGENT vs IMPOSED

### Genuinely EMERGENT (computed from mechanism, per cell)

| Quantity | Emerges from |
|---|---|
| **Tumour-conditionality of the inducible arms** | `p_eng,i` (EQ-5) is built from *this* T cell's local drug and *its own* neighbours' antigen. A T cell with no antigen-bearing neighbours has `Cb≈0 ⟹ a→0 ⟹ R≡R_rest`: **it never upregulates.** **No tumour flag, no organ label, no spatial mask is consulted anywhere.** This is the central emergence claim of T7 and it is real. |
| **Which T cells carry costim receptor, and how much** | Per-cell `a_i` (EQ-1) × per-cell `R_rest,i`. The induced population is a spatially-structured subset, not a compartment average. |
| **The induction↔occupancy feedback** | Induced `R_costim` (EQ-2) raises `occ` (EQ-6) on the **next** step, which raises the programs, which raise kill. The amplification is a closed loop, not a lookup. |
| **Prozone / hook in the costim response** | Inherited from `Cb` (EQ-0): at high drug the trimer collapses, `p_eng` falls, induction *relaxes*. Costim inherits the bell shape without being told about it. |
| **Effector-vs-exhaustion net benefit** | `exp(0.55·p_E − 0.30·p_X)` (EQ-10/11) — an arm can win on effector and still lose on exhaustion. |
| **Immunity to Treg suppression in Treg-free organs** | The product form of EQ-13: `n_treg=0 ⟹ supp_extra cannot act`. |
| **Cis vs trans coincidence** | `p_cis` from **span geometry** (EQ-8) × real per-cell CD3 binding — not a phenomenological switch. |

### IMPOSED (handed to the model as a constant)

| Quantity | Why it's imposed |
|---|---|
| **`FOLD`** | A **per-arm constant**, not computed. Even fully sourced, the *magnitude* of induction is a literature input. (Currently `None` ⟹ nothing runs.) |
| **`k_on`, `k_off` of induction** | Constants (currently the `1.0` / `0.05` fallbacks). The **timescale** of conditionality is imposed, only its **spatial pattern** is emergent. |
| **All six signaling programs** | Produced by `PerCellSignaling` from a **"lane drive table"** — i.e. a **per-arm lookup of program magnitudes** (:184–186). The programs are *not* derived from receptor biochemistry; they are tabulated per arm and then given per-cell dynamics. **This is the single biggest limit on the "emergence" claim** — and the table is **not in the tree** (BLOCKER-1), so its degree of imposition cannot even be audited. |
| **`kE_gain=0.55`, exhaustion `0.30`, cyto weights `0.45/0.32/0.18`** | Pure constants mapping programs → phenotype. Unsourced. |
| **`KD_costim = 1.0 nM` for every arm** | Imposed and **identical across all 11 arms** — erases affinity as a design variable. |
| **`supp_extra`** | **Mean-field**: one scalar per organ (EQ-12), applied to every synapse regardless of distance to any Treg. An imposed shortcut inside a per-cell model. |
| **`R_rest`** | An input from the scRNA-seq/HPA pipeline. Not computed here (correctly). |

### Honest summary of where emergence stops

T7's emergence claim is **narrow but real**: *the spatial/temporal pattern of who upregulates costim
receptor, and therefore where the costim arm fires, is genuinely computed from per-cell engagement.* That
is the claim worth defending, and it is defensible.

Everything about **how strongly** a given receptor signals once occupied — the program magnitudes, the
kill/cytokine coefficients, the fold-change, the affinity — is **tabulated per arm**. T7 is therefore best
described as **an emergent spatial gate on an imposed per-arm response**, and the doc should say so rather
than claim the response itself emerges. Overstating this in front of a committee would be the easiest
attack to land.

---

## 5. KNOWN LIMITATIONS / OPEN QUESTIONS

### The four blockers (restated — these are what a reviewer hits first)

1. **`signaling_dynamics.py` is missing ⟹ the subsystem cannot run at all.** Every costim arm raises
   `ModuleNotFoundError` at `OrganPD.__init__` (wholebody_pd.py:192). RUN-VERIFIED. Until this module
   lands, **there are no costim results, only costim code.** It is also, right now, **unaudited science**:
   the program magnitudes and timescales that decide the entire screen live in a file nobody can read.
2. **`FOLD` unsourced for all four inducible arms ⟹ they fail closed (`ValueError`).** RUN-VERIFIED. The
   fail-closed design is *correct*; the gap is *open*. **This is an OPEN BLOCKER, not a delivered feature.**
   Any statement of the form "the model now handles activation-induced costim" must be qualified with
   "…once the folds are sourced."
3. **Induction is silently inert on the kinetic path** (`_step_kinetic` never calls `costim_ind.step()`).
   Under `PD_KINETICS=True` the receptor stays at the resting density — **the original bug, back, silently.**
   This one is a genuine defect, not a gap.
4. **4 of the 11 screened arms (CD40, TNFRSF25, CD226, TNFRSF14) are not in `COSTIM_INDUCTION` and fail
   OPEN at fold=1.0.** RUN-VERIFIED. `strict=True` does not reach them. The model asserts, without data,
   that they are constitutive. **Fails open ⟹ worse than BLOCKER-2, which at least fails closed.**

### Attacks a reviewer will land

- **"Your induction fixes an under-rating bug, but your `clip(occ, 0, 5)` re-introduces it."** (EQ-7.)
  Because the occupancy anchor is the *resting* mean, a `FOLD` above ~5 is truncated by the clip. If the
  literature fold for 4-1BB is large, **the fix is capped at 5× and the inducible arms are still
  under-rated.** This must be resolved *simultaneously* with the FOLD sourcing, or the fix is cosmetic.

- **"`p_eng` means three different things in one file."** (EQ-5.) `Cb/(Cb+R_CD3)` drives induction;
  `clip(Cb/R_CD3,0,1)` drives myeloid IL-6; `clip(B2/RC,0,1)` drives it on the kinetic path. And
  `costim_induction.py:17` documents the drive as the *third* form while receiving the *first*. Pick one,
  justify it, and make the docstring match.

- **"Costim suppression damps kill but not cytokine — in one path only."** (EQ-13.) QSS omits
  `supp_extra` from `eng` (:364); kinetic includes it (:471). The two paths give different CRS for the
  same arm. Biologically, Treg suppression of cytokine release is well-described, so the QSS path is
  likely wrong.

- **"`kE_gain = 0.55` is a 'locked calib' with no named target."** (EQ-10.) Neither live file states what
  it was fit to. Same for `0.30` and `0.45/0.32/0.18`, which are bare literals with **no comment at all**.
  These four numbers jointly determine the efficacy/tox trade that *is* the counterscreen's output. **An
  unsourced constant that can flip the ranking is worse than a missing one**, because it produces a
  confident answer.

- **"Every costim arm has KD = 1.0 nM."** (§3.3.) Affinity is a primary lever in costim engager design
  (too-tight costim binding drives peripheral toxicity; too-weak loses the conditional window). A single
  default across 11 arms means the screen currently cannot see affinity at all.

- **"Your fail-closed guard has a hole you can drive four arms through."** (BLOCKER-4.) The screen has 11
  arms; the induction table has 7. CD40, TNFRSF25, CD226 and TNFRSF14 run at an **undeclared, unsourced
  fold=1.0** with a stderr warning that no one reads in a batch run. The module's header claims such arms are
  *"FLAGGED, never guessed"* — **fold=1.0 IS the guess.** This is the single easiest defect to demonstrate:
  one `CostimInduction('CD226', …, strict=True)` call returns rather than raises.

- **"`supp_extra` is mean-field inside a per-cell model."** (EQ-12.) Justify or make it spatial — the
  `n_treg` neighbourhood machinery to do so already exists (:275–277).

- **"`a(0) = 0` assumes a naive repertoire."** (§3.2.) Correct for first dose. For **cycle 2+** — the
  regime where costim engagers actually differentiate themselves, and where the step-up dosing lives —
  T cells arrive **pre-activated**, so `a(0)=0` **under-rates the inducible arms on every dose after the
  first.** The activation state is not carried across doses.

- **"The `t_peak`/`t_decay` branches are dead code."** (EQ-3/EQ-4.) The stated design ("kinetics from
  literature") is not what runs; every arm takes the `1.0`/`0.05` fallbacks. And EQ-4 does not specify
  whether `t_decay_h` is a 1/e time or a half-life — a **1.44× error** waiting to be made when the values land.

- **Numerical:** forward Euler on EQ-1 is monotone only for `dt_h·(k_on + k_off) ≤ 1`. At the step the
  screen actually runs (`dt = 0.02 d`, passed at run_tce_pd_reval.py:143 — **not** the class default, which
  is 0.01) this permits only `t_peak_h ≥ 3/(1/0.48 − 0.05) = 1.475 h`. A faster sourced `t_peak`
  **will overshoot** (RUN-DEMONSTRATED at `dt_h=1.0`: `a` hits exactly 1.0 on step 1, rescued only by the
  clip). Check this when the kinetics land — and note `PD_DT` in the environment can make the step *larger*
  (coupled_percell_pk.py:74), which widens the unsafe region.

- **Data-integrity:** `CYTO_IL6_CLINICAL_ANCHOR_PGML = 570.0` (:23) is a **known-contaminated value with no
  source** sitting in this file with a confident-looking citation attached. It is currently unread, which
  makes it a **trap**, not a bug. Delete it. Likewise the false *"measured IL-6 clearance… NOTHING here is
  fitted"* claim at :484.

### Open questions that need Max / the sourcing agents

1. **What are the fold-upregulations** for TNFRSF9, TNFRSF4, ICOS, TNFRSF18, with PMIDs? (BLOCKER-2.)
   And on what cells / what stimulus / what timepoint — a fold measured by flow on PBMC anti-CD3/CD28
   blasts is not the same quantity as a fold in a tumour-infiltrating T cell.
2. **What are `t_peak_h` and `t_decay_h`** for each, and **is `t_decay_h` a 1/e time or a half-life?**
3. **Where is `signaling_dynamics.py`,** and what is in the "lane drive table"? (BLOCKER-1.) Without it,
   §2 EQ-9…EQ-15 describe the *consumption* of numbers whose *production* is unauditable.
4. **Is the resting scRNA-seq density even the right `R_rest`** for a receptor that is ~absent at rest? If
   the measured resting copies are at the detection floor, then `R(t) = R_rest·(1+(FOLD−1)a)` **amplifies
   noise**, and the induced density inherits the floor's error. A receptor that is genuinely *zero* at rest
   stays *zero* forever under EQ-2 — **multiplying zero by any fold gives zero.** This may be the deepest
   problem in the subsystem: **the multiplicative form cannot induce a receptor from true zero.** An
   additive form (`R = R_rest + R_induced_max·a`) would not have this failure mode. **This needs a
   decision before the folds land.**
5. **Should the constitutive arms really be fold=1.0?** CD28 is constitutive but is *modulated* on
   activation; CD27 is downregulated on terminal effectors. "Constitutive" is doing real work here and is
   uncited (`source="TBD"`).
6. **Does the `clip(occ,0,5)` cap need to move** when the folds land? (EQ-7.)
7. **What are the induction status and fold for CD40, TNFRSF25 (DR3), CD226 (DNAM-1) and TNFRSF14 (HVEM)?**
   (BLOCKER-4.) They are being screened *right now* at an undeclared fold=1.0. Each needs an explicit
   `COSTIM_INDUCTION` row — `kind="constitutive", fold=1.0` **with a citation**, or `kind="inducible",
   fold=<sourced>` — and the `p is None` branch must be made to `raise` under `strict=True` so this class of
   hole cannot reopen when a 12th arm is added.

---

## 6. VERIFICATION RECORD (this task)

All executed with `/usr/local/Caskroom/miniconda/base/envs/claude-skills/bin/python`, `sys.path` = `engine/`.

| Check | Result |
|---|---|
| `import signaling_dynamics` (engine/ + rundir/handoff/ on path) | **ModuleNotFoundError** |
| `find` for `signaling_dynamics*` over the whole project tree | **zero hits** |
| `OrganPD(...)` **without** costim arm, 1 PD step | **OK**, `kill_frac = 0.0028471`, `sig=None`, `costim_ind=None` |
| `OrganPD(..., costim_arm='TNFRSF9', ...)` | **ModuleNotFoundError: signaling_dynamics** |
| `OrganPD(..., costim_arm='CD28', ...)` | **ModuleNotFoundError: signaling_dynamics** |
| `CostimInduction('TNFRSF9'/'TNFRSF4'/'ICOS'/'TNFRSF18', strict=True)` | **ValueError** (all four) |
| `CostimInduction('CD28'/'CD2'/'CD27').step(0.02, [1.0,0.5,0.0])` | `R(t) = [5000,5000,5000] ≡ R_rest`, `np.allclose → True` (**fold=1 identity CONFIRMED**) |
| `CostimInduction('CD40')` (not in table) | falls back `fold=1.0, kind='unknown'`, warns on stderr |
| `k_on`, `k_off` actually taken by every arm | `1.0 /hr`, `0.05 /hr` (fallbacks) |
| steady state `a*` at `p_eng=1` | `1.0/(1.0+0.05) = 0.952381` — integration converges to `0.9524` |
| forward-Euler overshoot at `dt_h = 1.0` | `a → 1.0000` on step 1, relaxes to `0.9524` (**clip-rescued**) |
| `costim_ind.step()` call sites in `wholebody_pd.py` | **exactly one: line 334 (QSS only)** — none in `_step_kinetic` |
| `CYTO_IL6_CLINICAL_ANCHOR_PGML` consumers | **none** (definition only) — dead constant |
| `RCAP_CD3` / `RCAP_TAA` consumers | **none** (definitions only) — dead constants |

### 6.1 ADVERSARIAL RE-VERIFICATION (second pass, independent of the pass above)

Every citation, line number and parameter value in this doc was re-checked against the live tree. Results:

| Check | Result |
|---|---|
| `PMID 31268236` (§3.5) — is it fabricated? | **NO — real in-code citation.** Present at `engine/coupled_percell_pd.py:277`, `engine/cytokine_pbpk.py:75`, `engine/il6_pbpk.py:50`, and audited in `docs/PROVENANCE_AND_VALIDATION.md:100–105`. The doc's characterization (modeling paper, no measured clearance) matches the repo's own audit. |
| IL-6 anchors 152 (mosun) / 21 (teclistamab), §3.5 | **CONFIRMED in code:** `engine/run_tce_pd_reval.py:74, :77, :81` — both tagged "population MEAN"; teclistamab 21 carries PMID 38831634. |
| "the code names tall 4-1BB CRD1, ~60 Å" (EQ-8) | **CONFIRMED:** `multiarm_binding.py:42` and `:116`; also `wholebody_pd.py:203`. |
| `_cis_feasibility` formula + `gap_match=12.5`, `tol=8.0` (EQ-8) | **CONFIRMED:** `multiarm_binding.py:39` (signature), `:44` (None→0.0), `:45–46` (Gaussian). |
| "11 costim arms" (§3.3, §4) | **CONFIRMED:** enumerated in `run_tce_pd_reval.py:150–151`. |
| Runner threading quote (§1) | **CONFIRMED verbatim:** `run_tce_pd_reval.py:146–149`; threaded at `:152–154`. |
| `attach_pd` has always accepted `costim_arm` (§1) | **CONFIRMED:** `coupled_percell_pd.py:17–18`; `'<ARM>_copies'` at `:26–27`; load at `:38–43`. |
| Any DEAD module documented (`cytokine_pbpk`, `il6_pbpk`, `unified_binding`, `multiarm_kinetic`, `biexact_solver`, `rna_to_receptor`, `convert_copies_ALL`, `calib_kdeath`)? | **NONE.** grep of this doc: zero hits. |
| Arm→slot mapping in `ternary_equilibrium` (EQ-0) | **DOC WAS WRONG — FIXED.** `R_A`=CD3/`KD1`=KD_CD3, `R_B`=TAA/`KD2`=KD_TAA (call site :313), so `R_Af`=free **CD3**, `R_Bf`=free **TAA**. The doc had the two labels swapped. |
| Euler monotonicity bound at `dt_h=0.48` | **DOC WAS WRONG (1.44 h) — FIXED to 1.475 h.** Computed: `3/(1/0.48 − 0.05) = 1.4754 h`. The 1.44 figure drops the `k_off` term. |
| "engine default `dt=0.02 d`" | **DOC WAS WRONG — FIXED.** Class default is `dt=0.01` (`coupled_percell_pk.py:65`); `0.02` is passed by the runner (`run_tce_pd_reval.py:143`); `PD_DT` env overrides (`:74`). |
| `R_costim` tagged `[MEASURED]` | **MISTAG — FIXED to [DERIVED — external RNA→protein pipeline; UNVERIFIED].** No measurement exists in either live file; the copies are scRNA-seq-derived. |
| 11 screened arms vs 7 table rows | **NEW DEFECT — BLOCKER-4 added.** `CostimInduction('CD40'/'TNFRSF25'/'CD226'/'TNFRSF14', strict=True)` **returns** (`fold=1.0, kind='unknown', source='NOT IN TABLE'`, stderr only) instead of raising. RUN-VERIFIED, all four. |
| Myeloid IL-6 drive reads the effector-gained `Cb` (QSS only) | **NEW DEFECT — added to EQ-5.** `Cb` is reassigned at `:350`; `:379` reads the gained value; the kinetic path's `:485` (`B2/RC`) is not gain-scaled. |
| All other RUN-VERIFIED claims in §6 (fold=1 identity, 4× ValueError, `a*`=0.952381, `dt_h=1.0` overshoot→1.0, `kill_frac=0.0028471`, `signaling_dynamics` absent, single `costim_ind.step()` call site) | **INDEPENDENTLY REPRODUCED — all hold exactly**, including `kill_frac = 0.0028471054464428436`. |
