---
title: "T3 — Rhoden bivalent binding core (per-cell kinetic, two solvers)"
subsystem: T3
model: costim_engager_counterscreen
source_file: engine/kinetic_rhoden_percell.py  (162 lines, mtime 2026-07-13 14:44, 10,569 bytes)
repo_state: NOT a git repo (no HEAD to pin) — line numbers are against the file as read on 2026-07-13
date: 2026-07-13
generated_by: workflow-subagent T3
adversarially_verified: 2026-07-13 — every equation re-checked at cited line; EQ-1/EQ-5/EQ-12/§5.2 re-run against
  live code (all reproduce); 8 commentary errors found and corrected (1 misattributed quote, 2 wrong line cites,
  1 false consistency claim, 1 unsupported number range, 1 inconsistent timing, 1 wrong provenance tag, 1
  non-monotonicity omission). ZERO fabricated citations. See §6.
---

# T3 — Rhoden Bivalent Binding Core (per-cell kinetic; explicit + backward-Euler solvers)

> **Scope.** This document covers `engine/kinetic_rhoden_percell.py` **only**, top to bottom, all 162 lines.
> Every equation, constant and guard in the file gets an entry. Sibling modules (`multiarm_binding.py`,
> `kinetic_synapse.py`, `coupled_percell_pk.py`, `coupled_percell_pd.py`, `wholebody_percell.py`) are
> referenced only where they are the **caller** or the **source of an input**.
>
> **Dead-code exclusion (enforced).** `unified_binding.py`, `multiarm_kinetic.py`, `biexact_solver.py`,
> `cytokine_pbpk.py`, `il6_pbpk.py`, `rna_to_receptor.py`, `convert_copies_ALL.py`, `calib_kdeath.py` are NOT
> in the live import graph and are NOT documented, cited, or treated as part of the model — even where they
> contain a near-identical `geo_ageff_nM`.

> [!warning]+ Three findings a reviewer will reach for first — stated up front, all RUN-verified this session
> 1. **`rhoden_bivalent_step` (the explicit two-antigen solver) has ZERO live call sites.** It is retained in
>    the live module but nothing in the live import graph calls it any more (§1.4). Every live binding solve
>    goes through the backward-Euler `rhoden_samecell_bivalent_step`.
> 2. **The BE solver's docstring claim "census exactly preserved (no manufacture) even at large AgEFF"
>    (`kinetic_rhoden_percell.py:96–97`) is TRUE of the linear solve and FALSE of the function as written.**
>    The linear solve is exactly census-preserving (I prove the left-invariant in EQ-11 and verify it to
>    machine precision), but the `np.maximum(y, 0)` clamp on line 160 fires whenever the frozen-coefficient
>    step drives a state negative, and *that* manufactures receptor census — up to **+450 % of `Ag_0`
>    cumulatively over 5 simulated days**, with the internalisation flux (the PK TMDD sink) overstated by
>    **+44 % to +124 %** on the affected rows, in the free-drug regime C ≳ 10 nM (§5.2). At C ≲ 5 nM the step
>    is exact (< 0.01 %). The error is **not monotone** in receptor density — see the §5.2 caveat.
> 3. **The singular-row guard's stated rationale is a non-sequitur, but its remedy is sound.** The comment
>    (`:135–136`) argues A = I − hM "cannot be singular" because its diagonal is ≥ 1; a ≥1 diagonal does not
>    imply non-singularity. The correct argument (which I verified numerically, §2 EQ-12) is that **M's
>    spectrum is strictly negative-real-part for all finite non-negative rates**, so 1/h > 0 is never an
>    eigenvalue and A is never singular for h > 0. The guard therefore fires **only** on non-finite input —
>    exactly as the code's own `bad` mask assumes — and the freeze-and-warn behaviour is correct.

---

## 1. PURPOSE & DATAFLOW POSITION

### 1.1 What this subsystem is

T3 is the **literal kinetic bivalent binding kernel**, evaluated **once per cell, per timestep**, for every
cell in the tissue graph (~1e5 cells per organ). It is the single place in the engine where free drug meets
free receptor. Everything downstream — plasma clearance, tumour kill, cytokine release — is a functional of
its outputs.

The module docstring (`:13–15`) states the design invariant explicitly:

> **"ONE binding solve, used IDENTICALLY by PK and PD:
>   PK : TMDD sink = kTMD*(every bound drug species) per cell, summed -> plasma clearance.
>   PD : kill ~ cross-cell bridge (CD3.drug.TAA); CRS ~ engaged dwell; costim ~ costim bridge."**

This matters scientifically: the TMDD sink is not a separately-parameterised clearance term bolted onto a PK
model. It is the *sum over cells* of the same internalisation flux that the PD module reads to compute kill.
A drug that binds more avidly is *automatically* cleared faster and kills more — the coupling is structural,
not fitted.

The second design invariant (`:17–19`) is that **receptor turnover is explicit**, not quasi-steady-state:

> **"Receptor turnover EXPLICIT (the piece the QSS sink lacked -> broke terminal TMDD):
>   dAg/dt has +KSYN (=Ag0*kDEG set-point) and -kDEG*Ag; every bound species internalizes at kTMD.
>   At C=0 steady state Ag -> Ag0 exactly. No QSS assumption anywhere."**

This is the difference between a target-mediated clearance that saturates and recovers on the receptor's own
resynthesis timescale (physiological) and one that snaps instantaneously back to the free-receptor set-point
(QSS — which is what produced the broken terminal phase the docstring is reacting to).

### 1.2 What feeds it (inputs — all owned elsewhere)

| Input | Symbol in code | Units | Supplied by |
|---|---|---|---|
| Free drug concentration, per cell | `C_free_nM` | nM | Organ path: `TissueGraph.C` (interstitial per-cell field, `wholebody_percell.py:179`). Heme/blood sinks: **plasma `C_pl` directly** (`coupled_percell_pd.py:320,334`) |
| Receptor set-point, per cell | `Ag1_0` / `Ag_0` | nM (synapse basis) | `R_copies × NM_PER_COPY` where `NM_PER_COPY = 6.0/257000 = 2.335e-5` nM/copy (`wholebody_pd.py:83`, mirrored `wholebody_percell.py:41,148`) |
| Live receptor pools | `Ag`, `BAg1`, `Bdbl` (`Ag2`,`BAg2`) | nM | **Carried across steps by the caller** — this module is a pure stepper, it owns no state |
| On-rate | `kon` | **/nM/day** | `coupled_percell_pk.py:70`: `kon1 = kon1_perM_pers/1e9*86400` (unit conversion from /M/s) |
| Off-rate | `koff` | **/day** | `coupled_percell_pk.py:71`: `koff1 = koff1_pers*86400`; or `kon*KD` if unset (`coupled_percell_pd.py:135,196`) |
| Receptor degradation | `kDEG` | /day | `coupled_percell_pk.py:72`; heme/blood fallback **0.5/day** (`coupled_percell_pd.py:133,194`) |
| Internalisation of bound complex | `kTMD` | /day | `kint_perday`; `pd_model_config.py:37` default `kint_bridge_perday = 0.9` |
| Effective 2nd-arm concentration | `Ag1EFF`/`AgEFF` | nM | Computed **by this module** — `geo_ageff_nM` (EQ-1) |
| Co-engagement span | `span_nm` | nm | `pd_model_config.py:38–39` (`span_bridge_nm`, `span_cis_nm` = 12.5); module default 12.5 (`:37`) |
| Timestep | `dt` | **days** | `wholebody_percell.TissueGraph.dt` (default 0.01 d = 14.4 min; `PD_DT` env override, `coupled_percell_pk.py:74`) |

### 1.3 What it feeds

- **PK / TMDD.** `intern_flux` (nM/day, synapse basis) → `intern_copies_cell` → `organ_sink` (nmol/day) →
  subtracted from the plasma amount (`wholebody_percell.py:184–185,208`; plasma ODE at
  `coupled_percell_pd.py:357`). The heme and blood sinks feed the same plasma ODE via
  `_heme_sink_nmol_day` / `_blood_sink_nmol_day` (`coupled_percell_pd.py:323,336`).
- **PD readouts.** `S = (BAg1 + Bdbl)/NM_PER_COPY` (bound copies/cell, monovalent + bivalent) and
  `D = Bdbl/NM_PER_COPY` (bivalently-crosslinked copies) — `wholebody_percell.py:186–187`. These are the
  per-cell bound-drug fields that the spatial snapshots and the PD module read.
- **Drug removal from the interstitium.** `loss_nM_day` enters the per-cell transport solve as a sink
  (`wholebody_percell.py:192,206–207`), so binding *depletes the local drug field* — the origin of the
  binding-site barrier / perivascular-shell behaviour.

### 1.4 Which solver runs where — and the dead one

**Live call sites (exhaustive, repo-wide grep this session):**

| Caller | Line | Solver | Compartment |
|---|---|---|---|
| `wholebody_percell.TissueGraph.step` | `:179` | `rhoden_samecell_bivalent_step` | every organ, every cell (the main PK/PD path) |
| `coupled_percell_pd.simulate_pd` | `:319` | `rhoden_samecell_bivalent_step` | heme (circulating malignant blast) TMDD sink |
| `coupled_percell_pd.simulate_pd` | `:333` | `rhoden_samecell_bivalent_step` | normal circulating blood TMDD sink |
| `wholebody_percell` / `coupled_percell_pd` | `:178`, `:318`, `:332` | `geo_ageff_nM` | effective 2nd-arm concentration for all three |

**`rhoden_bivalent_step` (`:47–80`) is called by NOTHING in the live graph.** The only surviving references
are (a) the comment at `wholebody_percell.py:177` that records its removal — *"Replaces the old two-pool
rhoden_bivalent_step call whose Bdbl was INERT for single-antigen (fed from BAg2=0)"* — and (b) archived
snapshots under `model_snapshots/pre_unification_20260713_050844/` and `model/_archived_2026-07-13/`, which
are not on the import path. It is documented here in full because it is **live code in a live module** (a
reader will find it, and the two-antigen scheme it implements is the reference form the whole subsystem is
derived from), but the honest statement is: **the heterotypic two-antigen path described in the module
docstring `:29–32` is currently NOT WIRED.** Every live solve is the reduced same-antigen 3-state system.

### 1.5 Why the BE solver exists (the 2026-07-13 stall)

Provenance is in-repo, not reconstructed: `reference_unified_binding/LIVE_ENGINE_FIXES.md:76–79` and
`CHANGELOG_2026-07-13.md:238–250`.

The explicit solver sizes its substep count from the **fastest pseudo-first-order rate**, and that rate is
`kon·AgEFF` — the avidity crosslink term, which is **linear in kon** (`:64–66`). While every molecule ran on
the generic fallback `kon = 1e5 /M/s`, `nsub` sat at a few hundred and the run was merely slow. When the
measured kinetics were wired in (`CHANGELOG_2026-07-13.md:236`, corroborated at `:242`: teclistamab
**kon_TAA = 1.28e6 /M/s**, 12.8× the generic default), `nsub` scaled by the same 12.8× and the run stalled.

**MEASURED, quoted verbatim from `LIVE_ENGINE_FIXES.md:76–79`** (isolated timing, 100 k cells):

> *"at kon 1e5→1.28e6, `rhoden_bivalent_step` nsub **442→5651**, wall **0.805 s→7.83 s**;
> `rhoden_samecell_bivalent_step` flat **0.026 s** (kon-independent) → **~300× faster** on the heme/blood step
> at teclistamab's kon."*

`py-spy` on the stuck run (PID 1767409) showed the main thread frozen inside `rhoden_bivalent_step` at
`kinetic_rhoden_percell.py:75` (`CHANGELOG_2026-07-13.md:239–240`); one PD step took minutes, and the run
never reached step 1.

**RUN-VERIFY (this session, independent):** the two reported `nsub` values are *exactly* consistent with the
code's own law (`:64–66`), which is strictly linear in `kon`:
`5651 / 442 = 12.785` vs `1.28e6 / 1e5 = 12.8` — agreement to 0.1 %. Inverting the law, `nsub = 5651` at
`dt = 0.01 d` implies `rmax = 0.1·5651/0.01 = 56,510 /day`, which at teclistamab's `kon = 110.6 /nM/day`
(= 1.28e6/1e9×86400) is dominated by a single `kon·AgEFF` term with `AgEFF ≈ 5.1e2 nM` — i.e. a BCMA density
of order 2×10³ copies/cell. **The `nsub` ratio reproduces** (it is a pure consequence of the `:64–66` law).
My own timing on an M-series Mac (N = 1e5 cells, `dt = 0.01 d`, generic kon, receptor lognormal about 5e4
copies; `T3_gaps.py`) gave the explicit step **29.8 s at nsub = 20,000 (the cap)** versus **0.039 s** for the
BE step — a **770×** ratio, the same order as the reported ~300×. (Wall-clock numbers are machine- and
load-dependent and are **not** a reproduction of the reported 0.805 s / 7.83 s / 0.026 s, which were measured
on different hardware; only the *ratio's order* is being corroborated here.)

**The crucial asymmetry:** the explicit solver's `nsub` is set by a **stability** requirement (h·λ_max ≲ 0.1),
which scales with the *stiffest* rate in the system — and `kon·AgEFF` is ~1e4–1e6 /day. The BE solver's `nsub`
is set by an **accuracy** requirement on the *slow* scales only (`kDEG + 2·kTMD` ≈ 2.3 /day) because
backward-Euler is A-stable — the stiff mode needs no resolution to remain bounded. That is the whole reason the
BE step is flat in kon: **kon does not appear in its substep law at all** (`:108`).

---

## 2. GOVERNING EQUATIONS

**Global unit convention (module-wide, `:51`):** rates are **/day**, concentrations are **nM**, `dt` is in
**days**. `AVO = 6.02214076e23 /mol` (`:35`). All state arrays are per-cell and vectorised.

---

### EQ-1 — Rhoden geometric effective 2nd-arm concentration `AgEFF` (`kinetic_rhoden_percell.py:37–45`)

```
rec_pc  = max(rec_pc, 0)                                             (:39)
SA_cell = 4·π·r_cell²                                    [µm²]       (:40)
r_Ab    = max(span_nm, 1e-3) · 1e-3                      [µm]        (:41)
SA_Ab   = π·r_Ab²                                        [µm²]       (:42)
V_Ab    = (2/3)·π·r_Ab³                                  [µm³]       (:42)
Ag_bulk = rec_pc · 1e9 / AVO · 1e9                                   (:43)
Am_cell = Ag_bulk/1e9 ;  Am_SA = Am_cell/SA_cell ;  Ag_SA = Am_SA·SA_Ab   (:44)
AgEFF   = Ag_SA / V_Ab · 1e15                            [nM]        (:45)
```

**Algebraically (the `1e9` factors cancel exactly; verified below), this collapses to:**

```
                       rec_pc          SA_Ab            rec_pc        3
    AgEFF [nM]  =  ─────────────── · ─────────  =  ─────────────── · ─────  · (1e15/AVO·1e9)
                      SA_cell          V_Ab           SA_cell        2·r_Ab

              =  σ_receptor · 3/(2·r_Ab)   [copies/µm³]   × 1.6606  [nM per copy·µm⁻³]
```

because `SA_Ab/V_Ab = π r² / ((2/3)π r³) = 3/(2r)` exactly.

**(b) Biological meaning.** An antibody that is already bound by one arm cannot wander. Its free arm is
tethered, and can only explore a **hemisphere of radius = the arm-to-arm span** centred on the anchored
paratope. The number of *partner* receptors it can possibly reach is the surface density σ times the **disc**
that hemisphere projects onto the membrane (`SA_Ab = π r²`). The *concentration* that free arm experiences is
that count divided by the **volume** it explores (`V_Ab = (2/3)π r³`, the hemisphere). So `AgEFF` is a genuine
molarity: *"how concentrated does the second antigen look, from the point of view of a tethered Fab?"* This is
the entire physical content of avidity in this model — the second binding event is not given a bonus rate
constant, it is given the **same kon** acting on a **much larger local concentration**.

**(c) Mechanistic rationale, and what it rejects.** The alternative — and the thing that is *not* done here —
is to hand the bivalent complex an empirical "avidity factor" or a folded-down apparent KD. That would be a
fitted constant. Here the avidity enhancement is *derived from geometry*: it depends only on the receptor
copy number, the cell radius, and the molecular span, all of which are independently specified. Note the
scaling that falls out and that a reviewer should check: **`AgEFF ∝ 1/span`, not `1/span³`** — because the
numerator (reachable partners) grows as `span²` while the denominator (explored volume) grows as `span³`. A
longer, floppier linker therefore *dilutes* avidity only linearly. The module docstring's sibling
(`multiarm_binding.py:30–31`) states this in words: *"Larger span -> larger explored shell -> LOWER c_eff
(dilution); shorter -> higher."*

**(d) Units.** `rec_pc` [copies/cell]; `r_cell_um` [µm]; `span_nm` [nm]; `SA_cell`,`SA_Ab` [µm²];
`V_Ab` [µm³]; `AgEFF` [nM]. The trailing `1e15` converts copies·µm⁻³ → copies·L⁻¹, and the `/AVO·1e9` pair
converts to nanomolar.

**Attribution — read the code, not the folklore.** The module header (`:3–5`) attributes the **bivalent
crosslink structure** (BAg1/BAg2/Bdbl with a geometric `AgEFF`) to *"Rhoden et al. 2016 (kinetic form: bioRxiv
10.1101/2022.09.12.507653)"*. That is the only citation in the file. **[UNVERIFIED CITATION]** — I did not
fetch the DOI in this task, and the header gives no PMID; the *form* of the equation is standard tethered-ligand
effective-concentration geometry, but the specific attribution is asserted by the comment, not proven here.

**RUN-VERIFY (this session, `T3_rhoden_verify.py` against the live import):**

| `rec_pc` (copies) | `geo_ageff_nM(rec, 8.0, 12.5)` | closed form `σ·3/(2r_Ab)·1.6606` | ratio |
|---|---|---|---|
| 1,000 | 247.7653 nM | 247.7653 nM | 1.0000000000 |
| 10,000 | 2,477.6531 nM | 2,477.6531 nM | 1.0000000000 |
| 100,000 | 24,776.5313 nM | 24,776.5313 nM | 1.0000000000 |
| 500,000 | 123,882.6566 nM | 123,882.6566 nM | 1.0000000000 |

Per-copy slope at (r_cell = 8 µm, span = 12.5 nm): **0.2477653 nM per receptor copy**. Span scaling verified:
`AgEFF(12.5)/AgEFF(25) = 2.000000` exactly — confirming the **1/span** law, not 1/span³.
**Magnitude sanity:** a 1e5-copy cell presents its own tethered second arm an effective **~25 µM** of partner
antigen. Against a KD of ~1 nM that is a ~2.5e4 avidity drive — which is precisely why the crosslink term is
the stiffest rate in the system and why the explicit solver died.

---

### EQ-2 — Receptor synthesis set-point `KSYN` (`:57` explicit; `:104` BE)

```
KSYN_k = Ag_k0 · kDEG                                    [nM/day]
```

**(b) Biological meaning.** Zero-order receptor synthesis, pinned so that the *drug-free* steady state of the
free-receptor pool is exactly the measured/imputed baseline: at `C = 0`, `dAg/dt = KSYN − kDEG·Ag = 0 ⇒
Ag → Ag_0`. The cell is not "given" its receptor number — it *maintains* it against first-order turnover.

**(c) Mechanistic rationale.** This is the piece the docstring (`:17–19`) says the QSS sink lacked. With a QSS
receptor pool, receptor lost to internalisation is replaced instantaneously, so the target sink never
saturates in the way real TMDD does and the terminal PK phase is wrong. Making synthesis explicit means
**target recovery has a timescale** (`1/kDEG`), so a second dose lands on a partially-depleted target — which
is the mechanism behind step-up dosing tolerance. The alternative (fitting a first-order `k_syn` independent of
`Ag_0`) was rejected because it would let the drug-free receptor density drift away from the measured value.

**(d) Units.** `Ag_0` [nM], `kDEG` [/day], `KSYN` [nM/day].

**RUN-VERIFY:** starting from `Ag = 0` with `C = 0` for 40 simulated days, the BE stepper recovers
`Ag/Ag_0 = [1. 1. 1. 1.]` across a 1e3–3e5 copy range (5 significant figures). The `C = 0` fixed point is exact.

---

### EQ-3 — Fractional-availability scaling of the crosslink `x` (`:58–59, 70`; `:105, 116`)

```
a_r  = 1/Ag_0   if Ag_0 > 1e-30   else 0                 (:58–59, :105)
x_k  = Ag_kEFF · (Ag_k · a_rk)  =  Ag_kEFF · (Ag_k / Ag_k0)      [nM]     (:70, :116)
```

**(b) Biological meaning.** The tethered second arm does not see the *baseline* partner density — it sees
whatever partner density is **still free right now**. As the drug consumes free receptor, the local partner
field thins out. `Ag/Ag_0 ∈ [0,1]` is the surviving free fraction, and it scales the geometric `AgEFF` down
proportionally.

**(c) Mechanistic rationale.** This makes the crosslink term **depletion-aware and self-limiting**: at
saturation the free-receptor fraction collapses, the crosslink on-rate collapses with it, and `Bdbl` formation
stops — automatically producing the bell-shaped (hook-effect) dependence of crosslink on drug concentration
without any explicit hook term. The rejected alternative — using the constant geometric `AgEFF` — would let a
saturated cell keep forming crosslinks out of receptors it no longer has, and (in a clamped explicit scheme)
manufacture mass.

**(d) Units.** `a_r` [1/nM]; `Ag_kEFF` [nM]; `x_k` [nM]. **Guard:** the `1e-30` threshold means a cell with
**zero receptors** gets `a_r = 0 ⇒ x = 0 ⇒ no crosslink` — fail-closed, and finite (this matters for EQ-12).

---

### EQ-4 — The reference 6-species two-antigen scheme (`:71–75`; docstring `:22–27`)

```
dAg1  = KSYN1 − kon1·C·Ag1 + koff1·BAg1 − kon1·x1·BAg2 + koff1·Bdbl − kDEG·Ag1       (:71)
dAg2  = KSYN2 − kon2·C·Ag2 + koff2·BAg2 − kon2·x2·BAg1 + koff2·Bdbl − kDEG·Ag2       (:72)
dBAg1 = kon1·C·Ag1 − koff1·BAg1 − kon2·x2·BAg1 + koff2·Bdbl − kTMD·BAg1              (:73)
dBAg2 = kon2·C·Ag2 − koff2·BAg2 − kon1·x1·BAg2 + koff1·Bdbl − kTMD·BAg2              (:74)
dBdbl = kon1·x1·BAg2 + kon2·x2·BAg1 − (koff1+koff2)·Bdbl − 2·kTMD·Bdbl                (:75)
```

**(b) Biological meaning — term by term.**
- `kon_k·C·Ag_k` : **first-bond capture from bulk.** Free drug in solution grabs a free receptor. Bimolecular,
  driven by the *bulk* drug concentration.
- `koff_k·BAg_k` : **first-bond release.** The singly-bound complex lets go and returns drug to bulk.
- `kon_k·x_k·BAg_j` : **second-bond crosslink.** A drug already anchored by arm *j* swings its free arm *k*
  onto a partner receptor. **Identical `kon`** — the *only* difference from the first bond is that the
  concentration is the tethered `x_k` (EQ-1/EQ-3) instead of the bulk `C`. This is the avidity mechanism, and
  it is emergent, not parameterised.
- `koff_k·Bdbl` (appearing in *both* `dAg_k` and `dBAg_j`) : **crosslink dissolution.** One arm of a
  doubly-bound complex releases → one receptor goes free AND the complex demotes to singly-bound. Both
  bookkeeping entries are required for census closure.
- `kDEG·Ag_k` : constitutive turnover of **free** receptor only (bound receptor is protected from the
  degradation route and instead follows the drug).
- `kTMD·BAg_k`, `2·kTMD·Bdbl` : **internalisation of the complex.** The doubly-engaged complex internalises at
  **twice** the rate — see (c).

**(c) Mechanistic rationale, and the assumption a reviewer should attack.** The factor **2** on `kTMD·Bdbl`
encodes *"each engaged arm contributes an independent internalisation hazard"* — crosslinking two receptors
doubles the endocytic drive. This is a **modelling assumption, not a measurement**; nothing in this file sources
it. It is not innocuous: it means the avidity advantage of a bivalent binder shows up *twice* — once in
increased occupancy (via `AgEFF`) and again in a doubled internalisation rate — and it is the direct driver of
the extra TMDD clearance of bivalent constructs. An equally defensible alternative (one complex = one
internalisation event, rate `kTMD`, carrying 2 receptors) would give a materially different TMDD magnitude.
See §5.4.

**Provenance of the scheme itself.** Docstring `:5–8`: *"Rhoden's own model assumes antigen kinetics NEGLIGIBLE
(no turnover); the receptor-turnover terms here (KSYN=Ag0*kDEG synthesis, -kDEG*Ag free-antigen degradation,
kTMD internalization of bound species) are ADDED on top -- standard TMDD receptor turnover, taken VERBATIM from
the user-supplied MATLAB scheme."* So: **crosslink structure ← Rhoden 2016 [UNVERIFIED CITATION]; turnover
terms ← a user-supplied MATLAB scheme, transcribed verbatim [UNSOURCED — TBD]** (no paper is named for the
turnover block).

**(d) Units.** All states [nM]; `kon` [/nM/day]; `koff`,`kDEG`,`kTMD` [/day]; `C`,`x` [nM]; derivatives [nM/day].

**Census (derived, this session).** Weighting the two conserved receptor pools:
`Σ1 = Ag1 + BAg1 + Bdbl` and `Σ2 = Ag2 + BAg2 + Bdbl` (the crosslinked complex holds exactly **one copy of
each** antigen). Summing (`:71`+`:73`+`:75`), every crosslink and koff term cancels identically, leaving
`dΣ1/dt = KSYN1 − kDEG·Ag1 − kTMD·BAg1 − 2·kTMD·Bdbl`. The ODE is exactly receptor-conserving up to turnover.
**The discretisation is not** — see EQ-6.

---

### EQ-5 — Explicit substep count `nsub` (stability-limited) (`:60–66`)

```
rmax = max over cells of [ kon1·C + koff1 + kon2·C + koff2
                           + kon1·max(Ag1EFF,0) + kon2·max(Ag2EFF,0)
                           + kDEG + 2·kTMD ]                              [/day]     (:64–65)
nsub = int( min( max(1, ceil(rmax·dt/0.1)), 20000 ) )                                (:66)
h    = dt / max(nsub, 1)                                                  [days]     (:67)
```

**(b) Meaning.** The sum of **all pseudo-first-order rates** in the system — an upper bound on the spectral
radius. The step is refined until `h·rmax ≤ 0.1`, i.e. no state may change by more than ~10 % of itself in one
substep. The `20000` ceiling is a wall-clock guard.

**(c) Mechanistic rationale, from the code's own comment (`:61–63`):**
> *"ALL pseudo-first-order rates incl. the AVIDITY CROSSLINK (kon*AgEFF ~ hundreds/day for real densities) —
> omitting it under-substeps -> explicit-Euler overshoot -> clamp injects mass."*

This comment is the scar tissue from an earlier bug: the crosslink rate had been left out of `rmax`, the step
overshot, the non-negativity clamp (EQ-6) caught the negative and thereby *created* receptors. Including
`kon·AgEFF` fixed the mass injection — **and simultaneously made the solver unusable**, because `kon·AgEFF` is
the largest rate in the problem by 3–5 orders of magnitude. **The fix and the failure are the same line.**

**(d) Units.** `rmax` [/day]; `dt`,`h` [days]; `0.1` = dimensionless relative-change tolerance; `nsub` = count.

**RUN-VERIFY (this session; kon = 8.64 /nM/day = the generic 1e5 /M/s, KD 1.45 nM, kDEG 0.5, kTMD 0.9,
dt = 0.01 d, C = 100 nM):**

| receptors/cell | `AgEFF` (nM) | `rmax` (/day) | `nsub` (explicit) | `nsub` (BE, EQ-9) |
|---|---|---|---|---|
| 1e4 | 2,477.7 | 44,569 | **4,457** | **1** |
| 1e5 | 24,776.5 | 429,894 | **20,000 (CAP HIT)** | **1** |
| 5e5 | 123,882.7 | 2,142,448 | **20,000 (CAP HIT)** | **1** |
| 1e6 | 247,765.3 | 4,283,140 | **20,000 (CAP HIT)** | **1** |

**Finding (§5.1):** at ≥1e5 receptors/cell the explicit solver **saturates its own cap**, so the stability
criterion it exists to enforce is silently violated (at 1e5 copies the realised `h·rmax = 0.215`, i.e. 2.15× the
0.1 target). The cap converts a slow-but-correct solver into a fast-but-wrong one with no warning. CD20 on a
B cell is ~1e5–3e5 copies — this is not an exotic corner.

---

### EQ-6 — Explicit Euler update + non-negativity clamp (`:76–78`)

```
Ag1 ← max(Ag1 + h·dAg1, 0) ;  Ag2 ← max(Ag2 + h·dAg2, 0)                 (:76)
BAg1 ← max(BAg1 + h·dB1, 0) ; BAg2 ← max(BAg2 + h·dB2, 0)                (:77)
Bdbl ← max(Bdbl + h·dBd, 0)                                              (:78)
```

**(b/c) Meaning and rationale.** Forward Euler with a hard floor at zero. The clamp is a **physical-realisability
guard** (a concentration cannot be negative) but it is **not conservative**: clamping a negative state to zero
*adds* the deficit to the census. It is precisely the mass-injection channel that `:62–63` warns about. In a
correctly-substepped explicit scheme the clamp should never fire; when `nsub` hits its cap (EQ-5 RUN-VERIFY) it
can and does.

**(d) Units.** All [nM]; `h` [days].

---

### EQ-7 — Explicit internalisation flux accumulator (`:68, 79, 80`)

```
intern_accum += kTMD·(BAg1 + BAg2 + 2·Bdbl)·h                            (:79)
intern_flux   = intern_accum / max(dt, 1e-30)                 [nM/day]   (:80)
```

**(b) Biological meaning.** The **drug** internalisation rate. One `BAg` complex = one drug molecule, taken in
at `kTMD`. One `Bdbl` = one drug molecule, but doubly engaged, taken in at `2·kTMD` — hence the factor 2. Note
the asymmetry with the *receptor* census: a `Bdbl` internalisation event removes **two** receptor copies but
only **one** drug molecule, and the two ledgers are separately consistent (verified in EQ-4/EQ-11).

**(c) Rationale.** Returning a **step-averaged** flux (rather than the end-of-step instantaneous rate) means the
caller's coarse `dt` sees the correctly time-integrated sink even though this module substeps internally — the
PK mass balance stays exact across the interface. Note the accumulation uses the **post-clamp, post-update**
pools (`:79` follows `:76–78`), so any clamp-injected mass is *also* counted into the sink.

**(d) Units.** `intern_accum` [nM]; `intern_flux` [nM/day].

---

### EQ-8 — Reduced same-antigen 3-state system (the LIVE system) (`:83–98` docstring; implemented `:123–129`)

Both arms bind the **same** antigen (glofitamab 2×CD20, alnuctamab 2×BCMA — docstring `:85–86`):

```
states: Ag (free)   BAg1 (one arm bound)   Bdbl (both arms -> 2 neighbour copies)

dAg   = KSYN − kon·C·Ag + koff·BAg1 − kon·x·BAg1 + 2·koff·Bdbl − kDEG·Ag     (:90)
dBAg1 = kon·C·Ag       − koff·BAg1 − kon·x·BAg1 + 2·koff·Bdbl − kTMD·BAg1     (:91)
dBdbl =                              kon·x·BAg1 − 2·koff·Bdbl  − 2·kTMD·Bdbl  (:92)
        x = AgEFF·Ag/Ag_0 ;  KSYN = Ag_0·kDEG                                 (:93)
        census Ag_tot = Ag + BAg1 + 2·Bdbl                                    (:94)
```

**(b) Biological meaning.** The **correct reduction when Ag1 ≡ Ag2**. `Bdbl` now consumes **two copies of the
same** receptor, so it carries **census weight 2** (contrast EQ-4, where it carried one copy of each of two
distinct antigens). The `2·koff` on the `Bdbl → BAg1` transition is the statistical factor for *"either of the
two identical bound arms may release"*.

**(c) Mechanistic rationale — why this replaced the two-pool call.** The comment at `wholebody_percell.py:177`
is explicit: the old code called the two-antigen `rhoden_bivalent_step` with `Ag2 = 0` for a single-antigen
target, and `dBdbl = kon1·x1·BAg2 + kon2·x2·BAg1` with `BAg2 ≡ 0`, `kon2 ≡ 0` is **identically zero**. So
`Bdbl` was **inert** — a bivalent same-antigen binder was silently being simulated as monovalent, with **no
avidity at all**. That is not a numerics bug, it is a *physics* bug, and it was live until 2026-07-13
(`reference_unified_binding/LIVE_ENGINE_FIXES.md:57–58` flags it as fixed in the same pass: *"Switching also
FIXES the inert-`Bdbl` avidity bug for glofitamab's heme binding"*; the independent audit record is
`reference_unified_binding/audit_confirmed.json` findings 2/3/7, and `CHANGELOG_2026-07-13.md:223` records the
consequence — glofitamab is retained as the *bivalent-avidity mechanism* test *"(does Bdbl form via
rhoden_samecell_bivalent_step)"*).

**(d) Units.** As EQ-4. `2·koff`, `2·kTMD` [/day].

---

### EQ-9 — BE substep count `nsub` (accuracy-limited, kon-INDEPENDENT) (`:107–109`)

```
nsub = int( min( max(1, ceil( (kDEG + 2·kTMD)·dt / 0.25 )), 16 ) )       (:108)
h    = dt / max(nsub, 1)                                                 (:109)
```

with the in-code justification `# SLOW scales only (BE is stable)` (`:108`) and, in the docstring (`:96–97`):
> *"BACKWARD-EULER, receptors-as-states -> unconditionally stable + census exactly preserved (no manufacture)
> even at large AgEFF. nsub set by SLOW scales (kDEG/kTMD), not the stiff crosslink."*

**(b/c) Meaning and rationale.** Because backward-Euler is A-stable, the stiff crosslink mode (`kon·x`, up to
~1e6 /day) is *damped* rather than amplified at any `h` — it needs no resolution for **stability**. Only the
slow physiological scales (`kDEG ≈ 0.5`, `kTMD ≈ 0.9 /day`) need resolving for **accuracy**. Hence `kon` is
absent from `:108`, and hence the flat 0.026 s.

**(d) Units.** `kDEG`,`kTMD` [/day]; `dt`,`h` [days]; `0.25` = dimensionless tolerance; cap 16.

**RUN-VERIFY (production values kDEG = 0.5, kTMD = 0.9, dt = 0.01 d):**
`(0.5 + 1.8)·0.01/0.25 = 0.092 → ceil = 1 → nsub = 1`. **In production this solver takes exactly ONE
backward-Euler solve per step, with h = dt = 0.01 d = 14.4 min.** That is the configuration whose accuracy is
audited in §5.2 — and it is where the residual problem lives.

---

### EQ-10 — The backward-Euler linear solve (`:123–133, 159–160`)

Linear generator `M` over `y = [Ag, BAg1, Bdbl]ᵀ`, with the **receptors carried as states** (`:119`):

```
        ⎡ −(kb + kDEG)      (koff − kx)        2·koff            ⎤
   M =  ⎢  kb              −(koff + kx + kTMD)  2·koff           ⎥        (:123–125)
        ⎣  0                 kx                −(2·koff + 2·kTMD)⎦

   kb = kon·C   [/day]   (:117)          kx = kon·x   [/day]   (:118)
   s  = [KSYN, 0, 0]ᵀ                                            (:130)
   A  = I − h·M                                                  (:132)
   rhs = y₀ + h·s                                                (:133)
   y  = A⁻¹ · rhs      (batched over cells, np.linalg.solve)     (:159)
   y ← max(y, 0)                                                 (:160)
```

**(b) Biological meaning of "receptors as states".** The free-receptor pool `Ag` is not eliminated by a
quasi-steady-state assumption and is not treated as a parameter — it is a **row of the generator**, so every
receptor that leaves the free pool arrives somewhere else *inside the same implicit solve*. Binding, release,
crosslinking and turnover are all resolved simultaneously and self-consistently at the end of the step.

**(c) Mechanistic rationale, and the one honest caveat.** The implicit solve is what buys A-stability. But note
`x` (hence `kx`) is computed **at the start of the substep from the old `Ag`** (`:116`, comment: *"effective
2nd-arm conc, frozen within the BE solve"*). So this is a **frozen-coefficient / IMEX backward-Euler**, not a
full Newton solve of the nonlinear system. The linear part is unconditionally stable; the *nonlinear* coupling
(`x` depends on `Ag`) is treated **explicitly**, and carries an O(h) splitting error that is **not** controlled
by the `:108` substep law. **This is the origin of the accuracy defect in §5.2.** A Newton iteration on the
3×3 system (2–3 iterations would suffice) would remove it at negligible cost.

**Structural note (important for EQ-12):** the entry `M[0,1] = koff − kx` is **negative whenever kx > koff**,
which at real densities means *always* (RUN: at 3e5 copies, `koff = 12.5 /day`, `kx = 642,208 /day`, so
`M[0,1] = −642,195`). A matrix with negative off-diagonals is **not a Metzler matrix**, so the flow is **not
guaranteed to preserve non-negativity** — which is exactly why the `max(y,0)` clamp on `:160` exists, and
exactly why it fires.

**(d) Units.** `M` [/day]; `h` [days]; `A` dimensionless; `y`,`rhs`,`s·h` [nM]; `s` [nM/day].

---

### EQ-11 — The census left-invariant (`:94` docstring; derived + verified here)

Let `w = [1, 1, 2]ᵀ` (a `Bdbl` holds two receptor copies). Then the **weighted column sums of M** are:

```
 wᵀM  =  [ −kDEG ,  −kTMD ,  −4·kTMD ]
```
Column-by-column (multiply the M entries above by `w` and add):
- col 0: `−(kb+kDEG) + kb + 0 = −kDEG` — binding moves receptors, only degradation destroys them.
- col 1: `(koff−kx) − (koff+kx+kTMD) + 2·kx = −kTMD` — release and crosslink cancel exactly; only
  internalisation removes census.
- col 2: `2·koff + 2·koff + 2·(−2koff−2kTMD) = −4·kTMD` — a `Bdbl` internalising at `2·kTMD` takes **2**
  receptor copies with it → `2 × 2·kTMD`.

**Consequence (the actual theorem the docstring is reaching for):** if `kDEG = kTMD = 0` then `wᵀM = 0`, so
`wᵀA = wᵀ(I − hM) = wᵀ`, hence `wᵀA⁻¹ = wᵀ`, hence `wᵀy = wᵀy₀` — **the backward-Euler solve conserves receptor
census EXACTLY, at any h, for any AgEFF.** With turnover on, it exactly integrates
`d(census)/dt = KSYN − kDEG·Ag − kTMD·BAg1 − 4·kTMD·Bdbl` (implicit form).

**RUN-VERIFY:** across receptor densities 1e3–3e5 and C = 0.1–50 nM, the **pre-clamp** census after one BE
step matches the exact implicit census law to **all printed digits** (12 test cells, `T3_census_diag.py`).
The theorem holds.

**AND YET — the function does not.** The clamp on `:160` is applied *after* the solve. Once `y` acquires a
negative component (which the non-Metzler `M` permits, see EQ-10), `max(y,0)` **breaks the invariant**. This is
quantified in §5.2. **The docstring claim at `:96–97` is a statement about the linear algebra, and it is
correct; it is not a statement about the function, and as a statement about the function it is false.**

---

### EQ-12 — The SINGULAR-ROW GUARD (added 2026-07-13) (`:134–158`)

```
bad  = ¬( isfinite(A).all(axis=(1,2)) ∧ isfinite(rhs).all(axis=1) )      (:142)
det  = det(A) on the finite rows only (bad rows get det := 1.0)          (:143–147)
sing = bad ∨ (|det| < 1e-12)                                             (:148)
if any(sing):  warn ONCE ; A[sing] := I ; rhs[sing] := nan_to_num(y₀[sing])   (:149–158)
```

**(b) What it does, biologically.** A cell whose linear system is degenerate has **no well-posed update**, so
the guard **freezes it** (`A := I`, `rhs := y₀` ⇒ `y = y₀`) rather than crashing the run or letting a NaN
propagate into the tissue field. The code's own justification (`:139–141`):
> *"A cell with a degenerate matrix has no well-posed update, so FREEZE it (y = y0) rather than crash or
> silently propagate NaN. Freezing is the physically-correct no-op for a fully depleted cell, but it is
> COUNTED and WARNED so it can never hide a real bug."*

**(c) Why it was added — the failure it prevents (`:137–138`, verbatim):**
> *"Observed 2026-07-13: elranatamab at TSIM=24 crashed the whole run with LinAlgError -- one bad cell out of
> ~1e5 killed a 40-minute simulation."*

`np.linalg.solve` on a batched `(n,3,3)` stack raises `LinAlgError` if **any** matrix in the batch is singular.
There is no partial-failure mode: one pathological cell out of ~100,000 terminates the process, and 40 minutes
of tissue transport, PD and cytokine integration are lost. The guard converts a **run-fatal exception** into a
**counted, warned, single-cell no-op** — a proportionate response to a defect whose blast radius was absurdly
disproportionate to its size.

**(c′) The stated rationale is wrong; the conclusion is right.** The comment (`:135–136`) argues:
> *"A = I - h*M has diagonal 1 + h*(non-negative rates) >= 1, so for FINITE inputs it cannot be singular."*

A diagonal ≥ 1 does **not** imply non-singularity (row 0 is not diagonally dominant once `h·kx ≫ 1`: at 3e5
copies, `|A[0,1]| = h·|koff − kx| ≈ 6422` versus `A[0,0] ≈ 5.3`). **The correct argument** — which I verified —
is spectral: `A = I − hM` is singular **iff `1/h` is an eigenvalue of M**. Since `h > 0`, singularity requires
`M` to have a **positive real eigenvalue**. It does not:

**RUN-VERIFY (this session):** sweeping `kx` over `[1e-2, 1e8] /day` (200 log-spaced points, `kb = 1`,
`koff = 1`, `kDEG = 0.5`, `kTMD = 0.9`): the number of `kx` values for which `M` has **any** eigenvalue with
positive real part is **0**. At the extreme `kx = 1e8` the spectrum is
`{−1.0000e+08, −2.5000, −1.8000}` — strictly negative. Minimum `|det(A)|` over the whole sweep at `h = 0.01`
was **1.0736** (nowhere near the 1e-12 trigger). **So: for finite non-negative rates A is never singular, and
the guard's `|det| < 1e-12` branch is unreachable — the guard fires only via the `bad` (non-finite) mask.**
The code's `bad` mask is therefore the operative one, and its diagnosis (*"a NON-FINITE rate leaked in (NaN/Inf
drug conc...)"*, `:136`) is the right one.

One sub-claim in the comment is nonetheless **incorrect**: `:136–137` also blames *"Ag_0=0 giving a degenerate
a_r"*. It cannot — `a_r` is explicitly zero-guarded at `:105` (`np.where(Ag_0 > 1e-30, 1/Ag_0, 0)`), so a
zero-receptor cell yields `x = 0`, `kx = 0`, a perfectly finite `M`, and a non-singular `A` (verified: the
`rec = 0` cell in my 5-cell test integrates cleanly and stays identically zero). **The only real trigger is a
non-finite drug concentration.**

**RUN-VERIFY of the guard itself (this session):** inject `C[7] = NaN` into a 100,000-cell array →
`rhoden_samecell_bivalent_step` **does not raise**, emits exactly the intended warning
(`[rhoden-BE] 1/100000 cells had a SINGULAR/non-finite matrix -> FROZEN this substep ... Non-finite inputs: 1`),
returns **all-finite** arrays, and the 99,999 healthy cells step normally (neighbour cell 8 advances to
`Ag = 0.0603, BAg1 = 6.11e-4, Bdbl = 0.2550` while cell 7 is neutralised). **The guard holds.**

**Contrast — the explicit solver has NO guard (verified):** feeding the same NaN to `rhoden_bivalent_step`
(i) makes `rmax = NaN` at `:64`, (ii) `max(1, NaN)` returns `1` in Python, so **`nsub` silently collapses to 1
for the ENTIRE array** — the stability control is destroyed for all 100,000 cells, not just the bad one — and
(iii) the NaN propagates into the output with **no exception and no warning**. Verified: output
`Ag1 = [2.10128, nan, 1.90256]`.

**(d) Units.** `det(A)` dimensionless; `1e-12` a dimensionless threshold (see §5.5 — it is **not**
scale-invariant, though on present evidence it is also never reached).

---

### EQ-13 — BE internalisation flux accumulator (`:112, 161, 162`)

```
intern += kTMD·(BAg1 + 2·Bdbl)·h            (post-update, post-clamp)    (:161)
return ... , intern / max(dt, 1e-30)         [nM/day]                    (:162)
```

Same construction and same step-averaging rationale as EQ-7, with the reduced species set. `BAg1` internalises
at `kTMD`, `Bdbl` (one drug, two receptors, doubly engaged) at `2·kTMD`.

**(d) Units.** `intern` [nM]; returned flux [nM/day]. The caller converts it to copies/cell/day by dividing by
`NM_PER_COPY` (`wholebody_percell.py:184`).

---

### EQ-14 — Defensive input handling (small but load-bearing)

| Guard | Line | What it does | Why |
|---|---|---|---|
| `C = max(C, 0)` | `:53`, `:100` | negative drug → 0 | transport solve can undershoot |
| `rec_pc = max(rec_pc, 0)` | `:39` | negative receptors → 0 | same |
| `r_Ab = max(span_nm, 1e-3)·1e-3` | `:41` | span floor 1e-3 nm | `V_Ab ∝ r³` in a denominator → div-by-zero |
| `a_r = 0 if Ag_0 ≤ 1e-30` | `:58–59`, `:105` | zero-receptor cells → no crosslink | fail-closed; keeps `M` finite |
| `max(nsub, 1)` | `:67`, `:109` | h is finite | — |
| `/max(dt, 1e-30)` | `:80`, `:162` | flux average is finite | dt=0 call |
| `np.array(..., copy=True)` | `:54–56`, `:101–102` | **no in-place mutation of caller state** | the caller owns the pools; a stepper that mutated them would double-count |
| scalar/array broadcast | `:110–114` | accepts scalar or per-cell `Ag`, `C`, `KSYN` | heme/blood pass a scalar `AgEFF = 0.0` for monovalent (`coupled_percell_pd.py:318,332`) |

---

## 3. PARAMETERS OWNED

**"Owned" = defined inside `kinetic_rhoden_percell.py`.** Kinetic rates (`kon`, `koff`, `kDEG`, `kTMD`) are
**arguments**, owned by `pd_model_config.py` / `eng_params_normalized.json` / the callers, and are documented in
their own subsystem docs — they are listed in §1.2 for traceability only, and are **not** tagged here.

| Symbol | Value | Units | Line | Provenance tag | Source | Mechanistic rationale (why this value) |
|---|---|---|---|---|---|---|
| `AVO` | 6.02214076e23 | /mol | `:35` | **[EXACT — SI defining constant, not measured]** | Avogadro constant, exact by definition (2019 SI redefinition). *Deliberately NOT tagged `[MEASURED]`: it is a defined value with zero uncertainty, not the result of a measurement.* | Not a fitted quantity. Consistent with `wholebody_percell.py:29`, `multiarm_binding.py:24`. |
| `r_cell_um` (default) | **8.0** | µm | `:37` | **[ASSUMED: generic cell radius]** | **No in-code citation.** Matches `wholebody_percell.R_CELL_UM = 8.0` (`:33`) and the hard-coded `8.0` at the heme/blood call sites (`coupled_percell_pd.py:318,332`) | Sets `SA_cell = 804.2 µm²`, the surface over which receptor copies are spread → the receptor **surface density** σ, which is what `AgEFF` actually depends on. A generic ~8 µm radius is a reasonable lymphocyte/blast/carcinoma-cell value but it is **one number applied to every cell type in the model** — see §5.6. |
| `span_nm` (default) | **12.5** | nm | `:37` | **[ASSUMED — no source in this file]** | **No in-code citation.** Sibling `multiarm_binding.py:26–27` calls it *"arm reach when a span is unset (compact within-module reach, nm)"* (`DEFAULT_ARM_REACH_NM = 12.5`); `pd_model_config.py:38–39` sets `span_bridge_nm = span_cis_nm = 12.5` with the comment *"AF3/format override per construct"*. | The arm-to-arm reach of the construct. `AgEFF ∝ 1/span`, so this value sets the **absolute scale of all avidity** in the model. 12.5 nm is the right order for an IgG Fab–Fab hemisphere radius, **but nothing in this codebase sources it** — do not let its plausibility launder it into a measured value. **[UNSOURCED — TBD]** for a citable origin. |
| span floor | 1e-3 | nm | `:41` | **[ASSUMED: numerical guard]** | code-internal | prevents `V_Ab = (2/3)πr³ → 0` in a denominator |
| explicit substep tolerance | **0.1** | — | `:66` | **[ASSUMED: numerical tolerance]** | code-internal | max fractional state change per explicit substep (`h·rmax ≤ 0.1`). A conventional forward-Euler accuracy/stability heuristic; not derived, not fitted to data. |
| explicit `nsub` cap | **20000** | count | `:66` | **[ASSUMED: wall-clock guard]** | code-internal | Ceiling to stop a run from hanging. **It is reached at ≥1e5 receptors/cell (§2 EQ-5 RUN-VERIFY), at which point the 0.1 tolerance is silently violated.** |
| BE substep tolerance | **0.25** | — | `:108` | **[ASSUMED: numerical tolerance]** | code-internal, comment *"SLOW scales only (BE is stable)"* | max fractional change per BE substep **on the slow scales only** — deliberately looser than 0.1 because A-stability removes the stiff-mode constraint. §5.2 shows 0.25 is **too loose for the frozen-`x` nonlinearity** at high free drug. |
| BE `nsub` cap | **16** | count | `:108` | **[ASSUMED: wall-clock guard]** | code-internal | Never approached in production (`nsub = 1` at `kDEG=0.5, kTMD=0.9, dt=0.01`). |
| singularity threshold | **1e-12** | — | `:148` | **[ASSUMED: numerical tolerance]** | code-internal | `|det(A)| < 1e-12 ⇒ treat as singular`. **RUN-VERIFY: unreachable for finite rates** (min `|det|` observed over a 10-decade `kx` sweep = 1.07). Not scale-invariant (§5.5). |
| zero-receptor threshold | **1e-30** | nM | `:58,:59,:105` | **[ASSUMED: numerical guard]** | code-internal | below this, `a_r := 0` → crosslink fail-closed |
| dt floor | **1e-30** | days | `:80,:162` | **[ASSUMED: numerical guard]** | code-internal | flux average |
| `2.0` on `kTMD·Bdbl` | 2 | — | `:75,:79,:92,:125,:161` | **[ASSUMED: mechanistic — "each engaged arm carries an independent internalisation hazard"]** | **No source in code.** | Doubly-engaged complexes internalise at `2·kTMD`. **This is a mechanistic assumption with a real numerical consequence** (it doubles the TMDD sink contribution of every crosslinked complex). See §5.4. |
| `2.0` on `koff·Bdbl` | 2 | — | `:90–92`, `:123–125` | **[DERIVED: combinatorics]** | statistical factor | in the **same-antigen** system the two bound arms are identical → either may release → the `Bdbl → BAg1` rate is `2·koff`. This one is not an assumption; it is a correct multiplicity count. (Correctly **absent** from the two-antigen EQ-4, where the arms are distinguishable and the rate is `koff1 + koff2`, `:75`.) |
| census weight on `Bdbl` | 2 (same-antigen) / 1+1 (two-antigen) | copies | `:94` / EQ-4 | **[DERIVED: stoichiometry]** | — | a same-antigen crosslink consumes **2 copies of one** antigen; a heterotypic crosslink consumes **1 copy of each of two**. |

**Parameters this module DOES NOT own but whose values determine its behaviour** (for the reader's line of
sight; provenance belongs to their own docs): `kon` (generic default `1e5 /M/s` → `8.64 /nM/day`,
`pd_model_config.py:33–34`; measured per-molecule values merged from `eng_params_normalized.json`, e.g.
teclistamab `kon_TAA = 1.28e6 /M/s`, `CHANGELOG_2026-07-13.md:236`), `koff` (`= kon·KD` when not measured,
`coupled_percell_pd.py:135,196`), `kDEG` (fallback **0.5 /day**, `coupled_percell_pd.py:133,194` —
**[UNSOURCED — TBD]** in this file's line of sight), `kTMD` (`kint_bridge_perday = 0.9 /day`,
`pd_model_config.py:37`), `NM_PER_COPY = 6.0/257000 = 2.335e-5 nM/copy` (`wholebody_pd.py:83`, a 71 pL
"synapse" reaction volume), `dt = 0.01 d`.

---

## 4. WHAT IS EMERGENT vs IMPOSED

### 4.1 Genuinely emergent (computed from mechanism, not handed in)

| Emergent quantity | Emerges from |
|---|---|
| **Avidity / bivalent enhancement** | Not a parameter anywhere. It is `kon × AgEFF`, where `AgEFF` is **pure geometry** (EQ-1: receptor copies, cell radius, arm span) and `kon` is the **same** on-rate as the first bond. No "avidity factor", no folded apparent-KD, no bivalency multiplier exists in this file. This is the strongest emergence claim the subsystem makes, and it holds. |
| **The hook effect / self-limiting crosslink** | `x = AgEFF·(Ag/Ag_0)` (EQ-3). At saturating drug the free-receptor fraction collapses and crosslink formation shuts down **by itself**. No explicit bell-shape term. |
| **TMDD clearance** | `intern_flux` summed over ~1e5 cells × organs (EQ-7/EQ-13 → `wholebody_percell.py:208`). The plasma sink is **not** a fitted `Vmax/Km`; it is the integral of a per-cell mechanism. A high-density organ clears more because it *binds* more. |
| **Target-recovery kinetics between doses** | `KSYN = Ag_0·kDEG` with explicit `−kDEG·Ag` (EQ-2). Receptor rebound has a real timescale; the C=0 fixed point `Ag → Ag_0` is exact (RUN-verified to 5 digits). |
| **Per-cell heterogeneity of everything above** | Every rate is a per-cell array. `AgEFF` varies cell-to-cell because `R_copies` does (real Xenium/IHC per-cell copies, `coupled_percell_pk.py:88–93`). A 10× receptor spread produces a 10× `AgEFF` spread and a **non-linear** spread in `Bdbl` — this is not averaged away. |
| **Bivalent > monovalent clearance** | Falls out of `2·kTMD·Bdbl` + higher occupancy. Nothing tells the model that bivalent binders clear faster. |

### 4.2 Imposed (constants handed to the mechanism)

| Imposed | Where | Honest status |
|---|---|---|
| **The arm span, 12.5 nm** | `:37` | The **single scale factor for all avidity** in the model (`AgEFF ∝ 1/span`). **[ASSUMED, UNSOURCED in code.]** Change it to 25 nm and every avidity in the model halves. This is the most load-bearing unsourced number in the subsystem. |
| **Cell radius, 8 µm** | `:37` | One value for T cells, myeloma blasts, hepatocytes, carcinoma cells. `AgEFF ∝ 1/r_cell²`. **[ASSUMED.]** |
| **The hemisphere geometry itself** | `:42` | `SA_Ab = πr²`, `V_Ab = (2/3)πr³` — a **uniformly-explored hemisphere**. The real tethered-arm distribution is not uniform (it is a worm-like-chain / Gaussian-tether density that peaks well inside the maximum reach). This is a deliberate simplification, and it is the *same* one Rhoden's construction makes. |
| **`kTMD` factor of 2 for `Bdbl`** | `:75,:92` | **[ASSUMED.]** Doubles the TMDD contribution of every crosslink. §5.4. |
| **The receptor is immobile / well-mixed on the surface** | (implicit) | `AgEFF` is a *mean-field* local concentration. There is **no** lateral diffusion, no receptor clustering, no lipid-raft partitioning. A clustered receptor would present a far higher local `AgEFF` than the mean. **[ASSUMED.]** |
| **Bound receptor does not degrade at `kDEG`** | `:71–75, :90–92` | Only **free** `Ag` carries `−kDEG·Ag`. Bound receptor's only exit is `kTMD`. A defensible convention (the drug protects/redirects the receptor) but it is a convention. |
| **`kon` is identical for the first bond (3D, from bulk) and the second bond (2D, tethered)** | `:71–75, :90–92, :118` | This is the load-bearing assumption *underneath* the emergence claim in §4.1. The second bond's on-rate is taken to be the same intrinsic `kon` acting on a higher concentration. If the tethered arm has a different intrinsic reactivity (orientation constraints, strain), that difference is **not modelled**. |

### 4.3 Where the emergence stops (the honest boundary)

The subsystem computes avidity from geometry — but the **geometry is imposed** (12.5 nm, 8 µm, hemisphere,
uniform density, immobile receptors). So the correct claim is: *avidity is emergent **given** a specified
molecular span and receptor density; it is not fitted to any binding or PK observable.* That is a strong and
defensible claim. The claim it is **not** entitled to is that the geometric inputs themselves are measured —
in this file, they are not sourced.

---

## 5. KNOWN LIMITATIONS / OPEN QUESTIONS

*(These are the attacks a reviewer will make. All five numbered items were RUN-verified this session; scripts in
`/Volumes/T7_SSD/claude-tmp/2026-07-13/`.)*

### 5.1 The explicit solver's `nsub` cap silently voids its own stability criterion — and it is dead code anyway

At ≥1e5 receptors/cell, `nsub` saturates the 20,000 cap (`:66`) and the realised `h·rmax = 0.215` — **2.15× the
0.1 tolerance the substep law exists to enforce**. At 5e5 and 1e6 copies it is 1.07 and 2.14 (RUN-VERIFY, EQ-5).
Beyond the cap, the solver is an under-substepped forward Euler with a mass-injecting clamp, which is exactly
the failure mode its own comment (`:61–63`) warns about.

**An in-repo reassurance on this point does NOT hold up, and I am flagging it rather than repeating it.**
`CHANGELOG_2026-07-13.md:191–196` records testing `rhoden_bivalent_step` at CD20 2.5e5 / 1e6 / EGFR-ceiling 2.3e6
copies at 50 nM and reports `nsub = 271 / 1073 / 2465`, *"ALL well under the 20,000 cap"*, concluding the
explicit PK bivalent path *"is safe for the submission panel."* **Those numbers are not reproducible from the
module's own substep law.** Driving `:64–66` directly with `geo_ageff_nM` at exactly those copy numbers, generic
`kon = 8.64 /nM/day`, `C = 50 nM`, `dt = 0.01 d` gives uncapped `nsub = 107,124 / 428,228 / 984,808` — **all three
cap out**, and are ~400× larger than the logged values. (The `271:1073:2465` ratios *are* linear in copy number,
so the logged run did obey a `kon·AgEFF`-dominated law — but with an `AgEFF` ~400× below the geometric one, i.e.
that test cannot have been passing the module's own `geo_ageff_nM` at those densities.) I could not reconstruct
what it *was* passing. **Do not rely on the "well under the cap" conclusion:** the doc's own EQ-5 RUN-VERIFY
shows the cap is already hit at **1e5** copies — a *lower* density than the lowest one that note declares safe.
This is unresolved, and it is flagged as such. **[UNRESOLVED — in-repo claim contradicted by RUN-VERIFY]**

**Mitigation: none needed — the function has no live caller (§1.4, exhaustive grep this session).** The correct
action is to **delete it or mark it `_reference_only`**; leaving an unguarded, cap-violating solver in a live
module is an invitation for someone to wire it back in.

### 5.2 ⚠ THE MATERIAL ONE — the BE solver's non-negativity clamp manufactures receptor census and inflates the TMDD sink at high free drug

The docstring (`:96–97`) promises *"census exactly preserved (no manufacture) even at large AgEFF."* The linear
solve delivers that (EQ-11, verified to machine precision). **The function does not**, because `:160` clamps.

**RUN-VERIFY — production configuration** (`nsub = 1`, `h = dt = 0.01 d`, `kon = 8.64 /nM/day`, `KD = 1.45 nM`,
`kDEG = 0.5`, `kTMD = 0.9`, marched 5 simulated days), benchmarked against a **tight-tolerance reference
integration of the identical 3-state ODE** (`scipy Radau`, `rtol 1e-10`, `atol 1e-14`):

| receptors/cell | free drug C (nM) | `Bdbl` (BE) | `Bdbl` (reference) | internalised drug (BE) | (reference) | **flux error** | **census injected** |
|---|---|---|---|---|---|---|---|
| 1e4 | 1 | 0.0310 | 0.0310 | 0.3620 | 0.3620 | **−0.00 %** | 0.0 % |
| 1e4 | 10 | 0.0680 | 0.0308 | 0.7889 | 0.3866 | **+104 %** | **+375 %** |
| 1e4 | 50 | 0.0393 | 0.0287 | 0.5873 | 0.4084 | **+44 %** | **+161 %** |
| 5e4 | 10 | 0.4341 | 0.1591 | 4.1378 | 1.8938 | **+118 %** | **+432 %** |
| 1e5 | 1 | 0.3146 | 0.3146 | 3.6232 | 3.6232 | **−0.00 %** | 0.0 % |
| 1e5 | 10 | 1.1378 | 0.3204 | 8.4323 | 3.7709 | **+124 %** | **+450 %** |
| 1e5 | 50 | 0.4860 | 0.3147 | 6.4100 | 3.8437 | **+67 %** | **+226 %** |
| 3e5 | 50 | 1.4805 | 0.9597 | 19.2468 | 11.3936 | **+69 %** | **+229 %** |

("census injected" = cumulative mass added by `max(y,0)` over 5 d, as % of `Ag_0`; computed by re-running the
identical linear solve without the clamp and differencing the weighted census.)

**Mechanism.** `M` is not Metzler (`M[0,1] = koff − kx < 0` whenever `kx > koff`, i.e. always at real densities
— EQ-10). At `h·kx ~ 1e3` the frozen-`x` BE step can therefore drive the free-`Ag` row negative; `max(y,0)`
zeroes it, and the zeroed deficit re-appears as receptor census. In the 1e4-copy / C=10 nM case the pre-clamp
state goes negative on **241 of 500 steps**.

**It is NOT a stability failure and NOT a stiffness failure — it is the frozen-coefficient splitting error, and
it is cheap to fix.** Refining the substep:

| forced `nsub` | h (days) | `Bdbl` | reference | flux error |
|---|---|---|---|---|
| **1 (production)** | 0.010000 | 1.13784 | 0.32041 | **+123.62 %** |
| 4 | 0.002500 | 0.32041 | 0.32041 | **−0.00 %** |
| 16 | 0.000625 | 0.32041 | 0.32041 | −0.00 % |
| 1024 | 0.000010 | 0.32041 | 0.32041 | +0.00 % |

**A 4× finer substep removes the entire error.** The `nsub` required for <1 % flux error across a
(receptor × free-drug) sweep is **1–128** — versus the **4,300–20,000** the *explicit* solver needs merely to be
stable. The BE architecture is right; **only its substep law (`:108`) is too lax.** Two clean fixes, in order of
preference:
1. **Newton-iterate `x`** within the BE step (2–3 iterations on the 3×3 system) — removes the splitting error at
   its source and keeps `nsub = 1`.
2. Add a term to `:108` that senses the nonlinearity (e.g. require `h·kx·(Ag/Ag_0) ≲ O(10²)`, or simply raise the
   cap and tighten the tolerance) — cheaper to implement, cruder.

**⚠ The error is NOT monotone in receptor density — do not use a simple threshold rule.** RUN-VERIFY (same
script): at `rec = 3e5`, `C = 10 nM` the production step is **exact (−0.00 %)**, while at `rec = 1e5`, `C = 10 nM`
it is **+123.6 %**. The frozen-`x` splitting error depends on where `h·kx` lands relative to the step's own
relaxation, and that is non-monotone in `Ag_0`. A low-`nsub` C-sweep at one receptor density therefore proves
nothing about another. The verified `C`-dependence at `rec = 1e5` is: `C ≤ 7 nM` → **|err| < 0.01 %**;
`C = 10 nM` → **+123.6 %** — the transition is essentially a cliff, not a ramp.

**Blast radius — who is actually in the bad regime?** The error is negligible below ~5 nM free drug and severe
above ~10 nM **at most (but not all) densities — see the non-monotonicity warning above**. The **organ path**
sees *interstitial* drug (a fraction of plasma) and is often safe. The **heme
and blood TMDD sinks are passed plasma `C_pl` DIRECTLY** (`coupled_percell_pd.py:320,334`), and TCE plasma Cmax
is routinely 10–100 nM — **so the two compartments most exposed to this error are exactly the ones that carry
the heme validation panel** (teclistamab, elranatamab, mosunetuzumab, epcoritamab, glofitamab). **Any TMDD /
clearance number from those compartments should be treated as provisional until the substep law is fixed and the
run repeated.** (Note the sign: the error **inflates** the sink, so it inflates target-mediated clearance and
`Bdbl`-driven readouts.)

### 5.3 The singular-row guard's stated justification is a non-sequitur (the guard itself is fine)

See EQ-12. The `:135–136` argument ("diagonal ≥ 1 ⇒ cannot be singular") does not follow. The **correct**
justification is spectral (M's eigenvalues are strictly negative-real ⇒ `1/h > 0` is never an eigenvalue ⇒ A is
never singular for `h > 0`), and I verified it over a 10-decade `kx` sweep. Consequences: (i) the
`|det| < 1e-12` branch (`:148`) is **dead** for finite input — the guard is, in practice, a **non-finite-input
guard**; (ii) the comment's secondary blame on *"Ag_0=0 giving a degenerate a_r"* (`:137`) is **wrong** — `a_r`
is zero-guarded at `:105`, and a zero-receptor cell integrates cleanly. Neither error changes the behaviour, but
both would mislead the next person to debug this path. **The comment should be corrected, not the code.**

Secondary: the warning fires **once per process** (`_warned` attribute, `:151–156`). If the first bad cell is
benign and a genuine divergence appears at hour 3, **it is frozen and never reported**. The count is not
accumulated, not returned, and not surfaced to the caller. Recommend: return `nbad` (or accumulate it on the
function/caller) so a run can be *audited* rather than merely *survived*. Also note `rhs[sing] =
nan_to_num(y0[sing])` (`:158`) **zeroes** a NaN state rather than freezing it — if a cell's *state* (not just its
input) has already gone NaN, its receptors are silently deleted, not preserved.

### 5.4 The `2·kTMD` for crosslinked complexes is an unsourced mechanistic assumption with a first-order effect

Nothing in the file sources it. It says: two engaged receptors ⇒ double the endocytic hazard. The alternative
convention (one complex, one internalisation event at `kTMD`, carrying two receptors) is equally arguable and
would **halve** the crosslink contribution to the TMDD sink. Because the sink is emergent, this single factor of
2 propagates directly into plasma clearance for every bivalent construct. A reviewer will ask for it; there is
currently no answer. **[UNSOURCED — TBD]** — it needs either a literature anchor (receptor-crosslinking →
internalisation-rate data) or an explicit sensitivity analysis.

### 5.5 Minor numerics

- **`|det| < 1e-12` is not scale-invariant** (`:148`). `det(A)` scales as the product of the three diagonal
  entries, which at large `h·kx` is ~1e3–1e4; a *relatively* near-singular matrix would never trip an
  *absolute* 1e-12 threshold. A condition-number or scaled-determinant test would be the correct construction —
  though on present evidence (§5.3) the branch is unreachable anyway.
- **The flux accumulators (`:79`, `:161`) integrate the post-clamp state**, so clamp-injected mass is counted
  into the TMDD sink rather than discarded. This is what converts the census artefact of §5.2 into a *PK*
  artefact.
- **`np.linalg.det` on the whole (n,3,3) batch every substep** (`:145`/`:147`) costs a second full batched
  factorisation purely to decide whether to do the first. For a 3×3 it is cheap in absolute terms, but it is
  ~30–40 % of the BE step's arithmetic and is entirely avoidable (the finiteness mask alone is sufficient,
  given §5.3).

### 5.6 Structural / biological open questions

1. **One cell radius (8 µm) for every cell type.** `AgEFF ∝ 1/r_cell²`. A 4 µm T cell presents its tethered arm
   **4× the effective partner concentration** of an 8 µm blast at equal copy number. The heme/blood call sites
   hard-code `8.0` (`coupled_percell_pd.py:318,332`) — so T-cell-side cis avidity, if ever routed through this
   module, would be systematically under-estimated. (`multiarm_binding.multiarm_bound` does carry separate
   `r_cell_T = 4.0` / `r_cell_tum = 8.0`, so the engine *knows* the distinction elsewhere.)
2. **Mean-field receptor distribution.** No clustering, no lateral diffusion, no co-localisation. Clustered
   receptors would give a local `AgEFF` far above the surface-average — the model cannot express this, and it is
   plausibly a major source of real bivalent avidity.
3. **The uniform hemisphere.** A tethered Fab's reachable-partner density is not uniform over the hemisphere.
   A worm-like-chain / Gaussian-tether kernel would concentrate reach at intermediate radii and change `AgEFF`
   by an O(1) factor.
4. **The two-antigen heterotypic path is written but not wired** (§1.4). The docstring advertises three variants
   (cis-T `CD3×costim`, cis-tumour `TAA×TAA`, trans-bridge `CD3×TAA`, `:29–32`) but the only live solver is the
   **same-antigen** reduction. The cross-cell bridge (the thing that actually kills tumour cells) is computed
   **elsewhere** (`kinetic_synapse.py` / `wholebody_pd.py`), on a different code path with a different solver.
   The docstring's *"ONE binding solve, used IDENTICALLY by PK and PD"* (`:13`) is therefore **aspirational, not
   current**: this module is the PK/TMDD binding solve and the same-antigen avidity solve. A reviewer reading the
   docstring will believe more unification than exists.
5. **`kDEG` fallback 0.5 /day** (`coupled_percell_pd.py:133,194`) is applied to whatever receptor the molecule
   targets — BCMA, CD20, CD19, GPRC5D all get the same turnover unless overridden. **[UNSOURCED — TBD]** from
   this module's line of sight; it belongs to the receptor-parameter subsystem, but it directly sets `KSYN` here.

### 5.7 What I could NOT verify

- The **Rhoden et al. 2016** attribution for the crosslink structure (`:3–5`, bioRxiv 10.1101/2022.09.12.507653).
  I did not fetch it. The equation *form* is standard tethered-effective-concentration geometry; the
  attribution is the comment's claim. **This is the ONLY external citation in this document, and it is present
  verbatim in the code at `:4–5` — it is not introduced by this doc.** **[UNVERIFIED CITATION]**
- The **"user-supplied MATLAB scheme"** (`:7–8`) that the turnover terms were transcribed from. It is not in
  this repository's line of sight. **[UNSOURCED — TBD]**
- The origin of **12.5 nm** and **8.0 µm**. No citation exists anywhere in the live module or its callers.
- **The `nsub = 271 / 1073 / 2465` figures logged at `CHANGELOG_2026-07-13.md:191–196`.** They are ~400× below
  what the module's own substep law produces at the stated copy numbers and are **not reproducible** (§5.1). I
  could not determine what `AgEFF` that test actually passed. The "safe for the submission panel" conclusion
  drawn there is **not supported** and should not be cited. **[UNRESOLVED]**

---

## 6. ADVERSARIAL VERIFICATION RECORD (2026-07-13)

This document was re-checked line-by-line against the live source under an explicit find-the-errors mandate.
What was **confirmed by re-running the code** (scripts: `T3_adversarial_verify.py`, `T3_52_verify.py`,
`T3_gaps.py`, all in `/Volumes/T7_SSD/claude-tmp/2026-07-13/`):

- Every equation EQ-1 … EQ-14 appears at the cited line of `kinetic_rhoden_percell.py` (file confirmed 162
  lines / 10,569 bytes / mtime 2026-07-13 14:44; not a git repo).
- **EQ-1 table reproduces exactly** (247.7653 / 2,477.6531 / 24,776.5313 / 123,882.6566 nM; slope 0.2477653
  nM/copy; span law exactly 2.000000).
- **EQ-5 table reproduces exactly** (rmax 44,569 / 429,894 / 2,142,448 / 4,283,140; `nsub` 4,457 / cap / cap /
  cap; realised `h·rmax` 0.1000 / 0.2149 / 1.0712 / 2.1416).
- **EQ-12 reproduces exactly** (0 positive-real eigenvalues over a 10-decade `kx` sweep; min `|det(A)|` =
  1.0736; spectrum at `kx = 1e8` = {−1.0000e+08, −2.5000, −1.8000}).
- **§5.2 — the headline defect — reproduces to every printed digit**, including the negative-step counts
  (241/500 at 1e4 copies / C = 10 nM) and the substep-refinement table (`nsub` 1 → +123.62 %, `nsub` 4 →
  −0.00 %). The `nsub`-for-<1 %-error range of **1–128** is confirmed (worst case 128 at 3e5 copies / C = 100 nM).
- **`rhoden_bivalent_step` has zero live call sites** — confirmed by exhaustive repo-wide grep.
- **No fabricated citations.** The doc introduces no PMID, DOI or paper that is not already in the source file.

**Errors found and corrected in this pass** (all were in the *commentary*, none in the equation transcriptions):
1. A quotation (*"Switching also FIXES the inert-`Bdbl` avidity bug…"*) was attributed to
   `CHANGELOG_2026-07-13.md:243–244`; it actually lives at `reference_unified_binding/LIVE_ENGINE_FIXES.md:57–58`.
   The quote was real but **the source was wrong** — corrected in EQ-8(c).
2. `CHANGELOG_2026-07-13.md:216–220` was cited (twice) for teclistamab's `kon_TAA = 1.28e6 /M/s`; the correct
   line is **`:236`** (corroborated at `:242`). Corrected in §1.5 and §3.
3. §5.1 asserted the changelog's `nsub = 271/1073/2465` were "consistent with the generic kon" — **they are not**,
   and that assertion contradicted this doc's own EQ-5 table. Rewritten as an open, flagged discrepancy.
4. The top-of-doc callout claimed a flux inflation of "+17 % to +130 %"; **neither bound is supported** by the
   §5.2 run. Corrected to the measured **+44 % to +124 %**.
5. The §1.5 timing (36.4 s / 0.045 s → "802×") was internally inconsistent (36.4/0.045 = 809). Replaced with
   this session's re-measurement (29.8 s / 0.039 s → 770×) and explicitly de-rated to a machine-dependent
   order-of-magnitude corroboration rather than a reproduction.
6. `AVO` was tagged **[MEASURED]**; it is a *defined* SI constant with zero uncertainty. Retagged **[EXACT]**.
7. Off-by-one line cites corrected: `coupled_percell_pk.py:73 → :74` (PD_DT), `multiarm_binding.py:27 → :26–27`,
   and the `"NON-FINITE rate leaked in"` comment `:137 → :136`.
8. §5.2's "negligible below 5 nM / severe above 10 nM" rule was stated as if monotone in receptor density. It
   is **not** (at 3e5 copies / C = 10 nM the step is exact). Caveat added.
