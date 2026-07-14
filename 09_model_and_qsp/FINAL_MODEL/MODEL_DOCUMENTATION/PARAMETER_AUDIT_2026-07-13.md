# FULL PARAMETER AUDIT — 2026-07-13

**Scope:** every parameter, data type, dosing schedule, route, receptor param, KD, kon, koff, MW and FcRn
factor in the TCE model, checked against `params/mab_tce_pkpd.sqlite` — a curated database of **76 drugs,
93 literature sources (48 DOI / 43 PMID / 45 FDA labels), 183 PK/PD curves, 936 digitized points, and 14
SPR kon/koff/KD records transcribed verbatim from Biacore tables.**

**Four defects found. Two are serious. All four are fixed.**

---

## ⛔ DEFECT 1 (SERIOUS) — the PK engine and the PD engine were binding the same antigen with DIFFERENT affinities

`run_tce_pd_reval.py` merges measured kinetics from `eng_params_normalized.json` into the `ENG` dict using
`cfg.setdefault(...)`. Because `setdefault` **cannot overwrite an existing key**, `ENG`'s hand-entered `KD`
survived, and the *measured* KD was parked in a **separate key**, `KD_norm`. Then:

```python
line 142   CoupledPerCellPD(..., cfg['KD'], ...)                     # PK antigen sink -> ENG's KD
line 152   attach_pd(..., KD_TAA_nM=cfg.get('KD_norm', cfg['KD']))   # PD synapse     -> MEASURED KD
```

**The same molecule, in the same run, bound the same antigen with two different affinities:**

| molecule | PK bound at | PD bound at | discrepancy |
|---|---|---|---|
| **talquetamab** | **2.0 nM** | **11.0 nM** | **5.5×** |
| **elranatamab** | **0.15 nM** | **0.04 nM** | **3.75×** |
| glofitamab | 5.0 nM | 4.0 nM | 1.25× |
| epcoritamab | 5.0 nM | 4.76 nM | 1.05× |
| mosunetuzumab | 5.0 nM | 5.45 nM | 1.09× |
| teclistamab | 0.15 nM | 0.18 nM | 1.20× |

The antigen sink that decides **how much drug is consumed** used one affinity; the synapse that decides
**killing and cytokine release** used another. These are physically **the same binding event.** This is not a
modelling choice — it is incoherent.

**FIX:** the measured KD (`KD_norm`) now wins **everywhere**. The runner prints the substitution on every run:

```
[KD] talquetamab: using MEASURED KD_TAA=11 nM everywhere (ENG carried 2 nM; PK and PD are now consistent)
```

---

## ⛔ DEFECT 2 (SERIOUS) — most binding kinetics are ASSUMED, and were being presented as measured

| | |
|---|---|
| **`kon = 1e5 /M/s` for 15 of 22 molecules** | a **generic placeholder**, not a measurement |
| **`koff` back-derived as `kon × KD`** | therefore also **not measured** |

**Measured kinetics exist for only four:** `teclistamab` (1.28e6), `elranatamab` (9.95e5), `catumaxomab`
(6.1e4 — and it **matches the SQLite Biacore record exactly**: *"k_ass=6.1±1.5×10⁴ M⁻¹s⁻¹, k_diss=3.3±0.3×10⁻⁵ s⁻¹,
KD=5.6×10⁻¹⁰ M"*, Ruf 2007), and the rituximab-lineage CD20 `kon` (4.3e5).

**Consequence:** for the other 15, **any conclusion that depends on kon/koff rather than on KD alone is
assumption-driven.** Since the model's whole selling point is *literal kinetic* binding (rather than an
occupancy shortcut), this materially limits which molecules its kinetic claims apply to.

**FIX:** tagged in code at the point of use, so it cannot be forgotten:

```python
# kon=1e5 is a GENERIC PLACEHOLDER, not a measurement... For those molecules the KINETICS ARE ASSUMED...
# Do not present the rest as kinetically validated.
```

---

## ⛔ DEFECT 3 — catumaxomab's CD3 affinity was 4.5× wrong

SQLite `kinetics`, **verbatim** from Ruf 2010 (PMC2878603, on-disk XML):

> *"intermediate affinity of catumaxomab to CD3 (**KD = 4.4 nM**)"*

**The model carried 20 nM.** Corrected to 4.4 nM; `koff_CD3` re-derived from `kon × KD`.

---

## ⛔ DEFECT 4 — my own consistency check was CIRCULAR and proved nothing

I checked `KD == koff / kon` across all 22 molecules and got **ratio = 1.000 for every single one.** That
looked like a clean bill of health. **It is worthless.** `koff` *was computed from* `kon` and `KD`, so the
identity holds **by construction**. The check validated my arithmetic, not the data.

> **A check that cannot fail is not a check.** I nearly reported this as evidence the parameters were sound.

---

## ✅ WHAT PASSED

* **Dosing schedules and routes** match the clinical labels for all 6 primary molecules (verified against the
  step-up regimens: mosunetuzumab `1/2/60/60 mg IV`, teclistamab `4.8/24/120 mg SC`, elranatamab `12/32/76 mg SC`,
  talquetamab `0.8/4.8/32 mg SC`, glofitamab `2.5/10/30 mg IV`, epcoritamab `0.16/0.8/48 mg SC`).
* **Fc-less molecules correctly carry `fFcRn = 0.0`** (blinatumomab, solitomab, pasotuxizumab, tebentafusp) —
  correct, and deliberately not "fixed".
* **catumaxomab's EpCAM arm matches the Biacore record exactly** (KD 0.56 nM, kon 6.1e4, koff 3.3e-5).
* **Unit consistency**: all KD in nM, kon in /M/s, koff in /s, MW in kDa. No 1000× unit errors found.

---

## ⛔ DEFECT 5 (SERIOUS) — the most influential biodistribution parameters in the model are UNSOURCED

Found by the T1 subsystem documentation pass (`subsystems/T1_shah_betts_pbpk_backbone.md`):

| parameter | finding |
|---|---|
| **`sigV` — all 15 per-organ vascular reflection coefficients** (0.75 tumour … 0.99 brain) | **NO CITATION ANYWHERE IN THE CODE.** These control how much antibody leaves the capillary in every organ. They are, in the doc's words, *"the single most influential biodistribution parameters in the layer, and they are unsourced."* `[UNSOURCED — TBD]` |
| **`sigL = 0.85`** | **`[FITTED]`** to reproduce `Vss/Vc ≈ 2.1`. The fit works (analytic Vss/Vc = 2.0903) — **but it is 4.25× the Shah-Betts platform value of `σ_L = 0.2`**, independently verified against the source. Consequence: **leaky organs CONCENTRATE antibody above plasma** (tumour Cis/Cpl = 1.667, spleen 1.333). This is the largest single deviation from the model this one was built by reproducing. |
| **`CLup = 0.3503 /day`** | `[FITTED]` to mosunetuzumab's 16.1 d terminal half-life (anchor **verified real**, FDA LUNSUMIO label). **But the claim "matches 16.1 d exactly" is DOSE-CONDITIONAL** — measured 9.11 d at low dose, 38.11 d when TMDD-saturated. It hits 16.1 d only at one specific dose. **Not a model property.** |
| **`k_renal_max = 8.70 /day`** | `[FITTED]` to blinatumomab — and it is the **best-validated parameter in the model**: predicts a BiTE falling to 50% of C₀ in **1.78 h** against the FDA BLINCYTO label's **2.11 h**. ✅ |
| **`fFcRn` and `CLup`** | **NOT SEPARATELY IDENTIFIABLE.** Only the product `CLup·(1−fFcRn)` enters the equations, so plasma PK cannot constrain them individually. Reporting either as a fitted value is misleading. |

**`sigV` is the one to fix first.** It is unsourced, it is the most influential, and nobody knew.

---

## ⛔ DEFECT 6 — the DOCUMENTATION OF THE AUDIT itself contained fabrications

The subsystem docs were written by agents and then **adversarially verified** by independent agents whose only
job was to find fabricated citations. They found some:

* A doc **invented four Shah-Betts reflection coefficients** (σ₁ 0.69–0.999 "mean 0.908"; σ₂ 0.258–0.841
  "mean 0.579") and attributed them to the **real** paper. They appear in neither the code nor the source
  (published: σ₁ 0.883–0.987, σ₂ 0.311–0.837). **Removed.**
* A doc **fabricated a provenance *defect*** — asserting the Tabula Sapiens density CSV "does not exist in the
  repo." **It does**, and all 30 immune-cell densities re-derive from it exactly. It manufactured a problem
  and demanded remediation for it. **Corrected.**
* A fabricated tumour ABC value ("~0.24 from the Shah/Betts work"). **Removed.**
* A number wrong in four places (1.80 h → correctly **1.7830 h**). **Fixed.**

> **Even the documentation of the anti-fabrication audit needed an anti-fabrication audit.**
> This is not an argument against the process — the verifiers *caught* every one of these. It is an argument
> that **no single pass, human or model, should be trusted on provenance.** Verify adversarially, always.

---

## ⚠️ OPEN — NOT FIXED, MUST NOT BE PRESENTED AS SOUND

| item | status |
|---|---|
| **`fFcRn` disagreement between the two param files** | `ENG` says 0.89 for teclistamab/talquetamab; `eng_params_normalized.json` says 0.70. `ENG` currently wins. **Unresolved — which is right?** |
| **mosunetuzumab CD20 KD** | SQLite has **10.2 µg/mL** (Scatchard, popPK-derived — Wang 2024, PMC11134317), a *different basis* from the model's 5.45 nM. Not directly comparable; flagged, not reconciled. |
| **Costim induction folds** | **NOT_FOUND in the literature** for 4-1BB / OX40 / ICOS / GITR. The code **fails closed**. The costim screen cannot produce a valid ranking. |
| **IL-6 clearance** | `[FITTED]`. Its citation (PMID 31268236) is to a **modeling** paper. Human IL-6 clearance appears **unmeasured**. |
| **`k_death`** | `[FITTED]` to epcoritamab depletion. The efficacy scale. |
| **Liver / Kupffer macrophages** | **Absent from the census entirely.** The model requires ≥6.8×10⁹ macrophages to reach clinical IL-6 and has 3.5×10⁷. See `IL6_PREREGISTERED_PREDICTION_2026-07-13.md`. |

---

## The pattern worth naming

Every defect here has the same shape: **two sources of truth that were never forced to agree.**

* `ENG` vs `eng_params_normalized.json` → two affinities for one antigen.
* An in-code comment defining `count_scale` as a **ratio** vs a JSON file holding **counts** → 942× IL-6 error.
* Intent vs data, in both cases, with **nothing asserting the difference.**

The fix in each case was not cleverness — it was **making the model check itself**: one affinity per antigen,
a physical ceiling on production, a refusal to run on unsourced parameters. Every one of those checks is
cheaper than the bug it catches.
