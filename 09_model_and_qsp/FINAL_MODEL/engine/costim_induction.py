"""costim_induction.py — ACTIVATION-INDUCED costim receptor expression.

THE PROBLEM THIS FIXES (found 2026-07-13, Max's catch):
  The model read costim receptor density ONCE from resting healthy-tissue scRNA-seq and never changed it.
  But 4-1BB (TNFRSF9), OX40 (TNFRSF4) and ICOS are ACTIVATION-INDUCED -- they are essentially ABSENT on
  resting T cells and appear only AFTER TCR engagement. That is precisely WHY the field targets them: the
  costim fires only on T cells that have ALREADY engaged a target, i.e. it is TUMOUR-CONDITIONAL.

  With a static resting density the model sees them at ~zero and systematically UNDER-RATES exactly the arms
  that matter, while OVER-RATING constitutive arms (CD2, CD28, CD27). A screen run on that model would have
  produced a confidently-wrong ranking ("CD2 beats 4-1BB") that any immunologist would reject on sight.

THE MECHANISM (this turns the bug into the actual biology):
    da_i/dt   = k_on * p_eng_i * (1 - a_i)  -  k_off * a_i        # per-T-cell activation memory
    R_i(t)    = R_rest_i * (1 + (FOLD - 1) * a_i)                 # induced receptor density

  p_eng_i = that T cell's OWN engaged-synapse fraction (B2/RC) -- so induction is DRIVEN BY ENGAGEMENT,
  spatially, per cell. A T cell that never engages never upregulates: the conditionality is EMERGENT.

  FOLD = 1  -> R(t) == R_rest  -> EXACTLY the previous static behaviour (constitutive arms are unchanged).
  This is strictly additive; the clinical validation molecules (no costim arm) are unaffected.

PARAMETERS: fold-upregulation + induction kinetics per arm are LITERATURE values. Any arm not in the table
falls back to FOLD=1 (constitutive / unknown) and is FLAGGED, never guessed.
"""
import os
import numpy as np

# ---- LITERATURE TABLE ------------------------------------------------------------------------------------
# arm -> dict(fold=<activation fold-upregulation>, t_peak_h=<hours to peak>, t_decay_h=<hours to decay>,
#             kind='inducible'|'constitutive', source='<PMID/DOI>')
#
# THE FOLD-UPREGULATION FOR EVERY INDUCIBLE ARM IS **NOT_FOUND** IN THE LITERATURE (audited 2026-07-13).
# A dedicated sourcing pass found KINETICS but NO clean surface-density fold-change for any of the four.
# What exists is mostly %-POSITIVE data, which is a DIFFERENT QUANTITY and must never be silently converted
# into a density fold. So `fold=None` stands, and the model REFUSES to run these arms rather than guess.
#
# WHY THE REFUSAL MATTERS: 4-1BB / OX40 / ICOS / GITR are ~ABSENT on resting T cells and appear only AFTER
# TCR engagement -- that conditionality is precisely WHY they are targeted. Running them at their RESTING
# density (i.e. quietly assuming fold=1) UNDER-RATES them by the entire induction factor and yields a
# confidently-wrong ranking such as "CD2 beats 4-1BB". A fabricated fold would be worse still.
COSTIM_INDUCTION = {
    # constitutive on resting T cells -> no induction; the static resting density IS correct for these
    "CD28":     dict(fold=1.0, t_peak_h=None, t_decay_h=None, kind="constitutive", source="TBD"),
    "CD2":      dict(fold=1.0, t_peak_h=None, t_decay_h=None, kind="constitutive", source="TBD"),
    "CD27":     dict(fold=1.0, t_peak_h=None, t_decay_h=None, kind="constitutive", source="TBD"),
    # ACTIVATION-INDUCED. fold = NOT_FOUND -> None -> hard refusal (never silently treated as 1.0).
    # Kinetics WERE sourced and are used; only the fold magnitude is missing.
    "TNFRSF9":  dict(fold=None, t_peak_h=24.0, t_decay_h=72.0,  kind="inducible",   # 4-1BB
                     source="kinetics PMID 24634374 (t_peak 24h, t_decay 72h); FOLD = NOT_FOUND"),
    "TNFRSF4":  dict(fold=None, t_peak_h=60.0, t_decay_h=120.0, kind="inducible",   # OX40 (t_peak 48-72h -> 60)
                     source="kinetics PMID 19538134 (t_peak 48-72h, t_decay ~120h); FOLD = NOT_FOUND"),
    "ICOS":     dict(fold=None, t_peak_h=None, t_decay_h=None,  kind="inducible",
                     source="FOLD = NOT_FOUND; kinetics = NOT_FOUND"),
    "TNFRSF18": dict(fold=None, t_peak_h=None, t_decay_h=None,  kind="inducible",   # GITR
                     source="FOLD = NOT_FOUND; kinetics = NOT_FOUND"),
}

# ---- SENSITIVITY SWEEP over the UNKNOWN fold -------------------------------------------------------------
# Since the fold cannot be sourced, the scientifically honest move is NOT to invent one -- it is to ask
# whether the RANKING even depends on it. Set COSTIM_FOLD to run every inducible arm at an ASSUMED fold, and
# sweep it (e.g. 3, 10, 30). If an arm wins at EVERY plausible fold, that conclusion is ROBUST TO THE UNKNOWN,
# which is a stronger claim than one made with a fabricated constant. If the ranking flips inside the range,
# we have identified exactly which measurement the field is missing -- also a result.
#
# The assumed fold is RECORDED on every run (`assumed_fold`, `fold_is_assumed=True`), so no output can ever
# silently inherit an unlabelled assumption. A result produced under this sweep must ALWAYS be reported as
# "conditional on an assumed fold of N", never as a bare ranking.
COSTIM_FOLD_ENV = "COSTIM_FOLD"


class CostimInduction:
    """Per-T-cell activation memory driving activation-induced costim receptor density."""

    def __init__(self, arm, n_T, R_rest, strict=True):
        self.arm = arm
        self.R_rest = np.asarray(R_rest, float)
        self.a = np.zeros(int(n_T))                 # per-cell activation state in [0,1]
        p = COSTIM_INDUCTION.get(arm)
        if p is None:
            self.fold, self.kind, self.source = 1.0, "unknown", "NOT IN TABLE"
            import sys; sys.stderr.write(
                f"[costim] '{arm}' not in the induction table -> treated as CONSTITUTIVE (fold=1). "
                f"If it is activation-induced, this UNDER-RATES it.\n")
        else:
            self.fold, self.kind, self.source = p["fold"], p["kind"], p["source"]
            self.fold_is_assumed = False
            if p["kind"] == "inducible" and p["fold"] is None:
                # SENSITIVITY SWEEP: an assumed fold must be EXPLICITLY opted into via COSTIM_FOLD, and is
                # then RECORDED on the object so no downstream result can inherit it silently. There is no
                # default. Absence of the env var is a REFUSAL, not a fallback to 1.0 -- because a silent
                # fold=1 is precisely the under-rating bug this module exists to fix.
                env = os.environ.get(COSTIM_FOLD_ENV)
                if env is not None:
                    self.fold = float(env)
                    self.fold_is_assumed = True
                    import sys; sys.stderr.write(
                        f"[costim] '{arm}': fold is NOT_FOUND in the literature. Running at an ASSUMED fold "
                        f"of {self.fold:g} ({COSTIM_FOLD_ENV}). THIS RESULT IS CONDITIONAL ON THAT ASSUMPTION "
                        f"and must be reported as such — sweep the fold and check the ranking is robust.\n")
                else:
                    msg = (f"[costim] '{arm}' is ACTIVATION-INDUCED but its fold-upregulation is NOT_FOUND in "
                           f"the literature (audited 2026-07-13: kinetics exist, the fold magnitude does not). "
                           f"Running it at the resting density would UNDER-RATE it by the whole induction "
                           f"factor — the exact bug this module exists to fix — and would produce a "
                           f"confidently-wrong ranking (e.g. 'CD2 beats 4-1BB').")
                    if strict:
                        raise ValueError(
                            msg + f" Either source the fold, or run a SENSITIVITY SWEEP by setting "
                                  f"{COSTIM_FOLD_ENV}=<fold> (e.g. 3, 10, 30) — the assumed value is then "
                                  f"recorded on every result. Do NOT pass strict=False to make this go away.")
                    import sys; sys.stderr.write(msg + " Proceeding at fold=1 (UNDER-RATED).\n")
                    self.fold = 1.0
        # kinetics: k_on from time-to-peak (reach ~95% of max at t_peak), k_off from decay
        tp = (p or {}).get("t_peak_h"); td = (p or {}).get("t_decay_h")
        self.k_on = (3.0 / tp) if tp else 1.0        # /hr
        self.k_off = (1.0 / td) if td else 0.05      # /hr

    def step(self, dt_days, p_eng):
        """Advance each T cell's activation memory from ITS OWN engaged fraction, return induced density."""
        dt_h = dt_days * 24.0
        e = np.clip(np.asarray(p_eng, float), 0.0, 1.0)
        if e.size != self.a.size:
            e = np.resize(e, self.a.size)
        self.a += dt_h * (self.k_on * e * (1.0 - self.a) - self.k_off * self.a)
        np.clip(self.a, 0.0, 1.0, out=self.a)
        return self.R_rest * (1.0 + (self.fold - 1.0) * self.a)

    def summary(self):
        # fold_is_assumed TRAVELS WITH THE RESULT. A ranking produced under an assumed fold must never be
        # reported as if the fold were known — the flag makes that impossible to forget downstream.
        return dict(arm=self.arm, kind=self.kind, fold=self.fold, source=self.source,
                    fold_is_assumed=bool(getattr(self, "fold_is_assumed", False)),
                    mean_activation=float(self.a.mean()) if self.a.size else 0.0,
                    frac_induced=float((self.a > 0.1).mean()) if self.a.size else 0.0)
