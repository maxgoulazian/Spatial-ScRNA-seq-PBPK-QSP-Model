"""pd_model_config.py — THE SINGLE SOURCE OF TRUTH for which PD kill law + parameters a run uses.

Created 2026-07-12 (user: "lets make sure to separate out and make a 1 final model so we dont get
confused on other runs later"). Every whole-body PD/PK-PD harness imports CANONICAL from here instead
of hard-coding engine flags. Change the engine in ONE place; all runs follow.

THREE PD kill laws exist in the codebase — they are NOT interchangeable, do NOT mix them in one run:

  1. 'abstract'  qsp_costim_window_v2.py, pd_driver_mode=None
        The bell-curve window-scoring used for the COSTIM NOMINATION. Byte-frozen; do not touch.
        Use ONLY for the abstract 12-arm nomination / window ranking, never for clinical validation.

  2. 'qss'       wholebody_pd.OrganPD(pd_kinetics=False)
        Schropp ternary EQUILIBRIUM (RC_AB = C*R_A*R_B/(alpha*KD1*KD2)) evaluated per cell per step.
        Prozone emergent, abundance-linear. The prior validated whole-body PD. Fast. No dwell-time /
        serial-killing kinetics (equilibrium magnitude only).

  3. 'kinetic'   wholebody_pd.OrganPD(pd_kinetics=True) -> kinetic_synapse.KineticSynapse   <-- CANONICAL
        Literal per-cell engage/hit/detach bond ODE (exact 2x2 matrix-exponential). Two-sided
        conservation (CD3 AND TAA), avidity-dependent synapse dwell time, explicit serial-killing cycle.
        Captures dwell-time / serial-killing / synapse-stability effects the QSS law structurally cannot
        (the axes an affinity+format design sweep tunes). Anchor = CLINICAL PD (mosunetuzumab IL-6 570
        pg/mL + near-complete B-cell depletion); k_hit fixed from literature, only k_death & IL6 scale set.

The CANONICAL engine for all NEW clinical PK-PD validation and design-sweep runs is 'kinetic'.
"""

# ---- CANONICAL engine selection (change here, everywhere follows) --------------------------------
PD_ENGINE = "kinetic"                 # one of {"kinetic","qss"} for whole-body validation runs

# ---- kinetic-synapse parameters (clinical-anchored; see handoff/kinetic_calib.json) --------------
KINETIC = dict(
    kon_CD3_perM_s = 1e5,             # CD3 arm on-rate (/M/s), standard mAb assoc; koff = kon*KD if not given
    kon_TAA_perM_s = 1e5,             # TAA arm on-rate (/M/s)
    koff_CD3_pers  = None,            # per-molecule koff (/s) override; None -> kon*KD_CD3 (KD-consistent)
    koff_TAA_pers  = None,
    kint_bridge_perday = 0.9,         # trimer internalization (/day)
    span_bridge_nm = 12.5,            # trans CD3<->TAA arm span (nm); AF3/format override per construct
    span_cis_nm    = 12.5,            # cis CD3<->costim arm span (nm) for trispecific avidity
    cis_avidity    = 0.0,             # costim co-engagement avidity [0,1); 0 for plain CD3xTAA validation
    k_hit_perday   = 12.0,            # /day (~1 hit/2h) FIXED from serial-killing literature, not fitted
)

# clinical-anchored calibration constants (mosunetuzumab IL-6 570 pg/mL + depletion)
IL6_SCALE_KIN = None                  # set below from calib
K_HIT         = None
K_DEATH       = None

def _load_calib():
    """Load the kinetic calibration (IL6 scale, k_hit, k_death) from handoff/kinetic_calib.json.
    Kept as a function so a harness can point KWS at its own handoff dir."""
    import json, os
    global IL6_SCALE_KIN, K_HIT, K_DEATH
    here = os.path.dirname(os.path.abspath(__file__))
    for cand in (os.path.join(here,"handoff","kinetic_calib.json"),
                 os.path.join(here,"kinetic_calib.json")):
        if os.path.exists(cand):
            d = json.load(open(cand))
            IL6_SCALE_KIN = d["IL6_SCALE_kin"]; K_HIT = d.get("k_hit",12.0); K_DEATH = d.get("k_death",1.0)
            KINETIC["k_hit_perday"] = K_HIT
            return d
    # fallback literals (documented) if calib file absent
    IL6_SCALE_KIN = 0.05473; K_HIT = 12.0; K_DEATH = 1.0
    KINETIC["k_hit_perday"] = K_HIT
    return dict(IL6_SCALE_kin=IL6_SCALE_KIN, k_hit=K_HIT, k_death=K_DEATH, source="fallback")

def kin_params():
    """kwargs dict for coupled_percell_pd.attach_pd(kin_params=...)."""
    if K_HIT is None: _load_calib()
    return dict(KINETIC)

def is_kinetic():
    return PD_ENGINE == "kinetic"

_load_calib()
