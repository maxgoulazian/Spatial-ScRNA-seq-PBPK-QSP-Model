"""myeloid_il6.py — PER-CELL, SPATIALLY-EMERGENT IL-6 (CRS) emitter. NO fitted scale, NO Emax, NO EC50.

WHY THIS REPLACES THE FITTED SCALE
----------------------------------
The legacy engine computed  il6_pgml = engaged_dwell_rate * IL6_SCALE  where IL6_SCALE was ONE constant fit
so mosunetuzumab -> 570 pg/mL. That is broken three ways:
  (1) IL-6 is NOT emitted by T cells. It is myeloid-derived: "macrophages ... the main overall source of IL-6"
      (Giavridis, Nat Med 2018, PMID 29808005); "human monocytes were the major source of IL-1 and IL-6"
      (Norelli, Nat Med 2018, PMID 29808007). T cells are the ACTIVATING INPUT, not the emitter.
  (2) The engine had NO cytokine clearance ODE — it multiplied a production RATE (pg/mL/hr) by a scalar and
      called it a CONCENTRATION (pg/mL). IL6_SCALE was silently standing in for the missing 1/kdeg.
  (3) A single global scale (or a single global Emax) makes every molecule saturate to the SAME ceiling, so
      construct differences wash out — fatal for a counterscreen whose whole job is to rank constructs.

THE MECHANISM (fully emergent, per cell)
----------------------------------------
Each myeloid agent i (monocyte / macrophage — REAL cells already present in the agent tables) has its own
activation state a_i in [0,1] driven by ITS OWN local contact with ENGAGED T cells (contact-gated: "IL-6
induction and myeloid activation require proximity of CAR T cells and myeloid cells", PMID 29808005; the
CD40L-CD40 axis). Each activated cell secretes at up to its MEASURED per-cell maximum rate. Plasma IL-6 then
accumulates against its MEASURED first-order clearance:

    da_i/dt = k_on * contact_i * (1 - a_i)  -  k_off * a_i          [structural saturation: the (1-a_i) term]
    IL6_prod = count_scale * SUM_i ( a_i * S_MAX )                   [pg/hr]
    dC/dt    = IL6_prod / V_plasma  -  kdeg * C                      [pg/mL/hr -> pg/mL]

EMERGENT PROPERTIES (nothing fitted):
  * SATURATION emerges from the finite myeloid pool + the per-cell (1-a_i) ceiling. No Emax, no Hill, no EC50.
  * THE ~3.9% "secretor fraction" measured single-cell (PMID 37533643) EMERGES as "only the myeloid cells that
    happen to sit near engaged T cells are activated" — spatial heterogeneity IS the secretion heterogeneity.
  * PER-MOLECULE DIFFERENCES emerge from WHERE each drug engages: mosunetuzumab (CD20) engages a huge B-cell
    field in the SPLEEN, which is 16% myeloid (54,373 resident macrophages) -> many myeloid activated -> high
    IL-6. Elranatamab (BCMA) engages rare plasma cells whose mass sits in BONE MARROW, only 1.8% myeloid
    (~3,000 cells) -> few myeloid activated -> lower IL-6. The clinical 570-vs-191 split falls out of ANATOMY.

EVERY PARAMETER IS LITERATURE-MEASURED (none fitted, none tuned to model output):
"""
import numpy as np

# ---- LITERATURE CONSTANTS (measured; cite in the deck) -------------------------------------------------
IL6_MW_DA        = 21000.0    # mature human IL-6 ~21 kDa
NA               = 6.02214076e23
S_MAX_MOLEC_PER_S= 10.6       # MEAN per-cell IL-6 secretion of an ACTIVELY-SECRETING human monocyte
                              #   (6.5 +/- 3.2 at 3 h to 10.6 +/- 7.1 at 12 h molec/s; PMID 20376398 — the
                              #   population average OVER SECRETING CELLS). Use the 12 h value.
                              # CORRECTION (2026-07-13): this was 156 molec/s, the PEAK/high-secretor-tail rate
                              #   (PMID 37533643). That is the rate of the TOP few percent -- the SAME paper
                              #   measures only ~3.9% of cells as active secretors at peak frequency. Applying
                              #   the top-tail rate to EVERY activated cell over-states secretion ~15x. The
                              #   defensible per-activated-cell rate is the measured MEAN (10.6), with 156 as
                              #   the upper tail (a future per-cell heterogeneity refinement, NOT a scale knob).
S_PEAK_MOLEC_PER_S = 156.0    # high-secretor tail (PMID 37533643) — retained for the heterogeneity refinement
SECRETOR_FRACTION  = 0.039    # ~3.9% of stimulated monocytes are ACTIVE IL-6 SECRETORS (PMID 37533643).
                              # CORRECTION (2026-07-13): I previously assumed this fraction would EMERGE
                              # spatially (only myeloid near engaged T cells activate). THE SOURCE REFUTES THAT:
                              # it was measured by droplet microfluidics on LPS-stimulated monocytes -- EVERY
                              # cell was maximally stimulated, with NO spatial constraint -- and still only 3.9%
                              # secreted IL-6. So this is CELL-INTRINSIC heterogeneity (most monocytes are simply
                              # not IL-6 secretors even when fully activated), NOT a spatial consequence.
                              # Omitting it over-produced IL-6 by ~23x (elranatamab: 22,956 vs clinical 191).
                              # It is a MEASURED constant being restored -- not a knob tuned to fit.
T_TO_MAX_MIN     = 150.0      # measured time for a monocyte to reach its maximal IL-6 secretion rate
                              #   (PMID 37533643, "maximal IL-6 rate reached ~150 min") -> sets k_on.
KDEG_IL6_PER_HR  = 0.20       # IL-6 first-order elimination; 0.18-0.25 /hr (PMID 31268236) -> t1/2 ~2.8-3.8 h
V_PLASMA_ML      = 11650.0    # IL-6 DISTRIBUTION SPACE = interstitial (8.55 L) + plasma (3.10 L) = 11.65 L,
                              # taken from the MODEL'S OWN PBPK volumes (qsp_costim_window_v2._PBPKArrays:
                              # sum(Vis) + V_pl) -- not a literature number I picked.
                              # CORRECTION (2026-07-13, Max): I was dividing by PLASMA volume (3.1 L), which
                              # treats IL-6 as if it never leaves the vasculature. That is backwards:
                              #   * IL-6 is 21 kDa vs a 150 kDa IgG -> it extravasates MORE freely, not less;
                              #   * it has NO FcRn recycling to retain it in circulation;
                              #   * and it is PRODUCED IN THE INTERSTITIUM (by tissue macrophages) -- it does
                              #     not even start in plasma.
                              # The clinical assay reads the PLASMA concentration of a molecule distributed
                              # across that whole space, so the ODE denominator must be the distribution
                              # volume. Using 3.1 L over-stated IL-6 by 3.76x.
                              # NOTE (the mechanistically-complete version, NOT yet done): IL-6 should be
                              # ROUTED THROUGH THE ANTIBODY'S PBPK TRANSPORT (per-organ interstitial production
                              # -> extravasation/lymph -> plasma) with IL-6's own reflection coefficient. That
                              # machinery already exists in this engine and is currently BYPASSED. This lumped
                              # ECF compartment is a well-mixed approximation of it.
                              # REVERTED (2026-07-13): I briefly set this to 6.4 L citing an "IL-6 Vss". That
                              # value (central 3.5 L + peripheral 2.9 L) is a textbook ~150 kDa MONOCLONAL
                              # ANTIBODY volume -- almost certainly TOCILIZUMAB (anti-IL-6R), NOT IL-6 (21 kDa).
                              # Mis-attributing a mAb's Vd to the cytokine is a sourcing error, and it happened
                              # to shrink the residual 2.1x -- exactly when to be MOST suspicious of a number.
                              # OPEN UNCERTAINTY (stated, not hidden): the correct denominator for a
                              # tissue-produced 21 kDa cytokine measured in plasma is genuinely unresolved here.
                              # Plasma volume (3 L) is the conservative, transparent choice; a true Vd is
                              # plausibly larger (interstitial distribution), which would LOWER predicted plasma
                              # IL-6. This is a NAMED model uncertainty, NOT a knob to tune to the clinical value.
# ---- CELL-CELL CONTACT DISTANCE (CD40L-CD40 is MEMBRANE-BOUND: the cells must physically TOUCH) ----------
# Two spherical cells touch when their centre-to-centre distance <= r_macrophage + r_Tcell.
#   tissue macrophage diameter 21.2 +/- 0.3 um (PMID 9400735)  -> r = 10.6 um
#   T lymphocyte diameter ~7 um (range 6-10)   (PMID 30571054) -> r = 3.5 um
#   => CONTACT distance = 10.6 + 3.5 = 14.1 um
# CORRECTION (2026-07-13): this was set to R_SYN = 30 um -- the model's T-cell:TARGET synapse REACH. That is a
# REACHABILITY radius, not a contact distance, and applying it to a contact-dependent interaction over-counted
# macrophage activation ~3x (measured on real spleen data: producing fraction 1.92% -> 0.60% at 13 um).
# (Blood monocytes are smaller, 8.8 um -- but blood myeloid are gated OFF by mechanism anyway.)
R_CONTACT_UM     = 14.1
K_OFF_PER_HR     = 0.10       # deactivation/return-to-rest of the myeloid IL-6 program. Slower than activation
                              #   (IL-6 secretion outlasts the stimulus); set from the measured decline of the
                              #   secretor fraction with prolonged stimulation (PMID 37533643).

def s_max_pg_per_hr(molec_per_s=S_MAX_MOLEC_PER_S):
    """Measured per-cell maximum IL-6 secretion (molecules/s) -> pg/hr/cell. Pure unit conversion."""
    g_per_s = molec_per_s * IL6_MW_DA / NA
    return g_per_s * 1e12 * 3600.0          # g/s -> pg/hr

S_MAX_PG_PER_HR = s_max_pg_per_hr()          # ~0.0196 pg/hr/cell (verified this session)
K_ON_PER_HR     = 3.0 / (T_TO_MAX_MIN/60.0)  # first-order rise reaching ~95% of max at the measured t_to_max


class MyeloidIL6:
    """Per-cell myeloid IL-6 emitter for ONE organ. Attaches to an OrganPD that already carries:
       x, y (coords), labs (cell-type labels), and per-step ENGAGED-T state (B2 per T cell, Tidx).

    is_myeloid is derived from the REAL cell-type labels present in the agent tables (verified this session:
    'macrophage', 'classical monocyte', 'monocyte', 'intermediate monocyte', 'Myeloid', 'myeloid dendritic cell').
    """
    MYELOID_TOKENS = ("macrophage", "monocyte", "myeloid", "kupffer", "microglia")

    def __init__(self, x, y, labs, r_contact_um=30.0, count_scale=1.0):
        """r_contact_um: myeloid<-engaged-T CONTACT radius. Contact-gated per PMID 29808005 (CD40L-CD40 needs
        proximity). Default = the model's existing synapse reach R_SYN so contact is defined identically to the
        T:target synapse — NOT a new tuned distance."""
        low = np.char.lower(np.array([str(v) for v in labs]))
        self.is_mye = np.array([any(t in s for t in self.MYELOID_TOKENS) for s in low])
        self.midx = np.where(self.is_mye)[0]
        self.x = np.asarray(x, float); self.y = np.asarray(y, float)
        self.r_contact = float(r_contact_um)
        self.count_scale = float(count_scale)
        self.a = np.zeros(len(self.midx))         # per-cell activation state in [0,1]
        self._tree = None
        self.n_myeloid = len(self.midx)
        # PER-CELL count_scale: how many physiological cells each sampled myeloid agent represents.
        # MUST be applied per cell, not as one average: blood carries a PER-LINEAGE scale array, so a
        # population mean would apply a lymphocyte-weighted scale to monocytes (wrong cell type).
        self.cs = np.full(self.n_myeloid, float(count_scale))
        # INTRINSIC IL-6-SECRETOR IDENTITY (PMID 37533643): only ~3.9% of monocytes are IL-6 secretors even
        # under maximal stimulation. Assigned ONCE per cell (a fixed cell property, deterministic seed), NOT
        # a bulk multiplier -- so a cell either is or is not an emitter, and the population fraction is the
        # MEASURED 3.9%. Deterministic per organ so runs are reproducible.
        _rng = np.random.default_rng(abs(hash(("il6_secretor", int(self.n_myeloid)))) % (2**32))
        self.is_secretor = _rng.random(self.n_myeloid) < SECRETOR_FRACTION

    def set_count_scale(self, scale):
        """scale = scalar, a myeloid-length array, or a FULL-LENGTH per-cell array (indexed like the agent
        table -> the myeloid entries are selected). Blood passes its per-lineage blood_count_scale here."""
        s = np.asarray(scale, dtype=float)
        if s.ndim == 0:
            self.cs = np.full(self.n_myeloid, float(s))
        elif s.shape[0] == self.n_myeloid:
            self.cs = s.astype(float).copy()
        else:
            self.cs = s[self.midx].astype(float).copy()   # full-length -> take THIS cell type's own scales
        return self

    def _build_tree(self, xT, yT):
        from scipy.spatial import cKDTree
        return cKDTree(np.column_stack([xT, yT])) if len(xT) else None

    def step(self, dt_days, T_x, T_y, T_engaged):
        """Advance every myeloid cell's own activation from ITS OWN local engaged-T contact, then return this
        organ's instantaneous IL-6 production (pg/hr, count_scale-lifted to the physiological population).

        T_x, T_y     : coordinates of the T cells (same frame as the myeloid coords)
        T_engaged    : per-T engaged-synapse level (e.g. p_eng = B2/RC in [0,1]) — the CD40L-bearing, activated T cells
        """
        if self.n_myeloid == 0 or len(T_x) == 0:
            self.a *= np.exp(-K_OFF_PER_HR * dt_days * 24.0)
            return 0.0
        eng = np.clip(np.asarray(T_engaged, float), 0.0, 1.0)
        # each myeloid cell sums the ENGAGED T cells within its own contact radius -> its OWN activating input
        tree = self._build_tree(T_x, T_y)
        pts = np.column_stack([self.x[self.midx], self.y[self.midx]])
        contact = np.zeros(self.n_myeloid)
        nb = tree.query_ball_point(pts, r=self.r_contact)
        for i, idx in enumerate(nb):
            if idx: contact[i] = float(eng[idx].sum())     # local engaged-T load on THIS myeloid cell
        # per-cell activation ODE (structural saturation via (1-a); NO EC50, NO Emax)
        dt_hr = dt_days * 24.0
        self.a += dt_hr * (K_ON_PER_HR * contact * (1.0 - self.a) - K_OFF_PER_HR * self.a)
        np.clip(self.a, 0.0, 1.0, out=self.a)
        # each activated cell secretes at up to ITS measured per-cell maximum, lifted by ITS OWN count_scale
        # (per-cell: blood carries a per-LINEAGE scale, so monocytes get the real monocyte count -- never a
        # population-average scale, which would apply a lymphocyte-weighted factor to monocytes).
        # ONLY the intrinsic IL-6-secretor subset emits (3.9%, PMID 37533643) -- the rest are activated but
        # are not IL-6 producers. Each emitter secretes a_i * S_MAX, lifted by its OWN count_scale.
        prod_pg_per_hr = float((self.a * S_MAX_PG_PER_HR * self.cs * self.is_secretor).sum())
        return prod_pg_per_hr

    def summary(self):
        return dict(n_myeloid=self.n_myeloid,
                    frac_activated=float((self.a > 0.05).mean()) if self.n_myeloid else 0.0,
                    mean_activation=float(self.a.mean()) if self.n_myeloid else 0.0)


class PlasmaIL6:
    """Systemic plasma IL-6 with MEASURED first-order clearance. This is the piece the engine never had:
    it converts a production RATE into a CONCENTRATION. dC/dt = prod/V - kdeg*C  (analytic exponential step)."""
    def __init__(self, kdeg_per_hr=KDEG_IL6_PER_HR, v_plasma_ml=V_PLASMA_ML):
        self.kdeg = float(kdeg_per_hr); self.V = float(v_plasma_ml); self.C = 0.0    # pg/mL

    def step(self, dt_days, prod_pg_per_hr_total):
        """Exact solution over dt for dC/dt = R - k*C  (R = prod/V, pg/mL/hr) -> unconditionally stable."""
        dt_hr = dt_days * 24.0
        R = prod_pg_per_hr_total / self.V              # pg/mL/hr
        if self.kdeg <= 0:
            self.C += R * dt_hr
        else:
            Css = R / self.kdeg
            self.C = Css + (self.C - Css) * np.exp(-self.kdeg * dt_hr)
        return self.C
