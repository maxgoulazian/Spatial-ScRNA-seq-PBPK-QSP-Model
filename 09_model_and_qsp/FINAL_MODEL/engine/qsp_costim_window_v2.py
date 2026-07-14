"""
qsp_costim_window_v2.py
=======================
Enhanced multi-cell PBPK / QSP therapeutic-window model for a CD3 x TAA T-cell
engager carrying a *costimulatory second arm* (TAA x costim-receptor).

Lane: science-qsp  (costim engager toxicity counter-screen fleet)
Supersedes qsp_costim_window.py (single-node v1). Built to the fleet's
CrossRef-verified 5-layer reference architecture (QSP_MODEL_REFERENCES.md).

WHY v2 EXISTS
-------------
v1 collapsed CD4 biology into scalar multipliers on a single CD8 activation
node. The project's core thesis is that CD3 engages BOTH lineages and the
costim arm amplifies whatever the engaged cell does, so the effector (CD8)
wiring must be separated from the CD4 sub-programs (helper/CRS vs Treg/
suppression). v2 resolves the cells into distinct compartments so each screen
axis enters its OWN mechanism, and adds a minimal-physiological PBPK layer so
distribution, TMDD, and a LIVER compartment (4-1BB myeloid hepatotox) are
explicit rather than lumped.

THE CD4-vs-CD8 PROBLEM, MADE EXPLICIT (project centerpiece)
----------------------------------------------------------
CD3 engagement is lineage-blind: it forms a trimer on any T cell whose TCR/CD3
is bridged to a TAA-bearing target. The engager therefore activates:
  - CD8 effectors        -> tumor killing            (BENEFIT, effector axis)
  - CD4 conventional     -> cytokines + CD8 help     (CRS LIABILITY + benefit)
  - CD4 regulatory (Treg)-> suppression              (SUPPRESSION LIABILITY)
A costim second arm multiplies signal-2 on whichever compartment expresses the
costim receptor. The therapeutic question is whether that multiplication lands
on the CD8 killers (widen window) or feeds the CD4 CRS/Treg programs (collapse
window). v2 encodes this as compartment-specific engagement weighted by the
receptor's expression differential across CD8 / CD4-conv / Treg.

5-LAYER ARCHITECTURE (QSP_MODEL_REFERENCES.md)
----------------------------------------------
  L1 PK/distribution : central + peripheral + tumor + liver + lymphoid,
                       IgG biodistribution + TMDD sinks.        [Schropp 2019; Hadjigeorgiou 2025]
  L2 Trimer/synapse  : 3-body binding; per-compartment ternary complex from
                       affinity x expression differential.      [Douglass 2013; Betts 2019]
  L3 Activation node : resolved CD8 / CD4-conv / Treg / myeloid activation;
                       killing + resolved cytokines.            [Ma 2020; Yang 2023]
  L3b Cytokine casc. : TNF/IFNg/IL2 (CD4-conv) -> IL-6 (myeloid). [Hosseini 2020; Teachey 2016]
  L4 Suppression fb  : Treg expansion discounts net CD8 killing. [Serrano 2024; Kasai 2026]
  L5 Window metric   : dose-to-kill / dose-to-CRS-severity;
                       + explicit liver-tox readout.            [Dudal 2016; Elmeliegy 2024]

CALIBRATION PHILOSOPHY (honest scope)
-------------------------------------
Screen z-scores are RELATIVE (empirical-null normalized), so cross-arm ranking
stays COMPARATIVE: anchored so (a) the CD3-only backbone is sub-curative and
(b) CD28 reproduces the pan-costim / TGN1412 window collapse. The PBPK/PD layer
is additionally shape-validated against literature PK/PD (Betts/Schropp/Hosseini)
so absolute time-courses are credible, but per-arm WINDOW numbers remain
comparative, not clinical predictions.
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass, field
from scipy.integrate import solve_ivp

# ============================================================================
# LAYER 1 — FULL-BODY PBPK (Shah & Betts 2012 platform; 2-pore + FcRn)
# ============================================================================
# Whole-body mAb PBPK: 15 tissues + tumor, each with a vascular (plasma) and an
# interstitial sub-compartment, connected by blood flow Q and lymph flow L
# (L = Q/500, Shah-Betts). Antibody extravasates by 2-pore CONVECTION governed
# by the vascular reflection coefficient sigma_V (tight vs leaky endothelium);
# it returns to blood via lymph (reflection sigma_L). FcRn salvage is modelled
# as an endosomal recycled fraction (quasi-steady-state reduction of catabolism)
# rather than an explicit on/off endosome ODE, for numerical robustness while
# preserving the platform's long IgG half-life. Portal organs drain to liver;
# liver + all organs drain to a venous pool -> lung (series) -> arterial pool.
#
# Physiological volumes/flows are reference-human (ICRP/Brown 71 kg); the
# 4 platform system parameters (FcRn salvage, pinocytosis, Kdeg, lymph return)
# are set to platform-typical values giving IgG t1/2 ~ 2-3 wk and are refined by
# the PK/PD validation track against the literature TCE-QSP PK curves.
#
# Tissue arrays are frozen at import from the tabulated physiology below.
# ---------------------------------------------------------------------------
_PBPK_TISSUES = {
    # name          V(L)   Qfrac  sigV   fV     fIS    portal
    "lung":       (0.50,  1.000, 0.95, 0.105, 0.188, False),
    "heart":      (0.33,  0.040, 0.95, 0.157, 0.320, False),
    "kidney":     (0.31,  0.190, 0.90, 0.105, 0.200, False),
    "brain":      (1.45,  0.120, 0.99, 0.037, 0.150, False),
    "muscle":     (30.0,  0.170, 0.95, 0.026, 0.120, False),
    "skin":       (3.40,  0.050, 0.95, 0.038, 0.302, False),
    "adipose":    (13.0,  0.050, 0.95, 0.020, 0.135, False),
    "bone":       (10.0,  0.050, 0.95, 0.041, 0.100, False),
    "stomach":    (0.15,  0.010, 0.90, 0.038, 0.170, True),
    "small_int":  (0.65,  0.100, 0.90, 0.038, 0.200, True),
    "large_int":  (0.37,  0.040, 0.90, 0.038, 0.200, True),
    "pancreas":   (0.18,  0.010, 0.90, 0.038, 0.180, True),
    "spleen":     (0.19,  0.030, 0.80, 0.110, 0.200, True),
    "liver":      (1.80,  0.065, 0.85, 0.115, 0.200, False),
    "tumor":      (0.10,  0.020, 0.75, 0.070, 0.300, False),
}
_PLASMA_CO = 5000.0        # L/day total plasma flow (cardiac output x plasma fraction)
_LYMPH_RATIO = 1.0/500.0   # Shah-Betts: lymph = plasma flow / 500

# Relative T-cell and myeloid density per tissue (spleen=1). Drives the CD3-arm
# receptor capacity (=> where TMDD is strongest) and the myeloid-mediated liver-
# tox term. DATA-GROUNDED: fractions of T cells and myeloid cells per tissue were
# computed from the Tabula Sapiens single-cell atlas (~1.14M primary cells across
# 24 tissues; CZ CELLxGENE Census stable 2025-11-08), then normalised to spleen.
# Brain (immune-privileged; absent in TS) and tumor (TIL prior from project brief)
# keep literature priors. See qsp_tabula_sapiens_densities.csv for provenance.
# Reference: Tabula Sapiens Consortium, Science 2022, DOI 10.1126/science.abl4896.
_TCELL_DENSITY = {
    "lung":0.309, "heart":0.117, "kidney":0.292, "brain":0.020, "muscle":0.183,
    "skin":1.510, "adipose":0.767, "bone":0.860, "stomach":1.028, "small_int":1.240,
    "large_int":2.359, "pancreas":0.182, "spleen":1.000, "liver":0.429, "tumor":0.400,
}
# Myeloid density per tissue (spleen=1), same Tabula Sapiens source. Liver is
# myeloid-rich (Kupffer cells, 0.84) — this is the substrate for 4-1BB / costim
# hepatotox, which the imposed-prior version understated.
_MYELOID_DENSITY = {
    "lung":0.924, "heart":0.306, "kidney":0.106, "brain":0.050, "muscle":0.178,
    "skin":0.572, "adipose":0.564, "bone":1.287, "stomach":0.305, "small_int":0.035,
    "large_int":0.017, "pancreas":0.394, "spleen":1.000, "liver":0.839, "tumor":0.300,
}

class _PBPKArrays:
    """Frozen physiological arrays for the full-body PBPK (built once at import)."""
    def __init__(self):
        n = list(_PBPK_TISSUES)
        self.names = n
        self.N = len(n)
        self.V   = np.array([_PBPK_TISSUES[k][0] for k in n])
        Qf       = np.array([_PBPK_TISSUES[k][1] for k in n])
        self.sigV= np.array([_PBPK_TISSUES[k][2] for k in n])
        fV       = np.array([_PBPK_TISSUES[k][3] for k in n])
        fIS      = np.array([_PBPK_TISSUES[k][4] for k in n])
        self.portal = np.array([_PBPK_TISSUES[k][5] for k in n])
        # normalize parallel-organ flow to conserve the circulatory loop
        mask = np.array([k != "lung" for k in n])
        Qf = Qf.copy(); Qf[mask] = Qf[mask]/Qf[mask].sum()
        self.Q = Qf*_PLASMA_CO; self.Q[n.index("lung")] = _PLASMA_CO
        self.L = self.Q*_LYMPH_RATIO
        self.Vv  = fV*self.V          # vascular sub-volume per organ
        self.Vis = fIS*self.V         # interstitial sub-volume per organ
        self.i_lung  = n.index("lung")
        self.i_liver = n.index("liver")
        self.i_tumor = n.index("tumor")
        self.i_portal = np.where(self.portal)[0]
        self.i_parallel = np.where(mask & (np.arange(self.N)!=self.i_lung))[0]
        # T-cell density -> CD3-arm receptor capacity per tissue (Tabula Sapiens)
        self.tcell = np.array([_TCELL_DENSITY[k] for k in n])
        # myeloid density -> local CRS amplification + liver-tox substrate (Tabula Sapiens)
        self.myeloid = np.array([_MYELOID_DENSITY[k] for k in n])
        # tumor-associated antigen present only in tumor (TAA-arm capacity)
        self.taa = np.zeros(self.N); self.taa[self.i_tumor] = 1.0
        # --- SPATIAL hooks (default = pure well-mixed) ---------------------------
        # spatial_exposure[organ] in (0,1]: scales drug conc reaching each organ's
        # interstitial node (reaction-diffusion penetration penalty from the
        # all-organ Xenium spatial layer). All 15 organs default to 1.0, so every
        # organ is pure well-mixed until the spatial lane delivers per-organ values
        # (spatial_exposure_perorgan.csv). Interface contract: 1.0 => unchanged.
        self.spatial_exposure = np.ones(self.N)

PB = _PBPKArrays()   # module-level frozen physiology


@dataclass
class PBPK:
    """Tunable full-body PBPK system parameters (physiology is frozen in PB).
    The 4 platform params (Shah-Betts) plus tumor-target TMDD."""
    sigL: float = 0.85          # lymph (return) reflection coefficient. CALIBRATED so the
                                # steady-state interstitial:plasma amount ratio gives
                                # Vss/Vc ~= 2.1 (class-typical mAb; pembrolizumab 2.17,
                                # trastuzumab 2.7, mosunetuzumab 2.1). SS ratio is
                                # (1-sigV)/(1-sigL), independent of k_dist, so this knob
                                # sets DISTRIBUTION EXTENT (Vss) without touching the rate.
    k_dist: float = 3.0         # 2-pore convective distribution-RATE multiplier on both
                                # extravasation and lymph return. Cancels from the SS ratio
                                # (rate-only knob): sets the alpha (distribution) phase
                                # DURATION/depth to match the class-typical ~1.3 d biexp
                                # distribution, calibrated vs pembrolizumab day-1/3/7 fall.
    fFcRn: float = 0.90         # FcRn salvage fraction of pinocytosed mAb (recycled)
    CLup: float = 0.3503        # /day plasma pinocytic uptake -> (1-fFcRn) catabolised.
                                # CALIBRATED (PK/PD validation track) so the backbone
                                # IgG-engager terminal t1/2 matches the mosunetuzumab
                                # clinical anchor 16.1 d exactly (FDA LUNSUMIO label);
                                # comparative window ranking unchanged (tumor/plasma
                                # AUC ratio 0.234 and prozone bell-shape preserved).
                                # CLup*(1-fFcRn) = 0.0385/day => IgG terminal t1/2 ~ 18-20 d
    Kdeg: float = 26.0          # /day catabolism of UNsalvaged endosomal mAb (fast)
    k_lymph_return: float = 24.0  # /day lymph pool -> plasma return (fast turnover)
    # ---- MOLECULAR PROPERTIES (drive real-unit conversion + format-dependent PK) ----
    # A single mechanistic PK layer must reproduce the whole T-cell-engager class,
    # whose plasma half-lives span TWO ORDERS OF MAGNITUDE, entirely from format:
    #   * canonical IgG mAb / IgG-based TCE  -> intact Fc -> FcRn salvage -> t1/2 ~ weeks
    #   * BiTE (tandem scFv, NO Fc)          -> NO FcRn recycling + small (<glomerular
    #                                           cutoff) so renally filtered -> t1/2 ~ hours
    # So Fc status and molecular weight are the two levers, NOT free per-molecule knobs.
    mw_kda: float = 146.9       # molar mass (kDa); IgG ~146.9. Sets ug/mL <-> nM and renal sieving.
    has_fc: bool = True         # intact Fc? If False, FcRn salvage is OFF (fFcRn forced 0).
    # Size-gated renal/catabolic clearance: a sigmoid in MW around the ~69 kDa
    # glomerular filtration cutoff (serum albumin). ~0 for a 150 kDa IgG (not
    # filtered); dominant for a ~55 kDa BiTE (filtered + proteolysed). This is the
    # SECOND half of the no-Fc story: losing FcRn alone drops t1/2 ~18 d -> ~2 d;
    # renal filtration of the small non-Fc scaffold is what takes it to ~2 h.
    k_renal_max: float = 8.70   # /day maximal size-based clearance (small-protein limit).
                                # CALIBRATED so a no-Fc BiTE at 54 kDa (blinatumomab) reaches
                                # the clinical ~2.1 h terminal t1/2 (Blincyto FDA label); the
                                # Hill sieving on MW then gives ~0/day for a 150 kDa IgG.
    mw50_renal: float = 69.0    # kDa: glomerular sieving midpoint (albumin ~66-69 kDa)
    hill_renal: float = 10.0    # steepness of the MW cutoff (sharp glomerular threshold)

    def eff_fFcRn(self) -> float:
        """FcRn salvage fraction actually in force: 0 when the construct has no Fc."""
        return self.fFcRn if self.has_fc else 0.0

    def k_renal(self) -> float:
        """Size-based plasma clearance (/day): ~0 for IgG, large for a small no-Fc BiTE."""
        return self.k_renal_max / (1.0 + (self.mw_kda / self.mw50_renal) ** self.hill_renal)
    # ---- EMERGENT TMDD via bivalent binding (Rhoden 2016) ----
    # TMDD is NOT imposed; it arises per tissue from the engager's two arms binding
    # their targets where those targets are expressed, then internalizing.
    #   CD3 arm: binds T cells (capacity ~ tissue T-cell density, PB.tcell)
    #   TAA arm: binds tumor cells (capacity ~ PB.taa, tumor only)
    KD_CD3: float = 3.0         # monovalent CD3-arm dissociation const (conc units)
    KD_TAA: float = 0.3         # monovalent TAA-arm dissociation const (higher-affinity TAA arm)
    Rcap_CD3: float = 2.0       # CD3 receptor capacity scale (x tissue T-cell density)
    Rcap_TAA: float = 6.0       # TAA receptor capacity scale in tumor (high antigen density)
    avidity: float = 8.0        # DEPRECATED heuristic cap — retained only for back-compat / legacy
                                # callers. The synapse now uses the EXACT geometric Rhoden avidity
                                # (Ag_eff from receptor surface density + arm reach); see bivalent_binding.
    # ---- GEOMETRIC AVIDITY (exact Rhoden 2016; preprint recipe 10.1101/2022.09.12.507653) ----
    # Ag_eff = local 2nd-arm antigen conc = sigma * 3/(2*r_Ab) * 1e15, sigma = surface density.
    # From a LOCAL density R_local (nM) and target-cell count n_cell in that voxel volume V, the
    # per-cell receptor number is rec_pc = R_local*V*NA/1e9/n_cell; sigma = rec_pc/(4*pi*r_cell^2).
    # Equivalent closed form used here: Ag_eff = R_local * (V/n_cell) * 3/(2*r_Ab*4*pi*r_cell^2) * 1e15/1e9.
    r_cell_um_CD3: float = 4.0   # T-cell radius (µm)
    r_cell_um_TAA: float = 8.0   # tumor-cell radius (µm)
    r_Ab_nm_avidity: float = 12.5 # antibody arm-to-arm reach (nm); AF3-derived per format
    rec_per_cell_CD3: float = 1.0e4  # CD3 copies per T cell (for Ag_eff geometry)
    rec_per_cell_TAA: float = 1.0e5  # TAA copies per tumor cell (for Ag_eff geometry)
    # ---- CONSTRUCT VALENCY (format geometry) ----
    # Number of binding arms against each cell target. Avidity (multivalent
    # affinity enhancement) applies to a cell ONLY when the construct engages it
    # with >=2 arms. Examples: 1+1 bispecific (n_CD3=1,n_TAA=1) -> no avidity on
    # either; 2+1 (n_CD3=2,n_TAA=1) -> avidity on the T cell (CD3) only; 1+2 ->
    # avidity on tumor (TAA) only; 2+2 -> avidity on both. This makes the binding
    # layer format-aware so affinity/valency can be optimised downstream.
    n_CD3: int = 1              # arms against CD3 (T cell)
    n_TAA: int = 1              # arms against TAA (tumor cell)
    n_costim: int = 0           # arms against the costim receptor (0 = pure CD3xTAA
                                # backbone). Reserved for the construct-format library
                                # (BiTE/IgG-scFv/CrossMab 2:1/DART-Fc/...): drives costim-
                                # arm avidity once per-format valency CSVs land. Default 0
                                # leaves the current backbone binding unchanged.
    hill_CD3: float = 1.0       # avidity Hill exponent per extra CD3 arm (set at build)
    hill_TAA: float = 1.0       # avidity Hill exponent per extra TAA arm
    # --- inter-arm REACH gate (construct-format hook) -----------------------------
    # A format's CD3<->TAA arm span must bridge the ~15 nm immune-synapse cleft to
    # form a productive trimer. reach_gate in (0,1] DOWN-WEIGHTS trimer formation
    # when the predicted span (from AF3 folds, supplied per-format) can't bridge.
    # Default 1.0 = fully reachable => trimer term unchanged (pure well-mixed).
    reach_gate: float = 1.0
    # ---- FORMAT-AWARE MECHANISTIC PD DRIVER (Schropp 2019 ternary_equilibrium) --------
    # When pd_driver_mode is set, the PD engagement (eng_tum/eng_sys) is driven by the
    # REAL complex species for the molecule's binding format, not the phenomenological
    # trimer_bell(C_opt). This removes the per-molecule fitted C_opt and the single
    # IL6_ANCHOR: cross-molecule PD magnitude then emerges from KD + antigen density.
    #   "trimer"    -> bridged ternary CD3xTAA (engagers); prozone emergent
    #   "occupancy" -> receptor occupancy C/(C+KD) (checkpoint / costim agonist)
    #   "soluble"   -> free-ligand fractional neutralization (soluble-target mAbs)
    #   "binary"    -> binary membrane complex fraction (depleting/blocking mAbs)
    #   None        -> LEGACY phenomenological bell (abstract window-scoring / nomination,
    #                  which is dimensionless and must stay untouched — DEFAULT).
    pd_driver_mode: str = None
    # ONE intrinsic normalization per model (NOT per molecule): the reference bridged-trimer
    # concentration (nM) that maps to unit engagement. Calibrated ONCE to the mosunetuzumab
    # operating point so the validated downstream cascade (which expects eng ~O(1)) is
    # preserved; every other molecule's engagement is trimer/TRIMER_REF, so its magnitude is
    # set purely by its own KD/density relative to the anchor. Replaces C_opt AND IL6_ANCHOR.
    trimer_ref_nM: float = 0.0455       # mosunetuzumab peak ternary (R_TAA=741, KD_TAA=5, KD_CD3=40)
    KD_CD3_pd: float = 40.0             # CD3-arm KD for the PD trimer (weak, T-cell-sparing anti-CD3 class)
    R_CD3_pd_nM: float = 0.05           # accessible CD3 pool driving the PD synapse (plasma T-cell CD3)
    occ_KD_nM: float = 1.0              # KD for occupancy-mode drivers (checkpoint/costim); set per molecule
    liquid_target: bool = False         # True for blood/marrow targets (B-cell/BCMA engagers): killing uses the
                                        # circulating CD8_b effector pool (no solid-tumour infiltration barrier),
                                        # not the tumour-infiltrating CD8_t. A real target-compartment property.
    kint_cplx: float = 0.90     # /day internalization of the target-BOUND (bridged) complex
    kint_mono: float = 0.25     # /day internalization of singly-bound drug
    # convenience mirrors used by the QSP layer (tumor/liver interstitial volumes)
    @property
    def V_tum(self): return PB.Vis[PB.i_tumor]
    @property
    def V_liv(self): return PB.Vis[PB.i_liver]
    @property
    def V_cen(self): return PB.Vv[PB.i_lung]  # plasma proxy (arterial pool volume)


# ============================================================================
# LAYER 3 (states) — full-body PBPK PK states + resolved cell populations
# ============================================================================
# State vector layout (QSS-vascular full-body PBPK):
#   The organ vascular sub-compartments equilibrate ~10^5x faster than the
#   antibody timescale (Q_i/Vv_i up to ~1e5/day vs pharmacology over days), so
#   they are solved ALGEBRAICALLY (quasi-steady-state) each RHS call rather than
#   integrated. This is physically exact for perfusion-limited vascular spaces
#   and removes the artificial stiffness. Dynamic PK states are therefore the
#   plasma pool, the per-organ INTERSTITIAL spaces (where targets/cells live),
#   and a lymph pool that returns extravasated mAb to plasma.
#   PK block:
#     0              A_pl      plasma pool amount
#     [1 : N+1)      A_is_i    interstitial amount per organ i
#     N+1            A_ly      lymph pool amount
#   QSP block (offset QOFF = N+2) — the 3 MOBILE T-cell types each have a BLOOD
#   (_b) pool and a TUMOR (_t) pool that TRAFFIC (blood<->tumor); myeloid is
#   resident. Killing comes ONLY from CD8_t (tumor); systemic CRS from CD4c_b:
#     QOFF+0  CD8_b   circulating CD8 (traffics to tumor; not a killer in blood)
#     QOFF+1  CD8_t   TUMOR CD8 effector (the KILLER; efficacy site)
#     QOFF+2  CD4c_b  circulating CD4-conv (systemic activation = CRS engine)
#     QOFF+3  CD4c_t  TUMOR CD4-conv (local help -> CD8 persistence)
#     QOFF+4  Treg_b  circulating Treg (systemic expansion liability)
#     QOFF+5  Treg_t  TUMOR Treg (local suppression of CD8 killing)
#     QOFF+6  Mye     activated myeloid/monocyte (RESIDENT; CRS amplifier + liver tox)
#     QOFF+7  Tumor   tumor burden (logistic)
#     QOFF+8  cTNF    systemic TNF   (T-cell storm-core)
#     QOFF+9  cIFN    systemic IFN-g (T-cell storm-core; drives myeloid IL-6)
#     QOFF+10 cIL2    systemic IL-2  (T-cell storm-core)
#     QOFF+11 cIL6    systemic IL-6  (MYELOID-derived; CRS-severity driver)
#     QOFF+12 Liv     liver-injury readout (integrated myeloid-in-liver activity)
#     QOFF+13 Dcyto   acute cytokine-desensitization state (self-limits the first-dose spike)
#     QOFF+14 Ecum    slow FRESH-DELIVERY / serial-engagement state in [0,1). Integrates
#                     the DRUG DELIVERY RATE (ongoing infusion), NOT standing concentration.
#                     Fresh drug arriving forms NEW CD3 synapses on cells recovered from
#                     desensitization; a continuous infusion (blina) charges Ecum to a high
#                     plateau -> sustained cytokine tail, while a single bolus (mosun)
#                     delivers once (input rate ->0 after t=0) so Ecum stays ~0 and the tail
#                     resolves. This is the mechanistic driver of the bolus-vs-infusion tail
#                     difference — instantaneous plasma conc is nearly identical (and even
#                     crosses over) for the two, so it cannot supply the ~4x clinical gap.
#                     Every window-SCORED arm is a bolus (no infusion) -> Ecum=0 -> the tail
#                     reduces to the arm-dependent CRS biology (gC), so the NOMINATION is
#                     structurally independent of this delivery term.
_N = PB.N
I_PL = 0
I_IS0 = 1                     # first interstitial index; organ i -> I_IS0+i
I_LY = _N + 1
QOFF = _N + 2
NSTATE = QOFF + 15            # 6 trafficking T-cell pools + Mye + Tumor + 4 cyto + Liv + Dcyto + Ecum
# QSP state layout — the three mobile T-cell types each have a BLOOD (_b) and a
# TUMOR (_t) pool (trafficking); myeloid is resident; then tumor, 4 cytokines, liver,
# and a systemic cytokine-DESENSITIZATION state (Dcyto) that reproduces the acute
# first-dose CRS spike-and-resolve (Hosseini 2020; step-up-dosing rationale): cytokine
# release is driven by ACUTE engagement gated by (1-Dcyto), and Dcyto builds under
# engagement (fast) and recovers slowly, so the burst self-limits within ~1-2 days
# even while drug (and slow T-cell expansion) persist for weeks.
IDX = dict(CD8_b=QOFF+0, CD8_t=QOFF+1, CD4c_b=QOFF+2, CD4c_t=QOFF+3,
           Treg_b=QOFF+4, Treg_t=QOFF+5, Mye=QOFF+6, Tumor=QOFF+7,
           cTNF=QOFF+8, cIFN=QOFF+9, cIL2=QOFF+10, cIL6=QOFF+11, Liv=QOFF+12,
           Dcyto=QOFF+13, Ecum=QOFF+14)
def _cis_amt(y): return y[I_IS0:I_IS0+_N]        # interstitial amounts per organ

# plasma & lymph physiological volumes (L)
V_PLASMA = 3.1
V_LYMPH  = 2.6


@dataclass
class CellParams:
    """Kinetics for the resolved cell compartments (relative units)."""
    # baseline pools (relative) — T-cell types now TRAFFIC between a BLOOD pool
    # and a TUMOR pool (see trafficking params). Blood baselines are the
    # circulating pool; tumor baselines are the resident TIL fraction at t=0.
    CD8_b0: float = 1.0       # CD8 blood baseline
    CD8_t0: float = 0.20      # CD8 tumor-resident (TIL) baseline
    CD4c_b0: float = 1.0      # CD4-conv blood baseline
    CD4c_t0: float = 0.15     # CD4-conv tumor baseline (fewer effectors infiltrate)
    Treg_b0: float = 0.10     # Treg blood baseline
    Treg_t0: float = 0.05     # Treg tumor baseline (Tregs are tumor-enriched rel. to their blood pool)
    Mye_0: float = 0.30       # myeloid stays RESIDENT (tumor/liver), not trafficked as a pool
    # activation rates (per unit compartment-engagement, /day)
    kact_CD8: float = 1.5
    kact_CD4c: float = 1.3
    kact_Treg: float = 0.9
    kact_Mye: float = 1.0     # myeloid activated by T-cell cytokines, not trimer directly
    # contraction / decay (/day)
    kdec_CD8: float = 0.16
    kdec_CD4c: float = 0.20
    kdec_Treg: float = 0.15
    kdec_Mye: float = 0.8     # myeloid resolves faster
    # ---- TRAFFICKING (blood <-> tumor), /day ----
    # T cells infiltrate the tumor from blood and egress back. Infiltration is
    # ENHANCED by local engagement/inflammation (activated T cells upregulate
    # adhesion + the inflamed tumor recruits more). Egress rate is derived per
    # type in the RHS so the untreated baseline is an exact fixed point
    # (k_infil*blood_0 == k_egress*tumor_0). CD8 is given the strongest
    # activation-driven infiltration (the efficacy lever); Treg the weakest.
    k_infil: float = 0.15     # base blood->tumor infiltration rate
    infil_gain_CD8: float = 4.0   # activation-driven CD8 infiltration boost (per unit tumor engagement)
    infil_gain_CD4c: float = 2.0
    infil_gain_Treg: float = 1.5
    # ---- SYSTEMIC engagement (CRS driver) ----
    # CD3 is lineage-blind and present on ALL T cells, so the drug engages the
    # circulating (blood) T-cell pool everywhere — but WITHOUT TAA there is no
    # full cytotoxic trimer, so systemic engagement is a fraction of the
    # tumor-localised trimer. This weak systemic CD4 engagement is the CRS origin
    # (why even TAA-targeted TCEs cause cytokine release; TGN1412 was worse
    # because a superagonist crosslinks CD3/CD28 with no TAA needed).
    sys_eng_frac: float = 0.35    # systemic (no-TAA) engagement efficiency vs tumor trimer
    # killing (calibrated so the CD3-only backbone is SUB-CURATIVE ~37% tumor
    # reduction at a co-tolerable dose -> costim's effector boost has real room)
    kkill: float = 35.0       # per TUMOR-CD8 effector, /day. Re-anchored for the
                              # trafficking + TRIMER-LICENSED-killing structure
                              # (killing comes only from tumor CD8, and only when
                              # the drug licenses the synapse) against a GROWING
                              # tumor (untreated ~6x/56d); backbone SUB-CURATIVE
                              # (~27% best control, holds but cannot clear ->
                              # costim's effector boost has real room to matter).
    # tumor
    tumor0: float = 100.0     # STARTING tumor burden at dosing
    tumor_K: float = 800.0    # logistic CARRYING CAPACITY (>> tumor0 so an
                              # untreated tumor has room to grow ~exponentially
                              # early, then decelerates near K)
    kgrow: float = 0.055      # logistic growth rate (verified: untreated tumor
                              # ~5.7x over the 56-day window; doubling ~11-12 d)
    # suppression: Treg blunts CD8 killing (gentle Hill, slope beta)
    supp_beta: float = 0.35
    # help: CD4-conv sustains CD8 (persistence) — reduces CD8 decay in proportion to CD4c help
    help_gain: float = 0.35   # max fractional CD8-decay slowdown from full CD4c help
    help_kill_gain: float = 0.30  # per-CD8 kill-EFFICIENCY penalty per unit help-erosion
                                  # (gH-1); how strongly eroded CD4 help degrades CD8
                                  # cytotoxic quality (the CD30 help-erosion efficacy cost)


@dataclass
class CytokineParams:
    """Resolved cytokine cascade: T-cell storm-core -> myeloid IL-6 (Hosseini 2020)."""
    # T-cell storm-core production per unit activated T-cell (CD4-conv dominant, CD8 minor)
    p_TNF_CD4c: float = 0.10
    p_IFN_CD4c: float = 0.10
    p_IL2_CD4c: float = 0.08
    p_TNF_CD8: float = 0.035   # CD8 contributes but CD4-conv dominates the storm
    p_IFN_CD8: float = 0.045
    # partial decoupling of storm from raw effector potency (cytotoxicity-cytokine split)
    cyto_eff_coupling: float = 0.5
    # myeloid amplification: TNF+IFN activate myeloid -> IL-6 (the severity driver)
    k_mye_act: float = 1.0    # T-cell TNF/IFN -> myeloid activation
    p_IL6_Mye: float = 2.0    # IL-6 per activated myeloid (amplification)
    # clearances (/day, cytokines fast)
    d_TNF: float = 2.0
    d_IFN: float = 1.8
    d_IL2: float = 2.0
    d_IL6: float = 6.0
    # ---- acute first-dose CRS: engagement-driven burst + desensitization ----
    # Cytokine release is driven by ACUTE systemic engagement (not the slow T-cell
    # pool size), gated by a desensitization state Dcyto in [0,1). Dcyto builds fast
    # under engagement and recovers slowly, so the storm peaks within ~1 day and
    # resolves within ~2-3 days even while drug persists (Hosseini 2020; the
    # clinical basis for step-up dosing). Tuned to mosun/blina IL-6 spike-resolve.
    k_cyto_burst:  float = 120.0   # amplitude of engagement-driven storm production
    k_cyto_desens: float = 30.0    # /day desensitization build rate per unit engagement
    k_cyto_recov:  float = 0.003   # /day recovery of cytokine responsiveness
    f_IL6_burst:   float = 15.0    # burst_drive multiplier feeding IL-6 directly (acute spike)
    # per-cytokine acute-burst multipliers (set to preserve the IL-6>IFN>TNF>IL-2
    # hierarchy at the acute peak; the burst is co-timed across cytokines and
    # IL-6-dominant, matching the clinical CRS hierarchy).
    f_TNF_burst:   float = 2.31
    f_IFN_burst:   float = 2.55
    f_IL2_burst:   float = 1.38
    sust_floor:    float = 0.01   # minimal tonic sustained-cytokine fraction
    sust_hill:     float = 1.0    # Hill exponent mapping cumulative engagement Ecum -> sustained tail
    k_cum:         float = 9.0    # /day charge rate of fresh-delivery state Ecum
    k_cum_decay:   float = 0.30   # /day decay rate of Ecum (tail-resolution time constant for a bolus)
    EC50_cum:      float = 2.0    # plasma conc (µg/mL) at half-max CD3 occupancy gating fresh delivery
    inf_ref:       float = 30.0   # reference infusion rate (mg/day) normalizing the delivery flux
    mye_floor:     float = 0.03   # minimal tonic myeloid activation (rest gated by desensitization resp)
    # liver injury: myeloid activity IN the liver compartment integrates to injury
    k_liv_injury: float = 1.0
    d_liv_repair: float = 0.3
    # --- SPATIAL liver-tox hook (per costim arm) ---------------------------------
    # spatial_kupffer_coloc in (0,1]: scales the costim->Kupffer engagement feeding
    # the Liv injury state, from the spatial layer's costim-arm/Kupffer co-
    # localization (spatial_kupffer_coloc_perarm.csv). MODULATE-ONLY: it can only
    # LOWER the liver-tox drive, never rescue the liver ceiling (the tox analogue of
    # the tumor efficacy penalty). Default 1.0 => liver-tox drive unchanged.
    spatial_kupffer_coloc: float = 1.0


# ============================================================================
# LAYER 2 — per-arm screen scores + costim-receptor EXPRESSION DIFFERENTIAL
# ============================================================================
@dataclass
class ArmScores:
    """Agonism-oriented screen z-scores for one costim arm (Stim48hr primary),
    plus the costim-receptor EXPRESSION differential across compartments.

    Screen axes (functional wiring, HOW STRONG the response is per engagement):
        effector_z : Axis1 Schmidt CD8 IFNG CRISPRa z (higher = more CD8 effector)
        crs_z      : Axis2 CD4 CRS empirical-null z    (higher = more CRS)
        supp_z     : Axis3 CD4 suppression z           (higher = more Treg/IL-10)
        help_z     : Axis4 CD4 help-preservation z     (higher = more help)

    Expression differential (WHERE the costim second-arm boost is delivered):
        expr_CD8, expr_CD4c, expr_Treg  — relative costim-receptor surface
        density on each compartment. This routes the costim boost. Defaults to
        uniform (1,1,1) = neutral; populated by set_expression_prior() from the
        RNA-to-surface-protein layer (here a literature/project-stated PRIOR,
        overridable). The CD3 arm itself is pan-T (uniform) — only the costim
        boost is expression-gated, which is the crux of the CD4-vs-CD8 problem.
    """
    name: str
    effector_z: float
    crs_z: float | None
    supp_z: float | None
    help_z: float | None
    effector_hit: bool = True
    tier: str = "core-11"
    intrinsic_tox: str = ""
    # expression differential (relative costim-receptor density per compartment)
    expr_CD8: float = 1.0
    expr_CD4c: float = 1.0
    expr_Treg: float = 1.0
    # TYPED suppression decomposition (composition vs state), populated by the tox
    # lane's composition-vs-state analysis; None => fall back to splitting supp_z.
    supp_size_z: float | None = None   # Treg-compartment-SIZE component (compositional)
    supp_rate_z: float | None = None   # per-cell suppressive-RATE component (state)
    frac_size: float = 0.5             # split of supp_z into size vs rate when un-decomposed
    # ACTIVATION-KINETICS scores (Rest -> Stim8hr -> Stim48hr), for PD onset:
    crs_z_8hr: float | None = None     # early (8hr) CRS z; None => use sustained crs_z
    crs_z_48hr: float | None = None    # sustained (48hr) CRS z
    supp_z_8hr: float | None = None
    supp_z_48hr: float | None = None
    # SINGLE-CELL help direction (science-latent D1 continuous-distribution test).
    # The pseudobulk help_z can have the WRONG SIGN: e.g. CD30 pseudobulk help_z
    # = +0.107 (mild BONUS) but the single-cell continuous-distribution test shows
    # agonism ERODES CD4 help (lineage-shutdown). Two scales for the SAME finding:
    #   - raw effect size (Cliff's-delta, panel_help_supp_crs_D1_Stim48hr.csv,
    #     verified against the pulled file): CD30 HELP_agonism = -0.2762,
    #     mwu_p = 5.0e-5, BH q = 5.7e-4
    #   - rescaled to THIS module's robust-z frame (qsp_singlecell_zframe, the
    #     values actually loaded as sc_help_z): CD30 help_z_singlecell = -1.425
    #     (help_q = 0.001)
    # Same direction (erosion), same significance; the z-frame value is the one
    # wired below. When present, sc_help_z OVERRIDES help_z for the persistence
    # gain so the model penalises help-eroding arms on the exact axis they fail
    # (CD8-durability-from-CD4-help). Sign: >0 = help preserved/raised, <0 = eroded.
    sc_help_z: float | None = None
    # cross-donor variability (SD of the axis z across donors D1-4), for CIs:
    effector_z_sd: float = 0.0
    crs_z_sd: float = 0.0
    supp_z_sd: float = 0.0
    # ---- NEW mechanistic axes (GRN/perturbation-response driven; COMPLEMENT the veto gate) ----
    # CD8 exhaustion drive (agonism frame): >0 => agonizing this arm PUSHES CD8 toward exhaustion
    #   (accelerates CD8 contraction). Source: Schmidt CD8 CRISPRa perturbation-response where the
    #   guide exists (dense, causal); else genome-scale hero CD4-LOF cross-lineage proxy. z across panel.
    exh_z: float | None = None
    exh_z_sd: float = 0.0
    # CD8 proliferation drive (agonism frame): >0 => agonism EXPANDS the CD8 effector pool
    #   (boosts activation/expansion term). Same source hierarchy as exh_z.
    prolif_z: float | None = None
    prolif_z_sd: float = 0.0

    def gains(self, gp: "GainParams"):
        """Map screen z-scores -> multiplicative mechanistic gains, referenced
        to the CD3-only backbone (z=0 => gain 1.0 on every axis).

        TYPED SUPPRESSION (per the deepening mandate): the single suppression
        axis is decomposed into TWO mechanistically distinct gains that enter
        the window differently:
          gS_size : modifies the Treg COMPARTMENT SIZE (expansion rate). Driven
                    by the compositional part of the suppression signal
                    (supp_size_z; falls back to supp_z until tox emits the
                    composition/state decomposition).
          gS_rate : modifies PER-CELL suppressive output (per-Treg potency).
                    Driven by the state/per-cell part (supp_rate_z; same
                    fallback). An arm can negate suppression by shrinking the
                    compartment (gS_size<1), by quieting each cell (gS_rate<1),
                    or both — and these are NOT interchangeable in the window.
        The fraction split (frac_size) is set by the tox composition-vs-state
        decomposition per receptor; default 0.5/0.5 until that lands.
        """
        def _nan(x): return x is None or (isinstance(x,float) and np.isnan(x))
        gE = float(np.clip(np.exp(gp.kE * self.effector_z), 0.2, None))
        # ---- CRS gain with ACTIVATION-KINETICS onset weighting ----
        # Acute CRS grade is an ONSET phenomenon: an arm whose cytokine liability
        # is already present at 8hr floods early (dangerous acute CRS), even if it
        # settles by 48hr. We blend the sustained (48hr) and early (8hr) CRS z with
        # weight w_onset on the early rate, so fast-onset arms (CD28/CD2) are
        # penalised for acute CRS beyond their steady-state value. Falls back to the
        # single 48hr snapshot when 8hr is unavailable.
        if _nan(self.crs_z):
            gC = 1.0
        else:
            z48 = self.crs_z_48hr if not _nan(self.crs_z_48hr) else self.crs_z
            z8  = self.crs_z_8hr  if not _nan(self.crs_z_8hr)  else z48
            z_eff = (1.0 - gp.w_onset)*z48 + gp.w_onset*max(z8, z48)
            gC = float(np.exp(gp.kC * z_eff))
        # persistence gain gH from CD4-help axis. PREFER the single-cell help
        # direction (sc_help_z) when present — the pseudobulk help_z can have the
        # wrong sign (CD30: pb +0.107 vs sc -1.425). gH = exp(-kH*z_help):
        #   help preserved (z>0) -> gH<1 -> SLOWS CD8 decay (durability bonus)
        #   help eroded  (z<0) -> gH>1 -> SPEEDS CD8 decay (durability PENALTY)
        z_help = self.sc_help_z if not _nan(self.sc_help_z) else self.help_z
        gH = 1.0 if _nan(z_help) else float(np.exp(-gp.kH * z_help))
        # typed suppression: use decomposed z's if present, else split supp_z
        if _nan(self.supp_z):
            gS_size = gS_rate = 1.0
        else:
            zs = self.supp_size_z if not _nan(self.supp_size_z) else self.frac_size*self.supp_z
            zr = self.supp_rate_z if not _nan(self.supp_rate_z) else (1-self.frac_size)*self.supp_z
            gS_size = float(np.exp(gp.kS_size * zs))
            gS_rate = float(np.exp(gp.kS_rate * zr))
        # NEW: exhaustion & proliferation gains (COMPLEMENT the veto gate; modest sensitivities).
        #   gX>1 (exh_z>0) => faster CD8 contraction; gP>1 (prolif_z>0) => larger CD8 effector pool.
        gX = 1.0 if _nan(self.exh_z)    else float(np.exp(gp.kX * self.exh_z))
        gP = 1.0 if _nan(self.prolif_z) else float(np.exp(gp.kP * self.prolif_z))
        return dict(gE=gE, gC=gC, gH=gH, gS_size=gS_size, gS_rate=gS_rate, gX=gX, gP=gP)

    def expr_vec(self):
        return np.array([self.expr_CD8, self.expr_CD4c, self.expr_Treg], float)


@dataclass
class GainParams:
    """z -> gain sensitivity constants (swept in the sensitivity analysis).
    Typed suppression uses two constants so the compartment-size and per-cell-
    rate pathways can be weighted independently in the typed-term sensitivity."""
    kE: float = 0.11
    kC: float = 0.45
    kS_size: float = 0.45      # suppression z -> Treg compartment-size gain
    kS_rate: float = 0.45      # suppression z -> per-Treg suppressive-rate gain
    kH: float = 0.20
    w_onset: float = 0.35      # weight on EARLY (8hr) CRS rate for acute-CRS liability
                               # (0 => 48hr snapshot only; 1 => early rate dominates)
    kX: float = 0.12           # exhaustion z -> CD8-decay-acceleration gain (gX=exp(+kX*exh_z));
                               #   >0 exhaustion SPEEDS CD8 contraction. Modest (0.12) so it complements
                               #   the binary EXH veto, never rescues/overrides it.
    kP: float = 0.10           # proliferation z -> CD8-expansion gain (gP=exp(+kP*prolif_z));
                               #   >0 expands the CD8 effector pool. Modest, complement-only.


# Costim-receptor expression PRIOR across {CD8, CD4-conv, Treg}.
# SOURCING: the project brief states Tregs express CD28 and 4-1BB at high levels
# and that Treg-mediated costim toxicity is "sharpest for CD28 and for 4-1BB";
# OX40/GITR are canonical Treg-high TNFRSF costim receptors. These set a COARSE,
# 3-level (low=0.5 / mid=1.0 / high=1.6) prior; everything else stays neutral.
# This is a LITERATURE/PROJECT-STATED PRIOR, explicitly overridable by the PI's
# RNA-to-surface-protein (QIFIKIT) layer. The nomination is checked with this
# prior ON and OFF (uniform) to report whether it depends on the prior.
_EXPR_PRIOR = {
    # gene:      (CD8,  CD4c, Treg)
    "CD28":      (1.0,  1.0,  1.6),   # broad; Treg-high (brief)
    "TNFRSF9":   (1.6,  1.0,  1.6),   # 4-1BB: activated-CD8 high AND Treg-high (brief)
    "TNFRSF4":   (0.8,  1.2,  1.6),   # OX40: Treg/activated-CD4 > CD8
    "TNFRSF18":  (0.8,  1.0,  1.6),   # GITR: Treg-high
    "ICOS":      (0.8,  1.4,  1.4),   # CD4/Treg biased
    "TNFRSF8":   (1.2,  1.0,  0.8),   # CD30: activation-restricted, not Treg-enriched
}

def set_expression_prior(arm: ArmScores, gene: str | None = None):
    """Populate arm.expr_* from the coarse literature/project prior (if present).
    Genes not in the prior stay uniform (neutral)."""
    g = gene or arm.name
    if g in _EXPR_PRIOR:
        arm.expr_CD8, arm.expr_CD4c, arm.expr_Treg = _EXPR_PRIOR[g]
    return arm


# ============================================================================
# LAYER 2 — trimer / synapse engagement (3-body, per compartment)
# ============================================================================
def trimer_bell(conc, C_opt, tmax=1.0):
    """Reduced ternary-complex abundance, bell-shaped (prozone) in concentration.
    Closed-form small-model of the Douglass 2013 three-body equilibrium:
    T(C) = tmax * (C/Copt) / (1 + C/Copt)^2, peaking at C = Copt (the prozone apex).
    """
    x = np.maximum(conc, 0.0) / C_opt
    return tmax * x / (1.0 + x) ** 2

def engagement(conc, C_opt, gE, prozone_rescue, tmax=1.0):
    """CD3-trimer engagement with costim anti-prozone rescue. The costim arm
    sustains activation at high concentration (killing plateaus rather than
    following the bell), scaled by the arm's effector gain."""
    bell = trimer_bell(conc, C_opt, tmax)
    x = np.maximum(conc, 0.0) / C_opt
    sat = tmax * x / (1.0 + x)                       # monotone-saturating (Emax) form
    r = prozone_rescue * np.clip((gE - 1.0) / 2.0, 0.0, 1.0)
    return (1.0 - r) * bell + r * (0.25 * sat)


def bivalent_binding(Cfree, R_CD3, R_TAA, p: "PBPK"):
    """Emergent target binding per tissue via a Rhoden-2016 bivalent-binding QSS.

    The engager's two arms (CD3, TAA) bind INDEPENDENTLY with monovalent
    dissociation constants KD_CD3, KD_TAA. After one arm binds, the second arm
    accesses its target within a constrained local volume -> an AVIDITY
    enhancement that scales with LOCAL RECEPTOR DENSITY (Rhoden's antigen-density
    dependence). We solve the quasi-steady-state and return the BOUND-COMPLEX
    CONCENTRATIONS (saturating in Cfree, capped by the finite receptor pool), so
    that internalization of bound complex is a SATURABLE elimination -> TMDD
    emerges with the correct sign (fractional clearance falls as dose rises).

    Returns (Cb_mono, Cb_bridge): local concentration of singly-bound and
    bridged (avidity-enhanced) drug-receptor complex. Both -> internalization.

    Reference: Rhoden JJ, Dyas GL, Wroblewski VJ (2016) J Biol Chem 291:11337-47,
    DOI 10.1074/jbc.M116.714287 (bispecific bivalent binding to cell-surface
    targets; antigen density / affinity / expression-ratio effects; PMID 27022022).
    Cross-cell TCE caveat: avidity cooperativity is a WITHIN-CELL effect (CD3
    density on T cells; TAA density on tumor cells), so it enters through the
    per-target local densities R_CD3, R_TAA separately, not as a single
    cross-cell bridge constant.
    """
    Cfree = np.maximum(Cfree, 0.0)
    # --- monovalent single-site binding of each arm, with VALENCY-dependent
    #     avidity. When the construct has >=2 arms against a given cell's target,
    #     multivalent binding raises the EFFECTIVE affinity of that arm to that
    #     cell (lower apparent KD), scaled by local receptor density (Rhoden:
    #     avidity gain grows with antigen density). A single arm (n=1) gets no
    #     enhancement. This makes binding format-aware (1+1 / 2+1 / 1+2 / 2+2).
    # EXACT geometric Rhoden avidity via the SATURATING receptor-limited QSS (NOT the linear
    # limit — Ag_eff here is 1e4-1e5 nM, so the 2nd-arm bridge is near-saturated and the
    # enhancement CANNOT run to 1e4x; it saturates as the receptor pool fills). The bound
    # fraction of a bivalent arm on a cell of density R_local (nM) is solved from
    #   (2*Ag_eff*C/KD^2) x^2 + (1+2C/KD) x - 1 = 0,  x = Rfree/R0,
    # bound = R_local*(2 x C/KD + Ag_eff x^2 C/KD^2). n_arm<=1 => monovalent R*C/(KD+C).
    def _geo_ageff(rec_pc, r_cell_um, r_Ab_nm):
        # EXACT validated preprint formula (rhoden_tmdd.Ag_eff_nM); cell-count cancels.
        _cc = 1e9
        SA_cell = 4.0*np.pi*r_cell_um**2; r_Ab_um = r_Ab_nm*1e-3
        SA_Ab = np.pi*r_Ab_um**2; V_Ab = (2.0/3.0)*np.pi*r_Ab_um**3
        Ag_bulk = rec_pc * _cc / 6.02214076e23 * 1e9             # receptors_to_nM
        Am_cell = Ag_bulk / _cc; Am_SA = Am_cell/SA_cell; Ag_SA = Am_SA*SA_Ab
        return Ag_SA / V_Ab * 1e15                               # nM
    def _bound_arm(n_arm, R_local, C, rec_pc, r_cell_um, KD):
        R_local = np.maximum(R_local, 0.0); C = np.maximum(C, 0.0)
        if n_arm <= 1:
            return R_local * C/(KD + C + 1e-12)                  # monovalent single-site
        ageff = _geo_ageff(rec_pc, r_cell_um, p.r_Ab_nm_avidity)
        a = 2.0*ageff*C/KD**2; b = 1.0 + 2.0*C/KD
        x = np.where(a > 1e-30, (-b + np.sqrt(b*b + 4.0*a))/(2.0*a + 1e-300), 1.0/np.maximum(b,1e-30))
        x = np.clip(x, 0.0, 1.0)
        return R_local*(2.0*x*C/KD + ageff*x*x*C/KD**2)          # S + D (saturating)
    Cb3 = _bound_arm(p.n_CD3, R_CD3, Cfree, p.rec_per_cell_CD3, p.r_cell_um_CD3, p.KD_CD3)
    CbT = _bound_arm(p.n_TAA, R_TAA, Cfree, p.rec_per_cell_TAA, p.r_cell_um_TAA, p.KD_TAA)
    # --- cross-cell bridge (the cytotoxic synapse): a drug molecule bound by its
    #     CD3 arm to a T cell AND its TAA arm to a tumor cell. This is the
    #     productive trimer; it is the co-localised (min) of the two bound pools,
    #     weighted by 2nd-arm occupancy. Distinct from within-cell avidity above.
    # AF3 REACH as an EFFECTIVE cross-cell AFFINITY shift on the 2nd (bridging) arm.
    # Parallel to the within-cell avidity above (KD_eff = KD/gain): a format whose
    # arms cannot span the ~15 nm synapse cleft has a WEAKER apparent cross-cell KD
    # for the bridging step -> KD_bridge_eff = KD_TAA / reach_gate. The bridge then
    # forms with lower occupancy at a given Cfree but is fully recoverable at high
    # Cfree (mass action) — the reduced bridging is EMERGENT from the shifted
    # occupancy, not an imposed cut. reach_gate=1.0 => identity (baseline). Applied
    # to the bridge only; within-cell mono binding (Cb3, CbT) crosses no cleft.
    f_reach = np.clip(getattr(p, "reach_gate", 1.0), 1e-3, 1.0)
    KD_bridge_eff = p.KD_TAA / f_reach                 # short reach => weaker apparent bridge affinity
    occ = Cfree/(KD_bridge_eff + Cfree + 1e-12)        # fractional 2nd-arm (cross-cell) engagement
    Cb_bridge = np.minimum(Cb3, CbT) * occ
    Cb_mono   = np.maximum(Cb3 + CbT - 2.0*Cb_bridge, 0.0)   # unbridged remainder (mass-conserved)
    return Cb_mono, Cb_bridge


def ternary_equilibrium(C, RtotA, RtotB, KD1, KD2, alpha=1.0):
    """EXACT bridged ternary-complex (trimer) abundance at binding equilibrium.

    Schropp, Khot, Shah & Koch (2019) CPT:PSP 8:177-187, DOI 10.1002/psp4.12369
    (equilibrium-binding / EB model, Eqs. 26-30 & 33). A bispecific of free
    concentration C bridges receptor A (e.g. CD3 on a T cell, total R_totA) and
    receptor B (e.g. TAA on a tumour cell, total R_totB); KD1/KD2 are the two
    monovalent dissociation constants and alpha the cross-arm cooperativity
    (alpha=1 for independent arms). Returns the ternary complex

        RC_AB = C * R_A(C) * R_B(C) / (alpha*KD1*KD2)                    (Eq. 33)

    with the free receptors R_A, R_B solved from the closed-form QE relations
    (Eqs. 26-30). This is the REAL immune-synapse trimer that drives engager PD:
    the prozone (hook) is EMERGENT — at high C both arms saturate as BINARY
    complexes (R_A,R_B -> 0) so the trimer collapses (Eq. 35, lim_{C->inf}=0),
    reproducing the bell WITHOUT any fitted optimum. Cross-molecule magnitude
    differences arise purely from KD1,KD2,R_totA,R_totB (affinity + density).

    Scalar or vectorised in C. Units: all concentrations nM; RC_AB in nM.
    """
    scalar = np.isscalar(C)
    C = np.atleast_1d(np.asarray(C, float))
    out = np.zeros_like(C)
    if RtotA > 0.0 and RtotB > 0.0:
        aKK = alpha * KD1 * KD2
        for i, c in enumerate(C):
            if c <= 0.0:
                continue
            a = (1.0 + c / KD2) * c / aKK                                  # Eq. 28
            b = c * (RtotA - RtotB) / aKK + (1.0 + c / KD2) * (1.0 + c / KD1)  # Eq. 29
            d = -RtotB * (1.0 + c / KD1)                                   # Eq. 30
            disc = b * b - 4.0 * a * d
            RB = (-b + np.sqrt(disc if disc > 0.0 else 0.0)) / (2.0 * a)   # Eq. 27
            RA = RtotA / (1.0 + c / KD1 + RB * c / aKK)                    # Eq. 26
            out[i] = c * RA * RB / aKK                                     # Eq. 33
    return float(out[0]) if scalar else out


# ============================================================================
# Master parameter bundle
# ============================================================================
@dataclass
class Params:
    pbpk: PBPK = field(default_factory=PBPK)
    cell: CellParams = field(default_factory=CellParams)
    cyto: CytokineParams = field(default_factory=CytokineParams)
    gain: GainParams = field(default_factory=GainParams)
    # dosing: IV bolus is dose->A_pl(0). A CONTINUOUS INFUSION (e.g. blinatumomab) is
    # specified here so it is native to the RHS: inf_rate (mg/day into the plasma pool)
    # for t <= inf_duration. Default 0 => pure bolus, which is what EVERY window-scored
    # arm uses, so the infusion-driven fresh-delivery tail term (Ecum) is inactive during
    # scoring and the nomination is structurally independent of it.
    inf_rate: float = 0.0        # continuous infusion rate (mg/day) into plasma
    inf_duration: float = 0.0    # infusion duration (days)
    # engagement
    C_opt: float = 1.0
    prozone_rescue: float = 0.6
    # ---- MECHANISTIC systemic-CRS target dependence (real-unit runs only) ----
    # TCE-driven CRS requires the CD3:drug:target trimer, so systemic cytokine
    # release scales with the ACCESSIBLE systemic target burden (splenic/circulating
    # target reachable by circulating T cells), NOT with a fixed fraction. This is
    # why B-cell targets (CD20/CD19: huge splenic burden) drive severe CRS while
    # marrow-restricted plasma-cell targets (BCMA/GPRC5D) drive mild CRS, and why
    # step-up dosing / debulking lower CRS. When sys_target_burden_nM is None
    # (default, abstract window score), behaviour is IDENTICAL to the fixed
    # sys_eng_frac form (backward-compatible). When set (per-molecule, by the
    # unified model), eng_sys uses an Emax in accessible burden with half-saturation
    # sys_burden_K_nM. burden_factor is normalised to the anchor so the mosunetuzumab
    # calibration is preserved.
    sys_target_burden_nM: float = None   # accessible systemic target burden (per-molecule); None => fixed sys_eng_frac
    sys_burden_K_nM: float = 50.0        # half-saturation of systemic trimer in accessible target burden.
                                         # Literature-order (not panel-fit): the accessible-target scale that
                                         # separates hematologic B-cell targets (severe CRS) from marrow-restricted
                                         # plasma-cell / solid-tumor targets (mild CRS). Encodes the established
                                         # "total accessible antigen load drives IL-6/CRS" relationship
                                         # (Synnott 2025 Cancer 10.1002/cncr.70069; AAPS J 2025 10.1208/s12248-025-01177-9).
    sys_burden_anchor_nM: float = 746.554767  # anchor burden = the UnifiedModel's exact class-computed
                                         # accessible burden for mosunetuzumab (CD20 splenic+marrow B-cell pool,
                                         # access_w spleen/bone/tumor=1.0, large_int=0.3). At this burden the
                                         # Emax reproduces the fixed sys_eng_frac EXACTLY (max|Δ| over all state
                                         # derivatives = 0 vs the None/abstract path), so the window/nomination
                                         # score is byte-identical to the committed model. Verified in-session.
    k_lic: float = 0.02       # trimer half-saturation for KILLING LICENSE: TCE
                              # cytotoxicity requires the drug-formed synapse, so
                              # kill scales with eng_tum/(eng_tum+k_lic) (no drug
                              # => no redirected lysis => untreated tumor grows).
    # costim boost strength: how much signal-2 multiplies activation where the
    # costim receptor is expressed (0 => CD3-only backbone behaviour)
    costim_boost: float = 0.7
    # tolerability + efficacy horizons
    efficacy_target: float = 0.35
    crs_severity_ceiling: float = None   # set by calibrate_ceiling (IL-6-based)
    t_read: float = 42.0
    tspan: tuple = (0.0, 56.0)
    crs_tol_factor: float = 1.15
    liver_ceiling: float = None      # set by calibrate_ceiling (liver-injury units)
    liver_tol_factor: float = 1.30   # tolerable liver injury = factor x backbone-at-MED
    # composite-window diagnostic weights (small; E_tol is the primary term)
    w_treg: float = 0.15     # residual Treg burden at operating dose (tie-breaker)
    w_liver: float = 0.15    # residual liver-tox at operating dose (tie-breaker)


# ============================================================================
# LAYER 3+4 — the multi-cell ODE system
# ============================================================================
def rhs(t, y, P: Params, arm: ArmScores, gains: dict):
    p, cp, cy = P.pbpk, P.cell, P.cyto
    gE  = gains["gE"];         gC     = gains["gC"]
    gSsize = gains["gS_size"]; gSrate = gains["gS_rate"]; gH = gains["gH"]
    gX = gains.get("gX", 1.0); gP = gains.get("gP", 1.0)   # exhaustion / proliferation gains
    y = np.asarray(y, float)

    # ---- unpack states ----
    A_pl = y[I_PL]                               # plasma pool (incl. fast-equilibrated organ vascular)
    A_is = _cis_amt(y)                           # per-organ interstitial amounts
    A_ly = y[I_LY]                               # lymph pool
    CD8_b = y[IDX["CD8_b"]]; CD8_t = y[IDX["CD8_t"]]
    CD4c_b= y[IDX["CD4c_b"]];CD4c_t= y[IDX["CD4c_t"]]
    Treg_b= y[IDX["Treg_b"]];Treg_t= y[IDX["Treg_t"]]
    Mye  = y[IDX["Mye"]];  Tumor= y[IDX["Tumor"]]
    cTNF = y[IDX["cTNF"]]; cIFN = y[IDX["cIFN"]]; cIL2 = y[IDX["cIL2"]]
    cIL6 = y[IDX["cIL6"]]; Liv  = y[IDX["Liv"]]
    Dcyto = y[IDX["Dcyto"]]
    Ecum  = y[IDX["Ecum"]]

    # ---- L1 FULL-BODY PBPK (lymph-flow-limited; QSS vascular lumped into plasma) ----
    # Antibody tissue distribution is governed by LYMPH flow, not blood flow
    # (tight endothelium; Shah DK & Betts AM 2012, J Pharmacokinet Pharmacodyn
    # 39:67-86, DOI 10.1007/s10928-011-9232-2, PMID 22143261). Organ vascular spaces
    # equilibrate with plasma ~instantly (Q/Vv up to ~1e5/day), so they are lumped
    # into the plasma pool and the DYNAMIC distribution is 2-pore convective
    # extravasation into each organ INTERSTITIUM (reflected by sigma_V) and lymph
    # return (sigma_L), all lymph collecting to a pool that drains back to plasma.
    # Vc = PHYSIOLOGICAL plasma volume (3.1 L). The antibody is measured as a plasma
    # concentration; organ vascular spaces hold drug at plasma concentration and are
    # part of that same circulating plasma, so they are NOT added as a separate central
    # volume (doing so put Vc at 5.16 L and depressed C0 ~40% below the clinical value).
    # This is the physiological Vc that 2-cpt clinical fits recover (~3.0-3.5 L).
    V_pl = V_PLASMA
    C_pl = A_pl / V_pl
    Cis  = A_is / PB.Vis                          # per-organ interstitial conc
    # 2-pore convective distribution, scaled by k_dist (distribution RATE knob; cancels
    # from the SS interstitial:plasma ratio so it sets the alpha-phase DURATION only).
    J_extrav = p.k_dist*PB.L*(1.0 - PB.sigV)*C_pl # plasma -> each interstitium (amount/day)
    J_return = p.k_dist*PB.L*(1.0 - p.sigL)*Cis   # each interstitium -> lymph (amount/day)

    # ---- EMERGENT TMDD (Rhoden 2016 bivalent binding, per tissue) ----
    # Interstitial drug engages targets where expressed: CD3 capacity ~ tissue
    # T-cell density, TAA capacity tumor-only. Internalization of the bound
    # species IS the target-mediated clearance -> nonlinear, saturable, largest
    # in tumor (TAA bridge) and T-cell-dense tissues. Nothing imposed.
    R_CD3 = p.Rcap_CD3*PB.tcell
    R_TAA = p.Rcap_TAA*PB.taa
    Cb_mono, Cb_bridge = bivalent_binding(Cis, R_CD3, R_TAA, p)
    # SATURABLE target-mediated elimination (amount/day per tissue): internalize
    # the bound complex. Because Cb saturates at the finite receptor pool, the
    # fractional clearance falls as dose rises => emergent TMDD nonlinearity.
    tmdd_sink = (p.kint_mono*Cb_mono + p.kint_cplx*Cb_bridge)*PB.Vis

    # Plasma elimination has TWO mechanistic routes, both format-dependent:
    #  (1) endosomal catabolism of the pinocytosed-but-NOT-FcRn-salvaged fraction.
    #      For a no-Fc construct eff_fFcRn()=0, so the ENTIRE pinocytosed pool is
    #      catabolised (no recycling) -> much faster than an Fc-bearing mAb.
    #  (2) size-gated renal filtration + proteolysis, ~0 for a 150 kDa IgG but
    #      dominant for a small (<~69 kDa) non-Fc scaffold like a BiTE.
    k_cat = p.CLup*(1.0 - p.eff_fFcRn()) + p.k_renal()
    infusion = P.inf_rate if t <= P.inf_duration else 0.0   # continuous IV infusion (mg/day)
    dA_pl = infusion - k_cat*A_pl - J_extrav.sum() + p.k_lymph_return*A_ly
    dA_is = J_extrav - J_return - tmdd_sink
    dA_ly = J_return.sum() - p.k_lymph_return*A_ly

    # ---- L2 engagement: two DISTINCT engagement sites ----
    # (a) TUMOR trimer: full CD3-TCE-TAA ternary complex in tumor interstitium
    #     (TAA present) -> drives CD8 KILLING. This is the EFFICACY site.
    # (b) SYSTEMIC engagement: CD3 is lineage-blind, so the drug engages the
    #     circulating (blood) T-cell pool everywhere, but WITHOUT TAA there is no
    #     full trimer -> weak (sys_eng_frac). This CD4-dominated systemic
    #     engagement is the CRS ORIGIN (the TOXICITY site). Separating the two is
    #     what lets the model dissect tox (CD4/systemic) from efficacy (CD8/tumor).
    # SPATIAL hook (a): per-organ exposure multiplier scales the interstitial conc
    # that drives engagement. Default spatial_exposure = 1.0 for all organs => the
    # tumor conc (and hence eng_tum) is exactly the well-mixed value. When the
    # spatial lane delivers a tumor penetration penalty <1, the tumor trimer (and
    # thus killing) is down-weighted; other organs' exposure enters via Cis below.
    C_tum = Cis[PB.i_tumor] * PB.spatial_exposure[PB.i_tumor]
    A_liv = A_is[PB.i_liver]                       # liver interstitial drug amount (hepatotox)
    # SPATIAL hook (c): inter-arm reach gate down-weights productive trimer when a
    # format's CD3<->TAA span can't bridge the synapse cleft. Default reach_gate=1.0
    # => trimer unchanged. Applied to the tumor trimer only (the killing license).
    # AF3 REACH as an EFFECTIVE-POTENCY shift (not an output multiplier). A format
    # whose CD3<->TAA arms cannot span the ~15 nm synapse cleft does not form ZERO
    # trimer — it needs MORE drug to form the same trimer (mass action). So reach
    # enters as a right-shift of the trimer dose-response: C_opt_eff = C_opt/reach_gate.
    # reach_gate=1.0 => identity (baseline). reach_gate<1 => potency drops, but a
    # saturating dose still recovers killing. The reduced-kill-at-fixed-dose result
    # is therefore EMERGENT from the shifted curve + the tumor-growth competition in
    # dTumor, not imposed. (Elimination side handled in bivalent_binding via KD.)
    C_opt_eff = P.C_opt / max(P.pbpk.reach_gate, 1e-3)
    # ---- FORMAT-AWARE MECHANISTIC PD DRIVER (efficacy site) --------------------------
    # When pd_driver_mode is set (real-unit validation runs), the tumour-site PD driver is
    # the REAL complex species for the molecule's binding format, normalised to trimer_ref_nM
    # so the validated downstream cascade is preserved. reach enters as an effective 2nd-arm
    # KD shift (mass action), exactly as in bivalent_binding. Legacy bell kept as default
    # (pd_driver_mode is None) so the abstract window-scoring / nomination is byte-unchanged.
    mode = getattr(p, "pd_driver_mode", None)
    if mode == "trimer":
        KD_TAA_eff = p.KD_TAA / max(p.reach_gate, 1e-3)          # short reach => weaker apparent bridge
        R_TAA_tum = R_TAA[PB.i_tumor]                            # real tumour TAA pool (nM)
        tri_tum = ternary_equilibrium(C_tum, p.R_CD3_pd_nM, R_TAA_tum, p.KD_CD3_pd, KD_TAA_eff)
        eng_tum = float(np.clip(tri_tum / max(p.trimer_ref_nM, 1e-9), 0.0, 4.0))
    elif mode == "occupancy":
        eng_tum = C_tum / (C_tum + p.occ_KD_nM + 1e-12)          # receptor occupancy (checkpoint/costim)
    elif mode == "soluble":
        eng_tum = C_tum / (C_tum + p.KD_TAA + 1e-12)             # fractional target neutralization
    elif mode == "binary":
        Cb_m, Cb_b = bivalent_binding(np.array([C_tum]), np.array([R_TAA[PB.i_tumor]]),
                                      np.array([R_TAA[PB.i_tumor]]), p)
        eng_tum = float(np.clip((Cb_m[0]+Cb_b[0]) / max(p.trimer_ref_nM, 1e-9), 0.0, 4.0))
    else:
        eng_tum = engagement(C_tum, C_opt_eff, gE, P.prozone_rescue)  # LEGACY bell (abstract)
    # Systemic engagement drives blood T cells via CD3 (no TAA). The instantaneous
    # engagement uses the same trimer/occupancy bell as the tumor, which reproduces both
    # clinical ACUTE peaks: blina's continuous-infusion Css (~1.1 ug/mL) sits near the
    # engagement optimum so it peaks HIGH (clinical 717 > mosun 570) despite lower conc
    # than mosun's bolus Cmax (the continuous-infusion effect). The SUSTAINED-tail
    # separation between a clearing bolus (mosun) and a maintained infusion (blina) is
    # carried by Ecum, a slow cumulative-engagement/adaptation state (see below), NOT by
    # instantaneous concentration (which is nearly identical for the two in the tail).
    # Systemic engagement efficiency: fixed fraction (abstract window score, default)
    # OR mechanistic accessible-target-burden Emax (per-molecule real-unit runs). The
    # burden form makes systemic CRS TARGET-DEPENDENT: the CD3:drug:target trimer needs
    # accessible systemic target, so a B-cell target (large splenic burden) saturates the
    # systemic trimer while a marrow-restricted plasma-cell target does not. Normalised so
    # the anchor burden reproduces sys_eng_frac exactly (mosunetuzumab calibration held).
    if P.sys_target_burden_nM is None:
        sys_eff = cp.sys_eng_frac
    else:
        b = max(P.sys_target_burden_nM, 0.0)
        emax_b   = b / (b + P.sys_burden_K_nM)
        emax_anc = P.sys_burden_anchor_nM / (P.sys_burden_anchor_nM + P.sys_burden_K_nM)
        sys_eff  = cp.sys_eng_frac * (emax_b / emax_anc)   # ==sys_eng_frac at anchor burden
    # ---- FORMAT-AWARE MECHANISTIC PD DRIVER (systemic CRS site) ----------------------
    # For a TCE the systemic storm is driven by the CD3xTAA trimer formed against the
    # ACCESSIBLE systemic target burden (circulating/splenic target) — so a B-cell target
    # (large splenic pool) seeds more CRS than a marrow-restricted plasma-cell target, and
    # a naked mAb (no CD3 arm) forms no systemic trimer at all. Normalised to the same
    # trimer_ref, scaled by sys_eff (which already carries the accessible-burden Emax).
    if mode == "trimer":
        b_sys = max(P.sys_target_burden_nM, 0.0) if P.sys_target_burden_nM is not None else 0.0
        tri_sys = ternary_equilibrium(C_pl, p.R_CD3_pd_nM, b_sys, p.KD_CD3_pd, p.KD_TAA)
        eng_sys = float(np.clip(tri_sys / max(p.trimer_ref_nM, 1e-9), 0.0, 4.0))
    elif mode in ("occupancy", "soluble", "binary"):
        eng_sys = sys_eff * (C_pl / (C_pl + p.occ_KD_nM + 1e-12))   # systemic occupancy (no bridged trimer)
    else:
        eng_sys = sys_eff * engagement(C_pl, P.C_opt, gE, P.prozone_rescue)  # LEGACY bell (abstract)
    e8,e4,er = arm.expr_vec()
    b8 = 1.0 + P.costim_boost * e8 * _costim_unit(gE)
    b4 = 1.0 + P.costim_boost * e4 * _costim_unit(gE)
    br = 1.0 + P.costim_boost * er * _costim_unit(gE)

    # ---- L3 activation, split by SITE ----
    # CD8: activated where the tumor trimer forms (needs TAA) -> tumor-pool CD8.
    act_CD8_t  = cp.kact_CD8  * eng_tum * b8 * gE * gP   # gP: agonism-driven CD8 proliferation/expansion
    # CD4-conv: activated BOTH at the tumor (trimer) AND systemically (blood pool,
    # lineage-blind CD3). The systemic CD4 activation is the CRS engine.
    act_CD4c_t = cp.kact_CD4c * eng_tum * b4
    act_CD4c_b = cp.kact_CD4c * eng_sys * b4
    # TYPED suppression term #1 (COMPARTMENT SIZE): gS_size scales Treg EXPANSION.
    act_Treg_t = cp.kact_Treg * eng_tum * br * gSsize
    act_Treg_b = cp.kact_Treg * eng_sys * br * gSsize

    # ---- trafficking: blood <-> tumor for each mobile T-cell type ----
    # infiltration is enhanced by local tumor engagement (inflamed tumor recruits
    # activated T cells). Egress rate per type is set so the UNTREATED baseline is
    # an exact fixed point: k_egress = k_infil*blood_0 / tumor_0.
    def _traffic(blood, tumor, b0, t0, infil_gain):
        k_in  = cp.k_infil * (1.0 + infil_gain * eng_tum)   # blood -> tumor (engagement-boosted)
        k_out = cp.k_infil * b0 / max(t0, 1e-6)             # tumor -> blood (baseline fixed point)
        flux  = k_in * blood - k_out * tumor
        return flux
    fl_CD8  = _traffic(CD8_b,  CD8_t,  cp.CD8_b0,  cp.CD8_t0,  cp.infil_gain_CD8)
    fl_CD4c = _traffic(CD4c_b, CD4c_t, cp.CD4c_b0, cp.CD4c_t0, cp.infil_gain_CD4c)
    fl_Treg = _traffic(Treg_b, Treg_t, cp.Treg_b0, cp.Treg_t0, cp.infil_gain_Treg)

    # help -> CD8 PERSISTENCE. Tumor CD4-conv provides help that sustains tumor
    # CD8 (slows its decay). The costim arm modulates the QUALITY of that help via
    # gH (single-cell-resolved): gH<1 preserves/raises help (slower CD8 decay,
    # durability bonus); gH>1 ERODES help (faster CD8 decay, durability PENALTY).
    # help_mult spans BELOW 1 (bonus) and ABOVE 1 (penalty) so an arm like CD30
    # that shuts down the CD4-help lineage pays a real persistence cost on CD8 --
    # the mechanism by which single-cell help-erosion lowers the therapeutic
    # window even though it also lowers suppression.
    help_avail = np.tanh(max(CD4c_t - cp.CD4c_t0, 0.0))    # available CD4 help (abundance)
    help_mult  = 1.0 - cp.help_gain*help_avail*(1.0 - gH)  # gH<1 -> <1 (bonus); gH>1 -> >1 (penalty)
    help_mult  = float(np.clip(help_mult, 0.4, 1.8))
    help_factor = help_mult
    # help also licenses CD8 CYTOTOXIC QUALITY: CD4 help (IL-2/IL-21, licensing)
    # sustains CD8 killing potency and functional avidity, not just cell number.
    # An arm that ERODES help (gH>1) therefore lowers per-CD8 kill EFFICIENCY even
    # when CD8 numbers hold — this is how help-erosion (CD30 lineage-shutdown)
    # costs efficacy in the sub-curative regime and separates it from a
    # help-sparing arm (4-1BB) that clears the same tumor and holds it. Bounded.
    help_quality = float(np.clip(1.0 - cp.help_kill_gain*(gH - 1.0), 0.55, 1.25))

    # CD8 pools: blood activates none (needs TAA), only traffics + decays; tumor
    # pool is activated by the trimer, receives infiltration, is helped, decays.
    dCD8_b  = -fl_CD8  - cp.kdec_CD8*(CD8_b - cp.CD8_b0)
    dCD8_t  = act_CD8_t + fl_CD8 - cp.kdec_CD8 * gX * help_factor * CD8_t   # gX: exhaustion-accelerated CD8 decay
    # CD4-conv pools: blood activated systemically (CRS), tumor activated by trimer
    dCD4c_b = act_CD4c_b - fl_CD4c - cp.kdec_CD4c*(CD4c_b - cp.CD4c_b0)
    dCD4c_t = act_CD4c_t + fl_CD4c - cp.kdec_CD4c*(CD4c_t - cp.CD4c_t0)
    # Treg pools
    dTreg_b = act_Treg_b - fl_Treg - cp.kdec_Treg*(Treg_b - cp.Treg_b0)
    dTreg_t = act_Treg_t + fl_Treg - cp.kdec_Treg*(Treg_t - cp.Treg_t0)

    # ---- L3b cytokines: SYSTEMIC storm from BLOOD T cells (CD4-conv dominant),
    # plus tumor-local contribution; IL-6 from myeloid. CRS is a SYSTEMIC readout,
    # so the blood CD4-conv pool is the dominant cytokine source (the CD4 problem).
    eff_pot = gE ** cy.cyto_eff_coupling
    CD4c_tot = CD4c_b + 0.5*CD4c_t     # systemic-weighted (blood dominant for CRS)
    CD8_tot  = CD8_b  + 0.5*CD8_t
    # ACUTE first-dose CRS driver: systemic engagement gated by desensitization.
    # Dcyto in [0,1); (1-Dcyto) is the fraction of cytokine-responsiveness remaining.
    # burst_drive is large at first engagement, self-limits as Dcyto builds -> the
    # storm spikes within ~1 day and resolves within ~2-3 days despite persistent drug.
    resp = max(1.0 - Dcyto, 0.0)
    burst_drive = cy.k_cyto_burst * eng_sys * resp * gC * eff_pot
    # SUSTAINED envelope (the slow term the window score integrates): per-cell output
    # from the expanding blood CD4-conv pool (the CD4 problem). Kept as before but
    # scaled so the plateau sits well below the acute peak (Hosseini: subsequent
    # peaks ~order of magnitude lower than the first).
    # Sustained production tracks ONGOING systemic engagement DIRECTLY (not the acute
    # desensitization resp): storm cytokines come from ACTIVELY-ENGAGED T-cells, not
    # resting expanded ones. Desensitization is an ACUTE first-dose phenomenon (it
    # self-limits the spike via burst_drive); it must NOT also suppress the chronic
    # low-level production from continuing engagement. Because eng_sys is monotone in
    # SUSTAINED tail gate. Two additive contributions:
    #  - sust_floor: the arm-dependent baseline sustained tail. Multiplied by gC below, so
    #    this is the term the CRS counter-screen (gC) acts on and the window score reads
    #    for a BOLUS-dosed arm (Ecum=0). A cleared bolus (mosun) resolves to this floor.
    #  - Ecum**hill: the FRESH-DELIVERY contribution, nonzero ONLY under a continuous
    #    infusion (blina), which sustains new-synapse formation and thus a higher tail.
    # This makes mosun (bolus) resolve low and blina (infusion) hold high — from the
    # dosing schedule, not instantaneous conc — while leaving the scored (bolus) arms'
    # tail = gC*sust_floor, so the nomination is structurally independent of the infusion.
    sust_gate = cy.sust_floor + (Ecum ** cy.sust_hill)
    sust_TNF = (cy.p_TNF_CD4c*CD4c_tot + cy.p_TNF_CD8*CD8_tot) * gC * eff_pot * sust_gate
    sust_IFN = (cy.p_IFN_CD4c*CD4c_tot + cy.p_IFN_CD8*CD8_tot) * gC * eff_pot * sust_gate
    sust_IL2 =  cy.p_IL2_CD4c*CD4c_tot * gC * eff_pot * sust_gate
    prod_TNF = cy.f_TNF_burst*burst_drive + sust_TNF
    prod_IFN = cy.f_IFN_burst*burst_drive + sust_IFN
    prod_IL2 = cy.f_IL2_burst*burst_drive + sust_IL2
    dcTNF = prod_TNF - cy.d_TNF*cTNF
    dcIFN = prod_IFN - cy.d_IFN*cIFN
    dcIL2 = prod_IL2 - cy.d_IL2*cIL2
    # Myeloid is activated by the T-cell cytokines it amplifies (cTNF + cIFN), which
    # now themselves track ongoing engagement — so the myeloid relay follows the same
    # emergent bolus-vs-infusion separation without any explicit resp/dosing gate: with
    # a clearing bolus (mosun) TNF/IFN fall and myeloid relaxes to baseline; with a
    # maintained infusion (blina) they persist and myeloid holds a higher residual.
    act_Mye = cy.k_mye_act * (cTNF + cIFN)
    dMye = act_Mye - cp.kdec_Mye*(Mye - cp.Mye_0)
    # IL-6 has BOTH an acute burst (direct engagement-driven monocyte/macrophage
    # activation — the first-dose CRS spike) and the sustained myeloid-relay term
    # (TNF/IFN -> myeloid -> IL-6). The burst carries the acute peak; the relay
    # carries the lower sustained envelope, so peak >> plateau (Hosseini: first peak
    # ~order of magnitude above subsequent levels).
    # IL-6 myeloid relay scales with myeloid activation ABOVE baseline (resting Kupffer/
    # monocytes make little storm IL-6) — so the tail reflects the ENGAGEMENT-DRIVEN
    # activation difference (mosun clears -> myeloid relaxes to baseline -> low IL-6;
    # blina infusion sustains activation -> higher IL-6) rather than being dominated by
    # a large baseline-myeloid offset that smooths the two molecules together.
    dcIL6 = cy.f_IL6_burst*burst_drive + cy.p_IL6_Mye*max(Mye - cp.Mye_0, 0.0) - cy.d_IL6*cIL6
    # desensitization state: builds under engagement, recovers slowly (both directions
    # bounded so Dcyto stays in [0,1)).
    dDcyto = cy.k_cyto_desens*eng_sys*resp - cy.k_cyto_recov*Dcyto
    # Fresh-delivery / serial-engagement state Ecum in [0,1). It integrates the DRUG
    # DELIVERY RATE (the ongoing infusion), gated by CD3 occupancy so delivered drug only
    # counts when there are receptors to engage. This is the ONE signal that genuinely
    # distinguishes the two molecules: instantaneous plasma conc is nearly identical for a
    # cleared bolus (mosun) and a maintained infusion (blina) in the tail (and mosun's is
    # actually HIGHER early), so no concentration/occupancy/exposure-integral can separate
    # them — but the DELIVERY RATE is categorically different (blina: continuous nonzero;
    # mosun: an impulse at t=0, zero thereafter). Fresh drug arriving forms NEW synapses on
    # cells recovered from desensitization -> sustained cytokine production under infusion,
    # none once a bolus has cleared. Normalized by a reference infusion rate so Ecum is O(1).
    # NOTE: every window-SCORED arm is a bolus (inf_rate=0) -> deliv=0 -> Ecum=0, so the
    # sustained tail reduces to the arm-dependent CRS term and the nomination is untouched.
    deliv = (infusion / cy.inf_ref) * (C_pl / (C_pl + cy.EC50_cum))   # normalized fresh-delivery flux
    dEcum = cy.k_cum * deliv * (1.0 - Ecum) - cy.k_cum_decay * Ecum

    # ---- L4 suppression: TUMOR Treg blunts TUMOR CD8 killing (gentle Hill) ----
    # TYPED suppression term #2 (PER-CELL RATE): gS_rate scales per-Treg potency.
    # Suppression is a LOCAL tumor effect -> uses the tumor Treg pool.
    supp = 1.0/(1.0 + cp.supp_beta*gSrate*max(Treg_t/cp.Treg_t0 - 1.0, 0.0))
    supp = float(np.clip(supp, 0.05, 1.0))

    # ---- tumor: logistic regrowth minus TUMOR-CD8 killing (Treg-suppressed) ----
    # TCE killing is TRIMER-LICENSED: a T-cell engager kills only via the
    # drug-formed cytotoxic synapse (CD3-TCE-TAA bridge). Without engager there is
    # no redirected lysis, so a resting/exhausted TIL pool does NOT control the
    # tumor (why the untreated tumor grows unchecked ~5x/56d) — the drug licenses
    # killing. Durability across the window comes from the long IgG half-life
    # keeping trimer present, not from bystander killing. lic in [0,1).
    lic = eng_tum / (eng_tum + P.k_lic)
    # Effector pool: liquid targets (B-cell/BCMA engagers) are lysed in blood/marrow by the
    # CIRCULATING CD8_b pool — there is no solid-tumour infiltration barrier — whereas solid
    # tumours are killed by the tumour-infiltrating CD8_t. A real target-compartment property,
    # set per molecule via p.liquid_target (not a fitted knob).
    eff_CD8 = CD8_b if getattr(p, "liquid_target", False) else CD8_t
    kill = cp.kkill * eff_CD8 * supp * lic * help_quality
    dTumor = cp.kgrow*Tumor*(1.0 - Tumor/cp.tumor_K) - kill*Tumor/(Tumor + 0.1*cp.tumor0)

    # ---- liver tox: myeloid activity IN liver interstitium (4-1BB hepatotox site) ----
    # liver-tox substrate scales with the LIVER'S myeloid density (Kupffer cells,
    # data-grounded from Tabula Sapiens) x local liver drug exposure. This is the
    # subset-independent hepatotox site (4-1BB / costim on liver myeloid).
    liver_myeloid = PB.myeloid[PB.i_liver] * Mye * (A_liv/p.V_liv) / (1.0 + A_liv/p.V_liv)
    # SPATIAL hook (b): costim->Kupffer co-localization MODULATES (only lowers) the
    # liver-tox drive. clip to (0,1] guarantees it can never AMPLIFY / rescue the
    # ceiling (tox analogue of the tumor efficacy penalty). Default 1.0 => unchanged.
    kupffer_coloc = float(np.clip(cy.spatial_kupffer_coloc, 0.0, 1.0))
    dLiv = cy.k_liv_injury*liver_myeloid*kupffer_coloc - cy.d_liv_repair*Liv

    # ---- assemble full derivative vector ----
    dydt = np.empty(NSTATE)
    dydt[I_PL] = dA_pl
    dydt[I_IS0:I_IS0+_N] = dA_is
    dydt[I_LY] = dA_ly
    dydt[IDX["CD8_b"]] = dCD8_b;  dydt[IDX["CD8_t"]]  = dCD8_t
    dydt[IDX["CD4c_b"]]= dCD4c_b; dydt[IDX["CD4c_t"]] = dCD4c_t
    dydt[IDX["Treg_b"]]= dTreg_b; dydt[IDX["Treg_t"]] = dTreg_t
    dydt[IDX["Mye"]]  = dMye;  dydt[IDX["Tumor"]]= dTumor
    dydt[IDX["cTNF"]] = dcTNF; dydt[IDX["cIFN"]] = dcIFN; dydt[IDX["cIL2"]] = dcIL2
    dydt[IDX["cIL6"]] = dcIL6; dydt[IDX["Liv"]]  = dLiv
    dydt[IDX["Dcyto"]] = dDcyto
    dydt[IDX["Ecum"]]  = dEcum
    return dydt


def _costim_unit(gE):
    """Costim signal-2 unit: rises with effector gain, bounded. A CD3-only
    backbone (gE=1) gives 0 (no costim); a strong costim arm gives up to ~1."""
    return float(np.clip((gE - 1.0), 0.0, 2.0) / 2.0)


# ============================================================================
# Simulation driver
# ============================================================================
def initial_state(P: Params):
    cp = P.cell
    y0 = [0.0]*NSTATE
    y0[IDX["CD8_b"]] = cp.CD8_b0;  y0[IDX["CD8_t"]]  = cp.CD8_t0
    y0[IDX["CD4c_b"]]= cp.CD4c_b0; y0[IDX["CD4c_t"]] = cp.CD4c_t0
    y0[IDX["Treg_b"]]= cp.Treg_b0; y0[IDX["Treg_t"]] = cp.Treg_t0
    y0[IDX["Mye"]]  = cp.Mye_0
    y0[IDX["Tumor"]]= cp.tumor0
    return y0

def simulate(arm: ArmScores, P: Params, dose: float, t_eval=400, dose_mg=None):
    """One arm, one dose. IV bolus -> dose enters the plasma pool."""
    gains = arm.gains(P.gain)
    y0 = initial_state(P)
    y0[I_PL] = dose                      # IV bolus into plasma pool
    sol = solve_ivp(rhs, P.tspan, y0, args=(P, arm, gains), method="LSODA",
                    dense_output=True, rtol=1e-6, atol=1e-8, max_step=0.5)
    t = np.linspace(P.tspan[0], P.tspan[1], t_eval)
    Y = sol.sol(t)
    d = {k: Y[i] for k,i in IDX.items()}
    d["t"] = t
    idx = np.argmin(np.abs(t - P.t_read))
    d["frac_reduction"] = float(1.0 - Y[IDX["Tumor"]][idx]/P.cell.tumor0)
    d["peak_IL6"]  = float(np.max(Y[IDX["cIL6"]]))     # acute spike peak (for PK/PD readout)
    d["peak_TNF"]  = float(np.max(Y[IDX["cTNF"]]))
    d["peak_IFN"]  = float(np.max(Y[IDX["cIFN"]]))
    d["peak_IL2"]  = float(np.max(Y[IDX["cIL2"]]))
    # SUSTAINED IL-6 (post-acute-spike): the acute first-dose spike is engagement-
    # driven and largely arm-independent (all arms share it via eng_sys), so it must
    # NOT drive the window CRS score — that would penalise CRS-neutral arms for a
    # generic spike. The CRS-axis biology (gC, from the CD4 counter-screen) lives in
    # the SUSTAINED envelope after the spike self-limits (~day 3). This is the
    # CRS-severity driver used by the tolerability ceiling.
    post = t >= 3.0
    d["sustained_IL6"] = float(np.max(Y[IDX["cIL6"]][post])) if post.any() else d["peak_IL6"]
    # Treg tox readout = TOTAL Treg (blood+tumor); tumor Treg drives local
    # suppression, blood Treg is the systemic expansion the costim arm risks.
    Treg_tot = Y[IDX["Treg_b"]] + Y[IDX["Treg_t"]]
    Treg_tot0 = P.cell.Treg_b0 + P.cell.Treg_t0
    d["max_Treg"]  = float(np.max(Treg_tot)/Treg_tot0*P.cell.Treg_b0)  # keep prior scaling ref
    d["max_Treg_tumor"] = float(np.max(Y[IDX["Treg_t"]]))
    d["max_Treg_blood"] = float(np.max(Y[IDX["Treg_b"]]))
    d["max_CD8_tumor"]  = float(np.max(Y[IDX["CD8_t"]]))    # the KILLERS (efficacy)
    d["max_CD8_blood"]  = float(np.max(Y[IDX["CD8_b"]]))
    d["max_CD8"]   = float(np.max(Y[IDX["CD8_t"]]))         # alias: tumor CD8
    d["max_CD4c_blood"] = float(np.max(Y[IDX["CD4c_b"]]))   # CRS engine (systemic CD4)
    d["max_CD4c_tumor"] = float(np.max(Y[IDX["CD4c_t"]]))
    d["max_CD4c"]  = float(np.max(Y[IDX["CD4c_b"]] + Y[IDX["CD4c_t"]]))
    V_pl = V_PLASMA                      # physiological plasma volume = Vc (see RHS note)
    C_pl_t = Y[I_PL]/V_pl                # plasma conc trajectory (model/abstract units)
    d["Cplasma"] = C_pl_t
    d["peak_Cplasma"] = float(np.max(C_pl_t))
    d["Cplasma_read"] = float(C_pl_t[idx])
    # ---- REAL-UNIT readout ----
    # The ODE state is an amount; dividing by plasma volume gives a concentration.
    # When the caller supplies a real IV dose in mg (dose_mg), the amount state IS
    # in mg (linear PK scaling holds because the model dose enters as A_pl(0)=dose),
    # so plasma concentration is mg/L == ug/mL directly, and nM = (ug/mL)/MW(kDa)*1000.
    if dose_mg is not None and dose:
        # Linear PK: the trajectory scales with the administered amount. The model
        # dose enters as A_pl(0)=dose, so multiplying the amount state by (dose_mg/dose)
        # gives the mg trajectory; /V_pl (L) => mg/L == ug/mL.
        A_mg = Y[I_PL] * (dose_mg / dose)
        C_ugml = A_mg / V_pl
        d["Cplasma_ugml"] = C_ugml                            # ug/mL == mg/L
        d["Cplasma_nM"]   = C_ugml / P.pbpk.mw_kda * 1000.0   # nM (MW in kDa)
        d["peak_Cplasma_ugml"] = float(np.max(C_ugml))
        d["peak_Cplasma_nM"]   = float(np.max(d["Cplasma_nM"]))
    d["C_tumor_is"] = Y[I_IS0+PB.i_tumor]/PB.Vis[PB.i_tumor]   # tumor interstitial conc
    d["max_liver_injury"] = float(np.max(Y[IDX["Liv"]]))
    d["gains"] = gains
    return d


def backbone_arm():
    """CD3-only engager backbone: no costim arm => z=0 everywhere, uniform expr."""
    return ArmScores(name="CD3-only backbone", effector_z=0.0, crs_z=0.0,
                     supp_z=0.0, help_z=0.0, tier="backbone")


def _g(row, k, default=None):
    v = row.get(k, default)
    try:
        import math
        return default if v is None or (isinstance(v,float) and math.isnan(v)) else v
    except Exception:
        return v if v is not None else default


def arm_from_row(row):
    """Build an ArmScores from an enriched-matrix row (dict-like), carrying the
    kinetics (8hr/48hr) and typed-suppression (size/rate) fields, then apply the
    per-compartment expression prior. `row` keys follow qsp_input_matrix_Xv2.csv."""
    arm = ArmScores(
        name        = row["receptor"],
        effector_z  = float(_g(row,"effector_z",0.0)),
        crs_z       = _g(row,"crs_z"),
        supp_z      = _g(row,"supp_z"),
        help_z      = _g(row,"help_z"),
        sc_help_z   = _g(row,"sc_help_z"),
        effector_hit= bool(_g(row,"effector_hit",True)),
        tier        = str(_g(row,"tier","core-11")),
        supp_size_z = _g(row,"supp_size_z"),
        supp_rate_z = _g(row,"supp_rate_z"),
        crs_z_8hr   = _g(row,"crs_z_8hr"),
        crs_z_48hr  = _g(row,"crs_z_48hr"),
        supp_z_8hr  = _g(row,"supp_z_8hr"),
        supp_z_48hr = _g(row,"supp_z_48hr"),
        effector_z_sd = float(_g(row,"effector_z_sd",0.0)),
        crs_z_sd      = float(_g(row,"crs_z_sd",0.0)),
        supp_z_sd     = float(_g(row,"supp_z_sd",0.0)),
    )
    set_expression_prior(arm)
    return arm


# ============================================================================
# WINDOW METRICS — dose-response, tolerability ceilings, composite score
# ============================================================================
# The screen axes map to DISTINCT mechanisms/compartments, and the therapeutic
# window is read off the resulting dose-response — NOT a linear z-sum:
#   Axis 1 effector    -> gE   -> CD8 activation + per-CD8 kill    (efficacy)
#   Axis 2 CRS         -> gC   -> T-cell storm cytokine -> IL-6    (acute tox ceiling)
#     (onset-weighted: early 8hr rate raises the ACUTE-CRS penalty)
#   Axis 3 suppression -> gS_size (Treg compartment) + gS_rate (per-Treg output)
#                                                       (efficacy erosion, TYPED)
#   Axis 4 help        -> gH   -> CD8 persistence (durable efficacy)
#   Liver-tox axis     -> emergent: myeloid-in-liver x local drug (subset-INDEP)
# ----------------------------------------------------------------------------

def dose_response(arm: ArmScores, P: Params, doses=None):
    if doses is None:
        doses = np.logspace(-0.5, 3.0, 34)
    out = []
    for d in doses:
        r = simulate(arm, P, float(d))
        out.append(dict(dose=float(d), eff=r["frac_reduction"], peak_IL6=r["peak_IL6"],
                        sustained_IL6=r["sustained_IL6"],
                        max_Treg=r["max_Treg"], max_CD8=r["max_CD8"],
                        liver=r["max_liver_injury"], Cmax=r["peak_Cplasma"]))
    return out


def calibrate_ceiling(P: Params, tol_factor=None, doses=None):
    """CRS-severity + liver tolerability ceilings, BOTH anchored to the CD3-only
    backbone at its operating dose (the dose giving its best tolerable efficacy).
    Tolerable peak IL-6 = crs_tol_factor x backbone IL-6; tolerable liver injury
    = liver_tol_factor x backbone liver injury, both at that anchor dose.
    Returns (crs_ceiling, liver_ceiling)."""
    tol_factor = P.crs_tol_factor if tol_factor is None else tol_factor
    bb = backbone_arm()
    dr = dose_response(bb, P, doses)
    eff = np.array([r["eff"] for r in dr]); il6 = np.array([r["sustained_IL6"] for r in dr])
    liv = np.array([r["liver"] for r in dr])
    ok = np.where(eff >= P.efficacy_target)[0]
    anchor = ok[0] if len(ok) else int(np.argmax(eff))   # backbone operating dose
    crs_ceiling   = float(tol_factor * il6[anchor])
    liver_ceiling = float(P.liver_tol_factor * liv[anchor])
    return crs_ceiling, liver_ceiling


def window_metrics(arm: ArmScores, P: Params, ceiling: float, liver_ceiling=None, doses=None):
    """Extract interpretable window metrics from the dose-response.

    A dose is TOLERABLE only if it clears BOTH tolerability ceilings:
      (1) peak IL-6 <= CRS-severity ceiling      (acute CRS)
      (2) liver-injury <= liver ceiling           (hepatotox; subset-independent)
    This makes the liver-tox axis a genuine SECOND constraint on the window
    rather than an additive score term — an arm with strong (cytokine-driven or
    intrinsic) hepatotox has its tolerable dose capped by liver, which lowers
    E_tol directly. Suppression is already inside E_tol (Treg discounts killing),
    with the TYPED terms (gS_size, gS_rate) shaping the Treg trajectory.

      E_tol   : best tumor reduction at a tolerable dose (co-limited)
      MTD     : highest co-tolerable dose
      MED     : lowest dose reaching efficacy_target
      cap     : which ceiling binds first at the operating dose ('CRS'|'liver'|'none')
    """
    if doses is None:
        doses = np.logspace(-0.5, 3.0, 40)
    if liver_ceiling is None:
        liver_ceiling = P.liver_ceiling
    dr = dose_response(arm, P, doses)
    dd  = np.array([r["dose"] for r in dr]);  eff = np.array([r["eff"] for r in dr])
    il6 = np.array([r["sustained_IL6"] for r in dr]); liv = np.array([r["liver"] for r in dr])
    treg= np.array([r["max_Treg"] for r in dr])
    crs_ok = il6 <= ceiling
    liv_ok = liv <= liver_ceiling
    safe = crs_ok & liv_ok
    E_tol = float(np.max(eff[safe])) if safe.any() else float(np.min(eff))
    itol  = np.where(safe)[0]
    MTD = float(dd[itol[-1]]) if len(itol) else np.nan
    ok  = np.where(eff >= P.efficacy_target)[0]
    MED = float(dd[ok[0]]) if len(ok) else np.nan
    TI  = float(MTD/MED) if (np.isfinite(MTD) and np.isfinite(MED) and MED>0) else np.nan
    # operating dose = SMALLEST co-tolerable dose reaching >=98% of E_tol.
    # (Not the argmax-efficacy dose: when an arm has a flat efficacy plateau over a
    # broad tolerable range, the argmax jumps between distant doses as the ceiling
    # boundary shifts, making the Treg/liver tie-breaker diagnostics jitter. The
    # lowest efficacious dose is both clinically sensible — you dose at the MED,
    # not the MTD — and a STABLE operating point for the diagnostics.)
    if len(itol):
        eff_tol = eff[itol]
        thresh = 0.98 * (E_tol if E_tol > 0 else eff_tol.max())
        reach = itol[eff_tol >= thresh]
        j = int(reach[0]) if len(reach) else int(itol[np.argmax(eff_tol)])
    else:
        j = int(np.argmax(eff))
    # what caps the window: the first ceiling breached just ABOVE the operating dose
    cap = "none"
    hi = min(j+1, len(dd)-1)
    if not crs_ok[hi] and liv_ok[hi]:   cap = "CRS"
    elif crs_ok[hi] and not liv_ok[hi]: cap = "liver"
    elif not crs_ok[hi] and not liv_ok[hi]:
        cap = "CRS" if (il6[hi]/ceiling) >= (liv[hi]/liver_ceiling) else "liver"
    liver_op = float(liv[j]); treg_op = float(treg[j]); dose_op = float(dd[j])
    return dict(name=arm.name, E_tol=E_tol, MTD=MTD, MED=MED, TI=TI, cap=cap,
                liver_tol=liver_op, max_Treg_tol=treg_op, dose_op=dose_op,
                max_eff=float(np.max(eff)), gains=arm.gains(P.gain))


def composite_window_score(m, bb_ref, P: Params):
    """Therapeutic-window score vs the CD3-only backbone.

    PRIMARY term: gain in efficacy-at-a-co-tolerable-dose over backbone. Because
    the tolerable dose is co-limited by BOTH the CRS and liver ceilings, and Treg
    suppression is inside the killing term, this SINGLE term already integrates
    every liability mechanistically (unlike a linear z-sum). The Treg/liver
    contributions below are reported as INTERPRETABILITY diagnostics — they
    decompose *why* the window moved — and enter the score only with small
    weights to break ties among arms with similar E_tol, never to override it.
    """
    eff_term  = m["E_tol"] - bb_ref["E_tol"]
    # diagnostics (backbone-relative), evaluated at each arm's operating dose.
    # BOUNDED tie-breakers: a tanh caps each contribution to +/- w, so they can
    # only ORDER arms with near-equal E_tol (e.g. the two arms that both fully
    # clear the tumor) and can NEVER override a real efficacy difference. This
    # also prevents a CRS-capped arm (forced to a low dose) from being spuriously
    # credited as Treg-/liver-clean just because it is dosed low.
    treg_term = -P.w_treg * np.tanh((m["max_Treg_tol"] - bb_ref["max_Treg_tol"]) / 1.0)
    liver_term= -P.w_liver* np.tanh((m["liver_tol"]    - bb_ref["liver_tol"])    / 1.0)
    score = eff_term + treg_term + liver_term
    return float(score), dict(eff_term=float(eff_term),
                              treg_term=float(treg_term), liver_term=float(liver_term))


def score_panel(Xrows, P: Params, ceiling=None, liver_ceiling=None, effector_gate=True,
                doses=None):
    """Score every arm in the enriched matrix. Returns (rows, (crs,liver), bb_ref).
    Pass a coarser `doses` grid for fast sensitivity sweeps (ranks are robust to grid)."""
    if ceiling is None or liver_ceiling is None:
        ceiling, liver_ceiling = calibrate_ceiling(P, doses=doses)
    bb = backbone_arm()
    bb_ref = window_metrics(bb, P, ceiling, liver_ceiling, doses=doses)
    rows = []
    for row in Xrows:
        arm = arm_from_row(row)
        m = window_metrics(arm, P, ceiling, liver_ceiling, doses=doses)
        qsp, terms = composite_window_score(m, bb_ref, P)
        gated_out = effector_gate and (not arm.effector_hit)
        rows.append(dict(receptor=arm.name, tier=arm.tier,
                         effector_hit=arm.effector_hit,
                         E_tol=m["E_tol"], MTD=m["MTD"], MED=m["MED"], TI=m["TI"],
                         cap=m["cap"], dose_op=m["dose_op"],
                         max_Treg_tol=m["max_Treg_tol"], liver_tol=m["liver_tol"],
                         qsp_window=qsp, **terms,
                         gE=m["gains"]["gE"], gC=m["gains"]["gC"],
                         gS_size=m["gains"]["gS_size"], gS_rate=m["gains"]["gS_rate"],
                         gH=m["gains"]["gH"],
                         nominated=(not gated_out)))
    return rows, (ceiling, liver_ceiling), bb_ref


# ============================================================================
# CROSS-DONOR CONFIDENCE INTERVALS (mandate #3)
# ============================================================================
# Once the tox lane emits per-donor axis SDs (by_donors.h5mu / D2-4), the window
# score inherits that uncertainty. We propagate it by Monte-Carlo: resample each
# arm's axis z from Normal(z, z_sd) across donors, re-score, and report the
# window-score distribution. SCAFFOLDED NOW against zero SDs (=> point estimate);
# swap in real per-donor SDs by populating effector_z_sd/crs_z_sd/supp_z_sd on the
# input rows. This makes CIs a drop-in, not a rebuild.

def window_score_ci(row, P: Params, ceiling, liver_ceiling, bb_ref,
                    n_mc=200, seed=0):
    """Monte-Carlo cross-donor CI for one arm's window score. Resamples effector/
    CRS/suppression z from their per-donor SDs, re-scores, returns (median, lo, hi,
    samples). With zero SDs (pre-tox-donor-data) this collapses to the point value."""
    rng = np.random.default_rng(seed)
    e_sd = float(_g(row,"effector_z_sd",0.0))
    c_sd = float(_g(row,"crs_z_sd",0.0))
    s_sd = float(_g(row,"supp_z_sd",0.0))
    base = arm_from_row(row)
    if max(e_sd, c_sd, s_sd) == 0.0:
        m = window_metrics(base, P, ceiling, liver_ceiling)
        q,_ = composite_window_score(m, bb_ref, P)
        return dict(median=q, lo=q, hi=q, sd=0.0, n=0, samples=np.array([q]))
    samples = np.empty(n_mc)
    for k in range(n_mc):
        r2 = dict(row)
        r2["effector_z"] = float(_g(row,"effector_z",0.0)) + rng.normal(0, e_sd)
        if _g(row,"crs_z")  is not None: r2["crs_z"]  = _g(row,"crs_z")  + rng.normal(0, c_sd)
        if _g(row,"supp_z") is not None: r2["supp_z"] = _g(row,"supp_z") + rng.normal(0, s_sd)
        arm = arm_from_row(r2)
        m = window_metrics(arm, P, ceiling, liver_ceiling)
        q,_ = composite_window_score(m, bb_ref, P)
        samples[k] = q
    return dict(median=float(np.median(samples)),
                lo=float(np.percentile(samples, 2.5)),
                hi=float(np.percentile(samples, 97.5)),
                sd=float(np.std(samples)), n=n_mc, samples=samples)


# =============================================================================
# INTEGRATION HOOKS — spatial layer + construct-format library (interface v1)
# -----------------------------------------------------------------------------
# The spatial lane (all-organ Xenium) and the construct-format lane deliver three
# CSVs. Every hook DEFAULTS to the identity (1.0 / no-op), so the well-mixed
# window is the exact result until the CSVs land (verified byte-identical). The
# loaders below are pure, side-effect-free readers: they take a Params (and its
# frozen PBPK arrays) and return them mutated in place, tolerating missing files
# and missing rows (unlisted organs/arms keep 1.0). Contract columns:
#   spatial_exposure_perorgan.csv   : organ,          spatial_exposure   in (0,1]
#   spatial_kupffer_coloc_perarm.csv: arm/receptor,   spatial_kupffer_coloc in (0,1]
#   construct_formats.csv           : format, mw_kda, has_fc, n_CD3, n_TAA,
#                                     n_costim, reach_gate  (+ optional k_renal_max)
# =============================================================================

def apply_spatial_exposure(P: "Params", csv_path: str) -> "Params":
    """Load per-organ spatial exposure multipliers into PB.spatial_exposure.
    Organs not present in the CSV keep their default 1.0. Missing file => no-op.
    Values are clipped to (0,1] (a spatial penalty can only REDUCE exposure)."""
    import os, csv
    if not os.path.exists(csv_path):
        return P
    with open(csv_path) as fh:
        for row in csv.DictReader(fh):
            organ = (row.get("organ") or row.get("tissue") or "").strip()
            if organ in PB.names:
                v = float(row.get("spatial_exposure", 1.0))
                PB.spatial_exposure[PB.names.index(organ)] = min(max(v, 0.0), 1.0)
    return P

def apply_kupffer_coloc(P: "Params", csv_path: str, arm_name: str) -> "Params":
    """Set P.cyto.spatial_kupffer_coloc for a given costim arm from the per-arm CSV.
    MODULATE-ONLY: clipped to (0,1] so it can never rescue the liver ceiling.
    Arm not present => keeps default 1.0. Missing file => no-op."""
    import os, csv
    if not os.path.exists(csv_path):
        return P
    with open(csv_path) as fh:
        for row in csv.DictReader(fh):
            nm = (row.get("arm") or row.get("receptor") or "").strip()
            if nm == arm_name:
                v = float(row.get("spatial_kupffer_coloc", 1.0))
                P.cyto.spatial_kupffer_coloc = min(max(v, 0.0), 1.0)
                break
    return P

def apply_construct_format(P: "Params", csv_path: str, fmt_name: str) -> "Params":
    """Load one construct format's PK + valency + reach into P.pbpk. Routes has_fc
    (IgG-based => FcRn salvage; no-Fc BiTE => renal), per-format MW, valency
    (n_CD3/n_TAA/n_costim -> avidity), and the inter-arm reach_gate. Format not
    found => Params unchanged (backbone). Missing file => no-op."""
    import os, csv
    if not os.path.exists(csv_path):
        return P
    with open(csv_path) as fh:
        for row in csv.DictReader(fh):
            if (row.get("format") or "").strip() != fmt_name:
                continue
            pb = P.pbpk
            if row.get("mw_kda"):     pb.mw_kda = float(row["mw_kda"])
            if row.get("has_fc") not in (None, ""):
                pb.has_fc = str(row["has_fc"]).strip().lower() in ("1","true","yes","t")
            if row.get("k_renal_max"): pb.k_renal_max = float(row["k_renal_max"])
            if row.get("n_CD3"):      pb.n_CD3 = int(float(row["n_CD3"]))
            if row.get("n_TAA"):      pb.n_TAA = int(float(row["n_TAA"]))
            if row.get("n_costim"):   pb.n_costim = int(float(row["n_costim"]))
            if row.get("reach_gate"): pb.reach_gate = min(max(float(row["reach_gate"]), 0.0), 1.0)
            break
    return P

def spatial_window(arm: "ArmScores", P: "Params", ceiling: float,
                   liver_ceiling=None, doses=None) -> dict:
    """Spatial therapeutic window = spatially-corrected tumor-kill / spatially-
    corrected liver-tox for one costim arm. Both corrections enter mechanistically
    through the RHS hooks already wired: spatial_exposure[tumor] + reach_gate
    down-weight the tumor trimer (=> lower E_tol), spatial_kupffer_coloc lowers the
    liver-tox drive (=> lower liver_tol). With ALL hooks at their 1.0 defaults this
    returns the well-mixed window exactly (spatial_ratio == well_mixed_ratio).

    Returns E_tol (kill proxy), liver_tol (tox proxy), their ratio, and the hook
    values in force, so a caller can confirm the well-mixed default is preserved.
    """
    if liver_ceiling is None:
        liver_ceiling = P.liver_ceiling
    m = window_metrics(arm, P, ceiling, liver_ceiling, doses)
    kill = max(m["E_tol"], 0.0)                      # spatially-corrected tumor-kill
    tox  = max(m["liver_tol"], 1e-6)                 # spatially-corrected liver-tox
    return dict(name=arm.name, E_tol=m["E_tol"], liver_tol=m["liver_tol"],
                spatial_window=kill / tox, cap=m["cap"], dose_op=m["dose_op"],
                spatial_exposure_tumor=float(PB.spatial_exposure[PB.i_tumor]),
                reach_gate=float(P.pbpk.reach_gate),
                spatial_kupffer_coloc=float(P.cyto.spatial_kupffer_coloc))
