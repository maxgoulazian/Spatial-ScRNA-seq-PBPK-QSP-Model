"""kinetic_synapse.py — LITERAL kinetic multivalent-Rhoden immune synapse (per cell).

Replaces the QSS-equilibrium trimer (Schropp closed-form) with a real, time-integrated
bond ODE per T cell, integrated by an EXACT 2x2 matrix-exponential scheme (unconditionally
stable at the hour-scale PD step). The QSS limit EMERGES automatically when koff is fast
vs the PD timescale; a slow (high-affinity) arm produces a genuine lag — no QSS assumption
is imposed anywhere.

WHY (user decision 2026-07-12, "go full literal ode kinetics so we dont miss important
potentially hidden by assumptions interactions"): synapse-lifetime effects (serial killing,
kinetic proofreading/selectivity, synapse-stability cytokine) are invisible to equilibrium
magnitude but are exactly the axes a CD3/TAA/costim affinity + format sweep tunes.

STATES per T cell i (nM, on the synapse reaction-volume basis, NM_PER_COPY):
  B1_i = [drug.CD3]           armed T cell (drug singly bound to its CD3), not yet bridged
  B2_i = [CD3.drug.TAA]       bridged trimer (summed over the T cell's alive target neighbors)
Per target cell j: surv_j in [0,1] survival (kill removes TAA -> serial re-bridging emerges).

GEOMETRY (the two distances the construct sweep tunes):
  span_bridge_nm : trans CD3<->TAA arm span. (a) FEASIBILITY gate vs the emergent synapse
                   cleft (cleft relaxes toward the bound-complex size); (b) sets c_eff,trans,
                   the effective neighbour-TAA concentration the armed CD3 arm samples.
  span_cis_nm    : cis CD3<->costim arm span on the SAME T cell. Sets c_eff,cis -> costim
                   co-engagement avidity, which slows the effective CD3 detachment (koff_eff)
                   = the avidity retention a costimxCD3 arm pair buys. (Inactive for the CD3xTAA
                   validation engagers; the hook the trispecific sweep uses.)

The bridge forward rate uses the SAME surface-density x reach geometry as percell_ageff_nM
(percell_binding.py): c_eff = Ag_eff(local TAA density, arm reach). One physics, both engines.
"""
import numpy as np
import scipy.sparse as sp

AVO = 6.02214076e23
# absolute copies->local synapse concentration (nM); pinned by the validated tumor
# (CEACAM5 257,000 copies -> Rcap_TAA 6.0 nM). SAME constant as wholebody_pd.NM_PER_COPY.
NM_PER_COPY = 6.0/257000.0
R_CELL_UM   = 8.0        # target-cell radius (µm), Rhoden default
SPAN_BRIDGE_DEFAULT_NM = 12.5   # trans CD3<->TAA arm span (nm); AF3/format override per construct
SPAN_CIS_DEFAULT_NM    = 12.5   # cis CD3<->costim arm span (nm)
CLEFT_MIN_NM = 13.0      # immune-synapse cleft floor (~TCR-pMHC dimension)
CLEFT_MAX_NM = 40.0      # cleft ceiling before bond is mechanically unfavourable
# serial-killing lethal-hit rate (/day) delivered by an ENGAGED T cell to its bound target.
# Literature: a CTL delivers a lethal hit in ~tens of min once engaged; serial killing ~ few
# targets/day at saturating engagement. K_HIT sets the MAX kill rate an engaged synapse can
# deliver; the ACTUAL serial throughput is min(k_hit, cycling rate set by CD3 detachment), so
# a slow-koff (high-affinity) T cell that cannot detach is throughput-limited by koff, not k_hit.
K_HIT_DEFAULT = 12.0    # /day (~1 lethal hit per 2 h engaged) — FIXED from serial-killing literature,
                        # NOT fitted. Engine anchor is CLINICAL PD (user decision 2026-07-12): only k_death
                        # is calibrated (to the clinical TCE IL-6 + depletion curves); tumor 28.7% is an
                        # OUTPUT of the calibrated engine, not a calibration target.

def ageff_nM(dens_copies, r_cell_um=R_CELL_UM, r_arm_nm=SPAN_BRIDGE_DEFAULT_NM):
    """Geometric effective 2nd-arm concentration (nM) from a per-cell surface copy number and
    the arm reach. IDENTICAL convention to percell_binding.percell_ageff_nM (surface density x
    reach-shell), so avidity/bridge geometry is consistent across the PK and PD engines.
    Larger arm reach -> larger explored shell -> LOWER c_eff (dilution); shorter -> higher."""
    r_Ab_um = max(r_arm_nm,1e-3)/1000.0
    N=1e9
    Ag_bulk = np.asarray(dens_copies,float)*N/AVO*1e9
    SA_cell = 4.0*np.pi*r_cell_um**2
    SA_Ab   = np.pi*r_Ab_um**2
    V_Ab    = (2.0/3.0)*np.pi*r_Ab_um**3
    Am_cell = Ag_bulk/N
    return (Am_cell/SA_cell)*SA_Ab/V_Ab*1e15

def cleft_feasibility(span_bridge_nm, cleft_nm):
    """Emergent-cleft bridge feasibility in [0,1]: an arm span must reach across the apposed
    membranes' cleft. Full reach (span>=cleft) -> 1; too short -> smoothly ->0. The cleft itself
    relaxes toward the bound-complex size (see KineticSynapse: cleft = clip(span_bridge))."""
    g = span_bridge_nm/max(cleft_nm,1e-6)
    return float(np.clip((g-0.6)/0.4, 0.0, 1.0))   # <0.6x cannot bridge, >=1.0x full

def _expm2x2_apply(m11,m12,m21,m22, x1,x2, b1, dt):
    """Vectorized EXACT solution of dX/dt = M X + b over dt for a 2x2 M (per element arrays).
    X(t+dt) = expm(M dt) X + M^{-1}(expm(M dt) - I) b, with b=[b1,0]. Eigen closed-form; stable
    for the stable M here (tr<0, det>0). Falls back to series where eigenvalues coincide."""
    tr = m11+m22
    det = m11*m22 - m12*m21
    disc = np.maximum(tr*tr - 4.0*det, 0.0)
    s = np.sqrt(disc)
    l1 = 0.5*(tr+s); l2 = 0.5*(tr-s)
    e1 = np.exp(np.clip(l1*dt,-50,50)); e2 = np.exp(np.clip(l2*dt,-50,50))
    dl = l1-l2
    near = np.abs(dl) < 1e-9
    dl_safe = np.where(near, 1.0, dl)
    # expm(Mt) = a0 I + a1 M   (Sylvester for 2x2)
    a1 = np.where(near, dt*e1, (e1-e2)/dl_safe)
    a0 = np.where(near, e1 - l1*dt*e1, (l1*e2 - l2*e1)/dl_safe)
    E11 = a0 + a1*m11; E12 = a1*m12; E21 = a1*m21; E22 = a0 + a1*m22
    # homogeneous part
    hx1 = E11*x1 + E12*x2
    hx2 = E21*x1 + E22*x2
    # particular part: M^{-1}(E-I) b, b=[b1,0]
    # M^{-1} = 1/det [[m22,-m12],[-m21,m11]]
    F11 = E11-1.0; F21 = E21
    # (E-I) b = [F11*b1, F21*b1]
    g1 = F11*b1; g2 = F21*b1
    idet = np.where(np.abs(det)<1e-30, 0.0, 1.0/np.where(np.abs(det)<1e-30,1.0,det))
    px1 = idet*( m22*g1 - m12*g2)
    px2 = idet*(-m21*g1 + m11*g2)
    return hx1+px1, hx2+px2


class KineticSynapse:
    """Per-T-cell kinetic bond states over one organ's synapse graph. Consumes the SAME
    Wt_norm apportionment + synapse neighbour structure the QSS OrganPD builds, so the two
    engines are drop-in comparable. Tracks per-target survival for emergent serial killing."""
    def __init__(self, R_CD3_nM, R_TAA_nM_target, Wt_norm, syn_TAA_dens_copies,
                 kon_CD3, koff_CD3, kon_TAA, koff_TAA, kint_bridge,
                 span_bridge_nm=SPAN_BRIDGE_DEFAULT_NM, span_cis_nm=SPAN_CIS_DEFAULT_NM,
                 cis_avidity=0.0, k_hit_perday=K_HIT_DEFAULT, W_incidence=None):
        self.RC   = np.asarray(R_CD3_nM,float)                 # per-T CD3 capacity (nM), evolves with turnover
        self.RC0  = self.RC.copy()                             # CD3 set-point (KSYN = RC0*kdeg_CD3)
        self.kdeg_CD3 = 0.0                                    # free-CD3 turnover (/day); set via set_turnover()
        self.Wt_norm = Wt_norm                                 # (nT x nTarget) row-normalized apportionment
        self.dens = np.asarray(syn_TAA_dens_copies,float)      # per-T mean neighbour TAA density (copies)
        self.RTAA = np.asarray(R_TAA_nM_target,float)          # per-target TAA capacity (nM)
        # raw T x target incidence (0/1) for the TWO-SIDED conservation cap: each T cell can bridge at most
        # the ALIVE TAA reachable in its synapse (sum_j W[i,j]*RTAA[j]*surv[j]). Without this cap the bridge
        # is limited only by CD3 (RC) and the huge geometric ceff saturates every synapse -> abundance signal
        # lost (low-copy BCMA killed like high-copy CD20). WITH it, high-abundance -> CD3-limited (full), low
        # -> TAA-limited (partial). This is the linear-in-R_TAA behaviour the Schropp QSS engine had.
        self.W_inc = W_incidence                               # csr (nT x nTarget), raw 0/1
        self.kon_CD3=kon_CD3; self.koff_CD3=koff_CD3
        self.kon_TAA=kon_TAA; self.koff_TAA=koff_TAA; self.kint=kint_bridge
        # emergent cleft: relaxes toward the bound-complex span (clamped to physical window)
        self.cleft_nm = float(np.clip(span_bridge_nm, CLEFT_MIN_NM, CLEFT_MAX_NM))
        self.feas = cleft_feasibility(span_bridge_nm, self.cleft_nm)
        # trans effective TAA conc (nM) the armed CD3 samples: geometry x feasibility
        self.ceff_trans = ageff_nM(self.dens, R_CELL_UM, span_bridge_nm) * self.feas
        # cis avidity: costim co-engagement (span_cis) slows effective CD3 detachment.
        # avidity factor = c_eff,cis/(c_eff,cis+KD_costim-ish); passed in as cis_avidity in [0,1).
        self.koff_CD3_eff = koff_CD3*(1.0 - np.clip(cis_avidity,0.0,0.95))
        self.k_hit = float(k_hit_perday)       # engaged lethal-hit rate (serial-kill ceiling)
        self.dwell_engaged = None              # per-T cumulative engaged dwell-time (for cytokine)
        nT=len(self.RC)
        self.B1=np.zeros(nT); self.B2=np.zeros(nT)             # armed, bridged (nM per T cell)
        self.dwell_engaged=np.zeros(nT)                        # cumulative engaged dwell (T-days), cytokine driver
        self.surv=None                                          # per-target survival, set by owner
        self._nTarget = Wt_norm.shape[1] if hasattr(Wt_norm,'shape') else 0

    def set_turnover(self, kdeg_CD3):
        """Enable free-CD3 receptor turnover: RC relaxes to its set-point RC0 with rate kdeg_CD3 (/day),
        KSYN = RC0*kdeg_CD3. Depletion by bound drug is already captured via the (RC-B1-B2) free pool."""
        self.kdeg_CD3 = float(kdeg_CD3)

    def step(self, C_percell_T, dt_days, k_death, surv_target):
        """Advance bond states by dt (days) at per-T local drug C (nM). surv_target = per-target
        survival array (0..1); alive TAA available for bridging shrinks as targets die -> the freed
        CD3 re-bridges surviving neighbours (SERIAL KILLING emerges). Returns per-target kill hazard
        increment (length nTarget)."""
        C=np.maximum(C_percell_T,0.0)
        # ---- free-CD3 receptor turnover (operator-split, slow vs the fast binding expm below) ----
        # dRC/dt = KSYN - kdeg*RC_free ; KSYN=RC0*kdeg. Bound CD3 (B1+B2) is not degraded here (its fate is
        # koff release or kint internalization in the binding solve). Backward-Euler on the free pool for stability.
        if self.kdeg_CD3 > 0.0:
            RC_free = np.maximum(self.RC - (self.B1 + self.B2), 0.0)
            KSYN = self.RC0*self.kdeg_CD3
            RC_free_new = (RC_free + dt_days*KSYN)/(1.0 + dt_days*self.kdeg_CD3)
            self.RC = RC_free_new + (self.B1 + self.B2)        # total capacity = free + currently bound
        # available (alive) neighbour TAA fraction per T cell: Wt_norm already antigen-weighted+row-norm,
        # multiply by per-target survival and sum -> fraction of this T cell's synapse TAA still alive.
        alive_frac = np.asarray(self.Wt_norm.multiply(surv_target[None,:]).sum(1)).ravel()
        alive_frac = np.clip(alive_frac,0.0,1.0)
        Tfree = self.ceff_trans*alive_frac                     # effective ALIVE bridgeable TAA (nM)
        rate_on = self.kon_CD3*C                               # CD3 loading rate (/day; kon in /nM/day)
        kf = self.kon_TAA*Tfree                                # bridge-forming (2nd-arm on) rate (/day)
        # Two-state reduction of the real 3-species bridge (free CD3 -> B1 armed -> B2 bridged trimer):
        #   dB1/dt = rate_on*(RC-B1-B2)  - koff_CD3_eff*B1  - kf*B1        + koff_TAA*B2
        #   dB2/dt =                        kf*B1            - koff_TAA*B2  - kint*B2
        # KEY avidity physics: from the BRIDGED trimer, drug can only leave the synapse by the TAA arm
        # releasing (koff_TAA -> back to armed B1) or by internalization (kint). The CD3 arm releasing
        # does NOT dissolve the trimer (drug still held by TAA) — so the trimer lifetime is set by the
        # SLOWER (TAA) arm, i.e. bivalent avidity. This is why B2 must NOT decay at koff_CD3.
        # rate_on*(RC-B1-B2): the (RC-B1-B2) free-CD3 pool makes it a linear-in-(B1,B2) system with a
        # constant source rate_on*RC and a -rate_on on the B1,B2 diagonal/off-diagonal (receptor depletion).
        m11 = -(rate_on + self.koff_CD3_eff + kf)
        m12 =  (self.koff_TAA - rate_on)                       # B2->B1 release (koff_TAA) minus CD3-pool depletion by B2
        m21 =  kf
        m22 = -(self.koff_TAA + self.kint)                     # trimer lost ONLY via TAA-arm release + internalization (avidity)
        b1  =  rate_on*self.RC
        self.B1, self.B2 = _expm2x2_apply(m11,m12,m21,m22, self.B1,self.B2, b1, dt_days)
        # conservation clamp per T cell: B1+B2 <= RC (CD3-side)
        tot=self.B1+self.B2; over=tot>self.RC
        if over.any():
            scale=np.where(over, self.RC/np.maximum(tot,1e-30), 1.0)
            self.B1*=scale; self.B2*=scale
        self.B1=np.maximum(self.B1,0.0); self.B2=np.maximum(self.B2,0.0)
        # TAA-SIDE conservation: the bridged trimer B2 cannot exceed the ALIVE TAA reachable in the synapse.
        # taa_cap[i] = sum_j W[i,j]*RTAA[j]*surv[j] (nM). CD3-limited when TAA plentiful (high abundance),
        # TAA-limited when scarce (low abundance) -> restores the abundance grading (linear-in-R_TAA regime).
        if self.W_inc is not None:
            taa_cap = self.W_inc @ (self.RTAA*surv_target)     # per-T reachable alive TAA (nM)
            over_taa = self.B2 > taa_cap
            if over_taa.any():
                excess=np.where(over_taa, self.B2-taa_cap, 0.0)
                self.B2=np.minimum(self.B2, taa_cap)           # cap trimer at available TAA
                self.B1=self.B1+excess                         # released drug returns to armed (CD3 still held)
        # kill hazard apportioned to ALIVE targets, RE-NORMALIZED onto survivors — this is the SERIAL
        # KILLING mechanism: as a T cell's neighbours die, its bridged capacity B2 concentrates on the
        # remaining live targets (freed CD3 re-bridges survivors) rather than being wasted on dead cells.
        # Build a survivor-weighted, per-T row-renormalized apportionment each step (cheap: reuse Wt_norm
        # sparsity). W_surv[i,j] = Wt_norm[i,j]*surv_j; renormalize rows to sum 1 over LIVE neighbours.
        Ws = self.Wt_norm.multiply(surv_target[None,:]).tocsr()
        rs = np.asarray(Ws.sum(1)).ravel(); rs[rs<1e-12]=1.0
        Ws = sp.diags(1.0/rs) @ Ws                              # rows renormalized onto survivors
        has_live = (np.asarray(Ws.sum(1)).ravel() > 1e-9).astype(float)   # 1 if any live neighbour
        # ---- EXPLICIT engage -> lethal-hit -> detach -> re-engage CYCLE (serial killing) --------------
        # p_eng = fraction of this T cell's CD3 committed to a bridged (engaged) synapse = B2/RC in [0,1].
        p_eng = np.clip(self.B2/np.maximum(self.RC,1e-30), 0.0, 1.0)
        # SERIAL THROUGHPUT via the engage/hit/detach RACE — the mechanism that makes CD3 affinity a
        # WINDOW optimum, not a monotone knob. Once bridged, two clocks race:
        #   k_hit         : deliver the lethal hit (needs sustained engagement)
        #   koff_CD3_eff  : detach and move on
        # P(productive hit before detach) = k_hit/(k_hit+koff_CD3_eff).  A T cell must then DETACH to
        # re-engage the next target; cycling rate ~ koff_CD3_eff. Productive serial rate (targets/day):
        #   serial_rate = koff_CD3_eff * [k_hit/(k_hit+koff_CD3_eff)] * p_eng
        #             = k_hit*koff_CD3_eff/(k_hit+koff_CD3_eff) * p_eng   (harmonic — a Michaelis-like hump)
        # LIMITS: slow koff -> can't cycle (rate -> koff, low); fast koff -> detaches before hitting
        # (rate -> k_hit, but hit prob ->0). The harmonic k_hit*koff/(k_hit+koff) rises with koff and
        # SATURATES toward k_hit for koff >> k_hit. EMPIRICAL sweep (this build, k_hit=12/day): endpoint
        # targets-killed peaks at koff ~ 346/day (KD_CD3 40 nM x kon) — i.e. the clinical CD3 affinity sits
        # on the high-koff plateau, NOT at koff~k_hit. The WINDOW metric (kill per engaged-dwell) improves
        # monotonically with koff across the tested range (detuning-widens-window on the CRS axis); the
        # efficacy endpoint is a broad plateau, not a sharp koff~k_hit optimum. Both emergent from the cycle.
        koff_e = self.koff_CD3_eff
        serial_rate = (self.k_hit*koff_e/(self.k_hit+koff_e+1e-30)) * p_eng * has_live   # targets/day per T
        # accumulate engaged dwell (T-days) — the SYNAPSE-STABILITY cytokine driver (slow koff -> longer
        # dwell per kill -> more cytokine per target killed = the over-stable-synapse CRS mechanism)
        self.dwell_engaged += dt_days*p_eng
        # apportion the serial kill onto surviving neighbours (renormalized), scale by k_death intrinsic pot.
        dkill_tgt = Ws.T @ (serial_rate)                       # length nTarget
        return dt_days*k_death*dkill_tgt

    def engaged_dwell_rate(self):
        """Instantaneous engaged-synapse level (sum p_eng) — the cytokine emission proxy. Over-stable
        synapses (slow koff, high p_eng sustained) emit more; fast-cycling arms emit less per kill."""
        return float(np.clip(self.B2/np.maximum(self.RC,1e-30),0,1).sum())

    def bridged_total(self):
        return float(self.B2.sum())
