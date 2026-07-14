"""kinetic_rhoden_percell.py - LITERAL kinetic bivalent Rhoden binding, per cell, spatially resolved.

6-species kinetic bivalent-binding ODE. The bivalent CROSSLINK structure (BAg1/BAg2/Bdbl with the
geometric effective 2nd-arm concentration AgEFF) is Rhoden et al. 2016 (kinetic form: bioRxiv
10.1101/2022.09.12.507653). Rhoden's own model assumes antigen kinetics NEGLIGIBLE (no turnover);
the receptor-turnover terms here (KSYN=Ag0*kDEG synthesis, -kDEG*Ag free-antigen degradation, kTMD
internalization of bound species) are ADDED on top -- standard TMDD receptor turnover, taken VERBATIM
from the user-supplied MATLAB scheme. Generalized to
(a) the multi-arm engager format space (CD3/costim/TAA, each arm valency n in {0,1,2}, total
<=4 = tetravalent max, three geometric spans) and (b) spatial per-cell resolution (the cross-cell
bridge couples a T cell to its tumour-cell neighbours via the tissue graph; per-cell contraction).

ONE binding solve, used IDENTICALLY by PK and PD:
  PK : TMDD sink = kTMD*(every bound drug species) per cell, summed -> plasma clearance.
  PD : kill ~ cross-cell bridge (CD3.drug.TAA); CRS ~ engaged dwell; costim ~ costim bridge.

Receptor turnover EXPLICIT (the piece the QSS sink lacked -> broke terminal TMDD):
  dAg/dt has +KSYN (=Ag0*kDEG set-point) and -kDEG*Ag; every bound species internalizes at kTMD.
  At C=0 steady state Ag -> Ag0 exactly. No QSS assumption anywhere.

REFERENCE SCHEME (same-cell, two co-engageable antigens Ag1/Ag2):
  dAg1  = KSYN1 - kon1*C*Ag1 + koff1*BAg1 - kon1*Ag1EFF*(Ag1/Ag1_0)*BAg2 + koff1*Bdbl - kDEG*Ag1
  dAg2  = KSYN2 - kon2*C*Ag2 + koff2*BAg2 - kon2*Ag2EFF*(Ag2/Ag2_0)*BAg1 + koff2*Bdbl - kDEG*Ag2
  dBAg1 = kon1*C*Ag1 - koff1*BAg1 - kon2*Ag2EFF*(Ag2/Ag2_0)*BAg1 + koff2*Bdbl - kTMD*BAg1
  dBAg2 = kon2*C*Ag2 - koff2*BAg2 - kon1*Ag1EFF*(Ag1/Ag1_0)*BAg2 + koff1*Bdbl - kTMD*BAg2
  dBdbl = kon1*BAg2*Ag1EFF*(Ag1/Ag1_0) + kon2*BAg1*Ag2EFF*(Ag2/Ag2_0) - (koff1+koff2)*Bdbl - 2*kTMD*Bdbl
  KSYNk = Agk_0*kDEG ;  AgkEFF = Rhoden geometric effective 2nd-arm conc (surface dens * IgG area / vol)

VARIANT = identity of the two arms + which span (all the SAME equations):
  cis T-side     : Ag1=CD3, Ag2=costim (or 2xCD3)            span=span_coeng_T
  cis tumour-side: Ag1=Ag2=TAA (or TAA1,TAA2 dual)           span=span_coeng_tumor
  trans bridge   : Ag1=CD3 (T cell), Ag2=TAA (tumour nbr)    span=span_bridge (Ag2EFF=graph-reachable TAA)
"""
import numpy as np
AVO = 6.02214076e23

def geo_ageff_nM(rec_pc, r_cell_um=8.0, span_nm=12.5):
    """Rhoden geometric effective 2nd-arm conc (nM). IDENTICAL to multiarm_binding.geo_ageff_nM."""
    rec_pc = np.maximum(np.asarray(rec_pc, float), 0.0)
    SA_cell = 4.0*np.pi*r_cell_um**2
    r_Ab_um = max(float(span_nm), 1e-3)*1e-3
    SA_Ab = np.pi*r_Ab_um**2; V_Ab = (2.0/3.0)*np.pi*r_Ab_um**3
    Ag_bulk = rec_pc*1e9/AVO*1e9
    Am_cell = Ag_bulk/1e9; Am_SA = Am_cell/SA_cell; Ag_SA = Am_SA*SA_Ab
    return Ag_SA/V_Ab*1e15

def rhoden_bivalent_step(C_free_nM, Ag1, Ag2, BAg1, BAg2, Bdbl,
                         Ag1_0, Ag2_0, kon1, koff1, kon2, koff2,
                         Ag1EFF, Ag2EFF, kDEG, kTMD, dt, nsub=None):
    """Advance the 6-species kinetic bivalent binding one step dt (days), VERBATIM the scheme,
    vectorized over per-cell arrays. Rates /day, conc nM. Sub-cycles for stability.
    Returns (Ag1,Ag2,BAg1,BAg2,Bdbl, intern_flux) ; intern_flux=kTMD*(BAg1+BAg2+2*Bdbl) step-avg (nM/day)."""
    C = np.maximum(np.asarray(C_free_nM, float), 0.0)
    Ag1 = np.array(Ag1, float, copy=True); Ag2 = np.array(Ag2, float, copy=True)
    BAg1 = np.array(BAg1, float, copy=True); BAg2 = np.array(BAg2, float, copy=True)
    Bdbl = np.array(Bdbl, float, copy=True)
    KSYN1 = Ag1_0*kDEG; KSYN2 = Ag2_0*kDEG
    a1r = np.where(np.asarray(Ag1_0) > 1e-30, 1.0/np.maximum(Ag1_0, 1e-30), 0.0)
    a2r = np.where(np.asarray(Ag2_0) > 1e-30, 1.0/np.maximum(Ag2_0, 1e-30), 0.0)
    if nsub is None:
        # ALL pseudo-first-order rates incl. the AVIDITY CROSSLINK (kon*AgEFF ~ hundreds/day for real
        # densities) — omitting it under-substeps -> explicit-Euler overshoot -> clamp injects mass.
        Cm=np.maximum(C,0.0)
        rmax = float(np.max(kon1*Cm+koff1+kon2*Cm+koff2+kon1*np.maximum(Ag1EFF,0.0)
                            +kon2*np.maximum(Ag2EFF,0.0)+kDEG+2.0*kTMD)) if np.size(C) else 0.0
        nsub = int(min(max(1, np.ceil(rmax*dt/0.1)), 20000))
    h = dt/max(nsub, 1)
    intern_accum = np.zeros_like(Ag1*1.0)
    for _ in range(int(nsub)):
        x1 = Ag1EFF*(Ag1*a1r); x2 = Ag2EFF*(Ag2*a2r)
        dAg1 = KSYN1 - kon1*C*Ag1 + koff1*BAg1 - kon1*x1*BAg2 + koff1*Bdbl - kDEG*Ag1
        dAg2 = KSYN2 - kon2*C*Ag2 + koff2*BAg2 - kon2*x2*BAg1 + koff2*Bdbl - kDEG*Ag2
        dB1  = kon1*C*Ag1 - koff1*BAg1 - kon2*x2*BAg1 + koff2*Bdbl - kTMD*BAg1
        dB2  = kon2*C*Ag2 - koff2*BAg2 - kon1*x1*BAg2 + koff1*Bdbl - kTMD*BAg2
        dBd  = kon1*x1*BAg2 + kon2*x2*BAg1 - (koff1+koff2)*Bdbl - 2.0*kTMD*Bdbl
        Ag1 = np.maximum(Ag1+h*dAg1, 0.0); Ag2 = np.maximum(Ag2+h*dAg2, 0.0)
        BAg1 = np.maximum(BAg1+h*dB1, 0.0); BAg2 = np.maximum(BAg2+h*dB2, 0.0)
        Bdbl = np.maximum(Bdbl+h*dBd, 0.0)
        intern_accum += (kTMD*(BAg1+BAg2+2.0*Bdbl))*h
    return Ag1, Ag2, BAg1, BAg2, Bdbl, intern_accum/max(dt,1e-30)


def rhoden_samecell_bivalent_step(C_free_nM, Ag, BAg1, Bdbl, Ag_0,
                                  kon, koff, AgEFF, kDEG, kTMD, dt, nsub=None):
    """SAME-antigen bivalent avidity (both arms bind the SAME antigen; e.g. glofitamab 2xCD20,
    alnuctamab 2xBCMA). Reduced Rhoden 3-state single-antigen system — the correct form when Ag1==Ag2,
    replacing the (Ag1,Ag2) two-pool call that leaves Bdbl inert when Ag2=0.

      states (nM, synapse basis): Ag (free), BAg1 (one arm bound), Bdbl (both arms -> 2 neighbour copies)
      dAg   = KSYN - kon*C*Ag + koff*BAg1  - kon*x*BAg1 + 2*koff*Bdbl - kDEG*Ag     (x=AgEFF*Ag/Ag0)
      dBAg1 = kon*C*Ag - koff*BAg1         - kon*x*BAg1 + 2*koff*Bdbl - kTMD*BAg1
      dBdbl = kon*x*BAg1                   - 2*koff*Bdbl              - 2*kTMD*Bdbl
      KSYN = Ag_0*kDEG.   crosslink feeds from BAg1 (second arm grabs a neighbour copy of the SAME Ag).
      census Ag_tot = Ag + BAg1 + 2*Bdbl conserved (turnover/internalisation aside).

    BACKWARD-EULER, receptors-as-states -> unconditionally stable + census exactly preserved (no
    manufacture) even at large AgEFF. nsub set by SLOW scales (kDEG/kTMD), not the stiff crosslink.
    Returns (Ag, BAg1, Bdbl, intern_flux) ; intern_flux = kTMD*(BAg1 + 2*Bdbl) step-avg (nM/day).
    """
    C = np.maximum(np.asarray(C_free_nM, float), 0.0)
    Ag = np.array(Ag, float, copy=True); BAg1 = np.array(BAg1, float, copy=True)
    Bdbl = np.array(Bdbl, float, copy=True)
    Ag_0 = np.asarray(Ag_0, float)
    KSYN = Ag_0*kDEG
    a_r = np.where(Ag_0 > 1e-30, 1.0/np.maximum(Ag_0, 1e-30), 0.0)
    AgEFF = np.asarray(AgEFF, float)
    if nsub is None:
        nsub = int(min(max(1, np.ceil((kDEG + 2.0*kTMD)*dt/0.25)), 16))   # SLOW scales only (BE is stable)
    h = dt/max(nsub, 1)
    n = np.size(Ag) if np.ndim(Ag) else 1
    Ag = np.broadcast_to(Ag, (n,)).astype(float).copy() if np.ndim(Ag) else np.full(n, float(Ag))
    intern = np.zeros(n)
    Cb = np.broadcast_to(C, (n,)).astype(float) if np.ndim(C) else np.full(n, float(C))
    KSYNb = np.broadcast_to(KSYN, (n,)).astype(float) if np.ndim(KSYN) else np.full(n, float(KSYN))
    for _ in range(int(nsub)):
        x = AgEFF*(Ag*a_r)                         # effective 2nd-arm conc, frozen within the BE solve
        kb = kon*Cb                                # first-arm on-rate (per cell)
        kx = kon*x                                 # crosslink on-rate (per cell)
        # linear generator M over y=[Ag, BAg1, Bdbl] with receptors as states (census left-invariant)
        # Ag'   = KSYN - kb*Ag + koff*BAg1 - kx*BAg1 + 2koff*Bdbl - kDEG*Ag
        # BAg1' = kb*Ag - koff*BAg1 - kx*BAg1 + 2koff*Bdbl - kTMD*BAg1
        # Bdbl' = kx*BAg1 - 2koff*Bdbl - 2kTMD*Bdbl
        M00 = -(kb + kDEG); M01 = (koff - kx); M02 = 2.0*koff
        M10 = kb;           M11 = -(koff + kx + kTMD); M12 = 2.0*koff
        M20 = np.zeros(n);  M21 = kx;          M22 = -(2.0*koff + 2.0*kTMD)
        M = np.empty((n,3,3))
        M[:,0,0]=M00; M[:,0,1]=M01; M[:,0,2]=M02
        M[:,1,0]=M10; M[:,1,1]=M11; M[:,1,2]=M12
        M[:,2,0]=M20; M[:,2,1]=M21; M[:,2,2]=M22
        s = np.zeros((n,3)); s[:,0]=KSYNb
        y0 = np.stack([Ag, BAg1, Bdbl], 1)
        A = (np.broadcast_to(np.eye(3),(n,3,3)) - h*M).copy()
        rhs = (y0 + h*s)
        # ---- SINGULAR-ROW GUARD ---------------------------------------------------------------------------
        # A = I - h*M has diagonal 1 + h*(non-negative rates) >= 1, so for FINITE inputs it cannot be singular.
        # A singular row therefore means a NON-FINITE rate leaked in (NaN/Inf drug conc, or Ag_0=0 giving a
        # degenerate a_r). Observed 2026-07-13: elranatamab at TSIM=24 crashed the whole run with
        # LinAlgError -- one bad cell out of ~1e5 killed a 40-minute simulation.
        # A cell with a degenerate matrix has no well-posed update, so FREEZE it (y = y0) rather than crash or
        # silently propagate NaN. Freezing is the physically-correct no-op for a fully depleted cell, but it is
        # COUNTED and WARNED so it can never hide a real bug.
        bad = ~(np.isfinite(A).all(axis=(1,2)) & np.isfinite(rhs).all(axis=1))
        if bad.any():
            det = np.full(n, 1.0)
            det[~bad] = np.linalg.det(A[~bad])
        else:
            det = np.linalg.det(A)
        sing = bad | (np.abs(det) < 1e-12)
        if sing.any():
            nbad = int(sing.sum())
            if not getattr(rhoden_samecell_bivalent_step, "_warned", False):
                import sys as _s
                _s.stderr.write(f"[rhoden-BE] {nbad}/{n} cells had a SINGULAR/non-finite matrix -> FROZEN this "
                                f"substep (no update). Non-finite inputs: {int(bad.sum())}. This is a guard, not "
                                f"a fix for the upstream cause -- if nbad is large, something is diverging.\n")
                rhoden_samecell_bivalent_step._warned = True
            A[sing] = np.eye(3)
            rhs[sing] = np.nan_to_num(y0[sing], nan=0.0, posinf=0.0, neginf=0.0)
        y = np.linalg.solve(A, rhs[:,:,None])[:,:,0]
        Ag = np.maximum(y[:,0],0.0); BAg1 = np.maximum(y[:,1],0.0); Bdbl = np.maximum(y[:,2],0.0)
        intern += (kTMD*(BAg1 + 2.0*Bdbl))*h
    return Ag, BAg1, Bdbl, intern/max(dt,1e-30)
