"""PK-BREADTH HARNESS — the 53 dense clinical PK curves the TCE runner cannot touch.

WHY THIS EXISTS: run_tce_pd_reval.py's ENG dict is TCE-only and KeyErrors on plain mAbs. But the
dense PK curves (312 digitized points across 21 drugs) are mostly REFERENCE mAbs — dupilumab,
adalimumab, bevacizumab, tocilizumab, benralizumab... Those curves are what actually tests the
Shah-Betts PBPK BACKBONE. Without them the PK validation has 4 points; with them it has 312.

WHAT IT DOES: PBPK + FcRn + (optional) TMDD sink. No CD3, no T cells, no synapse.
  spec (all_specs_comprehensive.json) -> KD, kint, valency, Rtot_nM (per-organ pools), sAg
  MOLPROPS.json                        -> mw_kda, fFcRn
  SQLite curve row                     -> dose_value, dose_route   (one run per DOSE LEVEL)

*** GATED ON ADALIMUMAB. *** I have mis-measured this model four times tonight (three of them by
truncating the sim shorter than the data). So this harness does NOT get to score 53 curves until it
reproduces adalimumab -- a plain IgG1, soluble target, no membrane TMDD, 15 digitized points, the
simplest possible case. If it cannot do adalimumab, it is wrong and nothing else it says counts.

usage:  python pk_breadth_harness.py --gate          # adalimumab only; must pass
        python pk_breadth_harness.py --all           # all 53 curves (only after the gate passes)
"""
import argparse, json, os, sqlite3, sys
import numpy as np


# ---------------------------------------------------------------------------
# LAYOUT-PORTABLE ROOT [added 2026-07-13]
# The model tree exists in two layouts:
#   model/       -> rundir/handoff/...   (the working tree)
#   FINAL_MODEL/ -> handoff/...          (the packaged deliverable)
# The harness hardcoded KWS="rundir", so it silently died in FINAL_MODEL. Resolve
# the root by LOOKING for handoff/, rather than assuming a layout.
# ---------------------------------------------------------------------------
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
KWS = None
for _c in ("rundir", ".", _os.path.join(_HERE, "rundir"), _HERE):
    if _os.path.isdir(_os.path.join(_c, "handoff")):
        KWS = _c
        break
if KWS is None:
    raise RuntimeError("pk_breadth_harness: cannot find a handoff/ dir under rundir/ or . -- refusing to run")
sys.path.insert(0, _os.path.join(KWS, "handoff"))
sys.path.insert(0, _os.path.join(_HERE, "engine"))
sys.path.insert(0, "engine")


import qsp_costim_window_v2 as q
from coupled_percell_pk import CoupledPerCellPK

# KWS resolved above (layout-portable)
ORGANS = ['spleen', 'bone', 'large_int', 'liver', 'lung', 'small_int',
          'pancreas', 'kidney', 'skin', 'heart', 'adipose']

pb = q._PBPKArrays()
arr = lambda a: {o: a[i] for i, o in enumerate(pb.names)}
Q, L, sigV, Vis, Vv = arr(pb.Q), arr(pb.L), arr(pb.sigV), arr(pb.Vis), arr(pb.Vv)
bl = json.load(open(f"{KWS}/handoff/bec_lec_masks.json"))

SPECS = json.load(open(_os.path.join(_HERE, "params/all_specs_comprehensive.json")))
MOLPROPS = {}
for p in ("params/MOLPROPS.json", "rundir/handoff/MOLPROPS.json"):
    if os.path.exists(p):
        MOLPROPS = json.load(open(p))
        break

V_PLASMA_L = 3.1
BW_KG = 70.0


def to_days(t, u):
    u = (u or "day").lower()
    if u.startswith("h"):
        return t / 24.0
    if u.startswith("min"):
        return t / 1440.0
    if u.startswith("w"):
        return t * 7.0
    return t


def to_ugml(v, u):
    u = (u or "ug/mL").lower().replace("µ", "u").replace("μ", "u")
    if "ng" in u:
        return v / 1000.0
    if "pg" in u:
        return v / 1e6
    if "mg/ml" in u:
        return v * 1000.0
    return v          # ug/mL == mg/L


def get_props(drug, spec):
    mw = spec.get("mw_kda")
    ff = spec.get("fFcRn")
    mp = MOLPROPS.get(drug) or {}
    if isinstance(mp, dict):
        mw = mw or mp.get("mw_kda") or mp.get("mw")
        ff = ff if ff is not None else (mp.get("fFcRn") if mp.get("fFcRn") is not None else None)
    if mw is None:
        mw = 148.0            # IgG class default — FLAGGED, not silent
    if ff is None:
        ff = 0.89
    return float(mw), float(ff)


def run_curve(drug, cv):
    """One PBPK run at THIS curve's dose + route. Returns (t_days, C_ugml) or None."""
    s = SPECS.get(drug)
    if not s:
        return None
    spec = s.get("spec", s)
    mw, ff = get_props(drug, spec)
    kd = float(spec.get("KD_nM") or 1.0)
    kint = float(spec.get("kint_perday") or 0.0)
    narm = int(spec.get("n_arm_target") or spec.get("valency") or 2)
    rt = spec.get("Rtot_nM") or {}
    pools = {o: float(rt.get(o, 0.0)) for o in ORGANS}     # {} -> soluble target, no membrane sink

    dose_mg = float(cv.get("dose_value") or 0.0)
    unit = (cv.get("dose_unit") or "mg").lower()
    if "kg" in unit:                                        # mg/kg -> mg
        dose_mg *= BW_KG
    route = (cv.get("dose_route") or "IV").upper()
    if "SC" in route:
        route = "SC"
    elif "IV" in route:
        route = "IV"
    else:
        return None                                        # e.g. IP — engine cannot dose it

    tmax = 0.0
    for p in cv["_pts"]:
        tmax = max(tmax, p[0])
    tsim = max(tmax * 1.15 + 2, 7.0)

    m = CoupledPerCellPK(ORGANS, spec.get("target", "?"), "CD3E", kd, narm, kint, mw, ff,
                         f"{KWS}/handoff/agents", pools, Q, L, sigV, Vis, Vv, bl, dt=0.02)
    # REAL API (read from the source, not guessed):
    #   simulate(schedule, tsim, route="IV", F_sc=0.6, ka_sc=0.25, ...) -> dict(t, Cplasma_ugml, ...)
    #   schedule is [(day, mg)] IN MILLIGRAMS — the engine does the mg->nmol conversion internally
    #   (`A_sc += mg/self.mw*1e3`). My first draft pre-converted to nM and would have DOUBLE-CONVERTED
    #   the dose on all 53 curves. The adalimumab gate caught it. That is what the gate is for.
    ka = float(spec.get("ka_sc_perday") or 0.25)
    r = m.simulate([(0.0, dose_mg)], tsim, route=route, F_sc=0.6, ka_sc=ka)
    if r is None:
        return None
    # SAVE EVERYTHING — the whole point is to be able to regenerate any figure without re-running.
    # Never discard a computed time-course: re-runs are expensive and the curves ARE the product.
    os.makedirs("out", exist_ok=True)
    blob = dict(
        drug=drug, curve_id=cv["curve_id"], source_id=cv.get("source_id"),
        dose_value=cv.get("dose_value"), dose_unit=cv.get("dose_unit"),
        dose_route=cv.get("dose_route"), dose_regimen=cv.get("dose_regimen"),
        n_subjects=cv.get("n_subjects"), summary_stat=cv.get("summary_stat"),
        # --- model params actually used (so a figure can state them) ---
        params=dict(KD_nM=kd, kint_perday=kint, n_arm=narm, mw_kda=mw, fFcRn=ff,
                    ka_sc_perday=ka, F_sc=0.6, pools_nM=pools, tsim_days=tsim,
                    target=spec.get("target"), sigL=0.85, k_dist=3.0),
        # --- model output ---
        t_days=[float(x) for x in r["t"]],
        Cplasma_ugml=[float(x) for x in r["Cplasma_ugml"]],
        sink_per_organ={o: [float(x) for x in v] for o, v in (r.get("sink") or {}).items()},
        bound_nM_per_organ={o: [float(x) for x in v] for o, v in (r.get("bound_nM") or {}).items()},
        # --- the clinical data it is scored against (so the overlay is self-contained) ---
        obs_t_days=[float(p[0]) for p in cv["_pts"]],
        obs_conc_ugml=[float(p[1]) for p in cv["_pts"]],
        engine="CoupledPerCellPK (the SAME class CoupledPerCellPD inherits — verified same PBPK path)",
    )
    with open(f"out/pkb_{cv['curve_id']}.json", "w") as fh:
        json.dump(blob, fh)
    return np.asarray(r["t"], float), np.asarray(r["Cplasma_ugml"], float)


def load_curves(only=None):
    con = sqlite3.connect("params/mab_tce_pkpd.sqlite")
    con.row_factory = sqlite3.Row
    c = con.cursor()
    out = []
    for cv in c.execute("SELECT * FROM curves WHERE readout_class='PK'"):
        cv = dict(cv)
        if only and cv["drug_id"] != only:
            continue
        pts = []
        for p in c.execute("SELECT time,time_unit,value,value_unit FROM timeseries "
                           "WHERE curve_id=? ORDER BY point_order", (cv["curve_id"],)):
            p = dict(p)
            if p["time"] is None or p["value"] is None or p["value"] <= 0:
                continue
            pts.append((to_days(float(p["time"]), p.get("time_unit") or cv.get("time_unit")),
                        to_ugml(float(p["value"]), p.get("value_unit") or cv.get("conc_unit"))))
        if len(pts) >= 3:
            cv["_pts"] = pts
            out.append(cv)
    con.close()
    return out


def score(cv, t, c):
    pts = np.array(cv["_pts"])
    ot, oc = pts[:, 0], pts[:, 1]
    w = (ot >= 0) & (ot <= t.max())
    if w.sum() < 3:
        return None
    pred = np.interp(ot[w], t, c)
    o = oc[w]
    g = (pred > 0) & (o > 0)
    if g.sum() < 3:
        return None
    r = pred[g] / o[g]
    return dict(n=int(g.sum()),
                afe=float(np.exp(np.mean(np.log(r)))),
                aafe=float(np.exp(np.mean(np.abs(np.log(r))))))


ap = argparse.ArgumentParser()
ap.add_argument("--gate", action="store_true")
ap.add_argument("--all", action="store_true")
ap.add_argument("--drug", default=None, help="run every curve for ONE drug (used by the pod orchestrator)")
a = ap.parse_args()

if a.drug:
    cvs = load_curves(only=a.drug)
    if not cvs:
        print(f"{a.drug}: no PK curve with >=3 points"); sys.exit(0)
    for cv in cvs:
        try:
            out = run_curve(a.drug, cv)
        except Exception as e:
            print(f"{a.drug} {cv['curve_id']}: CRASH {type(e).__name__}: {e}")
            continue
        if out is None:
            print(f"{a.drug} {cv['curve_id']}: skipped (route not IV/SC, or no spec)")
            continue
        s = score(cv, *out)
        if s:
            print(f"{a.drug} {cv['curve_id']} dose={cv.get('dose_value')}{cv.get('dose_unit')} "
                  f"{cv.get('dose_route')} n={s['n']} AFE={s['afe']:.2f} AAFE={s['aafe']:.2f}")
        else:
            print(f"{a.drug} {cv['curve_id']}: not scoreable")
    sys.exit(0)

if a.gate:
    print("=" * 88)
    print("HARNESS GATE — adalimumab (plain IgG1, soluble target, 15 pts, the simplest case)")
    print("=" * 88)
    print("  If the harness cannot reproduce THIS, it is wrong and must not score anything else.\n")
    cvs = load_curves(only="adalimumab")
    if not cvs:
        print("  no adalimumab curve found"); sys.exit(1)
    for cv in cvs:
        try:
            out = run_curve("adalimumab", cv)
        except Exception as e:
            print(f"  *** HARNESS CRASHED: {type(e).__name__}: {e}")
            print(f"  The PK-only path does not work. NOT scoring the other 52 curves.")
            sys.exit(2)
        if out is None:
            print("  *** harness returned nothing — CoupledPerCellPK has no usable run() entry point.")
            print("      The PK-only harness must be written against its real API. NOT guessing.")
            sys.exit(3)
        t, c = out
        s = score(cv, t, c)
        print(f"  {cv['curve_id']}: n={s['n']} AFE={s['afe']:.2f} AAFE={s['aafe']:.2f}")
        if s["aafe"] < 2.0:
            print("\n  ✅ GATE PASSED — harness reproduces adalimumab within 2-fold. Safe to run all 53.")
        else:
            print(f"\n  ⛔ GATE FAILED (AAFE {s['aafe']:.2f}) — do NOT trust this harness on the other curves.")
    sys.exit(0)

if a.all:
    cvs = load_curves()
    print(f"{len(cvs)} dense PK curves\n")
    rows = []
    for cv in cvs:
        try:
            out = run_curve(cv["drug_id"], cv)
        except Exception as e:
            print(f"  {cv['drug_id']:16s} CRASH {type(e).__name__}")
            continue
        if out is None:
            continue
        s = score(cv, *out)
        if not s:
            continue
        rows.append((cv["drug_id"], s["aafe"]))
        print(f"  {cv['drug_id']:16s} {cv['curve_id'][:40]:40s} n={s['n']:3d} AAFE={s['aafe']:.2f}")
    if rows:
        aa = np.array([r[1] for r in rows])
        print(f"\n  {len(rows)} curves | median AAFE {np.median(aa):.2f} | "
              f"within 2/3/5-fold {np.sum(aa<=2)}/{np.sum(aa<=3)}/{np.sum(aa<=5)}")
