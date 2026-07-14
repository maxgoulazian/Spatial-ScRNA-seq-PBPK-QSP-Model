#!/usr/bin/env python
"""
FINAL_MODEL entry point — run the whole-body per-cell PBPK-PD/QSP engine for one T-cell engager.

Usage:
    python run.py <engager_name> [--tsim DAYS] [--dt DT]

Examples:
    python run.py teclistamab
    python run.py mosunetuzumab --tsim 24
    python run.py glofitamab --tsim 24 --dt 0.04

Output:
    handoff/tce_pd_<engager>.json   (PK + PD + IL-6 + per-organ kill, one file per engager)

The engine runs PK and PD in a single coupled pass: each molecule's own plasma exposure drives its own
per-cell pharmacology on the spatial grid. All binding is true kon/koff Rhoden kinetics at single-cell
resolution. See MODEL_DOCUMENTATION/ for the full method and OPERATION.md for the parameter reference.
"""
import os, sys, argparse

HERE = os.path.dirname(os.path.abspath(__file__))

def main():
    ap = argparse.ArgumentParser(description="Run the per-cell PBPK-PD/QSP engine for one engager.")
    ap.add_argument("engager", help="engager name, e.g. teclistamab, mosunetuzumab, glofitamab")
    ap.add_argument("--tsim", type=float, default=None, help="simulated days (default: engine default)")
    ap.add_argument("--dt", type=float, default=None, help="timestep in days (default 0.01; 0.04 is ~4x faster, stable)")
    ap.add_argument("--graph-cache", default="0", help="WB_GRAPH_CACHE dir, or 0/off to rebuild each run (default 0)")
    args = ap.parse_args()

    # --- path + environment setup (the launcher's whole job) ---
    engine = os.path.join(HERE, "engine")
    handoff = os.path.join(HERE, "handoff")
    for p in (engine, handoff, HERE):
        if p not in sys.path:
            sys.path.insert(0, p)
    # local GPU driver is too old for CUDA init -> force CPU (harmless if no GPU).
    os.environ.setdefault("CUDA_VISIBLE_DEVICES", "")
    os.environ["WB_GRAPH_CACHE"] = args.graph_cache
    if args.tsim is not None: os.environ["TSIM"] = str(args.tsim)
    if args.dt is not None:   os.environ["PD_DT"] = str(args.dt)

    # the harness reads KWS=sys.argv[1] as the workspace root (this folder) and argv[2] as the engager.
    sys.argv = [os.path.join(engine, "run_tce_pd_reval.py"), HERE, args.engager]
    runpy_path = os.path.join(engine, "run_tce_pd_reval.py")
    code = compile(open(runpy_path).read(), runpy_path, "exec")
    exec(code, {"__name__": "__main__", "__file__": runpy_path})

if __name__ == "__main__":
    main()
