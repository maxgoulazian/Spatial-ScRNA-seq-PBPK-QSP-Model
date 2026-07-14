import sys, os, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_spatial_resolved as G
HB="model/rundir/handoff"
tag=sys.argv[1]  # e.g. full_teclistamab
base=tag
for p in ["full_","fin_","cal_","il6test_","mech_","pd2_","unif_","v3_","t7_","rd_","pk_","pkw_"]:
    if base.startswith(p): base=base[len(p):]; break
f=f"{HB}/tce_spatial_{tag}.pkl"
try:
    d=pickle.load(open(f,"rb"))
    tps=[k for k in d if isinstance(k,(int,float))]
    if not tps: print(f"SKIP {tag} no tp"); sys.exit()
    tp=max(tps); org=d[tp]
    done=0
    for o in list(org.keys()):
        outp=f"spatial_resolved_2026-07-13/spatialres_{base}_{o}.png"
        if os.path.exists(outp): done+=1; continue  # skip already-rendered
        try: G.render(base,tag,o,org[o],tp); done+=1
        except Exception as e: print(f"ERR {base}/{o}: {str(e)[:60]}")
    print(f"DONE {base}: {done} organs")
except Exception as e:
    print(f"FAIL {tag}: {str(e)[:80]}")
