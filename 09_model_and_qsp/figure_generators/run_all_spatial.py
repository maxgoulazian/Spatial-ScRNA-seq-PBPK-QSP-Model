"""Render ALL organs x ALL canonical molecules, 5-panel spatially-resolved. Loads each pkl once."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pickle, glob, os, datetime, traceback
import gen_spatial_resolved as G
HB="model/rundir/handoff"
# canonical variant per molecule
def rank(t):
    for i,p in enumerate(["full_","fin_","unif_","il6test_","mech_","pd2_","pk_"]):
        if t.startswith(p): return i
    return 50
bynm={}
for f in glob.glob(f"{HB}/tce_spatial_*.pkl"):
    tag=os.path.basename(f)[len("tce_spatial_"):-4]
    # derive molecule name
    base=tag
    for p in ["full_","fin_","cal_","il6test_","mech_","pd2_","unif_","v3_","t7_","rd_","pk_","pkw_"]:
        if base.startswith(p): base=base[len(p):]; break
    bynm.setdefault(base,[]).append(tag)
total=0
for nm in sorted(bynm):
    tag=sorted(bynm[nm],key=rank)[0]
    f=f"{HB}/tce_spatial_{tag}.pkl"
    try:
        d=pickle.load(open(f,"rb"))
        tps=[k for k in d if isinstance(k,(int,float))]
        if not tps: print(f"SKIP {tag}: no timepoints"); continue
        tp=max(tps); org=d[tp]
        for o in list(org.keys()):
            try:
                G.render(nm,tag,o,org[o],tp); total+=1
                print(f"OK {nm}/{o}",flush=True)
            except Exception as e:
                print(f"ERR {nm}/{o}: {str(e)[:80]}",flush=True)
        del d
    except Exception as e:
        print(f"FAIL {tag}: {str(e)[:100]}",flush=True)
print(f"TOTAL spatial figures: {total}")
