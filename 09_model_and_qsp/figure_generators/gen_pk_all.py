"""PK overlays for ALL molecules with clinical PK anchor points (digitized OR label-summary). Read-only."""
import json, glob, os, sqlite3, datetime
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT="."; HB="model/rundir/handoff"; DB="model/params/mab_tce_pkpd.sqlite"
DATE=datetime.date.today().isoformat(); OUT=f"overlays_{DATE}"; os.makedirs(OUT,exist_ok=True)
con=sqlite3.connect(DB); con.row_factory=sqlite3.Row
drugs={r["drug_id"]:r["inn"] for r in con.execute("SELECT drug_id,inn FROM drugs")}
name2id={v.lower():k for k,v in drugs.items() if v}
def pkcurves(nm):
    did=name2id.get(nm.lower())
    if not did: return []
    out=[]
    for r in con.execute("SELECT * FROM curves WHERE drug_id=? AND readout_class='PK'",(did,)):
        ts=con.execute("SELECT time,value,time_unit FROM timeseries WHERE curve_id=? ORDER BY point_order",(r["curve_id"],)).fetchall()
        if ts: out.append({"lab":r["dose_label"],"cu":r["conc_unit"],"tu":r["time_unit"],"dig":r["digitized"],
                           "t":[x["time"] for x in ts],"v":[x["value"] for x in ts]})
    return out
def rank(t):
    for i,p in enumerate(["fin_","full_","unif_","pk_","il6test_"]):
        if t.startswith(p): return i
    return 50
bynm={}
for f in glob.glob(f"{HB}/tce_pd_*.json"):
    if "_PREFIX" in f or "reval_result" in f: continue
    try:d=json.load(open(f))
    except:continue
    nm=d.get("name"); tag=os.path.basename(f)[7:-5]
    if not nm or not d.get("Cplasma_ugml"): continue
    if "_" in nm and nm.split("_")[-1] in ("sc","iv","low"): continue
    bynm.setdefault(nm,{})[tag]=d
made=[]
for nm,variants in bynm.items():
    cc=pkcurves(nm)
    if not cc: continue
    ct=sorted(variants,key=rank)[0]; d=variants[ct]
    t=np.array(d["t_pk"]); C=np.array(d["Cplasma_ugml"])
    fig,ax=plt.subplots(figsize=(7.5,5),dpi=190)
    ax.plot(t,C,'-',color="#1f4e79",lw=2.4,label=f"model {nm}",zorder=3)
    seen=set(); cols=plt.cm.autumn(np.linspace(0,.75,max(len(cc),1)))
    for i,c in enumerate(cc):
        key=(c["lab"],tuple(c["t"]))
        if key in seen: continue
        seen.add(key)
        tt=np.array(c["t"],float); vv=np.array(c["v"],float)
        u=(c["cu"] or "ug/mL").lower()
        if "ng" in u: vv/=1000.0
        elif "pg" in u: vv/=1e6
        style='o' if c["dig"] else 'D'
        lab=f"clin{'(fig)' if c['dig'] else '(label pts)'}: {(c['lab'] or '')[:26]}"
        ax.plot(tt,vv,style,ms=6,color=cols[i%len(cols)],label=lab,alpha=.9,zorder=4,mec='k',mew=0.4)
    ax.set_yscale("log"); ax.set_xlabel("time (day)"); ax.set_ylabel("plasma conc (µg/mL)")
    ax.set_title(f"PK overlay — {nm}  (model vs clinical)"); ax.legend(fontsize=6.5); ax.grid(alpha=.3,which="both")
    p=f"{OUT}/PK_overlay_{nm}.png"; fig.savefig(p,bbox_inches="tight"); plt.close(fig); made.append((nm,len(seen)))
con.close()
print(f"PK overlays: {len(made)}")
for nm,n in sorted(made): print(f"  {nm:16s} {n} clinical anchor(s)")
