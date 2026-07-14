"""ALL PK + PD overlays: model vs every digitized clinical curve (any analyte). Read-only handoff."""
import json, glob, os, sqlite3, datetime
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT="/media/balthasar-lab/RAID4/costim_engager_counterscreen"
HB=f"{ROOT}/model/rundir/handoff"; DB=f"{ROOT}/model/params/mab_tce_pkpd.sqlite"
DATE=datetime.date.today().isoformat(); OUT=f"{ROOT}/overlays_{DATE}"; os.makedirs(OUT,exist_ok=True)
con=sqlite3.connect(DB); con.row_factory=sqlite3.Row
drugs={r["drug_id"]:r["inn"] for r in con.execute("SELECT drug_id,inn FROM drugs")}
name2id={v.lower():k for k,v in drugs.items() if v}
def curves(drug,rc,analyte=None):
    did=name2id.get(drug.lower())
    if not did: return []
    out=[]
    for r in con.execute("SELECT * FROM curves WHERE drug_id=? AND readout_class=? AND digitized=1",(did,rc)):
        if analyte and analyte.lower() not in (r["analyte"] or "").lower(): continue
        ts=con.execute("SELECT time,value,time_unit FROM timeseries WHERE curve_id=? ORDER BY point_order",(r["curve_id"],)).fetchall()
        if ts: out.append({"lab":r["dose_label"] or r["analyte"],"an":r["analyte"],"cu":r["conc_unit"],"tu":r["time_unit"],
                           "t":[x["time"] for x in ts],"v":[x["value"] for x in ts]})
    return out
def rank(t):
    for i,p in enumerate(["fin_","full_","unif_","il6test_","pk_"]):
        if t.startswith(p): return i
    return 50
bynm={}
for f in glob.glob(f"{HB}/tce_pd_*.json"):
    if "_PREFIX" in f or "reval_result" in f: continue
    try:d=json.load(open(f))
    except:continue
    nm=d.get("name"); tag=os.path.basename(f)[7:-5]
    if not nm: continue
    if "_" in nm and nm.split("_")[-1] in ("sc","iv","low"): continue
    bynm.setdefault(nm,{})[tag]=d
made=[]
for nm,variants in bynm.items():
    # PK overlay
    ct=sorted(variants,key=rank)[0]; d=variants[ct]
    pk=curves(nm,"PK")
    if pk and d.get("Cplasma_ugml"):
        t=np.array(d["t_pk"]); C=np.array(d["Cplasma_ugml"])
        fig,ax=plt.subplots(figsize=(7,5),dpi=190)
        ax.plot(t,C,'-',color="#1f4e79",lw=2.4,label=f"model {nm}",zorder=3)
        cols=plt.cm.autumn(np.linspace(0,.75,len(pk)))
        for i,c in enumerate(pk):
            tt=np.array(c["t"],float); vv=np.array(c["v"],float)
            if c["cu"] and "ng" in c["cu"].lower(): vv/=1000.0
            ax.plot(tt,vv,'o',ms=5,color=cols[i],label=f"clin:{c['lab'][:32]}",alpha=.9,zorder=4)
        ax.set_yscale("log");ax.set_xlabel("day");ax.set_ylabel("plasma µg/mL")
        ax.set_title(f"PK — {nm}");ax.legend(fontsize=7);ax.grid(alpha=.3,which="both")
        p=f"{OUT}/PK_overlay_{nm}.png";fig.savefig(p,bbox_inches="tight");plt.close(fig);made.append(("PK",nm,ct))
    # PD IL-6
    ct2=sorted(variants,key=lambda t:(0 if t.startswith("il6test_") else rank(t)))[0]; d2=variants[ct2]
    pd=curves(nm,"PD","IL-6")
    if pd and d2.get("il6_pgml"):
        t=np.array(d2["t"]); il6=np.array(d2["il6_pgml"])
        fig,ax=plt.subplots(figsize=(7,5),dpi=190)
        ax.plot(t,il6,'-',color="#7b241c",lw=2.4,label=f"model {nm} IL-6",zorder=3)
        cols=plt.cm.winter(np.linspace(0,.8,len(pd)))
        for i,c in enumerate(pd):
            tt=np.array(c["t"],float);vv=np.array(c["v"],float)
            if (c["tu"] or "day").lower().startswith("h"): tt/=24.0
            ax.plot(tt,vv,'s',ms=6,color=cols[i],label=f"clin:{c['lab'][:30]}",alpha=.9,zorder=4)
        ax.set_xlabel("day");ax.set_ylabel("IL-6 pg/mL");ax.set_title(f"PD IL-6 — {nm}");ax.legend(fontsize=7);ax.grid(alpha=.3)
        p=f"{OUT}/PD_IL6_overlay_{nm}.png";fig.savefig(p,bbox_inches="tight");plt.close(fig);made.append(("PD-IL6",nm,ct2))
con.close()
print(f"generated {len(made)} overlays")
for k,nm,tg in sorted(made): print(f"  {k:8s} {nm} ({tg})")
