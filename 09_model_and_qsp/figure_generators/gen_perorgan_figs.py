"""Per-molecule, per-organ kill + cytokine figures from tce_pd_<name>.json (read-only).
One 3-panel figure per molecule: (A) per-organ kill trajectories, (B) final kill bar (target-cell-weighted),
(C) per-organ IL-6 time-course. Canonical result per molecule."""
import json, glob, os, datetime
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT="/media/balthasar-lab/RAID4/costim_engager_counterscreen"
HB=f"{ROOT}/model/rundir/handoff"
DATE=datetime.date.today().isoformat(); OUT=f"{ROOT}/perorgan_{DATE}"; os.makedirs(OUT,exist_ok=True)

# canonical result per molecule (prefer full/fin production run)
def rank(tag):
    for i,p in enumerate(["full_","fin_","unif_","il6test_","pk_"]):
        if tag.startswith(p): return i
    return 50
byname={}
for f in glob.glob(f"{HB}/tce_pd_*.json"):
    if "_PREFIX" in f or "reval_result" in f: continue
    try:d=json.load(open(f))
    except:continue
    nm=d.get("name"); tag=os.path.basename(f)[7:-5]
    if not nm or ("kill_traj" not in d): continue
    if "_" in nm and nm.split("_")[-1] in ("sc","iv","low"): continue
    byname.setdefault(nm,{})[tag]=d

# fixed organ colors
ORG=['tumor','spleen','bone','liver','lung','large_int','small_int','pancreas','kidney','skin','heart','adipose','blood']
CMAP={o:c for o,c in zip(ORG,plt.cm.tab20(np.linspace(0,1,len(ORG))))}
made=[]; summ=[]
for nm,variants in byname.items():
    tag=sorted(variants,key=rank)[0]; d=variants[tag]
    t=np.array(d["t"]); ct=d.get("cancer_type","?")
    kt=d["kill_traj"]; kf=d.get("kill_final",{}); knt=d.get("kill_ntgt",{}); co=d.get("cyto_organ",{})
    organs=[o for o in ORG if o in kt]+[o for o in kt if o not in ORG]
    fig,axs=plt.subplots(1,3,figsize=(17,5),dpi=190)
    # A: kill trajectories
    ax=axs[0]
    for o in organs:
        y=np.array(kt[o])
        if y.max()>0.005: ax.plot(t,y*100,lw=1.8,color=CMAP.get(o,'#888'),label=o)
    ax.set_xlabel("time (day)"); ax.set_ylabel("target-cell kill (%)"); ax.set_ylim(-2,102)
    ax.set_title(f"A. Per-organ kill trajectory"); ax.legend(fontsize=7,ncol=2); ax.grid(alpha=0.3)
    # B: final kill bar, sized by target-cell count
    ax=axs[1]
    bars=[(o,kf.get(o,0)*100,knt.get(o,0)) for o in organs if (kf.get(o,0)>0.005 or knt.get(o,0)>0)]
    bars.sort(key=lambda x:-x[1])
    if bars:
        names=[b[0] for b in bars]; vals=[b[1] for b in bars]; nt=[b[2] for b in bars]
        cols=[CMAP.get(n,'#888') for n in names]
        yb=ax.barh(range(len(names)),vals,color=cols); ax.set_yticks(range(len(names))); ax.set_yticklabels(names,fontsize=8)
        ax.invert_yaxis(); ax.set_xlabel("final kill (%)"); ax.set_xlim(0,105)
        for i,(v,c) in enumerate(zip(vals,nt)):
            ax.text(min(v+1,100),i,f"{v:.0f}%  (n={c:,})",va="center",fontsize=6.5)
    ax.set_title(f"B. Final kill by organ (n=target cells)"); ax.grid(alpha=0.3,axis="x")
    # C: per-organ IL-6
    ax=axs[2]
    for o in organs:
        if o in co and isinstance(co[o],dict) and "IL6" in co[o]:
            y=np.array(co[o]["IL6"])
            if y.max()>1: ax.plot(t[:len(y)],y,lw=1.8,color=CMAP.get(o,'#888'),label=o)
    ax.set_xlabel("time (day)"); ax.set_ylabel("IL-6 production rate (organ)")
    ax.set_title("C. Per-organ IL-6"); ax.legend(fontsize=7,ncol=2); ax.grid(alpha=0.3)
    dw=d.get("depletion_weighted"); tit=f"{nm}  —  {ct} tumor  |  depletion(wt)={dw:.2f}" if dw is not None else f"{nm} — {ct}"
    fig.suptitle(tit,fontsize=13,y=1.02); fig.tight_layout()
    p=f"{OUT}/perorgan_{nm}.png"; fig.savefig(p,bbox_inches="tight"); plt.close(fig); made.append(p)
    summ.append((nm,ct,tag,len([b for b in bars]),dw))
print(f"generated {len(made)} per-molecule per-organ figures in perorgan_{DATE}/")
for nm,ct,tag,nb,dw in sorted(summ): print(f"  {nm:18s} {ct:5s} organs={nb:2d} depl={dw if dw is None else round(dw,2)}  ({tag})")
