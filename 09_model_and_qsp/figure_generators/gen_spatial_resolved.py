"""Spatially-resolved per-cell overlay from tce_spatial_<name>.pkl (read-only).
Per molecule x per organ: 5 panels — cell types | interstitial drug C(nM) | bound antibody(nM) |
receptor R(copies) | cell kill(1-surv). Physical-size cells (EllipseCollection, units='xy') so cells are
visible. Endpoint timepoint. """
import pickle, glob, os, sys, datetime, json
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.collections import EllipseCollection
from matplotlib.lines import Line2D

ROOT="/media/balthasar-lab/RAID4/costim_engager_counterscreen"
HB=f"{ROOT}/model/rundir/handoff"
DATE=datetime.date.today().isoformat(); OUT=f"{ROOT}/spatial_resolved_{DATE}"; os.makedirs(OUT,exist_ok=True)
PAL=json.load(open(f"{ROOT}/SPATIAL_VIEWER_BUNDLE/cell_type_palette.json")) if os.path.exists(f"{ROOT}/SPATIAL_VIEWER_BUNDLE/cell_type_palette.json") else {}

def general(l):
    l=str(l).lower()
    if 'cd8' in l: return 'CD8 T'
    if 'cd4' in l or 'regulatory' in l or ('t cell' in l and 'nk' not in l): return 'CD4 T'
    if 'nk' in l or 'natural killer' in l: return 'NK'
    if 'plasma' in l: return 'Plasma cell'
    if 'b cell' in l: return 'B cell'
    if 'macrophage' in l or 'kupffer' in l: return 'Macrophage'
    if 'monocyte' in l or 'dendritic' in l: return 'Monocyte'
    if 'neutrophil' in l: return 'Neutrophil'
    if 'erythro' in l: return 'Erythrocyte'
    if 'endothel' in l: return 'Endothelial'
    if 'epitheli' in l or 'enterocyte' in l: return 'Epithelial'
    if 'hepatocyte' in l: return 'Hepatocyte'
    if 'malignant' in l or 'tumor' in l or 'cancer' in l: return 'Malignant'
    if 'fibro' in l or 'stromal' in l or 'stellate' in l: return 'Fibroblast'
    return 'Other'

def ecoll(ax,x,y,d,c,a,z):
    ax.add_collection(EllipseCollection(widths=d,heights=d,angles=0,units='xy',
        offsets=np.c_[x,y],transOffset=ax.transData,facecolors=c,edgecolors='none',alpha=a,zorder=z))

def panel_scalar(ax,x,y,val,d,cmap,title,vmin=None,vmax=None,log=False):
    v=np.array(val,float)
    if log:
        v=np.log10(np.clip(v,1e-6,None))
    vmin=np.nanpercentile(v,2) if vmin is None else vmin
    vmax=np.nanpercentile(v,99.5) if vmax is None else vmax
    if vmax<=vmin: vmax=vmin+1e-6
    import matplotlib.cm as cm; import matplotlib.colors as mc
    norm=mc.Normalize(vmin,vmax); sm=cm.ScalarMappable(norm=norm,cmap=cmap)
    cols=sm.to_rgba(v)
    ax.add_collection(EllipseCollection(widths=d,heights=d,angles=0,units='xy',
        offsets=np.c_[x,y],transOffset=ax.transData,facecolors=cols,edgecolors='none',alpha=0.9,zorder=2))
    ax.set_title(title,fontsize=10)
    cb=plt.colorbar(sm,ax=ax,fraction=0.046,pad=0.02); cb.ax.tick_params(labelsize=7)

def render(name, tag, organ, state, tp):
    x=state["x"]; y=state["y"]; C=state["C"]; b=state["bound_nM"]; R=state["R"]
    surv=state["surv"]; labs=state["labs"]; istgt=state.get("is_target")
    n=len(x)
    # physical cell size scaled to field so cells are visible but not overlapping too much
    span=max(x.max()-x.min(), y.max()-y.min())
    d=max(span/380.0, 8.0)   # µm diameter in data coords
    gen=np.array([general(l) for l in labs])
    kill=(1.0-np.array(surv))*100.0
    fig,axs=plt.subplots(1,5,figsize=(30,6.2),dpi=185)
    # 1 cell types
    ax=axs[0]; present=sorted(set(gen),key=lambda t:-(gen==t).sum())
    for t in present:
        m=gen==t; ecoll(ax,x[m],y[m],d,PAL.get(t,'#bbbbbb'),0.85,2)
    ax.set_title("Cell types",fontsize=10)
    counts={t:int((gen==t).sum()) for t in present}
    h=[Line2D([0],[0],marker='o',color='w',markerfacecolor=PAL.get(t,'#bbb'),markersize=7,label=f"{t} ({counts[t]:,})") for t in present[:12]]
    ax.legend(handles=h,fontsize=6,loc='center left',bbox_to_anchor=(0,-0.16),ncol=3,frameon=False)
    # 2 interstitial drug
    panel_scalar(axs[1],x,y,C,d,'viridis',"Interstitial drug C (nM)")
    # 3 bound antibody
    panel_scalar(axs[2],x,y,b,d,'magma',"Bound antibody (nM)")
    # 4 receptor
    panel_scalar(axs[3],x,y,R,d,'cividis',"Receptor R (copies)")
    # 5 kill (structure faded behind, kill in Reds)
    ax=axs[4]
    ecoll(ax,x,y,d,'#dddddd',0.35,1)
    mk=kill>2
    if mk.any():
        import matplotlib.cm as cm; import matplotlib.colors as mc
        norm=mc.Normalize(0,100); sm=cm.ScalarMappable(norm=norm,cmap='Reds')
        ecoll_c=sm.to_rgba(kill[mk])
        ax.add_collection(EllipseCollection(widths=d*1.1,heights=d*1.1,angles=0,units='xy',
            offsets=np.c_[x[mk],y[mk]],transOffset=ax.transData,facecolors=ecoll_c,edgecolors='none',alpha=0.95,zorder=3))
        cb=plt.colorbar(sm,ax=ax,fraction=0.046,pad=0.02); cb.ax.tick_params(labelsize=7)
    ax.set_title(f"Cell kill (%)  — mean {kill.mean():.0f}%",fontsize=10)
    for ax in axs:
        ax.set_xlim(x.min()-50,x.max()+50); ax.set_ylim(y.min()-50,y.max()+50)
        ax.set_aspect('equal'); ax.invert_yaxis(); ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle(f"{name} — {organ}  (day {tp:.0f}, {n:,} cells)  spatially-resolved per-cell state",fontsize=13,y=1.03)
    fig.tight_layout()
    p=f"{OUT}/spatialres_{name}_{organ}.png"; fig.savefig(p,dpi=185,bbox_inches="tight"); plt.close(fig)
    return p

if __name__=="__main__":
    mol=sys.argv[1] if len(sys.argv)>1 else "full_teclistamab"
    organs_arg=sys.argv[2].split(",") if len(sys.argv)>2 else None
    f=f"{HB}/tce_spatial_{mol}.pkl"
    d=pickle.load(open(f,"rb"))
    name=mol.split("_",1)[-1] if mol.split("_")[0] in ("full","fin","cal","il6test","mech","pd2","unif","v3","t7","rd","pk","pkw") else mol
    tps=[k for k in d if isinstance(k,(int,float))]; tp=max(tps)
    org=d[tp]
    organs=organs_arg or list(org.keys())
    made=[]
    for o in organs:
        if o not in org: continue
        try:
            p=render(name,mol,o,org[o],tp); made.append(p); print("OK",os.path.basename(p))
        except Exception as e:
            print("ERR",o,str(e)[:120])
    print("DONE",len(made))
