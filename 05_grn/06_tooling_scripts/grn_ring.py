"""GRN ring network generator — costim receptor circular layout.
Reconstructs grn_ring_cd8_effector.png: nodes = costim arms (fill = effector drive score,
size = out-edges), directed edges = GRNBoost2 links among arms (red=agonism raises target,
blue=lowers; solid=strong importance, dashed=weaker)."""
import numpy as np, pandas as pd, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle
from matplotlib.cm import ScalarMappable
from matplotlib.colors import TwoSlopeNorm

ARMS={'TNFRSF9':'4-1BB','TNFRSF4':'OX40','TNFRSF18':'GITR','TNFRSF25':'DR3','TNFRSF14':'HVEM',
      'TNFRSF8':'CD30','CD27':'CD27','CD28':'CD28','CD40':'CD40','ICOS':'ICOS','CD2':'CD2',
      'CD226':'DNAM1','LTBR':'LTBR'}

def grn_ring(edges_csv, drive_csv, out_png, title="", dpi=350):
    E=pd.read_csv(edges_csv)                    # TF,target,importance,corr_sign,signed_importance,validated
    D=pd.read_csv(drive_csv)                     # arm,effector_drive_score,...
    dmap=dict(zip(D['arm'],D['effector_drive_score']))
    arms=list(ARMS)
    # keep only arm->arm edges (ring is among costim receptors)
    EA=E[E['TF'].isin(arms) & E['target'].isin(arms)].copy()
    # node out-edge count -> size
    outdeg={a:int((EA['TF']==a).sum()) for a in arms}
    # circular layout
    n=len(arms); ang=np.linspace(90,90-360,n,endpoint=False)*np.pi/180
    pos={a:(np.cos(ang[i]),np.sin(ang[i])) for i,a in enumerate(arms)}
    fig,ax=plt.subplots(figsize=(11,11)); ax.set_aspect('equal'); ax.axis('off')
    # ---- node radii (data coords) computed FIRST so edges can terminate outside circles ----
    maxdeg=max(outdeg.values() or [1])
    radii={a: 0.055+0.045*(outdeg[a]/maxdeg) for a in arms}
    # ---- edges (trimmed to the target circle boundary + gap so arrowheads are never covered) ----
    imp=EA['importance'].abs(); imax=imp.max() if len(imp) else 1.0
    ARROW_GAP=0.030                                   # extra clearance beyond the circle edge
    for _,r in EA.iterrows():
        s=np.array(pos[r['TF']],float); t=np.array(pos[r['target']],float)
        d=t-s; L=np.hypot(*d)
        if L<1e-9: continue
        u=d/L
        s_trim=s+u*(radii[r['TF']]+0.006)
        t_trim=t-u*(radii[r['target']]+ARROW_GAP)
        sign=r.get('corr_sign', np.sign(r.get('signed_importance',1)))
        col='#c0392b' if sign>=0 else '#2c6fbb'          # red raises / blue lowers
        w=0.5+3.4*(abs(r['importance'])/imax)
        strong=abs(r['importance'])>=0.5*imax
        ax.add_patch(FancyArrowPatch(tuple(s_trim),tuple(t_trim),connectionstyle="arc3,rad=0.12",
            arrowstyle='-|>',mutation_scale=13+9*strong,lw=w,color=col,
            alpha=0.9 if strong else 0.45,ls='-' if strong else (0,(4,2)),
            zorder=2,shrinkA=0,shrinkB=0))
    # ---- nodes (drawn AFTER edges; arrowheads already stop outside, so nothing is covered) ----
    norm=TwoSlopeNorm(vmin=-9,vcenter=0,vmax=9); cmap=plt.cm.RdBu_r
    for a in arms:
        x,y=pos[a]; sc=dmap.get(a,0.0); rad=radii[a]
        ax.add_patch(Circle((x,y),rad,facecolor=cmap(norm(sc)),ec='#333',lw=1.3,zorder=3))
        ax.text(x,y,ARMS[a],ha='center',va='center',fontsize=9,fontweight='bold',
                color='white' if abs(sc)>4 else '#111',zorder=4)
    ax.set_xlim(-1.35,1.35); ax.set_ylim(-1.35,1.35)
    sm=ScalarMappable(norm=norm,cmap=cmap); sm.set_array([])
    cb=fig.colorbar(sm,ax=ax,fraction=0.035,pad=0.02); cb.set_label("effector drive score",fontsize=10)
    fig.suptitle(title,fontsize=12,y=0.97)
    # legend
    from matplotlib.lines import Line2D
    leg=[Line2D([0],[0],color='#c0392b',lw=2.5,label='agonism raises target'),
         Line2D([0],[0],color='#2c6fbb',lw=2.5,label='agonism lowers target'),
         Line2D([0],[0],color='#555',lw=2.5,ls='-',label='strong edge (importance)'),
         Line2D([0],[0],color='#555',lw=1.2,ls=(0,(4,2)),label='weaker edge'),
         Line2D([0],[0],marker='o',color='w',markerfacecolor='#ccc',markersize=13,label='size ∝ out-edges')]
    ax.legend(handles=leg,loc='lower center',bbox_to_anchor=(0.5,-0.08),ncol=5,frameon=False,fontsize=8)
    fig.savefig(out_png,dpi=dpi,bbox_inches='tight'); plt.close(fig)
    return out_png

if __name__=="__main__":
    print("grn_ring module OK")
