import os
os.environ["OMP_NUM_THREADS"]="1"; os.environ["OPENBLAS_NUM_THREADS"]="1"; os.environ["MKL_NUM_THREADS"]="1"
os.environ["NUMEXPR_MAX_THREADS"]="16"
import numpy as np, pandas as pd, time
import arboreto.core as ac
import dask.dataframe as dd
_orig = ac.from_delayed
def _patched(dfs, *a, **k):
    if isinstance(dfs,(list,tuple)) and len(dfs)==0:
        meta=k.get("meta");
        if meta is None and a: meta=a[0]
        return dd.from_pandas(meta.copy(), npartitions=1)
    return _orig(dfs,*a,**k)
ac.from_delayed=_patched

def main():
    from arboreto.algo import grnboost2
    from distributed import Client, LocalCluster
    cl=LocalCluster(n_workers=2,threads_per_worker=1,processes=True,memory_limit="2GB",dashboard_address=None)
    c=Client(cl)
    rng=np.random.default_rng(0)
    G=["r0","r1","r2","t0","t1"]
    Xd=rng.normal(size=(300,5)).astype("float32")
    Xd[:,3]=0.8*Xd[:,0]+0.2*rng.normal(size=300)
    ex=pd.DataFrame(Xd,columns=G)
    t=time.time()
    net=grnboost2(expression_data=ex, tf_names=["r0","r1","r2"], client_or_address=c, seed=0, verbose=False)
    print("grnboost2 OK", net.shape, "in %.1fs"%(time.time()-t))
    print(net.sort_values("importance",ascending=False).to_string())
    c.close(); cl.close()
    print("PATCHTEST PASS")

if __name__=="__main__":
    main()
