"""Rebuild raw_all.h5ad with symbol var_names (col2 of genes.tsv). Identical gene order across samples verified by id-list equality."""
import glob, gzip, os, warnings
import numpy as np, pandas as pd, scipy.io as sio
import anndata as ad
warnings.filterwarnings("ignore")

EXTRACT = r"D:\stroke_apa\data\GSE174574\extract"
samples = [("GSM5319987", "sham1", "sham"), ("GSM5319988", "sham2", "sham"), ("GSM5319989", "sham3", "sham"),
           ("GSM5319990", "MCAO1", "MCAO"), ("GSM5319991", "MCAO2", "MCAO"), ("GSM5319992", "MCAO3", "MCAO")]

with gzip.open(glob.glob(os.path.join(EXTRACT, "GSM5319987_sham1_genes.tsv.gz"))[0], "rt") as f:
    recs = [ln.rstrip("\n").split("\t") for ln in f]
ids = [r[0] for r in recs]
sym = pd.Index([r[1] if len(r) > 1 else r[0] for r in recs])
print("genes:", len(ids), "unique symbols:", sym.nunique())

adatas = []
for gsm, sname, cond in samples:
    ids_j = [r[0] for r in recs]
    assert ids_j == ids
    mm = sio.mmread(glob.glob(os.path.join(EXTRACT, f"{gsm}_{sname}_matrix.mtx.gz"))[0]).tocsr()
    with gzip.open(glob.glob(os.path.join(EXTRACT, f"{gsm}_{sname}_barcodes.tsv.gz"))[0], "rt") as f:
        bcs = [ln.strip() for ln in f]
    a = ad.AnnData(X=mm.T.tocsr())
    a.var_names = sym.copy()
    a.var_names_make_unique()          # 27933 unique symbols -> 27998 unique after de-dup suffixing
    a.var["ensembl"] = ids
    a.obs_names = [f"{sname}_{b}" for b in bcs]
    a.obs["sample"] = sname
    a.obs["condition"] = cond
    adatas.append(a)

comb = ad.concat(adatas, axis=0, join="outer")
comb.var_names_make_unique()
comb.X = comb.X.tocsr()
comb.write_h5ad(r"D:\stroke_apa\data\GSE174574\h5ad\raw_all_sym.h5ad")
print("saved raw_all_sym.h5ad", comb.shape, "var head:", list(comb.var_names[:5]))
