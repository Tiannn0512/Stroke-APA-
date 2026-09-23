"""Diagnose: astrocytes_annotated.h5ad var columns + GSE286075 DE index format."""
import anndata as ad, pandas as pd
a = ad.read_h5ad(r"D:\stroke_apa\data\GSE174574\h5ad\astrocytes_annotated.h5ad")
print("astro var columns:", list(a.var.columns))
print("astro var_names[:5]:", list(a.var_names[:5]))
r = pd.read_csv(r"D:\stroke_apa\results\p1_DE_paired_GSE286075.tsv", sep="\t", nrows=5)
print("286 index head:", list(r.iloc[:5, 0]))
