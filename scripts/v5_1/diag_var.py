"""Diagnose var_names content in raw_all.h5ad."""
import anndata as ad
a = ad.read_h5ad(r"D:\stroke_apa\data\GSE174574\h5ad\raw_all.h5ad")
print("shape:", a.shape)
print("var_names[:8]:", list(a.var_names[:8]))
print("var columns:", list(a.var.columns))
if "gene_names" in a.var:
    print("gene_names[:8]:", list(a.var["gene_names"][:8]))
# how many look like symbols (contain letters only, short)
import numpy as np
vn = np.array([str(x) for x in a.var_names])
like_symbol = np.char.str_len(vn) <= 12
print("first 30 raw:", vn[:30])
print("Ensembl-like count (starts ENS):", int(np.char.startswith(vn, 'ENS').sum()))
