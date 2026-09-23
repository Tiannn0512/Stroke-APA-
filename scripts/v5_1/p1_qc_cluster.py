"""P1-02/03: QC (MAD fallback, scrublet unavailable - logged), clustering, astrocyte calling, reactive score, pseudo-bulk."""
import os, re, warnings, time
import numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad, scanpy as sc
warnings.filterwarnings("ignore")
sc.settings.verbosity = 1

H5 = r"D:\stroke_apa\data\GSE174574\h5ad\raw_all_sym.h5ad"
OUT = r"D:\stroke_apa\results"
FIG = os.path.join(OUT, "p1_astrocyte_clusters")
os.makedirs(FIG, exist_ok=True)
t0 = time.time()

a = ad.read_h5ad(H5)
gnames = np.array([str(x) for x in a.var_names])
mt_mask = np.char.startswith(gnames, "mt-")
print("loaded", a.shape, f"{time.time()-t0:.0f}s", flush=True)

# ---------- P1-02 QC (frozen params; scrublet fallback -> MAD) ----------
counts = np.asarray(a.X.sum(axis=1)).ravel().astype(float)
ngenes = np.asarray((a.X > 0).sum(axis=1)).ravel().astype(float)
mt_frac = np.asarray(a.X[:, mt_mask].sum(axis=1)).ravel() / np.maximum(counts, 1)

keep_base = (ngenes >= 200) & (counts >= 500) & (mt_frac < 0.20)
# MAD outliers per sample on log10 scale (frozen fallback)
lc, lg = np.log10(counts + 1), np.log10(ngenes + 1)
keep_mad = np.ones(len(counts), bool)
for s in a.obs["sample"].unique():
    m = (a.obs["sample"] == s).values
    for v in (lc[m], lg[m]):
        med, mad = np.median(v), np.median(np.abs(v - np.median(v)))
        lo, hi = med - 3 * mad / 0.6745, med + 3 * mad / 0.6745
        idx = np.where(m)[0]
        keep_mad[idx] &= (v >= lo) & (v <= hi)
a.obs["passed_qc"] = keep_base & keep_mad
qc_tab = a.obs.groupby(["sample", "condition"])["passed_qc"].agg(["size", "sum"])
qc_tab.columns = ["total", "passed"]
qc_tab.to_csv(os.path.join(OUT, "p1_qc_summary.tsv"), sep="\t")
print(qc_tab.to_string(), flush=True)

adata = a[a.obs["passed_qc"].values].copy()
adata.layers["counts"] = adata.X.copy()
gnames = np.array([str(x) for x in adata.var_names])

# ---------- normalize / HVG / PCA / leiden ----------
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat", batch_key="sample")
sc.tl.pca(adata, n_comps=50, use_highly_variable=True, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.leiden(adata, resolution=1.0, key_added="leiden", flavor="igraph", n_iterations=2, directed=False)
sc.tl.umap(adata, min_dist=0.3)
print("clustered", adata.shape, f"leiden clusters={adata.obs['leiden'].nunique()}", f"{time.time()-t0:.0f}s", flush=True)

# ---------- cluster x marker table (log1p CPM=1e4 scale) ----------
markers = ["Aqp4","Slc1a3","Gfap","Aldh1l1","Serpina3n","C3","C1qa","Cx3cr1","Pecam1","Mbp","Rbfox3","Syt1"]
name_to_ix = {}
for i, g in enumerate(gnames):
    name_to_ix.setdefault(str(g), i)
rows = []
for cl in sorted(adata.obs["leiden"].unique(), key=int):
    m = (adata.obs["leiden"] == cl).values
    row = dict(cluster=cl, n_cells=int(m.sum()))
    for g in markers:
        ix = name_to_ix.get(g)
        row[g] = float(np.asarray(adata.X[m, ix].todense()).mean()) if ix is not None else np.nan
    rows.append(row)
mtab = pd.DataFrame(rows).set_index("cluster")
mtab.to_csv(os.path.join(OUT, "p1_cluster_markers.tsv"), sep="\t")
print(mtab.round(2).to_string(), flush=True)

# ---------- astrocyte cluster calling (frozen rules) ----------
astro_clusters, excluded = [], []
for cl, r in mtab.iterrows():
    is_astro = (r["Aqp4"] >= 1 and r["Slc1a3"] >= 1) or (r["Aqp4"] >= 0.5 and r["Gfap"] >= 1 and r["Slc1a3"] >= 0.5)
    immune = (r["C1qa"] >= 1 and r["Cx3cr1"] >= 1) or (r["Pecam1"] >= 1)
    if is_astro and not immune:
        astro_clusters.append(cl)
    elif is_astro and immune:
        excluded.append((cl, "astro-marker-high-but-immune-marker-high"))
print("astrocyte clusters:", astro_clusters, "excluded:", excluded, flush=True)

# ---------- cell-level reactive score on astrocyte cells (frozen gene set) ----------
astro_mask = adata.obs["leiden"].isin(astro_clusters).values
ad_astro = adata[astro_mask].copy()
sc.tl.score_genes(ad_astro, gene_list=["Serpina3n", "C3", "Gfap"], score_name="reactive_score", ctrl_size=50, n_bins=25)
ad_astro.obs["astro_subpop"] = np.where(ad_astro.obs["reactive_score"] > 0.5, "reactive", "homeostatic")
sub_tab = pd.crosstab(ad_astro.obs["sample"], ad_astro.obs["astro_subpop"])
sub_tab.to_csv(os.path.join(OUT, "p1_astrocyte_subpop_counts.tsv"), sep="\t")
print(sub_tab.to_string(), flush=True)

# save annotated astrocyte h5ad (raw counts layer intact)
ad_astro.write_h5ad(r"D:\stroke_apa\data\GSE174574\h5ad\astrocytes_annotated.h5ad")
# also save full clustered object (obs only essentials to keep size manageable)
adata.obs[["sample","condition","leiden","passed_qc"]].to_csv(os.path.join(OUT, "p1_all_cells_clusters.tsv"), sep="\t")

# figures
sc.settings.figdir = FIG
sc.pl.umap(adata, color=["leiden", "sample", "condition"], save="_all_clusters.png", show=False)
sc.pl.umap(ad_astro, color=["astro_subpop", "reactive_score", "sample"], save="_astro_subpop.png", show=False)

# ---------- pseudo-bulk (raw counts, sample x subpop) ----------
pb_rows, meta = [], []
Xc = ad_astro.layers["counts"]
obs = ad_astro.obs
for (s, sub), grp in obs.groupby(["sample", "astro_subpop"]):
    idx = np.where((obs["sample"] == s) & (obs["astro_subpop"] == sub))[0]
    if len(idx) < 20:
        print(f"unit {s}/{sub} n={len(idx)} <20 -> excluded from DE", flush=True)
        continue
    vec = np.asarray(Xc[idx].sum(axis=0)).ravel()
    pb_rows.append(vec)
    meta.append(dict(sample=s, subpop=sub, condition=("MCAO" if s.startswith("MCAO") else "sham"), n_cells=len(idx)))
pb = pd.DataFrame(pb_rows, index=[f"{m['sample']}__{m['subpop']}" for m in meta], columns=adata.var_names)
pb.to_csv(os.path.join(OUT, "p1_pseudobulk_astro_counts.tsv"), sep="\t")
pd.DataFrame(meta).to_csv(os.path.join(OUT, "p1_pseudobulk_meta.tsv"), sep="\t", index=False)
print("pseudo-bulk units:", len(meta), f"total {time.time()-t0:.0f}s", flush=True)
