# -*- coding: utf-8 -*-
"""P1-07: doublet + ambient sensitivity for GSE174574 (audit 7.3).
- Doublet: scanpy-native sc.pp.scrublet (batch_key=sample) on the full QC-passed object.
- Ambient: light proxy score (raw matrix unavailable from GEO -- filtered mtx only, limitation
  logged): per-sample ambient profile from bottom-5%-count cells; per-cell ambient score =
  fraction of counts from ambient-enriched genes; top-5% scorers flagged ambient_high.
- Variants: B = frozen QC minus scrublet doublets; C = B minus ambient_high.
- Each variant: same frozen pipeline (HVG2000/PCA50/leiden res1.0/astro call/reactive score/
  pseudo-bulk pydeseq2 overall) -> compare astrocyte counts, subpop, DE vs original (variant A).
Outputs: results/p1_07_sensitivity_summary.tsv + p1_07_variant_DE_<v>.tsv + md lines to log.
"""
import os, time, warnings
import numpy as np, pandas as pd
import anndata as ad, scanpy as sc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
warnings.filterwarnings("ignore")
sc.settings.verbosity = 1

H5 = r"D:\stroke_apa\data\GSE174574\h5ad\raw_all_sym.h5ad"
OUT = r"D:\stroke_apa\results"
ORIG_DE = os.path.join(OUT, "p1_DE_overall_astro.tsv")
ORIG_CAT = os.path.join(OUT, "p1_regulator_catalog.tsv")
ORIG_SUB = os.path.join(OUT, "p1_astrocyte_subpop_counts.tsv")
t0 = time.time()

REG = ["Pabpn1","Qk","Nudt21","Cpsf6","Cpsf4","Cpsf4l","Cpsf3","Cpsf1","Cpsf2","Cstf1","Cstf2",
       "Cstf3","Clp1","Wdr33","Fip1l1","Pcf11","Symplek","Elavl1","Elavl2","Ptbp1","Ptbp2",
       "Hnrnpa2b1","Srsf1","Nova1","Nova2","Fus","Tardbp","Khdrbs1","Mbnl1","Mbnl2","Rbfox1"]

a = ad.read_h5ad(H5)
gnames = np.array([str(x) for x in a.var_names])
mt_mask = np.char.startswith(gnames, "mt-")
counts = np.asarray(a.X.sum(axis=1)).ravel().astype(float)
ngenes = np.asarray((a.X > 0).sum(axis=1)).ravel().astype(float)
mt_frac = np.asarray(a.X[:, mt_mask].sum(axis=1)).ravel() / np.maximum(counts, 1)
keep_base = (ngenes >= 200) & (counts >= 500) & (mt_frac < 0.20)
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
adata = a[a.obs["passed_qc"].values].copy()
adata.layers["counts"] = adata.X.copy()
print(f"QC passed: {adata.n_obs} cells {time.time()-t0:.0f}s", flush=True)

# ---- doublet: scanpy-native scrublet ----
sc.pp.scrublet(adata, batch_key="sample")
n_dbl = int(adata.obs["predicted_doublet"].sum())
print(f"scrublet predicted doublets: {n_dbl} ({n_dbl/adata.n_obs*100:.2f}%) "
      f"mean score={adata.obs['doublet_score'].mean():.4f} {time.time()-t0:.0f}s", flush=True)

# ---- ambient proxy score ----
X = adata.X
tot = np.asarray(X.sum(axis=1)).ravel()
ambient_flag = np.zeros(adata.n_obs, bool)
amb_genes_all = set()
for s in adata.obs["sample"].unique():
    m = np.where((adata.obs["sample"] == s).values)[0]
    order = m[np.argsort(tot[m])]
    low = order[:max(20, int(len(order) * 0.05))]
    hi = order[-max(50, int(len(order) * 0.5)):]
    prof_low = np.asarray(X[low].sum(axis=0)).ravel()
    prof_hi = np.asarray(X[hi].sum(axis=0)).ravel()
    # genes enriched in low-count pool vs high-count pool (ambient signature candidates)
    frac_lo = prof_low / max(prof_low.sum(), 1)
    frac_hi = prof_hi / max(prof_hi.sum(), 1)
    amb = frac_lo > 3 * frac_hi
    amb_ix = np.where(amb)[0]
    for i in amb_ix:
        amb_genes_all.add(str(gnames[i]))
    sub_tot = np.asarray(X[m][:, amb_ix].sum(axis=1)).ravel() if len(amb_ix) else np.zeros(len(m))
    frac_amb = sub_tot / np.maximum(tot[m], 1)
    thr = np.quantile(frac_amb, 0.95)
    ambient_flag[m[frac_amb >= thr]] = True
adata.obs["ambient_high"] = ambient_flag
n_amb = int(ambient_flag.sum())
print(f"ambient proxy: {len(amb_genes_all)} ambient-enriched genes; flagged {n_amb} cells top-5% "
      f"({n_amb/adata.n_obs*100:.2f}%) {time.time()-t0:.0f}s", flush=True)

# ---- variants ----
variants = {
    "B_minus_doublet": ~adata.obs["predicted_doublet"].values,
    "C_minus_doublet_ambient": (~adata.obs["predicted_doublet"].values) & (~adata.obs["ambient_high"].values),
}
orig_de = pd.read_csv(ORIG_DE, sep="\t", index_col=0)
orig_cat = pd.read_csv(ORIG_CAT, sep="\t")
orig_sub = pd.read_csv(ORIG_SUB, sep="\t", index_col=0)
summary_rows = [dict(variant="A_original", n_cells_used=int(adata.n_obs), removed=0,
                     astro_cells=int(orig_sub.sum().sum()))]

def pipeline(variant, mask):
    ad_v = adata[mask].copy()
    ad_v.X = ad_v.layers["counts"].copy() if "counts" in ad_v.layers else ad_v.X.copy()
    # normalize on the variant's own cell set
    ad_v.X = ad_v.layers["counts"].copy()
    sc.pp.normalize_total(ad_v, target_sum=1e4)
    sc.pp.log1p(ad_v)
    sc.pp.highly_variable_genes(ad_v, n_top_genes=2000, flavor="seurat", batch_key="sample")
    sc.tl.pca(ad_v, n_comps=50, use_highly_variable=True, svd_solver="arpack")
    sc.pp.neighbors(ad_v, n_neighbors=15, n_pcs=30)
    sc.tl.leiden(ad_v, resolution=1.0, key_added="leiden", flavor="igraph", n_iterations=2, directed=False)
    gnames_v = np.array([str(x) for x in ad_v.var_names])
    ix = {}
    for i, g in enumerate(gnames_v):
        ix.setdefault(g, i)
    rows = []
    for cl in sorted(ad_v.obs["leiden"].unique(), key=int):
        m = (ad_v.obs["leiden"] == cl).values
        row = dict(cluster=cl, n_cells=int(m.sum()))
        for g in ["Aqp4", "Slc1a3", "Gfap", "C1qa", "Cx3cr1", "Pecam1"]:
            j = ix.get(g)
            row[g] = float(np.asarray(ad_v.X[m, j].todense()).mean()) if j is not None else np.nan
        rows.append(row)
    mtab = pd.DataFrame(rows).set_index("cluster")
    astro, _ = [], []
    for cl, r in mtab.iterrows():
        is_astro = (r["Aqp4"] >= 1 and r["Slc1a3"] >= 1) or (r["Aqp4"] >= 0.5 and r["Gfap"] >= 1 and r["Slc1a3"] >= 0.5)
        immune = (r["C1qa"] >= 1 and r["Cx3cr1"] >= 1) or (r["Pecam1"] >= 1)
        if is_astro and not immune:
            astro.append(cl)
    am = ad_v[ad_v.obs["leiden"].isin(astro).values].copy()
    sc.tl.score_genes(am, gene_list=["Serpina3n", "C3", "Gfap"], score_name="reactive_score",
                      ctrl_size=50, n_bins=25)
    am.obs["astro_subpop"] = np.where(am.obs["reactive_score"] > 0.5, "reactive", "homeostatic")
    sub = pd.crosstab(am.obs["sample"], am.obs["astro_subpop"])
    # pseudo-bulk DE (overall) on raw counts
    Xc = am.layers["counts"] if "counts" in am.layers else am.X
    pb_rows, meta = [], []
    for s, grp in am.obs.groupby("sample"):
        idx = np.where((am.obs["sample"] == s).values)[0]
        if len(idx) < 20:
            continue
        vec = np.asarray(Xc[idx].sum(axis=0)).ravel()
        pb_rows.append(vec)
        meta.append(dict(sample=s, condition=("MCAO" if s.startswith("MCAO") else "sham")))
    pb = pd.DataFrame(pb_rows, index=[m_["sample"] for m_ in meta], columns=ad_v.var_names)
    pb = pb.loc[:, pb.sum(axis=0) >= 10].astype(int)
    md = pd.DataFrame(index=pb.index)
    md["condition"] = ["MCAO" if s.startswith("MCAO") else "sham" for s in pb.index]
    dds = DeseqDataSet(counts=pb, metadata=md, design="~condition", quiet=True)
    dds.deseq2()
    sts = DeseqStats(dds, contrast=["condition", "MCAO", "sham"], quiet=True)
    sts.summary()
    res = sts.results_df
    res.to_csv(os.path.join(OUT, f"p1_07_variant_DE_{variant}.tsv"), sep="\t")
    sub.to_csv(os.path.join(OUT, f"p1_07_variant_subpop_{variant}.tsv"), sep="\t")
    n_sig = int((res["padj"] < 0.05).sum())
    sig0 = set(orig_de.index[orig_de["padj"] < 0.05]) if "padj" in orig_de.columns else set()
    sigv = set(res.index[res["padj"] < 0.05])
    # regulator deltas
    cat_rows = []
    for g in REG:
        if g in res.index and g in orig_cat["gene"].values:
            r = res.loc[g]
            oc = orig_cat[orig_cat["gene"] == g].iloc[0]
            cat_rows.append(dict(gene=g, lfc_new=r["log2FoldChange"], padj_new=r["padj"],
                                 lfc_orig=oc["layer1_log2FC"], padj_orig=oc["layer1_padj"]))
    pd.DataFrame(cat_rows).to_csv(os.path.join(OUT, f"p1_07_variant_regulators_{variant}.tsv"), sep="\t", index=False)
    # direction flips among orig-significant regulators
    flips = 0
    for c in cat_rows:
        if pd.notna(c["padj_orig"]) and c["padj_orig"] < 0.05 and pd.notna(c["padj_new"]):
            if np.sign(c["lfc_new"]) != np.sign(c["lfc_orig"]):
                flips += 1
    row = dict(variant=variant, n_cells_used=int(mask.sum()),
               removed=int((~mask).sum()),
               astro_cells=int(sub.values.sum()),
               astro_subpop_reactive=int(sub["reactive"].sum()) if "reactive" in sub.columns else 0,
               DE_sig=n_sig, DE_sig_overlap_with_A=len(sig0 & sigv) if sig0 else -1,
               regulator_direction_flips=flips)
    summary_rows.append(row)
    print(f"[{variant}] astro={row['astro_cells']} reactive={row['astro_subpop_reactive']} "
          f"DE_sig={n_sig} overlapA={row['DE_sig_overlap_with_A']} flips={flips} {time.time()-t0:.0f}s", flush=True)

for name, mask in variants.items():
    pipeline(name, mask)

pd.DataFrame(summary_rows).to_csv(os.path.join(OUT, "p1_07_sensitivity_summary.tsv"), sep="\t", index=False)
print("\nP1-07 DONE", flush=True)
