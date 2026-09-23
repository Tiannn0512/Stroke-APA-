"""P1-04/05/06: three DE layers (pydeseq2) + regulator catalog."""
import os, warnings, gzip, glob
import numpy as np, pandas as pd
import anndata as ad
warnings.filterwarnings("ignore")
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

OUT = r"D:\stroke_apa\results"
REG_CORE = ["Cpsf1","Cpsf2","Cpsf3","Cpsf4","Cpsf4l","Wdr33","Fip1l1","Cstf1","Cstf2","Cstf3",
            "Nudt21","Cpsf6","Cpsf7","Pabpn1","Pcf11","Clp1","Sympk"]
REG_RBP = ["Elavl1","Elavl2","Qk","Ptbp1","Ptbp2","Hnrnpa2b1","Srsf1","Nova1","Nova2",
           "Fus","Tardbp","Khdrbs1","Mbnl1","Mbnl2"]
REG = REG_CORE + REG_RBP

# ---------------- load astrocyte annotated ----------------
astro = ad.read_h5ad(r"D:\stroke_apa\data\GSE174574\h5ad\astrocytes_annotated.h5ad")
Xc = astro.layers["counts"]
sym = np.array([str(x) for x in astro.var_names])
# rebuild ensembl->symbol map directly from genes.tsv (var 'ensembl' column was dropped during concat)
import gzip as _gzip
def _load_map():
    with _gzip.open(r"D:\stroke_apa\data\GSE174574\extract\GSM5319987_sham1_genes.tsv.gz", "rt") as f:
        recs = [ln.rstrip("\n").split("\t") for ln in f]
    return recs
_recs = _load_map()
ens = np.array([r[0] for r in _recs])
sym_from_ens = {r[0]: (r[1] if len(r) > 1 else r[0]) for r in _recs}

def make_counts(mask_dict):
    rows, meta = [], []
    for name, mask in mask_dict.items():
        vec = np.asarray(Xc[np.asarray(mask)].sum(axis=0)).ravel()
        rows.append(vec)
        meta.append(dict(sample=name))
    df = pd.DataFrame(np.array(rows), index=[m["sample"] for m in meta], columns=sym)
    return df

obs = astro.obs
mask_overall = {s: (obs["sample"] == s).values for s in obs["sample"].unique()}
mask_homeo = {s: ((obs["sample"] == s) & (obs["astro_subpop"] == "homeostatic")).values for s in obs["sample"].unique()}
df_overall = make_counts(mask_overall)
df_homeo = make_counts(mask_homeo)

def run_deseq(counts_df, design_cols, contrast, out_prefix, keep_min=10):
    cnt = counts_df.loc[:, counts_df.sum(axis=0) >= keep_min].astype(int)
    meta = pd.DataFrame(index=cnt.index)
    for c, vals in design_cols.items():
        meta[c] = vals
    dds = DeseqDataSet(counts=cnt, metadata=meta, design="~" + "+".join(design_cols), quiet=True)
    dds.deseq2()
    ref, alt = contrast
    sts = DeseqStats(dds, contrast=["condition", alt, ref], quiet=True)
    sts.summary()
    res = sts.results_df.copy()
    res.to_csv(os.path.join(OUT, out_prefix), sep="\t")
    return res

cond_overall = ["MCAO" if s.startswith("MCAO") else "sham" for s in df_overall.index]
res_overall = run_deseq(df_overall, {"condition": cond_overall}, ("sham", "MCAO"), "p1_DE_overall_astro.tsv")
cond_homeo = ["MCAO" if s.startswith("MCAO") else "sham" for s in df_homeo.index]
res_homeo = run_deseq(df_homeo, {"condition": cond_homeo}, ("sham", "MCAO"), "p1_DE_homeostatic_astro.tsv")
print("overall DE done:", res_overall.shape, "| homeostatic DE done:", res_homeo.shape, flush=True)

# ---------------- GSE286075 paired DE ----------------
# Column correction 2026-09-16 (diag_286_strand.py): library is stranded 2nd-strand (dUTP) —
# col4 ≈ col2 (unstranded) while col3 ≈ 0 (e.g. Stroke1 col3=90,093 vs col2=21.0M).
# Original run used col3 -> false "extremely low depth" conclusion; use col4 (2nd strand).
files = {
    "CTR1": "GSM8717414_CTR1ReadsPerGene.out.tab.gz", "CTR2": "GSM8717415_CTR2ReadsPerGene.out.tab.gz",
    "CTR3": "GSM8717416_CTR3ReadsPerGene.out.tab.gz", "Stroke1": "GSM8717417_Stroke1ReadsPerGene.out.tab.gz",
    "Stroke2": "GSM8717418_Stroke2ReadsPerGene.out.tab.gz", "Stroke3": "GSM8717419_Stroke3ReadsPerGene.out.tab.gz",
}
series = {}
gids = None
for name, fn in files.items():
    with gzip.open(os.path.join(r"D:\stroke_apa\data\GSE286075\extract", fn), "rt") as f:
        rec = [ln.rstrip("\n").split("\t") for ln in f]
    rec = [r for r in rec if not r[0].startswith("N_")]
    if gids is None:
        gids = [r[0] for r in rec]
    else:
        assert [r[0] for r in rec] == gids
    series[name] = pd.Series([int(r[3]) for r in rec], index=gids)
cnt286 = pd.DataFrame(series).T
cnt286 = cnt286.loc[:, cnt286.sum(axis=0) >= 10].astype(int)
meta286 = pd.DataFrame(index=cnt286.index)
meta286["condition"] = ["CTR" if s.startswith("CTR") else "Stroke" for s in cnt286.index]
meta286["mouse"] = [s[-1] for s in cnt286.index]
dds2 = DeseqDataSet(counts=cnt286, metadata=meta286, design="~mouse+condition", quiet=True)
dds2.deseq2()
sts2 = DeseqStats(dds2, contrast=["condition", "Stroke", "CTR"], quiet=True)
sts2.summary()
res286 = sts2.results_df.copy()
res286.to_csv(os.path.join(OUT, "p1_DE_paired_GSE286075.tsv"), sep="\t")
print("paired GSE286075 DE done:", res286.shape, flush=True)

# map ensembl(versioned) -> symbol for GSE286075 results
res286 = res286.copy()
res286["ensembl_core"] = [g.split(".")[0] for g in res286.index]
res286["symbol"] = [sym_from_ens.get(e, e) for e in res286["ensembl_core"]]

# ---------------- regulator catalog ----------------
def grab(res, gene_col):
    out = {}
    idx = res.index
    for g in REG:
        if gene_col:
            hit = res[res["symbol"] == g]
            if len(hit):
                r = hit.iloc[0]
                out[g] = (r["log2FoldChange"], r["padj"])
                continue
        else:
            if g in idx:
                r = res.loc[g]
                out[g] = (r["log2FoldChange"], r["padj"])
                continue
        out[g] = (np.nan, np.nan)
    return out

g_overall = grab(res_overall, False)
g_homeo = grab(res_homeo, False)
g_286 = grab(res286, True)

rows = []
for g in REG:
    l1, p1 = g_overall[g]; l1h, p1h = g_homeo[g]; l2, p2 = g_286[g]
    dir1 = "up" if (pd.notna(l1) and l1 > 0) else ("down" if pd.notna(l1) and l1 < 0 else "NA")
    dir2 = "up" if (pd.notna(l2) and l2 > 0) else ("down" if pd.notna(l2) and l2 < 0 else "NA")
    consist = "concordant" if (dir1 == dir2 and dir1 != "NA") else ("discordant" if (dir1 != "NA" and dir2 != "NA") else "NA")
    rows.append(dict(gene=g, layer1_log2FC=l1, layer1_padj=p1, layer1_dir=dir1,
                     layer1_homeostatic_log2FC=l1h, layer1_homeostatic_padj=p1h,
                     layer2_paired_log2FC=l2, layer2_padj=p2, layer2_dir=dir2, direction_consistency=consist))
cat = pd.DataFrame(rows)
cat.to_csv(os.path.join(OUT, "p1_regulator_catalog.tsv"), sep="\t", index=False)
print("\n=== regulator catalog ===")
print(cat.round(3).to_string(index=False), flush=True)

n_consist = (cat["direction_consistency"] == "concordant").sum()
n_sig1 = (cat["layer1_padj"] < 0.05).sum()
print(f"\nconcordant: {n_consist}/{len(cat)} | layer1 padj<0.05: {n_sig1}/{len(cat)} | overall DE genes padj<0.05: {(res_overall['padj']<0.05).sum()}")
