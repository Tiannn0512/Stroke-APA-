"""P0-06/08/09 + P1 load: build AnnData from 6 samples, quickcheck regulators, audit cells."""
import gzip, os, re, warnings
import numpy as np, pandas as pd, scipy.io as sio
import anndata as ad
warnings.filterwarnings("ignore")

EXTRACT = r"D:\stroke_apa\data\GSE174574\extract"
OUT = r"D:\stroke_apa\results"
os.makedirs(OUT, exist_ok=True)
os.makedirs(r"D:\stroke_apa\data\GSE174574\h5ad", exist_ok=True)

samples = {
    "GSM5319987": ("sham1", "sham"), "GSM5319988": ("sham2", "sham"), "GSM5319999": ("sham3", "sham"),
    "GSM5319990": ("MCAO1", "MCAO"), "GSM5319991": ("MCAO2", "MCAO"), "GSM5319992": ("MCAO3", "MCAO"),
}
# fix sham3 gsm (typo guard): files are GSM5319989_sham3
samples["GSM5319989"] = samples.pop("GSM5319999")

def read_10x(gsm_prefix):
    """Return (counts_csr, genes_df, barcodes) from the extract dir."""
    m = os.path.join(EXTRACT, f"{gsm_prefix}_*")
    import glob
    mtx = glob.glob(m + "_matrix.mtx.gz")[0]
    genes = glob.glob(m + "_genes.tsv.gz")[0]
    bcs = glob.glob(m + "_barcodes.tsv.gz")[0]
    mm = sio.mmread(mtx).tocsr()          # genes x cells
    with gzip.open(genes, "rt") as f:
        gl = [ln.rstrip("\n").split("\t") for ln in f]
    gene_ids = [r[0] for r in gl]
    gene_names = [r[1] if len(r) > 1 else r[0] for r in gl]
    with gzip.open(bcs, "rt") as f:
        bcs_list = [ln.strip() for ln in f]
    assert mm.shape[1] == len(bcs_list), (mm.shape, len(bcs_list))
    assert mm.shape[0] == len(gene_ids), (mm.shape, len(gene_ids))
    return mm, gene_ids, gene_names, bcs_list

adatas = []
cell_audit_rows = []
qc_raw = {}
for gsm, (sname, cond) in sorted(samples.items()):
    mm, gids, gnames, bcs = read_10x(gsm)
    a = ad.AnnData(X=mm.T.tocsr())        # cells x genes
    a.var_names = pd.Index(gids)
    a.var["gene_ids"] = gids
    a.var["gene_names"] = pd.Index(gnames)
    a.obs_names = [f"{sname}_{b}" for b in bcs]
    a.obs["sample"] = sname
    a.obs["condition"] = cond
    adatas.append(a)

    # quick QC metrics on raw counts (pre-filter)
    mt_mask = np.array([bool(re.match(r"^mt-", str(g), re.I)) for g in gnames])
    counts = np.asarray(mm.sum(axis=0)).ravel().astype(float)
    ngenes = np.asarray((mm > 0).sum(axis=0)).ravel()
    mt_frac = np.asarray(mm[mt_mask].sum(axis=0)).ravel() / np.maximum(counts, 1)
    qc_raw[sname] = dict(total_cells=mm.shape[1],
        median_counts=float(np.median(counts)), median_genes=float(np.median(ngenes)),
        median_mt=float(np.median(mt_frac) * 100))

    # P0-09 pre-filter audit: standard 10x filters (genes>=200, counts>=500, mt<20%)
    keep = (ngenes >= 200) & (counts >= 500) & (mt_frac < 0.20)
    qc_raw[sname]["pass_qc"] = int(keep.sum())

    del a
combined = ad.concat(adatas, axis=0, join="outer", label=None)
# var alignment across samples: gene panels identical (same genes.tsv), so var columns survive concat.
# Rebuild a reliable gene_names column from var_names (gene ids) using the first sample's id->name map.
_id2name = dict(zip(adatas[0].var["gene_ids"], [str(x) for x in adatas[0].var_names if True]))
gnames = pd.Index([_id2name.get(g, str(g)) for g in combined.var_names]).astype(str)
combined.X = combined.X.tocsr()
combined.write_h5ad(r"D:\stroke_apa\data\GSE174574\h5ad\raw_all.h5ad")

audit = pd.DataFrame(qc_raw).T
audit.to_csv(os.path.join(OUT, "p0_cell_audit_raw.tsv"), sep="\t")

# P0-08 regulator quickcheck on raw counts (per-sample mean of raw counts, plus per-cell pct expressing)
REG = ["Nudt21", "Cpsf6", "Cstf2", "Qki", "Elavl1", "Aqp4", "Slc1a3", "Gfap", "Aldh1l1"]
gnames_idx = gnames
rows = []
for g in REG:
    hit = np.where(gnames_idx.str.lower() == g.lower())[0]
    if len(hit) == 0:
        hit = np.where(gnames_idx.str.startswith(g + "."))[0]
    for h in hit:
        vec = np.asarray(combined.X[:, h].todense()).ravel()
        for s in combined.obs["sample"].unique():
            m = (combined.obs["sample"] == s).values
            sub = vec[m]
            rows.append(dict(gene=str(gnames_idx[h]), sample=s,
                mean_raw=float(sub.mean()), pct_expr=float((sub > 0).mean() * 100)))
qc = pd.DataFrame(rows)
qc.to_csv(os.path.join(OUT, "p0_regulator_quickcheck.tsv"), sep="\t", index=False)

print("=== cell audit (raw, per sample) ===")
print(audit.to_string())
print("\n=== regulator quickcheck (pivot: mean raw counts) ===")
pv = qc.pivot_table(index="gene", columns="sample", values="mean_raw")
print(pv.to_string(float_format=lambda x: f"{x:.2f}"))
print("\n=== pct expressing ===")
pv2 = qc.pivot_table(index="gene", columns="sample", values="pct_expr")
print(pv2.to_string(float_format=lambda x: f"{x:.1f}"))
print("\ncombined AnnData:", combined.shape)
