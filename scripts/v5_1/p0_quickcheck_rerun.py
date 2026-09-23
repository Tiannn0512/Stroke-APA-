"""P0-08 rerun: regulator quickcheck with correct id->name map from genes.tsv col2."""
import glob, gzip, os, warnings
import numpy as np, pandas as pd, scipy.io as sio
import anndata as ad
warnings.filterwarnings("ignore")

EXTRACT = r"D:\stroke_apa\data\GSE174574\extract"
OUT = r"D:\stroke_apa\results"

samples = [("GSM5319987", "sham1", "sham"), ("GSM5319988", "sham2", "sham"), ("GSM5319989", "sham3", "sham"),
           ("GSM5319990", "MCAO1", "MCAO"), ("GSM5319991", "MCAO2", "MCAO"), ("GSM5319992", "MCAO3", "MCAO")]

REG = ["Nudt21", "Cpsf6", "Cstf2", "Qki", "Elavl1",
       "Aqp4", "Slc1a3", "Gfap", "Aldh1l1", "Serpina3n", "C3"]

# id->name map from sham1 genes.tsv (identical across samples; verified same size 217,843)
with gzip.open(glob.glob(os.path.join(EXTRACT, "GSM5319987_sham1_genes.tsv.gz"))[0], "rt") as f:
    recs = [ln.rstrip("\n").split("\t") for ln in f]
id2name = {r[0]: (r[1] if len(r) > 1 else r[0]) for r in recs}

rows = []
for gsm, sname, cond in samples:
    mm = sio.mmread(glob.glob(os.path.join(EXTRACT, f"{gsm}_{sname}_matrix.mtx.gz"))[0]).tocsr()
    gene_ids = [r[0] for r in recs]
    names = np.array([id2name.get(g, g) for g in gene_ids])
    counts = np.asarray(mm.sum(axis=0)).ravel().astype(float)
    for g in REG:
        hit = np.where(names == g)[0]
        if len(hit) == 0:
            hit = np.where(np.char.lower(names.astype(str)) == g.lower())[0]
        for h in hit:
            vec = np.asarray(mm[h].todense()).ravel().astype(float)
            rows.append(dict(gene=g, ensembl=gene_ids[h], sample=sname, condition=cond,
                             mean_raw=float(vec.mean()), median_raw=float(np.median(vec)),
                             pct_expr=float((vec > 0).mean() * 100),
                             sample_total_counts=float(counts.sum())))
qc = pd.DataFrame(rows)
qc.to_csv(os.path.join(OUT, "p0_regulator_quickcheck.tsv"), sep="\t", index=False)
print("rows:", len(qc))
pv = qc.pivot_table(index="gene", columns="sample", values="mean_raw")
print("=== mean raw counts ===")
print(pv.to_string(float_format=lambda x: f"{x:.3f}"))
pv2 = qc.pivot_table(index="gene", columns="sample", values="pct_expr")
print("\n=== % expressing ===")
print(pv2.to_string(float_format=lambda x: f"{x:.1f}"))
