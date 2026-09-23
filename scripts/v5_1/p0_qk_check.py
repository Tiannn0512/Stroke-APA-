"""P0-08 addendum: Qk (official MGI symbol for quaking) quickcheck."""
import glob, gzip, os, warnings
import numpy as np, pandas as pd, scipy.io as sio
warnings.filterwarnings("ignore")

EXTRACT = r"D:\stroke_apa\data\GSE174574\extract"
OUT = r"D:\stroke_apa\results"

samples = [("GSM5319987", "sham1"), ("GSM5319988", "sham2"), ("GSM5319989", "sham3"),
           ("GSM5319990", "MCAO1"), ("GSM5319991", "MCAO2"), ("GSM5319992", "MCAO3")]

with gzip.open(glob.glob(os.path.join(EXTRACT, "GSM5319987_sham1_genes.tsv.gz"))[0], "rt") as f:
    recs = [ln.rstrip("\n").split("\t") for ln in f]
id2name = {r[0]: (r[1] if len(r) > 1 else r[0]) for r in recs}
names = np.array([id2name.get(r[0], r[0]) for r in recs])

# find any gene whose name starts with Qk (Qk, Qki, Qk5 etc.)
cand = sorted(set(names[np.char.startswith(names.astype(str), "Qk")]))
print("Qk* gene names present:", cand)

rows = []
for gsm, sname in samples:
    mm = sio.mmread(glob.glob(os.path.join(EXTRACT, f"{gsm}_{sname}_matrix.mtx.gz"))[0]).tocsr()
    for g in cand:
        h = np.where(names == g)[0][0]
        vec = np.asarray(mm[h].todense()).ravel().astype(float)
        rows.append(dict(gene=g, sample=sname, mean_raw=float(vec.mean()), pct_expr=float((vec > 0).mean() * 100)))
qc = pd.DataFrame(rows)
qc.to_csv(os.path.join(OUT, "p0_regulator_quickcheck_Qk.tsv"), sep="\t", index=False)
print(qc.pivot_table(index="gene", columns="sample", values="mean_raw").to_string(float_format=lambda x: f"{x:.3f}"))
print(qc.pivot_table(index="gene", columns="sample", values="pct_expr").to_string(float_format=lambda x: f"{x:.1f}"))
