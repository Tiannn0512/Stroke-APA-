"""Identify the 6 significant genes in GSE286075 paired DE + add direction summary to catalog context."""
import pandas as pd
res = pd.read_csv(r"D:\stroke_apa\results\p1_DE_paired_GSE286075.tsv", sep="\t")
recs = {}
import gzip
with gzip.open(r"D:\stroke_apa\data\GSE174574\extract\GSM5319987_sham1_genes.tsv.gz", "rt") as f:
    for ln in f:
        p = ln.rstrip("\n").split("\t")
        recs[p[0]] = p[1] if len(p) > 1 else p[0]
res["symbol"] = [recs.get(g.split(".")[0], g) for g in res["Unnamed: 0"]]
sig = res[res["padj"] < 0.05].sort_values("padj")
print(sig[["symbol", "baseMean", "log2FoldChange", "padj"]].to_string(index=False))
res.to_csv(r"D:\stroke_apa\results\p1_DE_paired_GSE286075.tsv", sep="\t", index=False)
print("\nsaved with symbol column, rows:", len(res))
