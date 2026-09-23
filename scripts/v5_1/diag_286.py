"""Sanity check: why GSE286075 padj all ~1.0? Look at DE stats distribution + counts scale."""
import pandas as pd, gzip, os
res = pd.read_csv(r"D:\stroke_apa\results\p1_DE_paired_GSE286075.tsv", sep="\t")
print("padj<0.05:", (res['padj']<0.05).sum(), "/", len(res))
print("padj describe:\n", res['padj'].describe())
print("pvalue<0.05:", (res['pvalue']<0.05).sum())
print("top by pvalue:")
print(res.nsmallest(8, 'pvalue')[['baseMean','log2FoldChange','pvalue','padj']].to_string())
# raw counts scale for a highly expressed gene
fn = r"D:\stroke_apa\data\GSE286075\extract\GSM8717417_Stroke1ReadsPerGene.out.tab.gz"
with gzip.open(fn, 'rt') as f:
    recs = [l.split('\t') for l in f.read().splitlines() if not l.startswith('N_')]
tot = sum(int(r[2]) for r in recs)
print("\nStroke1 total N_unstranded reads:", tot)
