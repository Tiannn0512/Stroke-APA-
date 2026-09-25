# -*- coding: utf-8 -*-
"""Approximate recompute of Atp2a2 in GSE143531 (Mazare et al 2020 PAPome).
Design: 3 animals (C1-C3) x 2 compartments (Full Astrocytes GFP-IP vs PAP GFP-IP)
x 8 technical replicates. We sum tech reps per animal, apply median-of-ratios
size factors, log2(norm+1), paired Welch t-test PAP vs Full within animal,
BH FDR. Approximation of the authors' Eoulsan/DESeq2-style analysis."""
import gzip, numpy as np
from scipy import stats
from pathlib import Path

CACHE = Path("D:/stroke_apa_reaudit/results/07_followup_20260925/notes/gse_cache/GSE143531")
groups = {"2017545":"C1_Full","2017547":"C1_PAP","2017553":"C2_Full",
          "2017555":"C2_PAP","2017561":"C3_Full","2017563":"C3_PAP"}
import glob
samples = {}
order = None
for g, lab in groups.items():
    acc = []
    for f in sorted(glob.glob(str(CACHE/f"*expression_{g}[a-h].tsv.gz"))):
        ids, vals = [], []
        with gzip.open(f, "rt") as fh:
            next(fh)
            for line in fh:
                a, b = line.rstrip("\n").split("\t")
                ids.append(a); vals.append(int(b))
        if order is None: order = ids
        acc.append(np.array(vals, float))
    samples[lab] = np.sum(acc, axis=0)  # sum 8 tech reps
    print(lab, "files:", len(acc), "summed Atp2a2 count:", int(samples[lab][order.index("ENSMUSG00000029467")]))

lab_order = ["C1_Full","C1_PAP","C2_Full","C2_PAP","C3_Full","C3_PAP"]
mat = np.column_stack([samples[l] for l in lab_order])
names = np.array(order)
i = int(np.where(names == "ENSMUSG00000029467")[0][0])
pos = mat > 0
gm = np.exp(np.mean(np.log(np.where(pos, mat, np.nan)), axis=1))
valid = np.isfinite(gm)
sf = np.median(mat[valid] / gm[valid, None], axis=0)
norm = mat / sf
print("size factors:", dict(zip(lab_order, np.round(sf, 4))))
lg = np.log2(norm + 1)
full_idx = [0, 2, 4]; pap_idx = [1, 3, 5]
lfc = lg[:, pap_idx].mean(1) - lg[:, full_idx].mean(1)
# paired t-test across animals
ps = np.ones(mat.shape[0])
for g in range(mat.shape[0]):
    d = lg[g, pap_idx] - lg[g, full_idx]
    if d.std(ddof=1) == 0: continue
    ps[g] = stats.ttest_rel(lg[g, pap_idx], lg[g, full_idx]).pvalue
def bh(p):
    p = np.asarray(p); n = len(p)
    o = np.argsort(p); qs = p[o]*n/np.arange(1, n+1)
    qs = np.minimum.accumulate(qs[::-1])[::-1]
    q = np.empty(n); q[o] = np.minimum(qs, 1.0); return q
padj = bh(ps)
print("\n== contrast PAP GFP-IP vs Full Astrocyte GFP-IP (summed tech reps, n=3 animals, paired t, BH) ==")
print(f"Atp2a2 ({names[i]}) normalized counts Full: {np.round(norm[i, full_idx],1)}  PAP: {np.round(norm[i, pap_idx],1)}")
print(f"Atp2a2 log2FC = {lfc[i]:.4f}   p(paired) = {ps[i]:.4f}   padj(BH) = {padj[i]:.4f}")
print(f"rank by padj: {int((padj < padj[i]).sum())+1} of {mat.shape[0]}")
print(f"genes padj<0.05: {(padj<0.05).sum()}")
# unpaired version for reference
ps_u = np.ones(mat.shape[0])
for g in range(mat.shape[0]):
    a, b = lg[g, pap_idx], lg[g, full_idx]
    if a.std(ddof=1)==0 and b.std(ddof=1)==0: continue
    ps_u[g] = stats.ttest_ind(a, b, equal_var=False).pvalue
padj_u = bh(ps_u)
print(f"[unpaired Welch reference] Atp2a2 p = {ps_u[i]:.4f}  padj = {padj_u[i]:.4f}")
