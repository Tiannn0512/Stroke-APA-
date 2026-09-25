# -*- coding: utf-8 -*-
"""Approximate recompute of Atp2a2 enrichment in GSE74456 (PAPTRAP).
Method: DESeq2-style median-of-ratios size factors, log2(norm count+1),
Welch t-test per gene across n=3 vs n=3, BH FDR over all genes with
positive counts in the contrast. NOT identical to the authors' model.
"""
import gzip, numpy as np
from scipy import stats
from pathlib import Path

CACHE = Path("D:/stroke_apa_reaudit/results/07_followup_20260925/notes/gse_cache/GSE74456")
files = {
 "CtxPre1":"GSM1920988_CtxPre1.count.txt.gz","CtxPre2":"GSM1920989_CtxPre2.count.txt.gz",
 "CtxPre3":"GSM1920990_CtxPre3.count.txt.gz","CtxTRAP1":"GSM1920991_CtxTRAP1.count.txt.gz",
 "CtxTRAP2":"GSM1920992_CtxTRAP2.count.txt.gz","CtxTRAP3":"GSM1920993_CtxTRAP3.count.txt.gz",
 "PAPpre1":"GSM1920994_PAPpre1.count.txt.gz","PAPpre2":"GSM1920995_PAPpre2.count.txt.gz",
 "PAPpre3":"GSM1920996_PAPpre3.count.txt.gz","PAPTRAP1":"GSM1920997_PAPTRAP1.count.txt.gz",
 "PAPTRAP2":"GSM1920998_PAPTRAP2.count.txt.gz","PAPTRAP3":"GSM1920999_PAPTRAP3.count.txt.gz"}

counts, order = {}, None
for s, f in files.items():
    ids, vals = [], []
    with gzip.open(CACHE/f, "rt") as fh:
        for line in fh:
            a, b = line.rstrip("\n").split("\t")
            if a.startswith("__"): continue
            ids.append(a); vals.append(int(b))
    if order is None: order = ids
    assert ids == order, s
    counts[s] = np.array(vals, float)
mat = np.column_stack([counts[s] for s in files])   # genes x samples
names = np.array(order)
ATP2A2 = "ENSMUSG00000029467"  # Atp2a2, NCBI GeneID 11938 (NCBI eutils efetch 2026-09-25)
i = int(np.where(names == ATP2A2)[0][0])

# median-of-ratios size factors (only genes with all-nonzero proxy: use positive geometric mean)
pos = mat > 0
gm = np.exp(np.mean(np.log(np.where(pos, mat, np.nan)), axis=1))
valid = np.isfinite(gm)
sf = np.median(mat[valid] / gm[valid, None], axis=0)
norm = mat / sf
print("size factors:", dict(zip(files, np.round(sf, 4))))
print("library totals:", {s: int(mat[:, k].sum()) for k, s in enumerate(files)})

def bh(p):
    p = np.asarray(p); n = len(p)
    o = np.argsort(p); ps = p[o]
    qs = ps * n / np.arange(1, n+1)
    qs = np.minimum.accumulate(qs[::-1])[::-1]
    q = np.empty(n); q[o] = np.minimum(qs, 1.0); return q

groups = {"PAP": (["PAPTRAP1","PAPTRAP2","PAPTRAP3"], ["PAPpre1","PAPpre2","PAPpre3"]),
          "Ctx": (["CtxTRAP1","CtxTRAP2","CtxTRAP3"], ["CtxPre1","CtxPre2","CtxPre3"])}
for comp, (traps, pres) in groups.items():
    t_idx = [list(files).index(s) for s in traps]; p_idx = [list(files).index(s) for s in pres]
    lg = np.log2(norm + 1)
    lfc = lg[:, t_idx].mean(1) - lg[:, p_idx].mean(1)
    ts, ps = [], []
    for g in range(mat.shape[0]):
        a, b = lg[g, t_idx], lg[g, p_idx]
        if a.std(ddof=1) == 0 and b.std(ddof=1) == 0:
            ts.append(0.0); ps.append(1.0); continue
        r = stats.ttest_ind(a, b, equal_var=False)
        ts.append(r.statistic); ps.append(r.pvalue)
    padj = bh(np.array(ps))
    tested = (mat[:, t_idx].sum(1) + mat[:, p_idx].sum(1)) > 0
    print(f"\n== contrast {comp} TRAP vs {comp} pre (n=3 vs 3, Welch t on log2(norm+1), BH over {tested.sum()} tested genes) ==")
    print(f"Atp2a2 ({ATP2A2}) raw counts pre : {[int(mat[i, j]) for j in p_idx]}")
    print(f"Atp2a2 ({ATP2A2}) raw counts TRAP: {[int(mat[i, j]) for j in t_idx]}")
    print(f"Atp2a2 log2(norm+1) pre : {np.round(lg[i, p_idx], 3)}")
    print(f"Atp2a2 log2(norm+1) TRAP: {np.round(lg[i, t_idx], 3)}")
    print(f"Atp2a2 mean log2FC = {lfc[i]:.4f}   p = {ps[i]:.3e}   padj(BH) = {padj[i]:.3e}")
    print(f"rank by padj: {int((padj[tested] < padj[i]).sum())+1} of {tested.sum()}")
    # also simple mean-ratio version (DESeq2-like baseMean ratio on normalized counts)
    ratio = (norm[i, t_idx].mean() + 1) / (norm[i, p_idx].mean() + 1)
    print(f"Atp2a2 mean normalized counts: TRAP {norm[i, t_idx].mean():.1f} vs pre {norm[i, p_idx].mean():.1f}; ratio log2 = {np.log2(ratio):.4f}")

# Additional contrast: PAP TRAP vs Ctx TRAP (PAP-enriched translation, both GFP-IP)
lg = np.log2(norm + 1)
pap_t = [list(files).index(s) for s in ["PAPTRAP1","PAPTRAP2","PAPTRAP3"]]
ctx_t = [list(files).index(s) for s in ["CtxTRAP1","CtxTRAP2","CtxTRAP3"]]
lfc2 = lg[:, pap_t].mean(1) - lg[:, ctx_t].mean(1)
ts2, ps2 = [], []
for g in range(mat.shape[0]):
    a, b = lg[g, pap_t], lg[g, ctx_t]
    r = stats.ttest_ind(a, b, equal_var=False)
    ts2.append(r.statistic); ps2.append(r.pvalue)
padj2 = bh(np.array(ps2))
print("\n== contrast PAP TRAP vs Ctx TRAP (n=3 vs 3, Welch t on log2(norm+1), BH over all genes) ==")
print(f"Atp2a2 log2FC = {lfc2[i]:.4f}   p = {ps2[i]:.3e}   padj(BH) = {padj2[i]:.3e}")
print(f"rank by padj: {int((padj2 < padj2[i]).sum())+1} of {mat.shape[0]}")
print(f"Atp2a2 log2(norm+1) PAPTRAP: {np.round(lg[i, pap_t],3)}  CtxTRAP: {np.round(lg[i, ctx_t],3)}")
