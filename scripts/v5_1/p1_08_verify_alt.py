# -*- coding: utf-8 -*-
"""P1-08 verify: does pydeseq2 0.5.4 alt_hypothesis actually change padj?
Shared fitted dds; three DeseqStats variants on paper prefilter (col4)."""
import gzip, os
import numpy as np, pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

EXTRACT = r"D:\stroke_apa\data\GSE286075\extract"
files = {
    "CTR1": "GSM8717414_CTR1ReadsPerGene.out.tab.gz", "CTR2": "GSM8717415_CTR2ReadsPerGene.out.tab.gz",
    "CTR3": "GSM8717416_CTR3ReadsPerGene.out.tab.gz", "Stroke1": "GSM8717417_Stroke1ReadsPerGene.out.tab.gz",
    "Stroke2": "GSM8717418_Stroke2ReadsPerGene.out.tab.gz", "Stroke3": "GSM8717419_Stroke3ReadsPerGene.out.tab.gz",
}
series, gids = {}, None
for name, fn in files.items():
    with gzip.open(os.path.join(EXTRACT, fn), "rt") as f:
        rec = [ln.rstrip("\n").split("\t") for ln in f]
    rec = [r for r in rec if not r[0].startswith("N_")]
    if gids is None:
        gids = [r[0] for r in rec]
    series[name] = pd.Series([int(r[3]) for r in rec], index=gids)
cnt = pd.DataFrame(series).T
cnt = cnt.loc[:, (cnt >= 10).sum(axis=0) >= 3].astype(int)
meta = pd.DataFrame(index=cnt.index)
meta["condition"] = ["CTR" if s.startswith("CTR") else "Stroke" for s in cnt.index]
meta["mouse"] = [s[-1] for s in cnt.index]
dds = DeseqDataSet(counts=cnt, metadata=meta, design="~mouse+condition", quiet=True)
dds.deseq2()

out = {}
for tag, kw in [("plain", {}), ("greaterAbs", dict(lfc_null=0.322, alt_hypothesis="greaterAbs")),
                ("lessAbs", dict(lfc_null=0.322, alt_hypothesis="lessAbs"))]:
    sts = DeseqStats(dds, contrast=["condition", "Stroke", "CTR"], quiet=True, **kw)
    sts.summary()
    r = sts.results_df
    sig = r[(r["padj"] < 0.05)]
    out[tag] = r
    print(f"[{tag}] sig={len(sig)} up={int((sig['log2FoldChange']>0).sum())} "
          f"down={int((sig['log2FoldChange']<0).sum())} "
          f"min_pvalue={r['pvalue'].min():.3g}", flush=True)

common = out["plain"].index
print("pvalue identical plain vs greaterAbs:", bool(np.allclose(out['plain']['pvalue'].values, out['greaterAbs']['pvalue'].values, equal_nan=True)), flush=True)
print("pvalue identical plain vs lessAbs:", bool(np.allclose(out['plain']['pvalue'].values, out['lessAbs']['pvalue'].values, equal_nan=True)), flush=True)
g = out["greaterAbs"]; p = out["plain"]
s_g = set(g.index[g["padj"] < 0.05]); s_p = set(p.index[p["padj"] < 0.05])
print(f"greaterAbs sig={len(s_g)}, overlap with plain({len(s_p)})={len(s_g & s_p)}", flush=True)
