# -*- coding: utf-8 -*-
"""P1-08: GSE286075 205-vs-409 reproduction appendix (audit 14.4).
Original paper (Shim 2025 IJMS, PMC11720067) methods: DESeq2 v1.44.0, prefilter
'>=10 counts in >=3 samples' -> 10,862 genes, paired design (animal + condition),
ashr-shrunk LFC tested against interval null |log2FC| < 0.322 via s-values -> 205 DE (86 up/119 down).
Our P1-05 run: col4, sum>=10 loose filter -> 14,323 genes, ~mouse+condition, plain padj<0.05 -> 409.
Grid below isolates: (a) test criterion, (b) prefilter, (c) strand column.
"""
import gzip, os
import numpy as np, pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

EXTRACT = r"D:\stroke_apa\data\GSE286075\extract"
OUT = r"D:\stroke_apa\results"
LFC_NULL = 0.322  # paper: |log2FC| < 0.322 <=> <4/5 / >5/4 of control

files = {
    "CTR1": "GSM8717414_CTR1ReadsPerGene.out.tab.gz", "CTR2": "GSM8717415_CTR2ReadsPerGene.out.tab.gz",
    "CTR3": "GSM8717416_CTR3ReadsPerGene.out.tab.gz", "Stroke1": "GSM8717417_Stroke1ReadsPerGene.out.tab.gz",
    "Stroke2": "GSM8717418_Stroke2ReadsPerGene.out.tab.gz", "Stroke3": "GSM8717419_Stroke3ReadsPerGene.out.tab.gz",
}

def load(col):
    series, gids = {}, None
    for name, fn in files.items():
        with gzip.open(os.path.join(EXTRACT, fn), "rt") as f:
            rec = [ln.rstrip("\n").split("\t") for ln in f]
        rec = [r for r in rec if not r[0].startswith("N_")]
        if gids is None:
            gids = [r[0] for r in rec]
        else:
            assert [r[0] for r in rec] == gids
        series[name] = pd.Series([int(r[col]) for r in rec], index=gids)
    return pd.DataFrame(series).T

def run(cnt, label, alt=None):
    cnt = cnt.loc[:, cnt.sum(axis=0) >= 10].astype(int)
    meta = pd.DataFrame(index=cnt.index)
    meta["condition"] = ["CTR" if s.startswith("CTR") else "Stroke" for s in cnt.index]
    meta["mouse"] = [s[-1] for s in cnt.index]
    kw = dict(lfc_null=LFC_NULL, alt_hypothesis=alt) if alt else {}
    dds = DeseqDataSet(counts=cnt, metadata=meta, design="~mouse+condition", quiet=True)
    dds.deseq2()
    sts = DeseqStats(dds, contrast=["condition", "Stroke", "CTR"], quiet=True, **kw)
    sts.summary()
    res = sts.results_df
    n_padj = int((res["padj"] < 0.05).sum())
    up = int(((res["padj"] < 0.05) & (res["log2FoldChange"] > 0)).sum())
    dn = int(((res["padj"] < 0.05) & (res["log2FoldChange"] < 0)).sum())
    print(f"[{label}] genes_in_model={res.shape[0]} sig={n_padj} (up={up} down={dn})", flush=True)
    return res, n_padj, up, dn

rows = []

# V0: our original P1-05 run (col4, loose filter, plain padj) -- expect 409
res0, n, u, d = run(load(3), "V0 ours col4 loose plain")
rows.append(dict(variant="V0_ours_col4_loose_plain", genes=res0.shape[0], sig=n, up=u, down=d))

# V1: paper criteria on our data (col4, paper prefilter, interval-null |LFC|>0.322 test)
cnt4 = load(3)
mask = (cnt4 >= 10).sum(axis=0) >= 3
cnt_paper = cnt4.loc[:, mask]
print(f"paper-style prefilter genes = {cnt_paper.shape[1]} (paper: 10,862)", flush=True)
res1, n, u, d = run(cnt_paper, "V1 col4 paperfilter greaterAbs0.322", alt="greaterAbs")
rows.append(dict(variant="V1_col4_paperfilter_lfc0.322", genes=res1.shape[0], sig=n, up=u, down=d))

# V2: paper prefilter + plain padj (isolates test criterion)
res2, n, u, d = run(cnt_paper, "V2 col4 paperfilter plain")
rows.append(dict(variant="V2_col4_paperfilter_plain", genes=res2.shape[0], sig=n, up=u, down=d))

# V3: loose prefilter + interval-null test (isolates prefilter)
res3, n, u, d = run(cnt4, "V3 col4 loose greaterAbs0.322", alt="greaterAbs")
rows.append(dict(variant="V3_col4_loose_lfc0.322", genes=res3.shape[0], sig=n, up=u, down=d))

# V4: col2 (unstranded) with paper criteria (strand sensitivity)
cnt2 = load(1)
cnt2p = cnt2.loc[:, (cnt2 >= 10).sum(axis=0) >= 3]
res4, n, u, d = run(cnt2p, "V4 col2 paperfilter plain")
rows.append(dict(variant="V4_col2_paperfilter_plain", genes=res4.shape[0], sig=n, up=u, down=d))

grid = pd.DataFrame(rows)
grid.to_csv(os.path.join(OUT, "p1_08_grid_205v409.tsv"), sep="\t", index=False)
print(grid.to_string(index=False), flush=True)

# overlap of V1 (paper-criteria reproduction) with our V0 409 list
s0 = set(res0.index[res0["padj"] < 0.05])
s1 = set(res1.index[res1["padj"] < 0.05])
ov = s0 & s1
print(f"\noverlap V0(409) vs V1({len(s1)}): {len(ov)}", flush=True)

# top-gene anchor: paper abstract names Hspa1a/Hspa1b as top upregulated (HSF1-dependent)
for g in ["Hspa1a", "Hspa1b"]:
    for tag, res in (("V1", res1), ("V0", res0)):
        if g in res.index:
            r = res.loc[g]
            print(f"anchor {g} [{tag}]: log2FC={r['log2FoldChange']:.2f} padj={r['padj']:.2e}", flush=True)

# save V1 full table as the appendix's reproduction layer
res1.to_csv(os.path.join(OUT, "p1_286075_paper_criteria_repro.tsv"), sep="\t")
print("\nP1-08 grid DONE", flush=True)
