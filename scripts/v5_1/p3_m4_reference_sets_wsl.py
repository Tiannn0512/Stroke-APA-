"""P3/P4 workpackage 7: build M4 endfoot reference sets 1-3.

Set 1: GSE74456 PAPTRAP vs CtxTRAP (PAP-enriched set, pydeseq2, n=3+3)
Set 2: GSE286075 endfoot paired DE (corrected dUTP 2nd-strand version, 409 padj<0.05)
Set 3: GSE263986 stroke-responsive astrocyte (RiboTag IP) per zone, cortex & white matter
       Column map (decoded 2026-09-17 via snRNA enrichment + series design):
         L{1-5}=cortex stroke zones IP, Lb=cortex zones1-3 IP, Lip=cortex zones1-3 input,
         LSb=cortex uninjured IP, LSip=cortex uninjured input; M* = white matter mirror.
Outputs under results/p3_m4_*.tsv / *.txt
"""
import os, glob, gc
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

OUT = "/mnt/d/stroke_apa/results"
REFMAP = "/mnt/d/stroke_apa_data/ref_map_vM25.tsv"
refmap = {}
with open(REFMAP, encoding="utf-8") as f:
    for ln in f:
        a, b = ln.rstrip("\n").split("\t")
        refmap[a] = b
print("refmap:", len(refmap))

def syms(ids):
    return [refmap.get(i.split(".")[0], i) for i in ids]

# ---------------- Set 1: GSE74456 PAPTRAP vs CtxTRAP ----------------
files = {
    "CtxTRAP1": "GSM1920991_CtxTRAP1.count.txt.gz", "CtxTRAP2": "GSM1920992_CtxTRAP2.count.txt.gz",
    "CtxTRAP3": "GSM1920993_CtxTRAP3.count.txt.gz", "PAPTRAP1": "GSM1920997_PAPTRAP1.count.txt.gz",
    "PAPTRAP2": "GSM1920998_PAPTRAP2.count.txt.gz", "PAPTRAP3": "GSM1920999_PAPTRAP3.count.txt.gz",
}
series, gids = {}, None
for name, fn in files.items():
    tab = pd.read_csv(os.path.join("/mnt/d/stroke_apa/data/GSE74456", fn), sep="\t", header=None,
                      names=["gid", name])
    if gids is None:
        gids = tab["gid"]
    series[name] = tab.set_index("gid")[name]
cnt = pd.DataFrame(series)
cnt = cnt.loc[:, cnt.sum(axis=0) >= 10].astype(int)
meta = pd.DataFrame(index=cnt.index)
meta["condition"] = ["CtxTRAP" if s.startswith("Ctx") else "PAPTRAP" for s in cnt.index]
dds = DeseqDataSet(counts=cnt, metadata=meta, design="~condition", quiet=True)
dds.deseq2()
sts = DeseqStats(dds, contrast=["condition", "PAPTRAP", "CtxTRAP"], quiet=True)
sts.summary()
res = sts.results_df.copy()
del dds, sts
gc.collect()
res["symbol"] = syms(res.index.tolist())
res.to_csv(os.path.join(OUT, "p3_m4_set1_paptrap_de.tsv"), sep="\t")
pap = res[(res["padj"] < 0.05) & (res["log2FoldChange"] > 0) & (res["baseMean"] >= 10)]
pap["symbol"].drop_duplicates().to_csv(os.path.join(OUT, "p3_m4_set1_pap_enriched.txt"), index=False, header=False)
print(f"[SET1] tested={len(res)} PAP-enriched(padj<0.05,l2fc>0)={len(pap)}")

# ---------------- Set 2: GSE286075 endfoot (corrected DE) ----------------
de2 = pd.read_csv(os.path.join(OUT, "p1_DE_paired_GSE286075.tsv"), sep="\t")
sig2 = de2[de2["padj"] < 0.05]
de2["symbol"].dropna().drop_duplicates().to_csv(os.path.join(OUT, "p3_m4_set2_endfoot_background.txt"), index=False, header=False)
sig2["symbol"].dropna().drop_duplicates().to_csv(os.path.join(OUT, "p3_m4_set2_endfoot_de.txt"), index=False, header=False)
print(f"[SET2] background={de2['symbol'].nunique()} endfoot DE padj<0.05={len(sig2)}")

# ---------------- Set 3: GSE263986 stroke-responsive (IP only) ----------------
xl = pd.read_excel("/mnt/d/stroke_apa/data/GSE263986/GSE263986_ProcessedData.xlsx", sheet_name="count")
xl = xl.set_index("GeneID")
xl = xl.loc[:, xl.sum(axis=0) > 0]  # drop empty columns (WM uninjured input has 4 not 5)

def run_zone_contrasts(prefix, uninjured_col, tissue):
    """Per-zone stroke IP (n=5) vs uninjured IP (n=5)."""
    union = set()
    zone_rows = []
    for z in range(1, 6):
        cols = [f"{prefix}{z}_{i}" for i in range(1, 6) if f"{prefix}{z}_{i}" in xl.columns]
        ucols = [f"{uninjured_col}_{i}" for i in range(1, 6) if f"{uninjured_col}_{i}" in xl.columns]
        sub = xl[cols + ucols]
        sub = sub.loc[sub.sum(axis=1) >= 10]
        m = pd.DataFrame(index=sub.index)
        m["condition"] = ["stroke"] * len(cols) + ["uninjured"] * len(ucols)
        dds = DeseqDataSet(counts=sub.astype(int), metadata=m, design="~condition", quiet=True)
        dds.deseq2()
        st = DeseqStats(dds, contrast=["condition", "stroke", "uninjured"], quiet=True)
        st.summary()
        r = st.results_df.copy()
        del dds, st
        gc.collect()
        r["zone"] = z
        r["symbol"] = syms(r.index.tolist())
        zone_rows.append(r)
        hit = r[(r["padj"] < 0.05)]
        union |= set(hit["symbol"])
        print(f"  {tissue} zone{z}: tested={len(r)} sig={len(hit)}")
    allr = pd.concat(zone_rows)
    allr.to_csv(os.path.join(OUT, f"p3_m4_set3_zone_de_{tissue}.tsv"), sep="\t")
    pd.Series(sorted(union)).to_csv(os.path.join(OUT, f"p3_m4_set3_stroke_responsive_{tissue}.txt"),
                                    index=False, header=False)
    print(f"[SET3:{tissue}] stroke-responsive union={len(union)}")

print("[SET3] cortex zones (L*) vs LSb")
run_zone_contrasts("L", "LSb", "cortex")
print("[SET3] white matter zones (M*) vs MSb")
run_zone_contrasts("M", "MSb", "whitematter")
print("[DONE] M4 reference sets built")
