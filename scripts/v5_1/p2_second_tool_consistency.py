# -*- coding: utf-8 -*-
"""G2(d) second-tool consistency WITHOUT qapa annotation (fallback, registered).
Gene-level long-isoform fraction from salmon TPM (transcripts split by within-gene
median 3'UTR length); dLongFrac per contrast vs DaPars2 gene-mean dPDUI.
Expect POSITIVE correlation. Outputs: p2_second_tool_consistency.tsv + top40 signs.
Run in WSL: dapars2 env python.
"""
import re
import numpy as np, pandas as pd
from scipy.stats import spearmanr

GTF = "/home/taylor/reference_v2/gencode.vM25.annotation.gtf"
SALMON_DIR = "/home/taylor/stroke_APA_work/salmon"
OUT = "/mnt/d/stroke_apa/results"
SAMPLES = ["sham1","sham2","day1_rep1","day1_rep2","day3_rep1","day3_rep2",
           "day7_rep1","day7_rep2","day21_rep1","day21_rep2","day60_rep1","day60_rep2"]
TPS = ["day1","day3","day7","day21","day60"]

# ---- 1. GTF: tx -> (gene_name, utr3_len) ----
tx_gene, tx_strand = {}, {}
tx_3end, cds_edge = {}, {}
ex_edge = {}
for_needle = {}
for ln in open(GTF):
    if ln[0] == "#":
        continue
    p = ln.split("\t")
    if len(p) < 9:
        continue
    feat, strand, start, end = p[2], p[6], int(p[3]), int(p[4])
    attr = p[8]
    if feat == "transcript":
        m = re.search(r'transcript_id "([^"]+)"', attr)
        g = re.search(r'gene_name "([^"]+)"', attr)
        if m and g:
            tid = m.group(1)
            tx_gene[tid] = g.group(1)
            tx_strand[tid] = strand
            tx_3end[tid] = end if strand == "+" else -start
            cds_edge[tid] = None
            ex_edge[tid] = None
    elif feat == "exon":
        m = re.search(r'transcript_id "([^"]+)"', attr)
        if m and m.group(1) in tx_gene:
            tid = m.group(1)
            if strand == "+":
                ex_edge[tid] = end if ex_edge[tid] is None else max(ex_edge[tid], end)
            else:
                ex_edge[tid] = start if ex_edge[tid] is None else min(ex_edge[tid], start)
    elif feat == "CDS":
        m = re.search(r'transcript_id "([^"]+)"', attr)
        if m and m.group(1) in tx_gene:
            tid = m.group(1)
            if strand == "+":
                cds_edge[tid] = end if cds_edge[tid] is None else max(cds_edge[tid], end)
            else:
                cds_edge[tid] = start if cds_edge[tid] is None else min(cds_edge[tid], start)

tx_rows = []
for tid, g in tx_gene.items():
    ce = cds_edge.get(tid)
    if ce is None:
        continue
    strand = tx_strand[tid]
    e3 = ex_edge.get(tid)
    if e3 is None:
        continue
    utr = (e3 - ce) if strand == "+" else (ce - e3)
    tx_rows.append((tid, g, max(utr, 0)))
txmap_gene = dict((t, g) for t, g, _ in tx_rows)
txmap_utr = dict((t, u) for t, _, u in tx_rows)
print(f"GTF parsed: {len(tx_rows)} CDS+ transcripts", flush=True)

# ---- 2. per-sample gene long-isoform fraction ----
def gene_longfrac(sample):
    d = pd.read_csv(f"{SALMON_DIR}/{sample}/quant.sf", sep="\t", usecols=["Name", "TPM"])
    d["gene"] = d["Name"].map(txmap_gene)
    d["utr3"] = d["Name"].map(txmap_utr)
    d = d.dropna(subset=["gene"])
    out = {}
    for gene, g in d.groupby("gene", sort=False):
        if g["utr3"].nunique() < 2:
            continue
        tot = g["TPM"].sum()
        if tot < 1:
            continue
        med = g["utr3"].median()
        out[gene] = g[g["utr3"] > med]["TPM"].sum() / tot
    return pd.Series(out, name=sample)

LFS = pd.concat([gene_longfrac(s) for s in SAMPLES], axis=1)
LFS.columns = SAMPLES
print(f"genes with long/short isoforms: {len(LFS)}", flush=True)

# ---- 3. DaPars2 gene-level dPDUI ----
mat = pd.read_csv(f"{OUT}/p2_full_pdui_matrix.tsv", sep="\t")
mat["symbol"] = mat["Gene"].astype(str).str.split("|").str[2]
sham = mat[["sham1", "sham2"]].mean(axis=1)
dPD = {}
for tp in TPS:
    d = mat[[f"{tp}_rep1", f"{tp}_rep2"]].mean(axis=1) - sham
    d.index = mat["symbol"]
    dPD[tp] = d.groupby(level=0).mean()

# ---- 4. consistency across all testable genes ----
res, lines = [], []
for tp in TPS:
    dlf = (LFS[[f"{tp}_rep1", f"{tp}_rep2"]].mean(axis=1) - LFS[["sham1", "sham2"]].mean(axis=1)).dropna()
    j = pd.concat([dlf.rename("dLongFrac"), dPD[tp].rename("dPDUI")], axis=1).dropna()
    rho, pv = spearmanr(j["dLongFrac"], j["dPDUI"])
    res.append(dict(contrast=tp, n_genes=len(j), spearman_rho=round(rho, 3), p_value=f"{pv:.2e}"))
    lines.append(f"{tp}\tn={len(j)}\trho={rho:.3f}\tp={pv:.2e}")
res = pd.DataFrame(res)
res.to_csv(f"{OUT}/p2_second_tool_consistency.tsv", sep="\t", index=False)
print(res.to_string(index=False), flush=True)

# ---- 5. top40 sign concordance ----
top = pd.read_csv(f"{OUT}/p2_gate_g2_top_events.tsv", sep="\t")
mat1 = mat.copy(); mat1.index = mat1.index + 1  # layer2 event ids are 1-based
top["symbol"] = top["event"].astype(int).map(mat1["Gene"]).astype(str).str.split("|").str[2]
agree = tot_n = 0
detail = []
for _, r in top.iterrows():
    s, tp = r["symbol"], r["contrast"]
    if s not in LFS.index:
        continue
    dlf = LFS.loc[s, [f"{tp}_rep1", f"{tp}_rep2"]].mean() - LFS.loc[s, ["sham1", "sham2"]].mean()
    if pd.isna(dlf):
        continue
    tot_n += 1
    ok = np.sign(dlf) == np.sign(r["dPDUI"])
    agree += int(ok)
    detail.append((s, tp, round(float(dlf), 3), r["dPDUI"], "agree" if ok else "DISAGREE"))
print(f"\n[TOP40 sign concordance] {agree}/{tot_n}", flush=True)
with open(f"{OUT}/p2_second_tool_top40_signs.tsv", "w") as f:
    f.write("symbol\tcontrast\tdLongFrac\tdPDUI\tcall\n")
    for a, b, c, dd, e in detail:
        f.write(f"{a}\t{b}\t{c}\t{dd}\t{e}\n")
print("[DONE] second-tool consistency", flush=True)
