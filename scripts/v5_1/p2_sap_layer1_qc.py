# -*- coding: utf-8 -*-
"""SAP layer 1 (v1.1): sample-level global APA QC + technical confounders. Runs in WSL.
Inputs: merged PDUI matrix, fastp json x12 (WSL qc/), STAR per-sample dirs (bam/<s>_/Log.final.out),
        ReadsPerGene.out.tab x12 (gene expression, dUTP second-strand col4).
Outputs (to /mnt/d/stroke_apa/results/):
  p2_full_qc_sample_level.tsv  per-sample PDUI stats + dup/insert/unique-rate confounders
  p2_full_qc_correlation.tsv   12x12 Pearson
  p2_full_qc_pca.tsv           sample PC1-3
  p2_full_qc_stratified.tsv    d1-d60 shortening/lengthening by expression & UTR-length quartiles
Run: /home/taylor/miniconda3/envs/dapars2/bin/python p2_sap_layer1_qc.py
"""
import glob, json, os, re
import numpy as np, pandas as pd

BASE = "/home/taylor/stroke_APA_work"
OUT = "/mnt/d/stroke_apa/results"
MATRIX = os.path.join(OUT, "p2_full_pdui_matrix.tsv")
SAMPLES = ["sham1","sham2","day1_rep1","day1_rep2","day3_rep1","day3_rep2",
           "day7_rep1","day7_rep2","day21_rep1","day21_rep2","day60_rep1","day60_rep2"]

mat = pd.read_csv(MATRIX, sep="\t")
P = mat[SAMPLES].apply(lambda c: pd.to_numeric(c, errors="coerce"))

# ---- per-sample PDUI stats ----
rows = []
for s in SAMPLES:
    v = P[s].dropna()
    lo, hi = int(len(v) * 0.05), int(len(v) * 0.95)
    trim = v.sort_values().iloc[max(1, lo):max(2, hi)]
    rows.append(dict(sample=s, n_valid=len(v), median_PDUI=round(v.median(), 4),
                     trimmed_mean_PDUI=round(trim.mean(), 4)))

# ---- confounders: fastp + STAR ----
for i, s in enumerate(SAMPLES):
    dup = ins = uni = np.nan
    fp = f"{BASE}/qc/{s}_fastp.json"
    if os.path.exists(fp):
        j = json.load(open(fp))
        dup = j.get("duplication", {}).get("rate", np.nan)
        pk = j.get("insert_size", {}).get("peak", np.nan)
        ins = pk if isinstance(pk, (int, float)) else (pk or {}).get("mean", np.nan)
    logs = glob.glob(f"{BASE}/bam/{s}_*/Log.final.out")
    if logs:
        txt = open(logs[0], errors="ignore").read()
        m = re.search(r"Uniquely mapped reads number \|\s+(\d+)", txt)
        t = re.search(r"Number of input reads \|\s+(\d+)", txt)
        if m and t:
            uni = int(m.group(1)) / int(t.group(1))
    rows[i].update(dup_rate=dup, insert_mean=ins, unique_rate=uni)
qc = pd.DataFrame(rows)
qc.to_csv(os.path.join(OUT, "p2_full_qc_sample_level.tsv"), sep="\t", index=False)
print(qc.to_string(index=False), flush=True)

# ---- correlation ----
corr = P.corr()
corr.to_csv(os.path.join(OUT, "p2_full_qc_correlation.tsv"), sep="\t")
print("\ncorr min/max off-diagonal:", round(corr.values[np.triu_indices(12, 1)].min(), 3),
      round(corr.values[np.triu_indices(12, 1)].max(), 3), flush=True)

# ---- PCA (complete-case) ----
ok = P.dropna(axis=0)
X = ok.T
Xc = X - X.mean(axis=0)
sv, u = np.linalg.eigh(np.corrcoef(Xc.values) if False else np.cov(Xc.values))
order = np.argsort(sv)[::-1]
pc = pd.DataFrame(u[:, order][:, :3] * np.sqrt(sv[order][:3]),
                  index=ok.columns, columns=["PC1", "PC2", "PC3"])
pc.to_csv(os.path.join(OUT, "p2_full_qc_pca.tsv"), sep="\t")
print("PC1/PC2 var%:", (sv[order][:2] / sv.sum() * 100).round(1), flush=True)

# ---- gene expression from ReadsPerGene (col4 = second strand, dUTP) ----
expr = {}
for s in SAMPLES:
    f = glob.glob(f"{BASE}/bam/{s}_*/ReadsPerGene.out.tab")
    if not f:
        continue
    d = pd.read_csv(f[0], sep="\t", header=None)
    d = d[~d[0].astype(str).str.startswith("N_")]
    d[s] = d[3].astype(float)
    expr[s] = d.set_index(0)[s]
E = pd.DataFrame(expr)
E = E.loc[:, [s for s in SAMPLES if s in E.columns]]
gene_expr = E.mean(axis=1)  # per gene; NOTE: our STAR GeneCounts keys = gene SYMBOLS (verified empirically)

sym = mat["Gene"].astype(str).str.split("|").str[2]
ev_expr = sym.map(gene_expr)

# ---- stratified shortening/lengthening ----
loci = mat["Loci"].astype(str)
span = loci.str.extract(r":(\d+)-(\d+)").astype(float)
utr_len = span[1] - span[0]
sham_mean = P[["sham1", "sham2"]].mean(axis=1)
strat_rows = []
for tp in ["day1", "day3", "day7", "day21", "day60"]:
    cols = [c for c in SAMPLES if c.startswith(tp)]
    dd = (P[cols].mean(axis=1) - sham_mean).dropna()
    for axis_name, vals in [("expr_q", ev_expr.reindex(dd.index)), ("utr_len_q", utr_len.reindex(dd.index))]:
        good = vals.notna()
        ddv = dd[good]
        q = pd.qcut(vals[good].rank(method="first"), 4, labels=["Q1", "Q2", "Q3", "Q4"])
        for lab in ["Q1", "Q2", "Q3", "Q4"]:
            sub = ddv[q == lab]
            short = int((sub < -0.1).sum()); long_ = int((sub > 0.1).sum())
            strat_rows.append(dict(timepoint=tp, stratum=axis_name, quartile=lab, n=len(sub),
                                   shortening=short, lengthening=long_,
                                   ratio=round(short / max(long_, 1), 2),
                                   median_dPDUI=round(float(sub.median()), 4)))
pd.DataFrame(strat_rows).to_csv(os.path.join(OUT, "p2_full_qc_stratified.tsv"), sep="\t", index=False)
print("\n[d1 strata]")
print(pd.DataFrame(strat_rows)[lambda d: d.timepoint == "day1"].to_string(index=False), flush=True)
print("[DONE] layer1 QC", flush=True)
