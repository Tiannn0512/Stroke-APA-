# -*- coding: utf-8 -*-
"""P2-06 full merge + depth re-verification + descriptive QC (runs in WSL).
1. Depth check: every chr run.log's printed depth array must equal the canonical
   12-value table (lesson 2026-09-18: numbers must trace to logs).
2. Merge: glob dapars2_full/chr*_chr*/*_result_temp.chr*.txt -> single matrix on /mnt/d.
3. Descriptive QC (same metrics as pilot): valid PDUI per sample, sham r,
   per-timepoint shortening/lengthening at |dPDUI|>0.1 (SAP layer-1 descriptive input).
Run: /home/taylor/miniconda3/envs/dapars2/bin/python p2_full_merge_qc.py
"""
import glob, os, re, sys
import pandas as pd
import numpy as np

BASE = "/home/taylor/stroke_APA_work"
OUTDIR = "/mnt/d/stroke_apa/results"
DEPTH = {"sham1": 49652326, "sham2": 52273421, "day1_rep1": 39837975, "day1_rep2": 39038280,
         "day3_rep1": 38560545, "day3_rep2": 40906253, "day7_rep1": 39861605, "day7_rep2": 44513076,
         "day21_rep1": 38391433, "day21_rep2": 41344758, "day60_rep1": 40347218, "day60_rep2": 51568111}

# ---- 1. depth verification from every run.log ----
logs = sorted(glob.glob(f"{BASE}/dapars2_full/chr*/run.log"))
# canonical depth order = the depth FILE's row order (DaPars2 prints in file order)
exp_rows = []
with open(f"{BASE}/sequencing_depth_dapars2.tsv") as f:
    for ln in f:
        if ln.strip():
            name, val = ln.rstrip("\n").split("\t")
            exp_rows.append((name, int(val)))
            assert DEPTH[name] == int(val), f"depth file vs table mismatch: {name}"
exp = [v for _, v in exp_rows]
n_ok = n_bad = 0
for lg in logs:
    chrn = os.path.basename(os.path.dirname(lg))
    with open(lg) as f:
        txt = f.read()
    m = re.search(r"\[(\d+(?:\s+\d+)*)\s*\]", txt)
    if not m:
        print(f"[DEPTH] {chrn}: NO ARRAY FOUND")
        n_bad += 1
        continue
    got = [int(x) for x in m.group(0).strip("[]").split()]
    if got == exp:
        n_ok += 1
    else:
        print(f"[DEPTH] {chrn}: MISMATCH got={got}")
        n_bad += 1
print(f"[DEPTH] verified_ok={n_ok} mismatch={n_bad}", flush=True)
if n_bad:
    print("[DEPTH] FATAL: depth mismatch — investigate before trusting merge", flush=True)
    sys.exit(2)

# ---- 2. merge ----
files = sorted(glob.glob(f"{BASE}/dapars2_full/chr*_chr*/*_result_temp.chr*.txt"))
real = [f for f in files if not f.split("_result_temp")[0].endswith("chr22")]
done_logs = sum(1 for lg in logs if "exit=0" in open(lg).read())
print(f"[MERGE] {len(files)} result files ({len(real)} real chr); logs exit=0: {done_logs}/20", flush=True)
if len(real) < 20 or done_logs < 20:
    print(f"[MERGE] SMOKE MODE: run incomplete — no files written", flush=True)
    sys.exit(0)
frames = []
for fp in files:
    df = pd.read_csv(fp, sep="\t")
    if df.empty:
        print(f"  [skip empty] {os.path.basename(fp)}")
        continue
    frames.append(df)
mat = pd.concat(frames, ignore_index=True)
pdui_cols = [c for c in mat.columns if c.endswith("_PDUI")]
samples = [c[:-5].split("/")[-1] for c in pdui_cols]
mat = mat.rename(columns=dict(zip(pdui_cols, samples)))
mat["symbol"] = [str(g).split("|")[2] if "|" in str(g) else str(g) for g in mat["Gene"]]
out_fp = os.path.join(OUTDIR, "p2_full_pdui_matrix.tsv")
mat.to_csv(out_fp, sep="\t", index=False)
print(f"[MERGE] events={len(mat)} cols={mat.shape[1]} -> {out_fp}", flush=True)

# ---- 3. descriptive QC ----
rows = []
valid = {}
for s in samples:
    v = pd.to_numeric(mat[s], errors="coerce")
    valid[s] = v
    rows.append(dict(metric=f"pdui_valid_{s}", value=int(v.notna().sum())))
sham = pd.concat([valid["sham1"], valid["sham2"]], axis=1).dropna()
r = sham["sham1"].corr(sham["sham2"])
rows.append(dict(metric="sham1_sham2_pearson_r", value=round(float(r), 4)))
sham_mean = (valid["sham1"] + valid["sham2"]) / 2
for tp in ["day1", "day3", "day7", "day21", "day60"]:
    cols = [s for s in samples if s.startswith(tp)]
    d = pd.concat([valid[c] for c in cols], axis=1).dropna()
    dmean = d.mean(axis=1)
    common = dmean.dropna().index.intersection(sham_mean.dropna().index)
    dd = (dmean - sham_mean).loc[common]
    short = int((dd < -0.1).sum()); long_ = int((dd > 0.1).sum())
    rows.append(dict(metric=f"{tp}_n_vs_sham", value=len(common)))
    rows.append(dict(metric=f"{tp}_shortening_lt-0.1", value=short))
    rows.append(dict(metric=f"{tp}_lengthening_gt+0.1", value=long_))
    rows.append(dict(metric=f"{tp}_ratio", value=round(short / max(long_, 1), 2)))
    rows.append(dict(metric=f"{tp}_median_dPDUI", value=round(float(dd.median()), 4)))
qc = pd.DataFrame(rows)
qc.to_csv(os.path.join(OUTDIR, "p2_full_merge_qc.tsv"), sep="\t", index=False)
print(qc.to_string(index=False), flush=True)
print("[DONE] full merge + QC", flush=True)
