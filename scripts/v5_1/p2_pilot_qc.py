# -*- coding: utf-8 -*-
"""P2-05g pilot PDUI QC — descriptive stats only (no thresholds, no gates).
Reads DaPars2 subset PDUI output, reports: table shape, per-sample validity,
sham replicate correlation, day1-vs-sham deltaPDUI descriptive counts."""
import glob, os
import numpy as np
import pandas as pd

BASE = r"D:\stroke_apa_data\dapars2"
OUT = r"D:\stroke_apa\results"
SAMPLES = ["sham1", "sham2", "day1_rep1", "day1_rep2"]

files = sorted(glob.glob(os.path.join(BASE, "dapars2_subset*PDUI*.tsv")))
print("PDUI files found:", [os.path.basename(f) for f in files])
if not files:
    raise SystemExit("[QC FAIL] no PDUI output found in " + BASE)

dfs = []
for f in files:
    df = pd.read_csv(f, sep="\t")
    df["__src"] = os.path.basename(f)
    dfs.append(df)
tab = pd.concat(dfs, ignore_index=True)
idcol = tab.columns[0]
tab = tab.drop_duplicates(subset=[idcol])
print(f"total event rows={len(tab)} (from {len(files)} files), id col = '{idcol}'")

def col_for(s):
    exact = [c for c in tab.columns if str(c) == s]
    if exact:
        return exact[0]
    base = [c for c in tab.columns if s in os.path.basename(str(c))]
    if base:
        return base[0]
    return None

cols = {s: col_for(s) for s in SAMPLES}
missing = [s for s, c in cols.items() if c is None]
if missing:
    raise SystemExit(f"[QC FAIL] sample columns not found: {missing}; available={list(tab.columns)[:10]}")

sub = tab[[cols[s] for s in SAMPLES]].copy()
sub.columns = SAMPLES
sub = sub.apply(pd.to_numeric, errors="coerce")

qc = [("events_total", len(tab))] + [(f"pdui_valid_{s}", int(sub[s].notna().sum())) for s in SAMPLES]
r_sham = sub["sham1"].corr(sub["sham2"])
qc.append(("sham1_vs_sham2_pearson_r", round(float(r_sham), 4)))
sham_mean = sub[["sham1", "sham2"]].mean(axis=1)
for rep in ("day1_rep1", "day1_rep2"):
    d = sub[rep] - sham_mean
    d = d.dropna()
    qc.append((f"{rep}_vs_sham_n", int(len(d))))
    for th in (0.1, 0.2):
        qc.append((f"{rep}_deltaPDUI<-{th}", int((d < -th).sum())))
        qc.append((f"{rep}_deltaPDUI>+{th}", int((d > th).sum())))
    qc.append((f"{rep}_median_deltaPDUI", round(float(d.median()), 4)))

qcdf = pd.DataFrame(qc, columns=["metric", "value"])
qcdf.to_csv(os.path.join(OUT, "p2_pilot_pdui_qc.tsv"), sep="\t", index=False)
print(qcdf.to_string(index=False))
print("[QC DONE]")
