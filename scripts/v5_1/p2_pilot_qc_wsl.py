"""WSL version of pilot PDUI QC (reads /mnt/d paths, runs in WSL qapa python).
Descriptive stats only: PDUI counts, sham replicate correlation, day1 deltaPDUI overview."""
import glob, os
import numpy as np
import pandas as pd

BASE = "/mnt/d/stroke_apa_data"
OUT = "/mnt/d/stroke_apa/results"
SAMPLES = ["sham1", "sham2", "day1_rep1", "day1_rep2"]

files = sorted(glob.glob(os.path.join(BASE, "dapars2_single_chr*", "dapars2_slim4_result_temp.chr*.txt")))
print("PDUI files found:", len(files))
if not files:
    raise SystemExit("[QC FAIL] no PDUI output found under " + BASE)

dfs = []
for f in files:
    try:
        df = pd.read_csv(f, sep="\t")
        if len(df):
            dfs.append(df)
    except Exception as e:
        print("[WARN] skip", f, e)
if not dfs:
    raise SystemExit("[QC FAIL] all PDUI files empty/unreadable")
tab = pd.concat(dfs, ignore_index=True)
idcol = tab.columns[0]
tab = tab.drop_duplicates(subset=[idcol])
print(f"total event rows={len(tab)} (from {len(files)} chunk files), id col='{idcol}'")

def col_for(s):
    exact = [c for c in tab.columns if str(c) == s]
    if exact:
        return exact[0]
    cand = [c for c in tab.columns if s in os.path.basename(str(c))]
    return cand[0] if cand else None

cols = {s: col_for(s) for s in SAMPLES}
missing = [s for s, c in cols.items() if c is None]
if missing:
    raise SystemExit(f"[QC FAIL] sample columns not found: {missing}; have={list(tab.columns)[:8]}")

sub = tab[[cols[s] for s in SAMPLES]].copy()
sub.columns = SAMPLES
sub = sub.apply(pd.to_numeric, errors="coerce")

qc = [("events_total", len(tab))] + [(f"pdui_valid_{s}", int(sub[s].notna().sum())) for s in SAMPLES]
r_sham = sub["sham1"].corr(sub["sham2"])
qc.append(("sham1_vs_sham2_pearson_r", round(float(r_sham), 4)))
sham_mean = sub[["sham1", "sham2"]].mean(axis=1)
for rep in ("day1_rep1", "day1_rep2"):
    d = (sub[rep] - sham_mean).dropna()
    qc.append((f"{rep}_vs_sham_n", int(len(d))))
    for th in (0.1, 0.2):
        qc.append((f"{rep}_deltaPDUI<-{th}", int((d < -th).sum())))
        qc.append((f"{rep}_deltaPDUI>+{th}", int((d > th).sum())))
    qc.append((f"{rep}_median_deltaPDUI", round(float(d.median()), 4)))

qcdf = pd.DataFrame(qc, columns=["metric", "value"])
qcdf.to_csv(os.path.join(OUT, "p2_pilot_pdui_qc.tsv"), sep="\t", index=False)
print(qcdf.to_string(index=False))
print("[QC DONE]")
