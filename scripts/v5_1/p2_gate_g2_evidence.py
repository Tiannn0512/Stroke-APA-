# -*- coding: utf-8 -*-
"""Gate G2 evidence extras: LOSO global stability, d3 anomalous-replicate exclusion check,
per-contrast direction splits + union, top-event table for spot checks."""
import numpy as np, pandas as pd

OUT = "/mnt/d/stroke_apa/results"
mat = pd.read_csv(f"{OUT}/p2_full_pdui_matrix.tsv", sep="\t")
SAMPLES = ["sham1","sham2","day1_rep1","day1_rep2","day3_rep1","day3_rep2",
           "day7_rep1","day7_rep2","day21_rep1","day21_rep2","day60_rep1","day60_rep2"]
P = mat[SAMPLES].apply(lambda c: pd.to_numeric(c, errors="coerce"))
sham_mean = P[["sham1","sham2"]].mean(axis=1)
lines = []

def global_ratios(tag, dd, note=""):
    short = int((dd < -0.1).sum()); long_ = int((dd > 0.1).sum())
    lines.append(f"{tag}\tshort={short}\tlengthen={long_}\tratio={short/max(long_,1):.2f}\tmedian={dd.median():.4f}{note}")

# d1 short/long per replicate (LOSO-style sanity: both reps independently vs each single sham? keep simple: per-rep vs sham-mean)
lines.append("=== per-replicate d1 vs sham-mean (replicate-level consistency) ===")
for tp in ["day1","day7","day21","day60"]:
    for r in ["rep1","rep2"]:
        dd = (P[f"{tp}_{r}"] - sham_mean).dropna()
        global_ratios(f"{tp}_{r}", dd)

lines.append("=== LOSO: drop one sham, recompute d1/d7 global ratios ===")
for drop in ["sham1","sham2"]:
    sm = P[[s for s in ["sham1","sham2"] if s != drop]].mean(axis=1)
    for tp in ["day1","day7"]:
        dd = (P[[f"{tp}_rep1",f"{tp}_rep2"]].mean(axis=1) - sm).dropna()
        global_ratios(f"{tp}_vs_{drop}-excluded-sham", dd)

lines.append("=== LOSO: drop one timepoint replicate (d1) ===")
for keep in ["day1_rep1","day1_rep2"]:
    dd = (P[keep] - sham_mean).dropna()
    global_ratios(f"d1_{keep}-only", dd)

lines.append("=== d3 anomaly: exclude day3_rep1 (insert=150, dup=0.215) ===")
dd = (P["day3_rep2"] - sham_mean).dropna()
global_ratios("d3_rep2-only", dd)
dd = (P["day3_rep1"] - sham_mean).dropna()
global_ratios("d3_rep1-only", dd)
# drop the anomalous sample from sham side? no - anomaly is on day3 side.

ev = pd.read_csv(f"{OUT}/p2_sap_layer2_events.tsv", sep="\t")
om = pd.read_csv(f"{OUT}/p2_sap_layer2_omnibus.tsv", sep="\t")
sig = ev[ev["sig_joint"] == True]
lines.append("=== per-contrast significant direction split (sig_joint) ===")
for tp, g in sig.groupby("contrast"):
    s = int((g["dPDUI"] < 0).sum()); l = int((g["dPDUI"] > 0).sum())
    lines.append(f"{tp}\tsig={len(g)}\tshortening={s}\tlengthening={l}")
union = set(sig["event"])
lines.append(f"union_sig_events\t{len(union)}")
# cross-time concordance among union
wide = sig.pivot_table(index="event", columns="contrast", values="dPDUI", aggfunc="first")
concord2 = (wide.apply(lambda r: sum(1 for v in r if pd.notna(v)) >= 2 and
                       (r.dropna() < 0).all() or (r.dropna() > 0).all(), axis=1)).sum()
lines.append(f"union_events_significant_at_>=2_timepoints\t{int(concord2)}")

top = sig.reindex(sig["padj_joint"].abs().sort_values().index).head(40)
cols = ["event","contrast","logitFC","dPDUI","d_rep1","d_rep2","padj_joint","padj_within"]
top[cols].to_csv(f"{OUT}/p2_gate_g2_top_events.tsv", sep="\t", index=False)
lines.append(f"top40_written\t{len(top)}")
lines.append(f"omnibus_padj<0.1\t{int((om['padj_omnibus'] < 0.1).sum())}/{len(om)}")

with open(f"{OUT}/p2_gate_g2_evidence_extras.tsv", "w") as f:
    f.write("\n".join(lines) + "\n")
print("\n".join(lines))
print("[DONE] G2 evidence extras")
