# -*- coding: utf-8 -*-
"""Null BED v2: random fragments cut from NON-significant events' distal 3'UTR segments
(same UTR universe), length-matched to candidate lost segments. GC checked at the end."""
import numpy as np, pandas as pd

OUT = "/mnt/d/stroke_apa/results"
mat = pd.read_csv(f"{OUT}/p2_full_pdui_matrix.tsv", sep="\t")
mat1 = mat.copy(); mat1.index = mat1.index + 1
ev = pd.read_csv(f"{OUT}/p2_sap_layer2_events.tsv", sep="\t")
ev["Gene"] = ev["event"].astype(int).map(mat1["Gene"])
ev["symbol"] = ev["Gene"].astype(str).str.split("|").str[0]
ev["loci"] = ev["event"].astype(int).map(mat1["Loci"])
ev["prox"] = ev["event"].astype(int).map(mat1["Predicted_Proximal_APA"]).astype(float)
ev["chr"], ev["utr_end"] = ev["loci"].str.extract(r"(chr[^:]+):(\d+)-(\d+)")[0].values, \
    ev["loci"].str.extract(r":(\d+)-(\d+)")[1].astype(int).values
ev = ev.dropna(subset=["prox", "chr"])
ev["prox"] = np.floor(ev["prox"]).astype(int)
sig_ids = set(ev[ev["sig_joint"] == True]["event"].astype(str))
short_ids = set(ev[(ev["sig_joint"] == True) & (ev["dPDUI"] < 0)]["event"].astype(str))

# non-significant events with usable distal segments
bg = ev[~ev["event"].astype(str).isin(sig_ids)].copy()
bg = bg[bg["utr_end"] > bg["prox"] + 500]  # need >=500bp to cut from
print(f"null donor pool: {len(bg)} non-sig distal segments")

# candidate lost segment lengths
cand = pd.read_csv(f"{OUT}/fimo_seqs/lost_segments.bed", sep="\t", header=None,
                   names=["chr", "start", "end", "name"])
cand["len"] = cand["end"] - cand["start"]
# cap length at 1000 (very long lost segments -> cut null at same length from donor)
rng = np.random.default_rng(42)
donors = bg.sample(frac=1.0, random_state=42).reset_index(drop=True)
rows = []
di = 0
for i, L in enumerate(cand["len"].values):
    L = min(int(L), 2000)  # cap: donor segments must spare 500
    placed = False
    for _ in range(50):
        if di >= len(donors):
            donors = donors.sample(frac=1.0, random_state=i).reset_index(drop=True); di = 0
        d = donors.iloc[di]; di += 1
        span = d["utr_end"] - d["prox"]
        if span < L + 100:
            continue
        off = int(rng.integers(d["prox"], d["utr_end"] - L + 1))
        rows.append((d["chr"], off, off + L, f"unull_{i}"))
        placed = True
        break
    if not placed:
        continue
n = pd.DataFrame(rows, columns=["chr", "start", "end", "name"])
n.to_csv(f"{OUT}/fimo_seqs/null_utr.bed", sep="\t", header=False, index=False)
print(f"null v2: {len(n)} UTR-derived regions, total {(n['end']-n['start']).sum()/1e6:.2f} Mb "
      f"(candidates: {len(cand)} regions, {(cand['len']).sum()/1e6:.2f} Mb)")
