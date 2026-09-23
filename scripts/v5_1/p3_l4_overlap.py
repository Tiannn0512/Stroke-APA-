#!/usr/bin/env python3
"""Merge 870M-row per-base element records (streaming) is done in bash; this script
does the lost-segment overlap via per-chromosome sorted numpy sweep (memory-safe)."""
import numpy as np, pandas as pd, bisect, sys

MERGED = "/home/taylor/reference_v2/phastcons60way/elements_merged.bed"
LOST = "/mnt/d/stroke_apa/results/fimo_seqs/lost_segments.bed"
OUT = "/home/taylor/stroke_APA_work/fimo_seq/lost_segments_conservation.tsv"

# load merged elements per chromosome
by_chr = {}
for ln in open(MERGED):
    c, s, e = ln.split()
    by_chr.setdefault(c, ([], []))
    by_chr[c][0].append(int(s))
    by_chr[c][1].append(int(e))
for c in by_chr:
    s, e = by_chr[c]
    o = np.argsort(np.array(s))
    by_chr[c] = (np.array(s)[o], np.array(e)[o])
print(f"chroms: {len(by_chr)}; total merged elements: {sum(len(v[0]) for v in by_chr.values())}", flush=True)

rows = []
for ln in open(LOST):
    c, s, e, name = ln.rstrip("\n").split("\t")
    s, e = int(s), int(e)
    tot = 0
    if c in by_chr:
        S, E = by_chr[c]
        i = bisect.bisect_right(S, s) - 1
        j = i
        # walk overlapping elements
        if j < 0:
            j = 0
        while j < len(S) and S[j] < e:
            ov = min(e, E[j]) - max(s, S[j])
            if ov > 0:
                tot += ov
            j += 1
    rows.append((name, "CONSERVED" if tot >= 10 else "not", tot))
with open(OUT, "w") as f:
    for name, call, tot in rows:
        f.write(f"{name}\t{call}\t{tot}\n")
n_cons = sum(1 for _, c, _ in rows if c == "CONSERVED")
print(f"segments scored: {len(rows)}; CONSERVED(>=10bp): {n_cons}", flush=True)
print("[DONE] L4 overlap sweep", flush=True)
