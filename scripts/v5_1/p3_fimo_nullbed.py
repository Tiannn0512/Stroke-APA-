# -*- coding: utf-8 -*-
"""Null BED for lost-segment scan: matched on length distribution + chromosome usage, seed 42."""
import numpy as np, pandas as pd
bed = pd.read_csv("/mnt/d/stroke_apa/results/fimo_seqs/lost_segments.bed", sep="\t", header=None,
                  names=["chr", "start", "end", "name"])
lens = (bed["end"] - bed["start"]).values
chroms = {}
for ln in open("/home/taylor/reference_v2/chrom_sizes.txt"):
    c, s = ln.split()
    chroms[c] = int(s)
major = [c for c in bed["chr"].unique() if c in chroms]
rng = np.random.default_rng(42)
rows = []
for i, L in enumerate(lens):
    c = major[i % len(major)]
    start = int(rng.integers(0, max(chroms[c] - L - 1, 1)))
    rows.append((c, start, start + L, f"null_{i}"))
pd.DataFrame(rows, columns=["chr", "start", "end", "name"]).to_csv(
    "/mnt/d/stroke_apa/results/fimo_seqs/null_segments.bed", sep="\t", header=False, index=False)
print(f"[nullbed] {len(rows)} regions")
