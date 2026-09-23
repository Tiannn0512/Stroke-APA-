#!/usr/bin/env python3
"""Estimate chr6 grid-search cost vs chr2 (cost ~ sum of span^2 x 12 samples)."""
import sys
n = s = s2 = 0
max_span = 0
max_name = ""
spans6 = []
for ln in open("/home/taylor/reference_v2/gencode_M25_3UTR_slim.bed"):
    f = ln.rstrip("\n").split("\t")
    if f[0] != "chr6":
        continue
    span = int(f[2]) - int(f[1])
    n += 1; s += span; s2 += span * span
    spans6.append(span)
    if span > max_span:
        max_span, max_name = span, f[3]
print(f"chr6: n={n} sum_span={s} sum_span2={s2:.3g} max_span={max_span} {max_name}")
n2 = s2_ = 0
for ln in open("/home/taylor/reference_v2/gencode_M25_3UTR_slim.bed"):
    f = ln.rstrip("\n").split("\t")
    if f[0] != "chr2":
        continue
    span = int(f[2]) - int(f[1])
    n2 += 1; s2_ += span * span
print(f"chr2: n={n2} sum_span2={s2_:.3g}")
print(f"cost ratio chr6/chr2 (span^2): {s2/s2_:.2f}")
spans6.sort(reverse=True)
print("chr6 top-5 spans:", spans6[:5])
print("chr6 spans>50kb:", sum(1 for x in spans6 if x > 50000),
      " >100kb:", sum(1 for x in spans6 if x > 100000))
