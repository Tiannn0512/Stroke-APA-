#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Gm30970 与 L/S/M mRNA 的同源区精确作图（滑窗最小错配）。"""
BASE = "/mnt/d/stroke_apa_reaudit"
SEQS, cur = {}, None
with open(f"{BASE}/reference/gencode.vM25.transcripts.fa") as fh:
    for line in fh:
        if line.startswith(">"):
            cur = line[1:].split("|")[0].strip(); SEQS[cur] = []
        else:
            SEQS[cur].append(line.strip())
SEQS = {k: "".join(v).upper() for k, v in SEQS.items()}
L_, S_, M_, GMID = ("ENSMUST00000179939.7", "ENSMUST00000177974.7",
                    "ENSMUST00000031423.9", "ENSMUST00000227710.1")
GM = SEQS[GMID]
print(f"Gm30970 len={len(GM)}")

def rc(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))

for name, tid in (("L", L_), ("S", S_), ("M", M_)):
    s = SEQS[tid]
    # 每个 50nt GM 窗口 在 mRNA 上的最佳命中（步长 10）
    print(f"== {name} (cdna {len(s)})")
    spans = []
    for qstart in range(0, len(GM) - 49, 10):
        q = GM[qstart:qstart + 50]
        best = (99, -1)
        for i in range(len(s) - 49):
            mm = sum(1 for a, b in zip(q, s[i:i+50]) if a != b)
            if mm < best[0]:
                best = (mm, i)
        spans.append((qstart, best[0], best[1]))
    for qstart, mm, pos in spans:
        bar = "#" * max(0, 50 - mm * 5)
        print(f"  GM{qstart:4d}-{qstart+50:4d} -> {name}mRNA {pos+1:5d} mm={mm:2d} {bar}")
