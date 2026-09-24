#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""定位 Gm30970 假基因记录覆盖 Atp2a2 mRNA 的区间，为设计选址提供依据"""
import sys

sys.path.insert(0, "/mnt/d/stroke_apa_reaudit/scripts")
from g1_3_specificity import load_transcripts, rc

tx = load_transcripts()
pg = tx["ENSMUST00000227710.1"]
m = tx["ENSMUST00000179939.7"]
print(f"Gm30970 记录长: {len(pg)}, 179939 长: {len(m)}")

# 正反两方向找最长精确匹配
def longest_match(a, b):
    best = (0, 0, 0, 0)
    for i in range(0, len(a), 10):
        for L in range(len(b), 20, -20):
            if i + L > len(a):
                continue
            p = b.find(a[i:i + L])
            if p != -1:
                if L > best[0]:
                    best = (L, i, p, "+")
                break
    return best

L1, ai1, pi1, o1 = longest_match(pg, m)
print(f"pg→m 正向最长精确匹配: {L1}bp (mRNA 位置 {pi1}-{pi1+L1})")
L2, ai2, pi2, o2 = longest_match(rc(pg), m)
print(f"rc(pg)→m 最长精确匹配: {L2}bp (mRNA 位置 {pi2}-{pi2+L2})")

# 滑窗粗定位：用 pg 的每 30mer 在 m 上的命中分布
hits = []
for i in range(0, len(pg) - 30 + 1, 30):
    p = m.find(pg[i:i+30])
    pr = m.find(rc(pg[i:i+30]))
    hits.append((i, p if p != -1 else None, pr if pr != -1 else None))
for h in hits:
    print("  pg[%4d:%4d] → mRNA fwd=%s rc=%s" % h)
