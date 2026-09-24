#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""摊牌脚本：junction F 到底是什么、在哪"""
import sys

sys.path.insert(0, "/mnt/d/stroke_apa_reaudit/scripts")
from g1_3_seq_utils import fetch_pyfaidx, fetch_samtools

GENOME = "/mnt/d/stroke_apa_reaudit/reference/chr5.fa"
OLIGO = "GGAGTAACCGCTTCCTAAACCA"

def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]

seg1 = fetch_pyfaidx(GENOME, "chr5", 122453512, 122454304, "-")
seg2 = fetch_pyfaidx(GENOME, "chr5", 122455892, 122456500, "-")
cdna = seg1 + seg2

print("OLIGO 在 cdna 中的位置:", cdna.find(OLIGO), "| rc(OLIGO):", cdna.find(rc(OLIGO)))
print("cdna[6:28]  =", cdna[6:28])
print("cdna[28:50] =", cdna[28:50])
print()
print("终端外显子 DNA：")
gpos = fetch_pyfaidx(GENOME, "chr5", 122454250, 122454320, "+")
print("  genome+ 122454250-122454319:", gpos)
print("  其中 rc(OLIGO) 起点(0-based,相对窗口):", gpos.find(rc(OLIGO)))
print("  负链读向（revcomp 窗口）:", rc(gpos)[:40], "...")
print()
# samtools 独立验证 seg1 前 60
s1 = fetch_samtools(GENOME, "chr5", 122453512, 122454304, "-")
print("samtools seg1[0:60] =", s1[:60])
print("pyfaidx  seg1[0:60] =", seg1[:60])
print("一致:", s1[:60] == seg1[:60])

# 179939 mRNA 记录核对
name = None
buf = []
rec = {}
for line in open("/mnt/d/stroke_apa_reaudit/reference/gencode.vM25.transcripts.fa"):
    if line.startswith(">"):
        if name:
            rec[name] = "".join(buf)
        name = line[1:].split("|")[0]
        buf = []
    else:
        buf.append(line.strip())
rec[name] = "".join(buf)
t = rec["ENSMUST00000179939.7"]
print("\n179939 记录长:", len(t))
print("OLIGO 在 179939 中的位置:", t.find(OLIGO), "| rc:", t.find(rc(OLIGO)))
p = t.find(OLIGO)
if p >= 0:
    print("  上文[4860:4889]:", t[4860:4889])
    print("  命中[4889:4911]:", t[4889:4911])
    print("  下文[4911:4941]:", t[4911:4941])
q = t.find(rc(OLIGO))
if q >= 0:
    print("  rc 命中位置:", q)
