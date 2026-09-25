#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""调试 J_S 跨接头序列 vs S 真实转录本。"""
from pyfaidx import Fasta

BASE = "/mnt/d/stroke_apa_reaudit"
GENOME = Fasta(f"{BASE}/reference/chr5.fa")["chr5"]

SEQS, cur = {}, None
with open(f"{BASE}/reference/gencode.vM25.transcripts.fa") as fh:
    for line in fh:
        if line.startswith(">"):
            cur = line[1:].split("|")[0].strip(); SEQS[cur] = []
        else:
            SEQS[cur].append(line.strip())
SEQS = {k: "".join(v).upper() for k, v in SEQS.items()}
S_ = "ENSMUST00000177974.7"
s = SEQS[S_]

def rc(x):
    return x[::-1].translate(str.maketrans("ACGT", "TGCA"))

left = rc(str(GENOME[122457290:122457303]))
right = rc(str(GENOME[122454291:122454304]))
j = left + right
print("left :", left)
print("right:", right)
print("J_S  :", j)
i = s.find(left)
print("left in S at (0-based):", i)
if i >= 0:
    print("S context:", s[i-10:i+40])
    print("j  vs  S :", j[:13], "|", s[i:i+26])
# 在 S 中找 right
k = s.find(right)
print("right in S at:", k)
# 打印 S 接头区域: 应该在 penult(122457302-457426 124bp)->terminal
# S: penult mRNA 起点 = 找 exon 122457302-457426 在 S 的 cDNA 起点
# 用 genome 搜索：取 penult 前 20 碱基(mRNA sense)
pen20 = rc(str(GENOME[122457406:122457426]))
print("pen20 in S at:", s.find(pen20), "(penult exon 起)")
p = s.find(pen20)
if p >= 0:
    print("S [p-5:p+140]:", s[p-5:p+140])
