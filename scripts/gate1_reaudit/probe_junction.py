#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""实证 junction F 引物的序列来源与拷贝数"""
import subprocess
import sys

sys.path.insert(0, "/mnt/d/stroke_apa_reaudit/scripts")
from g1_3_seq_utils import fetch_pyfaidx, fetch_samtools

OLIGO = "GGAGTAACCGCTTCCTAAACCA"
GENOME = "/mnt/d/stroke_apa_reaudit/reference/chr5.fa"

def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]

# 1) 设计映射声称：cdna[28:50] 应为基因组 122454275-122454297（BED, 负链）的 revcomp
seg1 = fetch_pyfaidx(GENOME, "chr5", 122453512, 122454304, "-")
print("seg1 len:", len(seg1))
print("cdna[28:50] =", seg1[28:50])
print("预期 oligo  =", OLIGO, "→", "一致" if seg1[28:50] == OLIGO else "不一致!!")

# 2) 基因组 + 链上该 22-mer 的全部拷贝
seq = "".join(open(GENOME).read().split("\n")[1:]).upper()
hits = []
i = seq.find(OLIGO)
while i != -1:
    hits.append(i)
    i = seq.find(OLIGO, i + 1)
hits_rc = []
j = seq.find(rc(OLIGO))
while j != -1:
    hits_rc.append(j)
    j = seq.find(rc(OLIGO), j + 1)
print("oligo 在 chr5 + 链命中:", [h + 1 for h in hits], "(1-based)")
print("rc(oligo) 在 chr5 + 链命中:", [h + 1 for h in hits_rc], "(1-based)")
print("=> 目标位点(负链模板来源) rc 命中应含 122454276:",
      any(122454270 <= h + 1 <= 122454300 for h in hits_rc))
print("=> penultimate 外显子(122455892-122457153) 内的正向拷贝:",
      [h + 1 for h in hits if 122455892 <= h + 1 <= 122457153])

# 3) ENSMUST00000227710.1 是什么基因
import gzip
with gzip.open("/mnt/d/stroke_apa_reaudit/reference/gencode.vM25.annotation.gtf.gz", "rt") as f:
    for line in f:
        if "ENSMUST00000227710" in line and "\ttranscript\t" in line:
            print("\n227710 transcript 行:", line.strip()[:220])
            break
