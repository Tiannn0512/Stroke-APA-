#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""一次性对证：mRNA 记录段 vs 基因组正链 vs 引物寡核苷酸"""
import subprocess
import sys

sys.path.insert(0, "/mnt/d/stroke_apa_reaudit/scripts")
from g1_3_seq_utils import fetch_pyfaidx

GENOME = "/mnt/d/stroke_apa_reaudit/reference/chr5.fa"
FA_T = "/mnt/d/stroke_apa_reaudit/reference/gencode.vM25.transcripts.fa"

def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]

# 载入 179939 记录
name, buf, rec = None, [], {}
for line in open(FA_T):
    if line.startswith(">"):
        if name:
            rec[name] = "".join(buf)
        name = line[1:].split("|")[0]
        buf = []
    else:
        buf.append(line.strip())
rec[name] = "".join(buf)
m = rec["ENSMUST00000179939.7"]

print("=== junction R 引物（ATCCTACTATGCGGCAGAACAG）===")
print("mRNA 记录[4545:4567]      =", m[4545:4567])
gplus = fetch_pyfaidx(GENOME, "chr5", 122456206, 122456229, "+")
print("genome+ [122456206,122456229) =", gplus)
print("rc(genome+)                   =", rc(gplus))
print("→ R == mRNA 段?  ", m[4545:4567] == "ATCCTACTATGCGGCAGAACAG")
print("→ R == rc(genome+)?", rc(gplus) == "ATCCTACTATGCGGCAGAACAG")

print()
print("=== junction F 引物（GGAGTAACCGCTTCCTAAACCA）===")
print("mRNA 记录[4889:4911]      =", m[4889:4911])
g2 = fetch_pyfaidx(GENOME, "chr5", 122454275, 122454297, "+")
print("genome+ [122454275,122454297) =", g2)
print("rc(genome+)                   =", rc(g2))

print()
print("=== distal R 引物（TTGTAAGTGGCCAGATTGCTCT）===")
print("mRNA 记录[5270:5292]      =", m[5270:5292])
g3 = fetch_pyfaidx(GENOME, "chr5", 122453894, 122453916, "+")
print("genome+ [122453894,122453916) =", g3)
print("rc(genome+)                   =", rc(g3))

print()
print("=== common F 引物（TGAAGTGACATGGACAAGCAGA）===")
print("mRNA 记录[3295:3317] (=rc(F) 应在此) =", m[3295:3317])
print("rc(该段)                              =", rc(m[3295:3317]))
g4 = fetch_pyfaidx(GENOME, "chr5", 122457542, 122457564, "+")
print("genome+ [122457542,122457564) =", g4)
