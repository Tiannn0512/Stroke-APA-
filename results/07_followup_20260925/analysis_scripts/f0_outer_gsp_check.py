#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""核对两个 outer GSP 候选序列在 chr5 上的真实位置（pyfaidx, BED 0-based half-open）。"""
from pyfaidx import Fasta

GENOME = "/mnt/d/stroke_apa_reaudit/reference/chr5.fa"
fa = Fasta(GENOME)
ch = fa["chr5"]

def rc(s):
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))

def find_all(query, label):
    q = query.upper()
    hits = []
    start = 0
    while True:
        i = str(ch[start:]).find(q)  # chunked below; chr5 154Mb ok in slices
        break
    return hits

# 直接在目标窗口附近找（±500bp），再全染色体内精确找一次
candidates = {
    "table_seq_GACAGGAA": "GACAGGAATTAAAGCCTCAGGAAACC",
    "doc_seq_TCAAGAAG": "TCAAGAAGTGAAGTGACATGGACAAG",
}
lo, hi = 122459602 - 500, 122459626 + 500
win = str(ch[lo:hi]).upper()
for name, q in candidates.items():
    print(f"== {name} {q}")
    # 窗口内
    idx = win.find(q)
    while idx != -1:
        print(f"  near-window + hit at chr5:{lo+idx}-{lo+idx+len(q)} (0-based)")
        idx = win.find(q, idx + 1)
    idx = win.find(rc(q))
    while idx != -1:
        print(f"  near-window - hit at chr5:{lo+idx}-{lo+idx+len(q)} (0-based)")
        idx = win.find(rc(q), idx + 1)
    # 全染色体（+ 与 -）
    seq = str(ch[:]).upper()
    start = 0
    n_plus = n_minus = 0
    while True:
        i = seq.find(q, start)
        if i == -1:
            break
        n_plus += 1
        if n_plus <= 5:
            print(f"  chr5 whole + exact hit {i}-{i+len(q)}")
        start = i + 1
    start = 0
    qr = rc(q)
    while True:
        i = seq.find(qr, start)
        if i == -1:
            break
        n_minus += 1
        if n_minus <= 5:
            print(f"  chr5 whole - exact hit {i}-{i+len(q)}")
        start = i + 1
    print(f"  totals on chr5: + {n_plus}, - {n_minus}")
