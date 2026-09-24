#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""G1 复核 · 工作 1 附加：Atp2a2 三条转录本的 3' 端结构 + 窗口内全部基因。
（GTF 1-based closed；输出统一转 BED 0-based half-open）
"""
import gzip
from collections import defaultdict

GTF = "/mnt/d/stroke_apa_reaudit/reference/gencode.vM25.annotation.gtf.gz"
TARGETS = ["ENSMUST00000177974.7", "ENSMUST00000179939.7", "ENSMUST00000031423.9"]

exons = defaultdict(list)
utrs = defaultdict(list)
strand = {}
genes = {}

with gzip.open(GTF, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if p[0] != "chr5":
            continue
        s0, e0 = int(p[3]) - 1, int(p[4])   # BED 0-based
        feat = p[2]
        attr = p[8]
        if feat == "gene":
            if e0 >= 122452512 and s0 <= 122458153:
                gname = attr.split('gene_name "')[1].split('"')[0]
                genes[(s0, e0, gname)] = attr.split('gene_id "')[1].split('"')[0]
            continue
        if 'gene_name "Atp2a2"' not in attr:
            continue
        tid = attr.split('transcript_id "')[1].split('"')[0]
        if tid not in TARGETS:
            continue
        strand[tid] = p[6]
        if feat == "exon":
            exons[tid].append((s0, e0))
        elif feat == "UTR":
            utrs[tid].append((s0, e0))

print("===== 窗口 chr5:122452512-122458153 内的基因 =====")
for (s, e, n) in sorted(genes):
    print(f"  chr5:{s}-{e}  {n}  ({genes[(s,e,n)]})")

print("\n===== 三条转录本 3' 端结构（负链：坐标越小越靠 3' 端） =====")
for t in TARGETS:
    ex = sorted(exons[t])
    st = strand[t]
    term = ex[0]  # 负链终端外显子 = 坐标最小的外显子
    print(f"\n{t}  strand={st}  外显子数={len(ex)}")
    print(f"  终端外显子(3'末端): chr5:{term[0]}-{term[1]}  长 {term[1]-term[0]}bp")
    print(f"  其余低坐标侧外显子: {ex[1:4]}")
    print(f"  终端内含子(终端外显子-次末外显子): chr5:{term[1]}-{ex[1][0] if len(ex)>1 else '?'}")
    u = sorted(utrs[t])
    print(f"  UTR 块: {u if u else '无'}")
    # 终端外显子内的 CDS/UTR 划分
    if u:
        print(f"  3'UTR 范围(该转录本): chr5:{min(x[0] for x in u)}-{max(x[1] for x in u)}")

print("\n===== 事件四区段与各转录本终端外显子/UTR 的关系 =====")
segments = {
    "distal 全段":        (122453512, 122456498),
    "重叠段(177974 UTR)": (122453512, 122454293),
    "独有段":             (122454293, 122456498),
    "shared":             (122456498, 122457153),
}
for name, (s, e) in segments.items():
    print(f"  {name}: chr5:{s}-{e}")
    for t in TARGETS:
        term = sorted(exons[t])[0]
        u = sorted(utrs[t])
        us = (min(x[0] for x in u), max(x[1] for x in u)) if u else None
        in_term = "∈终端外显子" if term[0] <= s and e <= term[1] else "部分/不在"
        in_utr = ""
        if us and s >= us[0] and e <= us[1]:
            in_utr = "；∈该转录本3'UTR"
        elif us:
            ov = max(0, min(e, us[1]) - max(s, us[0]))
            in_utr = f"；与该转录本3'UTR重叠 {ov}bp"
        print(f"    {t}: {in_term}{in_utr}")
