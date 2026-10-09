#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""G1 复核 · 工作 1 自动测试：
A. 合成 CIGAR（含 N/S/H/D）断言剪接数与 soft-clipped 碱基数正确（H 不计）；
B. 接头坐标与 UTR 包含逻辑；
C. 真实 BAM 回归：sham1 审查窗内含 N 读段数 vs samtools 独立口径（-F 0x904, awk $6~/N/）。
全部通过退出码 0。
"""
import subprocess
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from g1_1_read_metrics import (cigar_op_bases, read_junctions,
                               REGIONS, EXCLUDE_MAIN)

FAIL = []


def check(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + ("" if cond else f"  [{detail}]"))
    if not cond:
        FAIL.append(name)


# ---------- A. 合成 CIGAR ----------
# 10M 5N 10M 3S 10M 2H 5M 4D 6M, reference_start = 122454290
ct = [(0, 10), (3, 5), (0, 10), (4, 3), (0, 10), (5, 2), (0, 5), (2, 4), (0, 6)]
cc = cigar_op_bases(ct)
check("A1 剪接 N 长度=5", cc["N"] == 5, str(cc))
check("A2 softclip S=3（H=2 不计入）", cc["S"] == 3 and cc["H"] == 2, str(cc))
check("A3 D=4", cc["D"] == 4, str(cc))
js = read_junctions(ct, 122454290)
# 10M 后进入 N：122454290+10=122454300 → 122454305
check("A4 唯一接头坐标 (122454300,122454305)", js == [(122454300, 122454305)], str(js))
check("A5 恰好 1 个接头", len(js) == 1, str(js))

# 多接头：5M 100N 5M 200N 5M
ct2 = [(0, 5), (3, 100), (0, 5), (3, 200), (0, 5)]
js2 = read_junctions(ct2, 1000)
check("A6 多接头坐标", js2 == [(1005, 1105), (1110, 1310)], str(js2))

# ---------- B. UTR 包含逻辑 ----------
u_start, u_end = REGIONS["annotated_UTR"][1], REGIONS["annotated_UTR"][2]
j_in = (122454304, 122455892)
j_out = (122454304, 122457302)   # 右端超出 UTR（122457153）
check("B1 122454304-122455892 两端均在 UTR 内",
      j_in[0] >= u_start and j_in[1] <= u_end)
check("B2 122454304-122457302 右端越界",
      not (j_out[0] >= u_start and j_out[1] <= u_end))

# ---------- C. 真实 BAM 回归（sham1） ----------
BAM = "/mnt/d/stroke_apa_reaudit/inputs/zip_2023/results/03_bam_review/slices/sham1.Atp2a2.cand.bam"
if os.path.exists(BAM):
    chrom, s, e = REGIONS["review_window"]
    reg = f"{chrom}:{s+1}-{e}"     # samtools 1-based
    sam_n = subprocess.run(
        ["samtools", "view", "-F", "0x904", BAM, reg],
        capture_output=True, text=True).stdout
    sam_spliced = sum(1 for line in sam_n.splitlines() if "N" in line.split("\t")[5])
    import pysam
    n_spliced = 0
    with pysam.AlignmentFile(BAM, "rb") as bam:
        for r in bam.fetch(chrom, s, e):
            if r.flag & EXCLUDE_MAIN:
                continue
            if any(op == 3 for op, _ in (r.cigartuples or [])):
                n_spliced += 1
    check("C1 sham1 剪接读段数 pysam==samtools", n_spliced == sam_spliced,
          f"pysam={n_spliced} samtools={sam_spliced}")
    print(f"    sham1 审查窗剪接读段数（修正口径）= {n_spliced}；旧脚本错误输出为 0")
else:
    check("C1 BAM 存在", False, BAM)

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("ALL TESTS PASSED")
