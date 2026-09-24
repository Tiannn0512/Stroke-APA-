#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G1 复核 · 工作 3 单元测试：FASTA 提取正确性（复现旧 bug 场景）。
构造行宽=10 和行宽=60 的合成 FASTA（含行内非零偏移起点与跨行区间），
断言 pyfaidx 提取 == samtools faidx 提取 == 已知序列；
再对真实 genome/transcripts 的设计窗口做逐碱基对拍。
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from g1_3_seq_utils import fetch_pyfaidx, fetch_samtools, verify_window

FAIL = []


def check(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + ("" if cond else f"  [{detail}]"))
    if not cond:
        FAIL.append(name)


def make_fasta(path, seq, width):
    with open(path, "w") as f:
        f.write(">chrT\n")
        for i in range(0, len(seq), width):
            f.write(seq[i:i + width] + "\n")


# 已知序列：401bp，含全部上下文
import random
random.seed(42)
KNOWN = "".join(random.choice("ACGT") for _ in range(401))

for width in (10, 60):
    fa = f"/tmp/test_w{width}.fa"
    make_fasta(fa, KNOWN, width)
    os.system(f"samtools faidx {fa}")
    # 场景1：行内非零偏移起点（旧 bug 触发条件）
    s = fetch_pyfaidx(fa, "chrT", 25, 90)
    check(f"W{width} 行内偏移起点 25-90", s == KNOWN[25:90], f"got {s[:20]}...")
    # 场景2：跨行整段
    s = fetch_pyfaidx(fa, "chrT", 0, 401)
    check(f"W{width} 全长跨行 0-401", s == KNOWN)
    # 场景3：两种实现逐碱基对拍（含负链）
    a = fetch_pyfaidx(fa, "chrT", 33, 77, "-")
    b = fetch_samtools(fa, "chrT", 33, 77, "-")
    comp = str.maketrans("ACGT", "TGCA")
    expect = KNOWN[33:77].translate(comp)[::-1]
    check(f"W{width} 负链反向互补", a == expect and b == expect)
    try:
        verify_window(fa, "chrT", 17, 200, "-")
        check(f"W{width} verify_window 通过", True)
    except AssertionError as e:
        check(f"W{width} verify_window 通过", False, str(e))

print()
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
print("ALL SYNTHETIC TESTS PASSED")

# ---------- 真实参考对拍 ----------
REAL = "/mnt/d/stroke_apa_reaudit/reference/chr5.fa"
if os.path.exists(REAL + ".fai"):
    wins = [
        ("chr5", 122453512, 122454304, "-"),   # 终端外显子 UTR
        ("chr5", 122455892, 122457153, "-"),   # 倒数第二外显子 UTR
        ("chr5", 122454262, 122454354, "-"),   # 跨 177974 UTR 边界 122454293
        ("chr5", 122456420, 122456620, "-"),   # 跨预测 PAS 122456498
        ("chr5", 122457240, 122457480, "-"),   # 跨三个转录本共享外显子边界
    ]
    ok = True
    for c, s, e, st in wins:
        try:
            verify_window(REAL, c, s, e, st)
            print(f"PASS genome 对拍 {c}:{s}-{e}({st})")
        except AssertionError as ex:
            ok = False
            print(f"FAIL genome 对拍 {c}:{s}-{e}: {ex}")
    sys.exit(0 if ok else 1)
else:
    print("genome.fa 尚未就绪，真实对拍待 genome 下载完成后运行")
