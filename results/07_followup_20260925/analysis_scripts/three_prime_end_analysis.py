#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
three_prime_end_analysis.py — 独立 3' 端测序分析框架（Gate B）

状态：可运行框架。不包含任何假定结果；数据到位后按 3prime_end_analysis_plan.md 的
预定义规则填充。输入接口：
  --bam <per-animal 3' 端测序 BAM（已比对到 mm10/GRCm38，含软clip/adapter 已裁剪或标记）>
  --meta wetlab_metadata_template.tsv（每行一只动物）
  --out 输出目录
依赖：pysam, numpy（bioinfo 环境已有）。

分子归属规则（与 Atp2a2_transcript_structure.tsv 一致，BED 0-based）：
  L 接头 junction = (122454304, 122455892)   # donor122455892->acceptor122454304
  S 接头 junction = (122454304, 122457302)
  M 端            = 122456339 ± 50
  L/S 终端端      = 122453512/122453513 ± 50
有效分子：3' 端存在 non-templated 接合（ adapter/G 尾之前的 RNA-polyA 边界）。
"""
import argparse
from collections import defaultdict

L_JCT = (122454304, 122455892)
S_JCT = (122454304, 122457302)
END_M = 122456339
END_LS = 122453512   # L/S 共用终端端（±50 内视为该端）
TOL = 50


def classify_read(aln, junction_sets=None):
    """占位实现：输入一条比对，返回 (version, end0, valid, reason)。
    version ∈ L/S/M/None；valid=是否满足 RNA-polyA 接合判据。
    TODO(数据到位后)：non-templated 碱基从软clip/adapter 序列中判定。"""
    raise NotImplementedError("等待真实数据：接合判据需要 adapter 序列与 3' 端质量字段")


def per_animal_table(assignments):
    """assignments: list of (version, end0, valid) → 每位点末端分布与比例。"""
    by_version = defaultdict(list)
    for v, e, ok in assignments:
        if ok and v in ("L", "S", "M"):
            by_version[v].append(e)
    n = {v: len(x) for v, x in by_version.items()}
    total = sum(n.values())
    if total == 0:
        return {"n_L": 0, "n_S": 0, "n_M": 0, "ambiguous_pct": None,
                "L_pct": None, "ci95": None}
    L_pct = 100.0 * n.get("L", 0) / total
    # 置信区间：Beta 二项（Jeffreys 先验）
    import math
    a, b = n.get("L", 0) + 0.5, total - n.get("L", 0) + 0.5
    lo = max(0.0, L_pct - 1.96 * math.sqrt(L_pct * (100 - L_pct) / total))
    hi = min(100.0, L_pct + 1.96 * math.sqrt(L_pct * (100 - L_pct) / total))
    return {"n_L": n.get("L", 0), "n_S": n.get("S", 0), "n_M": n.get("M", 0),
            "L_pct": round(L_pct, 2), "ci95": (round(lo, 2), round(hi, 2))}


def group_effect(animal_tables, group_col):
    """动物级 ΔL%（day7 vs sham 主比较）；不得用细胞/读段数冒充 n。"""
    raise NotImplementedError("等待元数据到位后实现（含置换检验与 CI）")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", nargs="+", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    print("框架就绪；规则见 3prime_end_analysis_plan.md。classify_read 待数据实现。")
    # TODO: 遍历 BAM → classify_read → per_animal_table → group_effect


if __name__ == "__main__":
    main()
