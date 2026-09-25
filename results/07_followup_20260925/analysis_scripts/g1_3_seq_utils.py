#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G1 复核 · 工作 3：修正版 FASTA 序列提取。

旧脚本 bug（勘误对象）：p4d_atp2a2_assay_design.py 自写的提取函数定位到目标所在行后，
没有跳过行内 start % line_bases 个碱基，导致行内非零偏移的窗口取到错误序列。

本实现：
  1) 主路径 pyfaidx（0-based half-open，与 BED 口径一致，库自身经过充分测试）；
  2) 独立对照 samtools faidx（1-based closed，换算后调用）；
  3) 逐碱基比对函数，供设计窗口全量对拍。
"""
import subprocess

from pyfaidx import Fasta

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def fetch_pyfaidx(fa_path, chrom, start0, end0, strand="+"):
    """BED 0-based half-open 提取；负链返回反向互补（转录方向 5'→3'）。"""
    fa = Fasta(fa_path)
    seq = str(fa[chrom][start0:end0]).upper()
    if strand == "-":
        seq = seq.translate(COMP)[::-1]
    return seq


def fetch_samtools(fa_path, chrom, start0, end0, strand="+"):
    """独立实现：samtools faidx（1-based closed）。"""
    r = subprocess.run(["samtools", "faidx", fa_path,
                        f"{chrom}:{start0+1}-{end0}"],
                       capture_output=True, text=True, check=True)
    seq = "".join(r.stdout.splitlines()[1:]).upper()
    if strand == "-":
        seq = seq.translate(COMP)[::-1]
    return seq


def verify_window(fa_path, chrom, start0, end0, strand="+"):
    """逐碱基对拍两种实现；不一致抛异常。返回统一序列。"""
    a = fetch_pyfaidx(fa_path, chrom, start0, end0, strand)
    b = fetch_samtools(fa_path, chrom, start0, end0, strand)
    if a != b:
        mism = [i for i, (x, y) in enumerate(zip(a, b)) if x != y][:10]
        raise AssertionError(f"序列不一致 {chrom}:{start0}-{end0}({strand}) "
                             f"len {len(a)} vs {len(b)}，首个错位 {mism}")
    return a
