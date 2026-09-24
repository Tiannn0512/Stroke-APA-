#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G1 复核 · 工作 1：修正版 CIGAR / 剪接读段审查（替代有 bug 的 p3c_read_metrics.sh）

旧脚本缺陷（勘误对象）：CIGAR 条件里写了 `op="="`（赋值而非比较），导致
spliced_read_pct 恒为 0。本脚本改用 pysam 的 cigartuples 显式枚举。

坐标口径：BED 0-based half-open（与任务书一致）。
过滤口径：
  main  = 排除 unmapped(0x4) + secondary(0x100) + supplementary(0x800)
  dedup = main 基础上再排除 duplicate(0x400)（敏感性口径）
计数单位：read（一条 primary alignment 记 1），不使用 fragment 计数。
soft-clipped：只计 S 的碱基；H 不作为测得的 soft-clipped 碱基。
"""
import argparse
import statistics
from collections import defaultdict

import pysam

REGIONS = {
    "annotated_UTR": ("chr5", 122453512, 122457153),
    "review_window": ("chr5", 122452512, 122458153),  # UTR 两侧各 1 kb
}
EXCLUDE_MAIN = 0x4 | 0x100 | 0x800
EXCLUDE_DEDUP = EXCLUDE_MAIN | 0x400
QUERY_CONSUMING = (0, 1, 4, 7, 8)   # M I S = X
REF_CONSUMING = (0, 2, 3, 7, 8)     # M D N = X


def cigar_op_bases(cigartuples):
    """按操作类型返回消耗碱基字典；H 单列，不计入 soft-clipped。"""
    d = {k: 0 for k in ("M", "I", "D", "N", "S", "H", "EQ", "X")}
    for op, ln in cigartuples:
        d[{0: "M", 1: "I", 2: "D", 3: "N", 4: "S", 5: "H", 7: "EQ", 8: "X"}[op]] += ln
    return d


def read_junctions(cigartuples, reference_start):
    """返回 [(start0, end0), ...]：每个 N 块的参考区间（即剪接内含子坐标）。"""
    out, ref = [], reference_start
    for op, ln in cigartuples:
        if op == 3:
            out.append((ref, ref + ln))
        if op in REF_CONSUMING:
            ref += ln
    return out


def sample_region_metrics(bam_path, chrom, start, end, exclude_flag):
    n = 0
    mapqs, nms = [], []
    sc_reads = sc_bases = query_bases = n_spliced = 0
    junctions = defaultdict(int)
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(chrom, start, end):
            if read.flag & exclude_flag:
                continue
            n += 1
            mapqs.append(read.mapping_quality)
            if read.has_tag("NM"):
                nms.append(int(read.get_tag("NM")))
            ct = read.cigartuples or []
            cc = cigar_op_bases(ct)
            sc_bases += cc["S"]
            query_bases += sum(cc[k] for k in ("M", "I", "S", "EQ", "X"))
            has_s = cc["S"] > 0
            if has_s:
                sc_reads += 1
            if cc["N"] > 0:
                n_spliced += 1
                for j in read_junctions(ct, read.reference_start):
                    junctions[j] += 1
    return {
        "n_reads": n,
        "mapq_median": statistics.median(mapqs) if mapqs else None,
        "mapq_lt20_pct": 100.0 * sum(1 for m in mapqs if m < 20) / n if n else None,
        "mapq_255_pct": 100.0 * sum(1 for m in mapqs if m == 255) / n if n else None,
        "mapq_ge30_pct": 100.0 * sum(1 for m in mapqs if m >= 30) / n if n else None,
        "softclip_read_pct": 100.0 * sc_reads / n if n else None,
        "softclip_base_pct": 100.0 * sc_bases / query_bases if query_bases else None,
        "spliced_read_pct": 100.0 * n_spliced / n if n else None,
        "nm_median": statistics.median(nms) if nms else None,
        "junctions": junctions,
    }


def collect(bam_paths, utr, exclude_flag):
    chrom, us, ue = utr["annotated_UTR"]
    per_sample = {}
    for name, path in bam_paths.items():
        per_sample[name] = {
            region: sample_region_metrics(path, c, s, e, exclude_flag)
            for region, (c, s, e) in
            list(utr.items()) + [("review_window", REGIONS["review_window"])]
        }
    return per_sample


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bams", nargs="+", required=True,
                    help="sample=path 形式，6 个 BAM 切片")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    bam_paths = dict(kv.split("=", 1) for kv in a.bams)
    utr = {"annotated_UTR": REGIONS["annotated_UTR"]}
    utr_start, utr_end = REGIONS["annotated_UTR"][1], REGIONS["annotated_UTR"][2]

    for tag, excl in (("main", EXCLUDE_MAIN), ("dedup", EXCLUDE_DEDUP)):
        res = collect(bam_paths, utr, excl)
        # ---- 指标表 ----
        out = f"{a.outdir}/corrected_read_metrics_{tag}.tsv"
        with open(out, "w") as f:
            f.write("sample\tregion\tfilter\tn_reads\tmapq_median\tmapq_lt20_pct\t"
                    "mapq_255_pct\tmapq_ge30_pct\tsoftclip_read_pct\tsoftclip_base_pct\t"
                    "spliced_read_pct\tnm_median\n")
            for s in sorted(res):
                for region in ("annotated_UTR", "review_window"):
                    m = res[s][region]
                    f.write(f"{s}\t{region}\t{tag}\t{m['n_reads']}\t{m['mapq_median']}\t"
                            f"{fmt(m['mapq_lt20_pct'])}\t{fmt(m['mapq_255_pct'])}\t"
                            f"{fmt(m['mapq_ge30_pct'])}\t{fmt(m['softclip_read_pct'])}\t"
                            f"{fmt(m['softclip_base_pct'])}\t{fmt(m['spliced_read_pct'])}\t"
                            f"{m['nm_median']}\n")
        # ---- 接头表（review_window 口径，含 UTR 内外拆分）----
        if tag == "main":
            jout = f"{a.outdir}/junction_counts_by_sample.tsv"
            all_j = set()
            per = {}
            for s in sorted(res):
                js = res[s]["review_window"]["junctions"]
                per[s] = js
                all_j |= set(js)
            samples = sorted(res)
            with open(jout, "w") as f:
                head = ["junction_start0", "junction_end0", "both_ends_in_UTR",
                        "total"] + [f"{s}_main" for s in samples] + [f"{s}_dedup" for s in samples]
                f.write("\t".join(head) + "\n")
                dedup_j = {
                    s: res2[s]["review_window"]["junctions"]
                    for s in samples
                    for res2 in [collect({s: bam_paths[s]}, utr, EXCLUDE_DEDUP)]
                }
                for j in sorted(all_j):
                    s0, e0 = j
                    both = "yes" if (s0 >= utr_start and e0 <= utr_end) else "no"
                    row = [s0, e0, both, sum(per[s].get(j, 0) for s in samples)]
                    row += [per[s].get(j, 0) for s in samples]
                    row += [dedup_j[s].get(j, 0) for s in samples]
                    f.write("\t".join(map(str, row)) + "\n")
    print("written:", a.outdir)


def fmt(x):
    return "" if x is None else f"{x:.3f}"


if __name__ == "__main__":
    main()
