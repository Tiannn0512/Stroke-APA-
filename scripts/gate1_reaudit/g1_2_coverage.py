#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G1 复核 · 工作 2：分区覆盖与替代解释（独立实现，不依赖旧 p3_coverage_stats）

坐标口径 BED 0-based half-open。深度口径：samtools depth -a -Q 20 -q 20（与旧管线一致），
另用 pysam.count_coverage 独立复算 sham1 交叉验证。

区段（两类）：
  canonical（预设分区）:
    distal_full            122453512-122456498
    overlap_177974_UTR     122453512-122454293   （与 177974 的 3'UTR 重叠）
    exclusive_3UTR_part    122454293-122456498   （避开重叠的部分）
    shared                 122456498-122457153
  structural（本轮 GTF 结构拆分，解释用）:
    terminal_exon          122453512-122454304   （177974+179939 共用终端外显子，含 UTR）
    spliced_UTR_intron     122454304-122455892   （179939 跨-UTR 内含子，成熟 mRNA 无覆盖）
    penult_UTR_5p_part     122455892-122456498   （179939 倒数第二外显子 UTR 前段）

敏感性：排除含主接头 122454304-122455892（=179939 成熟分子）的读段后重算各区段
平均深度 → 拆分"179939 成熟分子贡献"与"其余（pre-mRNA/177974/31423 等）"。
"""
import argparse
import subprocess
import json

import pysam

SEGMENTS = [
    ("distal_full",         122453512, 122456498),
    ("overlap_177974_UTR",  122453512, 122454293),
    ("exclusive_3UTR_part", 122454293, 122456498),
    ("shared",              122456498, 122457153),
    ("terminal_exon",       122453512, 122454304),
    ("spliced_UTR_intron",  122454304, 122455892),
    ("penult_UTR_5p_part",  122455892, 122456498),
]
KEY_JUNCTION = (122454304, 122455892)
CHROM, WSTART, WEND = "chr5", 122452512, 122458153


def depth_by_samtools(bam, chrom, start, end, mapq=20, bq=20):
    """samtools depth，返回 {pos0: depth}。"""
    out = subprocess.run(
        ["samtools", "depth", "-a", "-Q", str(mapq), "-q", str(bq),
         "-r", f"{chrom}:{start+1}-{end}", bam],
        capture_output=True, text=True, check=True).stdout
    d = {}
    for line in out.splitlines():
        p = line.split("\t")
        d[int(p[1]) - 1] = int(p[2])     # 1-based → 0-based
    return d


def depth_by_pysam(bam, chrom, start, end, mapq=20, bq=20):
    """pysam 独立实现（交叉验证用），返回 {pos0: depth}。"""
    with pysam.AlignmentFile(bam, "rb") as f:
        cov = f.count_coverage(chrom, start, end, read_callback="all",
                               quality_threshold=bq)
    total = [a + b + c + d for a, b, c, d in zip(*cov)]
    # count_coverage 无 MAPQ 过滤；仅用于与 samtools 口径对拍 MAPQ=0 人群可忽略的样本
    return {start + i: v for i, v in enumerate(total)}


def segment_stats(depth, s, e):
    vals = [depth.get(p, 0) for p in range(s, e)]
    mean = sum(vals) / len(vals) if vals else 0.0
    covered = sum(1 for v in vals if v >= 1)
    return mean, covered, (e - s)


def reads_with_key_junction(bam_path):
    """返回含主接头读段的 (query_name, flag, positions) 集合要素：这里用 read迭代
    时直接判断，不落盘名字集合（避免同名歧义），改为返回谓词所需的接头集合。"""
    marked = set()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for r in bam.fetch(CHROM, WSTART, WEND):
            if r.flag & (0x4 | 0x100 | 0x800):
                continue
            ref = r.reference_start
            for op, ln in (r.cigartuples or []):
                if op == 3 and (ref, ref + ln) == KEY_JUNCTION:
                    marked.add((r.query_name, r.flag, r.reference_start))
                if op in (0, 2, 3, 7, 8):
                    ref += ln
    return marked


def depth_excluding_marked(bam, chrom, start, end, marked, mapq=20, bq=20):
    """与 samtools depth 同口径（MAPQ/BQ 过滤），但排除 marked 读段的逐碱基深度。"""
    d = {}
    with pysam.AlignmentFile(bam, "rb") as f:
        for r in f.fetch(chrom, start, end):
            if r.flag & (0x4 | 0x100 | 0x800):
                continue
            if (r.query_name, r.flag, r.reference_start) in marked:
                continue
            ref = r.reference_start
            qseq = r.query_sequence
            qqual = r.query_qualities
            ref = r.reference_start
            qi = 0
            for op, ln in (r.cigartuples or []):
                if op in (0, 7, 8):          # M = X
                    for k in range(ln):
                        rp, qp = ref + k, qi + k
                        if start <= rp < end and qqual is not None and qqual[qp] >= bq:
                            d[rp] = d.get(rp, 0) + 1
                    ref += ln
                    qi += ln
                elif op == 1:                # I
                    qi += ln
                elif op == 4:                # S
                    qi += ln
                elif op in (2, 3):           # D N
                    ref += ln
                # H 不消耗
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bams", nargs="+", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    bams = dict(kv.split("=", 1) for kv in a.bams)
    samples = sorted(bams)

    # ---- 主表：samtools depth 口径 ----
    rows = []
    seg_depth = {}
    for s in samples:
        d = depth_by_samtools(bams[s], CHROM, WSTART, WEND)
        shared_mean, _, _ = segment_stats(d, 122456498, 122457153)
        seg_depth[s] = (d, shared_mean)
        for name, sg, se in SEGMENTS:
            m, cov, ln = segment_stats(d, sg, se)
            rows.append((s, name, sg, se, ln, f"{m:.3f}",
                         f"{100.0*cov/ln:.2f}", f"{m/shared_mean:.4f}"))
    with open(f"{a.outdir}/segment_coverage_by_sample.tsv", "w") as f:
        f.write("sample\tsegment\tstart0\tend0\tlength_bp\tmean_depth\t"
                "covered_bases_pct\tratio_vs_shared\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")

    # ---- 交叉验证：sham1 distal_full，pysam vs samtools ----
    d_s = seg_depth[samples[0]][0] if samples else {}
    d_p = depth_by_pysam(bams[samples[0]], CHROM, WSTART, WEND)
    m_s, _, _ = segment_stats(d_s, 122453512, 122456498)
    m_p, _, _ = segment_stats(d_p, 122453512, 122456498)
    xcheck = {"sample": samples[0], "samtools_mean_distal": round(m_s, 4),
              "pysam_mean_distal_no_mapq": round(m_p, 4),
              "note": "pysam count_coverage 无 MAPQ 过滤，仅作量级对拍；"
                      "MAPQ<20 比例见 corrected_read_metrics（<1%）"}

    # ---- 敏感性：排除含主接头（179939 成熟分子）的读段 ----
    sens = []
    for s in samples:
        marked = reads_with_key_junction(bams[s])
        d0, shared_mean0 = seg_depth[s]
        d1 = depth_excluding_marked(bams[s], CHROM, WSTART, WEND, marked)
        shared_mean1, _, _ = segment_stats(d1, 122456498, 122457153)
        row = {"sample": s, "marked_reads": len(marked),
               "shared_mean_full": round(shared_mean0, 3),
               "shared_mean_excl": round(shared_mean1, 3)}
        for name, sg, se in SEGMENTS:
            m0, _, _ = segment_stats(d0, sg, se)
            m1, _, _ = segment_stats(d1, sg, se)
            row["ratio_" + name] = round(m0 / shared_mean0, 4) if shared_mean0 else None
            row["ratio_exclJ_" + name] = round(m1 / shared_mean1, 4) if shared_mean1 else None
        sens.append(row)
    with open(f"{a.outdir}/coverage_sensitivity_data.tsv", "w") as f:
        keys = list(sens[0].keys())
        f.write("\t".join(keys) + "\n")
        for r in sens:
            f.write("\t".join(str(r[k]) for k in keys) + "\n")

    with open(f"{a.outdir}/crosscheck_sham1.json", "w") as f:
        json.dump(xcheck, f, ensure_ascii=False, indent=1)
    print(json.dumps(xcheck, ensure_ascii=False))


if __name__ == "__main__":
    main()
