#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
f2: 复算并整合 6 个 Atp2a2 BAM 切片的局部读段证据 →
  Atp2a2_local_read_evidence.tsv（长表，一行一样本×指标）

过滤规则（写入 filter_rule 列）：
  ruleA_primary : 排除 unmapped(0x4)/secondary(0x100)/supplementary(0x800)；保留重复读段与全部 MAPQ
  ruleA_dedup   : 同上再排除 duplicate(0x400)
  ruleB_depth20 : samtools depth -a -Q20 -q20（MAPQ>=20 且碱基Q>=20；含重复读段）
  junction计数  : CIGAR N 与目标内含子精确相等的 primary 读段（ruleA 口径）

分区（BED 0-based half-open）与类别：
  terminal_exon       122453512-122454304  exon(LS 共用终端外显子)
  spliced_UTR_intron  122454304-122455892  intron(仅 pre-mRNA；L 独有内含子)
  penult_UTR_5p_part  122455892-122456498  exon(L 独有段)
  LM_segment          122456339-122457302  exon(LM 共有)
  shared_penult       122457302-122457426  exon(LSM 共有)
  common_exon_511     122457511-122457629  exon(LSM 共有)
  intergenic_control  122452512-122453012  intergenic(基因 3' 端外侧 500bp)
敏感性：排除含 L 独有接头 122454304-122455892 的读段后重算（ruleB_depth20 口径）。
"""
import subprocess
import sys

PY = sys.executable
SCRIPTS = "/mnt/d/stroke_apa_reaudit/results/07_followup_20260925/analysis_scripts"  # 包内自洽(g1 脚本已随包)
OUT = "/mnt/d/stroke_apa_reaudit/results/07_followup_20260925"
SLICES = "/mnt/d/stroke_apa_reaudit/inputs/zip_2023/results/03_bam_review/slices"
SAMPLES = ["sham1", "sham2", "day3_rep1", "day3_rep2", "day7_rep1", "day7_rep2"]
GROUP = {"sham1": "sham", "sham2": "sham", "day3_rep1": "day3", "day3_rep2": "day3",
         "day7_rep1": "day7", "day7_rep2": "day7"}
BAMS = {s: f"{SLICES}/{s}.Atp2a2.cand.bam" for s in SAMPLES}

# ---- 1) 复跑 g1_1 / g1_2（与 06 包同一脚本，输出到新版本目录） ----
bam_arg = " ".join(f"{s}={BAMS[s]}" for s in SAMPLES)
for cmd in (
    f"{PY} {SCRIPTS}/g1_1_read_metrics.py --bams {bam_arg} --outdir {OUT}",
    f"{PY} {SCRIPTS}/g1_2_coverage.py --bams {bam_arg} --outdir {OUT}",
):
    print(">>", cmd[:120])
    subprocess.run(cmd, shell=True, check=True)

# ---- 2) 整合 ----
import pysam

KEY_JUNCTION = (122454304, 122455892)     # L 独有
SEGMENTS = [
    ("terminal_exon", 122453512, 122454304, "exon(LS)"),
    ("spliced_UTR_intron", 122454304, 122455892, "intron(L-unique;pre-mRNA)"),
    ("penult_UTR_5p_part", 122455892, 122456498, "exon(L-unique)"),
    ("LM_segment", 122456339, 122457302, "exon(LM)"),
    ("shared_penult", 122457302, 122457426, "exon(LSM)"),
    ("common_exon_511", 122457511, 122457629, "exon(LSM)"),
    ("intergenic_control", 122452512, 122453012, "intergenic"),
]

def depth_mapq20(bam, chrom, s, e):
    out = subprocess.run(["samtools", "depth", "-a", "-Q", "20", "-q", "20",
                          "-r", f"{chrom}:{s+1}-{e}", bam],
                         capture_output=True, text=True, check=True).stdout
    d = {}
    for line in out.splitlines():
        p = line.split("\t")
        d[int(p[1]) - 1] = int(p[2])
    return d

rows = []
# 2a. intergenic control 深度
for s in SAMPLES:
    d = depth_mapq20(BAMS[s], "chr5", 122452512, 122453012)
    m = sum(d.get(p, 0) for p in range(122452512, 122453012)) / 500.0
    rows.append([s, GROUP[s], "mean_depth", "intergenic_control", 122452512, 122453012,
                 "intergenic", "ruleB_depth20", f"{m:.3f}", "x"])

# 2b. 主指标合并
def read_tsv(path):
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            yield dict(zip(hdr, line.rstrip("\n").split("\t")))

for r in read_tsv(f"{OUT}/segment_coverage_by_sample.tsv"):
    s = r["sample"]
    seg = r["segment"]
    cls = next((c for n, _, _, c in SEGMENTS if n == seg), "composite(canonical)")
    st, en = r["start0"], r["end0"]
    rows.append([s, GROUP[s], "mean_depth_ratio_vs_shared", seg, st, en, cls,
                 "ruleB_depth20", r["ratio_vs_shared"], "ratio"])
    rows.append([s, GROUP[s], "mean_depth", seg, st, en, cls,
                 "ruleB_depth20", r["mean_depth"], "x"])
    rows.append([s, GROUP[s], "covered_bases_pct", seg, st, en, cls,
                 "ruleB_depth20", r["covered_bases_pct"], "%"])

for r in read_tsv(f"{OUT}/coverage_sensitivity_data.tsv"):
    s = r["sample"]
    rows.append([s, GROUP[s], "keyJ_read_count", "spliced_UTR_intron(L-junction)",
                 KEY_JUNCTION[0], KEY_JUNCTION[1], "intron(L-unique;pre-mRNA)",
                 "ruleA_primary", r["marked_reads"], "reads"])
    rows.append([s, GROUP[s], "distal_ratio_full_vs_exclJ",
                 "distal_full", 122453512, 122456498, "composite(canonical)",
                 "ruleB_depth20; excl-J sensitivity",
                 f"{r['ratio_distal_full']} -> {r['ratio_exclJ_distal_full']}", "ratio"])

for r in read_tsv(f"{OUT}/corrected_read_metrics_main.tsv"):
    s = r["sample"]
    rows.append([s, GROUP[s], f"reads_{r['region']}", r["region"], "", "",
                 "composite", "ruleA_primary", r["n_reads"], "reads"])
    rows.append([s, GROUP[s], f"spliced_pct_{r['region']}", r["region"], "", "",
                 "composite", "ruleA_primary", r["spliced_read_pct"], "%"])

# 2c. 关键接头逐样本计数（L 独有 / S 独有 / 共有内含子#2-#4）
S_JUNCTION = (122454304, 122457302)
COMMON_INTRONS = [(122455892, 122456498), (122456498, 122457153), (122457153, 122457302)]
def junction_reads(bam, j):
    n = 0
    with pysam.AlignmentFile(bam, "rb") as f:
        for r in f.fetch("chr5", j[0] - 2000, j[1] + 2000):
            if r.flag & (0x4 | 0x100 | 0x800):
                continue
            ref = r.reference_start
            for op, ln in (r.cigartuples or []):
                if op == 3 and (ref, ref + ln) == j:
                    n += 1
                if op in (0, 2, 3, 7, 8):
                    ref += ln
    return n

for s in SAMPLES:
    for j, lab, cls in [
        (KEY_JUNCTION, "junction_L_122454304_122455892", "junction(L-unique)"),
        (S_JUNCTION, "junction_S_122454304_122457302", "junction(S-unique)"),
    ]:
        n = junction_reads(BAMS[s], j)
        rows.append([s, GROUP[s], "spliced_reads", lab, j[0], j[1], cls,
                     "ruleA_primary", n, "reads"])

with open(f"{OUT}/Atp2a2_local_read_evidence.tsv", "w") as out:
    out.write("sample\tgroup\tmetric\tregion\tstart0\tend0\tregion_class\tfilter_rule\tvalue\tunit\n")
    for r in rows:
        out.write("\t".join(map(str, r)) + "\n")
print("written:", f"{OUT}/Atp2a2_local_read_evidence.tsv", len(rows), "rows")
