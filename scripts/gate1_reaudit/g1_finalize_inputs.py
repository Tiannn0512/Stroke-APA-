#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""收尾机械件：合并 metrics 表、修 edgeR CSV 冒充 TSV、生成 input_manifest + software_versions"""
import csv
import hashlib
import os

R = "/mnt/d/stroke_apa_reaudit/results/06_gate1_reaudit"
INP = "/mnt/d/stroke_apa_reaudit/inputs"
REF = "/mnt/d/stroke_apa_reaudit/reference"


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


# 1) 合并 metrics
with open(f"{R}/corrected_read_metrics.tsv", "w") as out:
    for tag in ("main", "dedup"):
        with open(f"{R}/corrected_read_metrics_{tag}.tsv") as f:
            head = f.readline()
            if tag == "main":
                out.write(head)
            for line in f:
                out.write(line)
print("corrected_read_metrics.tsv 合并完成")

# 2) edgeR 文件：确认分隔符并以真正 TSV 重新导出（保留旧件）
src = f"{INP}/zip_2023/results/04_public_pap/edgeR_TMM_PAP_vs_Full.tsv"
with open(src) as f:
    head = f.readline()
sep = "," if head.count(",") > head.count("\t") else "\t"
print(f"旧 edgeR 文件分隔符 = {sep!r}（任务书判定：{',实为 CSV' if sep==',' else '已是 TSV'}）")
with open(src) as f, open(f"{R}/edgeR_TMM_PAP_vs_Full.tsv", "w") as out:
    for line in f:
        out.write(line.replace(sep, "\t"))
print("已重新导出真正 TSV 至 results/06_gate1_reaudit/（旧件保留于 zip_2023 原位）")
# 顺带提取 Atp2a2 行供决策文书引用
with open(f"{R}/edgeR_TMM_PAP_vs_Full.tsv") as f:
    for line in f:
        if line.startswith("ENSMUSG00000029467"):
            print("Atp2a2 行:", line.strip())

# 3) input_manifest.tsv
items = [
    ("REANALYSIS_RESULTS_20260923.zip", f"{INP}/REANALYSIS_RESULTS_20260923.zip",
     "GitHub Release v2026-09-23 资产（字节数与远端登记一致）"),
    ("REANALYSIS_RESULTS_20260921.zip", f"{INP}/REANALYSIS_RESULTS_20260921.zip",
     "GitHub Release v2026-09-23 资产（含 120 个 BAM 切片）"),
    ("gencode.vM25.annotation.gtf.gz", f"{REF}/gencode.vM25.annotation.gtf.gz",
     "GENCODE 官网 release_M25"),
    ("gencode.vM25.transcripts.fa.gz", f"{REF}/gencode.vM25.transcripts.fa.gz",
     "GENCODE 官网 release_M25"),
    ("chr5.fa.gz", f"{REF}/chr5.fa.gz",
     "UCSC mm10 chromosomes/chr5.fa.gz（md5 53e996b8236720821f82a687423edd22 与官方一致；"
     "mm10=GRCm38 同一组装，chr 命名与 BAM 一致）"),
]
for s in ("sham1", "sham2", "day3_rep1", "day3_rep2", "day7_rep1", "day7_rep2"):
    items.append((f"{s}.Atp2a2.cand.bam",
                  f"{INP}/zip_2023/results/03_bam_review/slices/{s}.Atp2a2.cand.bam",
                  "0921/0923 双包内均有，字节数一致"))
items += [
    ("p2_full_pdui_matrix.tsv", f"{INP}/zip_2023/input_links/p2_full_pdui_matrix.tsv",
     "冻结输入（input_links sha256 清单在案）"),
    ("p2_sample_map.tsv", f"{INP}/zip_2023/input_links/p2_sample_map.tsv", "冻结输入"),
    ("Atp2a2_read_metrics.tsv(旧,含spliced=0勘误)",
     f"{INP}/zip_2023/results/03_bam_review/Atp2a2_read_metrics.tsv", "旧交付，本轮勘误对象"),
    ("event_review_table.v2.tsv(旧)",
     f"{INP}/zip_2023/results/03_bam_review/event_review_table.v2.tsv", "旧交付，本轮更新为 v3"),
    ("decision_summary_v1.md(旧)", f"{INP}/zip_2023/results/05_target_nomination/decision_summary_v1.md",
     "旧交付，本轮更新为 v2"),
    ("edgeR_TMM_PAP_vs_Full.tsv(旧,CSV冒名)", src, "旧交付，本轮重导出真 TSV"),
    ("Atp2a2_assay_design.tsv(旧)", f"{INP}/zip_2023/results/05_target_nomination/Atp2a2_assay_design.tsv",
     "旧交付（其序列提取有 bug），本轮重设计"),
]
with open(f"{R}/input_manifest.tsv", "w") as f:
    f.write("file\tsha256\tbytes\tsource_note\n")
    for name, p, note in items:
        if os.path.exists(p):
            f.write(f"{name}\t{sha256(p)}\t{os.path.getsize(p)}\t{note}\n")
        else:
            f.write(f"{name}\tMISSING\t\t{note}\n")
            print("!! 缺失:", name)
print("input_manifest.tsv 完成")

# 4) software_versions.txt
import subprocess
def run(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout.strip().split("\n")[0]
pyv = run("python --version")
pysam = run("python -c \"import pysam; print(pysam.__version__)\"")
p3 = run("python -c \"import primer3; print('primer3-py', primer3.__version__)\"")
pf = run("python -c \"import pyfaidx; print('pyfaidx', pyfaidx.__version__)\"")
mpl = run("python -c \"import matplotlib; print('matplotlib', matplotlib.__version__)\"")
sam = run("samtools --version | head -1")
with open(f"{R}/software_versions.txt", "w") as f:
    f.write("""# Gate1 复核环境（2026-09-24，WSL Ubuntu + conda env bioinfo）
OS: Windows 11 + WSL2 Ubuntu (WSL 2.4.13.0)
""")
    for line in (pyv, pysam, p3, pf, mpl, sam):
        f.write(line + "\n")
    f.write("IGV: 2.18.4（Windows，系统 OpenJP JDK：java 17.0.20.1 运行）\n")
    f.write("参考: UCSC mm10 chr5.fa (md5 53e996b8236720821f82a687423edd22) + GENCODE vM25 GTF\n")
print("software_versions.txt 完成")
