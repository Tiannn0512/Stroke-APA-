#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G1 复核 · 工作 2：六样本同尺度局部覆盖深度图。
注意：本图为程序绘制的 DEPTH PLOT（覆盖深度曲线），不是 IGV 比对截图；
IGV 比对截图见 igv/ 目录（由 IGV 软件生成）。
"""
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BAMS = {
    "sham1":     "/mnt/d/stroke_apa_reaudit/inputs/zip_2023/results/03_bam_review/slices/sham1.Atp2a2.cand.bam",
    "sham2":     "/mnt/d/stroke_apa_reaudit/inputs/zip_2023/results/03_bam_review/slices/sham2.Atp2a2.cand.bam",
    "day3_rep1": "/mnt/d/stroke_apa_reaudit/inputs/zip_2023/results/03_bam_review/slices/day3_rep1.Atp2a2.cand.bam",
    "day3_rep2": "/mnt/d/stroke_apa_reaudit/inputs/zip_2023/results/03_bam_review/slices/day3_rep2.Atp2a2.cand.bam",
    "day7_rep1": "/mnt/d/stroke_apa_reaudit/inputs/zip_2023/results/03_bam_review/slices/day7_rep1.Atp2a2.cand.bam",
    "day7_rep2": "/mnt/d/stroke_apa_reaudit/inputs/zip_2023/results/03_bam_review/slices/day7_rep2.Atp2a2.cand.bam",
}
COLORS = {"sham1": "#2166ac", "sham2": "#67a9cf", "day3_rep1": "#ef8a62",
          "day3_rep2": "#f4a582", "day7_rep1": "#b2182b", "day7_rep2": "#d6604d"}
OUT = "/mnt/d/stroke_apa_reaudit/results/06_gate1_reaudit"

W0, W1 = 122452512, 122458153
PAS = 122456498
JUN = (122454304, 122455892)


def depth(bam):
    out = subprocess.run(
        ["samtools", "depth", "-a", "-Q", "20", "-q", "20",
         "-r", f"chr5:{W0+1}-{W1}", bam],
        capture_output=True, text=True, check=True).stdout
    xs, ys = [], []
    for line in out.splitlines():
        p = line.split("\t")
        xs.append(int(p[1]) - 1)
        ys.append(int(p[2]))
    return xs, ys


data = {s: depth(b) for s, b in BAMS.items()}
ymax = max(max(y) for _, y in data.values()) * 1.05

def decorate(ax, title):
    ax.axvspan(122453512, 122456498, color="0.85", zorder=0, label="distal 122453512-122456498")
    ax.axvspan(JUN[0], JUN[1], color="#fddbc7", zorder=0,
               label="179939 spliced-UTR intron 122454304-122455892")
    ax.axvline(PAS, color="k", ls="--", lw=1, label="predicted PAS 122456498")
    ax.set_xlim(W0, W1)
    ax.set_title(title, fontsize=9)

# 图 A：线性六样本
fig, ax = plt.subplots(figsize=(11, 4.5))
for s in sorted(data):
    ax.plot(data[s][0], data[s][1], color=COLORS[s], lw=1, label=s)
decorate(ax, "Atp2a2 event 6953 region — 6-sample coverage DEPTH PLOT (samtools depth -Q20 -q20; NOT an IGV alignment screenshot)")
ax.legend(fontsize=7, ncol=2)
ax.set_xlabel("chr5 (mm10, 0-based position)")
ax.set_ylabel("depth")
fig.tight_layout()
fig.savefig(f"{OUT}/depth_plot_6samples_linear.png", dpi=150)

# 图 B：对数尺
fig, ax = plt.subplots(figsize=(11, 4.5))
for s in sorted(data):
    ax.plot(data[s][0], data[s][1], color=COLORS[s], lw=1, label=s)
decorate(ax, "same, log scale")
ax.set_yscale("symlog", linthresh=10)
ax.legend(fontsize=7, ncol=2)
ax.set_xlabel("chr5 (mm10, 0-based position)")
ax.set_ylabel("depth (symlog)")
fig.tight_layout()
fig.savefig(f"{OUT}/depth_plot_6samples_symlog.png", dpi=150)

# 图 C：均值 ± 组
import statistics
xs = data["sham1"][0]
sham = [statistics.mean(v) for v in zip(data["sham1"][1], data["sham2"][1])]
d3 = [statistics.mean(v) for v in zip(data["day3_rep1"][1], data["day3_rep2"][1])]
d7 = [statistics.mean(v) for v in zip(data["day7_rep1"][1], data["day7_rep2"][1])]
fig, ax = plt.subplots(figsize=(11, 4.5))
ax.plot(xs, sham, color="#2166ac", lw=1.2, label="sham mean (n=2)")
ax.plot(xs, d3, color="#ef8a62", lw=1.2, label="day3 mean (n=2)")
ax.plot(xs, d7, color="#b2182b", lw=1.2, label="day7 mean (n=2)")
decorate(ax, "group means — DEPTH PLOT")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(f"{OUT}/depth_plot_groupmeans.png", dpi=150)
print("depth plots written")
