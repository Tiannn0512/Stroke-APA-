#!/usr/bin/env python3
"""P3.2: fixed-scale per-gene coverage plots (6 samples) + IGV snapshot command sheet.

One PNG per gene: per-base depth (-Q20 -q20) for sham1/2, day3_rep1/2, day7_rep1/2 on a
common y-scale (max across samples over the region); UTR span shaded; predicted PAS and
distal/shared boundary marked; QKI peaks (if inside region) drawn as vspans.
Outputs: results/03_bam_review/plots/<gene>.png, igv_snapshot_commands.txt
"""
import csv, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = r"D:\stroke_apa_reanalysis"
BR = os.path.join(BASE, "results", "03_bam_review")
PLOTS = os.path.join(BR, "plots")
os.makedirs(PLOTS, exist_ok=True)
SAMPLES = ["sham1", "sham2", "day3_rep1", "day3_rep2", "day7_rep1", "day7_rep2"]
COLORS = {"sham1": "#4c78a8", "sham2": "#4c78a8", "day3_rep1": "#e45756", "day3_rep2": "#e45756",
          "day7_rep1": "#2ca02c", "day7_rep2": "#2ca02c"}

regions = {}
with open(os.path.join(BR, "regions.bed"), encoding="utf-8") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        regions[p[3].split("|")[0]] = (p[0], int(p[1]) + 5000, int(p[2]) - 5000)

seg = {}
for r in csv.DictReader(open(os.path.join(BASE, "results", "02_candidate_rebuild", "candidate_ranking_v2.tsv"), encoding="utf-8"), delimiter="\t"):
    if r["gene_symbol"] in regions:
        seg[r["gene_symbol"]] = r

qki = []
for line in open(os.path.join(BASE, "input_links", "p3_level2_qki_clip_peaks.bed"), encoding="utf-8"):
    p = line.rstrip("\n").split("\t")
    qki.append((p[0], int(p[1]), int(p[2]), p[3]))

depth = {g: {s: {} for s in SAMPLES} for g in regions}
for s in SAMPLES:
    with open(os.path.join(BR, f"depth_{s}.tsv"), encoding="utf-8") as f:
        for line in f:
            chrom, pos, v = line.rstrip("\n").split("\t")
            pos = int(pos)
            for g, (c, us, ue) in regions.items():
                if chrom == c and us <= pos < ue:
                    depth[g][s][pos] = int(v)

cmds = ["# IGV batch — paste into IGV (View > Preferences: mm10 genome loaded first) or run via Tools > Run Batch",
        "# one block per review gene; fixed scale = max across 6 samples over UTR+-5kb"]
for g, (chrom, us, ue) in sorted(regions.items()):
    r = seg[g]
    pas = int(r["proximal_pas_bed"])
    ymax = max((max(depth[g][s].values(), default=1) for s in SAMPLES), default=1)
    ymax = max(ymax, 10)
    fig, ax = plt.subplots(figsize=(12, 5))
    for s in SAMPLES:
        xs = sorted(depth[g][s])
        ax.plot(xs, [depth[g][s][x] for x in xs], color=COLORS[s], lw=0.9, alpha=0.85,
                label=s if s in ("sham1", "day3_rep1", "day7_rep1") else None)
    d = list(map(int, r["distal"].split("-")))
    sh = list(map(int, r["shared"].split("-")))
    ax.axvspan(d[0], d[1], color="#bbb", alpha=0.25)
    ax.axvline(pas, color="k", ls="--", lw=1)
    ax.axvline(us, color="gray", ls=":", lw=0.8); ax.axvline(ue, color="gray", ls=":", lw=0.8)
    for (pc, ps, pe, pg) in qki:
        if pc == chrom and pe > us and ps < ue:
            ax.axvspan(ps, pe, color="#f5a623", alpha=0.4)
            ax.text((ps+pe)/2, ymax*0.95, "QKI", ha="center", fontsize=8, color="#8a5b00")
    ax.set_xlim(us, ue); ax.set_ylim(0, ymax)
    ax.set_title(f"{g} ({r['chrom']}, {r['strand']}) class={r['class_v2']}  distal:{r['distal']}  PAS:{pas}\n"
                 f"scale 0-{ymax}; gray band=distal segment; dashed=PAS; orange=QKI peaks")
    ax.set_ylabel(f"depth (MQ20/BQ20)")
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, f"{g}.png"), dpi=150)
    plt.close(fig)
    cmds += [f"# ---- {g} ({chrom}) ----",
             f"goto {chrom}:{us}-{ue}",
             f"maxPanelHeight 400;  # set data range manually to 0-{ymax} (right-click track)",
             f"snapshot {PLOTS.replace(chr(92), '/')}/igv_{g}.png"]

with open(os.path.join(BR, "igv_snapshot_commands.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(cmds) + "\n")
print(f"plots written: {len(seg)}; commands: igv_snapshot_commands.txt")
