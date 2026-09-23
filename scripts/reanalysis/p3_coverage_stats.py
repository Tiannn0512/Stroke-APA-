#!/usr/bin/env python3
"""P3.3: per-sample distal/shared coverage statistics for the review set.

Reads samtools depth output (-a -Q20 -q20 over regions.bed = UTR+-5kb per event).
Per event x sample: mean depth over distal / shared segment, coverage_ratio
(distal_mean / shared_mean), covered-base fraction at >=1/5/10.
Output: results/03_bam_review/bam_coverage_review.tsv
"""
import csv, os

BASE = r"D:\stroke_apa_reanalysis"
BR = os.path.join(BASE, "results", "03_bam_review")
SAMPLES = ["sham1", "sham2", "day3_rep1", "day3_rep2", "day7_rep1", "day7_rep2"]

# regions: gene -> (chrom, utr_start, utr_end)
regions = {}
with open(os.path.join(BR, "regions.bed"), encoding="utf-8") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        gene = p[3].split("|")[0]
        regions[gene] = (p[0], int(p[1]) + 5000, int(p[2]) - 5000)

# segment bounds from ranking table
seg = {}
for r in csv.DictReader(open(os.path.join(BASE, "results", "02_candidate_rebuild", "candidate_ranking_v2.tsv"), encoding="utf-8"), delimiter="\t"):
    if r["gene_symbol"] in regions:
        d = list(map(int, r["distal"].split("-")))
        s = list(map(int, r["shared"].split("-")))
        seg[r["gene_symbol"]] = {"distal": (d[0], d[1]), "shared": (s[0], s[1]), "event_num": r["event_num"], "strand": r["strand"]}

# load depth per sample into per-gene arrays
depth = {g: {s: {} for s in SAMPLES} for g in regions}
for s in SAMPLES:
    with open(os.path.join(BR, f"depth_{s}.tsv"), encoding="utf-8") as f:
        for line in f:
            chrom, pos, v = line.rstrip("\n").split("\t")
            pos = int(pos)
            for g, (c, us, ue) in regions.items():
                if chrom == c and us <= pos < ue:
                    depth[g][s][pos] = int(v)

def stats(g, sample, bounds):
    a, b = bounds
    vals = [depth[g][sample].get(p, 0) for p in range(a, b)]
    n = len(vals)
    mean = sum(vals) / n if n else 0.0
    return mean, sum(1 for v in vals if v >= 1) / n, sum(1 for v in vals if v >= 5) / n, sum(1 for v in vals if v >= 10) / n

rows = []
for g, sg in seg.items():
    for s in SAMPLES:
        dm, d1, d5, d10 = stats(g, s, sg["distal"])
        sm_, s1, s5, s10 = stats(g, s, sg["shared"])
        rows.append({
            "gene_symbol": g, "event_num": sg["event_num"], "sample": s,
            "distal_mean_depth": round(dm, 2), "shared_mean_depth": round(sm_, 2),
            "coverage_ratio_distal_over_shared": round(dm / sm_, 3) if sm_ > 0 else "NA",
            "distal_cov_frac_ge1": round(d1, 3), "distal_cov_frac_ge5": round(d5, 3), "distal_cov_frac_ge10": round(d10, 3),
            "shared_cov_frac_ge1": round(s1, 3), "shared_cov_frac_ge5": round(s5, 3), "shared_cov_frac_ge10": round(s10, 3),
        })

with open(os.path.join(BR, "bam_coverage_review.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(rows)

# compact per-gene view
print(f"{'gene':<9} {'sample':<10} {'distal_mean':>11} {'shared_mean':>11} {'ratio':>7}")
for g in seg:
    for r in rows:
        if r["gene_symbol"] == g:
            print(f"{g:<9} {r['sample']:<10} {r['distal_mean_depth']:>11} {r['shared_mean_depth']:>11} {r['coverage_ratio_distal_over_shared']:>7}")
    print()
