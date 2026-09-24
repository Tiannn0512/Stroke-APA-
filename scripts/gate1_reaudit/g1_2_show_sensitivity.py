#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""打印排主接头敏感性表（工作 2 交付预览）"""
import csv

rows = list(csv.DictReader(
    open("/mnt/d/stroke_apa_reaudit/results/06_gate1_reaudit/coverage_sensitivity_data.tsv",
         encoding="utf-8"), delimiter="\t"))
segs = ["distal_full", "overlap_177974_UTR", "exclusive_3UTR_part", "shared",
        "terminal_exon", "spliced_UTR_intron", "penult_UTR_5p_part"]
print("sample     n_marked  " + "  ".join(s[:14].ljust(16) for s in segs))
print("                            " + "  ".join("full→exclJ".ljust(16) for _ in segs))
for r in rows:
    line = r["sample"].ljust(10) + str(r["marked_reads"]).rjust(8) + "  "
    for s in segs:
        line += (str(r["ratio_" + s]) + "→" + str(r["ratio_exclJ_" + s])).ljust(16)
    print(line)
