#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""打印接头×转录本注释摘要（工作 1 交付预览）"""
import csv

rows = list(csv.DictReader(
    open("/mnt/d/stroke_apa_reaudit/results/06_gate1_reaudit/junction_transcript_annotation.tsv",
         encoding="utf-8"), delimiter="\t"))
rows.sort(key=lambda r: -int(r["total"]))
hdr = f"{'junction':<26}{'total':>6}  {'vs_179939(长UTR)':<26}" \
      f"{'vs_177974(短UTR)':<24}{'vs_31423(中UTR)':<20}utr_rel_179939"
print(hdr)
for r in rows[:10]:
    j = r["junction_start0"] + "-" + r["junction_end0"]
    print(f"{j:<26}{r['total']:>6}  {r['vs_ENSMUST00000179939.7']:<26}"
          f"{r['vs_ENSMUST00000177974.7']:<24}{r['vs_ENSMUST00000031423.9']:<20}"
          f"{r['utr_rel_ENSMUST00000179939.7']}")
