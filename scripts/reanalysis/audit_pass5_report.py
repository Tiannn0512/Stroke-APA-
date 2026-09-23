#!/usr/bin/env python3
"""AUDIT PASS 5 — verify the final report's numbers against audited values.
Every load-bearing figure quoted in REANALYSIS_FINAL_REPORT.md must appear and
match the pass-1 recomputed values. Exit 1 on any FAIL."""
import os, sys

BASE = r"D:\stroke_apa_reanalysis"
rep = open(os.path.join(BASE, "REANALYSIS_FINAL_REPORT.md"), encoding="utf-8").read()
FAILS = []

def check(label, needle):
    ok = needle in rep
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if not ok:
        FAILS.append(label)

# P1
check("9693 events", "9,693 事件全部拿到链方向")
check("+4913/-4780", "+4,913/−4,780")
check("15/15 unit tests", "15/15")
check("1bp: proximal_pas_raw=proximal_pas_bed", "proximal_pas_raw = proximal_pas_bed")
# P2 QKI
check("Atp2a2 peak FDR 1.04e-05", "FDR 1.04e-05")
check("42/9,693 distal QKI", "42/9,693")
check("old bins FDR", "0.0055/0.0709")
# P2 motif
check("1541 pairs", "1,541")
check("824 unique", "824")
check("null med 21.5nt/0.03pp", "21.5nt、GC 差 0.03pp")
check("61488 hits", "61,488")
check("792/824", "792/824")
check("min p 0.054", "p=0.054")
check("q>=0.71", "q≥0.71")
# P2 conservation
check("9667/9693", "9,667/9,693")
check("99.7%", "99.7%")
check("median 0.299", "0.299")
check("v1 65% artifact", "bedGraph 中间管线伪影")
# class
check("class B=8 C=783 D=1 none=185", "B=8 / C=783 / D=1 / none=185")
check("346 flips", "346")
# P3 coverage
check("Atp2a2 0.302/0.084/0.085", "0.302 | 0.084 | 0.085")
check("-71.7%", "-71.7%")
check("Agpat3 -82.3%", "-82.3%")
check("shared rises 218->530", "218→530")
check("arich 10/10 low", "10/10 low")
check("8 advance 2 hold 0 drop", "8 advance（high 2）/ 2 hold / 0 drop")
# P4
check("6762 genes", "6,762")
check("median -3.42", "−3.42")
check("Atp2a2 54.0 pct", "54.0% 百分位（中位）")
check("Plec 95.8", "95.8%")
check("GSE330741 651 genes absent", "均不存在")
# Gate1
check("PAS 13bp", "13bp")
check("salmon traj", "0.719/0.541 → d3 0.130/0.112 → d7 0.039/0.039 → d60 0.275/0.127")
check("Aplp1 single-tx note", "1 条转录本")
check("Ndrg2 PAS 8bp", "PAS 8bp")
# extras
check("PAS atlas 301,006", "301,006")
import os as _os
_z = r"D:\stroke_apa\ATP2A2_IGV_CHECK_PACKAGE_20260920.zip"
_zok = _os.path.exists(_z) and _os.path.getsize(_z) == 29569120
print(f"[{'PASS' if _zok else 'FAIL'}] IGV zip on disk with exact size 29,569,120 (direct file check)")
if not _zok:
    FAILS.append("zip on disk")
check("Atp2a2 PDUI traj unchanged", "0.19/0.20 → d3 0.08/0.06 → d7 0.09/0.09")
check("five-pass audit stated", "五遍")

print()
print("AUDIT PASS 5:", "ALL PASS" if not FAILS else f"{len(FAILS)} FAILURES -> {FAILS}")
sys.exit(1 if FAILS else 0)
