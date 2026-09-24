#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G1 复核 · 工作 3：引物特异性筛查（落盘，修复旧交付"文件不存在"问题）。

口径：
  1) 逐引物全长度精确匹配：在 GENCODE vM25 全转录组（264MB，全部 chop 链）中
     搜索每条寡核苷酸的 5'→3' 序列及其反向互补链的完整命中，输出命中数与明细。
  2) 引物对水平：在转录组上扫描 F 引物正向命中与 R 引物反向命中之间的
     候选扩增子（80–3000bp），列出预期 on-target 与任何 off-target 扩增子。
注意：这只覆盖"全长精确匹配"；近似匹配（错配）需 Primer-BLAST（见
primer_specificity_review.md，无法访问 NCBI 时状态记 pending）。
"""
import csv
import sys

FA = "/mnt/d/stroke_apa_reaudit/reference/gencode.vM25.transcripts.fa"
DESIGN = "/mnt/d/stroke_apa_reaudit/results/06_gate1_reaudit/assay_design_corrected.tsv"
OUT = "/mnt/d/stroke_apa_reaudit/results/06_gate1_reaudit/primer_specificity_hits.tsv"
ATP2A2_TIDS = ("ENSMUST00000177974", "ENSMUST00000179939", "ENSMUST00000031423",
               "ENSMUST00000196490", "ENSMUST00000197415")  # Atp2a2 全部 5 转录本


def rc(s):
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def load_transcripts():
    tx = {}
    name = None
    buf = []
    with open(FA) as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    tx[name] = "".join(buf)
                name = line[1:].split("|")[0]
                buf = []
            else:
                buf.append(line.strip())
    if name:
        tx[name] = "".join(buf)
    return tx


def find_all(hay, needle):
    out = []
    i = hay.find(needle)
    while i != -1:
        out.append(i)
        i = hay.find(needle, i + 1)
    return out


def main():
    tx = load_transcripts()
    print(f"transcripts loaded: {len(tx)}")
    design = list(csv.DictReader(open(DESIGN), delimiter="\t"))

    # ---- 逐引物精确命中 ----
    hits = {}
    for r in design:
        oligo = r["seq5to3"]
        det = []
        n = 0
        for tid, seq in tx.items():
            for p in find_all(seq, oligo):
                n += 1
                det.append((tid, p, "+"))
            for p in find_all(seq, rc(oligo)):
                n += 1
                det.append((tid, p, "-"))
        hits[oligo] = (n, det)

    # ---- 引物对扫描 ----
    assays = {}
    for r in design:
        if r["role"] in ("F", "R") and "GSP" not in r["assay"]:
            assays.setdefault(r["assay"], {})[r["role"]] = r

    pair_rows = []
    for assay, pr in assays.items():
        if "F" not in pr or "R" not in pr:
            continue
        F, R = pr["F"]["seq5to3"], pr["R"]["seq5to3"]
        found = []
        for tid, seq in tx.items():
            # 方向1：F 正向命中，R 的反向互补在其下游（Atp2a2 为负链基因，
            # mRNA 记录上表现为 rc(F) 正向——两方向都要扫）
            for pF in find_all(seq, F):
                pR = seq.find(rc(R), pF + len(F), min(len(seq), pF + 3000 + len(F)))
                if pR != -1:
                    size = pR + len(R) - pF
                    on = tid.startswith(ATP2A2_TIDS)
                    found.append((tid, pF, size, "on-target" if on else "OFF-TARGET"))
            # 方向2：rc(F) 正向命中，R 正向在其上游
            for pF in find_all(seq, rc(F)):
                lo = max(0, pF - 3000)
                seg = seq[lo:pF]
                idx = seg.rfind(R)
                if idx != -1:
                    pR = lo + idx
                    size = pF + len(F) - pR
                    on = tid.startswith(ATP2A2_TIDS)
                    found.append((tid, pR, size, "on-target" if on else "OFF-TARGET"))
        pair_rows.append((assay, F, R, found))

    with open(OUT, "w") as f:
        f.write("assay\trole\toligo_seq5to3\ttranscriptome_exact_hits\thit_details(top6)\n")
        for r in design:
            n, det = hits[r["seq5to3"]]
            top = "; ".join(f"{t}:{p}({st})" for t, p, st in det[:6])
            f.write(f"{r['assay']}\t{r['role']}\t{r['seq5to3']}\t{n}\t{top}\n")
        f.write("\n# pair-level amplicon scan (transcriptome, 80-3000bp)\n")
        f.write("assay\tF_seq\tR_seq\tamplicons_found\tdetails\n")
        for assay, F, R, found in pair_rows:
            d = "; ".join(f"{t}:{p}:{s}bp:{tag}" for t, p, s, tag in found[:8])
            f.write(f"{assay}\t{F}\t{R}\t{len(found)}\t{d}\n")

    print(f"written {OUT}")
    for r in design:
        n, det = hits[r["seq5to3"]]
        print(f"  {r['assay']:<34}{r['role']:<3}hits={n:<4}"
              f"{('on-target' if n and det[0][0].startswith(ATP2A2_TIDS) else ('none' if not n else 'check'))}")
    for assay, F, R, found in pair_rows:
        on = sum(1 for t, p, s, tag in found if tag == "on-target")
        off = len(found) - on
        print(f"  PAIR {assay:<30}on-target amplicons={on}  OFF-TARGET={off}")


if __name__ == "__main__":
    main()


# ---------- 可导入的扫描接口（供设计脚本挑选零脱靶候选） ----------
ALLOWED = ATP2A2_TIDS


def oligo_all_hits(oligo, tx):
    """返回 [(tid, pos, orient)]：oligo 及其 rc 在转录组全部全长精确命中。"""
    det = []
    for tid, seq in tx.items():
        for p in find_all(seq, oligo):
            det.append((tid, p, "+"))
        for p in find_all(seq, rc(oligo)):
            det.append((tid, p, "-"))
    return det


def pair_amplicons(F, R, tx):
    """扩增子扫描（保守版）：dsDNA 下四种镜像排布均为真产物，故任何 F 变体位点
    （F 或 rc(F)）与任何 R 变体位点（R 或 rc(R)）距离 ≤3000bp 即计为一个候选
    扩增子。宁误报不漏报：此处 OFF>0 的候选一律淘汰。"""
    found = []
    for tid, seq in tx.items():
        f_sites = [(p, len(F)) for p in find_all(seq, F)] +                   [(p, len(F)) for p in find_all(seq, rc(F))]
        r_sites = [(p, len(R)) for p in find_all(seq, R)] +                   [(p, len(R)) for p in find_all(seq, rc(R))]
        for pf, lf in f_sites:
            for pr, lr in r_sites:
                if abs(pr - pf) <= 3000:
                    lo = min(pf, pr)
                    hi = max(pf + lf, pr + lr)
                    size = hi - lo
                    tag = "on" if tid.startswith(ALLOWED) else "OFF"
                    found.append((tid, lo, size, tag))
    return found
