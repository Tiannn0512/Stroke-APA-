#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G1 复核 · 工作 3 v5：引物设计（最终版）——mRNA 空间拼接模板 + 候选内零脱靶筛选。

结构事实（GENCODE vM25，Atp2a2 负链，BED 0-based half-open）：
  179939（长 UTR，事件载体）3'UTR 跨内含子 122454304-122455892
  177974（短 UTR）UTR 122453513-122454293；31423（中 UTR）UTR 122456339-122457153
  Gm30970（ENSMUST00000227710.1）= Atp2a2 加工型假基因（chr14，363bp 记录），
  精确匹配 Atp2a2 mRNA ~3156-3330 区域（覆盖部分共有外显子）。

选址逻辑（全部经内联转录组扫描，off-target==0 优先）：
  A. RACE_outer_GSP：拼接模板 [122458447 起的共有外显子]，候选须位于 122458221 上游
  B. RACE_nested_GSP：拼接模板 [122457511-457629 + 122458087-458221]，
     候选须位于 122457629 上游（=122458087-458221 外显子，假基因覆盖窗外）
  C. qPCR_common：拼接模板 [(122489237,122489376)+(122491680,122491785)]（远离假基因窗），
     跨内含子产物，on-target 须覆盖 3 转录本
  D. distal_terminal：终端外显子 UTR 122453513-122454293（假基因窗外，经验证干净）
  E. long_junction：cDNA 跨 spliced-UTR 内含子，仅 179939 可扩增
若 A/B 在清洁区无候选，允许唯一例外 = 额外命中仅为 Gm30970 假基因，并在 product
字段显式登记残余风险。
"""
import csv
import sys

import primer3

sys.path.insert(0, "/mnt/d/stroke_apa_reaudit/scripts")
from g1_3_seq_utils import verify_window
from g1_3_specificity import load_transcripts, oligo_all_hits, pair_amplicons
from g1_1_junction_annotation import parse_gtf

GENOME = "/mnt/d/stroke_apa_reaudit/reference/chr5.fa"
GTF = "/mnt/d/stroke_apa_reaudit/reference/gencode.vM25.annotation.gtf.gz"
OUT_TSV = "/mnt/d/stroke_apa_reaudit/results/06_gate1_reaudit/assay_design_corrected.tsv"
TARGETS = ["ENSMUST00000177974.7", "ENSMUST00000179939.7", "ENSMUST00000031423.9"]
ATP_ALL = ("ENSMUST00000031423", "ENSMUST00000177974", "ENSMUST00000179939",
          "ENSMUST00000196490", "ENSMUST00000197415")   # Atp2a2 全部 5 转录本
T3 = {"ENSMUST00000177974", "ENSMUST00000179939", "ENSMUST00000031423"}
PG = "ENSMUST00000227710"

PAIR_SETTINGS = {
    "PRIMER_TASK": "generic",
    "PRIMER_PICK_LEFT_PRIMER": 1,
    "PRIMER_PICK_RIGHT_PRIMER": 1,
    "PRIMER_PICK_INTERNAL_OLIGO": 0,
    "PRIMER_PRODUCT_SIZE_RANGE": [[80, 200]],
    "PRIMER_OPT_SIZE": 22, "PRIMER_MIN_SIZE": 20, "PRIMER_MAX_SIZE": 30,
    "PRIMER_OPT_TM": 60.0, "PRIMER_MIN_TM": 57.0, "PRIMER_MAX_TM": 63.0,
    "PRIMER_MIN_GC": 40.0, "PRIMER_MAX_GC": 60.0,
    "PRIMER_MAX_SELF_ANY_TH": 42.0, "PRIMER_MAX_SELF_END_TH": 32.0,
    "PRIMER_MAX_HAIRPIN_TH": 40.0,
    "PRIMER_NUM_RETURN": 20,
}


def common_exons():
    exons, _, _ = parse_gtf(GTF)
    sets = [set(map(tuple, exons[t])) for t in TARGETS]
    return sorted(sets[0] & sets[1] & sets[2])


def spliced_template(common, hi, lo):
    """[lo, hi] 范围内共有外显子按负链转录方向拼接；返回 (模板, gmap)。"""
    sel = [e for e in common if lo <= e[0] and e[1] <= hi]
    seq, gmap = [], []
    for ex in sorted(sel, reverse=True):
        s = verify_window(GENOME, "chr5", ex[0], ex[1], "-")
        seq.append(s)
        gmap += list(range(ex[1] - 1, ex[0] - 1, -1))
    return "".join(seq), gmap


def pick_gsp(name, tmpl, gmap, min_g, tm_settings, tx):
    """设计 RIGHT 候选并筛选：优退次序 = 清洁 > 仅 Gm30970 额外命中。"""
    res = primer3.bindings.design_primers(
        {"SEQUENCE_ID": name, "SEQUENCE_TEMPLATE": tmpl},
        {"PRIMER_TASK": "generic", "PRIMER_PICK_LEFT_PRIMER": 1,
         "PRIMER_PICK_RIGHT_PRIMER": 1, "PRIMER_PICK_INTERNAL_OLIGO": 0,
         "PRIMER_PRODUCT_SIZE_RANGE": [[70, len(tmpl)]], "PRIMER_NUM_RETURN": 20,
         "PRIMER_MIN_SIZE": 22, "PRIMER_OPT_SIZE": 26, "PRIMER_MAX_SIZE": 30,
         "PRIMER_MIN_GC": 40.0, "PRIMER_MAX_GC": 60.0,
         "PRIMER_MAX_SELF_ANY_TH": 42.0, "PRIMER_MAX_SELF_END_TH": 32.0,
         "PRIMER_MAX_HAIRPIN_TH": 40.0, **tm_settings})
    n = res.get("PRIMER_PAIR_NUM_RETURNED", 0)
    fallback = None
    for i in range(n):
        pos3, ln = res[f"PRIMER_RIGHT_{i}"]
        gl, gh = gmap[pos3 - ln + 1], gmap[pos3]
        if gl < min_g:
            continue
        seq = res[f"PRIMER_RIGHT_{i}_SEQUENCE"]
        hits = oligo_all_hits(seq, tx)
        bad = [h for h in hits if not h[0].startswith(tuple(t[:15] for t in TARGETS))]
        if not bad:
            return i, res, n, "清洁（精确命中全为 Atp2a2 转录本）"
        if all(h[0].startswith(PG) for h in bad) and fallback is None:
            fallback = (i, "含 Gm30970 假基因精确匹配——残余风险：假基因全长 363bp，"
                           "RACE 产物将短于 363bp 且测序可辨；建议 BLAST 复核")
    if fallback:
        return fallback[0], res, n, fallback[1]
    return None, res, n, ""


def main():
    tx = load_transcripts()
    print(f"transcripts loaded: {len(tx)}")
    common = common_exons()
    print("共有外显子（与本设计相关）:",
          [e for e in common if 122457000 <= e[0] <= 122500000])
    rows = []

    # ---------- A. RACE_outer_GSP ----------
    otmpl, ogmap = spliced_template(common, 122475000, 122458447)
    print(f"outer 模板长: {len(otmpl)}")
    i, res, n, note = pick_gsp("RACE_outer_GSP", otmpl, ogmap, 122458221,
                               {"PRIMER_OPT_TM": 62.0, "PRIMER_MIN_TM": 58.0,
                                "PRIMER_MAX_TM": 66.0}, tx)
    if i is None:
        print("!! RACE_outer_GSP: 无候选")
    else:
        pos3, ln = res[f"PRIMER_RIGHT_{i}"]
        gl, gh = ogmap[pos3 - ln + 1], ogmap[pos3]
        rows.append(dict(assay="RACE_outer_GSP", role="GSP(antisense)",
                         oligo=f"reverse-primer#cand{i}",
                         seq5to3=res[f"PRIMER_RIGHT_{i}_SEQUENCE"], chrom="chr5",
                         g_start0=gl, g_end0=gh + 1, strand="-",
                         tm=round(res[f"PRIMER_RIGHT_{i}_TM"], 2),
                         gc=round(res[f"PRIMER_RIGHT_{i}_GC_PERCENT"], 2),
                         product=f"3'RACE outer GSP；{note}（{i}/{n} 候选中选出）",
                         dimer_self_any=round(res[f"PRIMER_RIGHT_{i}_SELF_ANY_TH"], 1),
                         dimer_self_end=round(res[f"PRIMER_RIGHT_{i}_SELF_END_TH"], 1),
                         hairpin=round(res[f"PRIMER_RIGHT_{i}_HAIRPIN_TH"], 1)))

    # ---------- B. RACE_nested_GSP ----------
    ntmpl, ngmap = spliced_template(common, 122458221, 122457511)
    print(f"nested 模板长: {len(ntmpl)}")
    i, resB, n, note = pick_gsp("RACE_nested_GSP", ntmpl, ngmap, 122457629,
                                {"PRIMER_OPT_TM": 62.0, "PRIMER_MIN_TM": 58.0,
                                 "PRIMER_MAX_TM": 66.0}, tx)
    if i is None:
        print("!! RACE_nested_GSP: 无候选")
    else:
        pos3, ln = resB[f"PRIMER_RIGHT_{i}"]
        gl, gh = ngmap[pos3 - ln + 1], ngmap[pos3]
        rows.append(dict(assay="RACE_nested_GSP", role="GSP(antisense)",
                         oligo=f"reverse-primer#cand{i}",
                         seq5to3=resB[f"PRIMER_RIGHT_{i}_SEQUENCE"], chrom="chr5",
                         g_start0=gl, g_end0=gh + 1, strand="-",
                         tm=round(resB[f"PRIMER_RIGHT_{i}_TM"], 2),
                         gc=round(resB[f"PRIMER_RIGHT_{i}_GC_PERCENT"], 2),
                         product=f"3'RACE nested GSP；{note}（{i}/{n} 候选中选出）",
                         dimer_self_any=round(resB[f"PRIMER_RIGHT_{i}_SELF_ANY_TH"], 1),
                         dimer_self_end=round(resB[f"PRIMER_RIGHT_{i}_SELF_END_TH"], 1),
                         hairpin=round(resB[f"PRIMER_RIGHT_{i}_HAIRPIN_TH"], 1)))

    # ---------- C. qPCR_common（跨内含子，远离假基因窗） ----------
    ctmpl, cgmap = spliced_template(common, 122491785, 122489237)
    print(f"common 模板长: {len(ctmpl)}")
    resC = primer3.bindings.design_primers(
        {"SEQUENCE_ID": "qPCR_common", "SEQUENCE_TEMPLATE": ctmpl},
        dict(PAIR_SETTINGS, PRIMER_PRODUCT_SIZE_RANGE=[[120, len(ctmpl)]]))
    n = resC.get("PRIMER_PAIR_NUM_RETURNED", 0)
    chosen, note = None, ""
    T3 = {"ENSMUST00000177974", "ENSMUST00000179939", "ENSMUST00000031423"}
    for i in range(n):
        F = resC[f"PRIMER_LEFT_{i}_SEQUENCE"]
        R = resC[f"PRIMER_RIGHT_{i}_SEQUENCE"]
        amp = pair_amplicons(F, R, tx)
        on_tids = {a[0] for a in amp if a[3] == "on"}
        off = sum(1 for a in amp if a[3] == "OFF")
        if off == 0 and T3 <= {t.split(".")[0] for t in on_tids} and all(
                t.startswith(ATP_ALL) for t in on_tids):
            chosen, note = i, f"（{i}/{n} 候选中选出；3 目标转录本均扩增且命中限于"                               f"Atp2a2 基因，OFF-TARGET=0）"
            break
    if chosen is None:
        print("!! qPCR_common: 无零脱靶候选")
    else:
        i = chosen
        L = resC[f"PRIMER_LEFT_{i}"]
        R = resC[f"PRIMER_RIGHT_{i}"]
        fl, fh = cgmap[L[0]], cgmap[L[0] + L[1] - 1]
        rl, rh = cgmap[R[0] - R[1] + 1], cgmap[R[0]]
        rows.append(dict(assay="qPCR_common(3 transcripts)", role="F",
                         oligo=f"forward#cand{i}", seq5to3=resC[f"PRIMER_LEFT_{i}_SEQUENCE"],
                         chrom="chr5", g_start0=min(fl, fh), g_end0=max(fl, fh) + 1,
                         strand="-", tm=round(resC[f"PRIMER_LEFT_{i}_TM"], 2),
                         gc=round(resC[f"PRIMER_LEFT_{i}_GC_PERCENT"], 2),
                         product=f"跨内含子 cDNA 产物 {resC[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']}bp {note}",
                         dimer_self_any=round(resC[f"PRIMER_LEFT_{i}_SELF_ANY_TH"], 1),
                         dimer_self_end=round(resC[f"PRIMER_LEFT_{i}_SELF_END_TH"], 1),
                         hairpin=round(resC[f"PRIMER_LEFT_{i}_HAIRPIN_TH"], 1)))
        rows.append(dict(assay="qPCR_common(3 transcripts)", role="R",
                         oligo=f"reverse#cand{i}", seq5to3=resC[f"PRIMER_RIGHT_{i}_SEQUENCE"],
                         chrom="chr5", g_start0=min(rl, rh), g_end0=max(rl, rh) + 1,
                         strand="+", tm=round(resC[f"PRIMER_RIGHT_{i}_TM"], 2),
                         gc=round(resC[f"PRIMER_RIGHT_{i}_GC_PERCENT"], 2),
                         product=f"跨内含子 cDNA 产物 {resC[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']}bp",
                         dimer_self_any=round(resC[f"PRIMER_RIGHT_{i}_SELF_ANY_TH"], 1),
                         dimer_self_end=round(resC[f"PRIMER_RIGHT_{i}_SELF_END_TH"], 1),
                         hairpin=round(resC[f"PRIMER_RIGHT_{i}_HAIRPIN_TH"], 1)))

    # ---------- D. distal_terminal ----------
    tmplD = verify_window(GENOME, "chr5", 122453513, 122454293, "+")
    resD = primer3.bindings.design_primers(
        {"SEQUENCE_ID": "qPCR_distal_terminal", "SEQUENCE_TEMPLATE": tmplD},
        dict(PAIR_SETTINGS, PRIMER_PRODUCT_SIZE_RANGE=[[80, 180]]))
    n = resD.get("PRIMER_PAIR_NUM_RETURNED", 0)
    chosen, note = None, ""
    for i in range(n):
        F = resD[f"PRIMER_LEFT_{i}_SEQUENCE"]
        R = resD[f"PRIMER_RIGHT_{i}_SEQUENCE"]
        amp = pair_amplicons(F, R, tx)
        on_tids = {a[0] for a in amp if a[3] == "on"}
        off = sum(1 for a in amp if a[3] == "OFF")
        if off == 0 and {"ENSMUST00000177974", "ENSMUST00000179939"} <= {
                t.split(".")[0] for t in on_tids} and all(
                t.startswith(ATP_ALL) for t in on_tids):
            chosen, note = i, f"（{i}/{n} 候选中选出；177974+179939 均扩增且全部命中"                               f"限于 Atp2a2 基因，OFF-TARGET=0）"
            break
    if chosen is None:
        print("!! distal_terminal: 无零脱靶候选")
    else:
        i = chosen
        L = resD[f"PRIMER_LEFT_{i}"]
        R = resD[f"PRIMER_RIGHT_{i}"]
        rows.append(dict(assay="distal_terminal(177974+179939)", role="F",
                         oligo=f"forward#cand{i}", seq5to3=resD[f"PRIMER_LEFT_{i}_SEQUENCE"],
                         chrom="chr5",
                         g_start0=122453513 + L[0], g_end0=122453513 + L[0] + L[1],
                         strand="+", tm=round(resD[f"PRIMER_LEFT_{i}_TM"], 2),
                         gc=round(resD[f"PRIMER_LEFT_{i}_GC_PERCENT"], 2),
                         product=f"产物 {resD[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']}bp {note}",
                         dimer_self_any=round(resD[f"PRIMER_LEFT_{i}_SELF_ANY_TH"], 1),
                         dimer_self_end=round(resD[f"PRIMER_LEFT_{i}_SELF_END_TH"], 1),
                         hairpin=round(resD[f"PRIMER_LEFT_{i}_HAIRPIN_TH"], 1)))
        rows.append(dict(assay="distal_terminal(177974+179939)", role="R",
                         oligo=f"reverse#cand{i}", seq5to3=resD[f"PRIMER_RIGHT_{i}_SEQUENCE"],
                         chrom="chr5",
                         g_start0=122453513 + R[0] - R[1] + 1,
                         g_end0=122453513 + R[0] + 1,
                         strand="-", tm=round(resD[f"PRIMER_RIGHT_{i}_TM"], 2),
                         gc=round(resD[f"PRIMER_RIGHT_{i}_GC_PERCENT"], 2),
                         product=f"产物 {resD[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']}bp",
                         dimer_self_any=round(resD[f"PRIMER_RIGHT_{i}_SELF_ANY_TH"], 1),
                         dimer_self_end=round(resD[f"PRIMER_RIGHT_{i}_SELF_END_TH"], 1),
                         hairpin=round(resD[f"PRIMER_RIGHT_{i}_HAIRPIN_TH"], 1)))

    # ---------- E. long_junction ----------
    seg1 = verify_window(GENOME, "chr5", 122453512, 122454304, "-")
    seg2 = verify_window(GENOME, "chr5", 122455892, 122456500, "-")
    cdna = seg1 + seg2
    gmap = (list(range(122454303, 122453511, -1))
            + list(range(122456499, 122455891, -1)))
    resE = primer3.bindings.design_primers(
        {"SEQUENCE_ID": "junction_179939", "SEQUENCE_TEMPLATE": cdna},
        dict(PAIR_SETTINGS, PRIMER_PRODUCT_SIZE_RANGE=[[900, 1100]]))
    n = resE.get("PRIMER_PAIR_NUM_RETURNED", 0)
    chosen, note = None, ""
    for i in range(n):
        F = resE[f"PRIMER_LEFT_{i}_SEQUENCE"]
        R = resE[f"PRIMER_RIGHT_{i}_SEQUENCE"]
        amp = pair_amplicons(F, R, tx)
        on_tids = {a[0] for a in amp if a[3] == "on"}
        off = sum(1 for a in amp if a[3] == "OFF")
        if (off == 0 and on_tids == {"ENSMUST00000179939.7"}):
            chosen, note = i, f"（{i}/{n} 候选中选出；仅 179939 扩增，"                               f"177974/31423 不扩增，OFF=0）"
            break
    if chosen is None:
        print("!! long_junction: 无零脱靶候选")
    else:
        i = chosen
        Lp, Ll = resE[f"PRIMER_LEFT_{i}"]
        Rp, Rl = resE[f"PRIMER_RIGHT_{i}"]
        for role, p, ln, seq, tm, gc, sany, send, hp in (
                ("F", Lp, Ll, resE[f"PRIMER_LEFT_{i}_SEQUENCE"], resE[f"PRIMER_LEFT_{i}_TM"],
                 resE[f"PRIMER_LEFT_{i}_GC_PERCENT"], resE[f"PRIMER_LEFT_{i}_SELF_ANY_TH"],
                 resE[f"PRIMER_LEFT_{i}_SELF_END_TH"], resE[f"PRIMER_LEFT_{i}_HAIRPIN_TH"]),
                ("R", Rp, Rl, resE[f"PRIMER_RIGHT_{i}_SEQUENCE"], resE[f"PRIMER_RIGHT_{i}_TM"],
                 resE[f"PRIMER_RIGHT_{i}_GC_PERCENT"], resE[f"PRIMER_RIGHT_{i}_SELF_ANY_TH"],
                 resE[f"PRIMER_RIGHT_{i}_SELF_END_TH"], resE[f"PRIMER_RIGHT_{i}_HAIRPIN_TH"])):
            if role == "F":
                g0, g1 = gmap[p], gmap[p + ln - 1]
            else:
                g0, g1 = gmap[p - ln + 1], gmap[p]
            rows.append(dict(
                assay="long_junction(179939-specific)", role=role,
                oligo=f"{'forward' if role == 'F' else 'reverse'}#cand{i}", seq5to3=seq,
                chrom="chr5", g_start0=min(g0, g1), g_end0=max(g0, g1) + 1, strand="-",
                tm=round(tm, 2), gc=round(gc, 2),
                product=f"cDNA {resE[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']}bp across intron "
                        f"122454304-122455892 {note}",
                dimer_self_any=round(sany, 1), dimer_self_end=round(send, 1),
                hairpin=round(hp, 1)))

    cols = ["assay", "role", "oligo", "seq5to3", "chrom", "g_start0", "g_end0",
            "strand", "tm", "gc", "product", "dimer_self_any", "dimer_self_end", "hairpin"]
    with open(OUT_TSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"written {OUT_TSV}: {len(rows)} oligos")
    for r in rows:
        print(f"  {r['assay']:<32}{r['role']:<15}{r['seq5to3']:<28}"
              f"chr5:{r['g_start0']}-{r['g_end0']} Tm={r['tm']} GC={r['gc']}")


if __name__ == "__main__":
    main()
