#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G1 复核 · 工作 3：在修正的序列与修正的转录本结构上重做 Atp2a2 引物设计（v2）。

primer3 坐标约定（probe_primer3.py 实证）：
  LEFT_i  = (start, len)，占模板 [start, start+len)
  RIGHT_i = (pos, len)，pos 为引物 3' 端（模板最右）碱基索引，占 [pos-len+1, pos+1)

结构事实（GENCODE vM25，负链，BED 0-based half-open）：
  179939（长 UTR）UTR 跨内含子 122454304-122455892
  177974（短 UTR）UTR 122453513-122454293；31423（中 UTR）UTR 122456339-122457153
  三转录本共有外显子：122457302-122457426、122457511-122457629
"""
import csv
import sys

import primer3

sys.path.insert(0, "/mnt/d/stroke_apa_reaudit/scripts")
from g1_3_seq_utils import verify_window

GENOME = "/mnt/d/stroke_apa_reaudit/reference/chr5.fa"
OUT_TSV = "/mnt/d/stroke_apa_reaudit/results/06_gate1_reaudit/assay_design_corrected.tsv"

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
}


def map_left(start0, end0, strand, pos, ln):
    if strand == "+":
        return start0 + pos, start0 + pos + ln
    return end0 - (pos + ln), end0 - pos


def map_right(start0, end0, strand, pos3, ln):
    if strand == "+":
        return start0 + pos3 - ln + 1, start0 + pos3 + 1
    return end0 - pos3 - 1, end0 - pos3 + ln - 1


def oligo_rows(assay, res, start0, end0, strand):
    i = 0
    L = res[f"PRIMER_LEFT_{i}"]
    R = res[f"PRIMER_RIGHT_{i}"]
    gl = map_left(start0, end0, strand, L[0], L[1])
    gr = map_right(start0, end0, strand, R[0], R[1])
    return [
        dict(assay=assay, role="F", oligo="forward",
             seq5to3=res[f"PRIMER_LEFT_{i}_SEQUENCE"], chrom="chr5",
             g_start0=gl[0], g_end0=gl[1], strand=strand,
             tm=round(res[f"PRIMER_LEFT_{i}_TM"], 2),
             gc=round(res[f"PRIMER_LEFT_{i}_GC_PERCENT"], 2),
             product="amplicon %dbp penalty=%.2f" % (
                 res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"], res[f"PRIMER_PAIR_{i}_PENALTY"]),
             dimer_self_any=round(res[f"PRIMER_LEFT_{i}_SELF_ANY_TH"], 1),
             dimer_self_end=round(res[f"PRIMER_LEFT_{i}_SELF_END_TH"], 1),
             hairpin=round(res[f"PRIMER_LEFT_{i}_HAIRPIN_TH"], 1)),
        dict(assay=assay, role="R", oligo="reverse",
             seq5to3=res[f"PRIMER_RIGHT_{i}_SEQUENCE"], chrom="chr5",
             g_start0=gr[0], g_end0=gr[1], strand=strand,
             tm=round(res[f"PRIMER_RIGHT_{i}_TM"], 2),
             gc=round(res[f"PRIMER_RIGHT_{i}_GC_PERCENT"], 2),
             product="amplicon %dbp" % res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"],
             dimer_self_any=round(res[f"PRIMER_RIGHT_{i}_SELF_ANY_TH"], 1),
             dimer_self_end=round(res[f"PRIMER_RIGHT_{i}_SELF_END_TH"], 1),
             hairpin=round(res[f"PRIMER_RIGHT_{i}_HAIRPIN_TH"], 1)),
    ]


def main():
    rows = []

    # A/B: 3'RACE GSP（模板=转录方向；GSP=RIGHT 引物）
    for name, s0, e0 in (("RACE_outer_GSP", 122457511, 122457629),
                         ("RACE_nested_GSP", 122457302, 122457426)):
        tmpl = verify_window(GENOME, "chr5", s0, e0, "-")
        res = primer3.bindings.design_primers(
            {"SEQUENCE_ID": name, "SEQUENCE_TEMPLATE": tmpl},
            {"PRIMER_TASK": "generic", "PRIMER_PICK_LEFT_PRIMER": 0,
             "PRIMER_PICK_RIGHT_PRIMER": 1, "PRIMER_PICK_INTERNAL_OLIGO": 0,
             "PRIMER_PRODUCT_SIZE_RANGE": [[20, e0 - s0]], "PRIMER_NUM_RETURN": 1,
             "PRIMER_MIN_SIZE": 22, "PRIMER_OPT_SIZE": 26, "PRIMER_MAX_SIZE": 30,
             "PRIMER_OPT_TM": 62.0, "PRIMER_MIN_TM": 58.0, "PRIMER_MAX_TM": 66.0,
             "PRIMER_MIN_GC": 40.0, "PRIMER_MAX_GC": 60.0,
             "PRIMER_MAX_SELF_ANY_TH": 42.0, "PRIMER_MAX_SELF_END_TH": 32.0,
             "PRIMER_MAX_HAIRPIN_TH": 40.0})
        if not res.get("PRIMER_PAIR_NUM_RETURNED"):
            print("!! %s: no primer returned" % name)
            continue
        pos3, ln = res["PRIMER_RIGHT_0"]
        gl, gh = map_right(s0, e0, "-", pos3, ln)
        rows.append(dict(assay=name, role="GSP(antisense)", oligo="reverse-primer",
                         seq5to3=res["PRIMER_RIGHT_0_SEQUENCE"], chrom="chr5",
                         g_start0=gl, g_end0=gh, strand="-",
                         tm=round(res["PRIMER_RIGHT_0_TM"], 2),
                         gc=round(res["PRIMER_RIGHT_0_GC_PERCENT"], 2),
                         product="3'RACE GSP bind chr5:%d-%d (%d cands)" % (gl, gh, res["PRIMER_PAIR_NUM_RETURNED"]),
                         dimer_self_any=round(res["PRIMER_RIGHT_0_SELF_ANY_TH"], 1),
                         dimer_self_end=round(res["PRIMER_RIGHT_0_SELF_END_TH"], 1),
                         hairpin=round(res["PRIMER_RIGHT_0_HAIRPIN_TH"], 1)))

    # C: common 定量 assay
    tmplC = verify_window(GENOME, "chr5", 122457511, 122457629, "+")
    resC = primer3.bindings.design_primers(
        {"SEQUENCE_ID": "qPCR_common", "SEQUENCE_TEMPLATE": tmplC},
        dict(PAIR_SETTINGS, PRIMER_PRODUCT_SIZE_RANGE=[[70, 118]]))
    rows += oligo_rows("qPCR_common(3 transcripts)", resC, 122457511, 122457629, "+")

    # D: distal-terminal assay
    tmplD = verify_window(GENOME, "chr5", 122453513, 122454293, "+")
    resD = primer3.bindings.design_primers(
        {"SEQUENCE_ID": "qPCR_distal_terminal", "SEQUENCE_TEMPLATE": tmplD},
        dict(PAIR_SETTINGS, PRIMER_PRODUCT_SIZE_RANGE=[[80, 180]]))
    rows += oligo_rows("distal_terminal(177974+179939)", resD, 122453513, 122454293, "+")

    # E: long-isoform junction assay（cDNA 强制跨内含子）
    seg1 = verify_window(GENOME, "chr5", 122453512, 122454304, "-")
    seg2 = verify_window(GENOME, "chr5", 122455892, 122456500, "-")
    cdna = seg1 + seg2
    gmap = (list(range(122454303, 122453511, -1))
            + list(range(122456499, 122455891, -1)))
    assert len(gmap) == len(cdna), (len(gmap), len(cdna))
    resE = primer3.bindings.design_primers(
        {"SEQUENCE_ID": "junction_179939", "SEQUENCE_TEMPLATE": cdna},
        dict(PAIR_SETTINGS, PRIMER_PRODUCT_SIZE_RANGE=[[900, 1100]]))
    if not resE.get("PRIMER_PAIR_NUM_RETURNED"):
        print("!! junction assay: no pair returned")
    else:
        Lp, Ll = resE["PRIMER_LEFT_0"]
        Rp, Rl = resE["PRIMER_RIGHT_0"]
        specs = [
            ("F", Lp, Ll, resE["PRIMER_LEFT_0_SEQUENCE"], resE["PRIMER_LEFT_0_TM"],
             resE["PRIMER_LEFT_0_GC_PERCENT"], resE["PRIMER_LEFT_0_SELF_ANY_TH"],
             resE["PRIMER_LEFT_0_SELF_END_TH"], resE["PRIMER_LEFT_0_HAIRPIN_TH"]),
            ("R", Rp, Rl, resE["PRIMER_RIGHT_0_SEQUENCE"], resE["PRIMER_RIGHT_0_TM"],
             resE["PRIMER_RIGHT_0_GC_PERCENT"], resE["PRIMER_RIGHT_0_SELF_ANY_TH"],
             resE["PRIMER_RIGHT_0_SELF_END_TH"], resE["PRIMER_RIGHT_0_HAIRPIN_TH"]),
        ]
        for role, p, ln, seq, tm, gc, sany, send, hp in specs:
            if role == "F":
                g0, g1 = gmap[p], gmap[p + ln - 1]
            else:
                g0, g1 = gmap[p - ln + 1], gmap[p]
            rows.append(dict(
                assay="long_junction(179939-specific)", role=role,
                oligo="forward" if role == "F" else "reverse", seq5to3=seq,
                chrom="chr5", g_start0=min(g0, g1), g_end0=max(g0, g1) + 1, strand="-",
                tm=round(tm, 2), gc=round(gc, 2),
                product="cDNA %dbp across intron 122454304-122455892" % resE["PRIMER_PAIR_0_PRODUCT_SIZE"],
                dimer_self_any=round(sany, 1), dimer_self_end=round(send, 1),
                hairpin=round(hp, 1)))

    cols = ["assay", "role", "oligo", "seq5to3", "chrom", "g_start0", "g_end0",
            "strand", "tm", "gc", "pro
