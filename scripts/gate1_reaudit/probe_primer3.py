#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""探查 primer3 RIGHT 引物坐标约定与 pick right-only 的正确姿势"""
import sys
import primer3

sys.path.insert(0, "/mnt/d/stroke_apa_reaudit/scripts")
from g1_3_seq_utils import verify_window

GENOME = "/mnt/d/stroke_apa_reaudit/reference/chr5.fa"

# qPCR common 外显子（+ 链模板，119bp）
tmpl = verify_window(GENOME, "chr5", 122457511, 122457629, "+")
print("template len:", len(tmpl))
r = primer3.bindings.design_primers(
    {"SEQUENCE_ID": "t", "SEQUENCE_TEMPLATE": tmpl},
    {"PRIMER_PICK_LEFT_PRIMER": 1, "PRIMER_PICK_RIGHT_PRIMER": 1,
     "PRIMER_PRODUCT_SIZE_RANGE": [[70, 119]], "PRIMER_NUM_RETURN": 1})
print("LEFT_0:", r["PRIMER_LEFT_0"], r["PRIMER_LEFT_0_SEQUENCE"])
print("RIGHT_0:", r["PRIMER_RIGHT_0"], r["PRIMER_RIGHT_0_SEQUENCE"])
# 验证 RIGHT 序列与模板的对应关系
pos, ln = r["PRIMER_RIGHT_0"]
seg = tmpl[pos:pos + ln]
print("tmpl[pos:pos+ln] =", seg)
print("rc(tmpl[pos:pos+ln]) =", seg.translate(str.maketrans("ACGT", "TGCA"))[::-1])
seg2 = tmpl[pos + ln - len(r["PRIMER_RIGHT_0_SEQUENCE"]): pos + ln]
print("模板最右端比对：", tmpl[-ln:], "| right seq:", r["PRIMER_RIGHT_0_SEQUENCE"])

# right-only 任务
r2 = primer3.bindings.design_primers(
    {"SEQUENCE_ID": "t", "SEQUENCE_TEMPLATE": tmpl},
    {"PRIMER_PICK_LEFT_PRIMER": 0, "PRIMER_PICK_RIGHT_PRIMER": 1,
     "PRIMER_NUM_RETURN": 3})
print("right-only num:", r2.get("PRIMER_PAIR_NUM_RETURNED"),
      [k for k in r2 if "RIGHT_0" in k])
