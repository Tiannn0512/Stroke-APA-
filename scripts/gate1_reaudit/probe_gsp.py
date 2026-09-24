#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""探查 RACE GSP right-only 设计为何返回空：逐级放宽约束"""
import sys

import primer3

sys.path.insert(0, "/mnt/d/stroke_apa_reaudit/scripts")
from g1_3_seq_utils import verify_window

GENOME = "/mnt/d/stroke_apa_reaudit/reference/chr5.fa"
tmpl = verify_window(GENOME, "chr5", 122457511, 122457629, "-")
print("len:", len(tmpl), tmpl[:40])


def attempt(tag, settings):
    r = primer3.bindings.design_primers(
        {"SEQUENCE_ID": "gsp", "SEQUENCE_TEMPLATE": tmpl}, settings)
    n = r.get("PRIMER_PAIR_NUM_RETURNED", 0)
    seq = r.get("PRIMER_RIGHT_0_SEQUENCE", "-")
    print(f"{tag:<28} n={n} right0={seq}")


base = {"PRIMER_TASK": "generic", "PRIMER_PICK_LEFT_PRIMER": 0,
        "PRIMER_PICK_RIGHT_PRIMER": 1, "PRIMER_PICK_INTERNAL_OLIGO": 0,
        "PRIMER_PRODUCT_SIZE_RANGE": [[30, 118]]}
attempt("defaults+right-only", dict(base))
attempt("+Tm58-66", dict(base, PRIMER_OPT_TM=62.0, PRIMER_MIN_TM=58.0,
                         PRIMER_MAX_TM=66.0))
attempt("+size22-30", dict(base, PRIMER_OPT_TM=62.0, PRIMER_MIN_TM=58.0,
                           PRIMER_MAX_TM=66.0, PRIMER_MIN_SIZE=22,
                           PRIMER_OPT_SIZE=26, PRIMER_MAX_SIZE=30))
attempt("+GC40-60", dict(base, PRIMER_OPT_TM=62.0, PRIMER_MIN_TM=58.0,
                         PRIMER_MAX_TM=66.0, PRIMER_MIN_SIZE=22,
                         PRIMER_OPT_SIZE=26, PRIMER_MAX_SIZE=30,
                         PRIMER_MIN_GC=40.0, PRIMER_MAX_GC=60.0))
attempt("+TH42/32/40", dict(base, PRIMER_OPT_TM=62.0, PRIMER_MIN_TM=58.0,
                            PRIMER_MAX_TM=66.0, PRIMER_MIN_SIZE=22,
                            PRIMER_OPT_SIZE=26, PRIMER_MAX_SIZE=30,
                            PRIMER_MIN_GC=40.0, PRIMER_MAX_GC=60.0,
                            PRIMER_MAX_SELF_ANY_TH=42.0,
                            PRIMER_MAX_SELF_END_TH=32.0,
                            PRIMER_MAX_HAIRPIN_TH=40.0))
# 对照：双引物设计在同一模板
r = primer3.bindings.design_primers(
    {"SEQUENCE_ID": "gsp", "SEQUENCE_TEMPLATE": tmpl},
    {"PRIMER_PICK_LEFT_PRIMER": 1, "PRIMER_PICK_RIGHT_PRIMER": 1,
     "PRIMER_PRODUCT_SIZE_RANGE": [[70, 118]], "PRIMER_OPT_TM": 62.0,
     "PRIMER_MIN_TM": 58.0, "PRIMER_MAX_TM": 66.0, "PRIMER_MIN_SIZE": 22,
     "PRIMER_OPT_SIZE": 26, "PRIMER_MAX_SIZE": 30})
print("pair-design n=", r.get("PRIMER_PAIR_NUM_RETURNED"),
      "| RIGHT_0:", r.get("PRIMER_RIGHT_0"), r.get("PRIMER_RIGHT_0_SEQUENCE"))
