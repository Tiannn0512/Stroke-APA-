#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""调试 qPCR_common 候选"""
import sys

import primer3

sys.path.insert(0, "/mnt/d/stroke_apa_reaudit/scripts")
from g1_3_seq_utils import verify_window
from g1_3_specificity import load_transcripts, pair_amplicons

GENOME = "/mnt/d/stroke_apa_reaudit/reference/chr5.fa"
tx = load_transcripts()
ex1 = verify_window(GENOME, "chr5", 122489237, 122489376, "-")
ex2 = verify_window(GENOME, "chr5", 122491680, 122491785, "-")
tmpl = ex2 + ex1
print("模板长:", len(tmpl))
res = primer3.bindings.design_primers(
    {"SEQUENCE_ID": "c", "SEQUENCE_TEMPLATE": tmpl},
    {"PRIMER_TASK": "generic", "PRIMER_PICK_LEFT_PRIMER": 1, "PRIMER_PICK_RIGHT_PRIMER": 1,
     "PRIMER_PICK_INTERNAL_OLIGO": 0, "PRIMER_PRODUCT_SIZE_RANGE": [[120, 244]],
     "PRIMER_OPT_SIZE": 22, "PRIMER_MIN_SIZE": 20, "PRIMER_MAX_SIZE": 30,
     "PRIMER_OPT_TM": 60.0, "PRIMER_MIN_TM": 57.0, "PRIMER_MAX_TM": 63.0,
     "PRIMER_MIN_GC": 40.0, "PRIMER_MAX_GC": 60.0,
     "PRIMER_MAX_SELF_ANY_TH": 42.0, "PRIMER_MAX_SELF_END_TH": 32.0,
     "PRIMER_MAX_HAIRPIN_TH": 40.0, "PRIMER_NUM_RETURN": 20})
n = res.get("PRIMER_PAIR_NUM_RETURNED", 0)
print("候选数:", n)
for i in range(min(n, 4)):
    F = res[f"PRIMER_LEFT_{i}_SEQUENCE"]
    R = res[f"PRIMER_RIGHT_{i}_SEQUENCE"]
    amp = pair_amplicons(F, R, tx)
    print(f"cand{i} {F}/{R} product={res[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']}")
    for a in amp:
        print("   ", a)
