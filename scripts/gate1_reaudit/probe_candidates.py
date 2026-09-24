#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""调试：junction 与 common 各候选的扫描计数"""
import sys

import primer3

sys.path.insert(0, "/mnt/d/stroke_apa_reaudit/scripts")
from g1_3_seq_utils import verify_window
from g1_3_specificity import load_transcripts, pair_amplicons

GENOME = "/mnt/d/stroke_apa_reaudit/reference/chr5.fa"
tx = load_transcripts()
S = {"PRIMER_TASK": "generic", "PRIMER_PICK_LEFT_PRIMER": 1,
     "PRIMER_PICK_RIGHT_PRIMER": 1, "PRIMER_PICK_INTERNAL_OLIGO": 0,
     "PRIMER_OPT_SIZE": 22, "PRIMER_MIN_SIZE": 20, "PRIMER_MAX_SIZE": 30,
     "PRIMER_OPT_TM": 60.0, "PRIMER_MIN_TM": 57.0, "PRIMER_MAX_TM": 63.0,
     "PRIMER_MIN_GC": 40.0, "PRIMER_MAX_GC": 60.0,
     "PRIMER_MAX_SELF_ANY_TH": 42.0, "PRIMER_MAX_SELF_END_TH": 32.0,
     "PRIMER_MAX_HAIRPIN_TH": 40.0, "PRIMER_NUM_RETURN": 20}

# junction
seg1 = verify_window(GENOME, "chr5", 122453512, 122454304, "-")
seg2 = verify_window(GENOME, "chr5", 122455892, 122456500, "-")
cdna = seg1 + seg2
res = primer3.bindings.design_primers(
    {"SEQUENCE_ID": "j", "SEQUENCE_TEMPLATE": cdna},
    dict(S, PRIMER_PRODUCT_SIZE_RANGE=[[900, 1100]]))
n = res.get("PRIMER_PAIR_NUM_RETURNED", 0)
print("junction 候选数:", n)
for i in range(min(n, 5)):
    F = res[f"PRIMER_LEFT_{i}_SEQUENCE"]
    R = res[f"PRIMER_RIGHT_{i}_SEQUENCE"]
    amp = pair_amplicons(F, R, tx)
    print(f"  cand{i} {F}/{R} size={res[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']} amp={amp}")

# common
tmpl = verify_window(GENOME, "chr5", 122457511, 122457629, "+")
res2 = primer3.bindings.design_primers(
    {"SEQUENCE_ID": "c", "SEQUENCE_TEMPLATE": tmpl},
    dict(S, PRIMER_PRODUCT_SIZE_RANGE=[[70, 118]]))
n2 = res2.get("PRIMER_PAIR_NUM_RETURNED", 0)
print("common 候选数:", n2)
for i in range(min(n2, 6)):
    F = res2[f"PRIMER_LEFT_{i}_SEQUENCE"]
    R = res2[f"PRIMER_RIGHT_{i}_SEQUENCE"]
    amp = pair_amplicons(F, R, tx)
    on = sum(1 for a in amp if a[3] == "on")
    off = sum(1 for a in amp if a[3] == "OFF")
    print(f"  cand{i} {F}/{R} on={on} off={off} {[(a[0], a[2]) for a in amp if a[3]=='OFF']}")
