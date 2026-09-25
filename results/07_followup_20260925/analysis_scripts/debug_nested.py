#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""调试：nested GSP 候选被哪个过滤器卡死。"""
import json
import primer3
from pyfaidx import Fasta

BASE = "/mnt/d/stroke_apa_reaudit"
GENOME = Fasta(f"{BASE}/reference/chr5.fa")["chr5"]
LO, HI = 122400000, 122510000
LOCUS = str(GENOME[LO:HI]).upper()
SEQS, cur = {}, None
with open(f"{BASE}/reference/gencode.vM25.transcripts.fa") as fh:
    for line in fh:
        if line.startswith(">"):
            cur = line[1:].split("|")[0].strip(); SEQS[cur] = []
        else:
            SEQS[cur].append(line.strip())
SEQS = {k: "".join(v).upper() for k, v in SEQS.items()}
L_, S_, M_, T204, T205, GMID = ("ENSMUST00000179939.7", "ENSMUST00000177974.7",
                                "ENSMUST00000031423.9", "ENSMUST00000196490.1",
                                "ENSMUST00000197415.4", "ENSMUST00000227710.1")
GM = SEQS[GMID]

def rc(s):
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))

def mm_hits(q, text, maxmm=0):
    hits, Ln = [], len(q)
    for i in range(len(text) - Ln + 1):
        mm = 0
        for j in range(Ln):
            if text[i + j] != q[j]:
                mm += 1
                if mm > maxmm: break
        if mm <= maxmm:
            hits.append(i)
    return hits

g_hi1, g_lo1 = 122457629, 122457511
g_hi2, g_lo2 = 122457426, 122457302
tmpl = rc(str(GENOME[g_lo1:g_hi1])) + rc(str(GENOME[g_lo2:g_hi2]))
res = primer3.bindings.design_primers(
    {"SEQUENCE_ID": "x", "SEQUENCE_TEMPLATE": tmpl},
    dict(PRIMER_OPT_SIZE=23, PRIMER_MIN_SIZE=20, PRIMER_MAX_SIZE=28,
         PRIMER_OPT_TM=63.0, PRIMER_MIN_TM=60.0, PRIMER_MAX_TM=66.0,
         PRIMER_MIN_GC=35.0, PRIMER_MAX_GC=65.0,
         PRIMER_PICK_LEFT_PRIMER=1, PRIMER_PICK_RIGHT_PRIMER=1,
         PRIMER_PICK_INTERNAL_OLIGO=0, PRIMER_NUM_RETURN=200,
         PRIMER_PRODUCT_SIZE_RANGE=[[50, len(tmpl)]]))
n = res.get("PRIMER_PAIR_NUM_RETURNED", 0)
print("candidates:", n)
stats = dict(tm=0, dims=0, main_hits=0, minor_hits=0, gm=0, locus0=0, strand=0, far=0, ok=0)
for i in range(n):
    pos, ln = res[f"PRIMER_LEFT_{i}"]
    seq = res[f"PRIMER_LEFT_{i}_SEQUENCE"]
    tm = round(primer3.bindings.calc_tm(seq), 2)
    if not (60 <= tm <= 65):
        stats["tm"] += 1; continue
    sa = round(primer3.bindings.calc_homodimer(seq).tm, 1)
    se = 0
    hp = round(primer3.bindings.calc_hairpin(seq).tm, 1)
    if sa >= 45 or hp >= 45:
        stats["dims"] += 1; continue
    per = {}
    for t in (L_, S_, M_, T204, T205, GMID):
        per[t] = ([(p, "+") for p in mm_hits(seq, SEQS[t])]
                  + [(p, "-") for p in mm_hits(rc(seq), SEQS[t])])
    if any(len(per[t]) != 1 or per[t][0][1] != "+" for t in (L_, S_, M_)):
        stats["main_hits"] += 1
        print(f"  cand{i} {seq} main={[ (t.split('.')[0], per[t]) for t in (L_,S_,M_) if per[t] or True]}")
        continue
    if any(per[t] for t in (T204, T205, GMID)):
        stats["minor_hits"] += 1; continue
    best_mm = 99
    for q in (seq, rc(seq)):
        Ln = len(q)
        for j in range(len(GM) - Ln + 1):
            mm = sum(1 for a, b in zip(q, GM[j:j+Ln]) if a != b)
            best_mm = min(best_mm, mm)
    if best_mm < 2:
        stats["gm"] += 1
        print(f"  cand{i} {seq} gm_min_mm={best_mm}")
        continue
    lh0 = ([(LO + p, "+") for p in mm_hits(seq, LOCUS)]
           + [(LO + p, "-") for p in mm_hits(rc(seq), LOCUS)])
    if len(lh0) != 1:
        stats["locus0"] += 1
        print(f"  cand{i} {seq} locus0={lh0}")
        continue
    gpos, strand = lh0[0]
    if not (g_lo2 - 30 <= gpos <= g_hi1 + 30) or strand != "-":
        stats["strand"] += 1
        print(f"  cand{i} {seq} gpos={gpos} strand={strand}")
        continue
    far = [x for x in ([(LO + p, "+") for p in mm_hits(seq, LOCUS, 1)]
                       + [(LO + p, "-") for p in mm_hits(rc(seq), LOCUS, 1)])
          if not (gpos - 60 <= x[0] <= gpos + 60)]
    if far:
        stats["far"] += 1
        print(f"  cand{i} {seq} far_le1mm={far[:4]}")
        continue
    stats["ok"] += 1
    print(f"  OK cand{i} {seq} tm={tm} gpos={gpos} dist3={len(tmpl)-(pos+ln)}")
print(stats)
