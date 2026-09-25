#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
f4 (v2): Atp2a2 检测体系全量复核与重设计。

05 号方案指出的 3 个 HOLD 项 + 本轮新发现的 2 个结构缺陷：
  H1 long_junction 配对扫描 0 扩增子 → 根因：06 设计模板把终端外显子段拼在
     倒数第二外显子段之前（旋转模板），两引物在真实 mRNA 上背向。重设计 juncL_v2。
  H2 outer GSP 表/文档序列不一致 → canonical = GACAGGAA...（基因组唯一验证），
     文档 TCAAGAAG... 为早期泄漏（chr5:122457534），勘误。
  H3 nested GSP 命中 Gm30970（exact）→ 重设计，要求 minMM>=2。
  H4 (新) RACE outer/nested GSP 均为 mRNA 反义向（'-'位点，向 mRNA 5' 延伸）——
     3'RACE 的 GSP 必须正义向（'+'，向 poly(A) 延伸）。修复 = 取反向互补（Tm/GC 不变，
     脱靶重扫）。
  H5 (新) long_junction "366bp cDNA 产物" 在真实转录本上不成立（背向）。
     新对 juncL_v2：F 在 L 独有段（mRNA-sense 正义），R 在终端外显子（反义）。

输出：assay_design_verified.tsv + primer_pair_products.tsv + f4_design_log.json
方向约定：oligo 自身序列出现在 mRNA sense = '+'（右向延伸）；rc 后出现 = '-'（左向）。
isPCR 规则：'+'位点 pF 与 '-'位点 pR（rc(R) 匹配 sense 的窗口起点），
pF < pR 时产物 = pR+lenR-pF；镜像同理。
"""
import json

import primer3
from pyfaidx import Fasta

BASE = "/mnt/d/stroke_apa_reaudit"
OUT = f"{BASE}/results/07_followup_20260925"
FACTS = json.load(open(f"{OUT}/analysis_scripts/structure_facts.json"))
GENOME = Fasta(f"{BASE}/reference/chr5.fa")["chr5"]

LO, HI = 122400000, 122510000
LOCUS = str(GENOME[LO:HI]).upper()

SEQS, cur = {}, None
with open(f"{BASE}/reference/gencode.vM25.transcripts.fa") as fh:
    for line in fh:
        if line.startswith(">"):
            cur = line[1:].split("|")[0].strip()
            SEQS[cur] = []
        else:
            SEQS[cur].append(line.strip())
SEQS = {k: "".join(v).upper() for k, v in SEQS.items()}
L_, S_, M_, T204, T205, GMID = ("ENSMUST00000179939.7", "ENSMUST00000177974.7",
                                "ENSMUST00000031423.9", "ENSMUST00000196490.1",
                                "ENSMUST00000197415.4", "ENSMUST00000227710.1")
TIDS = [L_, S_, M_, T204, T205, GMID]
for t in TIDS:
    assert t in SEQS, t
GM = SEQS[GMID]
MAIN = [L_, S_, M_]

def rc(s):
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))

def mm_hits(q, text, maxmm=0):
    hits, Ln = [], len(q)
    for i in range(len(text) - Ln + 1):
        mm = 0
        for j in range(Ln):
            if text[i + j] != q[j]:
                mm += 1
                if mm > maxmm:
                    break
        if mm <= maxmm:
            hits.append(i)
    return hits

def all_hits(oligo, maxmm=0):
    res = {}
    for t in TIDS:
        res[t] = ([(p, "+") for p in mm_hits(oligo, SEQS[t], maxmm)]
                  + [(p, "-") for p in mm_hits(rc(oligo), SEQS[t], maxmm)])
    return res

def locus_hits(oligo, maxmm=0):
    return ([(LO + p, "+") for p in mm_hits(oligo, LOCUS, maxmm)]
            + [(LO + p, "-") for p in mm_hits(rc(oligo), LOCUS, maxmm)])

def gm_min_mm(oligo):
    best = 99
    for q in (oligo, rc(oligo)):
        Ln = len(q)
        for i in range(len(GM) - Ln + 1):
            mm = sum(1 for a, b in zip(q, GM[i:i + Ln]) if a != b)
            best = min(best, mm)
    return best

def products(f, r, max_len=8000):
    """isPCR 双向规则：返回 (space, id, prod_bp, detail)。space: cDNA:<tid> / genome。"""
    out = []
    for t in TIDS:
        s = SEQS[t]
        for a in mm_hits(f, s):
            for b in mm_hits(rc(r), s):
                pr = b + len(r) - a
                if 0 < pr <= max_len:
                    out.append((f"cDNA:{t}", pr, f"{a+1}-{b+len(r)}"))
        for b in mm_hits(r, s):
            for a in mm_hits(rc(f), s):
                pr = a + len(f) - b
                if 0 < pr <= max_len and b < a:
                    out.append((f"cDNA:{t}", pr, f"{b+1}-{a+len(f)}"))
    fa, rm = mm_hits(rc(f), LOCUS), mm_hits(rc(r), LOCUS)
    for a in mm_hits(f, LOCUS):
        for b in fa:
            pr = b + len(f) - a
            if 0 < pr <= max_len and a < b:
                out.append(("genome", pr, f"chr5:{LO+a}-{LO+b+len(f)} [F+/R-]"))
    for b in mm_hits(r, LOCUS):
        for a in fa:
            pr = a + len(f) - b
            if 0 < pr <= max_len and b < a:
                out.append(("genome", pr, f"chr5:{LO+b}-{LO+a+len(f)} [F-/R+]"))
    for a in mm_hits(f, GM):
        for b in mm_hits(rc(r), GM):
            pr = b + len(r) - a
            if 0 < pr <= max_len:
                out.append(("Gm30970", pr, f"{a+1}-{b+len(r)}"))
    for b in mm_hits(r, GM):
        for a in mm_hits(rc(f), GM):
            pr = a + len(f) - b
            if 0 < pr <= max_len and b < a:
                out.append(("Gm30970", pr, f"{b+1}-{a+len(f)}"))
    seen, res = set(), []
    for row in out:
        if row not in seen:
            seen.add(row)
            res.append(row)
    return res

def tm_of(seq):
    return round(primer3.bindings.calc_tm(seq), 2)

def gc_of(seq):
    seq = seq.upper()
    return round(100.0 * (seq.count("G") + seq.count("C")) / len(seq), 2)

def dims(seq):
    return (round(primer3.bindings.calc_homodimer(seq).tm, 1),
            round(primer3.bindings.calc_hairpin(seq).tm, 1))

def seed6_offhits(seq):
    seed = seq[-6:]
    n = LOCUS.count(seed) + LOCUS.count(rc(seed))
    for t in TIDS:
        n += SEQS[t].count(seed) + SEQS[t].count(rc(seed))
    return n

def end4_conflict(a, b):
    return a.endswith(rc(b[-4:])) or b.endswith(rc(a[-4:]))

def describe(oligo):
    h = all_hits(oligo)
    per = {t: h[t] for t in TIDS}
    gm = gm_min_mm(oligo)
    lh0, lh1 = locus_hits(oligo), locus_hits(oligo, 1)
    return per, gm, lh0, lh1

# ============ 交付 oligo 与修复 ============
DELIVERED = {
    "RACE_outer_GSP":  "GACAGGAATTAAAGCCTCAGGAAACC",
    "RACE_nested_GSP": "TGAGGGCATTACACATCTCTATGGTT",
    "qPCR_common_F":   "AGCCTTTGTAGAGCCGTTTGTA",
    "qPCR_common_R":   "CACACTCTTTCTGTCCTGTCGA",
    "distal_F":        "AGTTAGGACTGGAGGCCTATGT",
    "distal_R":        "TTGTAAGTGGCCAGATTGCTCT",
    "juncL_F":         "GGAGTAACCGCTTCCTAAACCA",
    "juncL_R":         "ATCCTACTATGCGGCAGAACAG",
}
print("=== 06 包 oligo 在真实 mRNA 上的方向（'+'=正义/右向）===")
for name, seq in DELIVERED.items():
    h = all_hits(seq)
    for t in MAIN:
        if h[t]:
            print(f"  {name:16s} {t.split('.')[0]}: "
                  f"{[f'{p+1}{s}' for p, s in h[t]]}")
    if h[GMID]:
        print(f"  {name:16s} Gm30970: {h[GMID]}")

# H4 修复：GSP 取反向互补
OUTER_V2 = rc(DELIVERED["RACE_outer_GSP"])
per, gm, lh0, lh1 = describe(OUTER_V2)
print(f"outer_GSP_v2(=rc): on-main={{{t.split('.')[0]}:[... for t in MAIN]}}"
      if False else
      f"outer_GSP_v2 hits: " + "; ".join(
          f"{t.split('.')[0]}:{[(p+1, s) for p, s in per[t]]}" for t in TIDS if per[t]))
print(f"  gm30970_minMM={gm} locus_exact={lh0} locus_le1mm={len(lh1)}")
assert len(per[L_]) == len(per[S_]) == len(per[M_]) == 1, "outer v2 主转录本命中异常"
assert all(per[t][0][1] == "+" for t in MAIN), "outer v2 应为正义向"
assert not any(per[t] for t in (T204, T205, GMID)), "outer v2 有额外命中"
assert gm >= 2, "outer v2 Gm30970 错配不足"
assert len(lh0) == 1, "outer v2 基因组命中不唯一"

# ============ H3/H4: nested GSP 重设计（LEFT primer，正义向） ============
print("=== nested GSP 重设计（共有块 122457511-457629 + 122457302-457426, mRNA 感知）===")
# 干净共有窗口（Gm30970 同源区 = L_mRNA 3157-3519 之外）：
#   exon 122458447-458533（86bp 全净）与 exon 122458087-122458159（72bp 净段）
g_hi1, g_lo1 = 122458533, 122458447
g_hi2, g_lo2 = 122458159, 122458087
tmpl = rc(str(GENOME[g_lo1:g_hi1])) + rc(str(GENOME[g_lo2:g_hi2]))

def tpos_to_g(i):
    return (g_hi1 - 1 - i) if i < (g_hi1 - g_lo1) else (g_hi2 - 1 - (i - (g_hi1 - g_lo1)))

NB1 = g_hi1 - g_lo1
res = primer3.bindings.design_primers(
    {"SEQUENCE_ID": "nested_gsp_common", "SEQUENCE_TEMPLATE": tmpl},
    dict(PRIMER_OPT_SIZE=23, PRIMER_MIN_SIZE=20, PRIMER_MAX_SIZE=28,
         PRIMER_OPT_TM=63.0, PRIMER_MIN_TM=60.0, PRIMER_MAX_TM=66.0,
         PRIMER_MIN_GC=35.0, PRIMER_MAX_GC=65.0,
         PRIMER_PICK_LEFT_PRIMER=1, PRIMER_PICK_RIGHT_PRIMER=1,
         PRIMER_PICK_INTERNAL_OLIGO=0, PRIMER_NUM_RETURN=200,
         PRIMER_PRODUCT_SIZE_RANGE=[[50, len(tmpl)]]))
n = res.get("PRIMER_PAIR_NUM_RETURNED", 0)
cands, best = [], None
for i in range(n):
    pos, ln = res[f"PRIMER_LEFT_{i}"]          # 占据 [pos, pos+ln)
    seq = res[f"PRIMER_LEFT_{i}_SEQUENCE"]
    if not (60 <= tm_of(seq) <= 65) or end4_conflict(seq, OUTER_V2):
        continue
    dm, hp = dims(seq)
    if dm >= 45 or hp >= 45:
        continue
    per, gm, lh0, lh1 = describe(seq)
    if any(len(per[t]) != 1 or per[t][0][1] != "+" for t in MAIN):
        continue
    if any(per[t] for t in (T204, T205, GMID)):
        continue
    if gm < 2 or len(lh0) != 1:
        continue
    gpos = lh0[0][0]
    if not (g_lo2 <= gpos <= g_hi1) or lh0[0][1] != "-":
        continue
    far = [x for x in locus_hits(seq, 1) if not (gpos - 60 <= x[0] <= gpos + 60)]
    if far:
        continue
    c = dict(seq=seq, tm=tm_of(seq), gc=gc_of(seq), gm=gm, dims=dims(seq),
             dist3=len(tmpl) - (pos + ln),
             g=f"{gpos}-{gpos+ln}")
    cands.append(c)
    score = c["dist3"] + abs(c["tm"] - 63) * 5
    c["score"] = score
    if best is None or score < best["score"]:
        best = c
assert best, "nested GSP 无合格候选"
NESTED_V2 = best["seq"]
print(f"candidates={len(cands)} chosen={best}")

# ============ H1/H5: long_junction v2 ============
print("=== long_junction v2（F 在 L 独有段，R 在终端外显子）===")
segL = rc(str(GENOME[122455892:122456339]))       # L 独有段 mRNA sense (447bp)
segT = rc(str(GENOME[122453512:122454600]))       # 终端外显子+下游 88bp mRNA sense
tmplL = segL + segT

def lpos_to_g(i):
    if i < len(segL):
        return 122456338 - i
    return 122454599 - (i - len(segL))

res = primer3.bindings.design_primers(
    {"SEQUENCE_ID": "juncL_v2", "SEQUENCE_TEMPLATE": tmplL},
    dict(PRIMER_OPT_SIZE=23, PRIMER_MIN_SIZE=20, PRIMER_MAX_SIZE=28,
         PRIMER_OPT_TM=63.0, PRIMER_MIN_TM=60.0, PRIMER_MAX_TM=66.0,
         PRIMER_MIN_GC=35.0, PRIMER_MAX_GC=65.0,
         PRIMER_PICK_LEFT_PRIMER=1, PRIMER_PICK_RIGHT_PRIMER=1,
         PRIMER_PICK_INTERNAL_OLIGO=0, PRIMER_NUM_RETURN=500,
         PRIMER_PRODUCT_SIZE_RANGE=[[120, 500]]))
n = res.get("PRIMER_PAIR_NUM_RETURNED", 0)
bestL, candsL = None, []
for i in range(n):
    fseq = res[f"PRIMER_LEFT_{i}_SEQUENCE"]
    rseq = res[f"PRIMER_RIGHT_{i}_SEQUENCE"]
    fpos, fln = res[f"PRIMER_LEFT_{i}"]
    rpos, rln = res[f"PRIMER_RIGHT_{i}"]
    # F 必须整体在 L 独有段内；R 必须整体在终端外显子段内
    if not (fpos + fln <= len(segL) and (rpos - rln + 1) >= len(segL)):
        continue
    if end4_conflict(fseq, rseq):
        continue
    # F: L 独有段内正义单命中；R: 终端外显子反义，L/S 各一命中，M 无
    def leg_ok(seq, want):
        # want = {tid: orientation} 精确命中要求；其余模板必须零命中
        per, gm, lh0, lh1 = describe(seq)
        for t in TIDS:
            if t in want:
                if len(per[t]) != 1 or per[t][0][1] != want[t]:
                    return False
            elif per[t]:
                return False
        if gm < 2 or len(lh0) != 1:
            return False
        return dims(seq)[0] < 45 and dims(seq)[1] < 45
    if not (leg_ok(fseq, {L_: "+"}) and leg_ok(rseq, {L_: "-", S_: "-"})):
        continue
    c = dict(F=fseq, R=rseq, prod=res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"],
             F_tm=tm_of(fseq), R_tm=tm_of(rseq))
    candsL.append(c)
    if bestL is None or abs(c["F_tm"] - c["R_tm"]) < abs(bestL["F_tm"] - bestL["R_tm"]):
        bestL = c
assert bestL, "long_junction v2 无合格候选"
JF_V2, JR_V2 = bestL["F"], bestL["R"]
print(f"cand={len(candsL)} chosen: F={JF_V2} R={JR_V2} prod={bestL['prod']}")

# ============ J_S / J_L 接头嵌合引物（+' 向，与 distal_F 配对） ============
def junction_oligo(exon_lo, term_hi, half=13, shift=0):
    """mRNA 感知跨接头 oligo。
    donor 侧：penult 外显子低坐标端 exon_lo 起的 a 个碱基（mRNA 末 a 位，基因组升序）；
    acceptor 侧：terminal 外显子高坐标端 term_hi 起的 b 个碱基（mRNA 首 b 位，基因组降序）。"""
    a, b = half + shift, half - shift
    left = rc(str(GENOME[exon_lo: exon_lo + a]).upper())
    right = rc(str(GENOME[term_hi - b + 1: term_hi + 1]).upper())
    return left + right

def junction_accept(seq, target):
    per, gm, lh0, lh1 = describe(seq)
    ok = (len(per[target]) == 1 and per[target][0][1] == "+"
          and all(not per[t] for t in TIDS if t != target)
          and gm >= 2 and len(lh0) == 0 and len(lh1) == 0
          and dims(seq)[0] < 45 and dims(seq)[1] <= 47
          and 58 <= tm_of(seq) <= 68)
    return ok, dict(per={t.split('.')[0]: len(per[t]) for t in TIDS}, gm=gm,
                    lh0=lh0, lh1=len(lh1), tm=tm_of(seq), dims=dims(seq))

def pick_junction(penult_lo, term_hi, target, label):
    seq = ok = d = None
    for half in (12, 11, 13, 14):
        for shift in (0, -2, 2, -3, 3, -1, 1, -4, 4, -5, 5):
            seq = junction_oligo(penult_lo, term_hi, half, shift)
            ok, d = junction_accept(seq, target)
            if ok:
                print(f"{label} chosen: {seq} half={half} shift={shift}")
                return seq, ok, d
    return seq, ok, d

J_S, okS, dS = pick_junction(122457302, 122454303, S_, "J_S")
J_L, okL, dL = pick_junction(122455892, 122454303, L_, "J_L")
print("J_S", okS, dS)
print("J_L", okL, dL)
assert okS and okL

# ============ 配对产物总表 ============
PAIRS = [
    ("qPCR_total(common)", DELIVERED["qPCR_common_F"], DELIVERED["qPCR_common_R"]),
    ("distal_terminal(LS)", DELIVERED["distal_F"], DELIVERED["distal_R"]),
    ("long_junction_v2(L)", JF_V2, JR_V2),
    ("qPCR_S_spliceform(J_S+distal_F)", J_S, DELIVERED["distal_F"]),
    ("qPCR_L_orthogonal(J_L+distal_F)", J_L, DELIVERED["distal_F"]),
]
print("=== 配对产物（cDNA×转录本 / 基因组 / Gm30970）===")
all_products = {}
for name, f, r in PAIRS:
    pr = products(f, r)
    all_products[name] = pr
    print(f"-- {name}")
    for sp, L, det in pr:
        print(f"   {sp}: {L}bp  {det}")

# ============ 输出表 ============
G2C = FACTS["transcripts"]

def mrna_pos(tid, g0):
    for lo, hi, start in G2C[tid]["g2c_blocks"]:
        if lo <= g0 < hi:
            return start + ((hi - 1) - g0)
    return None

FINAL = []
def add(assay, role, seq, status, note, g=None):
    lh0 = locus_hits(seq)
    if len(lh0) == 1:
        g0, strand = lh0[0]
    else:
        assert g is not None, (assay, lh0)
        g0, strand = g, "junction-spanning"
    per, gm, _, _ = describe(seq)
    mpos = "; ".join(
        f"{t.split('.')[0]}:" + (str(per[t][0][0] + 1) if per[t] else "NA")
        for t in MAIN)
    dm, hp = dims(seq)
    FINAL.append(dict(assay=assay, role=role, seq5to3=seq, chrom="chr5",
                      g_start0=g0, g_end0=g0 + len(seq), genome_strand_hit=strand,
                      extend_dir=("rightward->polyA(mRNA-sense)" if strand == "-"
                                  else "leftward(<-polyA,mRNA-antisense)" if strand == "+"
                                  else "rightward->polyA(junction-spanning)"),
                      on_LSM_mRNApos1=mpos, tm=tm_of(seq), gc=gc_of(seq),
                      homodimer_tm=dm, hairpin_tm=hp, gm30970_min_mm=gm,
                      seed6_offhits=seed6_offhits(seq), status=status, note=note))

add("RACE_outer_GSP_v2", "GSP(sense->polyA)", OUTER_V2, "verified(orientation-fixed)",
    "H4 修复：06 包 GACAGGAA... 为反义向(向mRNA5'), 弃用；v2=其反向互补, Tm/GC不变, 脱靶重扫通过")
add("RACE_outer_GSP_v1_deprecated", "GSP(antisense,错误方向)",
    DELIVERED["RACE_outer_GSP"], "deprecated(orientation)",
    "保留记录；不得订购")
add("RACE_nested_GSP_v2", "GSP(sense->polyA)", NESTED_V2, "verified(redesigned)",
    f"H3+H4 修复：Gm30970 minMM={best['gm']}；共有块 LSM；正义向")
add("RACE_nested_GSP_v1_deprecated", "GSP(antisense,错误方向+Gm30970)",
    DELIVERED["RACE_nested_GSP"], "deprecated(orientation+Gm30970)", "保留记录；不得订购")
add("qPCR_total_F", "F(sense)", DELIVERED["qPCR_common_F"], "verified",
    "5 条注释转录本全部扩增 154bp；由'qPCR_common(3 transcripts)'更名")
add("qPCR_total_R", "R(antisense)", DELIVERED["qPCR_common_R"], "verified", "同上")
add("distal_terminal_F(mRNA-space:R)", "R(antisense on mRNA)", DELIVERED["distal_F"],
    "verified", "终端外显子 LS 共有；设计模板为基因组+链,故 mRNA 空间角色为 reverse")
add("distal_terminal_R(mRNA-space:F)", "F(sense on mRNA)", DELIVERED["distal_R"],
    "verified", "同上；cDNA 产物 137bp(L+S), 基因组 136bp")
add("long_junction_v2_F", "F(sense,L-unique)", JF_V2, "verified(redesigned)",
    "H1/H5 修复：v1 对背向无产物；v2 F 位于 L 独有段(122455892-122456339)")
add("long_junction_v2_R", "R(antisense,terminal exon)", JR_V2, "verified(redesigned)",
    "v2 R 位于终端外显子；cDNA 产物仅 L, 基因组产物跨 1588bp 内含子可区分")
add("long_junction_v1_F_deprecated", "F(旧,背向)", DELIVERED["juncL_F"],
    "deprecated(pair orientation)", "v1 对在真实 mRNA 上无扩增子")
add("long_junction_v1_R_deprecated", "R(旧,背向)", DELIVERED["juncL_R"],
    "deprecated(pair orientation)", "同上")
add("qPCR_S_spliceform_J", "F(junction-spanning)", J_S, "verified(in-silico)",
    "跨 S 独有接头(donor122457302->acceptor122454304)；识别 S 剪接型, 3'端身份待 3'RACE",
    g=122457302)
add("qPCR_L_orthogonal_J", "F(junction-spanning)", J_L, "verified(in-silico)",
    "跨 L 独有接头；与 long_junction_v2 正交, 可作 smFISH 探针基础",
    g=122455892)

cols = ["assay", "role", "seq5to3", "chrom", "g_start0", "g_end0",
        "genome_strand_hit", "extend_dir", "on_LSM_mRNApos1", "tm", "gc",
        "homodimer_tm", "hairpin_tm", "gm30970_min_mm", "seed6_offhits",
        "status", "note"]
with open(f"{OUT}/assay_design_verified.tsv", "w") as out:
    out.write("\t".join(cols) + "\n")
    for r in FINAL:
        out.write("\t".join(str(r[c]) for c in cols) + "\n")

with open(f"{OUT}/primer_pair_products.tsv", "w") as out:
    out.write("pair\ttemplate\tproduct_bp\tdetail\n")
    for name, pr in all_products.items():
        for sp, L, det in pr:
            out.write(f"{name}\t{sp}\t{L}\t{det}\n")

with open(f"{OUT}/analysis_scripts/f4_design_log.json", "w") as out:
    json.dump(dict(nested_candidates=len(cands), nested_chosen=best,
                   juncL_v2_candidates=len(candsL), juncL_v2_chosen=bestL,
                   J_S=dS, J_L=dL,
                   doc_leak="chr5:122457534-122457560(+)",
                   defects=dict(H1="rotated design template -> pair faces outward",
                                H2="doc leaked old candidate seq",
                                H3="nested GSP exact Gm30970 hit",
                                H4="RACE GSPs antisense-oriented (should be sense)",
                                H5="claimed 366bp cDNA product not real"),
                   params=dict(size="20-28 opt23", tm="60-66 opt63", gc="35-65",
                               maxmm_local=1, min_gm30970_mm=2)), out, indent=1,
              ensure_ascii=False)

print("WRITTEN tables.")
print("OUTER_V2 =", OUTER_V2)
print("NESTED_V2 =", NESTED_V2)
print("J_S =", J_S)
print("J_L =", J_L)
