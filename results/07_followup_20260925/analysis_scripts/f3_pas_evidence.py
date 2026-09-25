#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
f3: Atp2a2 候选/注释 poly(A) 位点证据表 PAS_evidence_table.tsv

证据源：
  S1 GENCODE vM25 注释 3'端（主3转录本 exon 低坐标端）
  S2 PolyASite 2.0 GRCm38.96 atlas（3'端测序簇，下载 2026-09-25，
     https://polyasite.unibas.ch/download/atlas/2.0/GRCm38.96/atlas.clusters.2.0.GRCm38.96.bed.gz）
  S3 DaPars2 event 6953 预测近端 PAS（发现队列，计算推断）
对每个位点计算：
  - polyA 信号基序（mRNA 感知上游 -45..-10 内 AATAAA/ATTAAA 及单碱基变体）
  - 下游 20nt A 含量与连续 A（internal priming 风险）
  - 与三个注释 3' 端 / 预测 PAS 的距离
判定规则（写进 status）：
  annotated+external_support = GENCODE 端且有 PolyASite 簇重叠（signal>=1 或 n>=5）
  cluster_only / noise_only / predicted_window_boundary 等
坐标 BED 0-based half-open；负链基因，断裂点上游( mRNA 感知) = 基因组高坐标方向。
"""
import gzip
import json

FA = "/mnt/d/stroke_apa_reaudit/reference/chr5.fa"
BED = "/mnt/d/stroke_apa_reaudit/results/07_followup_20260925/notes/pas_cache/polyasite_atlas_GRCm38.bed.gz"
OUT = "/mnt/d/stroke_apa_reaudit/results/07_followup_20260925/PAS_evidence_table.tsv"

from pyfaidx import Fasta
fa = Fasta(FA)
CH = fa["chr5"]

def gseq(s0, e0):            # + 链序列
    return str(CH[s0:e0]).upper()

def rc(s):
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))

def mrna_downstream(cut0, n=20):     # 负链：mRNA 下游 = 基因组低坐标
    return rc(gseq(cut0 - n, cut0))

def mrna_upstream(cut0, n=45):       # 负链：mRNA 上游 = 基因组高坐标
    return rc(gseq(cut0, cut0 + n))

CANON = {"AATAAA": 0, "ATTAAA": 0}
VARIANTS = {"AATATA", "AATACA", "AATAGA", "AATGAA", "ACTAAA", "AAGAAA", "AATAAT",
            "CATAAA", "TATAAA", "GATAAA", "AATAAA", "AACAAA", "AATCAA", "AGTAAA",
            "AACAAG", "AATGTA", "AOTTAA"}

def signal_scan(cut0):
    """在 mRNA 感知上游 -45..-10（基因组高坐标 10..45）扫描信号。"""
    up = mrna_upstream(cut0, 45)          # index0 = cut+1 (mRNA -1)
    best = None
    for motif, exact in [(m, True) for m in CANON] + [(m, False) for m in VARIANTS]:
        # mRNA 坐标 -45..-10 → up 的 index: 位置 -k → index k-1
        for k in range(45, 9, -1):
            seg = up[45 - k:45 - k + 6]
            if len(seg) < 6:
                continue
            if seg == motif and (exact or best is None or best[1] != "canonical"):
                kind = "canonical" if motif in CANON else "variant"
                if best is None or (kind == "canonical" and best[1] != "canonical") or \
                   (kind == "canonical" and best[1] == "canonical" and abs(k - 20) < abs(best[0] - 20)):
                    best = (k, kind, motif)
                elif best is not None and kind == "variant" and best[1] == "variant" and \
                        abs(k - 20) < abs(best[0] - 20):
                    best = (k, kind, motif)
    return best

def internal_priming(cut0):
    ds = mrna_downstream(cut0, 30)
    a20 = ds[:20]
    pct = 100.0 * a20.count("A") / 20
    run = max((len(r) for r in a20.split("T") for r in r.split("C") for r in r.split("G")),
              default=0)
    maxrun = 0
    cur = 0
    for ch in a20:
        cur = cur + 1 if ch == "A" else 0
        maxrun = max(maxrun, cur)
    risk = "HIGH" if (maxrun >= 6 or pct >= 60) else ("MEDIUM" if maxrun >= 4 or pct >= 40 else "LOW")
    return pct, maxrun, risk, ds[:20]

# PolyASite 簇（chr5 目标窗）
clusters = []
with gzip.open(BED, "rt") as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if f[0] != "5":
            continue
        s0, e0 = int(f[1]), int(f[2])
        if not (122450000 <= s0 < 122462000):
            continue
        clusters.append(dict(start0=s0, end0=e0, id=f[3], signal=float(f[4]),
                             strand=f[5], quant=float(f[6]), n=int(f[7]),
                             loc=f[9], motif=f[10]))

def clusters_near(pos0, tol=80):
    out = [c for c in clusters if c["strand"] == "-" and
           c["start0"] - tol <= pos0 <= c["end0"] + tol]
    return sorted(out, key=lambda c: -c["signal"])

rows = []

def add(site_id, cut0, sources, cls, status, note):
    sig = signal_scan(cut0)
    pct, run, risk, ds20 = internal_priming(cut0)
    near = clusters_near(cut0)
    best = near[0] if near else None
    rows.append(dict(
        site_id=site_id, chrom="chr5", cut0=cut0, cut1=cut0 + 1, strand="-",
        sources=sources, evidence_class=cls, status=status,
        pas_signal=(f"{sig[2]}@mRNA{sig[0]}({sig[1]})" if sig else "none"),
        nearest_polyasite=(f"{best['id']} sig={best['signal']} n={best['n']} "
                           f"q={best['quant']} loc={best['loc']}" if best else "none"),
        downstream_A20_pct=round(pct, 1), downstream_max_Arun=run,
        internal_priming_risk=risk, note=note + f" downstream20={ds20}"))

# --- 注释 3' 端（GENCODE vM25 exon 低坐标端 = 负链 3' 端；断裂点取注释端点） ---
add("annotated_3end_L_ENSMUST00000179939", 122453512,
    "GENCODE_vM25+PolyASite2.0", "公共数据直接观察",
    "annotated+external_support",
    "L/S 共用终端外显子 3' 端；PolyASite 强簇 122453501-534(sig332,n9,q0.94)+AATAAA@-32@122453540")
add("annotated_3end_S_ENSMUST00000177974", 122453513,
    "GENCODE_vM25+PolyASite2.0", "公共数据直接观察",
    "annotated+external_support", "同 L 簇支持")
add("annotated_3end_M_ENSMUST00000031423", 122456339,
    "GENCODE_vM25+PolyASite2.0", "公共数据直接观察",
    "annotated+external_support",
    "M 3' 端；PolyASite 簇 122456301-355(sig48,n7,q0.74)+AATAAA@-31@122456357")

# --- DaPars2 预测近端 PAS ---
add("predicted_PAS_DaPars2_122456498", 122456498,
    "DaPars2_event6953+PolyASite2.0", "计算推断",
    "predicted_window_boundary(no signal support)",
    "DaPars2 近端窗边界；PolyASite 仅有噪声簇(sig0.03/0.01, 无基序)；建议重新锚定到 122456339")

# --- PolyASite 其余可用簇（signal>=1 或 n>=5，负链、事件窗内） ---
seen_ends = {122453512, 122453513, 122456339, 122456498}
for c in clusters:
    if c["strand"] != "-":
        continue
    if c["signal"] < 1.0 and c["n"] < 5:
        continue
    peak = c["start0"] + (c["end0"] - c["start0"]) // 2
    if any(abs(peak - e) < 60 for e in seen_ends):
        continue
    add(f"polyasite_cluster_{c['id'].split(':')[1]}", peak,
        "PolyASite2.0", "公共数据直接观察", "cluster_only(unmapped_to_annotation)",
        f"window={c['start0']}-{c['end0']} sig={c['signal']} n={c['n']} q={c['quant']} loc={c['loc']} motif_field={c['motif']}")

cols = ["site_id", "chrom", "cut0", "cut1", "strand", "sources", "evidence_class",
        "status", "pas_signal", "nearest_polyasite", "downstream_A20_pct",
        "downstream_max_Arun", "internal_priming_risk", "note"]
with open(OUT, "w") as out:
    out.write("\t".join(cols) + "\n")
    for r in rows:
        out.write("\t".join(str(r[c]) for c in cols) + "\n")

print(f"written {OUT}  rows={len(rows)}")
for r in rows:
    print(f"  {r['site_id']}: cut={r['cut0']} sig={r['pas_signal']} "
          f"ip={r['internal_priming_risk']}(A{r['downstream_A20_pct']}%/run{r['downstream_max_Arun']}) "
          f"status={r['status']}")
