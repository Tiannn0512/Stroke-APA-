#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
f1 (v2): GENCODE vM25 Atp2a2 全转录本结构 → 结构 TSV + facts JSON。

关键输出：
  1) 逐转录本 exon / 3p_UTR / 剪接接头(intron) / 注释3'端
  2) 122453000-122460000 逐碱基覆盖唯一性分段（哪些序列 179939 独有 / 177974+179939 共有 / 全共有…）
  3) 每个转录本的 基因组<->mRNA(cDNA,1-based) 坐标映射（供引物脚本用）
坐标统一 BED 0-based half-open。
先验位点（DaPars2 预测 PAS、QKI peak）单独行并标注状态。
"""
import gzip
import json

GTF = "/mnt/d/stroke_apa_reaudit/reference/gencode.vM25.annotation.gtf.gz"
OUTDIR = "/mnt/d/stroke_apa_reaudit/results/07_followup_20260925"
OUT_TSV = OUTDIR + "/Atp2a2_transcript_structure.tsv"
OUT_JSON = OUTDIR + "/analysis_scripts/structure_facts.json"
GENE = "Atp2a2"
MAIN3 = ["ENSMUST00000179939.7", "ENSMUST00000177974.7", "ENSMUST00000031423.9"]
SHORT = {"ENSMUST00000179939.7": "L(long-UTR,事件载体)",
         "ENSMUST00000177974.7": "S(short-UTR)",
         "ENSMUST00000031423.9": "M(medium-UTR)"}

def parse_attr(s):
    d = {}
    for kv in s.strip().rstrip(";").split("; "):
        if " " in kv:
            k, v = kv.split(" ", 1)
            d[k] = v.strip('"')
    return d

tx = {}
with gzip.open(GTF, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        chrom, feat, s1, e1, strand, attr = f[0], f[2], int(f[3]), int(f[4]), f[6], f[8]
        if chrom != "chr5" or not (122000000 <= s1 <= 123500000):
            continue
        d = parse_attr(attr)
        if d.get("gene_name") != GENE or "transcript_id" not in d:
            continue
        tid = d["transcript_id"]
        s0, e0 = s1 - 1, e1
        info = tx.setdefault(tid, dict(exons=[], utr_all=[], cds=[], strand=strand,
                                       name=d.get("transcript_name", "")))
        if feat == "exon":
            info["exons"].append((s0, e0))
        elif feat == "UTR":
            info["utr_all"].append((s0, e0))
        elif feat == "CDS":
            info["cds"].append((s0, e0))

for t in tx.values():
    t["exons"].sort()
    t["utr_all"] = t.pop("utr3", []) if "utr3" in t else t.get("utr_all", [])
    t.setdefault("utr_all", [])
    t.setdefault("cds", [])
    t["utr_all"].sort()
    t["cds"].sort()
    # GENCODE 的 UTR 不分 3'/5'：负链按“低于 CDS 最低坐标”划为 3'UTR
    cds_min = min((s for s, _ in t["cds"]), default=None)
    if cds_min is not None:
        t["utr3"] = [(s, e) for s, e in t["utr_all"] if e <= cds_min]
    else:
        t["utr3"] = list(t["utr_all"])

rows = []
for tid in sorted(tx):
    info = tx[tid]
    for i, (s0, e0) in enumerate(info["exons"]):
        rows.append((tid, "exon", s0, e0, f"exon {i+1}/{len(info['exons'])} genomic-asc"))
    for s0, e0 in info["utr3"]:
        rows.append((tid, "3p_UTR", s0, e0, "GENCODE vM25 three_prime_UTR"))
    # introns = junctions（mRNA 顺序）
    exs_hi2lo = sorted(info["exons"], reverse=True)
    for i in range(len(exs_hi2lo) - 1):
        hi2 = exs_hi2lo[i + 1][1]      # 下游（低坐标）exon 的 hi（0-based 不含）
        lo_up = exs_hi2lo[i][0]        # 上游（高坐标）exon 的 lo
        rows.append((tid, "splice_junction_intron", hi2, lo_up,
                     f"intron#{i+1}from3p donor_base={lo_up} acceptor_base={hi2-1} (BED)"))

# 先验位点
rows.append(("ENSMUST00000179939.7", "predicted_PAS_DaPars2", 122456498, 122456499,
             "predicted近端PAS; source=DaPars2 event6953(发现队列); 类别=计算推断; 状态=predicted(未测RNA-polyA接合)"))
rows.append(("ENSMUST00000179939.7", "QKI_CLIP_peak", 122454200, 122454250,
             "FDR=1.04e-5; source=发现队列eCLIP; 类别=计算推断"))

with open(OUT_TSV, "w") as out:
    out.write("transcript_id\tversion_label\tfeature\tchrom\tstart0\tend0\tstrand\tnote\n")
    for tid, feat, s0, e0, note in rows:
        out.write(f"{tid}\t{SHORT.get(tid,'')}\t{feat}\tchr5\t{s0}\t{e0}\t{tx[tid]['strand']}\t{note}\n")

# ---- facts ----
facts = dict(transcripts={}, coverage_segments=[], junctions={})
for tid in sorted(tx):
    info = tx[tid]
    exs = sorted(info["exons"], reverse=True)      # mRNA 顺序
    cmap = []                                       # (g_lo, g_hi, cdna_start1) mRNA 顺序
    pos = 1
    for (lo, hi) in exs:
        L = hi - lo
        cmap.append((lo, hi, pos))
        pos += L
    exs_asc = sorted(info["exons"])
    facts["transcripts"][tid] = dict(
        name=info["name"], strand=info["strand"], label=SHORT.get(tid, "minor"),
        exons_asc=exs_asc, exons_mrna_order=exs_asc[::-1] if info["strand"] == "-" else exs_asc,
        cdna_len=sum(h - l for l, h in exs_asc),
        utr3_blocks=info["utr3"],
        utr3_total=sum(e - s for s, e in info["utr3"]),
        annotated_3end_genomic=min(l for l, _ in exs_asc) if info["strand"] == "-" else None,
        g2c_blocks=cmap)

# 逐碱基覆盖唯一性分段（主3，窗口 122453000-122460000）
LO, HI = 122453000, 122460000
cov = {t: set() for t in MAIN3}
for t in MAIN3:
    for lo, hi in tx[t]["exons"]:
        a, b = max(lo, LO), min(hi, HI)
        cov[t].update(range(a, b))
segs, cur, curmem = [], LO, None
for p in range(LO, HI):
    mem = tuple(t for t in MAIN3 if p in cov[t])
    if mem != curmem:
        if curmem is not None:
            segs.append((cur, p, curmem))
        cur, curmem = p, mem
segs.append((cur, HI, curmem))
# 合并显示：只列 ≥1bp 且 membership 非空的段
merged = [s for s in segs if s[2]]
with open(OUT_TSV, "a") as out:
    out.write("#coverage_uniqueness_segments\twindow=chr5:122453000-122460000(BED)\t"
              "membership: L=179939 S=177974 M=31423\n")
    for s0, e0, mem in merged:
        lab = "".join({"ENSMUST00000179939.7": "L", "ENSMUST00000177974.7": "S",
                       "ENSMUST00000031423.9": "M"}[m] for m in mem)
        out.write(f"COVERAGE_SEG\t{lab}\t\tchr5\t{s0}\t{e0}\t-\t{len(mem)}转录本覆盖\n")
facts["coverage_segments"] = [list(s) for s in merged]

# 接头专属性（主3 各自的最后一两个 intron + 唯一性判定）
def introns_of(tid):
    exs = sorted(tx[tid]["exons"], reverse=True)
    return [(exs[i+1][1], exs[i][0]) for i in range(len(exs)-1)]

# 接头专属性：精确 (donor,acceptor) 对比——负链 mRNA 顺序里
# junction = (下游exon的hi, 上游exon的lo)，两个转录本只有坐标对完全相同才算同一接头
all_intr = {t: set(introns_of(t)) for t in MAIN3}
for tid in MAIN3:
    uniq = [intr for intr in all_intr[tid]
            if not any(intr in all_intr[t2] for t2 in MAIN3 if t2 != tid)]
    ins = sorted(all_intr[tid])
    facts["junctions"][tid] = dict(last_introns=ins[-2:], unique_introns=uniq)

with open(OUT_JSON, "w") as out:
    json.dump(facts, out, indent=1)

print("== transcripts ==")
for tid, f in facts["transcripts"].items():
    print(f"{tid} {f['name']} {f['label']} cdna={f['cdna_len']} utr3={f['utr3_total']}bp "
          f"3end={f['annotated_3end_genomic']}")
print("== coverage segments (L/S/M) ==")
for s0, e0, mem in merged:
    lab = "".join({"ENSMUST00000179939.7": "L", "ENSMUST00000177974.7": "S",
                   "ENSMUST00000031423.9": "M"}[m] for m in mem)
    print(f"  {s0}-{e0} ({e0-s0}bp) [{lab}]")
print("== unique junctions ==")
for tid, j in facts["junctions"].items():
    print(f"  {tid}: unique={j['unique_introns']}")
