#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G1 复核 · 工作 1：把观察到的剪接接头映射到 GENCODE vM25 的 Atp2a2 三条转录本结构。
输入：junction_counts_by_sample.tsv + gencode.vM25.annotation.gtf(.gz)
输出：junction_transcript_annotation.tsv
判定口径（对每条转录本独立判定）：
  exact_intron    接头 = 该转录本某相邻外显子对之间的内含子（供体+受体精确匹配）
  donor_only      仅供体端匹配某内含子起点
  acceptor_only   仅受体端匹配某内含子终点
  none            无内含子匹配
另标注接头是否落在该转录本 3'UTR 跨度内（three_prime_UTR 特征的最小起止）。
"""
import argparse
import gzip
from collections import defaultdict

TARGETS = ["ENSMUST00000177974.7", "ENSMUST00000179939.7", "ENSMUST00000031423.9"]


def parse_gtf(path, chrom="chr5", gene="Atp2a2"):
    exons = defaultdict(list)
    utrs = defaultdict(list)
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if p[0] != chrom or p[2] not in ("exon", "UTR"):
                continue
            attr = p[8]
            if f'gene_name "{gene}"' not in attr:
                continue
            tid = attr.split('transcript_id "')[1].split('"')[0]
            start0, end0 = int(p[3]) - 1, int(p[4])   # GTF 1-based closed → BED 0-based
            if p[2] == "exon":
                exons[tid].append((start0, end0))
            else:
                utrs[tid].append((start0, end0))
    introns = {}
    utr_relation = {}
    for tid, ex in exons.items():
        ex.sort()
        introns[tid] = [(ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)]
        # 负链基因：坐标最小的外显子是终端外显子；3' 端 UTR 块 = 低于第二低外显子上界的块
        u3 = sorted(b for b in utrs[tid] if b[1] <= ex[1][1]) if len(ex) > 1 else sorted(utrs[tid])
        utr_relation[tid] = u3
    return exons, introns, utr_relation


def utr_note(j, utr_blocks):
    """接头两端是否恰好连接该转录本的 3'UTR 块（负链：内含子左端=低坐标块末端，
    右端=高坐标块起点；donor/acceptor 为基因组左右侧坐标标签，非剪接方向）。"""
    s, e = j
    starts = {b[0] for b in utr_blocks}
    ends = {b[1] for b in utr_blocks}
    if s in ends and e in starts:
        return "connects_3UTR_block_to_3UTR_block(spliced_UTR)"
    if e in starts:
        return "right_end_abuts_3UTR_block"
    if s in ends:
        return "left_end_abuts_3UTR_block"
    return "-"


def classify(j, introns):
    s, e = j
    for k, (i0, i1) in enumerate(introns, 1):
        if (s, e) == (i0, i1):
            return f"exact_intron#{k}"
    for k, (i0, i1) in enumerate(introns, 1):
        if s == i0:
            return f"donor_only#intron{k}"
        if e == i1:
            return f"acceptor_only#intron{k}"
    return "none"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--junctions", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    exons, introns, utr_blocks = parse_gtf(a.gtf)
    for t in TARGETS:
        print(f"{t}: 外显子 {len(exons.get(t, []))} 个, "
              f"内含子 {len(introns.get(t, []))} 个, "
              f"3' 端 UTR 块 {utr_blocks.get(t, [])}")

    with open(a.out, "w") as out:
        out.write("junction_start0\tjunction_end0\tboth_ends_in_UTR\ttotal\t"
                  + "\t".join("vs_" + t for t in TARGETS) + "\t"
                  + "\t".join("utr_rel_" + t for t in TARGETS) + "\texon_counts\n")
        with open(a.junctions) as f:
            head = f.readline().rstrip("\n").split("\t")
            idx = {h: i for i, h in enumerate(head)}
            for line in f:
                p = line.rstrip("\n").split("\t")
                j = (int(p[idx["junction_start0"]]), int(p[idx["junction_end0"]]))
                ann = [classify(j, introns.get(t, [])) for t in TARGETS]
                rel = [utr_note(j, utr_blocks.get(t, [])) for t in TARGETS]
                exon_n = ";".join(f"{t}:{len(exons.get(t, []))}ex" for t in TARGETS)
                out.write(f"{j[0]}\t{j[1]}\t{p[idx['both_ends_in_UTR']]}\t{p[idx['total']]}\t"
                          + "\t".join(ann) + "\t" + "\t".join(rel) + "\t" + exon_n + "\n")
    print("written:", a.out)


if __name__ == "__main__":
    main()
