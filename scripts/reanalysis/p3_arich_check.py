#!/usr/bin/env python3
"""P3 A-rich internal-priming risk check around predicted PAS (+-50nt genomic window).

Flags per event: >=6 consecutive genomic A in window, or window A% >= 60.
Uses indexed genome.fa (.fai). Output: results/03_bam_review/arich_risk.tsv
"""
import csv, os

BASE = r"D:\stroke_apa_reanalysis"
FA = os.path.join(BASE, "data", "references", "genome.fa")
FAI = FA + ".fai"
BR = os.path.join(BASE, "results", "03_bam_review")
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

fai = {}
with open(FAI, encoding="utf-8") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        fai[p[0]] = (int(p[1]), int(p[2]), int(p[3]), int(p[4]))

def fetch(chrom, start, end):
    ln, off, lb, lw = fai[chrom]
    out, need, pos = [], end - start, off + (start // lb) * lw
    with open(FA, encoding="utf-8", errors="ignore") as f:
        f.seek(pos)
        while need > 0:
            line = f.readline().rstrip("\n")
            take = min(len(line), need)
            out.append(line[:take]); need -= take
    return "".join(out)

rows = []
for line in open(os.path.join(BR, "regions.bed"), encoding="utf-8"):
    p = line.rstrip("\n").split("\t")
    gene = p[3].split("|")[0]
    chrom, us, ue, strand = p[0], int(p[1]) + 5000, int(p[2]) - 5000, p[5]
    # predicted PAS from ranking table
    r = next(x for x in csv.DictReader(open(os.path.join(BASE, "results", "02_candidate_rebuild", "candidate_ranking_v2.tsv"), encoding="utf-8"), delimiter="\t") if x["gene_symbol"] == gene)
    pas = int(r["proximal_pas_bed"])
    w0, w1 = max(0, pas - 50), pas + 50
    gseq = fetch(chrom, w0, w1).upper()
    if strand == "-":
        gseq_t = gseq.translate(COMP)[::-1]  # transcript-direction window
    else:
        gseq_t = gseq
    max_run, run = 0, 0
    for ch in gseq_t:
        run = run + 1 if ch == "A" else 0
        max_run = max(max_run, run)
    apct = round(100.0 * gseq_t.count("A") / max(1, len(gseq_t)), 1)
    rows.append({"gene_symbol": gene, "event_num": r["event_num"], "strand": strand,
                 "pas_pas_bed": pas, "window": f"{w0}-{w1}",
                 "max_A_run": max_run, "A_pct": apct,
                 "arich_risk": "HIGH" if (max_run >= 6 or apct >= 60) else "low"})

with open(os.path.join(BR, "arich_risk.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(rows)
for r in rows:
    print(f"{r['gene_symbol']:<9} {r['strand']} maxArun={r['max_A_run']} A%={r['A_pct']} -> {r['arich_risk']}")
