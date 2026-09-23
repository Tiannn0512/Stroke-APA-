#!/usr/bin/env python3
"""Second-method (salmon TPM) long-isoform baseline for Gate-1 targets.

Atp2a2 and Agpat3 have multiple annotated 3'UTR transcripts -> transcript-level
long fraction computable. Aplp1 has a single annotated transcript in the tandem
BED -> not computable at transcript level (stays 'absent', honestly logged).

long_frac = TPM(long-3'UTR transcript = the DaPars2 event transcript)
            / sum(TPM(all gene transcripts in the DaPars2 annotation))
Output: results/02_candidate_rebuild/salmon_gate1_baseline.tsv
"""
import csv, os

BASE = r"D:\stroke_apa_reanalysis"
BED = os.path.join(BASE, "input_links", "gencode_M25_3UTR_for_DaPars2.bed")
S = ["sham1", "sham2", "day3_rep1", "day3_rep2", "day7_rep1", "day7_rep2",
     "day21_rep1", "day21_rep2", "day60_rep1", "day60_rep2", "day1_rep1", "day1_rep2"]

TARGETS = {
    "Atp2a2": ("ENSMUSG00000029467", "ENSMUST00000179939.7"),
    "Agpat3": ("ENSMUSG00000001211", "ENSMUST00000001240.11"),
}

gene_txs = {}
with open(BED, encoding="utf-8") as f:
    for line in f:
        p = line.split("\t")
        if len(p) > 3 and "|" in p[3]:
            tx, gid = p[3].split("|")[0], p[3].split("|")[1]
            gene_txs.setdefault(gid, set()).add(tx)

tpm = {s: {} for s in S}
for s in S:
    with open(os.path.join(BASE, "data", "salmon", s, "quant.sf"), encoding="utf-8") as f:
        hdr = f.readline().split("\t")
        ti, pi = hdr.index("Name"), hdr.index("TPM")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(hdr):
                continue
            tpm[s][p[ti]] = float(p[pi])

rows = []
for sym, (gid, longtx) in TARGETS.items():
    txs = gene_txs.get(gid, set())
    for s in S:
        tot = sum(tpm[s].get(t, 0.0) for t in txs)
        lon = tpm[s].get(longtx, 0.0)
        rows.append({"gene_symbol": sym, "sample": s,
                     "long_tx_tpm": round(lon, 2), "gene_total_tpm": round(tot, 2),
                     "long_frac_salmon": round(lon / tot, 3) if tot > 0 else "NA"})

with open(os.path.join(BASE, "results", "02_candidate_rebuild", "salmon_gate1_baseline.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(rows)

order = ["sham1", "sham2", "day3_rep1", "day3_rep2", "day7_rep1", "day7_rep2", "day21_rep1", "day21_rep2", "day60_rep1", "day60_rep2"]
by = {(r["gene_symbol"], r["sample"]): r for r in rows}
for sym in TARGETS:
    print(f"-- {sym} long_frac --")
    for s in order:
        print(f"  {s:<12} {by[(sym, s)]['long_frac_salmon']:>7}")
