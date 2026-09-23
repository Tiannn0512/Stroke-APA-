"""Rebuild DaPars2 tandem 3'UTR annotation from GENCODE GTF (rebuild after pilot deletion).
DaPars2-official semantics (mirrors DaPars_Extract_Anno bed logic):
  + strand: UTR = (max CDS end, max exon end]
  - strand: UTR = [min exon start, min CDS start)
per transcript, only protein_coding CDS transcripts. BED 0-based half-open:
  chr  start  end  transcript_id|gene_symbol  0  strand
"""
from collections import defaultdict

gtf = "/home/taylor/reference_v2/gencode.vM25.annotation.gtf"
out = "/home/taylor/reference_v2/gencode_M25_3UTR_for_DaPars2.bed"

ex = defaultdict(lambda: [None, None])   # tid -> [min_start, max_end]
cds = defaultdict(lambda: [None, None])  # tid -> [min_start, max_end]
info = {}                                # tid -> (chr, strand, gene_symbol)
for line in open(gtf):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    feat = f[2]
    if feat not in ("exon", "CDS"):
        continue
    attrs = f[8]
    tid = attrs.split('transcript_id "')[1].split('"')[0]
    s, e = int(f[3]), int(f[4])
    d = ex if feat == "exon" else cds
    rec = d[tid]
    rec[0] = s if rec[0] is None else min(rec[0], s)
    rec[1] = e if rec[1] is None else max(rec[1], e)
    if tid not in info:
        gname = attrs.split('gene_name "')[1].split('"')[0]
        gid = attrs.split('gene_id "')[1].split('"')[0].split(".")[0]
        info[tid] = (f[0], f[6], f"{gid}|{gname}")

n = 0
with open(out, "w") as fo:
    for tid, (c, st, name) in info.items():
        c_rec, d_rec = ex.get(tid), cds.get(tid)
        if not c_rec or not d_rec:
            continue
        if st == "+":
            u_start, u_end = d_rec[1], c_rec[1]      # CDS end -> exon end (0-based: [cds_end, exon_end))
        else:
            u_start, u_end = c_rec[0] - 1, d_rec[0] - 1  # exon start -> CDS start
        if u_end - u_start < 10:                      # skip <10 nt trivial UTRs
            continue
        fo.write(f"{c}\t{u_start}\t{u_end}\t{tid}|{name}\t0\t{st}\n")
        n += 1
print(f"3'UTR regions written: {n} transcripts -> {out}")
