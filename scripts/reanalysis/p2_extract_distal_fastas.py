#!/usr/bin/env python3
"""P2.2: strand-aware distal-segment FASTA extraction (indexed random access via .fai).

Sets:
  candidates = unique events with >=1 significant shortening contrast
               (sig_joint TRUE & dPDUI <= -0.1); 1,541 (event,contrast) pairs -> 824 unique events
               (v1 rule recovered exactly: old lost_segments.bed also had 1,541 rows)
  null pool  = distal segments of events with NO significant contrast (9,693 - 977)

Sequences are reverse-complemented for minus-strand segments (bedtools getfasta -s semantics).
Headers use simple ids (E####/N####) to avoid the FIMO '::' UCSC-parsing mangling seen in v1;
mapping tables are written alongside.

Outputs (results/02_candidate_rebuild/):
  distal_candidates.strand_aware.fa (+ _mapping.tsv with len/GC)
  distal_nullpool.fa (+ _mapping.tsv with len/GC, event_id)
"""
import csv, os

BASE = r"D:\stroke_apa_reanalysis"
FA = os.path.join(BASE, "data", "references", "genome.fa")
FAI = FA + ".fai"
SEG = os.path.join(BASE, "results", "01_strand_audit", "strand_aware_event_segments.tsv")
SAP = os.path.join(BASE, "input_links", "p2_sap_layer2_events.tsv")
OUT = os.path.join(BASE, "results", "02_candidate_rebuild")

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def load_fai():
    idx = {}
    with open(FAI, encoding="utf-8") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            idx[p[0]] = (int(p[1]), int(p[2]), int(p[3]), int(p[4]))
    return idx

def fetch(f, fai, chrom, start, end):
    ln, off, lb, lw = fai[chrom]
    seq = []
    pos, need = off + (start // lb) * lw, end - start
    f.seek(pos)
    while need > 0:
        line = f.readline().rstrip("\n")
        take = min(len(line), need)
        seq.append(line[:take])
        need -= take
    return "".join(seq)

def gc(seq):
    s = seq.upper()
    gc_n = s.count("G") + s.count("C") + s.count("S")
    acgt = sum(s.count(x) for x in "ACGT")
    return round(100.0 * gc_n / acgt, 2) if acgt else 0.0

fai = load_fai()

# event sets from SAP
sig_events, short_events = set(), set()
with open(SAP, encoding="utf-8") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["sig_joint"] == "TRUE":
            sig_events.add(r["event"])
            if float(r["dPDUI"]) <= -0.1:
                short_events.add(r["event"])

segs = {}  # keyed by event_num (string) — the integer `event` used by SAP tables
with open(SEG, encoding="utf-8") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        segs[r["event_num"]] = r

def extract(events, prefix, fa_name, map_name):
    meta, out = [], []
    with open(os.path.join(OUT, fa_name), "w", encoding="utf-8", newline="") as fa:
        for i, eid in enumerate(sorted(events), 1):
            r = segs[eid]
            sid = f"{prefix}{i:04d}"
            with open(FA, encoding="utf-8", errors="ignore") as g:
                seq = fetch(g, fai, r["chrom"], int(r["distal_start"]), int(r["distal_end"]))
            if r["strand"] == "-":
                seq = seq.translate(COMP)[::-1]
            fa.write(f">{sid}\n")
            for j in range(0, len(seq), 60):
                fa.write(seq[j:j+60] + "\n")
            meta.append({"seq_id": sid, "event_num": eid, "event_id": r["event_id"], "chrom": r["chrom"],
                         "start": r["distal_start"], "end": r["distal_end"],
                         "strand": r["strand"], "gene_symbol": r["gene_symbol"],
                         "len": len(seq), "gc_pct": gc(seq)})
    with open(os.path.join(OUT, map_name), "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(meta[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(meta)
    return meta

cand_events = sorted((e for e in short_events if e in segs), key=int)
nonsig_events = sorted((e for e in segs if e not in sig_events), key=int)
mc = extract(cand_events, "E", "distal_candidates.strand_aware.fa", "distal_candidates_mapping.tsv")
mn = extract(nonsig_events, "N", "distal_nullpool.fa", "distal_nullpool_mapping.tsv")
print(f"candidate segments: {len(mc)} (from {len(short_events)} shortening events; v1 pairs = 1541)")
print(f"null pool segments : {len(mn)} (events with no significant contrast)")
