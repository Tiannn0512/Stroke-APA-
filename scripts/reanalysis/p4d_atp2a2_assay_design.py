#!/usr/bin/env python3
"""Task D: Atp2a2 sequence-level assay design (3'RACE GSPs, common/distal qPCR pairs).

All coordinates mm10 BED 0-based half-open; gene on MINUS strand -> transcript-
direction sequence = reverse complement of the genomic slice. Specificity proxy:
exact-match screen of each primer (both orientations) against GENCODE vM25
transcriptome; web Primer-BLAST check remains pending with the experiment lead.

Outputs: results/05_target_nomination/Atp2a2_assay_design.tsv,
         Atp2a2_specificity_screen.txt
"""
import csv, os
import primer3

BASE = r"D:\stroke_apa_reanalysis"
FA = os.path.join(BASE, "data", "references", "genome.fa")
FAI = FA + ".fai"
TX = os.path.join(BASE, "data", "references", "gencode.vM25.transcripts.fa")
OUT = os.path.join(BASE, "results", "05_target_nomination")
os.makedirs(OUT, exist_ok=True)
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

fai = {}
with open(FAI, encoding="utf-8") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        fai[p[0]] = (int(p[1]), int(p[2]), int(p[3]), int(p[4]))

def fetch(chrom, start, end, strand="-"):
    ln, off, lb, lw = fai[chrom]
    seq, need, pos = [], end - start, off + (start // lb) * lw
    with open(FA, encoding="utf-8", errors="ignore") as f:
        f.seek(pos)
        while need > 0:
            line = f.readline().rstrip("\n")
            take = min(len(line), need)
            seq.append(line[:take]); need -= take
    s = "".join(seq).upper()
    return s.translate(COMP)[::-1] if strand == "-" else s

P3_COMMON = {
    "PRIMER_TASK": "generic",
    "PRIMER_PICK_LEFT_PRIMER": 1, "PRIMER_PICK_RIGHT_PRIMER": 1,
    "PRIMER_MIN_SIZE": 20, "PRIMER_OPT_SIZE": 22, "PRIMER_MAX_SIZE": 26,
    "PRIMER_MIN_TM": 58.0, "PRIMER_OPT_TM": 60.0, "PRIMER_MAX_TM": 62.0,
    "PRIMER_MIN_GC": 40.0, "PRIMER_MAX_GC": 60.0,
    "PRIMER_PRODUCT_SIZE_RANGE": [[80, 150]],
    "PRIMER_MAX_POLY_X": 4,
}
P3_GSP = {
    "PRIMER_TASK": "pick_primer_list", "PRIMER_PICK_LEFT_PRIMER": 1,
    "PRIMER_PICK_RIGHT_PRIMER": 0, "PRIMER_PICK_INTERNAL_OLIGO": 0,
    "PRIMER_MIN_SIZE": 22, "PRIMER_OPT_SIZE": 24, "PRIMER_MAX_SIZE": 28,
    "PRIMER_MIN_TM": 60.0, "PRIMER_OPT_TM": 62.0, "PRIMER_MAX_TM": 65.0,
    "PRIMER_MIN_GC": 40.0, "PRIMER_MAX_GC": 65.0, "PRIMER_MAX_POLY_X": 4,
}

def minus_seq(start, end):
    return fetch("chr5", start, end, strand="-")

print("loading transcriptome for specificity screen ...")
txseqs = []
sid, buf = None, []
with open(TX, encoding="utf-8", errors="ignore") as f:
    for line in f:
        if line.startswith(">"):
            if sid: txseqs.append((sid, "".join(buf)))
            sid = line[1:].split("|")[0].split(" ")[0]; buf = []
        else:
            buf.append(line.strip().upper())
    if sid: txseqs.append((sid, "".join(buf)))
print(f"transcripts loaded: {len(txseqs)}")

def screen(primer_seq):
    rc = primer_seq.translate(COMP)[::-1]
    hits = []
    for tid, s in txseqs:
        if primer_seq in s: hits.append((tid, "+"))
        elif rc in s: hits.append((tid, "-"))
    return hits

def design_gsp(zone, label, n=2):
    seq = minus_seq(*zone)  # transcript direction
    res = primer3.bindings.design_primers({"SEQUENCE_ID": label, "SEQUENCE_TEMPLATE": seq}, P3_GSP)
    out = []
    for i in range(n):
        key = f"PRIMER_LEFT_{i}_SEQUENCE"
        if key not in res: break
        st = res[f"PRIMER_LEFT_{i}"]
        out.append({"label": f"{label}_{i+1}",
                    "seq_transcript_dir": res[key],
                    "pos_in_zone": st[0],
                    "tm": round(res[f"PRIMER_LEFT_{i}_TM"], 1),
                    "gc": round(res[f"PRIMER_LEFT_{i}_GC_PERCENT"], 1)})
    return out

GSP_OUTER = (122456598, 122457098)
GSP_NESTED = (122456518, 122456698)
COMMON = (122456558, 122457103)
DISTAL = (122454500, 122456300)  # 179939 独有区段：避开 ENSMUST00000177974.7 的重叠 UTR(122453513-122454293)

results = []
for zone, label in ((GSP_OUTER, "3RACE_outer_GSP"), (GSP_NESTED, "3RACE_nested_GSP")):
    for g in design_gsp(zone, label):
        hits = screen(g["seq_transcript_dir"])
        off = len([h for h in hits if not h[0].startswith("ENSMUST00000179939")])
        g.update(assay=label, zone=f"{zone[0]}-{zone[1]}", exact_hits_offatp2a2=off)
        results.append(g)

for zname, zone in (("qPCR_common", COMMON), ("qPCR_distal", DISTAL)):
    zs = minus_seq(*zone)
    res_z = primer3.bindings.design_primers({"SEQUENCE_ID": zname, "SEQUENCE_TEMPLATE": zs}, P3_COMMON)
    for i in range(min(1, res_z["PRIMER_PAIR_NUM_RETURNED"])):
        lf, lr = res_z[f"PRIMER_LEFT_{i}"], res_z[f"PRIMER_RIGHT_{i}"]
        prod = lr[0] + lr[1] - lf[0]
        lseq, rseq = res_z[f"PRIMER_LEFT_{i}_SEQUENCE"], res_z[f"PRIMER_RIGHT_{i}_SEQUENCE"]
        offids = set()
        for half in (lseq, rseq):
            for tid, st in screen(half):
                if not tid.startswith("ENSMUST00000179939"): offids.add(tid)
        results.append(dict(assay=zname, seq_transcript_dir=f"{lseq} / {rseq}",
                            tm=f"{res_z[f'PRIMER_LEFT_{i}_TM']:.1f}/{res_z[f'PRIMER_RIGHT_{i}_TM']:.1f}",
                            gc=f"{res_z[f'PRIMER_LEFT_{i}_GC_PERCENT']:.0f}/{res_z[f'PRIMER_RIGHT_{i}_GC_PERCENT']:.0f}",
                            product_bp=prod, zone=f"{zone[0]}-{zone[1]}",
                            exact_hits_offatp2a2=len(offids), offatp2a2_ids=";".join(sorted(offids))))

fields = ["assay", "seq_transcript_dir", "tm", "gc", "product_bp", "pos_in_zone", "zone", "exact_hits_offatp2a2", "offatp2a2_ids"]
with open(os.path.join(OUT, "Atp2a2_assay_design.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    w.writeheader()
    for r in results:
        w.writerow({k: r.get(k, "") for k in fields})
print("designs:", len(results))
for r in results:
    print(" ", r["assay"], "|", str(r.get("seq_transcript_dir", ""))[:60], "| Tm", r.get("tm"), "| off-hits", r.get("exact_hits_offatp2a2", ""))
