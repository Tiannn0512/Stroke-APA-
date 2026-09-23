#!/usr/bin/env python3
"""P2.4: conservation per event from local mm10.60way.phastCons.bw (run in WSL bioinfo env).

For every event (all 9,693): distal + shared segments ->
  mean phastCons score over covered bases, bases with score>0, segment length,
  conserved_binary_v2 = (score>0 bases >= 10)  [frozen L4 rule, base-set identical to
  the v1 'merged element overlap >=10bp' because merging adjacent score>0 elements
  does not change the covered base set]
Coverage-none segments -> NA (marked).
Output: results/02_candidate_rebuild/distal_conservation.tsv
"""
import csv, os, sys
import pyBigWig

BASE = "/mnt/d/stroke_apa_reanalysis"
BW = os.path.join(BASE, "data", "references", "mm10.60way.phastCons.bw")
SEG = os.path.join(BASE, "results", "01_strand_audit", "strand_aware_event_segments.tsv")
OUT = os.path.join(BASE, "results", "02_candidate_rebuild", "distal_conservation.tsv")

bw = pyBigWig.open(BW)
chroms = bw.chroms()

def stats(chrom, s, e):
    s, e = int(s), int(e)
    if chrom not in chroms:
        return None
    iv = bw.intervals(chrom, s, e)
    if not iv:
        return {"mean": "", "pos_bases": 0, "covered_bases": 0}
    tot = pos = 0
    for b0, b1, v in iv:
        ov = min(e, b1) - max(s, b0)
        if ov <= 0:
            continue
        tot += ov
        if v > 0:
            pos += ov
    if tot == 0:
        return {"mean": "", "pos_bases": 0, "covered_bases": 0}
    mean = sum((b1v := (min(e, b1) - max(s, b0))) * v for b0, b1, v in iv if (b1v := min(e, b1) - max(s, b0)) > 0) / tot
    return {"mean": round(mean, 4), "pos_bases": pos, "covered_bases": tot}

rows_out = []
with open(SEG, encoding="utf-8") as f:
    for i, r in enumerate(csv.DictReader(f, delimiter="\t"), 1):
        d = stats(r["chrom"], r["distal_start"], r["distal_end"])
        s = stats(r["chrom"], r["shared_start"], r["shared_end"])
        dlen = int(r["distal_end"]) - int(r["distal_start"])
        rows_out.append({
            "event_num": r["event_num"], "event_id": r["event_id"], "gene_symbol": r["gene_symbol"],
            "chrom": r["chrom"], "strand": r["strand"],
            "distal_len": dlen,
            "distal_mean_phastcons": d["mean"] if d else "NA",
            "distal_pos_bases": d["pos_bases"] if d else "NA",
            "distal_covered_bases": d["covered_bases"] if d else "NA",
            "L4v2_conserved": ("yes" if d and d["pos_bases"] >= 10 else "no") if d else "NA",
            "shared_len": int(r["shared_end"]) - int(r["shared_start"]),
            "shared_mean_phastcons": s["mean"] if s else "NA",
            "shared_pos_bases": s["pos_bases"] if s else "NA",
        })
        if i % 1000 == 0:
            print(f"{i} events", flush=True)
bw.close()

with open(OUT, "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows_out[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(rows_out)

n_yes = sum(1 for r in rows_out if r["L4v2_conserved"] == "yes")
n_na = sum(1 for r in rows_out if r["L4v2_conserved"] == "NA")
print(f"written {len(rows_out)}; L4v2 conserved(>=10bp score>0): {n_yes}; NA(no coverage): {n_na}")
