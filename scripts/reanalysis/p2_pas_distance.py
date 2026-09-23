#!/usr/bin/env python3
"""P2.5 addendum: PolyASite 2.0 (GRCm38.96) nearest-PAS distance for the 977 significant events.

Distance = |predicted PAS (BED scale) - nearest same-strand atlas cluster position|,
computed via sorted per-chrom/per-strand arrays + binary search.
Also flags: <50bp = strong support; <200 = moderate; <500 = weak; else none.
Outputs pas_distance_977.tsv and appends pas_dist columns to
event_annotation_matrix.v2.tsv / candidate_ranking_v2.tsv in place.
"""
import bisect, csv, gzip, os

BASE = r"D:\stroke_apa_reanalysis"
PAS = os.path.join(BASE, "data", "references", "PAS", "atlas.clusters.2.0.GRCm38.96.bed.gz")
SEG = os.path.join(BASE, "results", "01_strand_audit", "strand_aware_event_segments.tsv")
R = os.path.join(BASE, "results", "02_candidate_rebuild")

idx = {}  # (chrom, strand) -> sorted list of positions
with gzip.open(PAS, "rt", encoding="utf-8") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < 6:
            continue
        chrom = p[0] if p[0].startswith("chr") else "chr" + p[0]
        idx.setdefault((chrom, p[5]), []).append(int(p[1]))
for k in idx:
    idx[k].sort()
print(f"atlas clusters loaded: {sum(len(v) for v in idx.values())}")

def nearest(chrom, strand, pos):
    arr = idx.get((chrom, strand))
    if not arr:
        return None
    i = bisect.bisect_left(arr, pos)
    best = None
    for j in (i - 1, i):
        if 0 <= j < len(arr):
            d = abs(arr[j] - pos)
            best = d if best is None else min(best, d)
    return best

out = []
with open(SEG, encoding="utf-8") as f:
    sap_sig = {r["event"] for r in csv.DictReader(
        open(os.path.join(BASE, "input_links", "p2_sap_layer2_events.tsv"), encoding="utf-8"), delimiter="\t")
        if r["sig_joint"] == "TRUE"}
    for r in csv.DictReader(f, delimiter="\t"):
        if r["event_num"] not in sap_sig:
            continue
        d = nearest(r["chrom"], r["strand"], int(r["proximal_pas_bed"]))
        flag = "none" if d is None else ("<50" if d < 50 else "<200" if d < 200 else "<500" if d < 500 else ">=500")
        out.append({"event_num": r["event_num"], "gene_symbol": r["gene_symbol"],
                    "pas_dist_bp": "" if d is None else d,
                    "pas_support": "NA" if d is None else flag})

with open(os.path.join(R, "pas_distance_977.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(out)

# append columns to matrix v2 and ranking v2
pd_map = {r["event_num"]: r for r in out}
for name in ("event_annotation_matrix.v2.tsv", "candidate_ranking_v2.tsv"):
    p = os.path.join(R, name)
    rows = list(csv.DictReader(open(p, encoding="utf-8"), delimiter="\t"))
    for r in rows:
        d = pd_map.get(r["event_num"], {})
        r["pas_dist_bp"] = d.get("pas_dist_bp", "")
        r["pas_support"] = d.get("pas_support", "NA")
    with open(p, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(rows)

import collections
c = collections.Counter(r["pas_support"] for r in out)
print("support buckets:", dict(c))
print("written:", len(out))
