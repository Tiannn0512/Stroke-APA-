#!/usr/bin/env python3
"""P2.3b: motif summary + family enrichment from FIMO dual scans.

Frozen L3 rule (P4-03): hit with p < 1e-4 in the distal segment AND regulator family
in the frozen 8 (QKI, ELAVL1, PTBP1, MBNL2, NOVA1, NOVA2, TARDBP, HNRNPA2B1).
Family enrichment: per regulator, Fisher exact (candidates vs matched null, segments
with >=1 hit at p<=1e-4), BH over the 16 regulator families.

Outputs: motif_event_summary.tsv, motif_family_enrichment.tsv, motif_null_matching_qc.tsv
"""
import csv, math, os

BASE = r"D:\stroke_apa_reanalysis"
D = os.path.join(BASE, "results", "02_candidate_rebuild")
FROZEN8 = {"QKI", "ELAVL1", "PTBP1", "MBNL2", "NOVA1", "NOVA2", "TARDBP", "HNRNPA2B1"}
PTHRESH = 1e-4

def load_fimo(path):
    hits = {}
    with open(path, encoding="utf-8") as f:
        hdr = f.readline().rstrip("#\n").strip().split("\t")
        col = {c: i for i, c in enumerate(hdr)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(hdr):
                continue  # FIMO 5.5.9 trailing footer line
            motif = p[col["motif_id"]]
            seq = p[col["sequence_name"]]
            pv = float(p[col["p-value"]])
            reg = motif.split("_")[0].upper()
            hits.setdefault(seq, []).append((reg, motif, pv))
    return hits

cand_hits = load_fimo(os.path.join(D, "fimo_distal", "fimo.tsv"))
null_hits = load_fimo(os.path.join(D, "fimo_null", "fimo.tsv"))

cmap = {r["seq_id"]: r for r in csv.DictReader(open(os.path.join(D, "distal_candidates_mapping.tsv"), encoding="utf-8"), delimiter="\t")}

rows = []
for sid, m in cmap.items():
    hs = cand_hits.get(sid, [])
    fam_p = {}
    for reg, motif, pv in hs:
        if pv < PTHRESH:
            fam_p.setdefault(reg, []).append(pv)
    frozen = {k: min(v) for k, v in fam_p.items() if k in FROZEN8}
    rows.append({
        "seq_id": sid, "event_num": m["event_num"], "event_id": m["event_id"],
        "gene_symbol": m["gene_symbol"], "n_hits_any": len(hs),
        "n_hits_lt1e4": sum(1 for _, _, pv in hs if pv < PTHRESH),
        "L3v2_frozen8": "yes" if frozen else "no",
        "frozen8_families_hit": ";".join(sorted(frozen)),
        "frozen8_best_p": ";".join(f"{k}:{v:.2e}" for k, v in sorted(frozen.items())),
        "all_families_lt1e4": ";".join(sorted(fam_p)),
    })

fields = list(rows[0].keys())
with open(os.path.join(D, "motif_event_summary.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    w.writeheader(); w.writerows(rows)

# --- family enrichment: Fisher exact, BH ---
def segs_with_family(hits, sids, reg):
    n = 0
    for s in sids:
        if any(r == reg and pv <= PTHRESH for r, _, pv in hits.get(s, [])):
            n += 1
    return n

cand_sids = list(cmap)
null_sids = list(null_hits)  # seq ids are "<E####>_null"
null_to_cand = {r["seq_id"] + "_null": r["seq_id"] for r in csv.DictReader(open(os.path.join(D, "distal_candidates_mapping.tsv"), encoding="utf-8"), delimiter="\t")}

regs = sorted({r["Gene_name"].upper() for r in csv.DictReader(open(os.path.join(BASE, "input_links", "frozen_motif_family_mapping.tsv"), encoding="utf-8"), delimiter="\t")})

def fisher(a, b, c, dd):
    # two-sided hypergeometric
    n_ = a + b + c + dd
    def logcomb(n, k):
        return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
    def p(k):
        return math.exp(logcomb(a + c, k) + logcomb(b + dd, (a + b) - k) - logcomb(n_, a + b))
    lo = max(0, (a + b) - (b + dd))
    hi = min(a + b, a + c)
    pobs = p(a)
    s = 0.0
    for k in range(lo, hi + 1):
        if p(k) <= pobs * (1 + 1e-9):
            s += p(k)
    return min(1.0, s)

enr = []
for reg in regs:
    a = segs_with_family(cand_hits, cand_sids, reg)
    b = len(cand_sids) - a
    c = segs_with_family(null_hits, null_sids, reg)
    dd = len(null_sids) - c
    p = fisher(a, b, c, dd)
    enr.append({"family": reg, "frozen8": "yes" if reg in FROZEN8 else "no",
                "cand_seg_hits": a, "cand_seg_no": b, "null_seg_hits": c, "null_seg_no": dd,
                "fisher_p": p})

# BH over 16 families
m = len(enr)
prev = 1.0
for i, r in enumerate(sorted(enr, key=lambda x: x["fisher_p"], reverse=True), 1):
    prev = min(prev, r["fisher_p"] * m / (m - i + 1))
    r["fisher_q_bh"] = prev
enr.sort(key=lambda x: x["fisher_p"])

with open(os.path.join(D, "motif_family_enrichment.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(enr[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(enr)

# null matching QC
mm = list(csv.DictReader(open(os.path.join(D, "null_matched_mapping.tsv"), encoding="utf-8"), delimiter="\t"))
import statistics
qc = {
    "pairs": len(mm),
    "median_len_diff": statistics.median(int(x["len_diff"]) for x in mm),
    "max_len_diff": max(int(x["len_diff"]) for x in mm),
    "median_gc_diff_pp": round(statistics.median(float(x["gc_diff"]) for x in mm), 3),
    "max_gc_diff_pp": max(float(x["gc_diff"]) for x in mm),
    "p_threshold": PTHRESH, "frozen8": ";".join(sorted(FROZEN8)),
}
with open(os.path.join(D, "motif_null_matching_qc.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    for k, v in qc.items():
        w.writerow([k, v])

n_l3 = sum(1 for r in rows if r["L3v2_frozen8"] == "yes")
print(f"candidate segs: {len(rows)}; L3v2 (frozen8, p<1e-4): {n_l3}")
print("top families by enrichment p:")
for r in enr[:6]:
    print(f"  {r['family']:<12} cand {r['cand_seg_hits']}/{r['cand_seg_hits']+r['cand_seg_no']}  null {r['null_seg_hits']}/{r['null_seg_hits']+r['null_seg_no']}  p={r['fisher_p']:.3g} q={r['fisher_q_bh']:.3g} {'FROZEN8' if r['frozen8']=='yes' else ''}")
