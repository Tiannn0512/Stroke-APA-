#!/usr/bin/env python3
"""P2.1 per-event QKI overlap summary + TODO completion-standard verification.

Reads bedtools -wa -wb intersections (distal & shared x QKI peaks).
Peak columns: chrom start end gene FDR strand; segment columns: chrom s e name score strand.
Writes qki_per_event_summary.tsv for all 9,693 events:
  event_id, gene_symbol, chrom, strand, distal_peak_count, distal_min_fdr,
  distal_overlap_bp, shared_peak_count, plus L2v2 flag (distal peak with FDR<0.05).
Verifies the three P2.1 completion standards and prints PASS/FAIL.
"""
import csv, os

BASE = r"D:\stroke_apa_reanalysis"
D = os.path.join(BASE, "results", "02_candidate_rebuild")

def load_intersections(path, side):
    hits = {}
    if not os.path.exists(path):
        return hits
    with open(path, encoding="utf-8") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            seg_chrom, s, e, name = p[0], int(p[1]), int(p[2]), p[3]
            peak_chrom, ps, pe, gene, fdr = p[6], int(p[7]), int(p[8]), p[9], float(p[10])
            ov = min(e, pe) - max(s, ps)
            hits.setdefault(name, []).append({
                "side": side, "peak": f"{peak_chrom}:{ps}-{pe}", "gene": gene,
                "fdr": fdr, "overlap_bp": max(0, ov),
            })
    return hits

distal = load_intersections(os.path.join(D, "qki_distal_intersections.tsv"), "distal")
shared = load_intersections(os.path.join(D, "qki_shared_intersections.tsv"), "shared")

events = {}
with open(os.path.join(BASE, "input_links", "dapars2_event_coordinates.tsv"), encoding="utf-8") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        events[r["event_id"]] = r

out_rows = []
for eid, ev in events.items():
    dh = distal.get(eid, [])
    sh = shared.get(eid, [])
    d_pass = [h for h in dh if h["fdr"] < 0.05]
    out_rows.append({
        "event_id": eid, "gene_symbol": ev["gene_symbol"], "chrom": ev["chrom"], "strand": ev["strand"],
        "distal_peak_count": len(dh), "distal_peak_count_fdr05": len(d_pass),
        "distal_min_fdr": min((h["fdr"] for h in dh), default=""),
        "distal_overlap_bp": sum(h["overlap_bp"] for h in dh),
        "distal_peaks": ";".join(f"{h['peak']}(fdr={h['fdr']:.2e})" for h in dh[:6]),
        "shared_peak_count": len(sh),
        "shared_peaks": ";".join(f"{h['peak']}(fdr={h['fdr']:.2e})" for h in sh[:6]),
        "L2v2_distal_qki": "yes" if d_pass else ("yes_rawfdr" if dh else "no"),
    })

fields = list(out_rows[0].keys())
outp = os.path.join(D, "qki_per_event_summary.tsv")
with open(outp, "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    w.writeheader(); w.writerows(out_rows)

# --- completion standards ---
print("== P2.1 completion standards ==")
ok1 = ok2 = ok3 = False
by_sym = {}
for r in out_rows:
    by_sym.setdefault(r["gene_symbol"], []).append(r)

a = by_sym.get("Atp2a2", [])[0]
d_hits = distal.get(a["event_id"], [])
ok1 = any(h["peak"] == "chr5:122454200-122454250" for h in d_hits)
print(f"  [1] Atp2a2 distal hits chr5:122454200-122454250 : {'PASS' if ok1 else 'FAIL'}")
print(f"      Atp2a2 distal peaks: {[h['peak'] for h in d_hits]}")

s_hits = shared.get(a["event_id"], [])
def covers_550_650(hits):
    # old-cited region chr5:122456550-122456650 = two adjacent 50nt bins; test by coordinate coverage
    cov = sorted((h["peak"].split(":")[1].split("-")) for h in hits)
    cov = [(int(s), int(e)) for s, e in cov]
    return any(s <= 122456550 and e >= 122456600 for s, e in cov)
ok2 = covers_550_650(s_hits) and not covers_550_650(d_hits)
print(f"  [2] old-cited region chr5:122456550-122456650 only in shared : {'PASS' if ok2 else 'FAIL'}")
print(f"      shared hits there: {[h['peak'] for h in s_hits if 122456500 <= int(h['peak'].split(':')[1].split('-')[0]) <= 122457000]}")

n = by_sym.get("Ndrg2", [])[0]
g = by_sym.get("Agpat3", [])[0]
ok3 = distal.get(n["event_id"], []) == [] and distal.get(g["event_id"], []) == []
print(f"  [3] Ndrg2 distal peaks = {distal.get(n['event_id'], []) and [h['peak'] for h in distal[n['event_id']]] or 'none'}; "
      f"Agpat3 distal peaks = {distal.get(g['event_id'], []) and [h['peak'] for h in distal[g['event_id']]] or 'none'} -> {'PASS (old evidence gone)' if ok3 else 'CHECK'}")

n_total = sum(1 for r in out_rows if r["L2v2_distal_qki"].startswith("yes"))
n_f05 = sum(1 for r in out_rows if r["L2v2_distal_qki"] == "yes")
print(f"events with any distal QKI peak: {n_total}; with FDR<0.05: {n_f05} (v1 claim was 11/977 on wrong-strand segments)")
print(f"summary -> {outp}")
print("STANDARDS:", "ALL PASS" if (ok1 and ok2 and ok3) else "SEE ABOVE")
