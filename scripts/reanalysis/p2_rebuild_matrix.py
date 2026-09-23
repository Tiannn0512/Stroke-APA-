#!/usr/bin/env python3
"""P2.5: rebuild event annotation matrix v2, old-vs-new diffs, candidate cards v2.

Scope: the 977 significant APA events (union of sig_joint TRUE across contrasts).
Strand-aware annotations: L2v2 (distal QKI peak, FDR<0.05), L3v2 (frozen-8-family
FIMO p<1e-4 in distal), L4v2 (>=10bp score>0 phastCons in distal).
Class rules (frozen P4-03): A unreachable (MPRA blocked); B=L2&L3&L4; C=L3&L4; D=L3; else none.

Ranking (TODO §5 P2.5): class -> n sig contrasts -> max|dPDUI| -> PAP membership.
Outputs: event_annotation_matrix.v2.tsv, candidate_ranking_v2.tsv,
         old_vs_new_annotation_diff.tsv, events_changed_qki_status.tsv,
         events_changed_class.tsv, candidate_cards_v2.md, rebuild_session_info.txt
"""
import csv, os, datetime

BASE = r"D:\stroke_apa_reanalysis"
IL = os.path.join(BASE, "input_links")
R = os.path.join(BASE, "results", "02_candidate_rebuild")
S1 = os.path.join(BASE, "results", "01_strand_audit")

FROZEN8 = {"QKI", "ELAVL1", "PTBP1", "MBNL2", "NOVA1", "NOVA2", "TARDBP", "HNRNPA2B1"}

def read_tsv(p):
    return list(csv.DictReader(open(p, encoding="utf-8"), delimiter="\t"))

coord = read_tsv(os.path.join(IL, "dapars2_event_coordinates.tsv"))
num2name = {r["event_num"]: r["event_id"] for r in coord}
name2num = {r["event_id"]: r["event_num"] for r in coord}

sap = read_tsv(os.path.join(IL, "p2_sap_layer2_events.tsv"))
sig = {}
for r in sap:
    if r["sig_joint"] != "TRUE":
        continue
    e = sig.setdefault(r["event"], {"contrasts": [], "dpduis": [], "pads": [], "conc": []})
    e["contrasts"].append(r["contrast"])
    e["dpduis"].append(float(r["dPDUI"]))
    e["pads"].append(float(r["padj_joint"]))
    e["conc"].append(float(r["d_rep1"]) * float(r["d_rep2"]) > 0)

segs = {r["event_num"]: r for r in read_tsv(os.path.join(S1, "strand_aware_event_segments.tsv"))}
qki = {r["event_id"]: r for r in read_tsv(os.path.join(R, "qki_per_event_summary.tsv"))}
mot = {r["event_num"]: r for r in read_tsv(os.path.join(R, "motif_event_summary.tsv"))}
cons = {r["event_num"]: r for r in read_tsv(os.path.join(R, "distal_conservation.tsv"))}
pap = {line.strip() for line in open(os.path.join(IL, "pap_gene_sets_set1_pap2729.txt"), encoding="utf-8") if line.strip()}
old = {r["event"]: r for r in read_tsv(os.path.join(IL, "p3_p4_event_annotation_matrix.tsv"))}
top40 = {r["symbol"] for r in read_tsv(os.path.join(IL, "p2_second_tool_top40_signs.tsv"))}

pdui = {r["Gene"]: r for r in read_tsv(os.path.join(IL, "p2_full_pdui_matrix.tsv"))}
SAMPLES = ["sham1", "sham2", "day3_rep1", "day3_rep2", "day7_rep1", "day7_rep2",
           "day21_rep1", "day21_rep2", "day60_rep1", "day60_rep2"]

def cls_v2(l2, l3, l4):
    if l3 == "yes" and l4 == "yes":
        return "B" if l2 == "yes" else "C"
    if l3 == "yes":
        return "D"
    return "none"

matrix, diffs, dqki, dcls, ranking = [], [], [], [], []
for num, e in sig.items():
    sg = segs.get(num)
    if sg is None:
        continue
    q = qki.get(num2name.get(num, ""))
    m = mot.get(num)
    c = cons.get(num)
    l2 = (q or {}).get("L2v2_distal_qki", "no")
    l2 = "yes" if l2 == "yes" else "no"
    l3 = (m or {}).get("L3v2_frozen8", "no")
    l4 = (c or {}).get("L4v2_conserved", "NA")
    l4eff = l4 if l4 in ("yes", "no") else "no"
    cls = cls_v2(l2, l3, l4eff)
    ds = e["dpduis"]
    row = {
        "event_num": num, "event_id": num2name.get(num, ""), "gene_symbol": sg["gene_symbol"],
        "chrom": sg["chrom"], "strand": sg["strand"],
        "distal": f"{sg['distal_start']}-{sg['distal_end']}", "distal_len": sg["distal_len"],
        "shared": f"{sg['shared_start']}-{sg['shared_end']}", "shared_len": sg["shared_len"],
        "proximal_pas_bed": sg["proximal_pas_bed"],
        "n_sig_contrasts": len(e["contrasts"]), "contrasts_sig": ";".join(e["contrasts"]),
        "replicates_concordant_all": "yes" if all(e["conc"]) else "no",
        "max_abs_dPDUI": round(max(abs(x) for x in ds), 3),
        "min_padj": round(min(e["pads"]), 4),
        "L2v2_distal_qki": l2,
        "L2v2_peaks": (q or {}).get("distal_peaks", ""),
        "L3v2_frozen8": l3, "L3v2_families": (m or {}).get("frozen8_families_hit", ""),
        "L3v2_best_p": (m or {}).get("frozen8_best_p", ""),
        "L4v2_conserved": l4,
        "L4v2_pos_bases": (c or {}).get("distal_pos_bases", "NA"),
        "distal_mean_phastcons": (c or {}).get("distal_mean_phastcons", "NA"),
        "class_v2": cls, "PAP_member": "yes" if sg["gene_symbol"] in pap else "no",
        "second_method": "direct" if sg["gene_symbol"] in top40 else "absent",
        "probe_designability": "ok_pending_P3" if int(sg["distal_len"]) >= 150 and int(sg["shared_len"]) >= 60 else "constrained",
    }
    matrix.append(row)
    o = old.get(num)
    if o:
        d = {"event_num": num, "gene_symbol": sg["gene_symbol"],
             "old_L2": o.get("L2_qki_lostseg", ""), "new_L2v2": l2,
             "old_L3": o.get("L3_motif_frozen", ""), "new_L3v2": l3,
             "old_L4": o.get("L4_conserved", ""), "new_L4v2": l4,
             "old_class": o.get("class", ""), "new_class_v2": cls}
        diffs.append(d)
        if (o.get("L2_qki_lostseg", "").lower() == "yes") != (l2 == "yes"):
            dqki.append({"event_num": num, "gene_symbol": sg["gene_symbol"],
                         "old_L2": o.get("L2_qki_lostseg", ""), "new_L2v2": l2,
                         "new_distal_peaks": (q or {}).get("distal_peaks", "")})
        if o.get("class", "") != cls:
            dcls.append({"event_num": num, "gene_symbol": sg["gene_symbol"],
                         "old_class": o.get("class", ""), "new_class_v2": cls,
                         "old_L2": o.get("L2_qki_lostseg", ""), "new_L2v2": l2,
                         "old_L3": o.get("L3_motif_frozen", ""), "new_L3v2": l3,
                         "old_L4": o.get("L4_conserved", ""), "new_L4v2": l4})
    cls_rank = {"B": 0, "C": 1, "D": 2, "none": 3}[cls]
    ranking.append((cls_rank, -len(e["contrasts"]), -max(abs(x) for x in ds),
                    0 if row["PAP_member"] == "yes" else 1, sg["gene_symbol"], row))

def dump(path, rows, fields=None):
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields or list(rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(rows)

dump(os.path.join(R, "event_annotation_matrix.v2.tsv"), matrix)
dump(os.path.join(R, "old_vs_new_annotation_diff.tsv"), diffs)
dump(os.path.join(R, "events_changed_qki_status.tsv"), dqki)
dump(os.path.join(R, "events_changed_class.tsv"), dcls)

ranking.sort(key=lambda x: x[:4])
dump(os.path.join(R, "candidate_ranking_v2.tsv"), [r[5] for r in ranking])

# summary counts
from collections import Counter
cc = Counter(r["class_v2"] for r in matrix)
oc = Counter(d["old_class"] for d in diffs)
lines = [
    "# Candidate cards v2 (strand-aware rebuild) — generated " + datetime.datetime.now().isoformat(timespec="seconds"),
    "",
    f"Scope: {len(matrix)} significant events (977 union; {len(segs)} events had segments).",
    f"Class v2 counts: " + ", ".join(f"{k}={cc.get(k,0)}" for k in ("B", "C", "D", "none")) +
    f"  (v1 counts on wrong-strand segments: B=5 C=475 D=250 none=247 over 977)",
    f"Annotation flips vs v1: L2 {sum(1 for d in diffs if (d['old_L2'].lower()=='yes')!=(d['new_L2v2']=='yes'))}, "
    f"L3 {sum(1 for d in diffs if (d['old_L3'].lower()=='yes')!=(d['new_L3v2']=='yes'))}, "
    f"L4 {sum(1 for d in diffs if (d['old_L4'].lower() in ('yes',))!=(d['new_L4v2']=='yes'))}, "
    f"class {len(dcls)} of {len(diffs)} events with v1 annotation.",
    "",
    "## Top candidates (class -> n contrasts -> max|dPDUI| -> PAP)",
    "",
]
for rk in ranking[:10]:
    r = rk[5]
    oldcls = next((d["old_class"] for d in diffs if d["event_num"] == r["event_num"]), "n/a")
    pd = pdui.get(r["event_id"], {})
    traj = " ".join(f"{s}={pd.get(s+'_PDUI', 'NA')}" for s in SAMPLES[:6])
    lines += [
        f"### {r['gene_symbol']} ({r['event_num']}) — class {r['class_v2']} (v1: {oldcls})",
        f"- {r['chrom']} {r['strand']} | PAS(pred) {r['proximal_pas_bed']} | distal {r['distal']} ({r['distal_len']}bp) | shared {r['shared']} ({r['shared_len']}bp)",
        f"- sig contrasts: {r['contrasts_sig']} | max|dPDUI| {r['max_abs_dPDUI']} | min padj {r['min_padj']} | reps concordant: {r['replicates_concordant_all']}",
        f"- PDUI (sham/d3/d7): {traj}",
        f"- L2v2: {r['L2v2_distal_qki']} {r['L2v2_peaks']}",
        f"- L3v2: {r['L3v2_frozen8']} [{r['L3v2_families']}] {r['L3v2_best_p']}",
        f"- L4v2: {r['L4v2_conserved']} (pos bases {r['L4v2_pos_bases']}, mean {r['distal_mean_phastcons']})",
        f"- PAP: {r['PAP_member']} | second method: {r['second_method']} | probe: {r['probe_designability']}",
        "",
    ]
open(os.path.join(R, "candidate_cards_v2.md"), "w", encoding="utf-8").write("\n".join(lines))

with open(os.path.join(R, "rebuild_session_info.txt"), "w", encoding="utf-8") as f:
    f.write(f"generated {datetime.datetime.now().isoformat(timespec='seconds')}\n")
    f.write("python: scripts/p2_rebuild_matrix.py (stdlib only)\n")
    f.write("inputs: input_links sha256 (results/00_inventory/input_files_sha256.tsv)\n")
    f.write("rules: config/coordinate_convention.yaml; class frozen P4-03 (input_links/P4_03_evidence_class_frozen.md)\n")
    f.write("note: PAS-distance ranking field pending (PolyASite site unreachable 2026-09-21)\n")

print(f"matrix v2: {len(matrix)}; diffs: {len(diffs)}; class changed: {len(dcls)}; L2 changed: {len(dqki)}")
print("class v2 counts:", dict(cc))
print("cards: candidate_cards_v2.md (top 10)")
