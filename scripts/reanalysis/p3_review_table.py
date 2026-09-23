#!/usr/bin/env python3
"""P3.4: assemble the event review table (TODO P3.4 fields) + proposed decisions.

Decision rule (objective, documented — final call is Gate 1 / user):
  advance : d7 ratio drop >=25% vs sham mean AND d7 replicate-concordant
            AND arich low AND probe ok
  hold    : everything else (weak ratio change, discordant, low coverage,
            or its significant timepoint not in the reviewed BAM set)
  drop    : artifacts that explain the event (none found so far)
"""
import csv, os

BASE = r"D:\stroke_apa_reanalysis"
BR = os.path.join(BASE, "results", "03_bam_review")

rank = {r["event_num"]: r for r in csv.DictReader(open(os.path.join(BASE, "results", "02_candidate_rebuild", "candidate_ranking_v2.tsv"), encoding="utf-8"), delimiter="\t")}
cov = list(csv.DictReader(open(os.path.join(BR, "bam_coverage_review.tsv"), encoding="utf-8"), delimiter="\t"))
arich = {r["gene_symbol"]: r for r in csv.DictReader(open(os.path.join(BR, "arich_risk.tsv"), encoding="utf-8"), delimiter="\t")}
pasm = {r["event_num"]: r for r in csv.DictReader(open(os.path.join(BASE, "results", "02_candidate_rebuild", "pas_distance_977.tsv"), encoding="utf-8"), delimiter="\t")}
pdui = {r["Gene"]: r for r in csv.DictReader(open(os.path.join(BASE, "input_links", "p2_full_pdui_matrix.tsv"), encoding="utf-8"), delimiter="\t")}

by_gene = {}
for r in cov:
    by_gene.setdefault(r["gene_symbol"], {})[r["sample"]] = r

def mean_ratio(g, samples):
    vals = [float(by_gene[g][s]["coverage_ratio_distal_over_shared"]) for s in samples if by_gene[g][s]["coverage_ratio_distal_over_shared"] != "NA"]
    return sum(vals) / len(vals) if vals else float("nan")

rows = []
for num, r in rank.items():
    g = r["gene_symbol"]
    if g not in by_gene:
        continue
    d3m = mean_ratio(g, ["day3_rep1", "day3_rep2"])
    d7m = mean_ratio(g, ["day7_rep1", "day7_rep2"])
    sham = mean_ratio(g, ["sham1", "sham2"])
    drop7 = (sham - d7m) / sham if sham else 0
    r7 = [float(by_gene[g][s]["coverage_ratio_distal_over_shared"]) for s in ("day7_rep1", "day7_rep2")]
    d3 = [float(by_gene[g][s]["coverage_ratio_distal_over_shared"]) for s in ("day3_rep1", "day3_rep2")]
    conc7 = (sham - r7[0]) * (sham - r7[1]) > 0
    conc3 = (sham - d3[0]) * (sham - d3[1]) > 0
    tot_sham = (float(by_gene[g]["sham1"]["distal_mean_depth"]) + float(by_gene[g]["sham1"]["shared_mean_depth"]))
    tot_d7 = (float(by_gene[g]["day7_rep1"]["distal_mean_depth"]) + float(by_gene[g]["day7_rep1"]["shared_mean_depth"]))
    ar = arich[g]
    ps = pasm.get(num, {})
    p = pdui.get(r["event_id"], {})
    traj = " ".join(f"{s}:P={p.get(s + '_PDUI', 'NA')}" for s in ("sham1", "sham2", "day3_rep1", "day3_rep2", "day7_rep1", "day7_rep2"))
    if drop7 >= 0.25 and conc7 and ar["arich_risk"] == "low" and r["probe_designability"] == "ok_pending_P3":
        decision = "advance"
    elif ar["arich_risk"] == "HIGH":
        decision = "drop (internal priming risk)"
    else:
        decision = "hold"
    # confidence: high = big drop + concordant at both d3 and d7; medium = d7-concordant with adequate d7 depth; low otherwise
    tot_d7_mean = (float(by_gene[g]["day7_rep1"]["distal_mean_depth"]) + float(by_gene[g]["day7_rep1"]["shared_mean_depth"])
                   + float(by_gene[g]["day7_rep2"]["distal_mean_depth"]) + float(by_gene[g]["day7_rep2"]["shared_mean_depth"])) / 2
    tot_sham_mean = (float(by_gene[g]["sham1"]["distal_mean_depth"]) + float(by_gene[g]["sham1"]["shared_mean_depth"])
                     + float(by_gene[g]["sham2"]["distal_mean_depth"]) + float(by_gene[g]["sham2"]["shared_mean_depth"])) / 2
    if decision == "advance":
        if drop7 >= 0.40 and conc3:
            conf = "high"
        elif tot_d7_mean >= 150 and tot_sham_mean >= 200:
            conf = "medium"
        else:
            conf = "low (thin depth)"
    else:
        conf = "-"
    rows.append({
        "confidence": conf,
        "sham_total_depth_mean": round(tot_sham_mean, 1), "d7_total_depth_mean": round(tot_d7_mean, 1),
        "gene_symbol": g, "event_num": num, "class_v2": r["class_v2"], "chrom": r["chrom"], "strand": r["strand"],
        "distal": r["distal"], "shared": r["shared"], "proximal_pas_bed": r["proximal_pas_bed"],
        "n_sig_contrasts": r["n_sig_contrasts"], "contrasts_sig": r["contrasts_sig"],
        "PDUI_traj": traj,
        "ratio_sham_mean": round(sham, 3), "ratio_d3_mean": round(d3m, 3), "ratio_d7_mean": round(d7m, 3),
        "d7_drop_pct": round(100 * drop7, 1), "d7_reps_concordant": "yes" if conc7 else "no",
        "d3_reps_concordant": "yes" if conc3 else "no",
        "total_depth_sham_vs_d7": f"{tot_sham:.0f} -> {tot_d7:.0f}",
        "PAS_dist_bp": ps.get("pas_dist_bp", "NA"), "PAS_support": ps.get("pas_support", "NA"),
        "A_rich_risk": ar["arich_risk"], "L2v2": r["L2v2_distal_qki"], "L3v2": r["L3v2_frozen8"],
        "L4v2": r["L4v2_conserved"], "PAP": r["PAP_member"], "second_method": r["second_method"],
        "probe": r["probe_designability"], "decision": decision,
    })

fields = list(rows[0].keys())
with open(os.path.join(BR, "event_review_table.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    w.writeheader(); w.writerows(rows)

print(f"{'gene':<9} {'cls':<4} {'sham':>5} {'d3':>5} {'d7':>5} {'drop%':>6} {'d7con':>5} {'PASbp':>6} {'Ari':<5} {'dec'}")
for r in sorted(rows, key=lambda x: (x["decision"] != "advance", -int(x["n_sig_contrasts"]))):
    print(f"{r['gene_symbol']:<9} {r['class_v2']:<4} {r['ratio_sham_mean']:>5} {r['ratio_d3_mean']:>5} {r['ratio_d7_mean']:>5} "
          f"{r['d7_drop_pct']:>6} {r['d7_reps_concordant']:>5} {r['PAS_dist_bp'] or 'NA':>6} {r['A_rich_risk']:<5} {r['decision']}")
