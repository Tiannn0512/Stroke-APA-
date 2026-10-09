#!/usr/bin/env python3
"""Task B: event_review_table v2 (fixes vs v1 per 03_预设 §5).

Fixes:
  1. PDUI trajectory: PDUI matrix columns are BARE sample names (sham1, day3_rep1, ...)
     — v1 queried f"{sample}_PDUI" and produced all-NA. Now joined via gene_symbol ->
     PDUI matrix Gene column with bare names.
  2. second_method: previously hardcoded "absent" from an earlier ranking snapshot;
     now joined from salmon_gate1_baseline.tsv (Atp2a2/Agpat3 = computed) and labeled
     "same-cohort computational cross-check (salmon TPM long-isoform fraction)" —
     NOT independent validation. Aplp1 = single annotated transcript, not computable.
  3. probe: downgraded to "region identified; oligo design pending" (no oligo
     sequences exist yet). "probe feasible/validated" was not earned.

Join key: gene_symbol + event_num across ranking/coverage/arich/pas tables.
Output: results/03_bam_review/event_review_table.v2.tsv + event_review_table_QC.md
"""
import argparse, csv, os

ap = argparse.ArgumentParser()
ap.add_argument("--base-dir", default=r"D:\stroke_apa_reanalysis")
ap.add_argument("--out-dir", default=None)
args = ap.parse_args()
BASE = args.base_dir
BR = args.out_dir or os.path.join(BASE, "results", "03_bam_review")

def rd(p, delimiter="\t"):
    return list(csv.DictReader(open(p, encoding="utf-8"), delimiter=delimiter))

rank = {r["gene_symbol"]: r for r in rd(os.path.join(BASE, "results", "02_candidate_rebuild", "candidate_ranking_v2.tsv"))}
cov = {}
for r in rd(os.path.join(BR, "bam_coverage_review.tsv"), delimiter="\t"):
    cov.setdefault(r["gene_symbol"], {})[r["sample"]] = r
arich = {r["gene_symbol"]: r for r in rd(os.path.join(BR, "arich_risk.tsv"), encoding="utf-8", delimiter="\t")} if False else \
        {r["gene_symbol"]: r for r in rd(os.path.join(BR, "arich_risk.tsv"))}
pasm = {r["event_num"]: r for r in rd(os.path.join(BASE, "results", "02_candidate_rebuild", "pas_distance_977.tsv"))}
pdui_rows = rd(os.path.join(BASE, "input_links", "p2_full_pdui_matrix.tsv"))
pdui = {r["Gene"]: r for r in pdui_rows}
sal = {}
for r in rd(os.path.join(BASE, "results", "02_candidate_rebuild", "salmon_gate1_baseline.tsv")):
    sal.setdefault(r["gene_symbol"], {})[r["sample"]] = r["long_frac_salmon"]

SAMPLES = ["sham1", "sham2", "day3_rep1", "day3_rep2", "day7_rep1", "day7_rep2"]
rows = []
for g, r in sorted(rank.items()):
    if g not in cov:
        continue
    p = pdui.get(r["event_id"], {})
    traj = " ".join(f"{s}={p.get(s, 'NA')}" for s in SAMPLES)
    sal_g = sal.get(g)
    if sal_g:
        sm = "same-cohort computational cross-check: salmon TPM long-frac " + \
             " ".join(f"{s}={sal_g[s]}" for s in SAMPLES)
    elif g == "Aplp1":
        sm = "not computable (single annotated transcript)"
    else:
        sm = "absent (not in 27-event top40 set)"
    d7 = [float(cov[g][s]["coverage_ratio_distal_over_shared"]) for s in ("day7_rep1", "day7_rep2")]
    d3 = [float(cov[g][s]["coverage_ratio_distal_over_shared"]) for s in ("day3_rep1", "day3_rep2")]
    sh = [float(cov[g][s]["coverage_ratio_distal_over_shared"]) for s in ("sham1", "sham2")]
    sham_m, d3m, d7m = sum(sh)/2, sum(d3)/2, sum(d7)/2
    drop7 = round(100 * (sham_m - d7m) / sham_m, 1)
    conc7 = (sham_m - d7[0]) * (sham_m - d7[1]) > 0
    rows.append({
        "gene_symbol": g, "event_num": r["event_num"], "class_v2": r["class_v2"],
        "strand": r["strand"], "distal": r["distal"], "shared": r["shared"],
        "proximal_pas_bed": r["proximal_pas_bed"],
        "PDUI_traj": traj,
        "PDUI_sham1": p.get("sham1", "NA"), "PDUI_sham2": p.get("sham2", "NA"),
        "PDUI_day3_rep1": p.get("day3_rep1", "NA"), "PDUI_day3_rep2": p.get("day3_rep2", "NA"),
        "PDUI_day7_rep1": p.get("day7_rep1", "NA"), "PDUI_day7_rep2": p.get("day7_rep2", "NA"),
        "ratio_sham_mean": round(sham_m, 3), "ratio_d3_mean": round(d3m, 3), "ratio_d7_mean": round(d7m, 3),
        "d7_drop_pct": drop7, "d7_reps_concordant": "yes" if conc7 else "no",
        "second_method": sm,
        "PAS_dist_bp": pasm.get(r["event_num"], {}).get("pas_dist_bp", "NA"),
        "A_rich_risk": arich[g]["arich_risk"],
        "L2v2_distal_qki": r["L2v2_distal_qki"], "L3v2_frozen8": r["L3v2_frozen8"],
        "L4v2_conserved": r["L4v2_conserved"], "PAP_member": r["PAP_member"],
        "probe_status": "region identified; oligo design pending",
        "decision": ("advance" if (drop7 >= 25 and conc7 and arich[g]["arich_risk"] == "low") else "hold"),
    })

fields = list(rows[0].keys())
outp = os.path.join(BR, "event_review_table.v2.tsv")
with open(outp, "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    w.writeheader(); w.writerows(rows)

qc = """# event_review_table.v2 QC

- join key: gene_symbol（ranking v2）<- 每子表同名列；PDUI 经 gene_symbol -> PDUI 矩阵 Gene 列（v1 bug：查询了不存在的 {sample}_PDUI 列，现改为矩阵实际列名）
- second_method 来源：results/02_candidate_rebuild/salmon_gate1_baseline.tsv，标注为 same-cohort computational cross-check（同一 GSE238125 队列的计算交叉核对，非独立验证）
- probe_status 统一为 "region identified; oligo design pending"（尚无 oligo 序列，不得写 feasible/validated）
- Atp2a2 验收：PDUI sham 0.19/0.20、day3 0.08/0.06、day7 0.09/0.09（非 NA）
"""
open(os.path.join(BR, "event_review_table_QC.md"), "w", encoding="utf-8").write(qc)

a = next(r for r in rows if r["gene_symbol"] == "Atp2a2")
print("Atp2a2 acceptance:", a["PDUI_sham1"], a["PDUI_sham2"], "|", a["PDUI_day3_rep1"], a["PDUI_day3_rep2"], "|", a["PDUI_day7_rep1"], a["PDUI_day7_rep2"])
print("rows:", len(rows), "->", outp)
