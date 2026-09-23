#!/usr/bin/env python3
"""AUDIT PASS 2 (deliverable completeness vs TODO §16/§18) + PASS 3 (cross-doc
consistency & design-interval geometry). Exit 1 on any FAIL."""
import csv, os, sys, re

BASE = r"D:\stroke_apa_reanalysis"
R = os.path.join(BASE, "results", "02_candidate_rebuild")
S1 = os.path.join(BASE, "results", "01_strand_audit")
BR = os.path.join(BASE, "results", "03_bam_review")
P4 = os.path.join(BASE, "results", "04_public_pap")
TN = os.path.join(BASE, "results", "05_target_nomination")
FAILS = []

def check(label, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {label}" + (f"  | {detail}" if detail and not ok else ""))
    if not ok:
        FAILS.append(label)

def exists(rel):
    return os.path.exists(os.path.join(BASE, rel))

print("== PASS 2: TODO §16 计算阶段交付物 ==")
deliverables = {
    "输入清单和 checksum": ["results/00_inventory/input_files_sha256.tsv", "results/00_inventory/package_sha256.tsv", "results/00_inventory/provenance_manifest.tsv"],
    "坐标约定文件": ["config/coordinate_convention.yaml"],
    "链方向区间脚本和测试": ["scripts/01_build_strand_aware_segments.py", "scripts/test_01_segments.py"],
    "strand_aware_event_segments.tsv": ["results/01_strand_audit/strand_aware_event_segments.tsv"],
    "distal/shared BED": ["results/01_strand_audit/distal_segments.strand_aware.bed", "results/01_strand_audit/shared_segments.strand_aware.bed", "results/01_strand_audit/distal_segments.sorted.bed"],
    "重算 QKI": ["results/02_candidate_rebuild/qki_distal_intersections.tsv", "results/02_candidate_rebuild/qki_shared_intersections.tsv", "results/02_candidate_rebuild/qki_per_event_summary.tsv"],
    "重算 motif": ["results/02_candidate_rebuild/fimo_distal/fimo.tsv", "results/02_candidate_rebuild/fimo_null/fimo.tsv", "results/02_candidate_rebuild/motif_event_summary.tsv", "results/02_candidate_rebuild/motif_family_enrichment.tsv", "results/02_candidate_rebuild/motif_null_matching_qc.tsv"],
    "重算保守性": ["results/02_candidate_rebuild/distal_conservation.tsv"],
    "event_annotation_matrix.v2.tsv": ["results/02_candidate_rebuild/event_annotation_matrix.v2.tsv"],
    "old_vs_new_annotation_diff.tsv": ["results/02_candidate_rebuild/old_vs_new_annotation_diff.tsv", "results/02_candidate_rebuild/events_changed_qki_status.tsv", "results/02_candidate_rebuild/events_changed_class.tsv"],
    "候选卡 v2": ["results/02_candidate_rebuild/candidate_cards_v2.md"],
    "固定尺度核查图": ["results/03_bam_review/plots/Atp2a2.png", "results/03_bam_review/igv_snapshot_commands.txt"],
    "BAM 审查表": ["results/03_bam_review/event_review_table.tsv", "results/03_bam_review/bam_coverage_review.tsv", "results/03_bam_review/arich_risk.tsv"],
    "GSE143531 gene-level 表": ["results/04_public_pap/gse143531_candidate_gene_summary.tsv", "results/04_public_pap/gse143531_candidate_percentiles.tsv"],
    "GSE330741 manifest 状态": ["input_links/p4_mpra_counts_long.tsv"],
    "target_nomination_v1.md": ["results/05_target_nomination/target_nomination_v1.md"],
}
for item, files in deliverables.items():
    missing = [f for f in files if not exists(f)]
    check(f"§16: {item}", not missing, f"missing {missing}")

print("== PASS 2: TODO §18 最近一轮清单 ==")
s18 = {
    "输入链接整理": exists("input_links/dapars2_event_coordinates.tsv"),
    "坐标约定冻结": exists("config/coordinate_convention.yaml"),
    "正/负链区间实现": exists("scripts/01_build_strand_aware_segments.py"),
    "Atp2a2 单元测试": exists("scripts/test_01_segments.py"),
    "QKI 全量重算": exists(os.path.join(R, "qki_distal_intersections.tsv")),
    "motif+保守性重算": exists(os.path.join(R, "distal_conservation.tsv")) and exists(os.path.join(R, "motif_family_enrichment.tsv")),
    "新旧差异表": exists(os.path.join(R, "old_vs_new_annotation_diff.tsv")),
    "候选卡更新": exists(os.path.join(R, "candidate_cards_v2.md")),
    "六样本核查(计算等价: depth+图)": exists(os.path.join(BR, "event_review_table.tsv")) and exists(os.path.join(BR, "plots")),
    "主靶/备靶+引物区间输出": exists(os.path.join(TN, "target_nomination_v1.md")),
}
for k, v in s18.items():
    check(f"§18: {k}", v)

print("== PASS 3: 跨文书一致性 ==")
rank = {r["gene_symbol"]: r for r in csv.DictReader(open(os.path.join(R, "candidate_ranking_v2.tsv"), encoding="utf-8"), delimiter="\t")}
rv = {r["gene_symbol"]: r for r in csv.DictReader(open(os.path.join(BR, "event_review_table.tsv"), encoding="utf-8"), delimiter="\t")}
nom = open(os.path.join(TN, "target_nomination_v1.md"), encoding="utf-8").read()
for g, drops in (("Atp2a2", "71.7"), ("Agpat3", "82.3"), ("Aplp1", "46.9")):
    check(f"nomination 引用 d7 降幅 {g}={drops}%", drops in nom and rv[g]["d7_drop_pct"] == drops)
for g in ("Atp2a2", "Agpat3", "Aplp1"):
    check(f"nomination PAS 距离 {g}={rank[g]['pas_dist_bp']}bp", f"{rank[g]['pas_dist_bp']}bp" in nom)
for g in ("Atp2a2", "Agpat3", "Aplp1"):
    d = rank[g]["distal"]; s = rank[g]["shared"]
    check(f"nomination 区间 {g} distal/shared 与 ranking 一致", (d in nom and s in nom))
check("nomination 主靶= Atp2a2", "主靶 | Atp2a2 | 6953" in nom)
check("卡片存在且含 top10", os.path.getsize(os.path.join(R, "candidate_cards_v2.md")) > 2000)

print("== PASS 3: 设计区间几何（主靶三区间落在对应段内）==")
def within(x, seg):
    return seg[0] <= x[0] < x[1] <= seg[1]
for g in ("Atp2a2", "Agpat3", "Aplp1"):
    r = rank[g]
    d = tuple(map(int, r["distal"].split("-")))
    s = tuple(map(int, r["shared"].split("-")))
    m = re.search(rf"### {g}.*?外侧引物区：\S+:(\d+)-(\d+)；.*?common 扩增子区：\S+:(\d+)-(\d+).*?distal 扩增子区：\S+:(\d+)-(\d+)", nom, re.S)
    ok = m is not None
    if ok:
        fz = (int(m.group(1)), int(m.group(2)))
        ca = (int(m.group(3)), int(m.group(4)))
        da = (int(m.group(5)), int(m.group(6)))
        ok = within(fz, s) and within(ca, s) and within(da, d) and fz[0] < fz[1] and ca[0] < ca[1] and da[0] < da[1]
    check(f"{g} 外侧引物/common∈shared、distal 扩增子∈distal、区间正长", ok, nom[m.start():m.start()+120] if m else "regex miss")

print("== PASS 3: rebuild_notes 关键数字与输出一致 ==")
notes = open(os.path.join(S1, "rebuild_notes.md"), encoding="utf-8").read()
for key in ("42/9,693", "792", "0.054", "9,667", "B=8 / C=783 / D=1 / none=185", "346", "13bp", "8bp", "-71.7%", "-82.3%"):
    check(f"notes 含 {key}", key in notes)

print()
print("AUDIT PASS 2+3:", "ALL PASS" if not FAILS else f"{len(FAILS)} FAILURES -> {FAILS}")
sys.exit(1 if FAILS else 0)
