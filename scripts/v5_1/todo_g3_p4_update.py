# -*- coding: utf-8 -*-
"""Gate G3 closure + P4 formalization in TODO (bulk edit)."""
import io

p = r"D:\stroke_apa\TODO_pipeline_v5.1.md"
t = io.open(p, encoding="utf-8").read()

# Gate G3 + P3 completion
t = t.replace(
    "- [ ] **P3-G1** 每个候选事件都有明确证据级别标注，无越级表述【记录：____】",
    "- [x] **P3-G1** 每个候选事件都有明确证据级别标注，无越级表述【2026-09-20 ✓ p3_p4_event_annotation_matrix.tsv 逐事件 L1/L2/L3/L4 + class 标注；全部表述按 L1/L2/L3 纪律】")
t = t.replace(
    "【记录：三轴各自=____】",
    "【2026-09-20 ✓ PAP 阳性（fgsea d7 padj≈1e-4、d3 0.0018、omnibus 0.0042、d1 0.044，多时点过 BH）；endfoot 阴性（day3 反向 +1.28 后全阴）；zone WM 弱同向 4 项、cortex 无——三轴分别表述，未合并】")
t = t.replace(
    "- [ ] **P3-G3** 【用户确认记录：____】",
    "- [x] **P3-G3** 【2026-09-20 13:4x 用户确认记录：**\"确认批准\"**——Gate G3 通过，P3 关闭，P4 执行】")
t = t.replace(
    "### ✅ P3 完成记录：____",
    """### ✅ P3 完成记录（2026-09-20 关闭）
- L1：12 对文献（P3-01）；L2：QKI 丢失段重叠 11/977（OR=2.71 p=0.013，P3-02）；L3：vs UTR 匹配 null 无富集（阴性照登，P3-03）；M4 三轴：PAP 阳性 / endfoot 阴性 / WM 弱同向（GSEA+Fisher 双层）
- P4-03 判据看交集前冻结（2a24dd7）；P3-04 矩阵 + 冻结规则 → **class B=5 / C=475 / D=250 / none=247**
- B 级 5 候选（全短缩、多时点、QKI 重叠）：Ndrg2 / Cnp / Fam107a / Pea15a / Agpat3""")

# P4 items formalization
t = t.replace(
    "→ `results/p4_evidence_matrix.tsv` + `results/p4_funnel_candidates.tsv`【记录：矩阵 n=____；描述性交集 n=____】",
    "→ `results/p3_p4_event_annotation_matrix.tsv`（977 事件，L1/L2/L3/L4+方向+三轴）+ `results/p4_top_b_class_candidates.tsv`【2026-09-20 ✓ 矩阵 n=977；描述性交集：APA∩PAP=161 / ∩endfoot=34 / ∩cortex=241 / ∩WM=429】")
t = t.replace(
    "- [ ] **P4-02 Δ3'UTR 区间提取**：每候选事件的 (proximal PAS, distal PAS) 区间，mm10【记录：____】",
    "- [x] **P4-02 Δ3'UTR 区间提取**：每候选事件的 (proximal PAS, distal PAS) 区间，mm10【2026-09-20 ✓ results/fimo_seqs/lost_segments.bed（1,541 缩短事件丢失段，1.68Mb）——即 FIMO 与保守性共用的区间层】")
t = t.replace(
    "- [ ] **P4-05 CLIP 层**：QKI/KHDRBS1/ELAVL1 等 peaks intersect【记录：____】",
    "- [x] **P4-05 CLIP 层**：QKI peaks intersect【2026-09-20 ✓ = P3-02：11/977（OR=2.71 p=0.013，10/11 缩短）；第二 CLIP 源（POSTAR3/KHDRBS1/ELAVL1）登记为局限未做】")
t = t.replace(
    "- [ ] **P4-06 Motif 层**：gain/loss 汇总【记录：____】",
    "- [x] **P4-06 Motif 层**：gain/loss 汇总【2026-09-20 ✓ **阴性照登**：vs UTR 匹配 null 无家族富集（0.89–1.07，无过 BH，5 家族轻度缺失）；逐事件描述表 p3_level3_motif.tsv；冻结规则下 L3=730/977 事件有 P1 家族强命中】")
t = t.replace(
    "- [ ] **P4-07 Conservation 层**：phastCons 60-way（mm10）平均保守分或保守性二值【记录：____】",
    "- [x] **P4-07 Conservation 层**：phastCons 60-way（mm10）【2026-09-20 ✓ 二值口径（冻结规则）：丢失段与保守元件重叠≥10bp = 1,004/1,541（65%）；管线：bigWig（Windows 代理 540KB/s 下载）→ bedGraph → 元件合并 6,788 万 → python 有序扫描（bedtools OOM×2 弃用）】")
t = t.replace(
    "- [ ] **P4-08 A/B/C/D 判定表**：每候选一张四层证据矩阵 → `results/p4_evidence_class.tsv`【记录：A=____ B=____ C=____ D=____】",
    "- [x] **P4-08 A/B/C/D 判定表**：每候选一张四层证据矩阵 → `results/p3_p4_event_annotation_matrix.tsv`（class 列）【2026-09-20 ✓ **A=0（MPRA 阻断，登记）B=5 C=475 D=250 none=247**；B 候选 = Ndrg2/Cnp/Fam107a/Pea15a/Agpat3（p4_top_b_class_candidates.tsv）】")

io.open(p, "w", encoding="utf-8", newline="").write(t)
print("P4 formalization done")
