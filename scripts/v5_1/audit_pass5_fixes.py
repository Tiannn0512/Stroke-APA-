# -*- coding: utf-8 -*-
"""Audit pass-5 fixes: stale G5 status, L1 terminology collision, misleading 三轴全中,
missing gene-name list + genome-wide consistency stat."""
import io

p = r"D:\stroke_apa\STROKE_APA_FINAL_REPORT.md"
t = io.open(p, encoding="utf-8").read()

# a. stale G5 status
t = t.replace("| P5 | 交付（本报告 + 候选卡 + 湿实验路线书） | ✅ 09-20 | **G5 待确认** |",
              "| P5 | 交付（本报告 + 候选卡 + 湿实验路线书） | ✅ 09-20 | **G5**（用户确认） |")

# b. L1 terminology collision (Layer-1 scRNA vs evidence-level L1)
t = t.replace("Pabpn1↓（L1 padj 3e-4，homeostatic 0.017 重复）",
              "Pabpn1↓（scRNA 层 padj 3e-4，homeostatic 层 0.017 重复）")

# c. Ndrg2 misleading 三轴全中
t = t.replace("保守重叠 447bp；三轴全中",
              "保守重叠 447bp；属于 endfoot/cortex/WM 三个响应集（注意：**不在**信号阳性的 PAP 轴清单内）")

# d. gene list reference + genome-wide consistency + 11 gene names
t = t.replace(
    "| 4 | 特异事件结合证据 | **11/977 事件丢失段含 QKI 真实结合位点**（OR=2.71，p=0.013，10/11 缩短） | 弱-中（单事件功效有限；top40 一致性 13/27 = 功效不足） | 事件级真实性需 3'RACE |",
    "| 4 | 特异事件结合证据 | **11/977 事件丢失段含 QKI 真实结合位点**（OR=2.71，p=0.013，10/11 缩短）：Agpat3/Cnp/Fam107a/Ndrg2/**Qk**/Pea15a/Igfbp7/Sparcl1/Atp2a2/Sirt2/Cpe。方法学一致性：salmon 与 DaPars2 两独立管线方向全时点一致（Spearman 0.041–0.139，5/5 时点 p≤0.0023） | 弱-中（单事件功效有限；top40 一致性 13/27 = 功效不足） | 事件级真实性需 3'RACE |")
t = t.replace(
    "候选卡全文（含坐标/motif p 值/保守 bp）= `results/p4_candidate_cards.md`。",
    "候选卡全文（含坐标/motif p 值/保守 bp）= `results/p4_candidate_cards.md`；**977 基因完整清单**（含 padj/方向/class/QKI/PAP 归属）= `results/final_answer_gene_list.tsv`。")

io.open(p, "w", encoding="utf-8", newline="").write(t)
print("pass-5 fixes applied")
