# -*- coding: utf-8 -*-
"""Gate G4 closure + P5 progress in TODO."""
import io
p = r"D:\stroke_apa\TODO_pipeline_v5.1.md"
t = io.open(p, encoding="utf-8").read()

t = t.replace(
    "- [ ] **P4-G1** 候选数 0–5 皆合法——0 个满足预设条件时**输出零候选并报告失败层级**，不得为凑数降标准；非零时每张候选卡有完整四要素（事件/regulator+级别/class/三轴富集统计）【记录：n=____】",
    "- [x] **P4-G1** 候选数 0–5 皆合法……【2026-09-20 ✓ n=5（B 级，冻结规则自动选出）；四要素齐：results/p4_candidate_cards.md】")
t = t.replace(
    "- [ ] **P4-G2** 无合成数值分数混入（红线 A4）【记录：____】",
    "- [x] **P4-G2** 无合成数值分数混入（红线 A4）【2026-09-20 ✓ class 为可枚举证据组合（L2/L3/L4 层的 yes/no），无加权打分】")
t = t.replace(
    "- [ ] **P4-G3** 【用户确认记录：____】",
    "- [x] **P4-G3** 【2026-09-20 13:5x 用户确认记录：**\"确认批准\"**——Gate G4 通过，P4 关闭，P5 执行】")

t = t.replace(
    "### ✅ P4 完成记录：____",
    """### ✅ P4 完成记录（2026-09-20 关闭）
- 事件级证据矩阵 977 事件（P4-01/P3-04）；Δ3'UTR 区间层（P4-02）；CLIP 层 11 事件（P4-05）；Motif 层阴性照登（P4-06）；保守层 65%（P4-07）
- 判定（冻结规则，P4-08）：**A=0（MPRA 阻断）/ B=5 / C=475 / D=250 / none=247**
- **Top 候选（G4 通过）**：Ndrg2（d1/21/60 + d7 边缘）、Agpat3（d1/3/7/21）、Cnp（d1，WM 轴）、Fam107a（d1）、Pea15a（d7/60）——多时点稳定性排序：Ndrg2 ≈ Agpat3 > Cnp ≈ Pea15a > Fam107a""")

t = t.replace(
    "进度总览：`P0 [x] / P1 [x] / P2 [x] / P3 [x]（09-20 G3 通过）/ P4 [~]（候选卡出，待 G4）/ P5 [ ]`",
    "进度总览：`P0 [x] / P1 [x] / P2 [x] / P3 [x] / P4 [x]（09-20 G4 通过）/ P5 [~]（交付中）`")

io.open(p, "w", encoding="utf-8", newline="").write(t)
print("G4 closure recorded")
