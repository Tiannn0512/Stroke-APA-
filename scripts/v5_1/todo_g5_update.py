# -*- coding: utf-8 -*-
"""Gate G5 closure: project completion records in TODO."""
import io
p = r"D:\stroke_apa\TODO_pipeline_v5.1.md"
t = io.open(p, encoding="utf-8").read()

t = t.replace(
    "- [ ] **P5-G1** 全链证据矩阵闭合、top candidates 定稿【记录：____】",
    "- [x] **P5-G1** 全链证据矩阵闭合、top candidates 定稿【2026-09-20 ✓ 七环节链表（报告 §二）+ 5 张候选卡定稿（p4_candidate_cards.md）】")
t = t.replace(
    "- [ ] **P5-G2** 【用户确认记录：____】",
    "- [x] **P5-G2** 【2026-09-20 14:0x 用户确认记录：**\"确认批准\"**——Gate G5 通过，课题干实验阶段关闭】")

t = t.replace(
    "### ✅ P5 完成记录：____",
    """### ✅ P5 完成记录（2026-09-20 关闭）
- 最终报告 STROKE_APA_FINAL_REPORT.md/.html（七环节链表 + 五段结论 + 九项局限）
- 候选卡 p4_candidate_cards.md（5 张 B 级）；湿实验路线书 STROKE_APA_WETLAB_ROUTE.md v2.0（PAP 修订版）
- 红线 grep：12 禁语 0 真实违规（p5_redline_check.tsv）
- **Gate G5 通过：课题干实验阶段关闭（2026-09-20）**；湿实验阶段按路线书移交""")

t = t.replace(
    "进度总览：`P0 [x] / P1 [x] / P2 [x] / P3 [x] / P4 [x]（09-20 G4 通过）/ P5 [~]（交付中）`",
    "进度总览：`P0 [x] / P1 [x] / P2 [x] / P3 [x] / P4 [x] / P5 [x] —— **全部关闭（2026-09-20，Gate G5 通过）**`")

io.open(p, "w", encoding="utf-8", newline="").write(t)
print("G5 closure recorded")
