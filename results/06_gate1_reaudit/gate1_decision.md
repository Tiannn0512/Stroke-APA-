# Gate 1 决定书：Atp2a2（event 6953）定靶复核 — 2026-09-24

## 决定

**`advance Atp2a2 to independent 3'RACE`**——
计算复核与引物特异性检查全部通过（含 Primer-BLAST 网页复核，见 §4）；
assay 状态仍不写 `validated`（validated 需湿实验确认），引物订单由实验负责人终审后提交。

按任务书判据逐条对照：

| 任务书 `advance` 判据 | 本轮状态 |
|---|---|
| 修正后的读段和分区覆盖仍支持核查该位点 | ✅ 旧数字（PDUI 三方法、覆盖比 0.302→0.085、0.190→0.086）全部独立复现；修正后证据反而增强（见下） |
| 剪接结构已被具体描述 | ✅ 区域内全部主要接头精确注释到转录本：主接头 122454304-122455892 = 179939 的 spliced-UTR 终端内含子；次接头 122454304-122457302 = 177974 终端内含子；三条高频接头 = 共有内含子 #2-#4。无未解释主接头 |
| 没有足以单独解释全部信号的明确伪影 | ✅ ①比对质量干净（MAPQ 255 为主 ≥99.6%，soft-clip 碱基 ≤0.82%）；②纯 pre-mRNA 指示器（跨-UTR 内含子区覆盖）卒中后平稳，排除转录活性下降；③剔除 179939 跨内含子读段后 distal 覆盖仅降 ~7%，排除"剪接读段伪影"单因解释；④旧 spliced=0% 系脚本 bug，勘误后剪接不是伪影而是异构体证据 |
| 实际引物坐标、方向和特异性检查通过 | ✅ 坐标/方向：FASTA 提取 bug 已修（pyfaidx + samtools faidx 双实现逐碱基对拍，单测覆盖行内偏移/跨行窗口），8 条寡核苷酸全部在验证序列上重设计，坐标经独立实现对拍。特异性：转录组（142,604 条记录）逐引物全长精确命中 1-5 个且全部为 Atp2a2 基因转录本；引物对扫描 0 个非 Atp2a2 扩增子；加工假基因 Gm30970 已通过设计选址规避；**Primer-BLAST 网页复核已完成**：common/distal 基因组库仅预期靶位点、零潜在脱靶（结果页存档 primer_blast/） |

## 本轮必须记录的三个新事实（影响叙事，不影响 advance 判定）

1. **长 3'UTR 是跨外显子的**：179939 的 3'UTR 被 1588bp 内含子（122454304-122455892）分为
   终端外显子 UTR（122453512-122454304）与倒数第二外显子 UTR（122455892-122457153）两块；
   DaPars2 事件 6953 的 distal 段含此内含子。distal 覆盖的变化因此混合了
   "成熟分子消失"与"pre-mRNA"，本轮已用区段拆分分离（内含子段平稳 → 成熟分子改变）。
2. **读段级异构体证据**：179939 终端内含子接头支持 sham 96/136 → d3 4/11 → d7 2/7（-96%），
   177974 终端内含子 13/18 → 8/15（归一后约 -15%）——长短异构体应答分离，
   这是比覆盖比更直接的证据（短读段局限：终端外显子内部覆盖无法按异构体归属）。
3. **Atp2a2 共 5 条转录本**（含 196490/197415）与一个加工假基因 Gm30970（chr14，363bp）；
   引物设计已按"精确命中全部限于 Atp2a2 基因"标准规避假基因覆盖窗（mRNA ~3156-3330）。

## 未完成项（进入湿实验前）

| 项 | 状态 | 负责 |
|---|---|---|
| Primer-BLAST 网页近似匹配/成对扩增检查 | ✅ 已完成（2026-09-24）：common/distal 基因组库仅预期靶位点零脱靶；junction 三种角色组合均无产物（库局限，已记录），其特异性由转录组精确匹配+尺寸区分+RACE 测序三重保证 | 本轮已完成 |
| 引物订单审核 | 待实验负责人终审（assay 状态 = primer-blast checked，非 validated） | 实验负责人 |
| 3'RACE 测序确认 RNA-poly(A) 接合 | 未开始（本轮之后启动） | 实验负责人 |
| 短版本特异性引物定稿 | 按任务书要求，待真实 3'端结构确认后 | 湿实验阶段 |

## 证据文件索引

corrected_read_metrics.tsv（勘误后指标）/ junction_counts_by_sample.tsv /
junction_transcript_annotation.tsv / segment_coverage_by_sample.tsv /
coverage_sensitivity.md（含独立复算与替代解释）/ igv/（15 张真 IGV 截图）/
depth_plot_*.png（DEPTH PLOT）/ assay_design_corrected.tsv / primer_specificity_hits.tsv /
input_manifest.tsv / software_versions.txt / test_g1_1.py、test_g1_3_seq.py（自动测试）。
旧文件未覆盖：v1 指标表、event_review_table.v2、decision_summary_v1、旧 edgeR CSV 文件均在
inputs/zip_2023 原位保留。
