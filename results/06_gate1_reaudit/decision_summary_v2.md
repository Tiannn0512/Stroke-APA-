# decision_summary_v2.md — Atp2a2 定靶决策摘要（v2，Gate1 复核后）

日期：2026-09-24 ｜ 本版替代 `decision_summary_v1.md`（v1 保留原位不覆盖）
前置文书：`gate1_decision.md`（正式决定书）、`Atp2a2_read_review_v2.md`、`coverage_sensitivity.md`、`primer_specificity_review.md`

## 决定（不变，证据增强）

**ADVANCE Atp2a2 → 独立 3'RACE。** Primer-BLAST 网页复核已完成（非 pending）；
assay 状态 = `primer-blast checked`（`validated` 留待湿实验）。

## 相对 v1 的五处实质更新

1. **读段剪接指标勘误并升级为证据**：v1 引用的 spliced=0% 系解析 bug；修正后
   剪接接头本身就是异构体分辨证据——179939（长 UTR）终端内含子支持
   sham 96/136 → d3 4/11 → d7 2/7（-96%），177974 约 -15%（经共有内含子归一）。
2. **结构新事实**：长 3'UTR 跨外显子（179939 UTR 被 1588bp 内含子分为两段），
   事件 6953 的 distal 段含该内含子；覆盖解释必须区隔成熟分子 vs pre-mRNA
   （内含子段平稳 → 成熟分子群体改变成立）。
3. **IGV 证据重做**：15 张 IGV 2.18.4 真截图替代 matplotlib 深度图；
   旧"IGV 已核查、无异常剪接"表述撤回。
4. **引物整体重设计**：旧引物基于 buggy 序列提取，作废。新 8 条寡核苷酸
   （5 个 assay：RACE outer/nested GSP、common qPCR、distal-terminal qPCR、
   179939 特异 junction qPCR）在双实现对拍序列上设计，
   转录组精确匹配零脱靶、Primer-BLAST 基因组库零潜在脱靶（common/distal）、
   junction 仅 179939 扩增（366bp cDNA 产物，备靶 177974/31423 均不扩增）。
   Gm30970 加工假基因经设计选址规避。
5. **edgeR 文件格式修正**：`edgeR_TMM_PAP_vs_Full` 旧文件实为 CSV，已重导出真 TSV；
   GSE143531 数值与表述不变。

## 证据一览（全部可溯源）

| 层 | 证据 | 数值 |
|---|---|---|
| 定量 | DaPars2 PDUI（不变） | sham 0.19/0.20 → d7 0.09/0.09 |
| 覆盖 | 独立复算 distal/shared | 0.3025 → 0.0849 → 0.0855（旧 0.302→0.084→0.085） |
| 剪接 | 179939 终端内含子接头 | sham 96/136 → d3 4/11 → d7 2/7 |
| 剪接 | 177974 终端内含子接头 | sham 13/18 → d3 14/10 → d7 8/15 |
| pre-mRNA | 跨-UTR 内含子区覆盖比 | 0.067/0.071 → 0.085/0.081（平稳） |
| 文献 | Sakers 2017 皮层 PAP | Atp2a2 log2FC +2.13（引用，非本轮计算） |
| 特异性 | 转录组精确匹配 + Primer-BLAST | 0 非本基因扩增子（8/8 寡核苷酸） |

## 未完成项 / 边界

- Primer-BLAST 对 junction 跨内含子对的检索局限已记录（三种角色组合均无产物；
  该接合不存在于 RefSeq mRNA 模型）；湿实验以产物尺寸（366 vs 1954）+ 测序终审。
- 短版本特异性引物：待真实 3' 端确认后定稿（任务书既定）。
- PAP 中长短版本空间分布、QKI 因果、BBB 功能：未测定，不在本轮主张范围。
