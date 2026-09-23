# Stroke-APA：卒中后星形胶质细胞 APA→突触旁突起定位课题（干实验全档案）

> 私有仓库 · 未发表课题档案 · 建库 2026-09-23
> 干实验阶段已关闭（v5.1 G1–G5 全过 → 链方向重建 → Gate 1 定靶 → ADVANCE **Atp2a2 → 3'RACE**），当前等待湿实验结果。

## 课题最终结论（一段话）

卒中后小鼠皮层星形胶质细胞发生广泛的 mRNA 3'UTR 变短（APA 位移，9,693 个定量事件中 977 个显著），候选调控子为 Pabpn1 与 Qk；主靶 **Atp2a2**（event 6953）满足五层证据：DaPars2 PDUI（sham 0.19/0.20→d7 0.09/0.09）、BAM 覆盖比（0.302→0.085，−71.7% 双重复一致）、salmon 长异构体占比（五时点同向）、远端丢失段含 QKI CLIP peak（FDR 1.04e-5）、预测断点距 PolyASite 已知 poly(A) 位点 13bp。文献（Sakers 2017，皮层）显示其野生型转录本突触旁富集 log2FC +2.13。**本课题自己的 PAP 数据核查（两个公共数据集）不支持"长短异构体分布不同"，该核心假说交由湿实验检验**（3'RACE → Gate 2A/2B → Gate 3）。结论上限表述为"altered localization potential"。

## 仓库地图

```
├── README.md                 ← 本文件
├── docs/
│   ├── 实验流程.md（在 docs/reanalysis/）★ 复现总手册：数据/软件版本/参数/预期结果锚点
│   ├── reanalysis/           ← Proposal、TODO、全程报告、重建报告、定靶决策
│   └── v5_1/                 ← v5.1 阶段全部文档、评审、工作日志（冻结参数在 results/v5_1/ 与 input_links/）
├── scripts/
│   ├── v5_1/                 ← P0–P6 全部脚本（QC/比对/DaPars2/limma/FIMO/定位层）
│   └── reanalysis/           ← R0–R4 + Task A–D 脚本（链修正/重建/定靶/引物）+ 单元测试
├── results/
│   ├── v5_1/                 ← P1 调控子层、P2 事件层+统计、P3 关联层、P4 定位层全部表与图
│   └── reanalysis/
│       ├── 00_inventory/     ← sha256 冻结清单、软件版本档案
│       ├── 01_strand_audit/  ← 链方向重建（9,693 事件区间）
│       ├── 02_candidate_rebuild/ ← QKI/motif/保守性/分类 v2（不含可再生大中间件）
│       ├── 03_bam_review/    ← 覆盖比审查、Atp2a2 读段核查（不含 IGV 切片 bam）
│       ├── 04_public_pap/    ← GSE143531 edgeR TMM 重算（Task A）
│       └── 05_target_nomination/ ← Atp2a2 引物设计与 ADVANCE 决策（Task D）
├── config/                   ← 冻结坐标约定（coordinate_convention.yaml）
├── input_links/              ← sha256 冻结输入（PDUI 矩阵、motif 库、QKI peaks 等 25 项）
├── data_frozen/              ← 小体积冻结数据（fastp QC、DaPars2 逐染色体输出、GEO 处理表）
├── wetlab/                   ← 湿实验工作目录骨架（3'RACE/异构体比例/成像/reporter）
├── envs/                     ← conda 环境导出
├── REGENERATE.md             ★ 数据再生手册：每一项未入库数据的确切重新获取/重建命令
└── scripts/audit_upload_coverage.py ← 上传覆盖审计脚本（可重跑验证）
```

## 交付包（GitHub Releases，不在仓库正文）

| Release 资产 | 内容 |
|---|---|
| `REANALYSIS_RESULTS_20260923.zip`（87MB） | 重建+定靶+Task A–D 最终交付（328 项） |
| `REANALYSIS_RESULTS_20260921.zip`（85MB） | 链方向重建期交付（313 项，含全部图与 IGV 切片） |
| `STROKE_APA_FINAL_DELIVERABLES_20260920.zip`（25MB） | v5.1 干实验收官交付 |
| `ATP2A2_IGV_CHECK_PACKAGE_20260920.zip`（29MB） | IGV 读段核查包 |

## 从哪里开始读

1. 只想知道结论 → 本 README + `docs/reanalysis/05…/decision_summary_v1.md`（在 `results/reanalysis/05_target_nomination/`）
2. 想复现 → `docs/reanalysis/实验流程.md`
3. 想查任何数字的来龙去脉 → `docs/reanalysis/PROJECT_FULL_REPORT.md`（全程报告，含代码详解与勘误总表）
4. 要重下数据/重建环境 → `REGENERATE.md`

## 纪律声明（沿用课题冻结规则）

一切 poly(A) 位点表述为"预测断点"；候选 ≠ 验证，结论上限为 altered localization potential；阴性结果（终足轴、L3 motif 富集、海马 PAP、d21/d60、GSE330741 manifest 矛盾）全部随档案保留；冻结文件（`results/v5_1/P1_frozen_params.md`、`input_links/P2_SAP_frozen.md`、`input_links/P4_03_evidence_class_frozen.md`、`config/coordinate_convention.yaml`）修改须登记理由。
