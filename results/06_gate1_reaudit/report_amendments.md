# report_amendments.md — 两份全程报告的修订清单（Gate1 复核，2026-09-24）

按要求更新 `PROJECT_FULL_REPORT.md` 与 `REANALYSIS_FINAL_REPORT.md` 中受本轮结果影响的段落。
本文件逐条给出"原说法 → 修正说法 → 证据文件 → 日期"；两报告本体的下一版本号
（v2）应按本清单合入（旧版本不覆盖，新副本另存）。凡未列出的段落不受本轮影响。

| # | 报告位置 | 原说法 | 修正说法 | 证据 | 日期 |
|---|---|---|---|---|---|
| 1 | 两报告所有涉及 Atp2a2 读段剪接处 | "spliced_read_pct = 0%" / "无剪接读段" | 脚本 bug（CIGAR 条件含赋值表达式）导致恒 0。修正后：注释 UTR 内剪接读段比例 sham 4.28%/4.32% → d3 0.72%/1.02% → d7 0.30%/0.67% | corrected_read_metrics.tsv；test_g1_1.py 回归（sham1 修正后 506 条 vs 旧 0） | 2026-09-24 |
| 2 | R3/Task C "IGV 读段已核查" | `igv/Atp2a2_*.png` 为 matplotlib 深度曲线，被表述为 IGV 读段核查 | 撤回。本轮出具 IGV 2.18.4 真截图 15 张（六样本同视野 + 单样本 + PAS/接头双放大区） | results/06_gate1_reaudit/igv/（combined_junction_zoom.png 可见 sham 剪接弧粗密、d3/d7 显著减弱） | 2026-09-24 |
| 3 | "无异常剪接" | 未经验证的断言 | 改为："区域内全部主要接头精确注释（无未解释主接头）；主接头 = 179939 的 spliced-UTR 终端内含子，其支持读段 sham 96/136 → d3 4/11 → d7 2/7，为长 UTR 异构体成熟分子特异衰减的读段级证据" | junction_transcript_annotation.tsv；junction_counts_by_sample.tsv | 2026-09-24 |
| 4 | 事件 6953 distal 段描述 | distal 122453512-122456498 为连续 UTR 区段 | 补充结构事实：该段内含 179939 的 1588bp 跨-UTR 内含子（122454304-122455892）；distal 覆盖 = 终端外显子成熟分子 + pre-mRNA + 倒数第二外显子 UTR 前段三者混合；内含子段（纯 pre-mRNA 指示器）卒中后平稳（0.067/0.071→0.085/0.081），排除转录活性下降的替代解释 | 01_strand_audit 结构 + segment_coverage_by_sample.tsv + coverage_sensitivity.md | 2026-09-24 |
| 5 | Task D "零脱靶" | "零脱靶" | 改为："转录组全长精确匹配零脱靶（含 Atp2a2 全部 5 条转录本白名单、Gm30970 假基因规避）；引物对扫描无非 Atp2a2 扩增子；Primer-BLAST 基因组库 common/distal 仅预期靶位点。旧引物序列系 buggy 提取产物，已整体重设计作废" | primer_specificity_hits.tsv；primer_specificity_review.md；primer_blast/*.html | 2026-09-24 |
| 6 | `edgeR_TMM_PAP_vs_Full.tsv` | 文件名为 TSV | 实为 CSV（逗号分隔）。已重新导出真 TSV（本轮目录），旧文件原位保留供追溯。GSE143531 表述本身不变："健康海马 PAP 数据中未检出相对 Full 的差异，估计 log2FC −0.273，FDR 0.683"（维持原表述，不得用于证明无富集） | 本轮 edgeR_TMM_PAP_vs_Full.tsv | 2026-09-24 |
| 7 | GSE74456 的 2,729 基因列表/Atp2a2 log2FC 2.13 | 若用于论文需交原始产物 | 维持既定要求：须交付产生该数值的原表、脚本、设计矩阵与样本单位，或清楚标为"未独立核验的历史结果"。本轮未重算该数据集 | paper_gene_set_provenance（后续待办） | 2026-09-24 |
| 8 | 两报告的"Top5 候选 / 定靶"段 | 基于 v5.1 错侧区间的候选排序（Ndrg2≈Agpat3>Cnp≈Pea15a>Fam107a）与后续定靶并存易误读 | 明确标注：该排序为**链修正前的历史结果**；链修正后当前主靶 = Atp2a2（B 级），备靶 Agpat3、Aplp1（见 gate1_decision.md 与 target_nomination_v1.md） | gate1_decision.md；target_nomination_v1.md | 2026-09-24 |

## 不受本轮影响、维持原状的段落

- P1 调控子层全部结果（31 调控子、Pabpn1/Qk）——未重算；
- P2 事件层与 SAP 统计（9,693/977）——未重算；
- GSE143531 库水平 edgeR TMM 结论——本轮仅修文件格式，未改数值；
- PolyASite 距离、保守性、motif 富集（阴性）——未重算；
- 中文结论与 Gate1 决定见 `gate1_decision.md`、`中文结论.md`。
