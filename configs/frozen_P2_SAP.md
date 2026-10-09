# P2 统计分析计划（SAP）v1.0 —— 全量 12 样本统计前冻结

> 依据：2026-09-19 论证逻辑审计 §3/§4/§6（reviews/2026-09-18_论证逻辑审计报告_20260918.md），用户采纳。
> 冻结时点：2026-09-19，**早于全量 DaPars2 完成与任何全量统计**（试点 4 样本描述性结果已知，全量结果未看——非 outcome-informed）。
> 铁律：本文件 commit 后方可执行层 2 统计；任何修改须显式登记理由与时间。
>
> **修订 v1.1（2026-09-19 晚，用户批准"用R"）**：层 2 的 moderated linear model 改用 **R limma 原生实现**（`lmFit` + `eBayes(robust=TRUE, trend=TRUE)`，与审计 §3.3 建议逐字对齐），替代原"纯 Python 重实现 eBayes"方案（该方案连同其回验要求一并作废，不执行）。理由：零实现风险、分钟级运行；属于 P0-07 原决策预留的"需要时装 R"分支。R 仅进入层 2 统计检验这一层，定量（DaPars2）与其余全部数据处理仍为 Python。环境 = WSL conda env `rlimma`（R 4.3.3 + limma 3.58.1，已冒烟验证）。

## 0. 设计事实与定位

- GSE238125：皮层 astrocyte bulk RNA-seq，Sham×2 + d1/d3/d7/d21/d60 各×2；**每组 n=2，全项目定位为 discovery 层，不做 confirmatory 推断**。
- 5 个 contrast 共享 sham1/sham2——跨时间点结果高度相关，**"≥2 时间点同向"仅为一致性标注，不是独立复现**。
- DaPars2 只负责事件定量（联合 breakpoint 推断 + PDUI），差异检验由本 SAP 定义。

## 1. 层 1：样本级全局 APA（QC + 全局位移）

每个样本计算：median PDUI、10% trimmed mean PDUI、相对 sham 均值的 ΔPDUI<-0.1 事件比例、ΔPDUI>0.1 事件比例、缩短/延长比。按 **基因表达量分位（4 档）、事件 coverage、UTR 长度分位** 分层重复上述指标。
样本级 QC：12×12 PDUI 相关热图 + PCA；每样本 PDUI 分布。
技术混杂排除（d1 缩短进正文的必要条件）：gene-body coverage（RSeQC 或 bedgraph 等价口径）、3′ bias、insert size/duplication/mapping 率（从 align + fastp 日志提取）；ΔPDUI 与基因表达 logFC、平均 coverage、UTR 长度的相关检查。
敏感性：leave-one-sample-out 全局比例重算。
产出：`results/p2_full_qc_sample_level.tsv/.pdf` + 分层表。**层 2 仅在层 1 无系统性技术伪迹时启动。**

## 2. 层 2：事件级探索性模型（moderated，纯 Python）

1. 事件进入某 contrast 的前提：**sham1+sham2+该时间点 rep1+rep2 共 4/4 有效 PDUI**（缺失 → `not_testable`，不插补）。量化层（DaPars2 ≥30% 规则）只决定事件是否被输出，不作为检验资格。
2. PDUI 边界修正 logit 变换：`logit((PDUI*(1-2ε)+ε))`，ε=0.01（预冻结）。
3. Moderated linear model：`~ 0 + timepoint`；方差借全部事件的经验贝叶斯收缩（limma-eBayes/Smyth 2004 算法纯 Python 重实现，robust+trend 语义；实现须对已发表基准数据回验，通过标准 = 与 limma 参考输出相关 r>0.99）。**维持免 R 红线**；若回验失败 → ⚠ 报批 R/limma 例外后再跑。
4. 检验顺序：先 omnibus time effect（联合 F，5 水平），再 5 个预设 contrast（d1/d3/d7/d21/d60 各 vs sham）。
5. 多重检验（预冻结）：**主口径 = 全部 event × 5 contrast 联合 BH，FDR<0.1**；各 contrast 内 BH 同步记录为辅助口径。
6. 效应阈值（预冻结）：**|ΔPDUI| ≥ 0.1**（ΔPDUI = 时间点均值 − sham 均值；负 = 缩短）。
7. 显著事件定义（主口径）= 联合 BH FDR<0.1 **且** |ΔPDUI|≥0.1 **且** 该 contrast 两重复方向一致。方向一致性标注沿用已批规则（≥2 时间点同向），明示非独立复现。

## 3. Gate G2 组合门槛（对照 TODO P2-G1）

(a) 层 1 QC 通过；(b) 全局位移分层后方向稳定 + LOSO 敏感；(c) 主口径显著事件满足上述三条件；(d) 显著事件数记录在案（原预注册 ≥100 为记录项）+ top 事件最低复核抽验：PolyASite mm10 PAS 距离注释 + 第二 APA 工具（QAPA）方向一致性抽验（≥10 事件或 10%，取大）。
任何一层失败 → 如实降级探索性叙事（kill-switch B 精神保留）；启用 PB-2 验证集前须预固定其规则（防 dataset shopping）。

## 4. 事件表述边界

全部事件表述为 "DaPars2-inferred proximal cleavage shift within the selected longest annotated terminal 3′UTR reference"；不笼统写"该基因发生 APA"。44.8% 试点产出率归因 = coverage threshold + 跨样本覆盖过滤 + 所选注释共同作用。

## 5. 产出清单

| 文件 | 内容 |
|---|---|
| results/p2_full_pdui_matrix.tsv | 全量 PDUI 矩阵（12 样本 × 事件） |
| results/p2_full_qc_sample_level.* | 层 1 全部 QC |
| results/p2_full_omnibus_time.tsv | omnibus 检验 |
| results/p2_apa_events_planB.tsv | 5 contrast 主口径事件表（含一致性标注） |
| results/p2_gate_g2_evidence.md | G2 证据包（组合门槛逐条） |
