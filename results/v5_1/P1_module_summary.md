# P1 模块小结 — 2026-09-16（Gate G1 证据包）

> ⚠ **2026-09-16 深夜勘误**：P1-05 的 GSE286075 配对 DE 原取 ReadsPerGene 错误链列（col3），实为 dUTP 第二链特异库（实证 scripts/diag_286_strand.py：col4≈col2≈12.5–29.5M reads、col3≈0）。**下方"P1-04/05 DE 结果"与"Gate G1 判定"两节中关于 L2 层的数字与"深度极低"结论已作废**，修正值见文末「勘误与修正重跑」节；P1-02/03/04（QC/聚类/L1 层 DE）不受影响。

## P1-02 QC 结果（冻结参数：genes≥200, counts≥500, mt<20%, MAD±3）

| 样本 | 总细胞 | 过 QC | 备注 |
|------|--------|-------|------|
| MCAO1/2/3 | 11,772/11,361/8,104 | 11,687/11,285/8,052 | MAD 过滤较松（MCAO 组 UMI 更高） |
| sham1/2/3 | 8,771/8,540/9,980 | 8,280/8,141/9,459 | |
| 合计 | 58,528 | **56,904** | 通过率 97.2% |

- Scrublet 依赖在本机安装失败（构建链缺 setuptools 隔离）→ 按冻结文件预登记降级为 MAD 离群过滤，偏差已记录。30 簇内无明显双胞专簇（无极小高复杂度异常簇），风险可控。

## P1-03 astrocyte 判定（30 个 Leiden 簇，分辨率 1.0）

- 按冻结规则命中：**簇 5（2,815 cells: Aqp4 1.38 / Slc1a3 2.95）+ 簇 21（230 cells: Gfap^high 1.39, Serpina3n 0.84）**，合计 3,045 astrocyte。
- 无簇被免疫/内皮标志排除（Aqp4+Slc1a3 高的簇 C1qa/Cx3cr1/Pecam1 均低于阈值）。
- 产出：results/p1_cluster_markers.tsv（30×13 marker 均值表）、p1_all_cells_clusters.tsv、两张 UAP 图。

## astrocyte 亚群（reactive_score = Serpina3n+C3+Gfap，阈值 0.5）

| 样本 | homeostatic | reactive |
|------|-------------|----------|
| MCAO1/2/3 | 529/511/362 | **216/194/157** |
| sham1/2/3 | 326/331/397 | 5/11/6 |

- MCAO 组 reactive 显著扩增（约 1/3 vs sham ~1.5%），与 P0-08 快查的 Gfap/Serpina3n 表达率跳升一致。
- **按冻结规则（单元 <20 细胞不入 DE）**：sham reactive 3 个单元（5/11/6 cells）排除 → reactive 口径 DE 无效（case 侧 3 vs control 侧 0），P1 DE 以 overall 口径为主证据，homeostatic 口径为辅（529/511/362 vs 326/331/397）。

## P1-04/05 DE 结果

- overall astrocyte（n=3+3 pseudo-bulk）：14,252 基因入模，**1,280 个 padj<0.05**——signal 充足。
- homeostatic 亚群 DE：13,837 基因入模（辅助口径）。
- GSE286075 配对 DE（~mouse+condition）：3,688 基因入模（≥10 reads 过滤后），**padj<0.05 仅 6 个**——深层原因：RiboTag IP 文库测序深度极低（Stroke1 总 reads ~9 万，约为常规 bulk 的 1/300），功效有限。6 个显著基因：Gm20481↑、Gm38287↓、Ccdc85b↑ 等，无 APA regulator。**解读纪律**：该层低功效不等于"regulator 未变"，方向一致性仍可用作弱证据（padj 全部≈1 的层不单独下结论）。

## P1-06 regulator 目录（31 个，Cpsf4l 未检出）

**两层方向一致的 4 个（G1 判定核心）**：

| gene | L1 astrocyte log2FC (padj) | L2 paired log2FC | 一致方向 | 生物学含义 |
|------|---------------------------|------------------|----------|-----------|
| **Qk** | −0.266 (0.077) | −0.826 | ↓↓ | RBP 调制型；边缘显著 |
| **Pabpn1** | **−0.630 (3e-4)** | −（L2 NA） | L1 显著↓ | 核心机器；唯一 padj<0.05 |
| Cstf3 | −0.246 (0.69) | −1.192 | ↓↓ | 核心机器 |
| Srsf1 | −0.246 (0.39) | −1.917 | ↓↓ | 剪接型调制 |
| Ptbp1 | +0.055 (0.94) | +2.978 | ↑↑ | 方向一致但 L1 无信号 |

- L2 层 Qk/Ptbp1/Srsf1/Cstf3 的 padj 全部 ≈1.0（深度限制），"concordant" 仅指方向；Pabpn1 是唯一的强信号（两层设计下 astrocyte L1 padj=0.0003，homeostatic 口径 padj=0.017 重复验证）。
- 快查表达率证据（P0-08）：Nudt21/Cpsf6/Cstf2/Elavl1 表达率在 MCAO 组上升 +4–6.5pp，Qk 下降 8.6pp——与 DE 的 Qk↓ 方向一致；Nudt21/Cpsf6 的 log2FC 在 L1 接近 0（−0.06/−0.07），表达率上升可能由 astrocyte 组成比例变化驱动（reactive 扩增），单看 regulator 分子层面的变化以 Pabpn1↓、Qk↓ 最实。

## Gate G1 判定（照 v5.1 规则）

- 规则要求：≥3 个 regulator 两层证据方向一致 → **形式上满足（4 个 concordant）**，但需如实说明：L2 层 padj 均无显著性（低深度），方向的先验价值有限；最强的独立信号是 Pabpn1↓（padj=3e-4）与 Qk↓（padj=0.077 边缘）。
- 对 M3 的影响：regulator 层证据强度定级为「弱-中」；M3 的 L1 文献池应重点收录 PABPN1 与 QKI 的 perturbation 证据（两者均有经典文献：PABPN1 KD→3'UTR 缩短；QKI 敲除→星形胶质细胞成熟与 3'UTR 结合谱改变）。
- 对 M2 的影响：不触发降级（M2 独立选源，事件数量判定在 Gate G2）。

## 产物清单（已登记 inventory）

p1_qc_summary.tsv / p1_cluster_markers.tsv / p1_all_cells_clusters.tsv / p1_astrocyte_subpop_counts.tsv / p1_pseudobulk_astro_counts.tsv+meta / p1_DE_overall_astro.tsv / p1_DE_homeostatic_astro.tsv / p1_DE_paired_GSE286075.tsv（含 symbol） / p1_regulator_catalog.tsv / p1_regulator_quickcheck.tsv(+Qk) / p1_astrocyte_clusters/*.png

## 勘误与修正重跑（2026-09-16 23:58，L2 链取列修正）

**根因**：p1_de_catalog.py 原 GSE286075 段取 `int(r[2])`（STAR ReadsPerGene col3=第一链计数），注释误写 "col3=N_unstranded"。diag_286_strand.py 实证：该库为 dUTP 第二链特异（col4≈col2，col3 仅 0.4–8%），"Stroke1 ~9 万 reads" 正是 Stroke1 col3 的值（90,093），真实深度为 col2/col4 ≈ 12.5–29.5M reads/样本——完全正常。

**修正动作**：仅改取列 r[2]→r[3]（design=~mouse+condition 等冻结参数未动），重跑 p1_de_catalog.py 全脚本；L1 两表逐位复现（overall 1,280 padj<0.05 不变），L2 得到修正值。

**修正后 L2（GSE286075 配对层）**：14,323 基因入模（≥10 reads），padj<0.05 = **409**，padj<0.1 = 576。该层由"方向性弱证据（原 6 个 sig）"升级为正常功效配对 DE 层。

**修正后 regulator 目录（P1-06）**：两层方向一致 **17/31**（原 4）。仍一致的：Qk、Ptbp1、Wdr33、Nudt21、Cpsf1、Cpsf3、Cpsf4、Fip1l1、Cstf2、Clp1、Elavl1、Nova1、Nova2、Fus、Tardbp、Hnrnpa2b1、Mbnl2；**Cstf3、Srsf1 翻转为 discordant**（L2 方向变号，如实登记）；Cpsf4l 在 L2 可检出（log2FC 2.34，padj NA）。

**Gate G1 复核**：结论不变且增强（17 ≥ 3），无需重开 Gate。regulator 层证据强度定级由「弱-中」升为「中」。最实信号：**Pabpn1↓**（L1 padj 3e-4，homeostatic 0.017 重复验证）与 **Qk 双层↓**（L1 padj 0.077 / L2 padj 0.122，log2FC −0.266/−0.558）。M3 L1 文献池优先级维持 PABPN1 + QKI。

**存档（负结果照登）**：错误版结果 = p1_DE_paired_GSE286075.wrongstrand.bak.tsv、p1_regulator_catalog.wrongstrand.bak.tsv（已登记 inventory）。
