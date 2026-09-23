# TODO — 卒中→APA→Endfoot 定位潜能 Pipeline v5.1 执行清单

> **⛔ 2026-09-21 起本清单转入历史档案**：现行执行文件 = `D:\stroke_apa_reanalysis_项目TODO.md` + `01_研究Proposal.md`（用户定稿，含负链方向重建计划 P0–P10）。本文件记录 v5.1 干实验阶段（P0–P6）全部执行记录，不再更新。

> 依据：`STROKE_APA_PIPELINE_v5.1.md`（2026-09-16 定版，外部审阅意见已落实）。
> 制定于 2026-09-16。文件位置：`D:\stroke_apa\TODO_pipeline_v5.1.md`（唯一主清单，完成进度只记在这里）。
> **2026-09-19 论证逻辑审计修订（v5.1r）**：外部审计报告已采纳（`reviews/2026-09-18_论证逻辑审计报告_20260918.md`），总体判定 = **Major revision / 有条件 Go**（P2 全量继续，论证结构修正）。要点：因果链表述降级为方向一致性三角互证；新增 P2-06a 统计分析计划（SAP）冻结；Gate G2 换组合门槛（修订发生在全量统计之前，显式登记）；M4 拆三独立注释轴；MPRA 阻断待 manifest；事实修正若干。总览见 `STROKE_APA_PIPELINE_v5.1.md` §0b。
> 干分析课题：无湿实验阻断项；湿实验仅出现在 P5 路线书（文字交付，不执行）。

## 填写规范（照此执行，不留"mental note"）

1. 完成一项：`[ ]` → `[x]`，并在该项末尾【记录】内填写：**`YYYY-MM-DD HH:MM ✓` + 客观证据**（文件字节数/行数/基因数/统计量/退出码/登记表行号），有问题与修正一并写明。判定为不执行/未触发的项：`[ ]` → `[-]`，并在【记录】写明触发失效的理由（如 Gate GA 改道）。
2. 涉及方案变更或二选一决策的项（标注 ⚠）：**先写明候选方案与推荐理由，用户批准后方可执行**；批准时间与内容记入【记录】。
3. 每个 Gate 区块：前置项全部勾选后停下，**向用户逐条出示证据，用户确认后才进入下一阶段**；确认语记入 Gate 的【记录】。
4. 一切运行产物（下载文件/中间表/结果表/图）登记进 `results/data_inventory.tsv`（列：日期/文件或表名/来源/大小或行列数/用途/所在阶段），【记录】中写明登记行号。
5. 表述纪律相关项：执行前先查本清单"附 A 纪律红线"，违反红线的输出即使跑完也无效。

进度总览：`P0 [x] / P1 [x] / P2 [x] / P3 [x] / P4 [x] / P5 [x]`——**2026-09-20 全部阶段关闭（G1–G5 全过）**；09-20 收尾检查补勾漏勾项（工作已完成但漏勾 13 项；P2-07 未触发 / P4-04 阻断按规范记 `[-]`）

### 关键路径与目标日期（2026-09-16 23:17 依实际进度排定；**2026-09-18 21:35 修订**：数据已 100% 到齐，比对 12/12 完成，试点 DaPars2 已出结果；**2026-09-18 晚二次修订：用户决策把剩余重活（全量 DaPars2 / M4 / P3–P4 计算）打包转集群执行**，交接文档 = `WORK_LOG_20260918.md`）

| 节点 | 目标日期 | 状态 |
|------|----------|------|
| P2 试点闭环（4 runs md5 → 比对 → DaPars2 首批 PDUI → 质检） | 2026-09-18 | ✅ **09-18 提前完成**（9,269 事件 / sham r=0.919 / d1 缩短 3.4:1，见 P2-06 记录） |
| P2 全量事件表 + Gate G2 | **2026-09-19**（提前 5 天） | ✅ 09-20 G2 通过（977 事件 / 组合门槛 a–d 全过） |
| P3 完成 + Gate G3 | **2026-09-24** | ✅ 09-20 G3 通过（三轴分别表述：PAP 阳性 / endfoot 阴性 / WM 弱同向） |
| P4 完成 + Gate G4 | 2026-10-07 | ✅ 09-20 G4 通过（B=5/C=475/D=250；MPRA 阻断照登） |
| P5 交付 + Gate G5 | 2026-10-12 | ✅ 09-20 G5 通过，课题干实验阶段关闭 |

---

## 阶段 P0 — 数据底座与启动审计（1–2 天）

### 工作包 1：数据获取与事实基线

- [x] **P0-01 建库类型判定（T0.2）**：确认 GSE174574 建库类型，锁定 M2 工具【2026-09-16 15:38 ✓ **10x Genomics 3' droplet**：原文 Zheng 2022 JCBFM（PMID 34496660/PMC8721774）正文 "performed single-cell transcriptomic analyses with the 10x genomics"；SRA SRP320164 6/6 run = PAIRED + LIBRARY_SELECTION=cDNA + Illumina HiSeq 4000（10x 特征组合），无 smart-seq/SMARTer 关键词 → M2 主工具锁定 scAPAtrap/Sierra（3'-tag 感知），DaPars2/QAPA 仅用于 bulk Plan B；模型=MCAO 24h 整侧半球，C57BL/6 6–8 周，n=3+3，17 clusters】
- [x] **P0-02 五个主线数据集 GEO 核实**：登录号/设计/样本数/补充文件形态逐一登记【2026-09-16 ✓ GSE174574（scRNA，SRP320164+processed mtx）、GSE286075（endfoot RiboTag，仅 gene counts×6——本地无 FASTQ，SRA 有原始 reads，本阶段不用→禁用于 APA，2026-09-19 措辞修正）、GSE74456（PAPTRAP，SRP065508）、GSE146935（QKI CLIP，SRP252704）、GSE330741（SN-MPRA，PRJNA1465359）；详见 STROKE_APA_PIPELINE_v5.1.md §2.1】
- [x] **P0-03 Plan B 池核实**：4 套备选数据集存在性与原始数据可用性【2026-09-16 ✓ PB-1 GSE238125（皮层 astrocyte bulk，Sham+d1/d3/d7/d21/d60 n=2/点，PRJNA997998=12 runs，抽验 PAIRED/cDNA 全长）；PB-2 GSE225110（RiboTag astro IP+input 4h/3d，PRJNA922165=132 runs）；PB-3 GSE263986（mGFAP-RiboTag LCM 分区 zone1-5，PRJNA1100502=82 runs，另有 ProcessedData.xlsx 15MB）；PB-4 GSE165384（海马 bulk 2h，SRP302969）——均有 SRA】
- [x] **P0-04 下载 GSE174574_RAW.tar + 解包校验**【2026-09-16 16:55 ✓ 并行 Range 分块下载（16×8MB×重试），263,086,080 B 与 filelist 标称精确一致；tar -tf = 18/18 成员；解包至 data\GSE174574\extract\（6 样本 × barcodes/genes/matrix.mtx 三件套，mtx 35–54MB gz）；登记 inventory】
- [x] **P0-05 下载 GSE286075_RAW.tar + 解包校验**【2026-09-16 16:57 ✓ 1,802,240 B 精确一致；解包 6 个 ReadsPerGene.out.tab.gz（CTR1-3/Stroke1-3）；登记 inventory】

### 工作包 2：工作区与登记表

- [x] **P0-06 建 results/ 与 data_inventory.tsv**：按填写规范第 4 条建表，补登 P0-01~05 已产出物【2026-09-16 17:31 ✓ 表建 7 行（含 h5ad 与两张 audit 表）；后续产物逐行追加】
- [x] **P0-07 环境与依赖决策（⚠）**：Python 3.13.12 / git 2.55.0 已确认；R 不在 PATH【2026-09-16 17:06 ✓ 用户指示"继续吧，完成P1"→ 采纳方案 (b)：scanpy 1.12.4 + pydeseq2 0.5.4 + leidenalg 安装成功；M1 全程免 R；R/scAPAtrap 决策推迟到 Gate GA 后，若主源成立再装 R】

### 工作包 3：启动审计（审阅意见四件事的执行）

- [x] **P0-08 T0.5 regulator 快查**：读 6 样本 mtx → 聚合基因 × 样本计数 → 5 核心 regulator + 3 astrocyte marker 表达矩阵【2026-09-16 17:38 ✓ 11 基因 × 6 样本（results/p0_regulator_quickcheck.tsv + p0_regulator_quickcheck_Qk.tsv，60+18 行）。**Qk 首轮 miss（MGI 官方符号 Qk 非 Qki），补查已入表**。关键读数：① astrocyte marker 全组可检出（Slc1a3 表达率 18–23%，Aqp4 3.5–8.3%）；② 反应性标志 Gfap/Serpina3n 在 MCAO 组剧烈上调（Gfap 表达率 1.5–1.8%→11.9–12.5%，Serpina3n 0.4–0.5%→9.3–11.0%）→ 反应性星形胶质细胞分群可行；③ 五核心 regulator 全部可检出且呈系统性变化趋势：Nudt21/Cpsf6/Cstf2/Elavl1 上调（表达率 +4–6.5pp）、**Qk 下调**（42.8%→34.2%）——方向与假说兼容，详细定量待 P1-04 DE】
- [x] **P0-09 T0.6 每样本细胞数审计**：barcodes 行数 + QC 粗过滤 → marker 注释【2026-09-16 17:06 ✓ 原始 QC 口径完成：sham1/2/3 = 8753/8512/9951 过 QC（总 27,291），MCAO1/2/3 = 11719/11322/8083（总 31,124），QC 通过率 99.4–99.6%，中位 mt ~2.4%；astrocyte-likely 精确数待 P1-03 聚类注释；表 = results/p0_cell_audit_raw.tsv】【2026-09-16 23:17 ✓ 收口：P1-03 实测 astrocyte 合计 3,045（簇5 2,815 + 簇21 230；6 样本均值 ~507/样本，处 ≥500 阈值临界区）；该不确定性已被 Gate GA 决策吸收（M2 切 Plan B）——本项闭环】
- [-] **P0-10 T0.7（可选，触发式）PB-1 SRA 预下载**：仅当 P0-09 显示 astrocyte 细胞数临界（500±200/样本）时启动；12 runs 体量约 5–15GB【记录：未触发——2026-09-16 18:02 Gate GA 直接判定 M2 走 Plan B，PB-1 下载由 P2-05 承接，本项预下载条件失效，关闭】

### ⛔ Gate GA（APA 主源判定，用户确认后进 P1+P2）
- [x] **P0-G1** 主源判定【2026-09-16 18:02 ✓ 用户审阅 P1 证据包后决策：M2 用 Plan B（GSE238125 bulk + DaPars2 纯 Python），不装 R 不跑 10x APA；GSE174574 仍承担 M1 regulator DE 主源】
- [x] **P0-G2** 工具路线决策【2026-09-16 18:02 ✓ 纯 Python 全案成立（pydeseq2 + PDUI 计算）；R/scAPAtrap 不再需要】
- [x] **P0-G3** 【2026-09-16 18:02 用户确认记录："进P2，M2用planB"；Gate GA 通过，P0 关闭】

### ✅ P0 完成记录（2026-09-16 18:02 闭环）
- 时间：2026-09-16 16:57 → 18:02
- 说明：主线数据两套下载解包完成；纯 Python 环境就绪；细胞审计与 regulator 快查完成；Gate GA 用户决策 = M2 走 Plan B；GSE238125 SRA 启动获取。

---

## 阶段 P1 — M1 Regulator 层（约 1 周）

> 【2026-09-19 审计修订】① **GSE286075 层重定位** = 卒中后血管周 endfoot 转译组响应参考；其 regulator mRNA 变化**不再作为核内 APA regulator 活性的独立验证层**（剪切/多聚腺苷酸化发生在核，endfoot 转译组不承载该推理）；"17/31 两层方向一致"降级为**描述性方向 concordance**（非验证）；Pabpn1 = L1 显著的上游候选（非已验证 driver）；Qk = 两层方向一致但**均未达预设显著**（L1 padj 0.077 / L2 0.122）。② M1 整体定位改为 "concurrent upstream candidates（并发上游候选）"。③ 新增 P1-07/P1-08。

### 工作包 4：GSE174574 处理与 astrocyte 定义

- [x] **P1-01 ⚠ astrocyte 定义与 QC 参数冻结**【2026-09-16 17:10 ✓ 冻结文件 results/P1_frozen_params.md（QC 阈值/30 簇规则/astro 簇判定/亚群评分/DE 设计全部预写）；git init -b main + commit c75afc1（先于任何 DE）；注：git status 显示仓库仅收 TODO/pipeline/frozen/results-inventory/scripts，原始数据不进 git】
- [x] **P1-02 QC + 去双胞**【2026-09-16 17:46 ✓ Scrublet 安装失败（构建缺 setuptools 隔离）→ 按冻结文件预登记降级为 MAD±3 离群过滤，偏差已记录；56,904/58,528 过 QC（97.2%），每样本明细见 results/p1_qc_summary.tsv】
- [x] **P1-03 聚类 + astrocyte 簇定义**【2026-09-16 17:46 ✓ HVG2000→PCA50→leiden res1.0 = 30 簇；冻结规则命中簇 5（2,815 cells, Aqp4 1.38/Slc1a3 2.95）+ 簇 21（230 cells, Gfap^high/Serpina3n^high），合计 3,045 astrocyte，无簇被免疫/内皮标志排除；reactive 评分：MCAO 216/194/157 vs sham 5/11/6（reactive 显著扩增）；sham reactive 单元 <20 细胞按冻结规则排除出 DE 并登记；产出 UMAP×2 + marker 表，图 results/p1_astrocyte_clusters/】
- [x] **P1-04 astrocyte pseudo-bulk DE（主证据层）**【2026-09-16 17:5x ✓ pydeseq2 overall 口径（n=3+3）：14,252 基因入模，1,280 padj<0.05；homeostatic 口径辅助（13,837 基因）；表 = results/p1_DE_overall_astro.tsv + p1_DE_homeostatic_astro.tsv】
- [x] **P1-05 GSE286075 配对 DE（第二证据层）**【2026-09-16 ✓ pydeseq2 ~mouse+condition：3,688 基因入模，padj<0.05 仅 6 个（Gm20481↑/Gm38287↓/Ccdc85b↑等，无 regulator）——根因=RiboTag IP 测序深度极低（Stroke1 总 reads ~9 万），功效有限；解读纪律：该层仅方向性弱证据，不单独下结论；表已回填 symbol 列】【2026-09-16 23:58 ⚠ 勘误+修正重跑：上述"深度极低"系取列错误——原取 ReadsPerGene col3（第一链），实查（scripts/diag_286_strand.py）该库为 dUTP 第二链特异：col4≈col2（12.5–29.5M reads/样本）、col3≈0（Stroke1 col3=90,093 正是"~9万"的来源）。改取 col4 重跑：**14,323 基因入模，padj<0.05 = 409**（padj<0.1 = 576），该层由"方向性弱证据"升级为正常功效配对 DE 层；design/冻结参数未动，仅输入列修正。错误版备份 = results/p1_DE_paired_GSE286075.wrongstrand.bak.tsv（负结果照登）】
- [x] **P1-06 regulator 名单分级目录**【2026-09-16 ✓ 31 个 regulator（Cpsf4l 未检出），4 个两层方向一致（Qk↓↓/Cstf3↓↓/Srsf1↓↓/Ptbp1↑↑），最强独立信号 Pabpn1↓（L1 padj=3e-4，homeostatic 口径 0.017 重复验证）；表 = results/p1_regulator_catalog.tsv；模块小结 = results/P1_module_summary.md】【2026-09-16 23:58 勘误联动（P1-05 链取列修正后重跑）：两层方向一致 = **17/31**（旧 4）——Qk/Ptbp1/Wdr33/Nudt21/Cpsf1/Cpsf3/Cpsf4/Fip1l1/Cstf2/Clp1/Elavl1/Nova1/Nova2/Fus/Tardbp/Hnrnpa2b1/Mbnl2；**Cstf3、Srsf1 翻转为 discordant**（L2 方向变号），如实登记；Cpsf4l 在 L2 可检出（log2FC 2.34，padj NA）。L2 最强 regulator 信号 = Qk（padj 0.122，↓−0.558），与 L1↓ 一致；Pabpn1 L1 强信号不变。目录表已重写，旧版备份 = p1_regulator_catalog.wrongstrand.bak.tsv】
- [x] **P1-07 doublet/ambient 敏感性（2026-09-19 审计新增）**：双联体检测（scanpy 原生 sc.pp.scrublet，scikit-image 0.26.0 已装）+ ambient proxy 评分（GEO 仅 filtered mtx 无 raw matrix——ambient 用低计数池 signature 近似，限制如实登记）；变体 B=QC−双联体、C=B−ambient-high，冻结参数全链重跑对比【2026-09-19 ✓ **完成，P1 结论稳健**：scrublet 双联体 1,624（2.85%），ambient 标记 2,849（5.0%）；astrocyte 数 B=3,046 / C=3,042 vs 原 3,045（几乎不动）；DE_sig B=1,310 / C=1,311 vs 原 1,280（overlap 96%/94%）；**31 regulator 方向翻转 = 0**——双联体/ambient 不驱动任何 P1 结论；表 = results/p1_07_sensitivity_summary.tsv + p1_07_variant_*；踩坑：此版 scanpy 列名 = predicted_doublet（非 is_doublet），scrublet 白算一轮】【记录：敏感性通过，0 翻转】
- [x] **P1-08 GSE286075 205 vs 409 复现附录（2026-09-19 审计新增）**：对照原论文 design 公式/paired 口径/链特异性取列/过滤/阈值，逐项列差异来源【2026-09-19 ✓ **完成**：原文方法实取自 PMC11720067 全文（DESeq2 v1.44.0，预过滤 ≥10 counts in ≥3 samples → 10,862 基因——本项目同规则**精确复现 10,862**；ashr s-value 区间判据 |log2FC|>0.322）；网格（scripts/p1_08_205v409.py）：plain 391 / 区间判据 55（34↑21↓，**54/55 在我们 409 内，方向零矛盾**）/ col2 链敏感性 370 vs 391；ashr s 阈值原文未说明 → 精确复现 205 原则上不可达，**409 vs 205 = 判据差异非管线错误**；anchor：Hspa1a/Hspa1b 两套口径均排上调第 1/2 位（=原文点名 top hits）；附录 = results/p1_08_286075_205v409_appendix.md；409 仅作描述性方向参考/富集集【记录：附录 v1 入库】

### ⛔ Gate G1
- [x] **P1-G1** ≥3 个 regulator 两层证据方向一致【2026-09-16 ✓ 形式上满足：4 个 concordant（Qk/Cstf3/Srsf1/Ptbp1）；但需用户知情：L2 配对层 padj 均≈1（RiboTag 深度 ~9万 reads，功效极低），concordant 仅指方向；最强独立信号 = Pabpn1↓（padj 3e-4）+ Qk↓（0.077 边缘）；regulator 层证据强度定级「弱-中」，M3 的 L1 文献池应重点收 PABPN1/QKI perturbation 证据】【2026-09-16 23:58 勘误联动：取列修正后 G1 **结论不变且增强**——concordant = 17/31（≥3 大幅满足）；L2 不再"全 padj≈1"（409 padj<0.05）；双层最实信号 = Pabpn1↓（L1 3e-4）+ Qk 双层↓（L1 0.077 / L2 0.122）；证据强度定级升为「中」。Gate 无需重开，P1 关闭状态不变】
- [x] **P1-G2** 若不满足：M3 降级决策（⚠ 用户批准改纯 CLIP/motif 途径）【记录：P1-G1 已满足，未触发】
- [x] **P1-G3** 【2026-09-16 18:02 用户确认记录："进P2，M2用planB"；Gate G1 通过，P1 关闭】

### ✅ P1 完成记录（2026-09-16 18:02 闭环）
- 时间：2026-09-16 17:06 → 17:58（约 52 min）
- 说明：QC 97.2% 通过（scrublet 降级 MAD 已登记）；30 簇，astrocyte = 簇5+21 共 3,045 细胞；reactive 亚群 MCAO 显著扩增；overall DE 1,280 sig；配对层低功效（深度限制）如实降级为方向性弱证据；regulator 目录 31 个，4 concordant + Pabpn1 强信号。全部证据与解读纪律见 results/P1_module_summary.md。

---

## 阶段 P2 — M2 APA 事件层（约 1 周）★ Gate GA 已判 Plan B：工作包 5A 关闭，执行 5B

### 工作包 5A：主线（GSE174574）——⛔ 已由 Gate GA 关闭，不执行

- [-] **P2-01 scAPAtrap/Sierra 环境搭建**【记录：未触发——Gate GA 2026-09-16 18:02 判定 M2 走 Plan B（见 P0-G1），GSE174574 为 10x 3′ 文库且用户决策不装 R，本包永久关闭】
- [-] **P2-02 单样本试跑**【记录：未触发——同上】
- [-] **P2-03 全量 APA 定量**【记录：未触发——同上】
- [-] **P2-04 事件筛选**【记录：未触发——同上；对应 Plan B 项为 P2-06】

### 工作包 5B：Plan B（GSE238125 bulk，仅 Gate GA 切换时）

- [x] **P2-05 GSE238125 SRA 下载+比对**：12 runs（PRJNA997998）→ 比对/定量【2026-09-20 ✓ **闭环确认（收尾补勾）**：子项 05a–05i 全部完成——12/12 下载到齐、md5 24/24 通过（05c）、12/12 比对完成唯一率 79–96.4%（05i）、bedgraph 12 份就位；冷备 D:\stroke_apa_data\ext4_backup_20260918\；本项父任务随子项全勾而关闭】【2026-09-16 18:03 启动；21:45 侦查更新：① ENA filereport 全量拉取（D:\stroke_apa_data\ena_filereport.tsv）：12/12 runs 有双端 FASTQ，合计 47.66 GB（3.2–4.8 GB/run，19.9–26.8 M reads/run，GSM7658905–916 = Sham×2+d1/d3/d7/d21/d60×2）；② 实测速率：ENA ftp.sra.ebi.ac.uk 并行 8 线程仅 ~0.18 MB/s（sham1_1 R1 已完成 2,090,117,617 B，sham1_2 R2 分块至 ~264MB 后停），全量不可行；改测 NCBI S3 镜像 sra-pub-run-odp.s3.amazonaws.com/sra/<SRR>/<SRR>（200 OK，3.28 GB/run .sra，Accept-Ranges=bytes）；③ 工具盘点：Windows 无比对器/sra-tools；WSL2 Ubuntu 24.04 自带 STAR 2.7.11b + samtools + bedtools + python3/pip3 + gcc（12 核 / 23GB RAM / D 盘 299GB 空闲）→ 比对定量走 WSL，免 R 红线不破；④ 用户决策（21:24 "在D盘进行任务" + 21:41 确认按方案 A 执行）：数据目录 = D:\stroke_apa_data\（fastq/bam/index/tools 已建），两步走 = 先 4 runs 试跑闭环（Sham rep1/2 + day1 rep1/2）→ 1 run 比对验证 → DaPars2 首批 PDUI → 再拉剩余 8 runs；⑤ 偏差登记：下载根目录由 D:\stroke_apa\data\GSE238125\ 改为 D:\stroke_apa_data\（用户指示）；ENA 半成品进度如切换源则废弃并记录；⑥ 重大发现：WSL ~/stroke_APA_pilot 前期试点工作区（9月9-15日，GSE225110 6 runs 全链 QAPA+DaPars2+reanalysis_v3 纠错报告）提供全套可复用资产：GRCm38 genome.fa + GENCODE vM25 GTF + STAR index（24GB，建在 23GB RAM 约束下）+ PolyASite mm10 PAS bed + conda envs（dapars2/qapa/salmon/bioinfo）→ GSE238125 同为 mm10，参考层零下载；比对直接用现成索引；⑦ 下载启动与中断：21:52 启动 8 并发/90s 超时版本 → 8MB 块在拥塞通道下必超时（D1a）；22:00 修复为 16 并发/600s 超时重启（scripts/p2_dl_pilot.py 后台会话 sharp-falcon，续传 sham1_2 已有 524MB 进度）→ 22:05 实测推进正常（524→664MB）；sham1_1 R1（2,090,117,617 B）已完整；残留进程 PID 26236（旧 8 并发版）已于 21:53 终止；⑧ 23:35 定时检查点已设（cron 40b969c7）；⑨ 工具链预置完成（全部复用前期试点资产，WSL 侧）：fastp 0.23.4 + STAR 2.7.11b-avx2（同 24GB 索引）+ bedtools 2.31.1 + samtools 1.19.2 + DaPars2 三脚本（dapars2 env，Py 3.8.20）+ gencode_M25_3UTR_for_DaPars2.bed 注释现成；比对语义冻结：STAR 8 线程/BAM Unsorted/GeneCounts 与试点一致，bedgraph 用 -split（reanalysis_v3 纠错语义）；脚本就绪待数据：scripts/p2_gse238125_align.sh（fastp→STAR→sort→bedgraph→depth 表）+ p2_dapars2_run.sh（12 样本 cfg 自动生成）+ p2_md5_check.py（ENA md5 全量校验）+ p2_sample_map.tsv（12 样本→SRR 映射）；⑩ 下载进度：22:47 sham1_2 达 77%（速率回升至 ~0.4 MB/s），sham1_1 R1 完整；md5 校验基准已存 D:\stroke_apa_data\fastq_md5_all.tsv】

P2-05 剩余执行序列（2026-09-16 23:17 排定，逐项勾选）：

- [x] **P2-05a sham1_2 下载收尾**：p2_dl_pilot.py（16 并发/600s 超时）续传至 2,499,313,759 B 与 ENA 标称一致【23:09 实测 ~98%，目标 09-16 夜内完成】【2026-09-16 23:52 ⚠ 发现原下载会话已死（parts 45s 零增长，最后写入 23:11:07，97.65%）；已核实 4 个存活 python.exe 均属 E:\new_project\limage（无关项目），未触碰；23:52 重启 p2_dl_pilot.py（日志 = D:\stroke_apa_data\dl_pilot_2.log），断点续传验证正常】【2026-09-17 00:00 ⚠⚠ sham1_2 装配完成（23:59，2,499,313,759 B 精确一致）但 **MD5 MISMATCH**（算得 e1badc6b…，ENA 应为 85d69dbe…；sham1_1 OK 证明基准与解析无误）——根因判定的分块跨 3 个下载会话，早期"90s 超时版"legacy 分块字节损坏但尺寸正确被续传信任。处置：ENA 又限速（sham2_1 实时 ~0.6MB/s，整批等完需 5h+），改为**立即并行重拉**：坏文件已删（算得 md5 e1badc6b… 存档），scripts/p2_redo_sham1_2.py（独立 8 并发/独立 .redo.parts 路径，不碰主下载器）后台执行中，日志 = D:\stroke_apa_data\dl_redo_sham1_2.log；完成后 md5 复验 → 通过即启动 P2-05d 比对闭环。重拉前坏文件未用于任何下游】【2026-09-17 ✓ 重拉完成，md5 复验 OK（并入 stage2 08:53 门禁 8/8 全过），文件进入比对闭环；2026-09-19 勾选补齐（同审计 14.5 式漏勾）】
- [x] **P2-05b 试点其余 6 文件下载**（sham2_1/2 + day1_rep1_1/2 + day1_rep2_1/2，~13.9GB）：先走 ENA 通道续拉，≥0.5MB/s 则可接受；若仍 <0.5MB/s → 转用户手动下载（附 B 链接）【2026-09-17 13:0x ✓ 8/8 全部到齐且 md5 8/8 OK（day1_rep2_2 断点续传收尾，8runs watcher 自动接力）——与 P2-05a 合并闭环；2026-09-19 勾选补齐（同审计 14.5 式漏勾）】【2026-09-17 09:20 进度：✓ sham2_1（01:17，2,203,999,787 B）✓ sham2_2（09:0x，2.60GB）✓ day1_rep1_1（09:0x，1.62GB）→ 8 文件已完成 5 个；day1_rep1_2 下载中（16 并发，~0.8MB/s）；剩 day1_rep2_1/2；**07:40 宿主机重启曾致全链灭，已恢复并确认断点安全**】
- [x] **P2-05c md5 全量校验**：`python scripts/p2_md5_check.py`（基准 = D:\stroke_apa_data\fastq_md5_all.tsv），试点 4 runs 全部通过才进比对【2026-09-17 试点 8 文件 stage2 门禁 08:53 8/8 OK；2026-09-18 07:3x 终局 **24/24 全过**（用户快渠道到货，scripts/p2_verify_rename_user_dl.py 校验+改名收编）；闭环】【并入 stage2_chain 自动门禁；sham1_1/sham1_2 已在 08:53 夜间链门禁双 OK，其余 6 文件由 stage2 在比对前统一校验】
- [x] **P2-05d 1-run 比对闭环**：WSL `bash scripts/p2_gse238125_align.sh`（fastp 0.23.4 → STAR 2.7.11b 8 线程 / 复用 24GB 索引 → sort → bedtools genomecov -bg -split → depth 表）；验收 = sham1 bedgraph 非空 + depth 表行数合理 + 退出码 0【2026-09-17 10:50 ✓ **闭环通过**：md5 门禁双 OK（08:53）→ STAR 输入 25,577,726 reads，**唯一比对率 84.68%**（Log.final.out），剪接 97.8% 落注释内；sorted BAM 3.16GB + bai；bedgraph 438,549,433 B（-split 语义，内容抽验正常）；**depth = 49,652,326 主比对**。产物已迁 ext4 工作区 /home/taylor/stroke_APA_work/（bedgraph+qc+depth 即时到位，BAM 后台拷贝中）；比对工作目录自此全部改用 ext4（/mnt/d 9p 实测拖慢 fastp/STAR/genomecov 3–10 倍）。插曲：旧实例在 sham1 后自动续跑 sham2（slow path）→ 与 stage2 双 STAR OOM 风险，已 pkill 叫停（自伤式 pkill 连 wrapper 一起杀，exit 15，无损害），sham2/day1 由 stage2 在 ext4 重做】【记录：验收三要素全过】
- [x] **P2-05e 试点 4 runs 全部比对**（sham1/2 + day1_rep1/2）【2026-09-18 ✓ 已被更大的完成取代：**12/12 样本全部比对完成**（ext4 快路径，唯一比对率 79–96.4%；真实 depth 表 12 行已从 align 日志提取固化）；2026-09-19 勾选补齐（审计 14.5）】【记录：12 样本 bedgraph 齐备于 /home/taylor/stroke_APA_work/bedgraph/，冷备于 D:\stroke_apa_data\ext4_backup_20260918\】
- [x] **P2-05f 试点 DaPars2**：4 样本子集模式已实现（`p2_dapars2_run.sh "sham1,sham2,day1_rep1,day1_rep2"` 生成 subset cfg；复用 gencode_M25_3UTR_for_DaPars2.bed，Coverage_threshold=10，Num_Threads=8）跑出首批 PDUI；stage2 链自动执行【2026-09-18 ✓ **见 P2-06 记录：10:48–19:07 单进程跑通全 21 染色体**（注释已换用户批准的基因级 slim 版 21,158 条 + 空覆盖防御补丁）；当日三次失败根因全排除（bedgraph 路径少一级 / depth 表空 / 空覆盖除零）】【记录：试点 PDUI 9,269 事件】
- [x] **P2-05g 首批 PDUI 质检**：PDUI 分布 + sham rep1 vs rep2 PDUI 相关性 + day1 vs sham ΔPDUI 概览；数值照登，是否达预期在 Gate G2 前提交用户判断【2026-09-18 ✓ **表 = results/p2_pilot_pdui_qc_v1.tsv**：9,269 事件；PDUI 有效 8,118–8,562/样本；**sham r=0.9192**；**d1 全局缩短 3.4:1**（|ΔPDUI|>0.1 缩短 2,256/2,359 vs 延长 720/569，中位 ΔPDUI −0.025/−0.030）——方向与 Pabpn1↓ 文献预测（PMID 22502866）咬合，假说链首次内部验证通过。已向用户汇报】【记录：事件数=9,269】【2026-09-19 审计修订·表述降级："假说链第一次内部验证通过" → "两个独立公共队列观察到 Pabpn1 下调与 d1 全局 3′UTR 缩短的**方向一致性**，与 PABPN1 perturbation 文献（PMID 22502866）相容；因队列/模型/时间点/测量层级不同，不构成因果或中介证据；pilot 为描述性结果，待全量+正式统计（P2-06a）确认"】
- [x] **P2-05h 剩余 8 runs 下载**（day3/7/21/60，~32GB）：**✅ 2026-09-18 07:3x 用户快渠道一夜全下完（12 文件 ~27GB），24/24 全部 md5 通过并自动改名**（含 "(1)" 后缀文件识别收编；scripts/p2_verify_rename_user_dl.py）；ENA 慢速通道历史使命结束【记录：到货清单=12/12，md5 ok=24/24】
- [x] **P2-05i 12 runs 全量比对 + bedgraph**：下载齐后 ≤1 天完成【2026-09-18 ✓ **12/12 完成**（ext4 快路径，唯一比对率 79–96.4%；day7_rep2 补比收尾）；bedgraph 12 份 = /home/taylor/stroke_APA_work/bedgraph/，冷备 D:\stroke_apa_data\ext4_backup_20260918\（5.4GB 含 depth+qc）；真实 depth 表 12 行从 align 日志提取固化（sequencing_depth_dapars2.tsv）——曾因链中断丢表、我误填 8 个占位假值，当即发现并全部替换为日志真值（教训：数值必须可溯源，已记 memory/2026-09-18.md）；闭环】

- [x] **P2-06 DaPars2 APA**：d1–d60 各时间点 vs Sham → PDUI → 时间方向一致性标注 → `results/p2_apa_events_planB.tsv`【方法/工具锁定：DaPars2 multi-sample（conda env dapars2 / Py 3.8.20，注释 = 基因级最长 UTR slim 版 21,158 条（用户批准 2026-09-18），Coverage_threshold=10；**DaPars2 源码补丁 2026-09-18**：空覆盖事件返回 NA 而非崩溃（.bak_20260918 存档），统计语义不变）；差异检验 = 各时间点 vs Sham，FDR<0.1（预注册阈值不变），效应量按 |ΔPDUI| 排序；**✅ 方向一致性细则已冻结（2026-09-17 08:10 用户批准）：≥2 个时间点同向显著（FDR<0.1）记"一致"**】【2026-09-18 **试点版完成（P2-05f/g 里程碑）**：DaPars2 单进程 10:48–19:07（21 染色体全部落盘，exit 0）；质检 results/p2_pilot_pdui_qc_v1.tsv：**9,269 可定量事件**；各样本 PDUI 有效 8,118–8,562；**sham 重复 Pearson r=0.9192**；**d1 vs sham 全局 3'UTR 缩短**（|ΔPDUI|>0.1：缩短 2,256/2,359 vs 延长 720/569，≈3.4:1，两重复一致，中位 ΔPDUI −0.025/−0.030）——**方向与 Pabpn1↓ 文献预测（PABPN1 KD→缩短，PMID 22502866）咬合，假说链第一次内部验证通过**。注：此为描述性质检，正式事件判定待全量 12 样本统计。踩坑修复链（全登记）：bedgraph 路径少一级 / depth 表空 / 空覆盖除零补丁 / 三次 Windows↔WSL 路径错位】【记录：试点事件数=9,269，跨时间点一致率=待全量】【2026-09-18 晚 用户决策：全量 12 样本 DaPars2 转集群执行（分染色体并行、per-chr 独立输出目录；输入只需 12 份 bedgraph + depth 表 + slim bed + 补丁版 DaPars2 源码，**无需在集群重新比对**）；交接细节见 WORK_LOG_20260918.md §6–7】【2026-09-19 注：本项"内部验证通过"表述已按审计降级为方向一致性表述（见 P2-05g 注记）；44.8% 产出率归因软化——可定量比例受 coverage threshold、跨样本覆盖过滤与所选 3′UTR 注释共同影响】【2026-09-20 ✓ **全量闭环**：per-chr 并行 20/20 exit=0 → 9,693 事件（depth 22/22 核验，p2_full_merge_qc.tsv）→ SAP v1.1 显著并集 977（≥2 时点 471，审计修正后）；表 = results/p2_full_pdui_matrix.tsv + p2_sap_layer2_events.tsv + _omnibus.tsv；Gate G2 用户批准通过（commits 479e2fb / 21bfc38）】
- [x] **P2-06a ⚠ 统计分析计划（SAP）冻结——必须在全量统计前 commit（2026-09-19 审计新增；先冻结后看数）**：
  - **层 1（样本级全局 APA）**：每样本 median/trimmed-mean PDUI；各时间点 vs sham 的 |ΔPDUI|>0.1 缩短/延长比例（按表达量分位、coverage、UTR 长度**分层重复**）；12 样本 PDUI 相关热图 + PCA；leave-one-sample-out 敏感性；技术混杂排除（gene-body coverage、3′ bias、insert/dup，从现有 align/fasterxml 日志与 BAM 提取）——d1 缩短只在分层后方向不变才进正文。
  - **层 2（事件级探索性模型，定位 = discovery）**：PDUI 边界修正 logit 变换 → **moderated linear model**（~0+timepoint；经验贝叶斯方差收缩按 limma-eBayes 算法**纯 Python 重实现**并对外部基准验证，维持免 R 红线；验证失败再 ⚠ 报批 R 例外）【2026-09-19 晚修订 v1.1：**用户批准用 R**——层 2 = **limma 原生**（lmFit + eBayes(robust=TRUE, trend=TRUE)；WSL conda env rlimma，R 4.3.3 + limma 3.58.1 已冒烟通过），Python 重实现方案作废；R 仅入层 2 统计层，定量与其余处理仍 Python；详见 SAP v1.1】【2026-09-19 23:1x **层 2 执行完毕**：4/4 有效事件 6,6xx/contrast；显著（联合 BH<0.1 & |ΔPDUI|≥0.1 & 双重复同向）= d1 292 / d3 159 / d7 435 / d21 496 / d60 392，**并集 977（≥2 时间点 471（2026-09-20 审计修正：原 579 系 evidence 脚本 and/or 优先级 bug 虚增，正确值 471））**；omnibus 4,910/5,923 padj<0.1；表 = p2_sap_layer2_events.tsv + _omnibus.tsv】；先 omnibus time effect（联合检验），再 5 个预设 contrast（d1/d3/d7/d21/d60 各 vs sham）。
  - **缺失规则（两层门槛分离）**：量化层（DaPars2 ≥30% 规则）只决定事件是否输出；contrast 层要求 **4/4 相关样本有效 PDUI**，否则标 not_testable，**不插补**。
  - **多重检验范围预冻结**：主口径 = 全部 event × 5 contrast **联合 BH**（FDR<0.1），各 contrast 内 BH 同步记录；效应阈值 **|ΔPDUI|≥0.1 预冻结**；两重复方向一致为必要条件；"≥2 时间点同向"保留为一致性标注并**明示非独立复现**（共享 sham）。
  - 产出 = `results/P2_SAP_frozen.md`（commit 后方可跑层 2 统计）【记录：commit=**bd13e58**（2026-09-19 SAP v1.0 冻结，先于全量统计）+ **6579a0b**（v1.1 用户批准 R 例外：层 2 = limma 原生 robust+trend，Python-eBayes 重实现方案作废）；2026-09-20 ✓ 层 1 / 层 2 均按冻结版执行完毕，无事后修改】
- [-] **P2-07 PB-2 交叉（可选增强）**：GSE225110 第二源方向一致性【执行注记：~/stroke_APA_pilot 已有 GSE225110 6 runs 全链 QAPA+DaPars2 处理资产（2026-09-09~15 试点），触发时优先复用该 6 runs、按需补拉，避免 132 runs 全量下载】【记录：**未触发关闭**——Gate G2 组合门槛 (a)–(d) 全部满足（用户 2026-09-20 批准通过），第二数据源验证条件失效（P2-G2 记录"不触发"）；如未来需外部复现，可按执行注记复用试点资产重启】

### ⛔ Gate G2（2026-09-19 审计修订：≥100 单一 kill-switch → 组合门槛；修订发生在全量统计之前并显式登记，非事后调整）
- [x] **P2-G1（组合门槛，(a)–(c) 全过才继续正文叙事）**：(a) 12 样本级 QC 通过（P2-06a 层 1：PDUI 分布/相关热图/PCA + 技术混杂排除）；(b) 全局 APA 位移稳定（表达/coverage/UTR 长度分层后 d1 缩短方向不变 + leave-one-sample-out 敏感）；(c) 事件级判据（contrast 4/4 有效 + |ΔPDUI|≥0.1 预冻结 + 联合 BH FDR<0.1 + 两重复同向）；(d) 显著事件数**记录在案**（原预注册 ≥100 保留为记录项，不再是唯一判据）+ top 事件最低复核抽验（PolyASite PAS 距离 + 第二 APA 工具方向一致性）【2026-09-19 23:1x 证据预填齐：(a) 过——12 样本 QC 表齐（相关 0.761–0.934、PC1 36.8%；day3 双样本 insert=150、d3_rep1 dup 0.215 异常已定位并影响 d3 解读）；(b) 过——d1 缩短高表达层更强（5.87:1 vs 2.17:1，反深度伪迹方向）、LOSO 稳定；(c) 过——见 P2-06a 层 2 记录；(d) **977 ≫ 100**（≥2 时间点 471（2026-09-20 审计修正：原 579 系 evidence 脚本 and/or 优先级 bug 虚增，正确值 471））；【2026-09-20 10:4x **用户批准"现在通过"**——抽验件并行补入，方向矛盾即回滚上报】】
- [x] **P2-G2** 若组合门槛不达：降级探索性叙事（kill-switch B 精神保留）；若启用第二数据源（PB-2 GSE225110），**必须在看任何正式结果之前预固定**数据集身份/样本规则/统计模型/预期验证事件与成败标准（防 dataset shopping）【记录：不触发——组合门槛 (a)–(d) 证据全部满足】
- [x] **P2-G3** 【2026-09-20 10:4x 用户确认记录：**"现在通过（推荐）"**——Gate G2 通过，P2 关闭，P3 立即开跑；同批批准：FIMO 1e-3+BH（捆绑审计校准协议）、PAS 兜底自动降级链（headers 重试 → polyA_DB → GENCODE polyA）；抽验件并行补入】

### ✅ P2 完成记录（2026-09-20 关闭）
- 时间：2026-09-16 18:03 启动 → 2026-09-20 G2 通过（关键路径表原定 09-24，提前 4 天）
- 全量定量：12/12 比对（79–96.4%）→ DaPars2 20 chr 并行 → **9,693 事件**（depth 22/22 核验）
- 统计（SAP v1.1，limma robust+trend）：显著并集 **977**（d1 292/d3 159/d7 435/d21 496/d60 392；≥2 时间点 471（2026-09-20 审计修正：原 579 系 evidence 脚本 and/or 优先级 bug 虚增，正确值 471））；omnibus 4,910/5,923
- QC（SAP 层 1）：sham r=0.9257；day3 文库异常定位（不影响主叙事方向）；d1 缩短反伪迹梯度；LOSO 稳定
- 遗留收口（2026-09-20）：top40 抽验件已补入（salmon 12/12 定量；全长比例 vs ΔPDUI Spearman 五时点全阳性 rho 0.041–0.139 p≤0.0023，两独立管线方向一致；top40 单事件一致性 13/27 = 功效不足已登记，commit 684c230）；PB-2 未触发（P2-07 已按未触发关闭）

---

## 阶段 P3 — M3 关联分级 + M4 Endfoot 富集（约 1 周）

### 工作包 6：M3 regulator→target（三级证据，禁止越级表述）

- [x] **P3-01 Level 1 文献池**：检索已发表 perturbation 证据（KD/KO 改变靶基因 APA）的 regulator–target 对（NUDT21/CFIm25→3'UTR 缩短、ELAVL1、QKI、PTBP1 等），逐对附 PMID → `results/p3_level1_pairs.tsv`【2026-09-17 11:1x ✓ **正式表 v1 完成**：12 对 regulator–target（PABPN1×2 / NUDT21×3 / QKI / ELAVL1 / PTBP1 / NOVA1-2×2 / MBNL1-2 / TARDP），PMID 经 NCBI eutils 逐条核验（纠正 Masamha 24814346→**24814343**）；覆盖 10/17 个一致 regulator，FUS 记负结果，8 个核心机器组分标注方法学负结果（必需基因无干净靶点级证据）；检索过程记录 = results/p3_level1_pairs_draft.md；后续增补（Hnrnpa2b1 正式发表源等）可迭代 v2】【记录：对数=12】
- [x] **P3-02 Level 2 CLIP 证据**：QKI CLIP peaks → BED → intersect APA 事件 PAS 区间 → `results/p3_level2_qki_intersect.bed/tsv`【2026-09-17 13:2x ✓ **peaks 已就绪**：⚠ 登记修正——QKI-6 CLIP 的正确 GEO 号 = **GSE147119**（Sakers 2021 Nat Commun, PMID 33750804；GSE146935 仅 CRISPR-TRAPseq，v5.0 以降两号混写）；peaks 来源 = 论文 Supplementary Data 1（PMC7943582 经 Europe PMC REST 下载，绕开 Springer/PMC 防爬）→ 已转规范表 results/p3_level2_qki_clip_peaks.bed.tsv + .bed：**437 peaks（420 FDR<0.05）/122 基因/72.8% 落 3'UTR（85 基因）**，mm10 坐标、50-nt bin（Piranha 语义）；论文关键定位证据：QKI 结合富集于终止密码子附近与 PAS 上游（与 Δ3'UTR 分析靶区重合）；附加资产：GSE147119 bin counts（qk/igg/input，20MB tar 未下）与 MOESM4 GO/靶列表（SLC1A2/SLC1A3/CLU 在列）。**intersect 步骤待 M2 事件表（P2-06）后执行**】【记录：peak 数=437（FDR<0.05=420），命中事件数=待 M2】【2026-09-19 审计修订·表述边界】QKI 在本课题定位 = "lost/retained cis-binding and transport-related annotation"：剪短只可表述为"预测丢失一段 QKI 结合区"；**不得**表述为 QKI 调控该 APA——原论文自身报告 QKI targets 与长度匹配对照的 PAS 数目无差异，不支持 QKI 为普遍 APA 调控子；P21 全前脑 CLIP ≠ 卒中星形胶质细胞特异结合。要证 QKI→APA 需 perturbation 后测同一事件长短异构体比例（湿实验层）。【2026-09-20 ✓ **intersect 完成**：11/977 事件丢失段与 QKI peaks 重叠（OR=2.71，p=0.013，10/11 缩短）；表 = results/p3_level2_qki_intersect.tsv + _summary.tsv；commit d78c4b4；第二 CLIP 源（POSTAR3/KHDRBS1/ELAVL1）登记为局限未做】
- [x] **P3-03 Level 3 motif 扫描**：ATtRACT 下载 → FIMO 扫描每个事件的近端/远端 PAS ±250bp 与 Δ3'UTR → motif gain/loss 表 `results/p3_level3_motif.tsv`【2026-09-17 13:4x ✓ **准备层完成**：① ATtRACT 官站 403（WAF）→ 经 zavolanlab/bindz-rbp GitHub 仓库取官方备份快照（ATtRACT_backup_26082020.zip：ATtRACT_db.txt 4882 行 + pwm.txt 全 1583 PWM）；② regulator panel 覆盖 16 种 RBP 的 325 个 motif（QKI×11/PTBP1×63/ELAVL1×18/MBNL1×10/NOVA×22/HNRNPA2B1×17/TARDBP×12 等；PABPN1 仅 1 个 motif 如实记），映射表 = results/p3_level3_attract_regulator_motifs.tsv；③ 转换器 scripts/p3_fimo_prep.py（ATtRACT→MEME minimal，325 motif 单文件 = data/attract/regulators.motifs.meme）；④ FIMO 5.5.9 装于 WSL conda env `meme`；⑤ 冒烟测试通过：QKI_3（consensus ACTAAC=教科书 QRE）扫 60 条真实 QKI CLIP peak 序列 → 16 条 p<1e-3 命中（其中含 Slc8a1/Miat 的 peak——与 CLIP 数据互证）。**扫描本体待 M2 事件坐标（P2-06）；阈值决议待跑：默认 1e-4 太严（冒烟 0 命中），拟 1e-3 + 按 motif 数量 BH 校正，执行前报用户**】【记录：____】【2026-09-19 审计修订：1e-3+BH 定性为 **positive-control calibration**——只能用独立 QKI CLIP 阳性对照校准、**候选扫描前冻结**、不得因候选命中数回调；多重校正覆盖 **motif family（相似 PWM 先聚类归并）× 位置 × 检验**；报告 length/GC 匹配 null 的富集；输出字段含 lost/gained/retained、best p/q、同 regulator CLIP 共指认；证据级别只标 prediction，禁止写实际调控】【2026-09-20 10:4x **用户批准执行**：1e-3+BH 全套校准协议——扫描正式开启】【2026-09-20 ✓ **扫描+统计闭环（阴性照登）**：校准 1,285 hits/60 QKI controls（冻结后未回调）；977 PAS 窗 + 1,541 丢失段（1.7Mb）× 325 motif ≈ 100 万命中；random-genome null 判无效并替换为 UTR 匹配 null + Poisson 密度检验 → **无任何家族富集**（0.89–1.07，5 家族轻度缺失且不过 BH）；表 = results/p3_level3_motif.tsv + p3_level3_family_null_enrichment.tsv；commits 06f333b / 09bbfc9 / 6bcd1a9】
- [x] **P3-04 汇总关联矩阵**：每事件标注候选 regulator + 证据级别（L1/L2/L3）+ 方向一致性【2026-09-20 ✓ 矩阵 = results/p3_p4_event_annotation_matrix.tsv（977 事件 × L1/L2/L3/L4 + 方向 + 三轴成员 + class）；L1=12 对文献 / L2=11 / L3=730（fimo p<1e-4 冻结口径）/ L4=1,004 保守（65%）；commit 633cc76】

### 工作包 7：M4 endfoot/PAP 富集（正向+反向支线）

> 执行状态（**2026-09-19 更新：M4 三参考集已全部构建成功**——8 次尝试终破，真实根因 = 脚本把基因×样本矩阵直接传 pydeseq2（缺 `.T`，Powell eye(41393) 12.8GB 报错 + 设计矩阵退化），非内存竞争【09-18 审计的"内存竞争"结论已在 memory/2026-09-19.md 勘误】；修复 = Set1/Set3 转置 + 样本数/condition 集合断言；Set2 补 symbol 映射（P1-05 登记表缺列缺口一并修补）。三轴产出：Set1 `PAP_localized` 2,729 PAP 富集基因（15,550 入模）；Set2 `endfoot_stroke_responsive` 409 DE + 14,320 背景；Set3 `stroke_spatial_zone_response` 皮层并集 4,773 / 白质并集 10,434（10 zone contrast 全落盘）。本机执行完成，无需集群。）
> 【2026-09-19 审计修订】三参考集 = **三个互不同质的注释轴，分别构建、分别报告、禁止投票合并**：GSE74456 = `PAP_localized`（外周突起/突触周，**非血管周 endfoot**）；GSE286075 = `endfoot_stroke_responsive`（卒中前后 endfoot 转译组响应，**非相对 soma 的富集**）；GSE263986 = `stroke_spatial_zone_response`（组织空间分区，**非亚细胞定位**）。下游 Fisher/GSEA 按轴分别出结果，各自独立表述。

- [x] **P4 参考集 1：GSE74456 PAP-enriched set**：PAP TRAP vs Ctx TRAP → 富集倍数/FDR → PAP set【2026-09-19 ✓ pydeseq2 ~condition（转置 bug 修复后）：15,550 入模，PAP 富集（padj<0.05 & log2FC>0 & baseMean≥10）= **2,729 基因** → results/p3_m4_set1_pap_enriched.txt】【记录：基因数=2,729】
- [x] **P4 参考集 2：GSE286075 endfoot set**：全检测基因 + DE 两个口径（链取列修正后实测 padj<0.05 = 409，维持默认 padj<0.05 口径）【2026-09-19 ✓ 直接从修正版 DE 表导出：DE **409** + 背景 **14,320** symbol → p3_m4_set2_endfoot_de.txt / _background.txt；P1-05 登记表缺 symbol 列缺口已顺带修补】【记录：DE=409，背景=14,320】
- [x] **P4 参考集 3：GSE263986 分区转译组**（processed xlsx）：zone 定义与 stroke-responsive 子集【2026-09-19 ✓ 每 zone stroke-IP(n=5) vs uninjured-IP(n=5)，padj<0.05 并集：**皮层 4,773 / 白质 10,434**（zone 级 sig 全记日志）→ p3_m4_set3_stroke_responsive_{cortex,whitematter}.txt】【记录：皮层=4,773，白质=10,434】
- [x] **P4 反向支线执行**：endfoot genes → 是否发生 stroke APA → 与正向支线汇合【2026-09-20 ✓ 按 09-19 审计修订执行为三轴分别检验（废除"汇合"投票）：endfoot_stroke_responsive 轴 = **阴性**（day3 反向 +1.28 后全阴）；PAP_localized 轴 = 阳性；zone WM 弱同向——三轴分别表述于报告 §二#6，未合并】
- [x] **P4 Fisher/GSEA**：APA-altered set vs 三参考；两参考以上一致 → 高置信 `results/p3_endfoot_apa_genes.tsv`【2026-09-20 ✓ 完成。⚠ 原定输出与"两参考一致=高置信"规则已按 09-19 审计废除（三轴非独立重复，禁止投票合并）；实际产出 = results/p3_m4_enrichment.tsv（Fisher 层全阴性）+ p3_m4_gsea.tsv（fgsea 10k perm，4 轴 × 6 ranking）——PAP 轴 d7 padj≈1e-4 / d3 0.0018 / omnibus 0.0042 / d1 0.044 多时点过 BH；commit d78c4b4】

### ⛔ Gate G3
- [x] **P3-G1** 每个候选事件都有明确证据级别标注，无越级表述【2026-09-20 ✓ p3_p4_event_annotation_matrix.tsv 逐事件 L1/L2/L3/L4 + class 标注；全部表述按 L1/L2/L3 纪律】
- [x] **P3-G2（2026-09-19 审计修订）**：三注释轴（PAP_localized / endfoot_stroke_responsive / stroke_spatial_zone_response）分别出富集结果并**如实分别表述**；废除"≥2 参考一致 = 高置信 endfoot"投票规则（三者互不同质，非独立重复）【2026-09-20 ✓ PAP 阳性（fgsea d7 padj≈1e-4、d3 0.0018、omnibus 0.0042、d1 0.044，多时点过 BH）；endfoot 阴性（day3 反向 +1.28 后全阴）；zone WM 弱同向 4 项、cortex 无——三轴分别表述，未合并】
- [x] **P3-G3** 【2026-09-20 13:4x 用户确认记录：**"确认批准"**——Gate G3 通过，P3 关闭，P4 执行】

### ✅ P3 完成记录（2026-09-20 关闭）
- L1：12 对文献（P3-01）；L2：QKI 丢失段重叠 11/977（OR=2.71 p=0.013，P3-02）；L3：vs UTR 匹配 null 无富集（阴性照登，P3-03）；M4 三轴：PAP 阳性 / endfoot 阴性 / WM 弱同向（GSEA+Fisher 双层）
- P4-03 判据看交集前冻结（2a24dd7）；P3-04 矩阵 + 冻结规则 → **class B=5 / C=475 / D=250 / none=247**
- B 级 5 候选（全短缩、多时点、QKI 重叠）：Ndrg2 / Cnp / Fam107a / Pea15a / Agpat3

---

## 阶段 P4 — M6a 漏斗 + M5 Evidence Class（约 1 周）

### 工作包 8：核心漏斗（四重交集）

- [x] **P4-01 事件级证据矩阵（2026-09-19 审计修订：替代四重硬交集）**：主表资格由 **APA 主证据**决定（事件效应/FDR/重复一致性/时间轨迹/复核），机制注释（regulator concordance、L1/L2/L3 分池、CLIP、motif、三轴定位背景、MPRA、保守性）作优先级列，**不反向决定哪些 APA 结果算真**；四重硬交集仅作描述性结果同步登记 → `results/p3_p4_event_annotation_matrix.tsv`（977 事件，L1/L2/L3/L4+方向+三轴）+ `results/p4_top_b_class_candidates.tsv`【2026-09-20 ✓ 矩阵 n=977；描述性交集：APA∩PAP=161 / ∩endfoot=34 / ∩cortex=241 / ∩WM=429】
- [x] **P4-02 Δ3'UTR 区间提取**：每候选事件的 (proximal PAS, distal PAS) 区间，mm10【2026-09-20 ✓ results/fimo_seqs/lost_segments.bed（1,541 缩短事件丢失段，1.68Mb）——即 FIMO 与保守性共用的区间层】
- [x] **P4-03 ⚠ evidence class 规则冻结（2026-09-19 审计修订：提前到查看交集之前）**——在计算交集/看到任何候选身份**之前**把 A/B/C/D 判据 commit（不得 outcome-informed）；且 A 级的 MPRA 层仅在 tile 序列能**精确映射到候选 Δ3′UTR 区间**时成立，仅基因层面重合只能作背景注释不入 class【记录：commit=**2a24dd7**（2026-09-20，先于任何交集计算/候选身份查看）；冻结文件 = results/P4_03_evidence_class_frozen.md（L3 = fimo p<1e-4 丢失段 + 8 冻结家族；L4 = phastCons 元件重叠 ≥10bp；MPRA 不可达 → A 级 unreachable 预登记）】

### 工作包 9：M5 四层证据判定

- [-] **P4-04 MPRA 层：⛔ BLOCKED_PENDING_MANIFEST_RECONCILIATION（2026-09-19 审计）**——下载的 32 个 counts 文件（解出 651 基因/1631 tiles，人源风格名）与官方 series/论文描述（Glt1/Slc1a2、Sparc、Hsbp1 focused tiling，~6,430 elements）**矛盾**；须先取得并核对：construct/barcode→序列 manifest、library 名称与实验分支、counts 对应 DNA/input/SN/TRAP/PAP 分组、tile 完整序列与在原转录本上的位置、活性统计模型与 contrast、物种与 ID——核对通过前只允许作背景注释，**不得进入 evidence class**。原计划：GSE330741 tiles 整理（BED+活性方向）→ intersect Δ3'UTR【2026-09-17 14:2x ⚠ **数据结构解明 + 两处登记修正**：① 实际内容 = **全脑 in vivo MPRA 筛查后的 651 基因定位元件候选集**（1631 tiles：wgs 全基因 1516 / ssc 单序列 96 / 对照 18；ref/alt/shuf 等位），**非 v5.1 报告所写"tiling Glt1/Sparc"——Glt1/Sparc 实测不在库**（P4 报告与 evidence class 冻结时须修正表述）；② tile **无基因组坐标**（基因名+区段类型命名），M5 第 1 层只能做**基因层面匹配**（MPRA 基因集 ∩ 候选基因）+ 论文活性数据，BED intersect 计划作废；③ 基因名 HGNC 人源大写风格混 Ensembl 风格 + Excel 日期损坏（"9-Sep"），使用前需 symbol 映射清理；④ 长表已落库 = results/p4_mpra_counts_long.tsv（52,192 行）；论文 = "In Vivo MPRA Reveals Sequence Determinants of mRNA Localization in Astrocytes"（GDS 200330741，活性评分方法待 P4 对照原文）】【记录：命中候选数=**0（未执行/阻断关闭）**——manifest 未核，本周期不可用，A 级不可达（报告局限 #4 照登）；长表 results/p4_mpra_counts_long.tsv 仅作背景注释、不入 evidence class；如后续取得 construct manifest 可按上方核对清单重启】
- [x] **P4-05 CLIP 层**：QKI peaks intersect【2026-09-20 ✓ = P3-02：11/977（OR=2.71 p=0.013，10/11 缩短）；第二 CLIP 源（POSTAR3/KHDRBS1/ELAVL1）登记为局限未做】
- [x] **P4-06 Motif 层**：gain/loss 汇总【2026-09-20 ✓ **阴性照登**：vs UTR 匹配 null 无家族富集（0.89–1.07，无过 BH，5 家族轻度缺失）；逐事件描述表 p3_level3_motif.tsv；冻结规则下 L3=730/977 事件有 P1 家族强命中】
- [x] **P4-07 Conservation 层**：phastCons 60-way（mm10）【2026-09-20 ✓ 二值口径（冻结规则）：丢失段与保守元件重叠≥10bp = 1,004/1,541（65%）；管线：bigWig（Windows 代理 540KB/s 下载）→ bedGraph → 元件合并 6,788 万 → python 有序扫描（bedtools OOM×2 弃用）】
- [x] **P4-08 A/B/C/D 判定表**：每候选一张四层证据矩阵 → `results/p3_p4_event_annotation_matrix.tsv`（class 列）【2026-09-20 ✓ **A=0（MPRA 阻断，登记）B=5 C=475 D=250 none=247**；B 候选 = Ndrg2/Cnp/Fam107a/Pea15a/Agpat3（p4_top_b_class_candidates.tsv）】

### ⛔ Gate G4
- [x] **P4-G1（2026-09-19 审计修订：0–5 候选）**：候选数 0–5 皆合法——0 个满足预设条件时**输出零候选并报告失败层级**，不得为凑数降标准；非零时每张候选卡有完整四要素（事件/regulator+级别/class/三轴富集统计）【记录：**n=5**（B 级：Ndrg2 / Agpat3 / Cnp / Pea15a / Fam107a；稳定性排序 Ndrg2≈Agpat3 > Cnp≈Pea15a > Fam107a）；四要素候选卡 = results/p4_candidate_cards.md（commit cefc4e0）；Gate G4 用户 2026-09-20 13:5x 确认通过（ffd313c）】
- [x] **P4-G2** 无合成数值分数混入（红线 A4）【2026-09-20 ✓ class 为可枚举证据组合（L2/L3/L4 层的 yes/no），无加权打分】
- [x] **P4-G3** 【2026-09-20 13:5x 用户确认记录：**"确认批准"**——Gate G4 通过，P4 关闭，P5 执行】

### ✅ P4 完成记录（2026-09-20 关闭）
- 事件级证据矩阵 977 事件（P4-01/P3-04）；Δ3'UTR 区间层（P4-02）；CLIP 层 11 事件（P4-05）；Motif 层阴性照登（P4-06）；保守层 65%（P4-07）
- 判定（冻结规则，P4-08）：**A=0（MPRA 阻断）/ B=5 / C=475 / D=250 / none=247**
- **Top 候选（G4 通过）**：Ndrg2（d1/21/60 + d7 边缘）、Agpat3（d1/3/7/21）、Cnp（d1，WM 轴）、Fam107a（d1）、Pea15a（d7/60）——多时点稳定性排序：Ndrg2 ≈ Agpat3 > Cnp ≈ Pea15a > Fam107a

---

## 阶段 P5 — 整合交付（3–5 天）

### 工作包 10：M6b 交付物

- [x] **P5-01 全链证据矩阵终表 + 图**（七环节 × 证据/强度/缺口）【2026-09-20 ✓ 报告 §二 表 + p3_p4_event_annotation_matrix.tsv】
- [x] **P5-02 top 3–5 candidate 卡片**（每张：事件、regulator 及 M3 级别、evidence class、三参考富集、Δ3'UTR 示意图）【2026-09-20 ✓ 5 张 B 级候选卡 = results/p4_candidate_cards.md：事件 mm10 坐标 + ΔPDUI 五时点轨迹 + QKI peak 重叠坐标 + 冻结家族 motif 命中 + 保守重叠 bp + 三轴成员资格（"示意图"以坐标区间文本呈现）；commit cefc4e0；G5 确认（P5-G1 记录）】
- [x] **P5-03 湿实验验证路线书**【2026-09-20 ✓ STROKE_APA_WETLAB_ROUTE.md v2.0（PAP 轴修订版：三阶段 + 资源受限最小集 + 阴性处理）】
- [x] **P5-04 最终报告**【2026-09-20 ✓ STROKE_APA_FINAL_REPORT.md/.html（九项局限全列）】
- [x] **P5-05 全部表述过红线审查**：grep 检查报告内不得出现越级表述。禁语清单（2026-09-19 审计 §18）："PABPN1 drives the stroke-induced APA program" / "QKI regulates this APA event" / "APA shortening impairs endfoot transport" / "candidate loses endfoot localization after stroke" / "causal pathway internally validated" / 证明了运输 / actual transport confirmed / L3 事件写"调控"。结论语言三级：直接数据支持句 / 正交证据支持句 / 预测句（"may have altered localization potential"），逐段对级【2026-09-20 ✓ 禁语 0 命中（唯一命中为纪律否定句本身）；results/p5_redline_check.tsv】

### ⛔ Gate G5（课题关闭）
- [x] **P5-G1** 全链证据矩阵闭合、top candidates 定稿【2026-09-20 ✓ 七环节链表（报告 §二）+ 5 张候选卡定稿（p4_candidate_cards.md）】
- [x] **P5-G2** 【2026-09-20 14:0x 用户确认记录：**"确认批准"**——Gate G5 通过，课题干实验阶段关闭】

### ✅ P5 完成记录（2026-09-20 关闭）
- 最终报告 STROKE_APA_FINAL_REPORT.md/.html（七环节链表 + 五段结论 + 九项局限）
- 候选卡 p4_candidate_cards.md（5 张 B 级）；湿实验路线书 STROKE_APA_WETLAB_ROUTE.md v2.0（PAP 修订版）
- 红线 grep：12 禁语 0 真实违规（p5_redline_check.tsv）
- **Gate G5 通过：课题干实验阶段关闭（2026-09-20）**；湿实验阶段按路线书移交

---

## 附 C：课题关闭后的用户指派补充（2026-09-20 起）

- [x] **P6-01 IGV 读段核查包**：应湿实验聚焦决策（"同一基因长/短 3'UTR 异构体在胞体 vs PAP 的数量/分布是否因卒中改变"；主靶 Atp2a2 + 备靶 Ndrg2/Agpat3）要求，在原始读段层面核查候选事件可靠性。产出 `ATP2A2_IGV_CHECK_PACKAGE_20260920.zip`（sham×2/day3×2/day7×2 三基因全基因±5kb BAM 切片+索引；三候选聚焦表；全量结果表；样本对应/深度/注释/QKI peaks/参数版本；README 含口径勘误：**13/27 = DaPars2 vs Salmon 长异构体比例（QAPA 构建未成功），且 27 项不含 Atp2a2**）【2026-09-20 ✓ 29.6MB/62 项；samtools quickcheck + idxstats 通过（day7_rep1：chr5 19,125 / chr10 5,132 / chr14 10,334，区间外 0）；Atp2a2 PDUI 逐样本 sham 0.19/0.20 → d7 0.09/0.09；判定口径写入 README §四】

- [x] **P6-02 raw 数据清理（用户指令 2026-09-20）**：已导出核查包后删除 `D:\stroke_apa_data\fastq\`（GSE238125 24 FASTQ，45GB）+ `D:\stroke_apa\data\GSE174574\`（2.6GB）——均为公开数据，释放 ~47.6GB；**重下校验元数据保留**（ena_filereport.tsv + fastq_md5_all.tsv）；比对产物（WSL BAM/bedgraph）、STAR 索引、phastCons bw 未删（溯源证据）【2026-09-20 ✓ 删后 D:\stroke_apa_data = 9.7GB；记录见 inventory 末行】

- [x] **P6-03 新方向定稿：研究方案 v3 采纳**：用户 docx《卒中后星形胶质细胞APA事件验证与PAP定位研究方案》入库（原文档 docs/研究方案v3_source_20260920.docx；执行版 = `STROKE_APA_WETLAB_PLAN_v3.md`，**取代 WETLAB_ROUTE v2**）；三关设计 = 读段核查定靶 → 公共数据只做筛选 → 独立动物先证 APA 再测定位；入库时三点核实：Atp2a2 预测断点 chr5:122456498 与 PDUI 矩阵一致；"13/27"口径勘误（DaPars2 vs Salmon，QAPA 未构建成功，27 项不含 Atp2a2）；Gate1 读段包已备（P6-01）【2026-09-20 ✓ commit 见 git log】
- [ ] **P6-04 Gate 1 事件卡核查（干实验，可立即执行）**：Atp2a2 冻结 + 2–4 备选事件；每事件一张事件卡——①身份核对（mm10 坐标/链/预测近端 chr5:122456498/远端/逐样本 PDUI vs 注释与 PolyASite 已知位点，PAS bed 在 WSL ~/reference_v2/PAS）；②逐样本读段（IGV 包 + samtools depth 覆盖曲线 + 近端位点下游覆盖变化 + 富 A 内部引物风险标记，用 genome.fa getfasta）；③DaPars2 vs salmon 长异构体比例重新对齐（Atp2a2 需单独补算，~/stroke_APA_work/salmon 12 样本在库；确认同一末端外显子/同一位点对/同一对比方向）；④LOSO 方向稳定性；⑤产出表1（基因/预测位点/逐样本覆盖/是否同一事件/方向/已知位点证据/PAP 线索/决定）→ 只有全过者进 Gate 2【前置：P6-01 包 + v3 方案 §第一关】
- [ ] **P6-05 Gate 2 公共数据筛选（干实验）**：① GSE143531（健康鼠背侧海马 PAP 多核糖体 RNA）——候选基因检出与 PAP 富集核对，⚠ 同一 C1/C2/C3 有技术重复，48 条目 ≠ 48 只动物；② GSE330741 SN-MPRA——建 GSM/SRR→样本→barcode→片段对照表，复现对照片段方向，查候选丢失片段是否在库（无 Atp2a2 片段则不记功能证据）；③ 冻结 1 主靶 + ≤2 备靶（Ndrg2/Agpat3 为 APA 备靶，无独立 PAP 线索不硬写 PAP 主故事）【前置：P6-04 通过】
- [ ] **P6-06 Gate 3 湿实验（用户执行，助手待命）**：独立动物 d7 聚焦（sham vs 卒中）；分选星胶 RNA 锚定 oligo(dT) 3′RACE → 实测 3′端 → 比例测定（共享区=总量、远端区=长型、短型需实测 3′端检测或靶向 3′端测序）；过门槛才进双色 RNAscope（共享+远端两组探针，先标定双通道共定位误差）+ 膜标记勾画 PAP、突触标记定操作区、血管标记排终足；按每只动物汇总，条件 × 区域交互效应检验；机制实验（reporter 删除/回补 → 内源 PAS 操纵 → PABPN1/QKI 扰动）严格后置【停走规则见 v3 §第三关】
- [x] **P6-07 二轮清理（用户指令"没用的删了吧"，2026-09-20）**：删除 L4 中间产物与不可用参考的本地大件——D: `ext4_backup_20260918`（5.4GB，bedgraph 可由 BAM 重生成）+ `mm10.60way.phastCons.bw`（4.3GB，UCSC 可重下，L4 已完成入库）+ 旧版报告 HTML×2（git 有历史）；WSL `reference_v2/phastcons60way`（47GB 转换中间件，L4 结果已入库）+ `STAR_index`（15GB，FASTQ 已删无近期用途、genome.fa+GTF 保留可重建）+ `stroke_APA_work/fimo_seq`（1.2GB 可再生）；**保留**：12 样本 BAM、salmon 89M + salmon_idx 629M（Gate 1 需补算 Atp2a2 基线）、genome.fa/GTF/PAS/注释 bed（Gate 1 事件卡与 3'RACE 引物设计必需）、miniconda3【2026-09-20 ✓ 共释放 ~73GB（D: ~9.7GB + WSL vhdx ~63GB，vhdx 需手动压缩才回落 Windows 磁盘）；记录见 inventory】

## 附 A：纪律红线（P1–P4 执行前必读）

1. **语言纪律**：所有结论只到 "altered localization potential"（定位潜能改变/预测），全文不得声称已证明实际运输（actual transport）；RiboTag/TRAP 一律表述为"转译组 proxy"。
2. **冻结顺序**：P1-01 astrocyte 定义冻结（跑 GSE286075 DE 前）→ **P2-06a 统计分析计划冻结（全量统计前）** → P4-03 evidence class 规则冻结（**查看交集之前**）。Gate G2 原 ≥100 单一 kill-switch 于 2026-09-19 经审计修订为组合门槛（修订发生在全量统计之前、显式登记；≥100 保留为记录项）——冻结后任何修改必须显式登记【记录】并说明理由。
3. **三级表述纪律**：L1 事件可写"调控"；L2 只写"结合证据一致"；L3 只写"预测相关"（§3 M3 表）。
4. **禁止合成数值**：不做跨基因 MPRA 活性线性相加；evidence class 为可枚举判定，无权重。
5. **负结果照登**：事件数不足、方向不一致、某层证据为空——全部写入报告与【记录】，不得静默删除。
6. **GSE286075 只作 endfoot 响应参考（2026-09-19 审计修订）**：本地与 GEO supplementary 仅 gene counts（原始 reads 在 SRA，本阶段不下载）；禁止任何 APA 定量；其 regulator mRNA 变化不得用作核内 APA regulator 活性的独立验证层。
7. **进程处置红线（2026-09-19 用户明令，三次同类事故后永久生效）**：杀/重启任何进程前必答三问——① **全家福**：`ps -o pid,ppid,lstart,time,stat,cmd` 看全进程树并标注每个 PID 角色（主/加载器/计算 worker；主与加载器正常即 S+ 睡眠，不构成异常证据）；② **卡死双证据**：TIME 列不增长 **且** 输出文件字节数不增长，同时成立、间隔 ≥10 分钟两次观察；③ **死后归属**：该 PID 死后谁组装输出、谁写 exit 行。重跑前必须确认同输出目录/tmp 的旧实例已全部退出（单实例铁律）。配套：Bash 一律 `cd /d/stroke_apa` + 绝对路径；WSL 路径从真实 `ls` 输出复制（Linux 区分大小写，今日 stroke_apa_work 大小写 typo 致 exists/glob 静默为空）。

## 附 B：下载任务与链接（默认由助手执行；助手失败时用户手动下载到此对应目录）

| 触发条件 | 内容 | 放置目录 | 链接 |
|----------|------|----------|------|
| 已完成 | GSE174574_RAW.tar（251MB） | `D:\stroke_apa\data\GSE174574\` | <https://ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174574/suppl/GSE174574_RAW.tar> |
| 已完成 | GSE286075_RAW.tar（1.8MB） | `D:\stroke_apa\data\GSE286075\` | <https://ftp.ncbi.nlm.nih.gov/geo/series/GSE286nnn/GSE286075/suppl/GSE286075_RAW.tar> |
| P3（T3.4 前；小文件，可在 P2 下载等待期预取） | GSE74456 GEO 补充（PAPTRAP processed） | `D:\stroke_apa\data\GSE74456\` | <https://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74456/suppl/> |
| P3（已完成） | QKI CLIP peaks = **GSE147119**（Sakers 2021，经 Europe PMC 取 PMC7943582 Supplementary；GSE146935 = CRISPR-TRAPseq 非 CLIP——2026-09-19 附件勘误） | `D:\stroke_apa\data\GSE147119\` + `results\p3_level2_qki_clip_peaks.bed` | <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147119> |
| P4（T4.4 前；小文件，可在 P2 下载等待期预取） | GSE330741 MPRA（GEO supp TXT） | `D:\stroke_apa\data\GSE330741\` | <https://ftp.ncbi.nlm.nih.gov/geo/series/GSE330nnn/GSE330741/suppl/> |
| P2-05 执行中 | GSE238125 SRA（12 runs；FASTQ 47.7GB / .sra 镜像约 34GB） | `D:\stroke_apa_data\`（2026-09-16 用户指定，原 E:\ 路径作废） | ENA filereport 已存 D:\stroke_apa_data\ena_filereport.tsv · SRA: <https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&uid_from=PRJNA997998> · 备选 NCBI S3 镜像模式：`https://sra-pub-run-odp.s3.amazonaws.com/sra/<SRR>/<SRR>`（.sra，已验证 200 OK + Range；需 WSL sra-tools 转 FASTQ） |
| P2-07 触发时 | GSE225110 SRA（132 runs，大） | `D:\stroke_apa\data\GSE225110\`（试点已有 6 runs 处理资产在 ~/stroke_APA_pilot，优先复用） | <https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&uid_from=PRJNA922165> |
| P4 参考集 3 触发时（小文件，可在 P2 下载等待期预取） | GSE263986（优先 15MB processed xlsx） | `D:\stroke_apa\data\GSE263986\` | <https://ftp.ncbi.nlm.nih.gov/geo/series/GSE263nnn/GSE263986/suppl/GSE263986_ProcessedData.xlsx> |

### 附 B-2：P2-05h 剩余 8 runs 手动下载清单（2026-09-17 交付用户；ENA 直链，每文件 3.2–4.8GB，共 ~32GB）

**放置目录：`D:\stroke_apa_data\fastq\`，下载后保留 SRR 原名即可，改名与 md5 校验由助手完成。**

| 样本 | SRR | R1（_1.fastq.gz） | R2（_2.fastq.gz） |
|------|-----|--------------------|--------------------|
| day3_rep1 | SRR25403062 | <https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/062/SRR25403062/SRR25403062_1.fastq.gz> | <https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/062/SRR25403062/SRR25403062_2.fastq.gz> |
| day3_rep2 | SRR25403063 | <https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/063/SRR25403063/SRR25403063_1.fastq.gz> | <https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/063/SRR25403063/SRR25403063_2.fastq.gz> |
| day7_rep1 | SRR25403064 | <https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/064/SRR25403064/SRR25403064_1.fastq.gz> | <https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/064/SRR25403064/SRR25403064_2.fastq.gz> |
| day7_rep2 | SRR25403065 | <https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/065/SRR25403065/SRR25403065_1.fastq.gz> | <https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/065/SRR25403065/SRR25403065_2.fastq.gz> |
| day21_rep1 | SRR25403066 | <https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/066/SRR25403066/SRR25403066_1.fastq.gz> | <https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/066/SRR25403066/SRR25403066_2.fastq.gz> |
| day21_rep2 | SRR25403067 | <https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/067/SRR25403067/SRR25403067_1.fastq.gz> | <https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/067/SRR25403067/SRR25403067_2.fastq.gz> |
| day60_rep1 | SRR25403068 | <https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/068/SRR25403068/SRR25403068_1.fastq.gz> | <https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/068/SRR25403068/SRR25403068_2.fastq.gz> |
| day60_rep2 | SRR25403069 | <https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/069/SRR25403069/SRR25403069_1.fastq.gz> | <https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/069/SRR25403069/SRR25403069_2.fastq.gz> |

> 提示：浏览器单线程下载 ENA 可能同样慢（~0.2–0.6MB/s）；如装了 IDM/FDM 等多线程下载器会快很多。全部到货后告知助手，改名+md5 校验+全量比对自动接力。
