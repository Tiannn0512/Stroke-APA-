# candidate_decision.md — Atp2a2 Gate 1 决策（v3，修订版）

日期：2026-09-25 ｜ 执行：干实验复核（05 号任务书"第一环节后续分析"） ｜ 参考版本：mm10/GRCm38 + GENCODE vM25
证据文件目录：`results/07_followup_20260925/`（本轮全部产物）；上一轮：`06_gate1_reaudit/`（GATE1_REAUDIT_20260924.zip）

---

## 决策：**advance to independent 3'RACE**

`advance` 仅表示：检测方案内部一致性已逐项机器验证、实验负责人可以审核并订购引物、启动独立 3'RACE。
**不表示** Atp2a2 APA 已成立（位点状态仍为 `predicted/annotated-supported`，未经 RNA–poly(A) 接合实测）。

---

## 1. 对 05 号任务书所列 HOLD 项的逐项复核结果

05 号任务书把 Gate 1 暂记 HOLD，理由是三项内部冲突。本轮逐项核实，**三项全部实锤**，并**新增发现两项更严重的缺陷**：

| # | 缺陷 | 核实结论 | 根因 | 处置 | 验证 |
|---|---|---|---|---|---|
| H1 | `long_junction` 配对扫描预测 0 扩增子 | 实锤（primer_specificity_hits.tsv 确为 0） | 06 版设计模板把终端外显子段拼在倒数第二外显子段**之前**（旋转模板）；两条引物在真实 mRNA 上**背向**，无任何扩增子 | 重设计为 juncL_v2：F 在 L 独有段(122455892-122456339)、R 在终端外显子 | cDNA 产物 **104bp 仅 179939**；基因组产物 1692bp（跨 1588bp 内含子）可区分（primer_pair_products.tsv） |
| H2 | outer GSP 表/文档序列不一致 | 实锤：表=`GACAGGAA...`；文档=`TCAAGAAG...` | 文档泄漏了设计早期候选（其真实位置 chr5:122457534-122457560，非声称的 122459602） | 勘误文档；canonical 以机器验证表为准 | f0_outer_gsp_check.py：`GACAGGAA...` 在 chr5 全染色体**仅 1 处精确命中** |
| H3 | nested GSP 有 Gm30970 命中 | 实锤：`TGAGGG...` 在 Gm30970 转录本 44 位**精确命中** | Gm30970(ENSMUST00000227710.1, 363bp) = L mRNA **3157-3519 的近乎完美拷贝**（滑窗作图，见 f4_design_log），旧引物位置(mRNA 3200)恰在其中 | 在共有干净块重设计 | nested_v2=`CGGTCCAAGAGTCTCCTTCTACCA`，Gm30970 最小错配 **12**，L/S/M 各精确单命中 |
| H4 | （本轮新发现）RACE outer/nested GSP **方向装反** | 实锤：两条 GSP 均为 mRNA 反义向（向 mRNA 5' 延伸） | 06 版在 mRNA-sense 模板上取了 primer3 的 RIGHT 引物；3'RACE 的 GSP 必须正义向（朝 poly(A) 延伸） | outer_v2 = v1 的反向互补（Tm/GC 不变，脱靶重扫通过）；nested 重设计 | 方向列（extend_dir）与逐碱基定位见 assay_design_verified.tsv |
| H5 | （本轮新发现）"366bp cDNA 产物"不成立 | 实锤：v1 对在真实 mRNA 上无产物；设计表/文档的产物长度互相矛盾 | H1 同根因 | 由 juncL_v2 取代（104bp cDNA / 1692bp 基因组） | primer_pair_products.tsv |

**处置后的 oligo 集**（assay_design_verified.tsv，14 行；4 条 v1 行仅留档标注 `deprecated`，不得订购）：

| 用途 | oligo | 序列(5'→3') | 关键验证 |
|---|---|---|---|
| 3'RACE outer GSP | RACE_outer_GSP_v2 | GGTTTCCTGAGGCTTTAATTCCTGTC | L/S/M 各 1 命中(cDNA 2831/2552/2882)；Gm30970 minMM=13；chr5 唯一 |
| 3'RACE nested GSP | RACE_nested_GSP_v2 | CGGTCCAAGAGTCTCCTTCTACCA | L/S/M 各 1 命中(3069/2790/3120)；Gm30970 minMM=12；位于共有干净块 122458451-475 |
| 总量 qPCR（5 转录本） | qPCR_total_F/R | AGCCTTTGTAGAGCCGTTTGTA / CACACTCTTTCTGTCCTGTCGA | 154bp × 5 条注释转录本；基因组产物 2458bp 可区分 |
| 终端外显子 qPCR（L+S） | distal_terminal_F/R | AGTTAGGACTGGAGGCCTATGT / TTGTAAGTGGCCAGATTGCTCT | cDNA 136bp(L+S)；基因组 136bp |
| L 特异 qPCR | long_junction_v2_F/R | TGTGGTGTTTTCCTCCAATGCCT / CACGCACCCGAACACCCTTATAT | cDNA 104bp **仅 L**；基因组 1692bp |
| S 剪接型 qPCR | qPCR_S_spliceform_J(+distal_F) | TGGAACAACCCGCAATACTGGAGT | 跨 S 独有接头；cDNA 536bp **仅 S** |
| L 正交 qPCR | qPCR_L_orthogonal_J(+distal_F) | CCATCAACTAACCAATACTGGAGT | 跨 L 独有接头；cDNA 536bp **仅 L**；可作 smFISH 探针基础 |

短版本特异检测**可以设计**（推翻 06 版"短版无独有序列"的前提——那是基于错误的"共有外显子"注解；实际上 S 有唯一的剪接接头 donor 122457302→acceptor 122454304，见下节）。但其确认对象是 **S 剪接型**，S 的 3' 端身份仍待 3'RACE 实测。Primer-BLAST 补跑（GRCm39 RefSeq，辅助性质）：juncL_v2 对在 intended targets 产物 1692bp 与本地一致；两个跨接头对 "No target templates" 属预期（接头序列参考中不存在）；NCBI 已不接受脚本化单引物提交（存证见 primer_blast_v2/README.md），两个 GSP 的等效近似匹配=本地三空间 ≤1 错配扫描，残余缺口（基因座外全基因组近似匹配）由实验负责人下单前手动补一次网页版 Primer-BLAST。**主判定为 mm10/vM25 本地逐碱基 + ≤1 错配扫描**（f4）。

---

## 2. RNA 版本定义（本项目口径，供 Gate B/C 使用）

逐碱基覆盖唯一性分段（Atp2a2_transcript_structure.tsv，chr5:122453000-122460000）：

| 版本 | 转录本 | 结构 | 可识别标识 | 3'端状态 |
|---|---|---|---|---|
| **L（长 UTR，事件载体）** | ENSMUST00000179939.7 (Atp2a2-203) | UTR 跨 1588bp 内含子（122454304-122455892），倒数第二外显子 UTR + 终端外显子 UTR | ① 独有序列 122455892-122456339（447bp）；② 独有接头 donor122455892→acceptor122454304 | PolyASite 强簇支持 + AATAAA@-32（外部支持，未实测接合） |
| **S（短 UTR）** | ENSMUST00000177974.7 (Atp2a2-202) | 终端外显子内 780bp UTR | **唯一接头** donor122457302→acceptor122454304（无独有序列，J_S 跨接头识别） | 同上（与 L 共簇） |
| **M（中 UTR）** | ENSMUST00000031423.9 (Atp2a2-201) | 倒数第二外显子内 814bp UTR，无终端外显子 | **无独有序列、无独有接头**；只能靠 3' 端位置(122456339)识别 | PolyASite 簇 122456301-355 (sig 48, n=7) + AATAAA@-31 |

另有两个次要转录本（ENSMUST00000196490.1/197415.4，698/615bp，3'端在 1224889xx），qPCR_total 会一并扩增——解读总量时注明。

**对 06 版结构描述的修正**：06 版称 "122457302-457426 为三转录本共有外显子"。GTF 实况：三个转录本共享 acceptor 122457426，差异在**供体/延伸**——177974 的倒数第二外显子止于 122457302，31423 止于 122456339（其末端即 3' 端），179939 延伸至 122455892 并经 1588bp 内含子接终端外显子。

---

## 3. 局部读段证据（6 个 BAM 切片复算，Atp2a2_local_read_evidence.tsv）

发现队列（GSE238125；sham 2 / day3 2 / day7 2；多个时间点共用 sham，非独立重复）：

- 远端/共有区覆盖比：sham 0.334/0.271 → day3 0.086/0.084 → day7 0.086/0.085（ruleB_depth20：MAPQ≥20、baseQ≥20）
- L 独有接头（122454304-122455892）支持读段：sham 96/136 → day3 4/11 → day7 2/7（≈-96%，归一化）
- S 独有接头（122454304-122457302）：sham 13/18 → day3 14/10 → day7 8/15（相对稳定）
- 敏感性：排除 L 接头读段后，sham 比值 0.334→0.311、0.271→0.251；day3/day7 基本不变（0.086→0.086 等）→ 覆盖差异**不能**被 L 接头读段单独解释（此操作不是 APA 定量分解）
- 基因 3' 端外侧 500bp intergenic 对照：0-1.5×（背景干净）
- 读段过滤规则逐条写入 TSV 的 filter_rule 列；图区分真 IGV 截图（06/igv/ 15 张）与程序深度图（06/depth_plot_*.png）

**这些全部来自同一发现队列，只能算一条证据线**（公共数据直接观察），不能互相充当独立验证。

## 4. PAS 证据与重锚定（PAS_evidence_table.tsv）

- **DaPars2 预测近端 PAS chr5:122456498 无信号基序支持**（PolyASite 该处仅噪声簇 sig 0.029/0.011、无 motif；本表 computed scan: none）→ 判定为 **预测窗口边界伪影**
- **重锚定**：近端 3' 端应为注释的 **122456339**（M 端）——PolyASite 簇 122456301-355（sig 48.06，7 数据集，分位 0.74）+ AATAAA@-31@122456357，预测断裂 ~122456327±，与注释端吻合
- L/S 终端端 122453512/513：PolyASite 强簇 122453501-534（**sig 332.5，9 数据集，分位 0.94**）+ AATAAA@-32@122453540 ✓
- 内部引发风险：L/S 端 LOW（下游 A20%=20%）；M 端 MEDIUM（A50%、run5）→ 写入 3'RACE 检查项
- 另发现终端外显子上游数个次级簇（122453580/122453643/122453690/122453873）与一处 ATTAAA 次级位点（122456938）——列为候选次要 PAS，待 3'RACE 判别

## 5. 公共数据究竟支持哪种 RNA 加工解释

三个候选解释——(i) APA（3' 端选择改变）、(ii) 剪接形式改变（L 型 spliced-UTR 剪接减少）、(iii) 终端外显子使用/加工改变：

- 支持"**L 型成熟分子减少**"的证据：L 独有接头读段 -96%（该接头只存在于 L 成熟 mRNA，pre-mRNA 不含）；L 独有段(447bp)覆盖下降；终端外显子覆盖下降
- **反对"转录塌缩"**：spliced-UTR 内含子区(pre-mRNA 指示)覆盖卒中后不降反稳（0.067/0.071→0.085/0.081）；共有外显子区稳定
- **区分(i)与(ii)是当前公共数据做不到的**：普通 RNA-seq 无法证明终端外显子 3' 端发生了 poly(A) 接合。故结论表述为：**"卒中伴随 Atp2a2 L 型（spliced-UTR 长版本）成熟分子比例下降"（计算推断），至于是 L 的 3' 端选择改变（APA）还是 L 的剪接/加工改变，待独立 3'RACE 判别**——这正是 Gate B 的目的。

## 6. 独立公共证据现状（本轮三个智能体交付，详见 notes/）

- **异构体层面无独立卒中队列**（public_dataset_suitability.tsv：18 条逐一核对；卒中+星形胶质细胞+3'端/长读段的组合不存在）→ 缺口如实记录，不制造复现
- GSE74456（Sakers 2017 PAPTRAP）：**PARTIAL**——历史值 log2FC 2.13/padj 8e-7 在 GEO counts、series matrix 与 PNAS S1-S4 中均不存在；我们复算 +2.03（padj 0.269）；且按作者定义 Atp2a2 **不是** PAP 富集（PAP TRAP vs SN input = -0.69）
- GSE143531（Mazaré 2020）：**未检出差异**——作者原表 mmc2 行 1175：log2FC +0.3913、padj 0.7304；Atp2a2 不在 PAP-enriched 名单
- GSE330741（Koester/Sakers/Dougherty 2026 AAV9-MPRA）：**NOT_APPLICABLE**——只 tile Glt1/Sparc，全部 6430 元素中 Atp2a2 命中为 0；**不得**解读为 Atp2a2 验证
- 基因层面最佳独立队列：GSE225110（Aldh1l1-RiboTag 卒中，n=3-12/组）等 3 条，均 GENE_LEVEL_ONLY
- PAS 佐证资源（非卒中）：GSE94054 (cTag-PAPERCLIP)、GSE171312（海马 direct RNA）等
- 新颖性（novelty_overlap.tsv，31 条）：四个关键环节（卒中→Atp2a2 APA、异构体 PAP 分布、卒中改变分布、QKI 调控 Atp2a2 APA）**均为空白**；最近邻 = Shim 2025（MCAO 改变星形胶质细胞 endfeet 翻译组）与 Sakers/Dougherty 组 SN-MPRA 平台（需跟踪预印本）

## 7. Gate 状态与条件

| 关口 | 本轮状态 |
|---|---|
| A 检测方案 | **通过（本轮）**：引物方向/产物/参考版本/脱靶核查一致；待实验负责人终审后可订购 |
| B 真实 3' 端 | 未开始（3'RACE 之后） |
| C 空间分布 | 未开始 |
| D 作用机制 | 未开始 |

`advance to independent 3'RACE` 的附带要求：① 3'RACE 测序确认 RNA–poly(A) 接合（凝胶/单 qPCR 不够）；② no-RT/无模板/内部引发对照（M 端 A50% 中风险必查）；③ S 版本身份以 3'RACE 结果为准。
