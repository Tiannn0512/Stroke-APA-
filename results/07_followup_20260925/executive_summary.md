# executive_summary.md — 一页中文结论（2026-09-25）

**任务**：05 号任务书（第一环节后续分析）。结论针对三问；证据类别标注在括号里。

## 决策

**Atp2a2 Gate 1 = advance to independent 3'RACE**（修订 06 包，此前暂记 HOLD）。
05 号任务书所列三项内部冲突全部核实为真，并新增发现两项更严重缺陷（RACE 引物方向装反、长接头对背向无产物）。五项缺陷本轮全部修复并逐碱基机器验证。`advance` 仅表示检测方案可以交实验负责人终审并启动独立 3'RACE，**不表示 APA 已成立**（计算推断）。

## 三问答案

**① Atp2a2 是否已具备开展独立 3'RACE 的检测条件？**
具备。最终 oligo 集（assay_design_verified.tsv）：outer GSP v2 = `GGTTTCCTGAGGCTTTAATTCCTGTC`、nested GSP v2 = `CGGTCCAAGAGTCTCCTTCTACCA`（Gm30970 最小错配 12、位于共有干净块）、L 特异对 juncL_v2（cDNA 104bp 仅 L）、S 剪接型对 J_S+distal_F（536bp 仅 S）、正交 J_L+distal_F（536bp 仅 L）、总量 qPCR_total（5 转录本 154bp）、终端外显子对 distal（136bp L+S）。全部在 mm10/vM25 转录本+基因组+假基因三空间重扫（≤1 错配），4 条 06 版问题引物留档 deprecated。订单需实验负责人终审（状态=primer-blast checked，非 validated）。

**② 目前公共数据究竟支持哪种 RNA 加工解释？**
支持"卒中伴随 **L 型（spliced-UTR 长版本 179939）成熟分子比例下降**"这一读段结构变化（公共数据直接观察）：L 独有接头读段 sham 96/136→day3 4/11→day7 2/7（≈-96%），远端/共有覆盖比 0.30→0.085，而 pre-mRNA 指示段稳定（排除转录塌缩）。**但无法区分**这是 3' 端选择改变（APA）还是 L 的剪接/加工改变——此外发现 DaPars2 的"近端 PAS 122456498"无 polyA 信号支持（预测窗口边界伪影），真实近端端应重锚定到注释端 122456339（PolyASite sig 48、7 数据集、AATAAA@-31）；L/S 终端端有强支持（sig 332、9 数据集）。普通 RNA-seq 的覆盖比、接头读段、Salmon 比例同属一个发现队列，不构成独立验证；独立异构体层面公共数据不存在（已系统检索，缺口如实记录）。

**③ 若第一环节实验确认，PAP 定位实验应检测哪个定义的 RNA 版本？**
主检测 **L 版本**（ENSMUST00000179939.7 结构）：探针用 L 独有段 122455892-122456339（447bp 多 Z 对）或跨 L 接头序列（J_L 已设计验证）；配对检测 **S 剪接型**（唯一接头 donor122457302→acceptor122454304，J_S 探针）；总量探针用共有外显子。M（31423）无独有序列/接头，只能以 common−L−S 推断（需检出效率校正）。主要检验：卒中(day7)×区室(PAP vs 胞体)对 L 占比的交互作用，动物为生物学单位（spatial_analysis_plan.md 预定义）。

## 交付与缺口

全部交付物在本目录（16 类文件 + analysis_scripts + figures + notes）：结构表、读段证据表（180 行，全数复现 06 包关键数字）、PAS 证据表、verified 引物表+产物表、决策文档、数据集筛查（18 条）、GSE74456/143531/330741 复核（PARTIAL / 无差异 / NOT_APPLICABLE——GSE330741 只 tile Glt1/Sparc，**零 Atp2a2 元素，不得解读为 Atp2a2 验证**）、新颖性（31 条，四个关键环节均空白）、三份预定义方案、元数据模板、分析框架代码。
Pending：独立 3'RACE 数据（Gate B）；GSE330741 manifest 重建待作者补充表；Primer-BLAST 为 GRCm39 辅助性质，mm10/vM25 本地核查为主。
