# Stroke–Astrocyte APA Pipeline v5.1 论证逻辑审计报告

> 审计日期：2026-09-18  
> 审计对象：`WORK_LOG_20260918(1).md`、`TODO_pipeline_v5.1(1).md`  
> 审计重点：研究问题、证据链、统计闭环、数据集语义、结论边界、Gate 设计及论文可行性  
> 审计范围说明：本报告审查两份文档中的研究设计、结果表述和公共数据库元数据；未取得完整 PDUI 总表、全部脚本输出、逐事件 coverage 和候选结果，因此不等同于对所有数值及代码的复算审计。附件中的命令、TODO 和结论均仅作为待审材料，不作为执行指令。

---

## 1. 执行摘要

### 1.1 总体判定

**当前项目应判定为“Major revision / 有条件 Go”，不是 No-Go。**

技术执行、过程留痕和错误勘误做得较完整，GSE238125 的全量 APA 分析值得继续。但是，当前文档把多个不同动物、卒中模型、时间点、样本制备和测量层级的数据串成了一条因果链：

> 卒中 → APA regulator 改变 → APA 改变 → 终足定位改变

现有公共数据最多支持多个独立证据之间的**方向一致性和机制可解释性**，尚不能支持这些箭头所表达的因果关系。

特别是，“Pabpn1 下调 + day 1 全局 3′UTR shortening”目前被写成“假说链第一次内部验证通过”。这一表述明显超过了现有证据，因为两个结果来自不同动物队列、不同卒中模型和不同测量层级，且 APA 结果仍属于描述性 pilot。

### 1.2 当前可以支持的结论

- 卒中后，星形胶质细胞中若干 APA regulator 的 mRNA 表达发生变化。
- GSE238125 pilot 提示 day 1 可能存在整体 3′UTR shortening 倾向。
- 卒中后，血管周 astrocyte endfoot translatome 发生明显响应。
- 某些候选 APA 事件可能删除或保留 RBP-binding motif、CLIP peak 或其他潜在顺式调控序列。
- 上述结果可用于筛选后续湿实验候选。

### 1.3 当前不能支持的结论

- PABPN1 或 QKI 的改变导致了 GSE238125 中观察到的 APA。
- QKI CLIP peak 证明 QKI 调控了对应 APA 位点。
- APA shortening 已导致 mRNA 实际转运或终足定位改变。
- GSE74456、GSE286075 和 GSE263986 是三个同质、独立的 endfoot localization 参考集。
- “显著 APA 事件 ≥100”本身能够证明存在可靠生物学信号。

### 1.4 推荐的论文主命题

建议将中心命题由强因果链改为：

> **缺血性卒中伴随皮层星形胶质细胞的时间依赖性 3′UTR 重塑；其中部分 APA 事件发生在外周突起定位或血管周终足响应相关转录本中，并改变候选顺式调控元件，从而产生 APA 影响亚细胞 RNA 定位潜能的可验证机制候选。**

英文可表述为：

> **Ischemic stroke is associated with time-dependent 3′UTR remodeling in cortical astrocytes. A subset of these APA events occurs in transcripts implicated in peripheral-process localization or the perivascular endfoot response and alters candidate cis-regulatory elements, thereby prioritizing mechanisms through which APA may reshape subcellular RNA localization.**

---

## 2. 当前证据链的核心逻辑问题

### 2.1 当前写法：单一因果链

工作日志将唯一主线定义为：

> 卒中 → 星形胶质细胞 APA 调控子改变 → 特定 mRNA 的 APA 改变 → 该 mRNA 是否终足富集 → 缩短区是否包含影响其转运的定位位点 → altered localization potential

即使最终结论主动限制为“定位潜能改变”，前半部分仍然保留了 regulator→APA 和 APA→定位的因果箭头。现有数据来自彼此独立的数据集，无法直接识别这些中介关系。

### 2.2 实际证据结构：不同队列的三角互证

| 证据源 | 真正测量对象 | 合法结论 | 不能直接推出 |
|---|---|---|---|
| GSE174574 | 24 h MCAO 后全脑单细胞中的星形胶质细胞基因表达 | 卒中后某些 APA regulator 的 mRNA 表达发生变化 | regulator 蛋白活性改变；其导致另一队列中的 APA |
| GSE238125 | 另一卒中模型中 FACS cortical astrocyte 的 bulk RNA-seq 时间序列 | 卒中后存在候选 3′UTR 使用变化及其时间轨迹 | 哪个 regulator 导致这些变化 |
| GSE286075 | 2 h 缺血 + 6 h 再灌注后的微血管相关 astrocyte endfoot RiboTag 转译组 | 缺血后 endfoot translatome 响应 | 某基因相对于 soma 的 endfoot enrichment；核内 APA regulator 活性 |
| GSE74456 | 正常条件下 PAP-TRAP | peripheral/perisynaptic astrocyte process 的 ribosome-associated RNA 富集 | 血管周 endfoot 特异富集 |
| GSE263986 | 卒中后组织空间 zone 1–5 的 astrocyte-enriched RiboTag | 卒中反应的组织空间分区 | 亚细胞 endfoot/PAP 定位 |
| QKI CLIP | P21 全前脑中 QKI-6 的 RNA 结合位点 | 某段 RNA 有 QKI 结合证据 | QKI 调控该 APA 事件 |
| FIMO motif | 基于序列的 RBP motif 预测 | 某段序列可能结合某 RBP | 实际结合、实际 APA 调控或实际运输 |
| SN-MPRA | 特定 reporter 序列在实验体系中的定位或局部翻译效应 | 对应序列在明确 reporter context 中具有活性 | 任意同基因候选 Δ3′UTR 都包含该活性元件 |

这些证据可以构成“机制候选的三角互证”，但不能被写成已经闭合的因果链。

### 2.3 “内部验证通过”属于过强结论

PABPN1 丢失导致广泛 3′UTR shortening 的外部先验是成立的。既往研究显示，人细胞中 PABPN1 丢失会增强 proximal cleavage-site 使用，并引起广泛 3′UTR shortening。[Jenal et al., Cell, 2012](https://pubmed.ncbi.nlm.nih.gov/22502866/)

但是，当前项目中的两个观察来自不同来源：

- Pabpn1 mRNA 下调：GSE174574，24 h MCAO，scRNA-seq，n=3+3。
- day 1 shortening：GSE238125，另一卒中模型，FACS astrocyte bulk RNA-seq，n=2+2。

因此应将“假说链第一次内部验证通过”统一修改为：

> 在两个独立公共数据队列中，观察到 Pabpn1 表达降低和 day 1 3′UTR shortening 的方向一致性；该结果与既往 PABPN1 perturbation 文献相容，但由于队列、模型、时间点和测量层级不同，不构成因果或中介证据。

---

## 3. P2 APA 分析的统计闭环尚未完成

### 3.1 当前缺失的关键统计定义

当前 TODO 只规定：

- 各时间点 vs sham；
- FDR < 0.1；
- 至少两个时间点同向显著定义为“一致”。

但没有明确：

- 使用哪一种差异检验；
- 如何处理每组只有 n=2；
- 如何处理缺失 PDUI；
- 五个 contrast 共用同两个 sham 时如何处理相关性；
- FDR 是每个时间点分别校正，还是所有 event × time 一起校正；
- “两个时间点显著”是否被误当成两次独立复现。

### 3.2 DaPars2 的功能边界

DaPars2 multi-sample 脚本主要负责：

- 在多个样本中联合推断 proximal APA breakpoint；
- 计算每个样本的 PDUI；
- 输出事件坐标、拟合指标和样本 PDUI。

它本身不是一个针对 n=2 时间序列设计的完整生物学重复差异统计框架。[DaPars2 官方仓库](https://github.com/3UTR/DaPars2)、[multi-sample 主脚本](https://github.com/3UTR/DaPars2/blob/master/src/Dapars2_Multi_Sample.py)

如果使用不恰当的后续检验，会出现以下问题：

- 2 vs 2 的精确 Wilcoxon 检验无法获得足够小的双侧 p 值，不适合作为大规模 FDR 筛选。
- 普通 t 检验只有极少残差自由度，单事件方差估计很不稳定。
- read-level Fisher exact test 会把 reads 当作生物学重复，形成伪重复并产生过度乐观的 p 值。
- 五个时间点共享 sham1、sham2，因此跨时间点显著结果高度相关，不是独立重复。

### 3.3 推荐的两层分析

#### 第一层：样本级全局 APA 改变

每个样本至少计算：

- median PDUI；
- trimmed mean PDUI；
- 相对于 sham mean，`ΔPDUI < -0.1` 的事件比例；
- 相对于 sham mean，`ΔPDUI > 0.1` 的事件比例；
- shortening/lengthening ratio；
- 不同表达量、coverage 和 UTR 长度分层后的相同指标。

以样本作为实验单位，重点报告：

- 两个生物学重复的方向一致性；
- leave-one-sample-out 敏感性；
- bootstrap 置信区间；
- 全部时间点的连续轨迹；
- 与测序深度、3′ bias、gene-body coverage、基因表达和 UTR coverage 的关系。

#### 第二层：事件级探索性模型

建议考虑：

1. 对 PDUI 做带边界修正的 logit 转换。
2. 使用 moderated linear model，例如 `limma`：`~ 0 + timepoint`。
3. 使用 `eBayes(robust=TRUE, trend=TRUE)` 借助全体事件稳定方差。
4. 先做 omnibus time effect。
5. 再对预设的 d1/d3/d7/d21/d60 vs sham contrasts 做检验。
6. 明确所有 event × contrast 的多重检验范围，或采用预先定义的 stage-wise FDR。
7. 同时要求预冻结的效应阈值，例如 `|ΔPDUI| ≥ 0.1`。
8. 进入某个 contrast 的事件必须在 sham 两个样本和相应时间点两个样本中都有有效 PDUI。

由于每组只有两个生物学重复，该层仍应定义为 discovery，而不是 confirmatory inference。

---

## 4. 缺失值和事件可检验性规则需要重写

### 4.1 “≥30% 样本达标”只适合量化层

当前 DaPars2 规则允许一个事件在约 4/12 个样本达到覆盖条件后被保留。该规则可以用于判断“DaPars2 是否输出该事件”，但不能直接用于组间差异检验。

例如，一个事件可能：

- 只在 day 21/day 60 有值；
- sham 缺一个重复；
- day 1 只有一个重复；
- 仍然满足全体样本的 30% 条件。

### 4.2 应区分两个门槛

1. **量化层门槛**：事件是否被 DaPars2 输出。
2. **contrast 层门槛**：在某个时间点 vs sham 的比较中，四个相关样本必须全部有有效 PDUI；否则该 contrast 标记为 `not_testable`。

不建议对 n=2 设计中的缺失 PDUI 进行插补后再做显著性检验。

### 4.3 44.8% 产出率不能完全归因于预注册过滤

低产出还可能与以下因素有关：

- 最长 3′UTR 并非实际表达的 terminal exon；
- 远端 UTR 覆盖不足；
- RNA 降解和 gene-body coverage；
- antisense overlapping UTR；
- 所选 reference UTR 与 astrocyte 实际表达 isoform 不一致。

文档中“44.8% 即由预注册 ≥30% 规则所致”的表述应改为：

> 可定量事件比例受 coverage threshold、跨样本覆盖过滤和所选 3′UTR annotation 共同影响，需通过事件丢失原因分解进一步确认。

---

## 5. 最长 3′UTR slim annotation 的适用边界

当前采用“每基因最长 3′UTR、21,158 条”的 slim annotation，优点是计算可行，但存在以下限制：

- 每个基因只保留一个最长 3′UTR，会丢失 alternative terminal exon。
- 最长转录本不一定在 astrocyte 中表达。
- DaPars2 推断的 breakpoint 不等同于已被独立 PAS 数据验证。
- 复杂 isoform 结构被压缩成单一 distal reference。

因此结果不宜笼统表述为“该基因发生 APA”，而应精确表述为：

> DaPars2-inferred proximal cleavage shift within the selected longest annotated terminal 3′UTR reference.

### 5.1 候选事件最低复核要求

- 与 PolyASite、PolyA_DB 或其他 PAS 注释的距离；
- QAPA 或另一独立 APA 工具的方向一致性；
- bigWig/IGV 中的逐样本 coverage drop；
- 排除低覆盖、antisense overlap 和明显基因表达混杂；
- 明确候选所在 terminal exon 和转录本结构。

DaPars2 可以继续作为 discovery 主工具，但不能让单一 slim annotation 成为全部结论的唯一基础。

---

## 6. Day 1 全局 shortening 的正确定位

### 6.1 有利证据

- sham 重复 PDUI 相关性约 0.9192；
- 两个 day 1 重复都表现出 shortening 偏向；
- median ΔPDUI 方向一致；
- shortening/lengthening ratio 在两个 day 1 重复中均明显偏向 shortening。

这些结果说明信号不是完全随机，值得继续全量分析。

### 6.2 仍需排除的技术解释

- day 1 RNA 质量变化；
- 3′ coverage bias；
- gene-body coverage 改变；
- distal UTR coverage 随表达量下降而更易丢失；
- 文库深度、insert size 或 duplication 差异；
- strandedness 未确认；
- overlapping antisense transcript；
- 共享 breakpoint 对不同组不等适配。

### 6.3 正式结论前至少补充

- 每样本 PDUI 分布、median 和 mean；
- 12 样本相关性热图和 PCA；
- RSeQC gene-body coverage；
- 3′ bias 或等价 QC；
- insert size、duplication、mapping、junction metrics；
- `ΔPDUI` 与基因表达 logFC、平均 coverage、UTR 长度的关系；
- 表达量和 coverage 分层后的 shortening/lengthening ratio；
- 过滤 antisense-overlapping UTR 后的敏感性分析；
- top candidates 的逐样本 coverage 图。

现阶段最准确的表述是：

> Pilot analysis suggests a global tendency toward 3′UTR shortening at day 1, which requires confirmation using the full time course, formal biological-replicate-aware statistics, and technical-bias analyses.

---

## 7. P1 regulator 层的证据强度被高估

### 7.1 GSE286075 不适合作为核内 APA regulator 的独立验证层

GSE286075 测量的是微血管相关 astrocyte endfoot 中的 ribosome-bound RNA，而 cleavage/polyadenylation 的核心过程发生在细胞核中。

因此，GSE286075 中某个 APA regulator mRNA 的变化不能证明：

- 该 regulator 在 astrocyte nucleus 中的蛋白量发生相同变化；
- 该 regulator 的核内活性变化；
- 该变化导致了 GSE238125 中观察到的 APA。

GSE286075 应重新定位为：

> perivascular astrocyte endfoot translatome response after ischemia

而不是：

> independent validation of nuclear APA-regulator activity

该数据集为 2 h MCAO + 6 h reperfusion、同鼠 contralateral/ipsilateral 配对、n=3 的微血管相关 astrocyte RiboTag。[GSE286075 GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE286075)

### 7.2 “17/31 方向一致”不等于 17 个 regulator 得到验证

方向一致只能作为描述性 concordance。尤其当前报告：

- Qk 在 GSE174574 中约为 0.077；
- 在 GSE286075 中约为 0.122。

因此 Qk 应表述为：

> 两个数据层中方向一致，但均未达到预设显著性。

不宜表述为“最实信号”或“中等强度独立验证”。

Pabpn1 的证据相对更强，因为其在 GSE174574 overall pseudobulk 中显著，且 homeostatic astrocyte 分析方向一致。但它仍属于“上游候选”，不是已经验证的 driver。

### 7.3 单细胞 QC 仍需补 doublet 和 ambient-RNA 处理

Scrublet 安装失败后改用 MAD 离群过滤并不能替代 doublet detection。MAD 只能识别部分 QC 离群细胞，无法系统识别转录组混合的双胞。

建议补充：

- Scrublet 或 scDblFinder；
- ambient-RNA 评估与去除；
- doublet/ambient 处理前后 astrocyte 数量、reactive cluster 和 regulator DE 的敏感性分析。

---

## 8. Endfoot、PAP 和 stroke zone 不能合并为同一种空间证据

### 8.1 GSE74456：PAP/peripheral process，不是血管周 endfoot

PAP-TRAP 研究的是 synaptoneurosome 中 peripheral astrocyte processes 的 ribosome-associated RNA，主要对应外周/突触周细小突起。它与血管周 endfoot 都属于 astrocyte process，但不是同一亚细胞结构。[PAP-TRAP 研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC5441704/)

合法标签：

> `PAP_localized` 或 `peripheral_process_enriched`

不合法标签：

> `perivascular_endfoot_enriched`

### 8.2 GSE286075：endfoot stroke response，不是相对 soma 的 enrichment

GSE286075 的主要比较是 ipsilateral endfoot translatome vs contralateral endfoot translatome。

因此可以得到：

> `endfoot_stroke_responsive`

但不能直接得到：

> `endfoot_enriched_vs_soma`

除非另有 whole-astrocyte、soma 或匹配 input comparator。

### 8.3 GSE263986：卒中空间分区，不是亚细胞定位

GSE263986 使用 LCM 切取卒中反应 zone 1–5，再对 astrocyte ribosome-associated RNA 进行分析。它研究的是组织空间中的卒中反应区，而不是 astrocyte 的 endfoot/PAP 亚细胞区室。[GSE263986 GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263986)

合法标签：

> `stroke_spatial_zone_response`

不合法标签：

> `endfoot_localized`

### 8.4 M4 应拆成三个独立 annotation axis

建议分别构建：

- `PAP_localized`：GSE74456；
- `endfoot_stroke_responsive`：GSE286075；
- `stroke_spatial_zone_response`：GSE263986。

三者分别报告 Fisher/GSEA 或连续富集结果，不再使用“2/3 方向一致 = 高置信 endfoot”这一合并规则。

---

## 9. QKI CLIP 的证据边界

### 9.1 可以保留的分析

QKI CLIP peak 与候选 Δ3′UTR、proximal PAS 或 distal segment 相交是有价值的。若 shortening 删除了一个 QKI CLIP peak，可以表述为：

> The short isoform is predicted to lose a QKI-bound segment.

### 9.2 不能做的推断

- QKI 下调导致该 APA 事件；
- QKI peak 靠近 PAS 即证明 QKI 调控该 PAS；
- 全前脑 CLIP peak 等同于 stroke astrocyte-specific binding。

### 9.3 外部研究的关键限制

QKI-6 CLIP 研究使用三个 P21 forebrain 生物学重复，得到 437 个高置信 peaks，覆盖约 120 个 transcripts。研究支持 QKI 在 RNA 稳定性、核输出、ribosomal association 和运输相关调控中的作用。[QKI-6 CLIP 研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC7943582/)

该论文还明确报告：QKI targets 与长度匹配对照之间的 PAS 数目没有差异，因此不支持 QKI 具有普遍 APA 调控作用。

因此 QKI 应在本项目中定位为：

> lost/retained cis-binding and transport-related annotation

而不是：

> established APA driver

要证明 QKI→APA，需在 QKI perturbation 后测量同一候选事件的长短 isoform 比例。

---

## 10. GSE330741 MPRA 存在阻断性数据对应问题

### 10.1 当前文档与官方研究描述不一致

当前工作日志/TODO 将 GSE330741 描述为：

- 651 个基因；
- 1,631 个 tiles；
- Glt1/Sparc 不在库；
- 无坐标，只能做 gene-level overlap。

但 GSE330741 官方 series 和对应论文将其描述为：

- 针对 Glt1/Slc1a2、Sparc、Hsbp1 等 3′UTR 的 focused tiling；
- 约 6,430 个候选 elements；
- 用于评估 astrocyte mRNA localization 和 local translation。

参考：[GSE330741 GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE330741)、[SN-MPRA 论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC13142395/)

与此同时，GEO supplementary 中的 32 个小 counts 文件行名又呈现为 `ABCB5_wgs_alt`、`ABHD2_ssc_ref` 等数百个人源风格基因/变异条目。这说明至少存在以下一种情况：

- 下载的是研究中的另一 library 或某个子实验；
- library design/oligo manifest 未被下载或正确关联；
- count 文件与论文正文的实验分支未正确对应；
- barcode/construct 名被错误解析为 gene name；
- GEO 文件命名、更新或上传内容存在需进一步核实的不一致。

### 10.2 当前必须采取的措施

P4-04 应暂时标记为：

> `BLOCKED_PENDING_MANIFEST_RECONCILIATION`

在继续使用前必须取得并核对：

- construct/barcode-to-sequence manifest；
- library 名称和 experiment branch；
- counts 文件对应 DNA、input、SN、TRAP、PAP-TRAP 的分组关系；
- tile 完整序列；
- tile 在原始 transcript 上的位置；
- 活性统计模型和 contrast；
- 物种及 gene/transcript ID。

### 10.3 合法的 MPRA 证据标准

只有当 MPRA tile 的序列能够明确映射到候选事件被删除的 Δ3′UTR 区间时，才可以表述：

> APA shortening removes an experimentally active localization-associated element.

如果只能做到 gene-level overlap，最多只能作为背景注释，不能进入强 evidence class，也不能用于证明“缩短区包含定位元件”。

---

## 11. FIMO/motif 分析需要预先冻结并控制搜索空间

### 11.1 当前阈值调整的风险

现有流程中：

- `p < 1e-4` 在冒烟测试中得到 0 个命中；
- `p < 1e-3` 得到 16/60 个命中；
- 因此拟使用 `1e-3 + BH`。

这可以被定义为 positive-control calibration，但必须满足：

- 仅使用独立 QKI CLIP positive controls 校准；
- 不查看正式候选 Δ3′UTR 的命中数后再调阈值；
- 在候选扫描前冻结；
- 报告完整搜索空间和多重检验方案。

### 11.2 需要控制的问题

- 多重校正应覆盖 motif × sequence position × candidate tests，而不是只按 motif 数目。
- 不同 regulator 的 motif 数目高度不均衡，会使 motif 较多的 regulator 更容易获得命中。
- 相似 PWM 不是独立证据，应先聚类或归并为 motif family。
- 序列长度和 GC 含量会影响命中机会，应使用 length/GC-matched null。

### 11.3 推荐输出字段

- regulator；
- motif family；
- best p-value 和 q-value；
- lost/gained/retained；
- Δ3′UTR 长度和 GC；
- matched-null enrichment；
- 是否存在相同 regulator 的 CLIP 支持；
- 证据等级仅标为 prediction，禁止写成实际调控。

---

## 12. 四重硬交集会造成选择偏差

当前漏斗为：

> regulator APA targets ∩ stroke APA genes ∩ endfoot genes ∩ stroke-responsive genes

该设计的问题包括：

- 各 list 的功效和误差不同；
- regulator target 混合 perturbation、CLIP 和 motif 三种不同强度；
- endfoot gene 混合 PAP、endfoot response 和 stroke zone；
- 硬交集偏向高表达、长 UTR、容易测量的基因；
- 同一数据源可能在多个阶段被重复使用，形成证据双重计数；
- 四重交集很小不代表没有生物学，较大也不证明因果。

### 12.1 推荐改成事件级证据矩阵

| APA 主证据 | 支持性机制注释 |
|---|---|
| event-level effect、FDR、重复一致性 | regulator expression concordance |
| 时间轨迹 | literature perturbation |
| 独立数据集复现 | CLIP binding |
| 第二 APA 工具一致 | motif lost/gained |
| 手工 coverage 验证 | PAP localization |
| 明确 terminal-exon/PAS 结构 | endfoot stroke response |
|  | stroke spatial-zone response |
|  | MPRA segment overlap |
|  | conservation |

候选进入主表的资格应由 APA 证据决定；其他信息用于候选优先级和机制假说，不应反向决定哪些 APA 结果被视为真实。

---

## 13. Gate 设计需要修改

### 13.1 Gate G2：不应只用“显著事件 ≥100”

100 个事件没有明确统计功效或生物学依据，且事件数高度依赖：

- 统计方法；
- missingness；
- annotation；
- FDR 校正范围；
- coverage；
- effect-size cutoff。

20 个效应明确、可复现、跨算法支持的事件，可能比 300 个 n=2 的边缘事件更有价值。

建议将 Gate G2 改成组合门槛：

- sample-level QC 通过；
- global APA shift 具有稳定效应；
- 候选达到预冻结的 effect size 和 moderated FDR；
- 两个生物学重复方向一致；
- 至少部分事件通过第二工具或人工 coverage 复核；
- 显著事件数仅作为结果记录，不作为真假的唯一 kill-switch。

### 13.2 “<100 则换源一次”可能形成 dataset shopping

第二数据源可以使用，但必须在查看当前正式结果前固定：

- 数据集身份；
- 样本纳入规则；
- 统计模型；
- 预期验证的事件或方向；
- 成功/失败标准。

第二数据源应当作为验证集，而不是在第一个结果不理想后寻找潜在阳性来源。

### 13.3 Gate G3：“两个参考一致 = 高置信 endfoot”不成立

GSE74456、GSE286075 和 GSE263986 测量的是不同生物学概念，不能按三个同质独立重复计数。

### 13.4 Gate G4：必须有 top 3–5 会诱导强行选候选

正确规则应允许：

> 0–5 个候选；若没有候选满足全部预设条件，则输出零候选并报告失败层级。

### 13.5 Evidence class 应在查看交集前冻结

当前设计是在“查看交集结果后、逐个判定前”冻结 evidence class。这仍然属于 outcome-informed rule setting。

应在以下时间点冻结：

- 查看候选身份之前；
- 最好连交集规模都不查看；
- 或者使用 discovery subset 制定规则，再应用到 untouched validation subset。

---

## 14. 文档内部需要立即修正的事实和状态冲突

### 14.1 “GSE286075 无 FASTQ”不准确

GEO supplementary 仅提供 gene-count 表，但官方记录显示原始 reads 在 SRA 中可用。

建议统一改为：

> 当前本地包和 GEO supplementary 只包含 gene counts；原始 reads 可从 SRA 获取。本阶段不下载、不用于 APA。

不能继续写成“无 FASTQ”。

### 14.2 P1 完成记录仍保留错误旧版本

P1 完成记录仍写“低深度、4 concordant、弱证据”，与修正链特异性取列后得到的 409 个 DE 和 17/31 方向一致冲突。

应更新主完成记录，将错误版本仅保留在 incident log 和 `.wrongstrand.bak` 审计材料中。

### 14.3 QKI CLIP accession 在附录仍写错

正文已改为 GSE147119，但附录仍存在 GSE146935 CLIP 的旧记录。应统一为：

- GSE147119：QKI CLIP；
- GSE146935：CRISPR-TRAPseq。

### 14.4 GSE286075 原论文 205 DE 与重分析 409 DE 需要解释

409 不一定错误，但与原研究报告的 205 相差较大。建议增加 reproduction appendix，列出：

- 原论文设计公式；
- 是否使用 paired design；
- strand-specific count column；
- filtering；
- DESeq2/PyDESeq2 版本；
- independent filtering；
- FDR 阈值；
- 原 205 个基因中有多少被重现。

在差异来源解释清楚前，不宜只把 409 作为“修正后的真实数字”。

### 14.5 P2-05e 已完成但未勾选

该问题不影响科学结论，但会降低交接和审计可信度。应将状态与完成记录同步。

---

## 15. 推荐的修订后研究结构

### 模块 A：Stroke-associated astrocyte APA discovery

主要数据：GSE238125。

主要问题：

> 卒中后 cortical astrocyte 是否存在可重复、时间依赖的 3′UTR 使用变化？

主要输出：

- 全局 PDUI 时间轨迹；
- event-level discovery；
- 跨算法一致性；
- coverage 和 PAS 验证；
- 技术偏差敏感性分析。

### 模块 B：Concurrent regulator-state annotation

主要数据：GSE174574。

主要问题：

> 卒中后 astrocyte 中哪些 APA/RNA-processing regulators 的表达发生变化？

输出应称为：

> concurrent upstream candidates

而不是：

> causal drivers

### 模块 C：Subcellular/spatial biological-context annotation

- GSE74456：PAP localization。
- GSE286075：endfoot stroke response。
- GSE263986：stroke spatial-zone response。

三者分别作为正交生物学背景，不进行“2/3 等票合并”。

### 模块 D：Cis-regulatory annotation

- 文献 perturbation：调控先验；
- CLIP：实际结合证据；
- motif：预测证据；
- conservation：进化约束；
- MPRA：仅在 tile sequence 可精确映射到 Δ3′UTR 时作为实验序列证据。

### 模块 E：Independent validation

第二公共数据集只验证预先固定的：

- 事件坐标；
- 方向；
- effect-size threshold；
- 候选顺序。

不得在验证集中重新调阈值、重选 breakpoint 或更改 candidate rule。

### 模块 F：Wet-lab causal closure

通过 regulator perturbation、isoform validation 和 localization reporter，把关联性结果升级为机制链。

---

## 16. 下一步执行优先级

### 16.1 可以立即继续

1. 跑完全部 12 样本 DaPars2。
2. 生成完整 PDUI matrix 和样本级 QC。
3. 检查 d1–d60 的全时间轨迹，不只关注 day 1。
4. 明确 P2 的差异模型、missingness 和 FDR 范围。
5. 对候选事件应用最低复核条件：
   - 四个对比样本均可量化；
   - 两个重复方向一致；
   - 达到预冻结 `|ΔPDUI|`；
   - 第二 APA 工具方向一致；
   - coverage 人工复核。

### 16.2 重写设计后再运行

1. M4 三参考集分析：先拆分 PAP、endfoot response 和 stroke zone。
2. Regulator→target 汇总：分开 perturbation、CLIP 和 motif。
3. 四重交集：改为事件级证据矩阵，硬交集仅作为描述性结果。

### 16.3 暂停

1. GSE330741 gene-level MPRA overlap，直至 library manifest 对应关系解决。
2. “假说链内部验证通过”的结论。
3. 用 QKI CLIP 证明 QKI 调控 APA。
4. 用 GSE286075 作为核内 APA regulator 的第二验证层。

---

## 17. 推荐湿实验闭环

若最终候选围绕 PABPN1，建议使用以下最短机制闭环。

### 17.1 候选 APA 验证

- astrocyte 3′RACE；
- proximal/distal isoform-specific qPCR；
- stroke 或 OGD/R 条件下验证 long/short isoform ratio；
- 必要时使用 targeted long-read sequencing 明确 terminal isoform。

### 17.2 上游因果

- Pabpn1 knockdown 是否重现 candidate shortening；
- Pabpn1 rescue 是否逆转 shortening；
- 测量 PABPN1 蛋白和核定位，而不仅是 mRNA；
- 同时检测候选基因总表达，区分 APA 与表达量改变。

### 17.3 下游定位

- long-UTR reporter；
- short-UTR reporter；
- smFISH/RNAscope；
- 用血管 lectin/CD31 和 astrocyte AQP4/GFAP 定义 perivascular endfoot；
- 定量 endfoot/soma RNA ratio，而不是只看总细胞荧光。

### 17.4 顺式元件因果

至少比较：

1. full long UTR；
2. short UTR；
3. long UTR 删除候选 Δ3′UTR element；
4. short UTR 补回候选 element。

比较各构建体的：

- RNA localization；
- local translation；
- total RNA abundance；
- RNA stability；
- reporter protein abundance。

这样才能真正闭合：

> PABPN1 → candidate APA → Δ3′UTR element → endfoot localization

如果资源有限，至少应完成以下三项中的两项：

- 3′RACE/isoform-specific qPCR；
- long-vs-short UTR reporter；
- vascular-endfoot smFISH/RNAscope。

只做总 RNA qPCR 可以验证表达或 APA，但不能证明亚细胞定位。

---

## 18. 推荐的结论语言分级

### Level 1：直接数据支持

- “Stroke was associated with altered PDUI in cortical astrocytes.”
- “The candidate event showed reproducible shortening in both biological replicates.”
- “The transcript was altered in the post-ischemic endfoot translatome.”

### Level 2：正交证据支持

- “The lost distal segment overlapped a QKI CLIP peak.”
- “The transcript was enriched in a PAP-TRAP reference.”
- “The observed direction was consistent with prior PABPN1 perturbation studies.”

### Level 3：预测性表述

- “The short isoform may have altered localization potential.”
- “Loss of the distal segment is predicted to remove candidate RBP-binding elements.”
- “PABPN1 is a candidate upstream regulator requiring perturbational validation.”

### 禁止在当前公共数据阶段使用

- “PABPN1 drives the stroke-induced APA program.”
- “QKI regulates this APA event.”
- “APA shortening impairs endfoot transport.”
- “The candidate loses endfoot localization after stroke.”
- “The causal pathway has been internally validated.”

---

## 19. 最终审计结论

### 19.1 Go/No-Go

**不是 No-Go。**

建议继续 GSE238125 全量 APA，因为它仍可能形成文章最重要的 discovery layer。当前需要修改的是论证结构，而不是放弃整个课题。

### 19.2 推荐决策

- P2 全量 APA：继续。
- P1 regulator：降级为 concurrent upstream candidate layer。
- P3：改为正交注释和候选优先化，不再声称机制调控已成立。
- endfoot/PAP/zone：拆开解释。
- QKI：定位为 cis-binding/transport-related annotation，不作为已建立 APA driver。
- MPRA：解决 manifest 与官方研究描述的矛盾前暂停使用。
- 最终论文：聚焦 `stroke-associated astrocyte APA remodeling and altered subcellular localization potential`。
- 真正的机制箭头：交给 PABPN1/QKI perturbation、isoform validation 和 localization wet experiments。

### 19.3 一句话评价

> **执行日志、失败记录和可追溯性做得很好；当前最需要的不是继续增加数据层，而是删除证据无法支持的因果箭头，让每一种数据只回答它真正能够测量的问题。**

---

## 20. 关键公共来源

- [GSE174574：Single-cell RNA-seq reveals the transcriptional landscape in ischemic stroke](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174574)
- [GSE238125 示例样本：Astrocyte day 1 rep 1](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7658907)
- [GSE286075：Astrocyte endfoot RiboTag after cerebral ischemia](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE286075)
- [GSE263986：Spatially defined astrocyte response zones after stroke](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263986)
- [GSE330741：In vivo MPRA of astrocyte mRNA localization](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE330741)
- [PABPN1 suppresses alternative cleavage and polyadenylation sites](https://pubmed.ncbi.nlm.nih.gov/22502866/)
- [QKI-6 CLIP and astrocyte maturation study](https://pmc.ncbi.nlm.nih.gov/articles/PMC7943582/)
- [Astrocytes locally translate transcripts in peripheral processes](https://pmc.ncbi.nlm.nih.gov/articles/PMC5441704/)
- [In vivo SN-MPRA study](https://pmc.ncbi.nlm.nih.gov/articles/PMC13142395/)
- [DaPars2 official repository](https://github.com/3UTR/DaPars2)

