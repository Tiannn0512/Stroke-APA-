# 新颖性/相似研究检索日志 — 07_followup_20260925

- 检索日期：2026-09-25
- 检索人：ZCode agent（新颖性复audit follow-up）
- 项目链条：卒中(MCAO, 小鼠) → 星形胶质细胞 Atp2a2(SERCA2) 3'UTR 长短异构体比例改变(候选 APA, 预测 PAS mm10 chr5:122456498) → 异构体在胞体 vs PAP 的分布 → 卒中改变该分布 → QKI 等调控因子参与
- 交付：`../novelty_overlap.tsv`（31 条外部文献线索，evidence_class 全部为"外部文献线索"）
- 纪律声明：TSV 中 `what_it_demonstrated` 仅写摘要/检索记录中实际读到的内容；未读到摘要的条目标注"仅标题/记录线索"并降为 low/medium；综述一律标注"综述背景"，不计为原始证明。

---

## 1. 检索过程与命中数

### 方向 1：Atp2a2/SERCA2 与卒中

| # | 日期 | 引擎/接口 | 检索式 | 命中 | 结果 |
|---|------|-----------|--------|------|------|
| 1 | 2026-09-25 | WebSearch | `SERCA2 Atp2a2 cerebral ischemia stroke astrocyte` | — | 命中 Salman(SERCA2激活抗卒中)、RIPostC-SERCA2、多篇星形胶质细胞钙信号综述；无星形胶质细胞特异 SERCA2 原始研究 |
| 2 | 2026-09-25 | PubMed eutils esearch | `(Atp2a2[MeSH Major Topic] OR Atp2a2)` | 1145 | 取 relevance 前 50 逐条 esummary 筛读：无脑卒中/星形胶质细胞主文章（最近者为 Shiga 2024 视网膜神经元 SERCA2 基因增强） |
| 3 | 2026-09-25 | PubMed eutils esearch | `(SERCA2[A/B] OR Atp2a2) AND (ischemi* OR stroke OR MCAO)` | 94 | 前 40 逐条 esummary：卒中相关仅 5 条（42161231, 37897375, 41208593 肌肉, 37132383 ERO1α, 29961273 电针）；全脑组织/药理层面，无细胞类型 |
| 4 | 2026-09-25 | Europe PMC REST | `ABSTRACT:Atp2a2 AND (ABSTRACT:MCAO OR "focal ischemia" OR photothrombosis OR "cerebral ischemia")` | 1 | 仅电针-海马-钙信号基因（PMID 29961273, 2017），bulk 组织 |
| 5 | 2026-09-25 | PubMed eutils esearch | `(SERCA2[A/B] OR Atp2a2) AND astrocyte` (Title/Abstract) | 1 | 唯一命中为晶状体 MANF/白内障（PMID 42419372），与星形胶质细胞-SERCA2 无关 |
| 6 | 2026-09-25 | Europe PMC REST | `ABSTRACT:Atp2a2 AND ABSTRACT:astrocyte` | **0** | **无任何摘要同时提及 Atp2a2 与 astrocyte 的论文** |

### 方向 2：Atp2a2 与 APA/3'UTR 异构体

| # | 日期 | 引擎/接口 | 检索式 | 命中 | 结果 |
|---|------|-----------|--------|------|------|
| 7 | 2026-09-25 | WebSearch | `ATP2A2 SERCA2 3'UTR alternative polyadenylation isoforms` | — | 命中 UniProt P16615（注释：SERCA2 各转录本仅 3'-UTR 不同、组织特异性表达）及经典肌肉文献 |
| 8 | 2026-09-25 | PubMed eutils esearch | `(ATP2A2[A/B] OR SERCA2[A/B]) AND (polyadenylation OR "3' UTR" OR "3' untranslated" OR "untranslated region")` | 16 | 全部 16 条 efetch 摘要：APA/3'末端加工证据集中在猪/丰年虾/哺乳动物肌肉（1991-1998）；3'UTR 稳定性/翻译为心肌-平滑肌与心脏（2002-2025）；**无一涉及脑、卒中或星形胶质细胞** |

### 方向 3：PAP 定位 / 卒中改变定位 / 卒中+APA+星形胶质细胞

| # | 日期 | 引擎/接口 | 检索式 | 命中 | 结果 |
|---|------|-----------|--------|------|------|
| 9 | 2026-09-25 | WebSearch | `astrocyte perisynaptic processes mRNA localization transcriptome Sakers` | — | 定位到 Sakers 2017 PNAS、Mazaré 2020 Cell Rep、Boulay 2017 Cell Discovery |
| 10 | 2026-09-25 | WebSearch | `Sakers "perisynaptic astrocytic processes" 2017 sensory experience PMID` | — | 搜索引擎自称存在"Sakers 2017 J Neurosci sensory deprivation"论文 PMID 28794256；经 Europe PMC 验证**不存在**（见 #11/#12），已弃用该线索 |
| 11 | 2026-09-25 | Europe PMC REST | `"synapse-associated astrocyte transcriptome"` | **0** | 证实 #10 中该标题为搜索引擎虚构/混淆 |
| 12 | 2026-09-25 | Europe PMC REST | `TITLE:"sensory deprivation" AND ABSTRACT:astrocyte` | **0** | 同上 |
| 13 | 2026-09-25 | Europe PMC REST | `AUTHOR:"Sakers K"` | 29 | 完整作者列表核对：PAP 线为 PNAS 2017 / Cell Rep 2020(合著) / Cell Rep 2022 / Nat Commun 2021(QKI)；另有 Koester 2026 bioRxiv SN-MPRA |
| 14 | 2026-09-25 | Europe PMC REST | `TITLE:"perisynaptic astrocytic processes"` | 7(5 正式发表) | 全部筛读：Mazaré 2020、Denizot 2026(PAP 内 ER)、Hussain 2023、Bellesi 2018、Puletti 2026(综述) |
| 15 | 2026-09-25 | WebSearch | `alternative polyadenylation cerebral ischemia stroke 3'UTR isoform brain` | — | 命中"rNLS8 APA"线索，核对后为 ALS/FTLD-TDP（Eck 2025, Mol Psychiatry），**非卒中** |
| 16 | 2026-09-25 | Europe PMC REST | `TITLE:"alternative polyadenylation" AND (ABSTRACT:stroke OR "cerebral ischemia" OR MCAO)` | **0** | **没有任何论文标题含 APA 且摘要涉及卒中/脑缺血** |
| 17 | 2026-09-25 | Europe PMC REST | `"alternative polyadenylation" AND (ABSTRACT:"middle cerebral artery occlusion" OR ABSTRACT:"focal ischemia")` | 1(松散) | 唯一命中为 FTO/m6A 去甲基化-脑缺血（Xu 2023），APA 仅以背景词出现，非 APA 研究 |
| 18 | 2026-09-25 | Europe PMC REST | `("3'UTR" OR "3' untranslated" OR "alternative polyadenylation") AND astrocyte AND (stroke OR "cerebral ischemia" OR MCAO)` | 37 | 37 条全部为 miRNA/exosome 靶向 3'UTR（AQP4 等）或星形胶质细胞-卒中综述；**无 APA 异构体研究**；唯一亮点 Shim 2025 endfeet RiboTag(PMID 39796165) |
| 19 | 2026-09-25 | WebSearch | `quaking QKI RNA binding protein astrocyte stroke cerebral ischemia` | — | 命中 Sakers 2021 Nat Commun(QKI-星形胶质细胞)、QKI-APA(EMT)、QKI6-CIRI(神经元)；检索引擎确认 QKI-卒中交叉是研究空白 |
| 20 | 2026-09-25 | Europe PMC REST | `(ABSTRACT:quaking OR ABSTRACT:QKI) AND (ABSTRACT:"alternative polyadenylation" OR "cerebral ischemia" OR "ischemic stroke")` | 7 | Neumann 2024 RNA Biol(QKI 调控 EMT APA)、Liu 2021 Brain Behav(QKI6-CIRI 神经元)、circFOXP1 等 |
| 21 | 2026-09-25 | Europe PMC REST | `TITLE:"dendritic localization" AND ABSTRACT:BDNF` + eutils 补充 | — | BDNF 长 3'UTR-树突定位先例线：Liao GY 2012 Nat Med (PMID 22426422) |
| 22 | 2026-09-25 | WebSearch | `astrocyte distal processes transcriptome Boulay 2017` | — | Boulay 2017 Cell Discovery, PMID 28377822（endfeetome） |

数据库检索均通过 eutils(eutils.ncbi.nlm.nih.gov) 与 Europe PMC REST(ebi.ac.uk/europepmc/webservices/rest)；PubMed 网页因 reCAPTCHA 不可读，全部改用 API。

---

## 2. 关键相邻工作（最接近我们链条的文献）

1. **Shim B et al., Int J Mol Sci 2025 (PMID 39796165)** — 小鼠 MCAO/再灌注 + 星形胶质细胞 RiboTag 于分离脑微血管：缺血引起 endfeet 翻译组 205 个变化、HSP70 在 endfeet 上调。**这是"卒中改变星形胶质细胞亚细胞翻译/定位"的最近邻证据**，但：定位在血管端 endfeet（非突触旁 PAP）、测的是全长转录本丰度/翻译（非 3'UTR 异构体/APA）、未涉及 Atp2a2。
2. **Mazaré 2020 / Sapkota 2022（Cell Rep）** — 生理活动（恐惧条件）与癫痫样活动改变 PAP 内核糖体结合 mRNA/蛋白组：证明"PAP 翻译组可被活动改写"，但非卒中、非异构体层面。
3. **Sakers 2017 PNAS (PMID 28439016)** — 星形胶质细胞外周突起局部翻译且转录本富集 **quaking 结合基序**：把 QKI 与星形胶质细胞突起 mRNA 定位在机制上连接起来，但无卒中、无 APA。
4. **Sakers 2021 Nat Commun (PMID 33750804)** — QKI-6 结合星形胶质细胞 mRNA 3'UTR 并稳定靶标。
5. **Neumann 2024 RNA Biol (PMID 38112323)** — QKI 是 APA 景观的调控因子（EMT/癌细胞）。
6. **Eggermont 1991 / Misquitta 2002, 2005 / Nie 2004 / Umar 2025** — Atp2a2 的 3'末端多聚腺苷酸化位点、3'UTR 决定 mRNA 稳定性、poly(A) 尾长度与 RBFOX1-3'UTR 翻译调控：**Atp2a2 的 3'UTR 异构体在肌肉/心脏确实存在并有功能**，但全部在非脑组织。
7. **Salman 2026 (PMID 42161231) / Chen 2024 (PMID 37897375)** — 卒中-SERCA2 直接证据，但为药理激活/组织水平表达，无细胞类型、无异构体。
8. **Denizot 2026 Glia (PMID 41098061)** — 75% 的 PAP 含内质网：SERCA2 蛋白在 PAP 发挥作用的结构前提，但没有 RNA 证据。

## 3. 明确的阴性结果（对新颖性判断最重要）

- `TITLE:"alternative polyadenylation" AND (stroke OR cerebral ischemia OR MCAO)`：**0 篇**。
- `"alternative polyadenylation" AND (MCAO OR focal ischemia)`：仅 1 篇松散命中（FTO/m6A，非 APA 研究）。
- `ABSTRACT:Atp2a2 AND ABSTRACT:astrocyte`：**0 篇**。
- `(SERCA2|Atp2a2) AND astrocyte`（Title/Abstract）：1 篇，核对为无关（晶状体/白内障）。
- `Atp2a2 AND 脑缺血`：仅 1 篇 2017 电针-海马 bulk 基因表。
- "Sakers 2017 J Neurosci sensory deprivation PAP" 检索引擎说法经验证不存在（0 命中）；真实 Sakers 2017 为 PNAS 局部翻译论文。
- 卒中方向已检索到 APA 的研究均不在卒中：Eck 2025 Mol Psychiatry（rNLS8 ALS/FTLD-TDP 的 APA 紊乱）说明"神经疾病+APA"已有先例，但疾病模型是 TDP-43 蛋白病而非缺血。

## 4. 结论：新颖性判定

**问题：是否已有研究覆盖"卒中改变星形胶质细胞 Atp2a2 3'UTR 异构体组成及其亚细胞定位"？**

**没有。** 截至 2026-09-25，没有任何已发表或预印本研究同时覆盖该链条中的任意两环以上，更无全链条研究：

- "卒中 → Atp2a2 丰度"：存在（药理激活/组织水平，2024/2026；电针 bulk 表 2017），但**无细胞类型分辨、无异构体分辨**。
- "Atp2a2 APA 存在"：存在（猪/肌肉 1991-1998、心脏 3'UTR 功能 2002-2025、UniProt 注释人 ATP2A2 仅 3'UTR 不同的组织特异转录本），但**从未在脑/星形胶质细胞中报道**。
- "星形胶质细胞 mRNA 异构体定位"：PAP 局部翻译与 mRNA 定位领域成熟（2017-2026），且 2026 SN-MPRA 预印本刚证明星形胶质细胞 mRNA 定位由顺式元件决定，但**无任何 3'UTR 异构体（APA）层面的 PAP 定位研究，未涉及 Atp2a2**。
- "卒中改变定位"：唯一近邻是 Shim 2025（缺血改变 endfeet 翻译组），但 endfeet≠PAP、丰度≠异构体。
- "QKI 参与"：三个环节各有独立文献（QKI=APA 因子；QKI=星形胶质细胞 mRNA 3'UTR 结合蛋白且 PAP 转录本富集 QKI 基序；QKI6 在缺血再灌注中于神经元起作用），**但 QKI→星形胶质细胞 Atp2a2 APA→PAP 定位这一通路无人连接**。

**我们链条中仍新颖（无任何外部覆盖）的环节**：
1. 卒中改变星形胶质细胞 **Atp2a2 3'UTR（APA）异构体组成**——完全空白；
2. Atp2a2 长短 3'UTR 异构体在**胞体 vs PAP 的差异分布**——完全空白；
3. 卒中**改变**该异构体空间分布——完全空白（Shim 2025 提示卒中可改变星形胶质细胞亚细胞翻译组，是我们的概念最近邻，可作为可行性引文而非竞争）；
4. QKI 对 Atp2a2 APA 的调控——完全空白（QKI-APA 与 QKI-星形胶质细胞 3'UTR 各有先例）。

**风险与保留**：
- Koester/Sakers 2026 bioRxiv SN-MPRA 表明该组正系统推进"星形胶质细胞 mRNA 定位顺式元件"，其平台随时可能延伸到 APA/疾病模型，建议持续跟踪该组预印本；
- "卒中+APA+脑"若发生，最可能先出现在全脑 3'Seq 资源性论文而非机制论文（截至检索日未见）；
- 检索盲区：中文文献、非英文数据库未系统覆盖；PubMed 1145 条 Atp2a2 记录仅逐条筛读 relevance 前 50，Atp2a2+缺血的重叠文献理论上可能存在于长尾（但方向性 esearch 已另行覆盖）。
