# P3-01 Level 1 文献池（草案 DRAFT — 2026-09-17，主会话检索）

> **2026-09-17 11:0x 更新：正式表已建 = `results/p3_level1_pairs.tsv`（12 对，PMID 经 NCBI eutils 逐条核验；QKI PMCID 行与 ELAVL1 经典 JBC 行标注 partial/n）。本草案保留作检索过程记录。**
> L1 定义（冻结）：已发表 perturbation（KD/KO/OE）实验证明该 regulator 改变靶基因 APA/3'UTR 亚型。唯一允许写"调控"的级别。
> 方向约定：我们的 M1 数据 = Pabpn1↓（padj 3e-4）、Qk↓（L1 0.077/L2 0.122）、Nudt21 微降、Elavl1↑、Ptbp1 微升（均两层一致）。

## 检索到的 regulator–target 对

| regulator | target | perturbation | APA 方向 | 物种/体系 | PMID | 年 | 证据摘要 | 置信 |
|-----------|--------|--------------|----------|-----------|------|----|----------|------|
| PABPN1 | 全基因组（数百基因；OPMD 肌肉） | 功能缺失/表达降低 | **3'UTR 缩短**（近端 PAS 使用↑） | 人/OPMD 肌肉+细胞 | 22772983 | 2012 | PABPN1 减少致替代性 polyA 位点选择改变、近端位点使用增加 | 高 |
| PABPN1 | 全局（遗传筛选） | KD | **全局 3'UTR 缩短** | 人细胞（HeLa 等） | Jenal 2012 Cell（PMID 待核验，标题 "The Poly(A)-Binding Protein Nuclear 1 Suppresses Alternative Polyadenylation"） | 2012 | 遗传筛选鉴定 PABPN1 为 APA 抑制子；KD 后近端 PAS 启用 | 高（PMID 待核） |
| NUDT21/CFIm25 | 全局；CCND1 等数百转录本 | KD | **广泛 3'UTR 缩短** | 人 HeLa + 胶质瘤 | 23103948 | 2012 | CFIm25 KD → 近端 PAS 启用、3'UTR 缩短；CCND1 短 3'UTR 逃逸 miR-342-3p；GBM 中缺失 | 高 |
| NUDT21/CFIm25 | 全局 | KD | 3'UTR 缩短 | 人 293T | Gruber 2012（PMID 待核验） | 2012 | RNA-seq 显示 CFIm25 KD 后全局 3'UTR 缩短 | 中高（PMID 待核） |
| NUDT21/CFIm25 | 神经分化程序 | 部分缺失（ hypomorph） | APA 改变 | 小鼠（在体，脑） | Alcott 2020 eLife（PMID 待核验） | 2020 | 部分 CFIm25 缺失 → 学习缺陷+神经分化异常——**脑在体上下文** | 中高（PMID 待核） |
| QKI | 自身（Qki-5 自调控）+ EMT 细胞 APA 子集 | 扰动（EMT 诱导/QKI-5） | APA 事件改变（含自身 3'UTR 长度） | 人细胞（EMT） | PMC10732628（RNA Biology） | 2023 | QKI 调控 3'UTR 内 APA 子集，CLIP 支持自调控 | 中 |
| ELAVL1/HuR | U-rich PAS 位点（位点封闭机制） | 体外重组/封闭实验 | 抑制 U-rich 位点的切割与多聚腺苷化（近端偏移机制） | 体外+果蝇 ewg 经典 | JBC 282:2203-2210（PMID 待核验）；Frontiers Genet 2022 综述 848626 | 2007/2022 | Hu 蛋白选择性阻断含 U-rich 序列的 PAS；nELAVL 神经特异 | 中（机制确凿，全基因组 KD-APA 表未直接检索到） |
| PTBP1 | 多基因（含 3'UTR miRNA 竞争） | KD | APA 事件改变（PTB 缺失诱导） | 人成纤维细胞→神经元 | 23260894 | 2013 | PTB KD 诱导 APA 改变；PTB 与 miRNA 在 3'UTR 位点直接竞争 | 中高 |
| PTBP1 | Pbx1 等 | KD/iCLIP | 剪接+表达程序（含 UTR 调控） | 小鼠神经元 | 23260894 关联文献 PMC4755740 | 2015 | PTBP1 控制神经元基因表达程序 | 中（非 APA 直接证据，作辅助） |
| NOVA1/2 | 全脑转录组（3'UTR CLIP+测序） | Nova KO（脑在体） | **Nova 调控脑内 APA**（3'UTR 结合全基因组证据） | 小鼠脑 | 18978773 | 2008 | HITS-CLIP 揭示 Nova 大量结合 3'UTR 并调控脑 APA（Nature） | 高 |
| NOVA2 | 跨神经元细胞类型差异 APA | NOVA2（小脑发育） | APA 差异加工 | 小鼠小脑在体 | Jereb 2018 eLife 7:e34042（PMID 待核验） | 2018 | NOVA2 驱动跨神经元类型的 3' 端差异加工 | 中高（PMID 待核） |
| MBNL1/2 | 发育性 APA 程序 | MBNL 隔离（CUGexp）+敲低 | **发育性 3'UTR 缩短/APA 选择破坏** | 人+小鼠（DM1 模型） | 25454948 | 2014 | MBNL 失活破坏发育调控的 APA（Mol Cell）；MBNL2 为脑型 | 高 |
| TARDBP/TDP-43 | 3'UTR polyadenylation 全局 | TDP-43 核缺失（rNLS8 在体） | 3'UTR polyadenylation 紊乱 | 小鼠（ALS/FTLD 模型，脑） | Mol Brain 2025, 10.1186/s13041-025-01174-1（PMID 待核验） | 2025 | TDP-43 失位致 3'UTR polyA 紊乱并与表达改变相关 | 中（模型为失位而非 KD，PMID 待核） |

## 上下文（非 L1，不计对数）

- QKI 结合（非 APA perturbation）：Hnrnpa1 3'UTR QRE（PLOS Genet 2010, pgen.1001269）；Sirt2 3'UTR（PMC5392665）；**星形胶质细胞 mRNA 3'UTR CLIP（= 我们 GSE146935 的原始文献，PMID 33750804，属 L2）**；qkv 小鼠少突转译组（Sci Rep 2017）。
- Brumbaugh 2018 Cell（NUDT21–染色质–APA）；Masamha 2014 综述/方法（PMC3912697）；DaPars 方法学（PMC4467577）。
- 3'aTWAS 2025（含缺血性卒中）——P4 讨论用。

## 覆盖度（17 个两层一致 regulator 的 L1 缺口，负结果照登）

已覆盖（有 L1 级证据）：**Pabpn1、Nudt21、Qk、Elavl1、Ptbp1、Nova1/2、Mbnl1/Mbnl2、Tardbp**（10/17 个 regulator 覆盖）。
未检索到直接 L1 证据（第二批检索 FUS 已做、无果；如实标"无 L1"）：**Fus**（FUS 文献集中于剪接/保留，未见 KD→APA 直接研究，2026-09-17 检索）；Wdr33、Cpsf1、Cpsf3、Cpsf4、Fip1l1、Cstf2、Clp1、Hnrnpa2b1（8 个待第三批或标无）。
  - 备注：Hnrnpa2b1 有 UCSD 学位论文级证据（depletion 后结合位点附近 3'UTR 选择改变，待正式发表源核验）；核心机器组分（Wdr33/Cpsf1/3/4/Fip1l1/Cstf2/Clp1）多为必需基因，KD 常致全局加工缺陷而非干净靶点级 APA 证据——这本身是可报告的负结果。

## 对 P4 漏斗的初步含义（不改变设计，仅备忘）

- Pabpn1↓ + 其 L1"KD→3'UTR 缩短"组合预测：卒中 astrocyte 中 PABPN1 靶基因整体偏缩短——可在 M2 结果中做方向性 sanity check（不许事后调阈值，仅描述）。
- Qk↓ + QKI 星形胶质细胞 3'UTR CLIP（GSE146935）= L2 主力组合；M3-02 intersect 时靶基因池应包含 QRE 载量高的 endfoot 基因。
