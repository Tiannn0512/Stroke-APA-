# 卒中→APA→星形胶质细胞 endfoot 定位潜能：Pipeline v5.1 与 Todolist

> 生成日期：2026-09-16 ｜ 版本：v5.1（v5.0 经导师审阅修订 + T0.1/T0.2 实查证据更新）；**v5.1r（2026-09-19）：论证逻辑审计修订——当前有效主张以 §0b 为准**，审计原文 = `reviews/2026-09-18_论证逻辑审计报告_20260918.md`
> 核心假说不变：Stroke → APA regulator 改变 → APA 改变 → 特定 astrocyte mRNA 3′UTR 亚型改变 → endfoot-localized mRNA 中的 APA 事件 → 3′UTR 序列（element）删除/获得 → localization-related cis-element / RBP-binding element 改变 → 预测 endfoot 定位潜能改变
> 定版原则（引用审阅意见）：**所有最终定位结论限定为 "localization potential"，不声称已证明实际运输。**

---

## 0. v5.0 → v5.1 变更记录

| # | 审阅意见 | v5.1 落实 |
|---|----------|-----------|
| R1 | M1 与 M2 绑定太死：不要让 GSE174574 单点承担 regulator DE + APA 两层证据；两个模块允许用不同数据集 | M1/M2 正式解耦：M1 = regulator 检测层（多数据集任选），M2 = APA 检测层（主源 + 已核实备选池）。T0.2 实查确认 GSE174574 为 10x 3′ 后，M2 工具链改写并预置 3 套已核实 Plan B |
| R2 | M3 有过度解释倾向：motif 存在 ≠ regulator 真的结合；Δexpr+ΔPDUI ≠ 因果 | M3 改为三级证据分级（Level 1 perturbation 文献证据 > Level 2 CLIP+方向一致 > Level 3 motif+相关），输出语言纪律随分级锁定 |
| R3 | M5 Δscore 模型锁死过早：跨 3′UTR 的 MPRA 活性不能线性相加；"有类似 motif" ≠ "经实验验证的 localization activity" | 删除 Δscore，改为 Localization evidence class（MPRA > CLIP > motif > conservation 四层 → A/B/C/D 级），每个 candidate 一张证据矩阵表 |
| R4 | 应保留从 endfoot 反推的平行支线，四重交集作为核心漏斗 | 反向支线升级为核心漏斗的一支；新增 GSE263986（分区反应性星形胶质细胞转译组）补强 endfoot 侧参考 |
| R5 | 停止继续设计，直接执行 T0.1/T0.2 四件事：① 文库类型 ② 每样本 astrocyte 细胞数 ③ 3′ 端信息可用性 ④ 核心 regulator 表达快查 | T0.2 已完成（结论：10x 3′）；③④ 所需的 count 矩阵已确认存在，tar 下载中；T0.5（regulator 快查，无需 FASTQ）加入清单 |

---

## 0b. v5.1 → v5.1r 修订记录（2026-09-19 论证逻辑审计，全文采纳；判定 = Major revision / 有条件 Go）

**修订后中心主张（替代原因果链表述）**：
> 缺血性卒中伴随皮层星形胶质细胞**时间依赖的 3′UTR 重塑**；其中一部分 APA 事件发生在与外周突起定位或血管周终足响应相关的转录本中，并改变候选顺式调控元件，从而产生"APA 影响亚细胞 RNA 定位潜能"的**可验证机制候选**。
> EN: Ischemic stroke is associated with time-dependent 3′UTR remodeling in cortical astrocytes. A subset of these APA events occurs in transcripts implicated in peripheral-process localization or the perivascular endfoot response and alters candidate cis-regulatory elements, thereby prioritizing mechanisms through which APA may reshape subcellular RNA localization.

**证据结构定性**：各数据集来自不同动物队列、卒中模型、时间点与测量层级，构成的是**三角互证（方向一致性 + 机制可解释性）**，不构成闭合因果链；regulator→APA 与 APA→定位两个因果箭头交由湿实验闭环（模块 F）检验。

| # | 审计发现 | v5.1r 落实 |
|---|----------|-----------|
| A1 | 单一因果链超出证据 | 改为模块化结构：A 发现层（GSE238125）/ B regulator 并发注释（GSE174574）/ C 空间背景三轴 / D 顺式元件注释 / E 预固定独立验证 / F 湿实验闭环 |
| A2 | "Pabpn1↓+d1 缩短 = 内部验证通过"过强 | 降级为方向一致性表述（TODO P2-05g 注记）；pilot = 描述性 |
| A3 | P2 统计闭环缺失（n=2 检验、FDR 范围、共享 sham 相关性未定义） | 新增 P2-06a SAP 冻结：两层分析（样本级全局 + 事件级 moderated model）；n=2 明示为 discovery 层 |
| A4 | 缺失值规则未分层；44.8% 产出率归因过单一 | 量化层（≥30% 规则）≠ contrast 层（4/4 有效 PDUI，不插补）；产出率归因软化（含注释选择因素） |
| A5 | slim 注释边界 | 事件表述 = "DaPars2-inferred proximal cleavage shift within the selected longest annotated terminal 3′UTR reference"；候选最低复核（PAS 距离/第二工具/coverage 手检） |
| A6 | d1 缩短需排技术混杂 | 纳入 SAP 层 1：表达/coverage/UTR 长度分层 + gene-body coverage + 3′ bias + 深度指标 + LOSO |
| A7 | GSE286075 层错位（endfoot 转译组 ≠ 核内 regulator 验证） | 重定位 = endfoot_stroke_responsive 参考；17/31 = 描述性 concordance；Qk 表述改为"方向一致但均未达预设显著" |
| A8 | 三参考集异质却被合并投票 | M4 拆三独立注释轴（PAP_localized / endfoot_stroke_responsive / stroke_spatial_zone_response），废除"2/3 一致 = 高置信" |
| A9 | QKI 证据越界 | 只作 cis-binding/transport-related annotation；不作已确立 APA driver（原论文 PAS 阴性结果引注） |
| A10 | MPRA 数据对应矛盾（651 基因 vs 官方 Glt1/Sparc focused tiling ~6,430） | P4-04 = BLOCKED_PENDING_MANIFEST_RECONCILIATION；解决前仅背景注释，不入 evidence class |
| A11 | FIMO 阈值调整风险 | 1e-3+BH 定性为阳性对照校准：仅用独立 CLIP 阳性对照、候选扫描前冻结、motif family 归并、length/GC-matched null |
| A12 | 四重硬交集选择偏差 | 改事件级证据矩阵（APA 主证据定资格，机制注释定优先级）；硬交集降为描述性结果 |
| A13 | Gate 设计缺陷 | G2 → 组合门槛（全量统计前显式修订）；G3 → 三轴分别表述；G4 → 0–5 候选（允许零）；P4-03 冻结提前到看交集之前 |
| A14 | 事实错位 | "GSE286075 无 FASTQ"→"本地仅 gene counts（SRA 有原始 reads，本阶段不用）"；GSE147119/GSE146935 区分；205 vs 409 复现附录（P1-08）；P2-05e 勾选补齐 |

---

## 1. T0 实查证据日志（2026-09-16）

| 问题 | 结论 | 证据 |
|------|------|------|
| GSE174574 建库类型？ | **10x Genomics 3′ droplet**（非 Smart-seq2 全长） | 原文 Zheng 2022 JCBFM（PMID 34496660 / PMC8721774）正文："performed single-cell transcriptomic analyses with the **10x genomics**"；SRA 元数据 6/6 run 均为 PAIRED + LIBRARY_SELECTION=cDNA + Illumina HiSeq 4000（10x 特征组合），无 smart-seq/SMARTer 关键词 |
| 模型与取样？ | MCAO 24h，取整侧缺血半球 vs sham，C57BL/6，6–8 周，n=3+3 | 原文 + SRA sample attributes（source_name=Brain, treatment=sham/ischemic） |
| GEO 有 processed 矩阵吗？ | **有**：6 样本 × 10x 三件套（barcodes/genes/matrix.mtx），tar 共 251MB | ftp filelist 实查；barcodes 32–46KB gzip/样本 |
| 3′ 端信息可用性？ | 10x 3′ 文库 → 只能做 3′-tag 感知的 APA 工具（scAPAtrap / Sierra）；QAPA/DaPars2 常规用法不适用 | 文库类型推断（10x 3′ 只测序靠近 polyA 的一端） |
| astrocyte 细胞数？ | 待矩阵解包后精确计数（原文 17 clusters，含 astrocyte 簇；未给出每样本 astrocyte 数） | **T0.6 待办**：barcodes 行数 + marker 注释即可关闭 |
| Plan B 数据集？ | 已核实 4 套（见 §2 Plan B 池），全部为星形胶质细胞特异或脑缺血 bulk 型 RNA-seq，均有 SRA 原始数据 | eutils + ftp 逐一实查 |

**M1 额外红利**：GSE174574 的 processed 10x 矩阵可直接做 regulator DE（Seurat/Scanpy 读 mtx，无需 FASTQ 比对）——审阅意见第 4 项（5 个核心 regulator 快查）不需要等比对流水线，T0.5 就能做。

---

## 2. 数据资产 v5.1

### 2.1 主线数据集（同 v5.0，均已在 GEO 核实）

| 登录号 | 内容 | 关键参数 | 角色 | 原始数据 |
|--------|------|----------|------|----------|
| GSE174574 | MCAO 24h 脑 scRNA-seq（Zheng 2022） | **10x 3′**，MCAO vs Sham n=3+3，整侧半球 | M1 regulator DE 主源；M2 APA 主源（条件成立时，scAPAtrap/Sierra） | ✅ SRA SRP320164 + processed mtx |
| GSE286075 | Astrocyte endfoot RiboTag（Shim 2025） | MCAO 2h+6h 再灌注；ipsi/contra 配对 n=3+3 | **endfoot_stroke_responsive 转译组参考**（v5.1r：不再作核内 regulator 活性验证层；配对 DE 保留为描述性方向参考） | 本地仅 gene counts（SRA 有原始 reads，本阶段不用） |
| GSE74456 | PAPTRAP（Sakers 2017） | PAP/Ctx × TRAP/input | **PAP_localized** 参考（外周突起/突触周——v5.1r 明示非血管周 endfoot） | ✅ SRP065508 |
| GSE147119 | QKI-6 CLIP（Sakers 2021；v5.1r 更正——GSE146935 仅 CRISPR-TRAPseq） | 437 peaks，富集于终止密码子附近/PAS 上游 | RBP element 参考（Level 2 = cis-binding annotation，非 APA 调控证明） | ✅ 已取（Europe PMC） |
| GSE330741 | SN-MPRA（2026） | 官方描述：Glt1/Slc1a2/Sparc/Hsbp1 focused tiling ~6,430 elements；与本地 32 counts（651 基因/1631 tiles）矛盾 | localization element 金标准 → **⛔ manifest 核对中，P4-04 阻断（v5.1r）** | ✅ PRJNA1465359（32 counts 已下，待核） |

### 2.2 M2 Plan B 池（本次新核实，按优先级排序）

| 优先级 | 登录号 | 内容 | 设计 | 对 M2 的适配性 | 原始数据 |
|--------|--------|------|------|----------------|----------|
| PB-1 | **GSE238125** | 卒中后星形胶质细胞动态表达谱 | bulk RNA-seq，皮层 astrocyte，Sham + d1/d3/d7/d21/d60，n=2/时间点 | ✅ 首选：纯 bulk 全长建库（PAIRED cDNA 实查）→ DaPars2/QAPA 直接可用；多时间点方向一致性弥补 n=2 | ✅ PRJNA997998（12 runs） |
| PB-2 | **GSE225110** | 卒中胶质转译组（RiboTag） | Astrocyte IP + input，超急性 4h + 急性 3d，sham/stroke，雌雄，n=132 runs | ✅ 次选：重复数最多；RiboTag IP 是 mRNA 层，APA 可做（读长覆盖 3′ 端） | ✅ PRJNA922165 |
| PB-3 | **GSE263986** | 分区反应性星形胶质细胞转译组 | mGFAP-RiboTag + LCM 分区（皮层/白质 zone1–5，uninjured vs stroke） | ✅ 兼用：既是 M2 备选，又补强 endfoot 侧反应性 astrocyte 参考 | ✅ PRJNA1100502（82 runs） |
| PB-4 | GSE165384 | 脑缺血 bulk 转录组 | 海马，缺血 2h vs sham，n=3+3 | ⚠️ 末位备选：非皮层、非 astrocyte 特异 | ✅ SRP302969 |

> Plan B 决策规则：GSE174574 astrocyte 细胞数不足或 scAPAtrap/Sierra 输出不可靠 → M2 主源切 PB-1（GSE238125），PB-2 作为第二 APA 源做跨数据集方向一致性；PB-3 同时供 M4 做 endfoot 侧第二参考。M1 永远不依赖 M2 的数据集选择。

---

## 3. Pipeline v5.1（核心漏斗结构）

```
                    STROKE
                       │
             ┌─────────┴─────────┐
             ↓                   ↓
   Astrocyte RNA-seq        Endfoot datasets
   (M1: GSE174574 scRNA     (M4: GSE286075 endfoot
    + GSE286075 paired)      RiboTag + GSE74456 PAPTRAP
             │               + GSE263986 zones)
             ↓                   │
       Regulator DE              │
             │                   │
             ↓                   │
   Stroke-associated            │
   APA regulators               │
             │                   │
             ↓                   │
   M2: APA events               │
   (primary: GSE174574 astrocyte│
    pseudo-bulk, scAPAtrap/     │
    Sierra; PB-1/2 fallback)    │
             │                   │
             └─────────┬─────────┘
                       ↓
   M6a: CORE FUNNEL（四重交集）
   Regulator APA targets
        ∩ Stroke APA genes
        ∩ Astrocyte endfoot genes
        ∩ Stroke-responsive genes
                       ↓
              Endfoot-APA candidates
                       ↓
                Δ3′UTR region
                       ↓
        ┌──────────────┼──────────────┐
        ↓              ↓              ↓
   M5: MPRA      M5: CLIP       M5: Motif
   (GSE330741)  (GSE146935)    (ATtRACT+FIMO)
        └──────────────┼──────────────┘
                       ↓
                 Conservation
                       ↓
        M5: Localization evidence class (A/B/C/D)
                       ↓
              TOP 3–5 candidates
                       ↓
      "altered localization potential"（仅预测）
```

### M0 数据与工具底座（进行中）
- [x] T0.2 文库类型判定 → 10x 3′，APA 工具锁定 scAPAtrap/Sierra（+ DaPars2 仅用于 bulk Plan B）
- [~] T0.1 下载：GSE174574_RAW.tar（251MB，下载中）、GSE286075 tar、其余按需
- 坐标基准 mm10；GSE330741 如为 mm39 则 liftOver
- 工具栈：Seurat/Scanpy（10x mtx 直读）；STAR（bulk Plan B 比对）；scAPAtrap / Sierra（10x APA）；DaPars2（bulk APA）；bedtools；FIMO+ATtRACT；R（DESeq2/edgeR/clusterProfiler）

### M1 Regulator 层（谁变了）——与 M2 解耦【v5.1r：定位 = concurrent upstream candidates（并发上游候选），非因果 driver】
- 主源：GSE174574 processed mtx → QC（Scrublet）→ 聚类 → astrocyte 簇（Aqp4/Slc1a3/Gfap/Aldh1l1）→ homeostatic vs reactive（C3/Serpina3n）→ pseudo-bulk DE（MCAO vs Sham）
- 第二参考：GSE286075 配对 DE——v5.1r 重定位：测的是 endfoot 转译组响应，其 regulator mRNA 变化只作**描述性方向 concordance**，不作为核内 APA regulator 活性的独立验证层
- Regulator 名单两级：核心机器（CPSF1/2/3/4、WDR33、FIP1L1、CSTF1/2/3、NUDT21/CFIm25、CPSF6/7、PABPN1、PCF11、CLP1、SYMPLEK）+ RBP 调制型（ELAVL1/2、QKI、PTBP1/2、HNRNPA2B1、SRSF1、NOVA1/2、FUS、TARDBP、KHDRBS1、MBNL1/2）
- **T0.5 快查**（无需等比对）：5 个核心 regulator（NUDT21、CPSF6、CSTF2、QKI、ELAVL1）在 mtx 中 astrocyte 簇的表达量先跑出来

### M2 APA 事件层（哪些 mRNA 的 APA 变了）——独立选源【v5.1r：统计闭环见 P2-06a SAP；n=2 定位为 discovery 层】
- 主源（条件成立时）：GSE174574 astrocyte pseudo-bulk → **scAPAtrap 或 Sierra**（3′-tag 感知）→ PDUI/近端 PAS 使用率，MCAO vs Sham，FDR<0.1
- Plan B（已触发执行中）：PB-1 GSE238125（bulk → DaPars2，多时间点）+ PB-2 GSE225110（**预固定验证集**——启用前必须冻结数据集身份/样本规则/模型/预期验证事件与成败标准，防 dataset shopping）
- 统计（v5.1r）：两层分析——样本级全局 APA（分层 + 敏感性 + 技术混杂排除）+ 事件级 moderated model（|ΔPDUI|≥0.1 预冻结、contrast 4/4 有效、联合 BH FDR<0.1），细则以 P2-06a 冻结版为准
- 判定（v5.1r）：**组合门槛**（样本 QC + 全局位移稳定 + 事件级判据 + top 事件抽验）；≥100 事件数保留为记录项

### M3 Regulator→target 关联（三级证据分级，禁止越级表述）
| 级别 | 证据要求 | 允许的表述 |
|------|----------|-----------|
| Level 1 | 已发表 perturbation 证据（KD/KO 实验证明该 regulator 改变靶基因 APA） | "regulator X 调控 gene Y 的 APA" |
| Level 2 | CLIP/binding 证据 + APA 方向一致性 | "与 regulator X 的结合证据一致" |
| Level 3 | motif 存在 + 表达/PDUI 相关 | "预测与 regulator X 相关" |
- 实操：每事件标注证据级别；文章中 Level 3 一律不写成因果。Level 1 文献池预检索（NUDT21/CFIm25→3′UTR 缩短、ELAVL1、QKI、PTBP 等已有经典 perturbation 工作）

### M4 Endfoot/PAP 富集层【v5.1r：三个互不同质注释轴，分别构建、分别报告、禁止投票合并】
- 轴 1 `PAP_localized`：GSE74456（外周突起/突触周过程，非血管周 endfoot）
- 轴 2 `endfoot_stroke_responsive`：GSE286075（卒中前后 endfoot 转译组响应，非相对 soma 的富集）
- 轴 3 `stroke_spatial_zone_response`：GSE263986（组织空间卒中分区，非亚细胞定位）
- 统计：Fisher exact + GSEA（按 ΔPDUI 排序）按轴分别输出与表述；"≥2 参考一致 = 高置信 endfoot"规则废除

### M6a 事件级证据矩阵【v5.1r：替代四重硬交集】
- 候选资格由 **APA 主证据**决定（事件效应/FDR/重复一致性/时间轨迹/复核）；机制注释（regulator concordance、L1/L2/L3 分池、CLIP、motif、三轴定位背景、MPRA、保守性）作优先级列，**不反向决定哪些 APA 结果被视为真实**
- 四重硬交集保留为描述性结果登记，不再作选择器

### M5 Localization evidence class（替代 Δscore）
- 对每个候选的 Δ3′UTR 区间做四层证据判定：

| 层 | 证据 | 数据 |
|----|------|------|
| 1. Experimental | Δ3′UTR 区域是否直接被 MPRA 验证为 localization element（**tile 序列精确映射到候选 Δ3′UTR 区间** + 活性方向；仅基因层面重合 = 背景注释，不入 class） | GSE330741（⛔ manifest 核对中） |
| 2. RBP | 是否有 CLIP 证据（QKI/KHDRBS1/ELAVL1…） | GSE146935、POSTAR3 |
| 3. Motif | 相关 RBP motif gain/loss | ATtRACT + FIMO |
| 4. Conservation | 改变区域/元件是否保守 | phastCons / UCSC 60-way |
- 输出类别：**A 级**（1+2+3+4 全有）、**B 级**（2+3+4，无直接 MPRA）、**C 级**（3+4）、**D 级**（仅 motif）
- 禁止：跨基因线性相加 MPRA 活性、发明合成数值分数

### M6b 整合交付
- 全链证据矩阵 + candidate 卡片 **0–5 张**（0 个满足预设条件 = 如实报零并报告失败层级，不得凑数；每张：事件、regulator 及其 M3 级别、evidence class、三轴富集统计）
- 湿实验延伸路线（v5.1r 按审计 §17 细化，围绕最上游候选如 PABPN1 的最短机制闭环）：3′RACE / 异构体特异 qPCR（长短比）→ regulator KD/rescue（测蛋白与核定位，非仅 mRNA）→ long/short UTR reporter + 元件删除/回补四构建体 → smFISH/RNAscope（血管 lectin/CD31 + AQP4/GFAP 定义 endfoot，报 endfoot/soma 比值）；资源受限时至少完成 3'RACE-qPCR / reporter / smFISH 三者之二
- 语言纪律全文统一：结论只到 "altered localization potential"；逐段对齐三级语言（直接数据支持 / 正交证据支持 / 预测），禁语清单见 TODO P5-05

---

## 4. Todolist v5.1（状态更新至 2026-09-16 16:00）

### 阶段 0：启动
- [~] T0.1 下载全部数据：GSE174574_RAW.tar 下载中（后台）；待续 GSE286075 tar、Plan B 池 SRA（按需）
- [x] T0.2 文库类型确认 → **10x 3′ droplet**（原文 + SRA 6/6 PAIRED/HiSeq4000/cDNA）；APA 工具锁定 scAPAtrap/Sierra；QAPA 仅用于 bulk Plan B —— **已关闭**
- [x] T0.3 坐标方案：mm10 基准（GENCODE mm10 + M2 PAS 注释同版本）
- [x] T0.4 项目结构：D:\stroke_apa（data/ 已建，pipeline 文档 v5.1 就位）
- [ ] T0.5 **regulator 快查**（新增，无需 FASTQ）：tar 解包 → 读 mtx → 5 核心 regulator（NUDT21/CPSF6/CSTF2/QKI/ELAVL1）在 astrocyte 簇的表达快查（审阅意见第 4 项）
- [ ] T0.6 **细胞数审计**（审阅意见第 2 项）：每样本 barcodes 数 + astrocyte marker 粗注释 → 关闭 kill-switch A 的判定输入
- [ ] T0.7（可选）Plan B 池 GSE238125 SRA 预下载（12 runs，体量小），作为保险

### 阶段 1：M1 regulator 层（1 周）
- [ ] T1.1 GSE174574 QC + 聚类 + astrocyte 簇定义（交付：UMAP + marker 图）
- [ ] T1.2 homeostatic/reactive 拆分 + pseudo-bulk DE（交付：regulator DE 表）
- [ ] T1.3 GSE286075 配对 DE（交付：paired DE 表）
- [ ] T1.4 双层交叉 + 分级目录（交付：M1 模块报告）
- **GATE 1**：≥3 个 regulator 两层方向一致；否则 M3 改纯 CLIP/motif 途径并明示

### 阶段 2：M2 APA 事件层（1 周）
- [ ] T2.1 触发判定：T0.6 astrocyte 细胞数 ≥500/样本？scAPAtrap 试跑 1 样本评估 3′ 端信息质量
- [ ] T2.2 主源 APA（GSE174574 pseudo-bulk + scAPAtrap/Sierra）或 Plan B（GSE238125 bulk + DaPars2）
- [ ] T2.3 FDR<0.1 + 方向分类 + （若双源）跨数据集方向一致性
- **GATE B**：显著事件 ≥100 继续；<100 → 探索性叙事或换 PB-2 复试一次

### 阶段 3：M3+M4 关联与富集（1 周）
- [ ] T3.1 QKI CLIP peaks → BED + intersect（Level 2 证据池）
- [ ] T3.2 Level 1 文献池整理（perturbation 证据的 regulator-target 对）
- [ ] T3.3 ATtRACT + FIMO motif 扫描（Level 3 证据池）
- [ ] T3.4 endfoot 三参考构建（GSE74456 + GSE286075 + GSE263986）+ Fisher/GSEA

### 阶段 4：M6a 漏斗 + M5 evidence class（1 周）
- [ ] T4.1 四重交集计算 → endfoot-APA candidates 池
- [ ] T4.2 Δ3′UTR 区间提取 → MPRA tile intersect → CLIP intersect → motif gain/loss → conservation
- [ ] T4.3 A/B/C/D evidence class 判定表（每 candidate 一张证据矩阵）

### 阶段 5：M6b 整合交付（3–5 天）
- [ ] T5.1 全链证据矩阵 + 图
- [ ] T5.2 top 3–5 candidate 卡片（含 M3 级别 + M5 class）
- [ ] T5.3 湿实验验证路线书
- [ ] T5.4 报告撰写（局限：n=3、mRNA≠蛋白、RiboTag proxy、GSE238125 n=2）

**周期估计：5–6 周（纯干）不变**；kill-switch A 已由"单一条件"细化为"细胞数 + 工具输出质量"双条件，且 Plan B 已从概念变成 4 套已核实数据集。

---

## 5. 风险表 v5.1（增量）

| 风险 | 概率 | 影响 | v5.1 缓解 |
|------|------|------|-----------|
| GSE174574 astrocyte 细胞不足 / 10x 3′ APA 信号弱 | **已确认为真风险**（10x 3′ 文库） | 高 | T0.5/T0.6 先行审计；PB-1 GSE238125（bulk，多时间点）为主切换目标；PB-2 交叉 |
| GSE238125 n=2/时间点功效弱 | 中 | 中 | 效应量 + d1–d60 方向一致性 + PB-2 交叉验证；预注册判定规则 |
| M3 过度解释 | 高（写作风险） | 中 | 三级证据分级 + 语言纪律表（§3 M3） |
| M5 合成分数被审稿质疑 | 中 | 中 | Δscore 已废除；evidence class 为可枚举证据，不引入权重 |
| 10x 矩阵无 QC 状态说明 | 低 | 低 | tar 内含 mtx 原始计数，自行 QC（Scrublet + 线粒体比例） |
| （继承 v5.0）n=3 伪重复 / RiboTag proxy 误读 / 坐标漂移 / regulator mRNA≠活性 | — | — | 措施同 v5.0（§6 风险表、双参考、liftOver、验证建议） |

---

## 6. 迭代轨迹（v5.0 五轮 + v5.1 两轮）

| 轮次 | 发现 | 处置 |
|------|------|------|
| v5.0 Loop 1–5 | （见 v5.0 报告 §6：数据核查→链条闭合→流程落地→统计强度→输出收敛） | 生成 v5.0 |
| **v5.1 Loop 6：外部审阅** | R1 M1/M2 耦合单点故障；R2 M3 过度解释；R3 M5 Δscore 不成立；R4 缺反向支线；R5 应立即执行 T0 | 解耦、分级、废 Δscore、立漏斗、转执行 |
| **v5.1 Loop 7：T0 实查** | GSE174574 实为 10x 3′（QAPA 出局）；10x 三件套矩阵存在（regulator 快查无需比对）；4 套 Plan B 全部核实有 SRA；单文件 supplements 不存在（在 tar 内） | 工具链改写；Plan B 池落位；T0.5/T0.6 增补 |
| **v5.1r Loop 8：论证逻辑审计（2026-09-18 报告，09-19 用户采纳全文）** | 因果链过强；P2 统计未闭环；三参考异质被投票合并；QKI/MPRA 证据越界；Gate 设计缺陷；多处事实错位 | §0b A1–A14 落实；判定 = Major revision / 有条件 Go，P2 全量继续 |

---

## 7. 关键引用（v5.1 增补）

1–6 同 v5.0（GSE286075/GSE74456/GSE146935/GSE330741/Nat Commun 2018/3'aTWAS 2025）
7. Zheng K, et al. Single-cell RNA-seq reveals the transcriptional landscape in ischemic stroke. *J Cereb Blood Flow Metab.* 2022;42(1):56-73. PMID 34496660（= GSE174574，**10x 3′ 建库证据源**）
8. GSE238125：Dynamic astrocytic gene expression profiles after ischemic stroke（PRJNA997998）
9. GSE225110：Glial translatomes in stroke（PRJNA922165）
10. GSE263986：Astrocyte-enriched murine cortical and white matter transcriptional profiles after stroke in spatially-defined zones（PRJNA1100502）
11. GSE165384：Quantitative Analysis of sham and cerebral ischemia Transcriptomes in mouse hippocampus（SRP302969）
