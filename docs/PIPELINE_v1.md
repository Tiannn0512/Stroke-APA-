# 卒中→APA→星形胶质细胞 endfoot 定位潜能 小课题：Pipeline v5.0 与 Todolist

> 生成日期：2026-09-16 ｜ 状态：可行性核查完成，五轮自检收敛
> 核心假说：Stroke → APA regulator 改变 → APA 改变 → 特定 astrocyte mRNA 3′UTR 亚型改变 → endfoot-localized mRNA 中的 APA 事件 → 3′UTR 序列（element）删除/获得 → localization-related cis-element / RBP-binding element 改变 → endfoot 定位潜能改变（预测）

---

## 0. 可行性判定（结论先行）

**判定：方向可行（纯干限定位），但原始设想有 4 处必须修正。**

| # | 原设想 | 核查结果 | 修正 |
|---|--------|----------|------|
| C1 | 用 GSE286075 检测/确认 APA | GSE286075 = Shim 2025 (IJMS, PMID 39796165)，Astrocyte endfoot **RiboTag IP RNA-seq**（MCAO 2h缺血+6h再灌注，ipsi vs contra 各 n=3），GEO 补充文件**只有基因水平 counts（ReadsPerGene.out.tab），无 FASTQ/BAM** | GSE286075 重定位为「endfoot translatome 参考 + regulator 第二层 DE」；APA 检测改用带原始数据的卒中 scRNA（GSE174574，SRP320164 有 FASTQ） |
| C2 | 只用 GSE286075 闭合全链 | 单数据集既不能测 APA，也没有 localization element 参考 | 引入三套互补数据：GSE74456（PAPTRAP，PAP 定位富集参考）、GSE146935（QKI CLIP-seq，RBP element 参考）、GSE330741（SN-MPRA，in vivo localization element 金标准） |
| C3 | "确认" endfoot 定位改变 | 手头只有 IP/转译组 proxy，无原位定位数据 | 措辞降级为「预测 endfoot localization 潜能」（与用户假说原意一致）；原位验证（RNAscope/smFISH）列入湿实验延伸 |
| C4 | regulator 改变 = 活性改变 | 手头证据是 mRNA 丰度（RiboTag 测的也是 mRNA），不是蛋白量/翻译后修饰 | 报告中明确标注局限；top candidates 给 Western/IF 验证建议 |

**查新**：stress→APA（ arsenite 全局 3′UTR 缩短，Nat Commun 2018）；3'aTWAS 已把缺血性卒中列为 APA 相关脑病（2025）；但「卒中→astrocyte APA→endfoot 定位潜能」链无发表。新颖性成立。

---

## 1. 数据资产（已逐一在 GEO 核实）

| 登录号 | 内容 | 关键参数 | 角色 | 原始数据 |
|--------|------|----------|------|----------|
| GSE286075 | Astrocyte endfoot RiboTag RNA-seq（Shim 2025） | MCAO 2h+6h 再灌注；ipsi/contra 同鼠配对 n=3+3；Mus musculus；GPL24247 | endfoot translatome 参考；regulator 第二层 DE（配对检验） | ❌ 仅 gene counts |
| GSE174574 | 卒中脑 scRNA-seq（2021） | MCAO vs Sham n=3+3 | astrocyte 分群 → regulator DE 第一层 + **APA 检测源（pseudo-bulk）** | ✅ SRA: SRP320164 |
| GSE74456 | PAPTRAP（Sakers 2017 PNAS, PMID 28439016） | 外周星形胶质细胞突起（PAP）核糖体结合 mRNA；PAP/Ctx × TRAP/input | **PAP/endfoot 定位富集参考**（双参考之一） | ✅ SRA: SRP065508 |
| GSE146935 | QKI CLIP-seq + CRISPR-TRAPseq（PMID 33750804） | QKI-6 结合 astrocytic mRNA 3′UTR，peak 富集于终止密码子附近 | **RBP-binding element 参考**（QKI 通路） | ✅ SRA: SRP252704 |
| GSE330741 | SN-MPRA（2026） | in vivo MPRA tiling Glt1/Sparc 3′UTR；定位 element 单碱基突变解析；KHDRBS1 等 motif | **localization element 金标准参考** | ✅ Bioproject PRJNA1465359 |

辅助数据库：GENCODE mouse 注释（统一坐标）；PolyASite 2.0 / APAatlas（PAS 目录）；POSTAR3、DoRINA（RBP–target）；ATtRACT（motif + FIMO）。

---

## 2. 因果链逐环证据矩阵（最终版）

| 环节 | 证据来源 | 数据支撑 | 强度 | 缺口 |
|------|----------|----------|------|------|
| Stroke → regulator 改变 | GSE174574 astrocyte 簇 pseudo-bulk DE；GSE286075 配对 DE | 直接 | ★★☆（mRNA 层） | 蛋白/活性未测 |
| regulator → APA 改变 | GSE174574 astrocyte APA 事件 × regulator motif/CLIP/文献 | 关联推断 | ★★☆ | 无因果操作（KD） |
| APA → 3′UTR 亚型改变 | PDUI/PAS-usage 定量 | 直接 | ★★★ | — |
| 3′UTR 改变 × endfoot 富集 | GSE74456 PAP set + GSE286075 endfoot set 双参考富集 | 富集统计 | ★★☆ | 原位定位未测 |
| Δ3′UTR × localization element | GSE330741 MPRA tiles + GSE146935 QKI CLIP + ATtRACT motif gain/loss | 序列证据 | ★★☆ | element 功能未单独验证 |
| → 定位潜能改变 | Δscore 模型（MPRA activity + motif 加权） | 预测 | ★☆☆（模型输出） | 需湿实验闭环 |

---

## 3. 最终 Pipeline（v5.0，六模块）

### M0 数据与工具底座
- 下载：GSE174574 FASTQ（SRP320164）、GSE74456（SRP065508）、GSE146935 CLIP（SRP252704）、GSE286075 counts tar。
- 注释：GENCODE mouse（M2 检测用 Gencode PAS；统一 liftOver 到 mm10 或 mm39，PAPTRAP/GSE146935 时代较早，注意坐标版本）。
- 工具链：Cell Ranger / STARsolo（scRNA 比对定量）；STAR + Samtools（pseudo-bulk BAM）；QAPA 或 DaPars2（bulk APA 主工具）+ scAPAtrap / Sierra（10x 专用交叉验证，视文库类型启用）；bedtools；FIMO + ATtRACT；R（DESeq2/edgeR/clusterProfiler）。
- GSE174574 文库类型确认（Smart-seq2 full-length vs 10x 3′）：决定 APA 主工具。**若 full-length → QAPA 为主；若 10x → scAPAtrap/Sierra 为主 + 谨慎解释。**

### M1 APA regulator 层（谁变了）
- GSE174574：Scrublet 去双胞 → 聚类 → 识别 astrocyte（Aqp4/Slc1a3/Gfap/Aldh1l1）→ 区分 homeostatic vs reactive（C3/Serpina3n/Gfap^high）→ pseudo-bulk（按样本×亚群）→ DESeq2 DE（MCAO vs Sham）。
- 检索范围分级 APA regulator 名单：
  - 核心剪切复合物：CPSF1/2/3/4、WDR33、FIP1L1、CSTF1/2/3、SYMPLEK、NUDT21(CFIm25)、CPSF6/7、PABPN1、PCF11、CLP1；
  - RBP 型 modulator：ELAVL1/2、QKI、PTBP1/2、HNRNPA2B1、SRSF1、NOVA1/2、FUS、TARDBP、KHDRBS1、MBNL1/2。
- GSE286075：同名单 ipsi vs contra 配对 DE（limma/DESeq2 paired design）→ 双层证据交叉。
- 产物：regulator 变化目录（表达变化 + 分级 + 双数据集方向一致性）。

### M2 APA 事件层（哪些 mRNA 的 APA 变了）
- GSE174574 astrocyte 细胞按样本 pseudo-bulk → BAM → QAPA/DaPars2（Gencode PAS 注释）→ 每基因 PDUI / distal PAS usage，MCAO vs Sham，FDR<0.1。
- 交叉验证：文库允许则 scAPAtrap/Sierra 单细胞级复核；方向一致事件标为高置信。
- 产物：astrocyte APA 事件目录（3′UTR 缩短/延长；ΔPDUI 效应量优先报告）。

### M3 regulator→target 关联层（谁能解释这些 APA）
- 对每个 APA-altered 基因的近端/远端 PAS 区间 + Δ3′UTR 区域：
  - intersect QKI CLIP peaks（GSE146935）；
  - FIMO 扫描 M1 中显著变化 regulator 的 motif（ATtRACT）；
  - POSTAR3 / DoRINA / 文献已知调控对；
  - 统计：regulator Δexpr × ΔPDUI 方向一致性（如 NUDT21↓ ↔ 3′UTR 缩短族）。
- 产物：regulator–target 关联矩阵（每事件标注候选 regulator 与证据类型）。

### M4 endfoot/PAP 富集层（这些 mRNA 在 endfoot 富集吗）
- 双参考定义：GSE74456 PAP-enriched set（PAP TRAP vs Ctx TRAP，按 Sakers 2017 定义复算或取文中列表）；GSE286075 endfoot translatome（全部检测基因 + 205 DE 基因）。
- 检验：APA-altered gene set vs 两套参考 → Fisher exact + GSEA（按 ΔPDUI 排序打分）；两参考一致者为高置信 endfoot-APA 交集。
- 产物：stroke-APA × endfoot 基因列表（附富集统计与效应方向）。

### M5 element 与定位潜能层（改变的序列里有什么）
- 对每个高置信 endfoot-APA 事件提取 Δ3′UTR（获得的/丢失的 [proximal PAS, distal PAS] 区间）：
  - intersect GSE330741 MPRA localization tiles（Glt1/Sparc tiling，及可推广的 motif 群）；
  - QKI/KHDRBS1/ELAVL1 等 motif 的 gain/loss 计数（FIMO，统一坐标）；
  - 构 localization potential Δscore：Δscore = Σ(获得的 MPRA tile 活性/motif 权重) − Σ(丢失的)。
- 产物：事件级 element gain/loss 表 + localization potential 变化排序。

### M6 整合与优先级（交付）
- 全链证据矩阵（每环节证据/强度/缺口）+ top 3–5 candidate mRNA（APA 改变 + 双参考 endfoot 富集 + element gain/loss 明确 + 对应 regulator 变化）。
- 湿实验验证路线（延伸，不属本小课题）：3′RACE 验证 PAS switch；RNAscope/smFISH 定位 endfoot；候选 element 单独 MPRA；regulator 干预（AAV-shRNA / 星形胶质细胞特异）。

---

## 4. Todolist（按模块，含判定标准与 kill-switch）

### 阶段 0：启动（1–2 天）
- [ ] T0.1 申请/下载全部 4 个 GSE 数据（含 SRA）与 Gencode 注释
- [ ] T0.2 确认 GSE174574 文库类型（读 SRA metadata / 文献方法学）→ 锁定 APA 主工具
- [ ] T0.3 统一坐标方案（建议 mm10 为基准；GSE330741 若 mm39 则 liftOver）
- [ ] T0.4 建项目 repo（数据/脚本/结果/笔记四目录）与环境（conda/Renv）
- **Kill-switch A**：GSE174574 若为纯 10x 且 astrocyte 细胞数 < 500/样本 → 切换备选 APA 源（预检索 2–3 个 bulk MCAO RNA-seq with FASTQ 作为 Plan B）

### 阶段 1：regulator 层（1 周）
- [ ] T1.1 GSE174574 QC + 聚类 + astrocyte 簇定义（交付：UMAP + marker 图）
- [ ] T1.2 astrocyte homeostatic/reactive 拆分 + pseudo-bulk DE（交付：regulator DE 表）
- [ ] T1.3 GSE286075 配对 DE（交付：paired DE 表）
- [ ] T1.4 双层交叉 + regulator 分级目录（交付：M1 模块报告）
- **判定**：至少 3 个核心/调制型 regulator 在两层证据中方向一致，否则 M3 改用纯 motif/CLIP 途径并在报告中说明

### 阶段 2：APA 事件层（1 周）
- [ ] T2.1 pseudo-bulk BAM 生成与 QC（比对率/3′覆盖评估）
- [ ] T2.2 QAPA/DaPars2（或 scAPAtrap/Sierra）跑通 → PDUI 全表
- [ ] T2.3 FDR<0.1 事件筛选 + 方向分类（缩短/延长）
- [ ] T2.4 可行性复盘点：显著事件数 ≥ 100 才继续全链；< 100 → 降级为探索性报告（kill-switch B）
- **产物**：astrocyte APA 事件目录

### 阶段 3：关联与富集（1 周）
- [ ] T3.1 QKI CLIP peaks 处理为 BED + intersect（交付：CLIP 支持的 APA 事件子集）
- [ ] T3.2 ATtRACT motif + FIMO 扫描（交付：motif gain/loss 表）
- [ ] T3.3 PAPTRAP PAP set 构建 + GSE286075 endfoot set 构建
- [ ] T3.4 Fisher/GSEA 双参考富集（交付：endfoot-APA 基因列表 + 富集统计）

### 阶段 4：element 与定位潜能（1 周）
- [ ] T4.1 GSE330741 MPRA tiles/活性表整理（BED + score）
- [ ] T4.2 Δ3′UTR 区间提取 + tile intersect + motif gain/loss
- [ ] T4.3 Δscore 模型构建与排序（交付：事件级 localization potential 表）

### 阶段 5：整合交付（3–5 天）
- [ ] T5.1 全链证据矩阵 + 图（每环节）
- [ ] T5.2 top 3–5 candidate mRNA 卡片（事件、regulator、element、富集四要素）
- [ ] T5.3 湿实验验证路线书（延伸部分）
- [ ] T5.4 报告撰写（含局限：n=3、mRNA≠蛋白、proxy 定位）

**总周期估计：5–6 周（纯干）**；两处 kill-switch（A 数据源、B 事件数量）保证不沉没。

---

## 5. 五轮自检迭代轨迹（收敛记录）

| 轮次 | 发现的问题 | 修正动作 |
|------|------------|----------|
| Loop 1 数据核查 | ① GSE286075 无 FASTQ，不能测 APA；② 因果链需要"regulator 检测层"与"APA 检测层"两个数据层 | GSE286075 → 参考层；引入 GSE174574 作 APA 检测源（有 FASTQ） |
| Loop 2 链条闭合 | ③ 定位只能预测不能确认；④ TRAP/RiboTag ≠ total localization；⑤ regulator 名单需分级 | 全链措辞降级为"潜能预测"；引入 PAPTRAP + CLIP + MPRA 三套数据闭合；regulator 名单分核心/调制两级 |
| Loop 3 流程落地 | ⑥ GSE174574 文库类型待定 → 工具分叉；⑦ astrocyte 需分 reactive/homeostatic；⑧ 坐标系跨数据集不一致 | M2 主备双工具策略；M1 分层 DE；M0 统一 liftOver |
| Loop 4 统计强度 | ⑨ n=3 功效有限；⑩ 富集单一参考易偏；⑪ 干湿边界模糊 | 效应量+跨数据集方向一致性优先；M4 双参考；M6 明确湿实验为延伸 |
| Loop 5 输出收敛 | ⑫ 缺可执行判定标准；⑬ 无风险预案 | Todolist 每阶段加判定标准；设 kill-switch A/B；锁定交付物清单 |

---

## 6. 风险表

| 风险 | 概率 | 影响 | 缓解 |
|------|------|------|------|
| GSE174574 文库不适合 APA | 中 | 高 | kill-switch A；Plan B：bulk MCAO RNA-seq（预检索备选）或 3′tag 数据 |
| astrocyte APA 显著事件过少 | 中 | 高 | kill-switch B（≥100）；或放宽至探索性报告 |
| n=3 伪重复/功效 | 高 | 中 | 效应量优先；scRNA 细胞级交叉；方向一致性报告 |
| TRAP proxy 误读为定位 | 中 | 中 | 双参考交叉 + 措辞纪律（"预测潜能"） |
| 坐标/注释版本漂移 | 中 | 低 | liftOver + 事件级坐标复核 |
| regulator mRNA≠活性 | 高 | 中 | 明示局限；top 候选附验证建议 |

---

## 7. 关键引用

1. Shim B, et al. RiboTag RNA Sequencing Identifies Local Translation of HSP70 in Astrocyte Endfeet After Cerebral Ischemia. *Int J Mol Sci.* 2025;26(1):309. PMID 39796165.（GSE286075）
2. Sakers K, et al. Astrocytes locally translate transcripts in their peripheral processes. *PNAS.* 2017;114(22):E3830-E3838. PMID 28439016.（GSE74456）
3. QKI CLIP/TRAP: PMID 33750804.（GSE146935）
4. SN-MPRA: In vivo MPRA reveals sequence determinants of mRNA localization in astrocytes.（GSE330741, 2026）
5. Cellular stress alters 3′UTR landscape through alternative polyadenylation and isoform-specific degradation. *Nat Commun.* 2018.（stress→APA 依据）
6. Single-cell APA TWAS (3'aTWAS) for brain disorders incl. ischemic stroke. *PLOS Genet/bioRxiv.* 2025.（卒中-APA 关联依据）
