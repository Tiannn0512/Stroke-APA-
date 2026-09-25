# GSE330741 补充文件清单核验（manifest check）

- 核验日期：2026-09-25（本机时间 Sep 25 11:08–11:50）
- 核验人：自动复核（stroke_apa 项目 07_followup）
- 下载缓存目录：`/d/stroke_apa_reaudit/results/07_followup_20260925/notes/gse_cache/`
- 证据类别约定：**[直接观察]** = 公共文件逐字/逐行核对；**[推断]** = 由公开文件计算或按外部文献对照得出；**[待独立实验]** = 公共数据不足以完成，需自行实验或联系作者。

---

## 1. Series 基本信息 [直接观察]

来源：`GSE330741_series_matrix.txt.gz`（下载自 `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE330nnn/GSE330741/matrix/GSE330741_series_matrix.txt.gz`，2026-09-25）

| 字段 | 值（逐字摘录） |
|---|---|
| Series_title | "In Vivo Massively Parallel Reporter Assay Reveals Sequence Determinants of mRNA Localization in Astrocytes" |
| Series_type | "Expression profiling by high throughput sequencing" |
| Series_status | "Public on May 17 2026"（submission 2026-05-12，last update 2026-09-09） |
| Series_contributor | "Sarah,K,Koester" / "Kristina,,Sakers" / "Gareth,M,Rurak" / "Stephen,P,Plassmeyer" / "Kelli,,McFarland White" / "Sebastian,,Alves Ferreira Dias" / "Peter,,Baird" / "Dina,J,Kornbluth" / "Joseph,D,Dougherty"（联系人 Joseph D. Dougherty, Genetics, Washington University in St. Louis） |
| Platform | GPL34290 = "Illumina NovaSeq X Plus (Mus musculus)"（soft 文件核对） |
| Assembly | "Assembly: GRCm38 (mm10)" [直接观察，series matrix data_processing 行] |
| 关联 | "Reanalysis of: GSM9188124 … GSM9188139"（16 条）；BioProject PRJNA1465359 |
| 更新说明 | "September 08, 2026: Samples GSM9822133-GSM9822156 were moved to Series GSE295080, per submitter request (Reason: The Samples were included by mistake in Series GSE330741)." |

**论文判定**：疑为 Koester 等 2026 年 SN-MPRA 论文属实（submitter = Koester/Sakers/Dougherty，WashU）。摘要仅提到 **Glt1 (Slc1a2) 与 Sparc 两个 3'UTR 的 tiling**，未提 Atp2a2。尚无 DOI/PMID（GEO 页未挂引用，2026-09-25 查询）。

**样本设计** [直接观察]：AAV9 递送报告文库感染 Aldh1l1-EGFP/RPL10a（JD130Htz）小鼠，测：
- `aav_1..3`（treatment: "Library AAV DNA"）= 文库 DNA 对照；
- `ctxin_*` / `snin_*`（treatment: "MPRA Library AAV Infection"，tissue "brain"）= 输入 RNA；
- `ctxtrap_*` / `sntrap_*` = GFP(RPL10a)-IP（TRAP）。
唯一文库编号（按文件名去重）：ctxin n=16、snin n=16、ctxtrap n=9、sntrap n=9、aav n=3；同编号的 `a`/`b` 后缀文件为重复测序（重跑）。文件总数 139（RAW.tar）。
**样本单位**：按文库编号推断为单只动物；GEO 未提供 sample sheet，in/trap 配对仅靠编号对应 [推断]。
**SN 含义**：series matrix data_processing 定义 localization=log2(sn input/ctx input)、local translation=log2(SN TRAP/ctx TRAP)；结合 Sakers 2017 PNAS 同组用法（SN=synaptoneurosome，突触神经体组分），判 SN=突触神经体组分 [推断，GEO 元数据只写 tissue "brain"]。

## 2. 补充文件清单与格式 [直接观察]

suppl 目录：`https://ftp.ncbi.nlm.nih.gov/geo/series/GSE330nnn/GSE330741/suppl/`（2026-09-25 列目录 + `filelist.txt`）

### 2.1 Pre/Post/DNA 计数文件（32 个，GSE330741_{DNA,Pre,Post}_{1-4/1-6}_{L1,L2}.counts.txt.gz）
- 格式：**无表头** 4 列，tab 分隔，1631 行/文件。已核对 `GSE330741_DNA_1_L1.counts.txt.gz`、`GSE330741_Pre_1_L1.counts.txt.gz`、`GSE330741_Post_1_L1.counts.txt.gz` 前 5 行，形如：
  ```
  9-Sep_wgs_alt	150	2032	0
  9-Sep_wgs_ref	150	2386	0
  ABCB5_wgs_alt	150	470	0
  ```
  推断列 = elementID、寡核苷酸长度(150)、计数、常数 0 [推断]。
- 元素名为**人类基因符号 + wgs/ssc + alt/ref/shuf**（ABCB5、ABHD2、"9-Sep"（即 SEPT9 的 Excel 日期污染写法）、AC010536.1 等，去重后约 979 个基因）。
- **与主文库完全不相交**：与 tiling 计数文件的 6430 个元素 ID 交集 = 0 [直接观察，comm 校验]。
- **判定**：这是一套来源不明的外来（人类变异型 MPRA 风格）计数文件，与 series 的 astrocyte tiling 文库无法互相对应；疑似与被移出/Reanalysis 的 GSM9188124-9188139 或 GSE295080 相关，但 GEO 未作任何说明 → **manifest 不一致，待作者澄清**。

### 2.2 RAW.tar（139 个 GSM 文件，GSM9731591–GSM9731729）
- 格式：**2 列带表头**（"ID\tcount"），6431 行/文件（6430 元素）。已核对 `GSM9731603_ctxtrap_1.txt.gz`、`GSM9731591_aav_1.txt.gz`、`GSM9731594_ctxin_1.txt.gz`。
- 元素 ID 家族及片段数（按"去掉末尾 barcode 序号的 ID"去重）：
  | 家族 | 元素数 | 片段数 | 说明 |
  |---|---|---|---|
  | Slc1a2_* | 4570 | 457 | Glt1/Slc1a2 3'UTR tiles（ID 形如 Slc1a2_1_1001_10：gene_?_tile_barcodeIdx） |
  | sparc_* | 660 | 66 | Sparc 3'UTR tiles（ID 形如 sparc_1001_1：gene_tile_barcodeIdx） |
  | scram_* | 660 | 66 | scrambled 对照 |
  | hsbp1_* | 400 | 40 | Hsbp1 相关 tiles |
  | Basal_* | 110 | 11 | 基础启动子对照 |
  | mbp_loc_*/mbp_qk_* | 20 | 2 | 定位阳性对照（mbp） |
  | act_zip_* | 10 | 1 | zipcode 对照 |
  | 合计 | **6430** | **643** | 每片段 ~10 个 barcode 序号 |

## 3. barcode → insert → sample 连接键核验

| 键 | 是否公开 | 说明 |
|---|---|---|
| sample → count 文件 | ✔ [直接观察] | filelist.txt 的 GSM 号 ↔ 文件名一一对应；series matrix 无 per-sample 样本表（无动物号/年龄/批次） |
| elementID → count | ✔ [直接观察] | 每个 GSM 文件 2 列 |
| elementID → barcode 序列 | ✘ | ID 只含"barcode 序号"（每片段 10 个），序列未公开 |
| elementID → insert 片段序列 / 基因组坐标 | ✘ | data_processing 提到 "aligned to reference libraries using custom scripts"，reference library 文件未随 GEO 存档 |
| tile 编号 → Glt1/Sparc 3'UTR 坐标 | ✘ | 无 data dictionary；tile_1001 等编号无法映射到核苷酸区间 |
| sample → 动物/批次配对表 | ✘ | 仅能按文件名编号推断 |

**结论 [直接观察+推断]**：当前公开文件可以重建"片段（fragment-level，聚合 10 个 barcode 序号）→ 各样本计数 → 按作者公式算 expression/localization/local translation 的相对排序"，**但**：(a) 无法做独立的质量控制（barcode 饱和度、跨 barcode 一致性）；(b) 无法把任何片段对应到具体核苷酸序列或 3'UTR 坐标；(c) 无法把结果映射到 Atp2a2 —— **文库中根本没有 Atp2a2 元素**（对全部 6430 个元素 ID 与所有已下载 txt.gz/counts.txt.gz grep "atp2a2|ENSMUSG00000029467" = 0 命中）。

## 4. 判定与最小缺失输入（status: pending）

- **Atp2a2 reporter 验证：不可能由本数据集提供**（非"数据不全"，而是"文库不含该基因"）。
- 若要完整重建"片段 → reporter 结果"映射（针对 Glt1/Sparc），最小缺失输入：
  1. **参考文库表**（每寡核苷酸的 ID、完整 insert 序列、来源基因/3'UTR 坐标、barcode 序列，含 BGH/attB 等骨架信息）——用途：把 elementID 映射到核苷酸片段并独立重对齐 reads；
  2. ** authors 的 "custom scripts"**（exact-match 对齐、CPM 归一、低计数过滤、离群剔除的代码）——用途：复现 expression/localization/local translation 指标；
  3. **样本表**（动物号 ↔ 文库名 ↔ region/trap 配对、批 次、测序 lane L1/L2 归属）——用途：确定生物学重复与配对结构；
  4. **Pre/Post/DNA 32 个文件的归属说明**（series 与主文库 0 交集）——用途：澄清 manifest 是否受 2026-09-08 移样事件影响。
- 获取途径：邮件联系 Dougherty lab（WashU Genetics）；或等待论文正式发表后的 supplementary。
- **不得把本数据集的 Glt1/Sparc 片段结果解释为 Atp2a2 的验证**；如需 Atp2a2 3'UTR 片段证据，属**待独立实验**（例如把 Atp2a2 3'UTR tiles 克隆进同一 SN-MPRA 骨架后自测）。
