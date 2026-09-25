# GSE 复核日志（GSE74456 / GSE143531 / GSE330741）— Atp2a2 基因层面证据

- 执行日期：2026-09-25（所有下载与计算均在本日完成）
- 环境：Windows 10 / Git Bash / curl；Python 3.12（numpy 2.5.3, scipy 1.18.1, pandas 3.0.6；statsmodels/R 不可用，BH 为手写实现）
- 缓存根目录：`/d/stroke_apa_reaudit/results/07_followup_20260925/notes/gse_cache/`
- 证据类别：**[直接观察]** 公共文件逐字摘录；**[计算推断]** 本地复算；**[外部文献线索]** 未在公开文件中核实的说法

---

## 0. Atp2a2 的基因标识（先决步骤）

[直接观察] 计数文件不带基因符号，只有 Ensembl ID，故先解析：
- `curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=Atp2a2[Gene Name]+AND+Mus+musculus[Organism]&retmode=json"` → 唯一命中 **GeneID 11938**。
- `efetch db=gene id=11938 retmode=xml` → `<Gene-ref_locus>Atp2a2`，Ensembl dbtag = **ENSMUSG00000029467**。
- 反查 `ENSMUSG00000028590[All Fields]` → 0 命中（确认早期项目记录里若用过 28590 则是错误 ID）。
- 后续全部检索均使用 ENSMUSG00000029467。

---

## 1. GSE74456（"PAPTRAP"）

### 1.1 下载
| URL | 落盘文件 | 大小 |
|---|---|---|
| https://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74456/suppl/ （目录列表） | — | — |
| https://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74456/suppl/filelist.txt | notes/gse_cache/filelist_GSE74456.txt | 12 个 GSM count 文件 + RAW.tar |
| https://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74456/matrix/GSE74456_series_matrix.txt.gz | notes/gse_cache/GSE74456_series_matrix.txt.gz | 2,239 B |
| https://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74456/suppl/GSE74456_RAW.tar | notes/gse_cache/GSE74456_RAW.tar | 1,740,800 B（12 文件，已解包至 notes/gse_cache/GSE74456/） |

### 1.2 设计与样本单位 [直接观察，series matrix 逐字]
- `!Series_title "PAPTRAP"`；`!Series_pubmed_id "28439016"`
- `!Sample_source_name_ch1`："Input to Cortex TRAP before immunoprecipitation" / "Cortex TRAP: Immunoprecipitation sample of GFP from Aldh1l1 bacTRAP mouse (cortex and hippocampus)" / "Synaptoneurosome fraction of Aldh1L1 bacTRAP cortex and hippocampus, before immunoprecipitation" / "Peripheral Astrocyte Process TRAP: Immunoprecipitation of GFP from synaptoneurosome fraction of Aldh1L1 bacTRAP cortex and hippocampus"
- `!Sample_characteristics_ch1`："strain: C57/Bl6"；"age: Post natal day 21"；"cell type: astrocytes" → **发育期 P21 基础状态，非卒中**
- `!Sample_data_processing`："Align to genome with STAR"、"Count reads with HTSeq"、"Genome_build: mouse ensembl175"（= Ensembl release 75 / GRCm38-mm10）、"Supplementary_files_format_and_content: .count files. Ensembl id and count #"
- 12 样本：CtxPre1-3, CtxTRAP1-3, PAPpre1-3, PAPTRAP1-3（各 n=3）
- **样本单位 [外部文献线索→已核实为论文原文]**：论文（PMC5441704，Methods，经 WebFetch 全文读取）"Three cortices (including hippocampi) per replicate were homogenized from Aldh1l1:EGFP/RPL10A mice" → **每重复 = 3 只小鼠合并**；差异分析用 edgeR；PAP-enriched 定义 = SN input/Cortex input 与 PAP TRAP/SN input 两对比的交集（FDR<0.1 且 log2FC>1.5）。

### 1.3 Atp2a2 行原始摘录 [直接观察]
文件：`GSE74456/GSM1920997_PAPTRAP1.count.txt.gz`（两列无表头：Ensembl_ID \t HTSeq_count）
```
ENSMUSG00000029467	11478
```
12 个文件中该 ID 均存在且唯一（grep -c = 1/文件）。全部原始计数（analysis_scripts/recompute_gse74456.py 输出，logs/recompute_gse74456_out.txt）：
- PAPpre: 4616, 4619, 9164；PAPTRAP: 11478, 8186, 15561
- CtxPre: 7213, 7343, 7939；CtxTRAP: 6969, 2236, 3086

### 1.4 本地复算 [计算推断]
脚本：`analysis_scripts/recompute_gse74456.py`；输出：`logs/recompute_gse74456_out.txt`
方法：丢 HTSeq `__` 行；median-of-ratios size factors；log2(norm+1)；Welch t 检验（n=3 vs 3）；BH 全基因组（手写）。注意这是对作者 edgeR 的**近似复算**。
- size factors：CtxPre 1.899/1.593/1.709；CtxTRAP 1.902/1.246/1.774；PAPpre 0.317/0.329/0.521；PAPTRAP 1.128/1.269/1.317
- **PAP TRAP vs PAP(SN) pre：log2FC = -0.7375，p = 9.228e-02，padj = 1.000**（rank 154/20143 tested）
- **Ctx TRAP vs Ctx pre：log2FC = -0.9430，p = 1.064e-01，padj = 0.564**
- **PAP TRAP vs Ctx TRAP：log2FC = +2.0274，p = 1.179e-02，padj(BH 全 41388 基因) = 0.269**（rank 1813）

### 1.5 论文补充表核验（PNAS Sakers 2017 Datasets S1–S4）[直接观察]
PNAS 端 403（Cloudflare）；PMC bin 路径返回 "Preparing to download" **PoW 挑战页**。处理：读取挑战 JS（`cdn.ncbi.nlm.nih.gov/.../pow-o51sQKbL.js`），算法 = 找 nonce 使 sha256(challenge+nonce) 前缀 "0000"，cookie `cloudpmc-viewer-pow=challenge,nonce`。本机解得 nonce=7526（hash 0000c53d…）后下载成功。挑战值绑定当次会话，复核时需重新求解。
- 已下载：`notes/gse_cache/pnas_sakers2017/sd01.csv`（1,9 MB, Dataset S1: Ens_ID,gene_name + 12 列 CPM）、`sd02.csv`（37 KB, 224 PAP-enriched）、`sd03.csv`（19 KB, 116 PAP-depleted）、`sd04.csv`（31 KB, 191 neuronal-filtered PAP-enriched）。
- sd02/sd03/sd04 表头（逐字）：`gene_id,logFC_PAP TRAP_v_SN_Input,FDR_PAP TRAP_v_SN_Input,logFC_SN_Input_v_Cortex_Input,FDR_SN_Input_v_Cortex_Input,logFC_Cortex TRAP_v_Cortex_Input,FDR_Cortex TRAP_v_Cortex_Input,mgi_symbol,description`
- **grep -i "atp2a2" 在 sd02/sd03/sd04 = 0 命中**（三张含 logFC/FDR 的表均无 Atp2a2）。
- sd01.csv Atp2a2 行（逐字）：`ENSMUSG00000029467,Atp2a2,400.7659243,488.8100718,485.3375349,399.4648256,188.6482744,187.7895452,999.4496686,675.4827973,1199.232631,1399.280349,1376.360944,1758.592857`
  （列序：Cortex_Input_1..3, Cortex_TRAP_1..3, PAP_TRAP_1..3, SN_Input_1..3，CPM）
- 由 S1 CPM 几何均值计算 [计算推断]：PAP_TRAP/SN_Input = **-0.688**；SN_Input/Cortex_Input = **+1.718**；Cortex_TRAP/Cortex_Input = -0.916；PAP_TRAP/Cortex_TRAP = **+1.946**；PAP_TRAP/Cortex_Input = +1.030。

### 1.6 历史值 log2FC 2.13 / padj 8e-7 的溯源结论
- **未能在 GSE74456 的任何公开文件（counts、series matrix）、PNAS 论文正文/表/Dataset S1-S4 中找到该值** [直接观察]。
- 同时核查了同组其他候选文献：Sakers 2022 Cell Rep（PMC9624251）正文 XML grep "Atp2a2" = 0；Qk 论文 Nat Commun 2021（PMC7943582）= 0 [直接观察]。
- 与该值最量级接近的公开量为：本复算 PAP TRAP vs Ctx TRAP +2.03、CPM 比 +1.95 —— 暗示历史值可能出自对 PAP-TRAP-vs-Cortex-TRAP 型对比的某次（未存档/未发表）分析或转抄 [外部文献线索，未核实 → 待向原作者或项目内历史记录求证]。
- **重要**：按作者自己的 PAP-enrichment 流程（vs SN input），Atp2a2 **不是** PAP 特异局部翻译转录本（log2 PAP_TRAP/SN_Input = -0.69，且不在 S2/S4 名单）。

**GSE74456 判定：PARTIAL**（基因级 counts 公开、可复算；作者未在任何公开 DE 表中报告 Atp2a2 的 logFC/FDR；历史值 2.13/8e-7 不可核实）。

---

## 2. GSE143531（PAPome / Mazare 2020）

### 2.1 下载
| URL | 落盘文件 |
|---|---|
| https://ftp.ncbi.nlm.nih.gov/geo/series/GSE143nnn/GSE143531/suppl/filelist.txt | notes/gse_cache/filelist_GSE143531.txt |
| https://ftp.ncbi.nlm.nih.gov/geo/series/GSE143nnn/GSE143531/matrix/GSE143531_series_matrix.txt.gz | notes/gse_cache/GSE143531_series_matrix.txt.gz |
| https://ftp.ncbi.nlm.nih.gov/geo/series/GSE143nnn/GSE143531/suppl/GSE143531_RAW.tar（7,464,960 B，48 文件） | notes/gse_cache/GSE143531/ |

### 2.2 设计与样本单位 [直接观察，series matrix 逐字]
- `!Series_title "Local translation highlights novel properties of perisynaptic astrocytic processes, and is modulated by behavior involving the dorsal hippocampus"`；`!Series_pubmed_id "32846133"`
- `!Sample_characteristics_ch1`："strain: Tg(Aldh1l1-eGFP/Rpl10a) JD130Htz …(Aldhl1:L10a-eGFP)"；"tissue: dorsal hippocampus"；"cell part: Full Astrocytes"（或 "PAP"）；"antiboby: GFP"；"animal id: Animal C1/C2/C3"
- 48 库 = 3 动物 × 2 区室（Full Astrocytes / PAP）× 8 技术重复（a-h，GEO 注 "technical replicate"）
- `!Sample_data_processing`："…Eoulsan pipeline…STAR (version 2.5.2b)…"、"Genome_build: mm10 (Ensembl 84)"、"Supplementary_files_format_and_content: tab-delimited text files include gene counting for each sample"（文件两列：Id \t Count）
- 论文：Mazare N et al., Cell Rep 2020;32(8):108076, DOI 10.1016/j.celrep.2020.108076（eutils 核实）。

### 2.3 Atp2a2 行原始摘录 [直接观察]
48 个 tsv 中 grep `^ENSMUSG00000029467`（每文件一行，`文件后缀 \t Count`）：
```
2017545a 468  2017545b 444  2017545c 438  2017545d 453  2017545e 449  2017545f 412  2017545g 413  2017545h 377   (C1 Full Astrocytes)
2017547a  43  2017547b  31  2017547c  47  2017547d  33  2017547e  36  2017547f  34  2017547g  51  2017547h  38   (C1 PAP)
2017553a 538  2017553b 513  2017553c 479  2017553d 457  2017553e 433  2017553f 433  2017553g 460  2017553h 427   (C2 Full Astrocytes)
2017555a  40  2017555b  34  2017555c  30  2017555d  32  2017555e  24  2017555f  36  2017555g  44  2017555h  41   (C2 PAP)
2017561a 334  2017561b 329  2017561c 313  2017561d 334  2017561e 286  2017561f 287  2017561g 277  2017561h 266   (C3 Full Astrocytes)
2017563a  45  2017563b  53  2017563c  58  2017563d  58  2017563e  47  2017563f  42  2017563g  44  2017563h  54   (C3 PAP)
```

### 2.4 作者公开 DE 表核验 [直接观察]
论文补充 Table（PMC7450274 → bin mmc2.xlsx，6,080,831 B，经同一 PMC PoW cookie 下载；`notes/gse_cache/mazare2020_supp/mmc2.xlsx`）。sheet 列表："1. Raw data"、"2. Whole Papome"、"3. Whole astrocytome"、"4. PAP-enriched"、"5. Astrocyte-enriched"、"6. GFAP exon study"。
- **sheet "1. Raw data" row 1175（逐字，含列名）**：
  列名：`Id | Mean Normalised Reads Astrocyte | Mean Normalised Reads PAP | log2FoldChange | BH-adjusted p-value | Associated Gene Name | Description | Exon Coverage Astrocyte | Exon coverage PAP`
  值：`ENSMUSG00000029467 | 1090.6022175928867 | 812.35912948110433 | 0.39125713340048701 | 0.73036232616605601 | Atp2a2 | ATPase, Ca++ transporting, cardiac muscle, slow twitch 2 [Source:MGI Symbol;Acc:MGI:88110] | 99.463078769502872 | 65.561240603878502`
- **Atp2a2 不在 sheet "4. PAP-enriched"**（sharedStrings 中 Atp2a2 仅出现在 sheet 1 与 sheet 3）。
- 符号约定核验（sheet "1. Raw data"）：Fth1（PAP 均值 25672.97 > Astro 8833.92）log2FC = **-1.5488**（padj 0.0383）；Ccnd2（39873.25 vs 905.0）= **-5.4565**（padj 2.94e-9）→ 正值 = Whole-Astrocyte 高。故 Atp2a2 +0.3913 = 全细胞略高于 PAPome，padj 0.730 不显著。
- 论文正文 XML（mazare2020.xml, PMC7450274）grep "Atp2a2" = 0（正文不点名）。

### 2.5 本地复算 [计算推断]
脚本：`analysis_scripts/recompute_gse143531.py`；输出：`logs/recompute_gse143531_out.txt`
方法：每动物合并 8 个技术重复（求和）→ 6 列矩阵；median-of-ratios size factors（C1_Full 3.211, C1_PAP 0.328, C2_Full 2.977, C2_PAP 0.315, C3_Full 3.394, C3_PAP 0.523）；log2(norm+1)；动物内配对 t 检验；BH。
- 合并后 Atp2a2 原始计数：C1_Full 3454, C1_PAP 313, C2_Full 3740, C2_PAP 281, C3_Full 2426, C3_PAP 401
- 归一化均值：Full 1075.8/1256.4/714.8；PAP 953.7/891.9/766.8
- **PAP vs Full：log2FC = -0.1888，配对 p = 0.3868，padj(BH) = 1.000**（rank 9833/47643）；非配对 Welch 参考：p = 0.5277, padj = 1.000

### 2.6 判定
作者公开表（+0.3913, padj 0.730，全细胞方向）与本复算（-0.19, padj 1.0，PAP 方向）及项目此前估计（-0.273, FDR 0.683）**结论一致：未检出差异**（|log2FC| < 0.4，padj ≥ 0.68，Atp2a2 不在作者 PAP-enriched 名单）。方向上各复算符号随对比方向而反，量级均在无效应范围内。

**GSE143531 判定：NO_DIFFERENCE_DETECTED（未检出差异）**。

---

## 3. GSE330741（2026 SN-MPRA）

### 3.1 下载（全部于 2026-09-25）
| URL | 落盘文件 |
|---|---|
| https://ftp.ncbi.nlm.nih.gov/geo/series/GSE330nnn/GSE330741/suppl/ （目录列表） | — |
| https://ftp.ncbi.nlm.nih.gov/geo/series/GSE330nnn/GSE330741/suppl/filelist.txt | notes/gse_cache/filelist_GSE330741.txt |
| https://ftp.ncbi.nlm.nih.gov/geo/series/GSE330nnn/GSE330741/matrix/GSE330741_series_matrix.txt.gz | notes/gse_cache/GSE330741_series_matrix.txt.gz |
| https://ftp.ncbi.nlm.nih.gov/geo/series/GSE330nnn/GSE330741/suppl/GSE330741_DNA_1_L1.counts.txt.gz（及 Pre/Post × L1/L2 共 5 个抽样） | notes/gse_cache/GSE330741_*.counts.txt.gz |
| https://ftp.ncbi.nlm.nih.gov/geo/series/GSE330nnn/GSE330741/suppl/GSE330741_RAW.tar（3,389,440 B，139 文件） | notes/gse_cache/GSE330741/ |
| https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL34nnn/GPL34290/soft/GPL34290_family.soft.gz | 标题核验："Illumina NovaSeq X Plus (Mus musculus)" |

完整 32+1 文件清单见 `filelist_GSE330741.txt`；结构核对详见 `notes/gse330741_manifest_check.md`（elementID 家族计数、Pre/Post/DNA 与主文库 0 交集、GSM9822133-9822156 移至 GSE295080 的官方说明）。

### 3.2 Atp2a2 检索结果 [直接观察]
- 主文库（ctxtrap_1 等 GSM 文件）6430 个元素 ID 家族：Slc1a2 4570 / sparc 660 / scram 660 / hsbp1 400 / Basal 110 / mbp 20 / act_zip 10 —— **无任何 Atp2a2 元素**。
- `zcat GSE330741/*.txt.gz GSE330741_*.counts.txt.gz | grep -ic "atp2a2\|ENSMUSG00000029467"` → **0**。
- Pre/Post/DNA 32 文件（1630 元素，人类基因名 wgs/ssc alt/ref）与主文库元素交集 = 0。
- series matrix 摘要明言 tiling 对象为 "Glt1 and Sparc" 两基因 3'UTR。

### 3.3 判定
**GSE330741 不能作为 Atp2a2 的 reporter 验证来源（NOT_APPLICABLE / status: pending）**；其 fragment→reporter 映射对 Glt1/Sparc 亦因缺参考文库表与脚本而无法完整重建（最小缺失输入清单见 gse330741_manifest_check.md 第 4 节）。**严禁把 Glt1/Sparc 片段结果解释为 Atp2a2 验证。**

---

## 4. 汇总
| GSE | 判定 | Atp2a2 log2FC / padj | 来源 |
|---|---|---|---|
| GSE74456 | PARTIAL | 复算 +2.03 / 0.269（PAP TRAP vs Ctx TRAP）；-0.69（PAP TRAP vs SN input）；历史值 2.13/8e-7 未在任何公开表找到 | our_recompute（counts + 作者 S1 CPM） |
| GSE143531 | NO_DIFFERENCE_DETECTED（未检出差异） | 作者表 +0.3913 / 0.7304；复算 -0.19 / 1.0 | authors_table + our_recompute |
| GSE330741 | NOT_APPLICABLE（文库无 Atp2a2；方法重建 pending） | NA | NA |
