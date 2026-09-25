# 独立公共数据集筛查日志 — Atp2a2 3'UTR/APA 定点核查用卒中队列

- 项目: 小鼠卒中-星形胶质细胞-APA 复核 (发现队列 GSE238125, 预测 PAS chr5:122456498)
- 检索与核对执行日期: **2026-09-25** (所有检索式与记录核对均于当日完成)
- 输出: `../public_dataset_suitability.tsv` (18 个数据集条目)

## 0. 访问路径说明

本会话中 `www.ncbi.nlm.nih.gov` 的 GEO HTML 页面(WebFetch 与 curl 两种方式)均返回 reCAPTCHA 拦截页。
因此所有记录改经以下两条可复现路径核实:

1. NCBI E-utilities: `esearch`/`esummary`(db=gds, db=pubmed), `efetch`(db=pubmed, rettype=abstract)。
2. GEO FTP 的 SOFT family 文件: `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE{前3位}nnn/{ACCESSION}/soft/{ACCESSION}_family.soft.gz`。

凡写入 TSV 的 accession, 均于 2026-09-25 打开过其 SOFT 记录(eutils 或 FTP), 逐字段核对物种/脑区/模型/时间点/细胞方法/重复数/读长与比对策略。本日志中不收录任何未能打开核实的 accession。

## 1. 检索式清单 (全部执行于 2026-09-25)

| # | 引擎/接口 | 检索式 | URL | 命中与去向 |
|---|---|---|---|---|
| Q1 | NCBI esearch db=gds | `astrocyte[All Fields] AND MCAO[All Fields] AND gse[ETYP]` | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=astrocyte%5BAll+Fields%5D+AND+MCAO%5BAll+Fields%5D+AND+gse%5BETYP%5D | 12 hits; 核对出 GSE244576(snRNA,药物), GSE268852(外泌体circRNA), GSE247474(scRNA 12h), GSE210674(ACSA-2-array), GSE197711(体外球状体), GSE120565(全脑array), GSE28731(Hsp72 array), GSE35338(FACS星形胶质细胞array), **GSE103783/GSE103782/GSE103781(Cx43-RiboTag tMCAO)** |
| Q2 | NCBI esearch db=gds | `astrocyte[All Fields] AND ischemic stroke[All Fields] AND RNA-Seq[All Fields] AND mus musculus[ORGN] AND gse[ETYP]` | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=astrocyte%5BAll+Fields%5D+AND+ischemic+stroke%5BAll+Fields%5D+AND+RNA-Seq%5BAll+Fields%5D+AND+mus+musculus%5BORGN%5D+AND+gse%5BETYP%5D | 6 hits; **GSE225110(星形胶质细胞/小胶质细胞 RiboTag, 4h+3d)**, GSE276260(原代星形胶质细胞 H2O2 体外), GSE279665(sc/sn 流程), GSE253799(Matn3 KO 皮层), GSE203050(Treg 共培养), GSE167593(scRNA 缺血vs出血) |
| Q3 | NCBI esearch db=gds | `ACSA2[All Fields] AND gse[ETYP]` | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&retmax=40&term=ACSA2%5BAll+Fields%5D+AND+gse%5BETYP%5D | 28 hits; 逐一核对标题: 全部为 MPTP 帕金森(GSE244652/GSE191131)、EAE(GSE260982)、ALS(GSE275841)、肿瘤、饮食(GSE205313)、癫痫、弓形虫等; **无卒中队列** |
| Q4 | NCBI esearch db=gds | `Aldh1l1[All Fields] AND stroke[All Fields] AND gse[ETYP]` | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&retmax=40&term=Aldh1l1%5BAll+Fields%5D+AND+stroke%5BAll+Fields%5D+AND+gse%5BETYP%5D | 4 hits; GSE277064(脊髓 RiboTag), **GSE286075(endfoot RiboTag 缺血)**, **GSE225110**, **GSE35338** |
| Q5 | NCBI esearch db=gds | `neurotoxic[All Fields] AND reactive astrocytes[All Fields] AND gse[ETYP]` | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&retmax=30&term=neurotoxic%5BAll+Fields%5D+AND+reactive+astrocytes%5BAll+Fields%5D+AND+gse%5BETYP%5D | 14 hits; 均非 Liddelow 2017 pMCAO 集散数据(GBA1, 线粒体复合物I, 眼压, 斑块微环境, C3, 癫痫, ALS 等) |
| Q6 | NCBI esearch db=gds | `RiboTag[All Fields] AND gse[ETYP]` | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&retmax=20&term=RiboTag%5BAll+Fields%5D+AND+gse%5BETYP%5D | 20 hits; GSE336013/GSE305388(MiRA-AAV 诱导,无卒中), GSE344979(生物材料注射), GSE344344(PV神经元), 其余为神经元/心肌/视网膜 |
| Q7 | NCBI esearch db=gds | `"3'READS"[All Fields] AND gse[ETYP]` | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&retmax=30&term=%223%27READS%22%5BAll+Fields%5D+AND+gse%5BETYP%5D | 30 hits; 以酵母/细胞系为主; GSE150607 为 3'READS+ 但仅 C2C12/NIH3T3(非脑); **无卒中脑 3'READS** |
| Q8 | NCBI esearch db=gds | `"alternative polyadenylation"[All Fields] AND brain[All Fields] AND mus musculus[ORGN] AND gse[ETYP]` | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&retmax=40&term=%22alternative+polyadenylation%22%5BAll+Fields%5D+AND+brain%5BAll+Fields%5D+AND+mus+musculus%5BORGN%5D+AND+gse%5BETYP%5D | 22 hits; **GSE94054(cTag-PAPERCLIP 小胶质细胞 LPS)**, **GSE108480(小脑神经元 PAPERCLIP)**, **GSE142683(Nudt21/CFIm25 海马)**, **GSE310389(酒精, QuantSeq 3')**, GSE347392/347390(催产素受体), GSE185662(U1 IPA), GSE119073(Elavl3), GSE169023(HuD KO), GSE86227/86224/86043(hnRNPA2B1) |
| Q9 | NCBI esearch db=gds | `"direct RNA"[All Fields] AND nanopore[All Fields] AND mus musculus[ORGN] AND gse[ETYP]` | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&retmax=25&term=%22direct+RNA%22%5BAll+Fields%5D+AND+nanopore%5BAll+Fields%5D+AND+mus+musculus%5BORGN%5D+AND+gse%5BETYP%5D | 25 hits; **GSE171312(海马 Nanopore direct RNA, APP NL-G-F vs WT)** 入选; 其余为巨噬细胞/mESC/非脑 |
| Q10 | NCBI esearch db=gds | `("Oxford Nanopore" OR "Iso-Seq") AND brain AND mus musculus[ORGN] AND gse[ETYP]` | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&retmax=25&term=%28%22Oxford+Nanopore%22%5BAll+Fields%5D+OR+%22Iso-Seq%22%5BAll+Fields%5D%29+AND+brain%5BAll+Fields%5D+AND+mus+musculus%5BORGN%5D+AND+gse%5BETYP%5D | 25 hits; GSE315051(成年皮层长+短读段, TRP通道), GSE271020(circRNA/inosine), 垂体/耳蜗等; **无卒中长读段队列** |
| Q11 | NCBI esearch db=gds | `astrocyte AND "permanent middle cerebral" AND gse[ETYP]` | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&retmax=30&term=astrocyte%5BAll+Fields%5D+AND+%22permanent+middle+cerebral%22%5BAll+Fields%5D+AND+gse%5BETYP%5D | 1 hit: GSE234052(周细胞, 非星形胶质细胞) |
| Q12 | NCBI esearch db=gds | `Liddelow[All Fields]` | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&retmax=30&term=Liddelow%5BAll+Fields%5D | 30 hits(Barres 实验室遗留系列为主: GSE148610/11/12 神经炎症星形胶质细胞亚型, GSE109352-4 IL-33, GSE52564 皮层细胞图谱); **未找到 2017 Nature pMCAO RiboTag 系列** |
| Q13 | Web 检索(工具) | `GSE238125`, `GSE74456`, `GSE143531`, `GSE330741` 身份确认; `Liddelow 2017 "Neurotoxic reactive astrocytes" GEO accession` | (GEO 镜像与 PMC 线索) | 四个已知 accession 身份确认; Liddelow 检索提示 GSE87497, 经 SOFT 核实为间充质干细胞数据(排除), **Liddelow 2017 数据登录号未能核实, 该资源未入表** |
| Q14 | GEO FTP | SOFT family 下载核对 | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE{p}nnn/{ACC}/soft/{ACC}_family.soft.gz` | GSE238125, GSE74456, GSE143531, GSE330741, GSE103782, GSE103781, GSE225110, GSE35338, GSE210674, GSE286075, GSE137482, GSE276260, GSE247474, GSE94054, GSE108480, GSE142683, GSE310389, GSE150607, GSE87497, GSE171312 |
| Q15 | NCBI efetch/esummary db=pubmed | PMIDs: 38713438, 28439016, 32846133, 30585358, 22553043, 32553170, 36880055, 39796165, 42049021, 37067534 | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id=... | 引文核对: Wang S Mol Neurobiol 2024; Sakers PNAS 2017; Mazaré Cell Rep 2020; Rakers Glia 2019; Zamanian J Neurosci 2012; Androvic Cell Rep 2020; Yamaguchi IBRO Neurosci Rep 2023; Shim Int J Mol Sci 2025; Lee JH Cell Metab 2026; Hernandez Glia 2023 |

## 2. 各检索方向结论

### 方向1: 卒中 + FACS/ACSA2/急性分离星形胶质细胞 bulk RNA-seq
除发现队列 GSE238125 外, GEO 中未检索到第二个"卒中 + FACS/ACSA2 星形胶质细胞 + bulk RNA-seq"数据集。ACSA2 检索式(Q3)的 28 个命中全部为非卒中疾病模型。最接近的是 GSE35338(Zamanian 2012, Aldh1l1-eGFP FACS 星形胶质细胞, MCAO 1d/3d/7d), 但其为 Affymetrix 430 2.0 微阵列, 只能做基因层面复核; GSE210674(ACSA-2+, Clariom S array)因无 sham、n=2、运动混杂判 INSUFFICIENT。

### 方向2: Aldh1l1-TRAP / RiboTag 卒中翻译组
本方向收获最大: **GSE225110**(Hernandez/Buckwalter, Glia 2023; 星形胶质细胞 IP, sham vs 卒中 4h 与 3d, 双性别, 132 样本, RSEM 定量, 原始数据开放)为最佳独立基因层面队列; GSE103782/103781(Rakers 2019, Cx43-RiboTag tMCAO, n=3+3)与 GSE286075(星形胶质细胞 endfoot RiboTag, tMCAO 2h+6h, 配对 n=3)为补充。均为 RiboTag IP RNA-seq, 按约定判 GENE_LEVEL_ONLY, 无异构体分辨力。

### 方向3: 卒中星形胶质细胞 scRNA-seq
GSE247474(12h, 2 个 pooled 文库)、GSE167593(缺血 vs 出血, 3 文库)、GSE279665(病程 sc/sn)均为 10x 类 3' 计数, 判 NOT_APPLICABLE。除非出现全长读段 scRNA, 此类数据一律不得用于 Atp2a2 APA 验证。

### 方向4: 小鼠脑 poly(A) 位点 / 3'端测序资源
存在可用的非卒中 PAS 佐证资源: GSE94054(cTag-PAPERCLIP, 小胶质细胞 LPS 激活, 内源 PAS 图谱, 与卒中神经炎症最接近)、GSE108480(小脑神经元 PAPERCLIP 基线)、GSE142683(Nudt21/CFIm25 海马 APA)、GSE171312(海马 Nanopore direct RNA, 唯一长读段全转录本脑资源, 可读出 Atp2a2 3'UTR 异构体基线)。GSE310389(QuantSeq 3', 酒精模型)仅作方法参考。**不存在卒中脑的 3'端/PAS 或长读段公共数据。**

### 已知条目复核
GSE238125(发现队列, sham/d1/d3/d7/d21/d60 各 n=2, FACS 星形胶质细胞, mm10, 基因计数矩阵+RAW.tar)、GSE74456(PAPTRAP, 健康P21)、GSE143531(PAPome, 恐惧条件化)、GSE330741(SN-MPRA, Glt1/Sparc 3'UTR tile; 2026-09-08 起 24 个 GSM 移至 GSE295080)均经 SOFT 复核, 结论与既有评估一致, 详见表。

## 3. 总结论 (3-5 句)

1. **最有希望的独立队列是 GSE225110**(Hernandez VG et al., Glia 2023; PMID 37067534): Aldh1l1-RiboTag 星形胶质细胞 IP、sham vs 卒中 4h/3d、双性别、每组 3-12 只、原始 FASTQ 开放, 可对 Atp2a2 做独立基因层面(及 RSEM 转录本层面参考性)复核; 其次为 GSE103782/103781(tMCAO Cx43-RiboTag)与 GSE35338(MCAO d1/d3/d7 FACS 星形胶质细胞微阵列), 组织层面可加 GSE137482(pMCAO 3d 皮层 QuantSeq 3')。
2. **异构体层面存在明确缺口**: 检索范围内不存在"卒中 + 星形胶质细胞富集 + 3'端/PAS 测序或长读段"的公共数据集, 也没有卒中脑的任何 3'READS/PAPERCLIP/DRS/nanopore 数据; GSE238125 原始 FASTQ 是唯一能对 chr5:122456498 PAS 做覆盖位移观察的卒中数据, 但 n=2 且即发现队列本身, 不构成独立验证。
3. PAS 佐证可在非卒中资源中先行完成: 用 GSE94054(小胶质细胞炎症)、GSE108480(神经元基线)、GSE142683(CFIm25)与 GSE171312(direct RNA)确认 Atp2a2 该 PAS 在鼠脑的基线使用与注释; 卒中特异性变化只能诉诸新实验(FACS/ACSA2 分选 MCAO 皮层星形胶质细胞 + 3'READS+ 或 Nanopore direct RNA)。
4. 卒中星形胶质细胞 scRNA-seq(10x 3' 计数)按约定全部判 NOT_APPLICABLE, 不得充当异构体证据。
5. Liddelow et al. 2017(Nature 541:481)pMCAO/LPS 反应性星形胶质细胞 RiboTag 资源的 GEO 登录号在本次会话无法核实(GSE87497 已查证为无关数据), 为避免编造 accession 未纳入表; 若后续人工确认, 可作为方向2的又一基因层面独立队列。

---
检索与记录核对: 自动化会话, 2026-09-25; 工具: NCBI E-utilities, NCBI FTP(SOFT), Web 检索。所有 accession 均以 SOFT/eutils 记录为准。
