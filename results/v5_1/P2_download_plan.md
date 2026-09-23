# P2-05 下载与比对方案（⚠ 用户决策项）

## 现状
- 本机无 SRA Toolkit（fastq-dump/fasterq-dump/prefetch）、无 STAR/hisat2/minimap2/salmon/kallisto。
- GSE238125 = 12 runs PE（SRR25403062–73），总 80.8 G bases ≈ 40 GB FASTQ（gzip 后），单 run 3–4 GB。
- 下载走 NCBI 单流限速约 0.02 MB/s（下午实测），并行 Range 分块对 SRA 对象是否有效待验证——但 SRA 原始文件（.sra）无法用 Range 分块工具直接下；需走 https://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR254SRXXX/SRR25403062/SRR25403062.sra 每文件整体下载，或用 AWS/EBI 镜像。
- 更关键的是比对：装 STAR（需 ~32GB RAM 对小鼠全基因组 + 32GB 索引磁盘）+ 12 runs 比对在笔记本上约需 6–12 小时 CPU；且 APA 定量还需要 DaPars2（Python 版 DaPars2_py 有 PyPI 包，但需要 BAM + Gencode 注释）。

## 候选方案

### 方案 A：FASTQ 全流程（传统）
- 用户手动下载 12 个 .sra 文件（或更快q直接 FASTQ：ENA 的 fastq_1/2.gz 直接可下，EBI 无 NCBI 式单流限速）→ 装 sra-tools（conda/winget）转 fastq → 装 STAR 比对 → featureCounts/salmon 定量 → DaPars2_py。
- 优点：最标准，后续所有分析可复用。缺点：用户要下 40GB；本机装比对器 + 索引构建 + 比对总耗时可能 1–2 天（笔记本 RAM 未知）。

### 方案 B：改用 Salmon/selective-alignment 或 kallisto 伪比对
- 免去全基因组比对，用 transcriptome 索引（~4GB），速度快 5–10 倍。但对 APA 分析不友好：APA 需要**基因组坐标层面的 3'UTR 覆盖**（DaPars2 输入是 BAM），transcriptome 伪比对不能直接产出 BAM。需换工具链（如 kallisto+bustools 的 3'UTR 定制或 tail-tools 类），风险高。

### 方案 C：跳过比对——用 SRA 的比对产物（推荐评估）
- SRA 高级格式/Elite 区域有 .sra 原始文件，但 NCBI 已下线大部分预比对 BAM。
- **检查点：GSE238125 的 GEO 补充只有 gene counts（408KB csv），无 per-transcript/PAS 信息 → 无法做 APA，必须回到 FASTQ。**

### 方案 D：降维打击——只做 M3 层面可复用的部分，M2 事件层改用现成 APA 目录
- 用 APAatlas/PolyASite 的已有 stroke-relevant 数据？——查新已确认无 stroke astrocyte APA 目录，此路不通。

## 推荐：方案 A（分两步走，减轻用户负担）
1. **第一步（本轮）**：先只下载 **Sham rep1/rep2 + day1 rep1/rep2 四个 run（~13.9 GB FASTQ）**做试跑闭环（比对 1 run → 定量 → DaPars2 单时间点），确认工具链在笔记本上跑得动、输出正确，再决定是否全量 12 runs。
2. **第二步**：试跑通过后，剩余 8 runs 由用户手动下载（附 B 链接表）或我按同法分批拉。

## 用户需要做的事（方案 A 第一步）
- 无需下载软件：sra-tools/STAR 我可通过 conda（若有）或 winget/pip 装；**若本机 RAM < 32GB，STAR 改用 hisat2（索引 ~8GB，内存需求低）或 minimap2 -ax sr（索引 ~4GB，速度最快）**——先测 RAM 再定。
- 若走 ENA 直链 FASTQ（免 sra-tools），四个 run 的链接我会生成到附 B。

## 请用户确认
1. 批准方案 A 两步走（先 4 runs 试跑闭环）？
2. 本机可用磁盘空间（FASTQ 40GB + 比对临时文件 ~50GB + 索引）是否充足？数据盘当前可用 ____ GB。
