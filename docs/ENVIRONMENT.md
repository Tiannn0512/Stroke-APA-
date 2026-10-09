# Stroke&APA 工作环境清单

> 生成：2026-09-16 23:45（版本均为当日实查，非记忆值）
> 课题：卒中→APA 调控子重塑→星形胶质细胞 mRNA 3'UTR 同工型→足突定位潜能（pipeline v5.1，纯干实验）
> 纪律：M2 主分析全程纯 Python（不新增 R 依赖）；R/salmon 仅存在于前期试点资产中

## 1. 主机与操作系统

| 项 | 值 |
|---|---|
| 宿主 OS | Windows 10.0.26200 x64（Windows_NT），主机名 Taylor |
| 内存 | 31.3 GB 物理（WSL2 可见 23 GB） |
| CPU | 12 逻辑核 |
| Linux 子系统 | WSL2 Ubuntu 24.04.2 LTS（发行版用户 taylor，无免密 sudo） |
| 默认 Shell | PowerShell（Windows 侧）；bash（WSL 侧） |

## 2. 工作目录（三个根，职责分离）

| 根目录 | 职责 | 关键内容 |
|---|---|---|
| `D:\stroke_apa\` | **主工作区**（代码/文档/小结果） | `TODO_pipeline_v5.1.md`（唯一主清单）、`STROKE_APA_PIPELINE_v5.1.md`、`results\`（P1 冻结参数/DE 表/regulator 目录/data_inventory.tsv）、`scripts\`（全部流程脚本）、`memory\`（日记） |
| `D:\stroke_apa\data\` | 小体量已下载数据 | `GSE174574\`（scRAW.tar 263,086,080 B + 解包 6 样本 10x 三件套 + h5ad）、`GSE286075\`（1,802,240 B，6×ReadsPerGene.out.tab） |
| `D:\stroke_apa_data\` | **大数据区**（用户指定，2026-09-16 生效） | `fastq\`（GSE238125 FASTQ）、`bam\`、`index\`、`tools\`、`ena_filereport.tsv`（12 runs 元数据）、`fastq_md5_all.tsv`（校验基准） |
| WSL `~/reference_v2/` | **⚠ 关键生产资产（2026-09-17 起，禁止删除！）** | STAR 索引、GRCm38 genome.fa、GENCODE vM25 GTF、DaPars2 3'UTR 注释、chromosomes.txt——**比对链与 DaPars2 的运行依赖，删除=全链瘫痪**（前身 `~/stroke_APA_pilot/` 于 09-17 16:4x 被误删，教训见 memory/2026-09-17.md） |
| WSL `~/stroke_APA_work/` | 比对工作区（ext4，性能关键） | bam/bedgraph/depth/qc；产物定期备份至 `D:\stroke_apa_data\ext4_backup_*\` |
| WSL `~/miniconda3/`、`~/DaPars2/` | ⚠ 关键依赖 | dapars2/meme 等 conda env 与 DaPars2 源码 |

## 3. Windows 侧软件与 Python 包

系统级：

| 软件 | 版本 | 路径 |
|---|---|---|
| Python | 3.13.12 | `D:\Auto Claw\AutoClaw\resources\python\python.exe` |
| Git | 2.55.0.windows.5 | `C:\Users\Taylor\Desktop\lxPC\Git\cmd\git.exe` |
| curl | 8.18.0 (Schannel) | `C:\WINDOWS\system32\curl.exe` |
| conda / sra-tools / 比对器 | **缺失**（本机不装；比对走 WSL） | — |

Python 包（P0/P1 实际使用，pip 实测版本）：

| 包 | 版本 | 用途 |
|---|---|---|
| scanpy | 1.12.4 | scRNA QC、聚类、UMAP |
| anndata | 0.13.3.post0 | h5ad 读写 |
| pydeseq2 | 0.5.4 | 伪批量 DE（P1 主证据层） |
| leidenalg | 0.12.0 | Leiden 聚类 |
| igraph | 1.0.0 | 图后端 |
| numpy | 2.5.3 | 数值基础 |
| pandas | 3.0.5 | 表格处理 |
| scipy | 1.18.1 | 稀疏矩阵/统计 |
| h5py | 3.16.0 | HDF5 底层 |
| matplotlib | 3.11.2 | 出图 |

注：scrublet 在本环境（py3.13）安装失败 → P1-02 按冻结文件降级为 MAD±3 双联体回退，偏差已登记。

## 4. WSL 侧软件（比对 + APA 全链）

系统路径（`/usr/bin`，Ubuntu 源）：

| 软件 | 版本 | 用途 |
|---|---|---|
| STAR（STAR-avx2） | 2.7.11b | 比对（试点实测 95.75% 唯一比对，~8 min/样本，8 线程） |
| fastp | 0.23.4 | 接头裁剪（PE，--detect_adapter_for_pe） |
| bedtools | v2.31.1 | genomecov → bedgraph（**用 -split**，reanalysis_v3 纠错语义） |
| samtools | 1.19.2 | sort/index/depth（深度口径 `-F 0x904` 主比对） |
| git / gcc / make | 系统版 | 编译后备 |

Miniconda3（`~/miniconda3/`）四个 env：

| env | Python | 关键工具 | 状态 |
|---|---|---|---|
| bioinfo | 3.10.21 | 通用生物信息 | 前期试点用 |
| dapars2 | 3.8.20 | **DaPars2 三脚本**（DaPars2_Multi_Sample_Multi_Chr.py 等，源码 `~/DaPars2/src/`） | P2-06 主力 |
| qapa | 3.9.19 | salmon 1.10.3、R 4.3.3（Rscript）、QAPA | 前期试点用；v5.1 不新增依赖 |
| salmon | —（salmon 独立） | salmon 1.10.3 | 前期试点用 |

## 5. 参考数据资产（WSL，GSE238125 直接复用，零下载）

| 资产 | 位置（WSL 路径） | 规格 |
|---|---|---|
| 基因组 FASTA（GRCm38/mm10） | `~/stroke_APA_pilot/reference/genome/genome.fa`（+.fai） | 2.6 GB，primary assembly |
| 基因注释 | `.../gencode.vM25.annotation.gtf`、`genes.gtf` | GENCODE vM25 |
| STAR 索引 | `.../STAR_index/` | 24 GB（23GB RAM 约束下建成）；另有 `STAR_index_lowmem`（8KB，占位未完成，**勿用**） |
| RSEM 索引 | `~/stroke_APA_pilot/reference/RSEM_mm10_M25_bt2/` | bowtie2 版（试点 QAPA 路线用） |
| 3'UTR 注释（DaPars2 输入） | `~/stroke_APA_pilot/dapars2/annotation/gencode_M25_3UTR_for_DaPars2.bed` | 现成，P2-06 直接引用 |
| PolyASite PAS 床文件 | `~/stroke_APA_pilot/reference/PAS/PolyASite_mouse_mm10_chr.bed` | 30,100 行（chr 坐标版 + sorted 版） |
| 染色体列表 | `~/stroke_APA_pilot/scripts_package/configs/chromosomes.txt` | DaPars2 多染色体循环用 |

## 6. 数据集与下载状态（截至 2026-09-16 23:45）

| 数据集 | 内容 | 状态 |
|---|---|---|
| GSE174574 | MCAO scRAW（10x 3'，6 样本） | ✅ 已下载解包校验（E 盘），P1 用毕 |
| GSE286075 | endfoot RiboTag gene counts ×6 | ✅ 已下载（E 盘），P1 用毕 |
| GSE238125 | 皮层 astrocyte bulk RNA-seq，12 runs 双端 FASTQ，**47.66 GB** | 🔄 下载中（D 盘；试点 4 runs = Sham×2 + d1×2，sham1 完整 4.59GB，sham1_2 ~98%，速率 ~0.4 MB/s，后台续传） |
| GSE225110 / GSE263986 / GSE74456 / GSE146935 / GSE330741 | P2-07/P3/P4 参考集 | ⏳ 未触发 |

## 7. 关键路径速查

```
主清单      D:\stroke_apa\TODO_pipeline_v5.1.md
结果登记表  D:\stroke_apa\results\data_inventory.tsv
比对脚本    D:\stroke_apa\scripts\p2_gse238125_align.sh   （WSL 内执行）
DaPars2    D:\stroke_apa\scripts\p2_dapars2_run.sh        （WSL 内执行）
下载器      D:\stroke_apa\scripts\p2_dl_pilot.py           （Windows python，后台）
MD5 校验    D:\stroke_apa\scripts\p2_md5_check.py
FASTQ      D:\stroke_apa_data\fastq\<sample>_<R>.fastq.gz
深度表      D:\stroke_apa_data\sequencing_depth.tsv        （比对产物）
样本映射    D:\stroke_apa\scripts\p2_sample_map.tsv         （12 样本↔SRR）
```

## 8. 前期试点（~/stroke_APA_pilot）一句话背景

9月9–15日已用 GSE225110 的 6 个 run 跑通 QAPA+DaPars2 全链并完成 reanalysis_v3 纠错：事件级重建（2,231 eligible tandem 事件）、-split 修正、置换检验后聚合信号不超零假设、旧候选全部降级、留 7 个探索级候选。该轮结论是 v5.1 方法学红线（PAS 背景匹配随机化、表述纪律）的来源；本轮 GSE238125 是独立的新数据层，不依赖其结论，只复用其参考资产与教训。
