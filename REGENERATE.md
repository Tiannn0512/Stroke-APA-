# REGENERATE —— 数据再生手册（本地删除后，一切数据如何找回）

> 本仓库刻意**不**收录任何可用网络重新获取或可由脚本确定性重新生成的大文件（原始 FASTQ、BAM、参考基因组/索引、FIMO 原始 XML、null 序列池、IGV 切片）。
> 本文件逐项写明每类数据的重新获取命令与重建顺序。**验证锚点**：全部重建后应复现 `docs/reanalysis/实验流程.md` §0 的数字锚点。

## 0. 重建顺序总览

```
① 环境还原（envs/）→ ② 参考数据下载（GENCODE/PolyASite/phastCons）→ ③ 索引构建（STAR/salmon/slim bed）
→ ④ FASTQ 重下（ENA）→ ⑤ 比对（p2_gse238125_align.sh）→ ⑥ DaPars2/SAP（按 实验流程.md §5.2–5.3）
→ ⑦ 其余 GEO 数据按需重下
```

## 1. 环境还原

| env | 用途 | 还原 |
|---|---|---|
| `dapars2` | DaPars2 运行环境 | `conda env create -f envs/dapars2.yml`；**打过补丁的 DaPars2 源码已随仓库入库**：`tools/DaPars2_patched/`（src/*.py + `.bak_20260918` 为打补丁前原件，可 diff 还原补丁语义=空覆盖事件返回 NA） |
| `rlimma` | R 4.3.3 + limma 3.58.1 | `conda env create -f envs/rlimma.yml`；edgeR 4.0.16 另装：先 sed 修 Makeconf 编译器（详见 实验流程.md §2.2），`install.packages('edgeR', lib='~/Rlibs')` |
| `bioinfo` | pysam/pyBigWig 等 | `conda env create -f envs/bioinfo.yml` |
| `salmon` | salmon 1.10.3 | `conda env create -f envs/salmon.yml` |
| `meme` | FIMO 5.5.9（MEME suite） | `conda env create -f envs/meme.yml` |
| scRNA 环境 | scanpy/pydeseq2/scrublet | 若 yml 导出失败，按 `envs/README.md` 手工清单重建 |

## 2. 参考数据（全部公开、免费）

| 数据 | 获取命令 |
|---|---|
| GRCm38 genome.fa | `wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz` |
| GENCODE vM25 GTF / transcripts | `wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz` 及 `gencode.vM25.transcripts.fa.gz` |
| PolyASite 2.0 mm10 atlas（301,006 clusters） | 新址 **polyasite.unibas.ch**（旧 .ethz.ch 已废弃）→ 下载页 `/download/atlas/2.0/GRCm38.96/`，文件名 `atlas.clusters.2.0.GRCm38.96.bed.gz`（以站点实际列表为准） |
| mm10.60way.phastCons.bw（4.58GB） | `python scripts/reanalysis/p2_dl_phastcons.py`（分块并行下载脚本已入库） |
| ATtRACT motif 库 | **无需重下**：冻结副本已在 `input_links/frozen_rbp_motifs.meme` + `input_links/frozen_motif_family_mapping.tsv`（原站 403 时的 zavolanlab GitHub 快照） |
| QKI CLIP peaks（GSE147119，437 peaks） | **无需重下**：冻结副本 `input_links/p3_level2_qki_clip_peaks.bed`；原始补充表可从 GEO GSE147119 重取 |

## 3. 索引构建

```bash
# STAR 索引（试点期为 23GB 内存约束下建成的 24GB 索引；内存不足时降 SAindex）
STAR --runMode genomeGenerate --genomeDir STAR_index \
     --genomeFastaFiles GRCm38.primary_assembly.genome.fa \
     --sjdbGTFfile gencode.vM25.annotation.gtf \
     --runThreadN 8 --limitGenomeGenerateRAM 55000000000
# 若 <32GB 内存：追加 --genomeSAindexNbases 12

# salmon 索引
salmon index -t gencode.vM25.transcripts.fa -i salmon_index -k 31 --gencode

# DaPars2 3'UTR 注释（两步）：
# 第一步 GTF→for_DaPars2.bed：用 tools/DaPars2_patched/src/DaPars_Extract_Anno.py
#   （⚠ 无入库存档配方，曾因 Generate_Annotation.py 路径报错；可跳过——冻结副本已在
#    input_links/gencode_M25_3UTR_for_DaPars2.bed，sha256 见 results/reanalysis/00_inventory/）
# 第二步 for_DaPars2.bed→slim（每基因最长 3'UTR，21,158 条）：
python scripts/v5_1/make_slim_3utr_bed.py   # 输入=gencode_M25_3UTR_for_DaPars2.bed；产物=slim bed
```

⚠ **DaPars2 运行脚本的硬编码路径**：`p2_dapars2_run.sh` 等写死 `/home/taylor/DaPars2/src`、`/home/taylor/reference_v2/...`、`/mnt/d/stroke_apa_data/...`。重建时把 `tools/DaPars2_patched` 软链到 `~/DaPars2`、按 REGENERATE §3 重建索引到 `~/reference_v2/STAR_index`，或直接改脚本内路径。

## 4. GSE238125 FASTQ 重下与比对（主数据）

- 样本映射：`scripts/v5_1/p2_sample_map.tsv`（12 样本 ↔ ENA run ID）
- 预检：`python scripts/v5_1/p2_ena_preflight.py` → 下载：`python scripts/v5_1/p2_dl_pilot.py`（16 并发断点续传）或 `p2_download_fastq.py`
- 校验：`python scripts/v5_1/p2_md5_check.py`（对 ENA md5 全量核对，**24/24 必须全过**）
- 比对+覆盖+深度：`bash scripts/v5_1/p2_gse238125_align.sh`（fastp→STAR→sort→`genomecov -bg -split`→`view -c -F 0x904`；工作目录放 ext4）
- DaPars2：`bash scripts/v5_1/p2_dapars2_run.sh` → 全量 `p2_dapars2_full_parallel.sh`（应得 **9,693 事件**）
- 统计：`python scripts/v5_1/p2_sap_layer1_qc.py` → `Rscript scripts/v5_1/p2_sap_layer2_limma.R`（应得 **977 显著**）

## 5. 其余 GEO 数据集

| 数据集 | 用途 | 重取方式 |
|---|---|---|
| GSE174574（scRNA，P1 主数据） | 调控子层 | GEO 页面 supplementary → `https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE174574&format=file`；已删可重下 |
| GSE143531（海马 PAP，Task A） | 库水平 edgeR | `https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE143531&format=file`（RAW tar，48 文件）+ series matrix（冻结副本已在 `data_frozen/reanalysis_data/`） |
| GSE74456（皮层 PAPTRAP，Sakers 2017） | PAP 注释 Set1 | GEO supplementary（处理表已冻结在 `data_frozen/v5_1_data/GSE74456/`，无需重下即可跑 P4） |
| GSE286075（终足 RiboTag） | Set2 + 配对 DE | `data_frozen/v5_1_data/GSE286075/` 已含处理表 |
| GSE263986（皮层/白质 RiboTag） | Set3 | xlsx 处理表未入库（16MB，重取：GEO supplementary）；Set3 的 DE 结果表已在 `results/v5_1/p3_m4_set3_*` |
| GSE147119 原始补充（100MB） | QKI CLIP | 无需重下（peaks 冻结副本在 input_links） |
| GSE330741（SN-MPRA） | reporter 背景 | ⚠ manifest 矛盾未解，仅背景，不必重下 |

## 6. 未入库的分析中间件再生（全部确定性可重建）

| 被剔除文件 | 再生方式 |
|---|---|
| `results/reanalysis/02_candidate_rebuild/distal_nullpool.fa`（37M） | `python scripts/reanalysis/p2_extract_distal_fastas.py` + `p2_build_null_set.py`（seed 42，重跑结果逐字节一致） |
| `fimo_distal/`、`fimo_null/` 的 `cisml.xml`、`fimo.gff`（~79M） | meme 环境 `fimo --thresh 1e-4 --motif <frozen_rbp_motifs.meme 全库> <fa>`（fimo.tsv 主表已入库）；tsv 末尾 footer 行需跳过 |
| `results/reanalysis/03_bam_review/slices/`（60 bam+60 bai，IGV 便利件） | BAM 重建后 `samtools view -b regions.bed`（`results/reanalysis/03_bam_review/regions.bed` 已入库）或 `bash scripts/v5_1/p6_igv_slices.sh` |
| `data/bam/`（12 个排序 BAM，29G） | 第 4 节管线完整重跑（约 1–2 天） |
| `data/references/`（8.6G：genome.fa/STAR 索引/salmon 索引等） | 第 2–3 节 |
| `data/salmon/` | salmon index/quant 重跑 |
| `data/GSE143531_RAW.tar` + 解包目录 | 第 5 节 GEO 直链 |
| `data/tools/DaPars2/.git`（上游克隆历史）、`Dapars2_Test_Dataset.zip` | `git clone https://github.com/3UTR/DaPars2`（上游地址见 tools/DaPars2_patched/README.md）——**打过补丁的 src 本体已入库，勿用上游覆盖** |
| GSE330741 manifest 原件 | **从未获取**（即"manifest 矛盾"本身：651 基因库缺全部候选）。32 个 counts 文件已入库 `data_frozen/v5_1_data/GSE330741/`，reporter 对照表 `input_links/p4_mpra_counts_long.tsv`。待 Koester 2026 supplementary 到手后重建对照，无历史数据可丢 |
| 两个 `.git` 历史 | 不迁移（含 1GB 早期误提交的大对象冗余）；本仓库即干净历史 |

## 7. 验证锚点（重建完成的验收标准）

全部跑通后核对 `docs/reanalysis/实验流程.md` §0 表：9,693 事件 / 977 显著 / 17:17 单元测试 / Atp2a2 坐标三件套（distal 122,453,512–122,456,498；PAS 122,456,498；QKI FDR 1.04e-5）/ edgeR TMM 后 Atp2a2 log2FC −0.273 & FDR 0.683 / 引物零脱靶。
另有机械审计脚本可随时复跑：`python scripts/audit_upload_coverage.py`（对本地两工作区跑覆盖审计；本地删除后无意义，仅存档）。

## 8. Releases 依赖说明

四个交付 zip 只存在于 GitHub Releases（不在仓库正文）。若 Release 资产丢失：`audit_manifest.tsv` 中每个 RELEASE_ZIP 行的 note 写明了该文件在哪个 zip 内、当年内容是什么；其中分析核心（表/图/报告）都已同时入库 results/ 与 docs/，zip 主要增加 IGV 切片、PAS atlas、GSE143531 原表等大快照——均可按本手册重建。
