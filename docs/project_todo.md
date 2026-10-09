# 卒中—星形胶质细胞 APA—PAP 项目 TODO

版本：2026-09-21

## 1. 当前起点

### 1.1 已有结果

- GSE238125 的 DaPars2 逐样本 PDUI 矩阵、事件统计表、候选注释矩阵和候选卡；
- sham×2、day 3×2、day 7×2 的 Atp2a2、Ndrg2、Agpat3 区域 BAM/BAI；
- GENCODE vM25 的 DaPars2 3′UTR BED、QKI peak、样本对应表、QC 和运行日志；
- Atp2a2 event 6953 的时间轨迹：sham 0.19/0.20，day 3 0.08/0.06，day 7 0.09/0.09；
- 第二路线在 27 个可比较重点结果中有 13 个同向、14 个反向，Atp2a2 不在这 27 个结果中；
- GSE143531 可用于健康条件下的 PAP gene-level 优先级；
- GSE330741 可用于 reporter 方法、样本 manifest 和阳性对照设计。

### 1.2 已识别的坐标问题

Atp2a2、Ndrg2、Agpat3 均位于负链。现有 `lost_segments.bed` 使用了预测位点到高坐标端的区间。负链转录本因近端 PAS 使用而缺失的远端 3′UTR 位于低坐标端。该问题影响下列派生结果：

- lost segment 序列；
- QKI 与其他 RBP peak 重叠；
- motif 扫描和 motif 家族富集；
- 保守性评分；
- 基于上述字段生成的 candidate class 和候选卡。

逐样本 PDUI、ΔPDUI 和事件时间轨迹暂不因该区间方向问题改变。它们仍受样本量、短读长推断和位点真实性限制。

Atp2a2 的手工校正结果用于后续单元测试：

- 3′UTR：`chr5:122453512–122457153`，负链；
- predicted proximal PAS：`chr5:122456498`；
- 暂定远端特异区：`chr5:122453512–122456498`；
- 真正远端特异区内的 QKI peak：`chr5:122454200–122454250`；
- 旧结果中的 QKI peaks `chr5:122456550–122456650` 位于共有或保留侧。

### 1.3 当前结论边界

Atp2a2 保留为暂定主候选，依据是多个时间点的缩短预测、组内重复方向一致和 PAP gene-level 注释。QKI 证据需要使用校正后的区间全量重建。Ndrg2 与 Agpat3 保留为 APA 备选，其旧 QKI-lost-segment 证据停止使用。

## 2. 文件保留与新增输入

### 2.1 立即保留

| 文件或目录 | 用途 | 当前是否进入计算 |
|---|---|---|
| `p2_full_pdui_matrix.tsv` | 逐样本 PDUI 和时间轨迹 | 是 |
| `p2_sap_layer2_events.tsv` | 事件—对比统计和筛选状态 | 是 |
| `p2_second_tool_top40_signs.tsv`、`p2_second_tool_consistency.tsv` | 第二路线方向一致性 | 是 |
| `p3_p4_event_annotation_matrix.tsv` | 旧注释审计和重建前后比较 | 是，旧的链相关字段不可直接选靶 |
| `p4_candidate_cards.md` | 旧候选叙述审计 | 是，需重写 |
| `fimo_seqs/lost_segments.bed` | 定位坐标错误来源 | 仅审计 |
| DaPars2 3′UTR BED | 事件边界和链方向 | 是 |
| QKI peak BED | 重算远端区间重叠 | 是 |
| 六个区域 BAM/BAI | 主候选覆盖核查 | 是 |
| 样本表、参数、软件版本和日志 | 可追溯性 | 是 |
| 原始 FASTQ、全量 BAM/CRAM | 结果追溯和必要时局部重导出 | 归档，当前不重跑 |
| mm10/GRCm38 FASTA、GENCODE vM25 GTF/BED 的路径、版本和 checksum | 序列提取、链方向序列、IGV 参考一致性 | 是；集中存储时只保留路径和校验值 |

### 2.2 需要取得

| 新输入 | 最低要求 | 用途 |
|---|---|---|
| GSE143531 处理后表达表与样本信息 | PAP polysome、whole-astrocyte polysome、三个生物学样本和技术重复关系 | 评估候选基因在健康 PAP 中是否可测，生成 gene-level 证据 |
| GSE330741 完整 sample/barcode/tile manifest 与处理后 counts | sample、fraction、barcode、tile sequence、gene、replicate 对应明确 | 选择 reporter 阳性对照，复核分析设计 |
| 公共 PAS 注释 | PolyASite 或同类小鼠 mm10 位点表，带版本 | 判断预测位点与已知 cleavage site 的距离 |
| 独立动物样本 | sham 与 stroke，day 7 主时间点；先导每组 4 只 | 3′RACE、异构体比例和方差估计 |
| 湿实验元数据模板 | mouse、sex、age、batch、surgery、region、RNA QC、assay batch | 随机化、盲法、排除记录和统计建模 |

## 3. P0：冻结输入和建立可追溯清单

### 任务 P0.1：解包并生成文件清单

**工具：** `unzip`、`tar`、`shasum`、`samtools`。

**调用：**

```bash
unzip -l ATP2A2_IGV_CHECK_PACKAGE_20260920\(1\).zip > inventory_igv_zip.txt
unzip -l STROKE_APA_FINAL_DELIVERABLES_20260920.zip > inventory_final_zip.txt
shasum -a 256 ATP2A2_IGV_CHECK_PACKAGE_20260920\(1\).zip \
  STROKE_APA_FINAL_DELIVERABLES_20260920.zip > package_sha256.tsv
samtools --version > samtools_version.txt
```

**输出：** `inventory_igv_zip.txt`、`inventory_final_zip.txt`、`package_sha256.tsv`、`software_versions.tsv`。

**完成标准：** 文件清单、大小、checksum、来源路径和取得日期齐全。

**目的：** 固定分析输入，保证后续重建使用同一版本。

### 任务 P0.2：建立分析目录

```text
stroke_apa_reanalysis/
├── config/
├── input_links/
├── metadata/
├── scripts/
├── results/00_inventory/
├── results/01_strand_audit/
├── results/02_candidate_rebuild/
├── results/03_bam_review/
├── results/04_public_pap/
├── results/05_target_nomination/
├── wetlab/01_race/
├── wetlab/02_isoform_ratio/
├── wetlab/03_imaging/
├── wetlab/04_reporter/
└── logs/
```

原始文件通过只读链接进入 `input_links/`。结果文件写入新的版本目录。旧结果不覆盖。

## 4. P1：链方向和坐标体系审计

### 任务 P1.1：冻结坐标约定

建立 `config/coordinate_convention.yaml`：

```yaml
genome: GRCm38_mm10
annotation: GENCODE_vM25
bed_coordinates: 0_based_half_open
event_coordinate_source: DaPars2_output
plus_distal_rule: "[proximal_pas, utr_end)"
minus_distal_rule: "[utr_start, proximal_pas)"
plus_shared_rule: "[utr_start, proximal_pas)"
minus_shared_rule: "[proximal_pas, utr_end)"
```

同时核对 DaPars2 输出的 breakpoint 是否需要 1 bp 转换。判断依据包括 DaPars2 源代码、原始 BED 的长度关系和人工事件。任何转换写入独立字段：`proximal_pas_raw`、`proximal_pas_bed`。

### 任务 P1.2：编写链方向区间脚本

**脚本：** `scripts/01_build_strand_aware_segments.py`

**输入字段：** `event_id`、`chrom`、`utr_start`、`utr_end`、`proximal_pas_raw`、`strand`、`gene_symbol`。

**算法：**

```text
assert utr_start < proximal_pas_bed < utr_end
if strand == "+":
    shared  = [utr_start, proximal_pas_bed)
    distal = [proximal_pas_bed, utr_end)
if strand == "-":
    distal = [utr_start, proximal_pas_bed)
    shared  = [proximal_pas_bed, utr_end)
assert length(shared) + length(distal) == utr_end - utr_start
```

**调用：**

```bash
python scripts/01_build_strand_aware_segments.py \
  --events input_links/dapars2_event_coordinates.tsv \
  --utr-bed input_links/gencode.vM25.dapars2_3utr.bed \
  --coordinate-config config/coordinate_convention.yaml \
  --out-table results/01_strand_audit/strand_aware_event_segments.tsv \
  --out-distal-bed results/01_strand_audit/distal_segments.strand_aware.bed \
  --out-shared-bed results/01_strand_audit/shared_segments.strand_aware.bed
```

**输出：** `strand_aware_event_segments.tsv`、`distal_segments.strand_aware.bed`、`shared_segments.strand_aware.bed`、`invalid_event_coordinates.tsv`。

### 任务 P1.3：单元测试

测试集包含一个正链人工事件、一个负链人工事件、Atp2a2、Ndrg2、Agpat3。

Atp2a2 必须满足：

```text
distal_start = 122453512
distal_end   = 122456498
shared_start = 122456498
shared_end   = 122457153
```

**全量完成标准：**

- 每个区间长度为正；
- shared 与 distal 无内部重叠；
- 二者长度之和等于原 3′UTR 区间长度；
- event ID 唯一映射；
- strand 仅为 `+` 或 `-`；
- 无法匹配注释的事件单独输出；
- Atp2a2 单元测试通过。

**目的：** 将坐标方向写成可重复检查的程序规则。

## 5. P2：重建 QKI、motif、保守性和候选分级

### 任务 P2.1：重算 QKI peak 重叠

**工具：** `bedtools intersect`。

```bash
bedtools sort \
  -i results/01_strand_audit/distal_segments.strand_aware.bed \
  > results/01_strand_audit/distal_segments.sorted.bed

bedtools sort \
  -i input_links/p3_level2_qki_clip_peaks.bed \
  > results/02_candidate_rebuild/qki_peaks.sorted.bed

bedtools intersect -wa -wb \
  -a results/01_strand_audit/distal_segments.sorted.bed \
  -b results/02_candidate_rebuild/qki_peaks.sorted.bed \
  > results/02_candidate_rebuild/qki_distal_intersections.tsv

bedtools intersect -wa -wb \
  -a results/01_strand_audit/shared_segments.strand_aware.bed \
  -b results/02_candidate_rebuild/qki_peaks.sorted.bed \
  > results/02_candidate_rebuild/qki_shared_intersections.tsv
```

**完成标准：**

- Atp2a2 distal 区间命中 `chr5:122454200–122454250`；
- `chr5:122456550–122456650` 只进入 Atp2a2 shared 结果；
- Ndrg2、Agpat3 的旧 QKI 证据重新判定；
- 输出每个事件的 distal peak 数、最小 FDR、总覆盖 bp 和 shared-side peak 数。

### 任务 P2.2：提取链方向序列

**工具：** `bedtools getfasta`。

```bash
bedtools getfasta \
  -fi references/mm10.fa \
  -bed results/01_strand_audit/distal_segments.sorted.bed \
  -s -name \
  -fo results/02_candidate_rebuild/distal_segments.strand_aware.fa
```

`-s` 用于负链区间的反向互补。参考 FASTA 在此步骤承担序列提取和参考一致性作用。参考位于公共目录时，保留绝对路径、`.fai`、版本和 checksum，无需复制进交付 ZIP。

### 任务 P2.3：重跑 motif 分析

**工具：** MEME Suite `fimo`。沿用原先冻结的 motif 数据库、阈值和匹配背景。原背景受到同一坐标问题影响时，同步重建长度与 GC 匹配的 null segments。

```bash
fimo --thresh 1e-4 \
  --oc results/02_candidate_rebuild/fimo_distal \
  input_links/frozen_rbp_motifs.meme \
  results/02_candidate_rebuild/distal_segments.strand_aware.fa
```

对每个 motif family 计算候选 distal segments 与长度/GC 匹配 null segments 的命中比例，沿用冻结的原统计方法或 Fisher exact test，并对 motif families 使用 Benjamini–Hochberg 校正。

**输出：** `fimo.tsv`、`motif_event_summary.tsv`、`motif_family_enrichment.tsv`、`motif_null_matching_qc.tsv`。

### 任务 P2.4：重算保守性

**工具：** `bigWigAverageOverBed`，输入使用与 mm10 匹配的 phyloP 或 phastCons bigWig。

```bash
bigWigAverageOverBed \
  references/mm10.phyloP_or_phastCons.bw \
  results/01_strand_audit/distal_segments.sorted.bed \
  results/02_candidate_rebuild/distal_conservation.tsv
```

同时报告均值、覆盖碱基数和区间长度。覆盖不足的事件标记 `NA`。

### 任务 P2.5：重建事件注释矩阵和候选卡

**脚本：** `scripts/02_rebuild_candidate_matrix.R`

```bash
Rscript scripts/02_rebuild_candidate_matrix.R \
  --pdui input_links/p2_full_pdui_matrix.tsv \
  --sap input_links/p2_sap_layer2_events.tsv \
  --segments results/01_strand_audit/strand_aware_event_segments.tsv \
  --qki results/02_candidate_rebuild/qki_distal_intersections.tsv \
  --motif results/02_candidate_rebuild/motif_event_summary.tsv \
  --conservation results/02_candidate_rebuild/distal_conservation.tsv \
  --pap input_links/pap_gene_sets.tsv \
  --out results/02_candidate_rebuild/event_annotation_matrix.v2.tsv
```

**必须同时输出：** `event_annotation_matrix.v2.tsv`、`candidate_cards.v2.md`、`old_vs_new_annotation_diff.tsv`、`events_changed_qki_status.tsv`、`events_changed_class.tsv`、`rebuild_session_info.txt`。

**候选排序字段：**

1. 事件的重复一致性和效应量；
2. 多时间点稳定性；
3. 预测 PAS 与公共 PAS 的距离；
4. 区域覆盖和探针可设计性；
5. PAP gene-level 可测性；
6. 校正后的 distal RBP/motif/保守性；
7. 第二方法是否直接覆盖该事件。

RBP、motif 和保守性承担优先级信息，不单独承担定位机制结论。

## 6. P3：BAM 与 IGV 事件核查

### 任务 P3.1：BAM 完整性

```bash
samtools quickcheck -v results/03_bam_review/*.bam \
  > results/03_bam_review/samtools_quickcheck.txt

for bam_file in results/03_bam_review/*.bam; do
  samtools idxstats "${bam_file}" >> results/03_bam_review/idxstats_all.tsv
done
```

**完成标准：** 六个 BAM 可读取，BAI 与 BAM 匹配，目标染色体存在。失败样本从全量 BAM 按同一坐标重新导出。

### 任务 P3.2：建立 IGV session

**IGV 设置：**

- genome：mm10/GRCm38，与 BAM 和 BED 一致；
- track：GENCODE vM25、predicted PAS、校正后的 distal/shared BED、QKI peaks；
- coverage scale：同一候选的所有样本固定相同；
- read filter：mapping quality、duplicate、secondary/supplementary 规则固定；
- display：按 sham、day 3、day 7 排列，每组两个生物学重复；
- snapshot：全 3′UTR、predicted PAS ±500 bp、distal/shared 边界各一张。

**每个候选记录：**

- 覆盖变化是否位于预测位置附近；
- 两个重复是否同向；
- 是否存在单样本局部尖峰；
- 是否有重叠基因、可变剪接或末端外显子切换；
- PAS 附近是否存在基因组富 A；
- 注释末端与实际覆盖末端是否一致。

### 任务 P3.3：定量共有区与远端区覆盖

**工具：** `samtools depth` 或 `mosdepth`。所有样本使用相同 MAPQ、base quality 和 duplicate 规则。

```bash
samtools depth -aa -Q 20 -q 20 \
  -b results/01_strand_audit/distal_segments.top_candidates.bed \
  sample.bam > sample.distal.depth.tsv

samtools depth -aa -Q 20 -q 20 \
  -b results/01_strand_audit/shared_segments.top_candidates.bed \
  sample.bam > sample.shared.depth.tsv
```

对每个样本计算：

```text
mean_distal_coverage
mean_shared_coverage
coverage_ratio = mean_distal_coverage / mean_shared_coverage
covered_base_fraction at depth >= 1, 5, 10
```

覆盖比用于读段证据描述。正式 APA 比例以 3′RACE 后的异构体检测为准。

### 任务 P3.4：事件审查表

为 Atp2a2、Ndrg2、Agpat3 和重建后新增的前五候选生成一行一事件的审查表：

| 字段 | 内容 |
|---|---|
| event/gene | 唯一标识 |
| strand-aware distal | 校正后的区间 |
| PDUI trajectory | sham 和各时间点 |
| replicate concordance | 两只动物方向 |
| BAM coverage | distal/shared 比和图 |
| PAS support | 已知 PAS 距离、富 A 风险 |
| second-tool status | direct、indirect、absent |
| PAP evidence | gene-level 数据来源 |
| probe feasibility | common/distal 可设计长度和特异性 |
| decision | advance、hold、drop |

**P3 完成标准：** 至少一个候选通过读段和探针可设计性核查；Atp2a2 的 advance、hold 或 drop 有明确证据记录。

## 7. P4：公共 PAP 与 SN-MPRA 辅助分析

### 任务 P4.1：GSE143531

**输入：** 处理后表达矩阵、样本信息、技术重复关系。

**分析：**

1. 核对 C1/C2/C3 生物学样本及技术重复；
2. 技术重复在同一生物学样本内合并，或在模型中保留重复层级；
3. 计算候选基因在 PAP polysome 与 whole-astrocyte polysome 中的表达和配对差；
4. 输出可测性、效应方向和三个生物学样本的一致性；
5. 结果标注为 gene-level，不推断长短异构体。

**输出：** `gse143531_candidate_gene_summary.tsv`、QC 图、样本关系图。

### 任务 P4.2：GSE330741

**输入：** 最新 GEO manifest、barcode/tile map、processed counts。

**分析：**

1. 核对 GEO 更新后的 sample—sequence 对应；
2. 建立 `sample_id → fraction → replicate → barcode → tile → gene` 长表；
3. 复算论文定义的 reporter 定位分数或使用作者处理后分数；
4. 选择 Glt1/Slc1a2 或 Sparc 中重复一致、效应稳定的片段作为阳性对照；
5. 保留 matched neutral tile 和 barcode-level 阴性对照。

**输出：** `gse330741_manifest_verified.tsv`、`reporter_positive_controls.tsv`、`reporter_negative_controls.tsv`。

## 8. Gate 1：湿实验靶标确定

### 8.1 必须满足

- 链方向重建和单元测试通过；
- BAM 未显示明显比对或注释伪影；
- 两个发现样本的方向一致；
- 校正后的 distal 区间可设计特异探针或引物；
- 预测 PAS 附近不存在足以解释 3′RACE 的严重内部引物结合风险；
- day 7 作为主时间点具备同向效应。

### 8.2 加分项

- 多时间点同向；
- 公共 PAS 支持；
- GSE143531 中 gene-level PAP 可测；
- 校正后的 distal 区间有 RBP、motif 或保守性证据；
- 第二方法直接覆盖并同向。

### 8.3 Gate 1 输出

- `target_nomination_v1.md`；
- 1 个主靶；
- 2 个备靶；
- 每个靶标的 common 区、distal 区、预测 PAS、引物区和探针区；
- advance、hold 或 drop 的理由。

## 9. P5：3′RACE 确认真实 RNA 3′端

### 任务 P5.1：样本设计

- 主时间点：day 7；
- 分组：sham、stroke；
- 先导：每组 4 只独立小鼠；
- 样本：同侧预设皮层区域的星形胶质细胞 RNA；
- 记录：性别、年龄、手术批次、梗死或神经功能指标、取材时间、RNA integrity；
- 随机化：按手术批次分层；
- 盲法：样本编码后完成 PCR、测序和峰图判读。

### 任务 P5.2：3′RACE

**实验设计：**

- 锚定 oligo(dT) 接头引物；
- 候选基因特异正向引物位于预测 proximal PAS 上游；
- 第二套 nested primer；
- no-RT、no-template 和已知 3′端阳性对照；
- 产物用 Sanger 或靶向 amplicon sequencing；
- 判定 cleavage site 时要求基因特异序列后紧接非模板 poly(A) 尾。

**保留文件：**

- primer table 和批号；
- 原始电泳图；
- `.ab1` 或 FASTQ；
- base-called FASTA；
- 比对文件；
- 每只动物的 cleavage site 计数；
- PAS 上下游 50 nt 序列及内部富 A 评估。

**输出：** `validated_cleavage_sites.tsv`、`race_read_counts_per_mouse.tsv`、位点示意图。

### Gate 2A

主靶进入比例测定需确认至少两个 RNA 3′端，其中至少一个近端和一个远端位点具有独立测序读段支持。缺少真实 poly(A) 接合的条带不计入位点。

## 10. P6：独立动物异构体比例

### 任务 P6.1：检测设计

- `common assay`：位于 proximal PAS 上游，测 total candidate RNA；
- `distal assay`：位于 distal-only 区，测 long isoform；
- `proximal junction assay`：仅在真实 cleavage junction 允许特异设计时测 short isoform；
- 平台优先级：ddPCR，其次为完成扩增效率、线性范围和标准曲线验证的 qPCR；
- 每只动物输出 long、total、long/total 和检测下限状态。

### 任务 P6.2：统计

主比较为 day 7 stroke 与 sham 的每只小鼠 `long/total`。报告组间差、95% CI 和原始点。手术批次作为固定或随机项依据批次数量预先确定。技术孔先在每只动物内汇总。

先导每组 4 只用于估计方差和失败率。正式样本量使用下式或等价软件计算：

```text
n_per_group = 2 × (z_0.975 + z_0.80)^2 × sigma^2 / delta^2
```

其中 `sigma` 来自盲态或先导方差，`delta` 为预先定义的最小重要差异。最终增加预登记的手术死亡和 assay failure 比例。

### Gate 2B

进入空间实验需满足：

- stroke 效应与发现数据同向；
- 结果未由单只动物或单一批次驱动；
- 长异构体位于可可靠定量范围；
- 预先冻结的效应量和不确定性门槛达到进入条件。

## 11. P7：长、短异构体的组织空间分布

### 任务 P7.1：探针验证

- common probe：测全部候选 RNA；
- distal probe：测 long isoform；
- 已知长、短构建体或合成 RNA：估计检出率、串色、交叉反应和双通道共定位率；
- 双通道性能不足时转用 BaseScope、padlock 或 cleavage-junction 检测；
- 探针批次在 sham 和 stroke 间平衡。

### 任务 P7.2：组织和区室

- 星形胶质细胞：膜标记或遗传 reporter，用于细突起分割；
- 胞体：核周和细胞体预设区域；
- 突触邻近突起：膜标记与预设突触标记的距离规则；
- 血管周终足：膜标记、CD31/lectin 和 AQP4 辅助定义；
- 读取组别前冻结分辨率、距离阈值、分割参数和排除规则；
- 部分样本使用 Airyscan、SIM、扩增显微或同等级平台确认空间定义。

### 任务 P7.3：主要终点

每只小鼠先汇总：

```text
long_fraction_PAP  = long_PAP / total_PAP
long_fraction_soma = long_soma / total_soma
delta_spatial      = long_fraction_PAP - long_fraction_soma
```

主模型：

```text
long_fraction ~ condition * compartment + batch + (1 | mouse)
```

计数型数据采用与分布匹配的广义线性混合模型。主要检验为 `condition × compartment`。切片、视野和细胞属于小鼠内嵌套观测。

### 任务 P7.4：同步控制

- 单位体积或突起长度的 total RNA；
- 突起长度、分支、体积和突触邻近面积；
- 总候选 RNA；
- 成像深度和背景；
- 血管周终足的独立区域结果。

### Gate 3

| 结果 | 处理 |
|---|---|
| APA 比例复现，`condition × compartment` 有明确效应 | 进入 reporter 和内源 PAS 操纵 |
| APA 比例复现，空间效应接近零且 CI 排除预设效应 | 论文主线收敛为卒中相关 APA，停止定位机制 |
| 空间变化与 total RNA 或形态同步 | 增加归一化和形态匹配，限制空间分配结论 |
| 信号低于检测能力 | 改进检测平台；一次优化后仍失败则换备靶 |

## 12. P8：顺式 3′UTR reporter

### 任务 P8.1：构建体

1. 3′RACE 确认的长 3′UTR；
2. 3′RACE 确认的短 3′UTR；
3. 长 3′UTR 删除候选元件；
4. 短 3′UTR 回补同一元件；
5. GSE330741 复核后的 PAP 定位阳性对照；
6. 长度和 GC 匹配的中性对照。

全部构建使用相同启动子、编码区、载体骨架、标签和终止设计。构建序列进行 Sanger 或全质粒测序。

### 任务 P8.2：实验顺序

1. 原代星形胶质细胞完成表达、稳定性和定位预筛；
2. 器官型切片或体内 astrocyte-specific 系统确认组织背景；
3. 使用 smFISH 或 barcode 检测 reporter RNA；
4. 同步测定 reporter 总 RNA、RNA 半衰期和蛋白；
5. 比较长与短、删除与长、回补与短三个预设对比。

### Gate 4

候选元件进入内源实验需在控制总 RNA 后改变空间分布，并在删除与回补方向上形成一致证据。

## 13. P9：内源 PAS 因果操纵

### 任务 P9.1：策略筛选

优先测试 proximal PAS 阻断寡核苷酸。编辑策略仅在星形胶质细胞递送和脱靶评估可行时启动。每项操作均设置 scrambled 对照和同区域不影响 PAS 的序列对照。

### 任务 P9.2：验证顺序

1. 3′RACE 确认 cleavage site 未产生异常新位点；
2. ddPCR/qPCR 确认 long/total 改变；
3. 原位检测 `Δspatial`；
4. rescue 恢复 PAS 使用或回补顺式元件；
5. 同步记录 total RNA、细胞活性和星形胶质细胞形态。

### Gate 5

内源 long/short 比例改变并引起同方向空间变化，rescue 使两者恢复，构成 APA 参与 RNA 空间分配的因果证据。

## 14. P10：条件性机制和功能

### 14.1 QKI

启动条件：校正后的 distal 区间含 QKI peak，目标区在星形胶质细胞中有 QKI 结合证据，顺式元件影响定位。

任务包括 QKI 蛋白检测、候选 RNA 结合验证、QKI 扰动、异构体比例、空间分布和 rescue。QKI 对 PAS 选择与对 RNA 运输的作用分别报告。

### 14.2 PABPN1

启动条件：独立样本中存在近端 PAS 增加，且 PABPN1 蛋白量或核定位提示相关变化。任务包括 PABPN1 扰动、PAS 使用、rescue 和细胞状态控制。

### 14.3 Atp2a2 功能分支

启动条件：Atp2a2 通过事件、空间和内源因果三层验证。终点包括局部 SERCA2 分布和星形胶质细胞突起的 cytosolic 或 ER Ca²⁺ 动态。功能实验比较恢复长异构体与维持短异构体，保留总 SERCA2、细胞活性和形态指标。

## 15. 数据分析和报告规则

### 15.1 实验单位

- 公共 GSE238125：小鼠样本；
- 独立动物实验：小鼠；
- 技术孔：先在小鼠内汇总；
- 切片、视野、细胞：小鼠内嵌套观测；
- reporter barcode：构建体内技术层级，生物学重复单独建模。

### 15.2 报告内容

- 原始点、效应量、95% CI、确切 `n`；
- 随机化、盲法、排除、死亡和 assay failure；
- 主要终点和预设对比；
- 多重比较的检验家族和校正方法；
- 缺失数据与低于检测限的处理；
- 软件、版本、参数、脚本 commit 和 session info。

### 15.3 结论层级

| 数据层 | 可以支持的结论 |
|---|---|
| DaPars2 + BAM | 卒中相关 APA 候选及其读段基础 |
| 3′RACE + 独立比例 | 真实 RNA 3′端和卒中相关异构体比例变化 |
| 双探针组织成像 | 长异构体在胞体、PAP 和终足的相对分布 |
| Reporter 删除/回补 | 远端片段中的顺式作用 |
| 内源 PAS 操纵 + rescue | APA 事件参与空间分配的因果作用 |
| QKI/PABPN1 扰动 | 特定调控因子的机制作用 |
| 局部 SERCA2/Ca²⁺ | Atp2a2 路线的局部功能后果 |

## 16. 交付物清单

### 计算阶段

- [ ] 输入清单和 checksum；
- [ ] 坐标约定文件；
- [ ] 链方向区间脚本和测试；
- [ ] `strand_aware_event_segments.tsv`；
- [ ] distal/shared BED；
- [ ] 重算 QKI、motif、保守性结果；
- [ ] `event_annotation_matrix.v2.tsv`；
- [ ] `old_vs_new_annotation_diff.tsv`；
- [ ] 候选卡 v2；
- [ ] 固定尺度 IGV 截图；
- [ ] BAM 审查表；
- [ ] GSE143531 gene-level 表；
- [ ] GSE330741 manifest 和 reporter 对照表；
- [ ] `target_nomination_v1.md`。

### 湿实验阶段

- [ ] 3′RACE 原始文件和 cleavage-site 表；
- [ ] 每只动物的 long/total；
- [ ] 样本量计算和预登记统计计划；
- [ ] 双探针验证结果；
- [ ] 盲法组织成像原始图和分割对象；
- [ ] mouse-level 空间终点表；
- [ ] reporter 序列、质控和定位结果；
- [ ] 内源 PAS 操纵和 rescue；
- [ ] 条件性 QKI/PABPN1 与功能结果。

## 17. 建议执行顺序

| 顺序 | 工作 | 预计用时 | 进入下一步的条件 |
|---:|---|---:|---|
| 1 | P0 输入冻结 | 0.5–1 天 | 清单和 checksum 完成 |
| 2 | P1 链方向重建 | 1–2 天 | 单元测试和全量 QC 通过 |
| 3 | P2 注释与候选分级重建 | 2–4 天 | 差异表和候选卡 v2 完成 |
| 4 | P3 BAM/IGV 核查 | 1–2 天 | 至少一个候选通过 |
| 5 | P4 公共 PAP/MPRA 辅助分析 | 2–5 天 | target nomination 完成 |
| 6 | P5 3′RACE | 3–6 周 | 真实近端和远端位点确认 |
| 7 | P6 独立动物比例 | 3–6 周，可与 P5 部分并行 | 同向效应达到 Gate 2B |
| 8 | P7 组织空间成像 | 2–4 月 | Gate 3 决定论文路线 |
| 9 | P8 reporter | 2–4 月 | 顺式元件通过 Gate 4 |
| 10 | P9 内源 PAS | 3–6 月 | 因果证据达到 Gate 5 |
| 11 | P10 机制和功能 | 3–6 月 | 相应启动条件满足 |

前五步完成前不订购大规模探针、AAV 或 reporter 文库。3′RACE 的少量引物可在主靶通过 BAM 核查后设计。

## 18. 最近一轮应完成的工作

- [ ] 从最终交付包和 IGV 包整理输入链接；
- [ ] 冻结 0-based、half-open 坐标约定；
- [ ] 实现正链和负链 distal/shared 区间生成；
- [ ] 用 Atp2a2 的两个 QKI 区域完成单元测试；
- [ ] 全量重算 QKI overlap；
- [ ] 重算 motif 和保守性；
- [ ] 生成旧结果与新结果差异表；
- [ ] 更新候选卡；
- [ ] 在 IGV 中完成 Atp2a2、Ndrg2、Agpat3 六样本核查；
- [ ] 输出主靶、备靶和 3′RACE 引物设计区间。

这一轮结束时应得到一个可执行决定：Atp2a2 进入 3′RACE、进入保留观察，或由链方向重建后的候选替换。
