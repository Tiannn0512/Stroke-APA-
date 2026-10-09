# 课题全程报告：从原始设计到干实验收官

> 生成：2026-09-21 ｜ 覆盖：v5.1 干实验阶段（09-16 至 09-20）＋ 链方向重建阶段（09-21）＋ Gate 1 定靶
> 工作区：`D:\stroke_apa\`（v5.1 档案，git）＋ `D:\stroke_apa_reanalysis\`（重建，git）
> 本报告初稿完成后经五轮自校验，修正问题后发布（校验记录见第 9 节）

---

## 1. 课题设计与科学问题（原始假设 → 四幕演化）

### 1.1 原始设计（pipeline v5.0/v5.1，2026-09-16 定版）

科学问题：**卒中后，星形胶质细胞内 RNA 加尾机器（APA 调控子）改变 → 特定 mRNA 的 3'UTR（尾巴）被剪短 → 这些 RNA 的远端隔室定位潜能改变？**（全程核心指标 PDUI = 远端加尾位点的使用比例，越高代表长尾巴版本占比越高）

原始设计为六环节链：M1 调控子层（scRNA 找差异表达的调控子）→ M2 事件层（全转录组测 APA 事件）→ M3 关联层（调控子↔事件的三级证据）→ M4 定位层（剪短事件是否偏向某个隔室的基因集）→ M5 证据分级 → M6 交付。当时关注的隔室是**血管终足**。

设计原则（贯穿全程）：**先冻结后看数**——每层的关键参数在运行前写入冻结文件并 commit；结论分级（直接数据/正交证据/预测）；一切负结果照登；不越级表述。

### 1.2 演化过程中的四次方向修正（全部由数据或用户追问触发）

1. **M2 主数据源切换（Gate GA，09-16）**：原计划的单细胞 APA 工具因细胞数临界与工具风险，切换为 Plan B（GSE238125 分选星形胶质细胞 bulk，12 样本时序）。
2. **定位轴切换（09-20）**：终足轴全阴性，PAP（突触旁突起）轴阳性——定位假设从终足修正为突触旁。
3. **问题磨尖（09-20–21，用户追问驱动）**：从"这批基因偏向哪里"收窄为"**同一基因的长、短两种 mRNA 在胞体与突触旁的分配是否因卒中改变**"，主案例 Atp2a2。
4. **链方向错误重建（09-21）**：发现 v5.1 对负链基因的"丢失段"定义在错误一侧，全部序列级证据推倒重建；重建后 **Atp2a2 升级为主靶**（用户批准）。

4. **收官期附 C 事件（v5.1 关闭后、重建启动前）**：P6-01 产出 IGV 读段核查包（29,569,120 B）；P6-02/07 两轮数据清理（删 FASTQ 45GB+GSE174574 2.6GB、bedGraph 备份/phastCons 中间件/STAR 索引等约 73GB→搬迁至重建区，WSL 侧原件清除）；P6-03 采纳用户定稿的研究方案 v3（三关设计，即本轮重建+湿实验的总框架）；v5.1 版计划的 P6-04/P6-05（原 Gate 1/Gate 2 计划）由重建阶段 R3/R4 实际承接。


---

## 2. 数据资源

| 数据集 | 内容 | 用途 | 状态 |
|---|---|---|---|
| GSE174574 | 小鼠 MCAO（大脑中动脉阻塞中风模型）24h 皮层 scRNA（10x 3'，Sham3+MCAO3） | M1 调控子层 | 已下载（原始文件 09-20 已删，公开可重下） |
| GSE238125 | 皮层星形胶质细胞 bulk RNA-seq，Sham×2 + d1/d3/d7/d21/d60×2，永久性远端 MCA 结扎+7min 双侧 CCA 阻断 | M2 事件层主数据 | FASTQ 已删（md5 元数据保留）；BAM 12 个在重建区 |
| GSE74456 | Sakers 2017 PAPTRAP：健康皮层星形胶质细胞**突触旁突起 vs 全皮层**翻译组 | M4 PAP 轴参考 | 处理表已入库 |
| GSE286075 | 血管终足 RiboTag 转译组（CTR/Stroke 各 3） | M1-L2 配对 DE 第二证据层 ＋ M4 终足轴参考 | 已入库 |
| GSE263986 | mGFAP-RiboTag 皮层/白质分区转译组（每区 n=5+5） | M4 空间分区轴 | xlsx 已入库 |
| GSE147119 | QKI CLIP（Sakers 2021，P21 全前脑），437 peaks | M3-L2 结合证据 | peaks bed 已入库 |
| GSE330741 | SN-MPRA（定位序列筛查） | reporter 对照 | ⚠ manifest 矛盾未解，仅背景 |
| GSE143531 | Mazaré 2020 健康**海马**PAP vs 全星胶翻译组（3 动物×2 分区×8 技术重复） | P4.1 可测性核对 | 已下载分析 |
| PolyASite 2.0 | mm10 已知 poly(A) 位点图谱，301,006 clusters | PAS 支持/距离 | 已下载（新址 unibas.ch） |
| mm10.60way.phastCons.bw | 60 物种保守性 bigWig（4.58GB） | L4 保守层 | 并行分块重下 |

参考与工具：GRCm38/mm10 genome.fa、GENCODE vM25 GTF/genePred、DaPars2 3'UTR 注释（tandem 版 3.9MB）、salmon 索引、STAR 索引（v5.1 期 24GB，重建期已删可重建）、WSL conda envs（dapars2/meme/rlimma/bioinfo/salmon）。

---

## 3. 阶段一：v5.1 干实验（P0–P5，09-16 至 09-20）

### 3.1 P0 数据底座与启动审计

**P0-01 建库判定**：查原文+SRA 元数据确认 GSE174574 为 10x 3' droplet（6/6 runs PAIRED+cDNA+HiSeq 4000）→ M2 原工具锁定 scAPAtrap/Sierra，为后续切 Plan B 埋点。
**P0-04/05 下载**：并行 Range 分块下载（16×8MB×重试），字节数与 filelist 精确核对；tar 成员逐一验证。
**P0-08/09 启动审计**：读 mtx 聚合基因×样本矩阵，确认 5 个核心调控子全部可检出且 Qk 下调（42.8%→34.2% 表达率）；astrocyte 细胞数 3,045（接近每样本 500 的预设下限）→ 触发 Gate GA（该阈值为 v5.0 设计期预设，非 P1 冻结文件内容）。

**Gate GA（用户决策）**：M2 切 Plan B（GSE238125 bulk + DaPars2），不装 R 不跑单细胞 APA。

### 3.2 P1 调控子层（scRNA）

**P1-01 冻结参数**（`results/P1_frozen_params.md`，DE 计算前 commit）：

| 参数 | 值 | 理由 |
|---|---|---|
| QC 过滤 | genes≥200, counts≥500, mt%≤20 | 10x 常规底线，防低质量细胞/空滴 |
| 双联体 | Scrublet（sim_ratio=2, n_neighbors=30）；不可用时降级 MAD±3 | 冻结时就写明降级路径——后来 scrublet 在 py3.13 装不上，按预案降级并登记 |
| 标准化 | normalize_total(1e4)→log1p→HVG2000→PCA50→neighbors(15,30 PCs)→leiden(1.0) | scanpy seurat 惯例；分辨率 1.0 取中等粒度 |
| astro 簇判定 | Aqp4≥1 且 Slc1a3≥1；或 Aqp4≥0.5 且 Gfap≥1 且 Slc1a3≥0.5 | 双规则覆盖稳态簇与反应性簇；微胶/内皮标志共高则排除 |
| reactive 评分 | score_genes(Serpina3n,C3,Gfap)，>0.5 为 reactive | 最小基因集防过拟合 |
| 伪bulk 下限 | 每单元<20 细胞不进 DE | 防小单元方差爆炸 |
| DE | pydeseq2 ~condition（n=3+3）；GSE286075 用 ~mouse+condition 配对 | 配对设计利用同鼠对照 |

**结果**：QC 通过 97.2%；30 簇中命中第 5、21 两个簇为 astrocyte（合计 3,045 细胞）；overall DE 1,280 sig；GSE286075 配对 DE 一度算出"深度极低仅 6 个 sig"，**勘误**后发现取错链特异性列（col3→col4，dUTP 第二链），重跑得 409 sig；regulator 目录 31 个，17/31 两层方向一致，最强信号 Pabpn1↓（padj 3e-4）。
**P1-07 敏感性**：双联体（2.85%）与 ambient（5.0%）过滤后 astro 数 3,046/3,042（vs 3,045），31 调控子方向翻转 0 → 结论稳健。
**P1-08 复现附录**：按原论文判据精确复现 10,862 基因入模、54/55 方向零矛盾 → 205 vs 409 = 判据差异非管线错误。
**Gate G1（用户确认）**：≥3 个调控子两层方向一致（17/31 大幅满足）→ 通过；G1 证据指向后续 L1 文献池应重点收 PABPN1/QKI 扰动证据。

### 3.3 P2 数据获取与比对

**下载**：ENA 慢速（0.2–0.6MB/s，`p2_dl_pilot.py` 16 并发/600s 超时断点续传）→ 用户快渠道一夜到齐 24 文件；`p2_md5_check.py` 对 ENA md5 基准全量校验 **24/24 通过**（中途 sham1_2 因跨会话分块损坏 md5 不符，立即重拉，坏文件未进下游）。

**比对管线**（`p2_gse238125_align.sh`，语义 09-16 冻结）：

```bash
fastp --thread 8 --detect_adapter_for_pe      # PE 接头自动识别；8 线程=12核留余量
STAR --runThreadN 8 --genomeDir <pilot 24GB index>
     --outSAMtype BAM Unsorted --quantMode GeneCounts
     --outTmpDir /home/taylor/STARtmp/...      # 临时目录必须在 ext4：NTFS 上无法建 FIFO
samtools sort -@ 4; samtools index
bedtools genomecov -ibam X -bg -split          # 关键：-split 按剪接拆分覆盖，防内含子读段污染外显子覆盖
depth = samtools view -c -F 0x904              # 0x904 = 去除未比对(0x4)+次要比对(0x100)+补充(0x800)：主比对口径
```

参数理由：复用试点期已验证的 24GB STAR 索引（GRCm38+GENCODE vM25，23GB 内存约束下建成）；`-split` 语义来自试点期 reanalysis_v3 纠错（不加会高估内含子覆盖、扭曲 3'UTR 覆盖剖面）；工作目录全部放 ext4（/mnt/d 的 9p 协议拖慢 fastp/STAR 3–10 倍，实测 sham1 比对 27min→快路径）。**结果**：12/12 比对，唯一比对率 79–96.4%，bedgraph 12 份，深度表 12 行（曾因链中断丢表误填 8 个占位值，当即发现并全部替换为日志真值）。

**DaPars2 定量**（`p2_dapars2_run.sh` + `p2_dapars2_full_parallel.sh`）：

| 参数 | 值 | 理由 |
|---|---|---|
| 注释 | gencode_M25_3UTR_for_DaPars2.bed（tandem 3'UTR）→ 后换 slim 版 21,158 条（用户批准） | DaPars2 需要 3'UTR 级注释；slim=每基因最长 3'UTR，事件定义收窄为"最长注释末端内的断点移位"（表述边界见 SAP §4） |
| Coverage_threshold | 10 | DaPars2 默认；低于此覆盖的 PDUI 估计不稳 |
| Num_Threads | 8（单进程）/ per-chr 8 进程并行 | 12 核机器留余量；逐染色体独立实例规避 multiprocessing Manager 字典不稳 |
| 源码补丁 | 空覆盖事件返回 NA（.bak 存档） | 修复空覆盖除零崩溃；统计语义不变 |
| 深度表口径 | -F 0x904 主比对 | 与比对深度统计同一语义，保证 PDUI 归一一致 |

试点（4 样本）先行：9,269 事件、sham 重复相关 0.919、d1 缩短:延长 = 3.4:1（显著事件数比，\|ΔPDUI\|>0.1 口径）——方向与 Pabpn1↓ 文献预测咬合。全量 12 样本：**9,693 事件**（depth 22/22 日志核验：合并脚本对"12 份比对深度日志 + 逐染色体定量日志"与深度表行的逐条对账计数，脚本 p2_full_merge_qc.py）。

### 3.4 P2 统计（SAP v1.0→v1.1，全量统计前冻结）

冻结文件 `results/P2_SAP_frozen.md`（09-19 commit，早于全量统计；v1.1 经用户批准改用 R limma）：

| 参数 | 值 | 理由 |
|---|---|---|
| 事件资格 | 每 contrast 需 sham1+2+时间点 rep1+2 **4/4 有效 PDUI**，缺失=not_testable 不插补 | 插补会制造假信号 |
| 变换 | logit((PDUI·(1−2ε)+ε))，**ε=0.01** | PDUI 有 0/1 边界值，logit 发散；ε 内缩边界 |
| 模型 | limma `~0+timepoint`，**eBayes(robust=TRUE, trend=TRUE)** | robust=抗离群事件方差收缩；trend=均值-方差趋势校正（RNA-seq 必要）；v1.1 弃纯 Python 重实现方案改 R 原生（零实现风险，用户批准） |
| 检验顺序 | 先 omnibus 联合 F（5 水平），再 5 个预设 contrast（各 vs sham） | omnibus 控家族、contrast 出方向 |
| 多重校正 | **全部事件×5 contrast 联合 BH，FDR<0.1**（各 contrast 内 BH 作辅助） | 联合校正最保守，防跨时间点选择性报告 |
| 效应阈值 | **\|ΔPDUI\|≥0.1**（预冻结） | 生物学可解释的最小移位 |
| 显著定义 | 联合 FDR<0.1 且 \|ΔPDUI\|≥0.1 且两重复同向 | 三条件缺一不可 |
| 一致性标注 | ≥2 时间点同向（共享 sham，明示非独立复现） | 防止把相关时点当独立验证 |
| 层 1 QC | median/trimmed PDUI、分层（表达四分位/coverage/UTR 长度）、12×12 相关热图+PCA、LOSO、混杂排查 | d1 缩短进正文的必要条件：分层后方向不变。⚠已定位 day3_rep1 文库异常（insert=150、dup 0.215），d3 0.73:1 假象由此而来，剔除后 d3 亦缩短（2.13:1） |

实现：层 1 = `p2_sap_layer1_qc.py`，层 2 = `p2_sap_layer2_limma.R`（limma lmFit+eBayes），合并质检 = `p2_full_merge_qc.py`（depth 22/22 日志核验）。**结果**：显著并集 **977**（d1 292/d3 159/d7 435/d21 496/d60 392；≥2 时点 471——原 579 系证据脚本 `A and B or C` 优先级 bug 虚增，审计修正）；sham 重复相关 0.9257；d1 缩短:延长在高表达层更强（5.87:1）反深度伪迹；LOSO 稳定。Gate G2 组合门槛 (a) 层1 QC 通过 (b) 分层后 d1 缩短方向稳定+LOSO 敏感 (c) 事件级三条件成立 (d) 977 记录在案+PolyASite 距离与第二工具抽验——全过（用户批准"现在通过"）。P2-07（GSE225110 第二数据源交叉）因门槛全过而**未触发关闭**；GSE225110 已有试点资产可复用。
**第二工具链**：salmon 1.10.3 12 样本定量 + PolyASite 注释，全基因组 dLongFrac vs ΔPDUI Spearman 五时点全阳性（rho 0.041–0.139，p≤0.0023）；top40 单事件一致性 13/27 = 功效不足如实登记（⚠ 27 项不含 Atp2a2；口径=DaPars2 vs **Salmon**，QAPA 未构建成功）。

### 3.5 P3 关联层（M3 三级证据）

- **L1 文献池**：12 对 regulator–target（PABPN1×2/NUDT21×3/QKI/ELAVL1/PTBP1/NOVA×2/MBNL/TARDBP），PMID 逐条经 NCBI 核验（纠正一处错引 24814346→24814343）；FUS 记阴性。
- **L2 CLIP**：QKI peaks 取自 GSE147119 论文补充表（经 Europe PMC REST 绕开防爬），437 peaks（FDR<0.05 420），mm10 50-nt bin。表述纪律：只可写"剪短预测丢失一段 QKI 结合区"，不得表述为"该事件受 QKI 的调控作用"（原论文自身显示 QKI 靶基因 PAS 数目无差异）。
- **L3 motif**：ATtRACT 官站 403 → 经 zavolanlab GitHub 备份快照取库；16 个 regulator 共 325 motif 转 MEME 格式；FIMO 5.5.9 扫描；**阳性对照校准先行**（60 条真实 QKI CLIP peak 序列，1,285 hits，冻结后不回调）；扫描阈值 1e-3；**null 两轮修正**——random-genome null 判无效（任何 UTR vs 随机基因组=平凡 73 倍），改用"非显著事件的远端段、长度匹配、seed 42"的 UTR 匹配 null + 逐家族 Poisson 密度检验。结论：**无任何家族富集（阴性照登）**。

### 3.6 P4 定位层（M4 三轴）与分级

- **三轴构建**（审计后废除"两参考一致=高置信"投票）：Set1 = GSE74456 PAP vs Ctx TRAP，pydeseq2 ~condition，PAP 富集 = padj<0.05 & log2FC>0 & baseMean≥10 → **2,729 基因**（8 次失败后根因=计数矩阵未转置传给 pydeseq2，导致其为 41,393×41,393 单位矩阵申请 12.8GB 内存而报错、设计矩阵退化——勘误记录在案）；Set2 = GSE286075 padj<0.05 → 409 + 背景 14,320；Set3 = GSE263986 每 zone n=5+5 padj<0.05 并集（皮层 4,773/白质 10,434）。
- **检验**：Fisher 交集全阴性 → fgsea 10,000 次置换、4 轴×6 排序：**PAP 轴一致向缩短偏移**（d7 padj≈1e-4、d3 0.0018、omnibus 0.0042、d1 0.044；d21/d60 阴性；后续按用户追问如实降级为"弱、急性期、rank 级"）；终足轴阴性。
- **L4 保守性**：v1 = bedGraph 中间管线合并元件 ≥10bp 交叠 = 65%（**后被重建阶段判定为管线伪影**）；v2（见重建阶段）直接逐碱基统计。
- **分级（P4-03，看交集前冻结）**：L3 = fimo p<1e-4 且家族∈冻结 8 家族（QKI/ELAVL1/PTBP1/MBNL2/NOVA1/NOVA2/TARDBP/HNRNPA2B1）；L4 = phastCons 元件 ≥10bp；A 级因 MPRA 阻断不可达；B=L2&L3&L4，C=L3&L4，D=L3。判定 **B=5/C=475/D=250/none=247**；Top 候选 Ndrg2≈Agpat3>Cnp≈Pea15a>Fam107a（**注意：此排序基于错侧区间，后被重建阶段全面修订**）。
- Gate G3 用户确认通过（三轴分别表述、废除投票合并）；Gate G4 用户确认通过；P5 交付（终报告 md+html、路线书、红线 grep 12 条禁语 0 违规；重建阶段另用 9 条精简清单复检，两套口径各自全过）；Gate G5 通过，v5.1 干实验阶段关闭。

---

## 4. 转折点记录（v5.1 收官 → 重建启动）

1. **用户追问 PAP 信号强度** → 如实复核：fgsea NES≈1.2（弱）、仅 d1–d7、无二值富集、d21/d60 阴性——结论降级为"弱、急性期、rank 级偏好"，主证据回归事件层。
2. **用户确认问题边界** → "不能直接说富集在 PAP；变短前/变短后两问"；聚焦为"同基因长/短异构体在胞体 vs 突触旁的分配"2×2 设计，主案例 Atp2a2（三要素：缩短+QKI 重叠+PAP 注释）。
3. **外部论证（用户提供的 docx/截图）发现负链方向错误** → Atp2a2/Ndrg2/Agpat3 均负链，旧 lost_segments.bed 把丢失段定义在 [PAS→高坐标]（错侧）；人工校正：Atp2a2 校正后丢失段（预测断点口径）= chr5:122,453,512–122,456,498，含 QKI peak 122,454,200–250（FDR 1.04e-5）；旧引用 peak 122,456,550–650 在保留侧。用户定稿 Proposal+TODO（P0–P10），进入重建阶段。

---

## 5. 阶段二：链方向重建（R0–R4 + Gate 1，09-21；R0–R4 即 02_项目TODO.md 的 P0–P4）

### 5.1 R0 输入冻结

`p0_freeze_inputs.py`：对 `input_links/` 17+4 个输入逐个 sha256+字节数（`input_files_sha256.tsv`）；两个交付 zip 清单+校验和；软件版本存档；provenance 清单（来源路径+取得日期）。工作区骨架按 TODO §3 建立（config/input_links/metadata/scripts/results/{00–05}/wetlab/logs），git init，`.gitignore` 排除大文件（bam/salmon/references/tools/*.fa/*.gtf/*.bw）。

### 5.2 R1 链方向与坐标审计

**1bp 转换裁定（P1.1）——无需转换**。源码证据：DaPars2 输出 `Loci` 在 20% 修剪与 ±1 收缩**之前**由原始 BED 坐标构造（Load_Target_Wig_files ~L325）；`Predicted_Proximal_APA` = 修剪窗（搜索边距 search_point_start=150、search_point_end=5%·len）内的基因组坐标，负链降序遍历但输出仍为基因组坐标（De_Novo ~L255–278）。两者同为 BED 0-based 尺度 → `proximal_pas_raw = proximal_pas_bed`。人工校验：Atp2a2 PAS 122,456,498 ∈ 负链搜索区 [122,453,693, 122,457,002] ✓。结构限制如实登记：PAS 永远在修剪窗内（不可能贴搜索侧两端），分辨率 ±1nt，表述一律为"预测断点"。

**事件表构建**（`p1_build_event_coordinates.py`）：从 PDUI 矩阵（9,693 行）提取 event_id（DaPars2 名）、event_num（1-based 数据行号——SAP 表的整数 event 的来源，经 Atp2a2 行 6954−1=6953 验证）、Loci 解析 utr_start/utr_end、Predicted_Proximal_APA、链方向（genePred 按转录本 ID 联接，142,604 条）。结果 0 重复/0 缺链/0 解析失败。

**区间切分**（`01_build_strand_aware_segments.py`）：

```text
规则（config/coordinate_convention.yaml 冻结）：
  正链：shared=[utr_start,PAS)，distal=[PAS,utr_end)
  负链：distal=[utr_start,PAS)，shared=[PAS,utr_end)
校验：utr_start < PAS < utr_end（严格在内）；两段长度>0；
      和 = UTR 长度；无重叠；event_id 唯一；strand∈{+,-}
```

**单元测试**（`test_01_segments.py`）**17/17 通过**：正/负链合成事件（1000-2000，PAS=1500）方向互验；4 种边界拒绝（PAS 贴边/出界/链非法）；Atp2a2 锚点精确复现（distal=122,453,512-122,456,498；shared=122,456,498-122,457,153）；断言新 QKI peak ∈ distal、旧 peak ∈ shared；Ndrg2/Agpat3 负链低坐标侧+长度和恒等。
**全量**：9,693/9,693（+4,913/−4,780），0 invalid，五项标准逐行全过。

### 5.3 R2 注释重建

- **QKI 重算**（bedtools sort → `intersect -wa -wb` distal/shared × peaks；`p2_qki_summary.py` 汇总逐事件 distal peak 数/min FDR/重叠 bp/shared 计数 + L2v2 标记）：**三条完成标准全过**（Atp2a2 distal 命中 122,454,200–250 FDR 1.04e-5；旧 550–650 区域只在 shared；Ndrg2/Agpat3 distal 零命中）。（v1 错侧口径的对照基线：11/977、OR=2.71 p=0.013）全景 42/9,693（FDR<0.05）；977 显著事件中 **9 个** distal 有 QKI peak（FDR<0.05）——其中 8 个进入 class B，第 9 个是 **Qk 自身的事件**（event 3015，distal 最小 FDR 4.13e-4，因 L3v2=no 落 none；构成 Qk 自调控回环的观察线索）。
- **取序**（`p2_extract_distal_fastas.py`）：`.fai` 索引随机访问 genome.fa，负链反向互补（bedtools getfasta -s 等价语义）；候选集 = sig & dPDUI≤−0.1 的 (事件,contrast) 对 **1,541**（与旧 bed 行数精确一致，规则复原）→ 824 唯一事件；简单序列 ID（E####/N####）+映射表——防 v1 踩过的 FIMO `::` UCSC 解析坑。
- **null 构建**（`p2_build_null_set.py`）：非显著事件（9,693−977=8,716）distal 段为池；逐候选按 |ΔGC| 在长度 ±10% 窗内选最优（seed 42 洗牌，两档回退 ±25%/最近长度）；824/824 配对，长度差中位 21.5nt、GC 差 0.03pp。
- **FIMO 双扫描**：`--thresh 1e-4`（5.5.9 中按 p 值截断；tsv 末尾有 footer 行需跳过）；候选 61,488 命中全 p<1e-4；792/824 有冻结 8 家族强命中（**L3 饱和**）。
- **家族富集**（`p2_motif_summary.py`）：逐家族 Fisher 精确（候选 vs null 有命中段数）+ BH 校正 16 家族：**全部无富集**（最小 p=0.054 TARDBP，q≥0.71）——v1 阴性结论（v1 口径：730/977 有冻结家族命中但 vs 匹配 null 无富集）在正确坐标下复现。
- **保守性**（`p2_conservation.py`，bioinfo env pyBigWig 本地读 bw）：逐事件 distal/shared 逐碱基统计 mean score 与 score>0 碱基数；**L4v2 = score>0 ≥10bp**（与 v1"合并元件≥10bp"按碱基集合等价——合并不改变覆盖碱基集）→ 9,667/9,693（99.7%），无区分度；**勘误 v1 65% = bedGraph 中间管线伪影**。连续量（mean score，候选中位 0.299）更有信息量。
- **矩阵 v2**（`p2_rebuild_matrix.py`）：977 事件 × 校正后 L2v2/L3v2/L4v2 + class v2（同冻结规则）+ 排序 7 字段 + PAP 成员 + 第二方法 + 探针可行性。**class v2：B=8/C=783/D=1/none=185**；346 翻转；Atp2a2 D→B、Ndrg2/Agpat3 B→C。新旧差异三表齐备。
- **PAS 距离**（`p2_pas_distance.py`）：PolyASite 同链最近 cluster 二分查找；977 事件 <50bp 302/<200 386/<500 197/≥500 92；Atp2a2 **13bp**。

### 5.4 R3 六样本 BAM 核查

- **P3.1**：quickcheck 6/6 OK；idxstats 全染色体在位。
- **切片**：审查集 10 基因（提案三基因 Atp2a2/Ndrg2/Agpat3 + 重建 B 级 8 个；Atp2a2 两者重叠）UTR±5kb × 6 样本，120 个 cand.bam+bai（IGV 便利件）。
- **P3.3 深度**：`samtools depth -a -Q 20 -q 20 -b regions.bed`（MQ20+BQ20）；逐事件逐样本 distal/shared 平均深度、覆盖比、≥1/5/10 覆盖分数。
- **P3.2 图**（`p3_plots.py`）：10 基因固定尺度 PNG（公共 y 轴=六样本最大深度；PAS 虚线、distal 灰带、QKI 橙带）+ IGV 快照指令单。
- **富 A 检查**（`p3_arich_check.py`）：PAS ±50nt 基因组窗（负链取转录方向窗口），max A-run ≥6 或 A%≥60 判 HIGH——**10/10 low**（内部引物风险排除）。
- **P3.4 审查表**（`p3_review_table.py`）：机械规则 advance = d7 覆盖比降幅 ≥25% 且双重复一致且富 A low 且探针可行；置信度 high（降幅≥40% 且 d3+d7 全一致）/medium（d7 一致且深度足）/low（深度薄）。逐样本覆盖比（distal/shared 平均深度比）：**Atp2a2 sham 0.302 → d3 0.084 → d7 0.085（-71.7%，d3+d7 双重复全一致；shared 侧组均值深度 291→565 反升，排除总 RNA 下降假阳性）**；**Agpat3 0.471 → 0.201 → 0.083（-82.3% 全一致）**；Aplp1 0.299 → 0.159（-46.9%，⚠总 RNA 同步暴跌 4,823→242 需控制）；Sirt2 0.381→0.240（-36.9%）；Ndrg2 0.542→0.365（-32.7%）；Kazn/Plec/Fam107a 方向一致但深度薄（d7 总深 71–204，判 low 置信）；Cnp -19.3%、Pea15a -22.4% 未过 25% 阈。审查表结果 **8 advance（high 2：Atp2a2、Agpat3；medium 3：Aplp1、Sirt2、Ndrg2；low 3：Kazn、Plec、Fam107a）/ 2 hold / 0 drop**。Aplp1 附总 RNA 塌陷警示（4,823→242）。

### 5.5 R4 公共数据辅助 + Gate 1

- **GSE143531**（`p4_1_gse143531.py`）：48 样本=3 动物×Full/PAP×8 技术重复（series matrix 解析）；技术重复按 (动物,分区) **求和合并为 6 组**（48 条目≠48 只动物）；主靶与两备靶全部检出 ✓（审查集中 Kazn 2/3 动物不可检出，如实登记）；PAP/Full 全基因组背景（6,762 可用基因，中位 log2 = −3.42）→ **候选均无富集**（Atp2a2 54.0% 百分位；仅 Plec 95.8%）——v1 PAP 富集先验降级为区域依赖待验证。
- **GSE330741**：本地 651 基因库缺 Glt1/Sparc/Hsbp1 及全部候选（可见部分基因符号已被 Excel 转成日期，如 "9-Sep"）——manifest 矛盾未解，对照待 Koester 2026 supplementary；不阻断。
- **第二方法基线**（`p3b_salmon_baseline.py`）：salmon TPM 长异构体占比——Atp2a2（数值=重复1/重复2 的 TPM 长异构体占比）sham 0.719/0.541 → d1 0.059/0.105 → d3 0.130/0.112 → d7 0.039/0.039 → d21 0.066/0.059 → d60 0.275/0.127，与 PDUI 五时点逐一同向；Agpat3 同向；Aplp1 单转录本不可算（如实）。
- **Gate 1 定靶**（主靶 Atp2a2 = event 6953；重建阶段的靶标冻结关口，区别于 v5.1 的 Gate G1；`p4_target_nomination.py` → `target_nomination_v1.md`，用户批准）：主靶 **Atp2a2**（B 级；QKI FDR 1.04e-5；PAS 13bp；覆盖比 -71.7% 全一致；PAP+（皮层 GSE74456 注释成员，海马数据不富集，见 §5.5/§5.6）；富 A low；第二方法 ✓（salmon 长异构体占比五时点同向，d1 0.059/0.105、d21 0.066/0.059 亦同向））；备靶 **Agpat3**（event 216；读段最强 -82.3%；class C 如实；PAS 距已知 171bp；**未声称 PAP——两套数据核对均不富集，矩阵 PAP_member=no**）、**Aplp1**（event 7821；B 级自有 distal QKI peak FDR 0.02；PAS 182bp；⚠总 RNA 塌陷需控）；Sirt2 替补、**Ndrg2 保留观察（event 2072；PAS 距已知 8bp（审查集 10 基因内最近；977 全景另有 0–1bp 事件），但校正后 QKI 证据消失、d3 覆盖比不一致、健康海马百分位 13.9%）**。三靶 3'RACE 外/内引物区、ddPCR common/distal 扩增子区、smFISH 铺瓦区全部给出（几何校验：引物区∈shared、distal 扩增子∈distal）。


### 5.6 文献佐证（2026-09-21 补充，事后描述性、不参与冻结分级）

- **Atp2a2 正常状态的突触旁定位有发表数据直接支持**：Sakers 2017（PNAS，GSE74456，即本课题 PAP 注释来源）中 Atp2a2 突触旁富集 **log2FC = +2.13（约 4.4 倍）、padj = 8×10⁻⁷**、baseMean 8,288（表达充足）——该统计量为候选确定后应用户问询补充调取（该基因行的原始数据自 09-19 建表起即在库），不改变任何冻结规则（B 级仅依赖成员资格）。
- **两套独立 PAP 数据对比**：皮层（Sakers）强富集 vs 健康海马（Mazaré 2020）中位（54% 百分位）——区域依赖，湿实验应在皮层直接测量。终足数据集无 Atp2a2。
- **概念先例**：Grzejda et al. 2025（EMBO reports 26:1792–1815，PMID 39984683）在果蝇脑突触区域证明"APA 异构体差异定位 + RBP（Pumilio）差异结合"机制成立——哺乳动物星形胶质细胞+卒中场景仍为空白（检索无同类文献）。
- 三靶标+审查集两数据集并排统计：`results/04_public_pap/targets_two_dataset_pap_stats.tsv`。

---

## 6. 关键参数速查总表

| 环节 | 参数 | 值 | 为什么 |
|---|---|---|---|
| scRNA QC | genes/counts/mt | 200/500/20% | 10x 底线过滤 |
| scRNA 双联体 | Scrublet→MAD±3 降级 | — | 预案化降级，防环境问题卡死 |
| 聚类 | HVG2000/PCA50/leiden 1.0 | — | 中等粒度标准流程 |
| astro 定义 | Aqp4/Slc1a3/Gfap 双规则 | — | 覆盖稳态+反应性簇 |
| 比对 | fastp 线程 8；STAR 线程 8；sort -@4 | — | 12 核留余量 |
| 覆盖 | genomecov **-bg -split** | — | 防剪接读段污染外显子覆盖 |
| 深度 | view -c **-F 0x904** | — | 主比对口径，与 DaPars2 归一一致 |
| DaPars2 | Coverage_threshold=10 | 默认 | 低覆盖 PDUI 不稳 |
| DaPars2 | 注释=slim 最长 3'UTR | 21,158 | 用户批准；事件定义收窄 |
| SAP | logit ε=0.01 | 预冻结 | 防 0/1 边界发散 |
| SAP | eBayes robust+trend | — | 抗离群+均值方差趋势 |
| SAP | 联合 BH FDR<0.1；\|ΔPDUI\|≥0.1；4/4 有效；双重复同向 | 预冻结 | 最保守主口径 |
| fgsea | 10,000 置换；4 轴×6 排序 | — | 置换下限≈1e-4 |
| FIMO | 校准 1e-3→L3 收紧 p<1e-4 | 冻结 | 阳性对照校准后不回调 |
| L3 null | UTR 匹配（非显著事件远端段）seed 42 | — | 随机基因组 null 无效 |
| v2 null | 长度 ±10% 内取 \|ΔGC\| 最小，seed 42 | — | 长度+GC 双匹配 |
| 保守性 | score>0 碱基 ≥10 | 冻结 | 与 v1 元件定义碱基等价 |
| R3 深度 | -Q20 -q20 | TODO 规定 | 双 20 过滤低质量比对/碱基 |
| R3 富 A | ±50nt，A-run≥6 或 A%≥60 | 常规口径 | 内部引物风险 |
| R3 advance | d7 降幅≥25% 且一致且无风险且探针可行 | 机械规则 | Gate 1 前置 |
| PAS 桶 | <50/<200/<500/≥500 | — | 支持强度分级 |
| 定靶区间 | 引物区∈shared、distal 扩增子避开边界±50–100bp | 几何校验 | 防跨段设计 |

## 7. 诚实登记与勘误总表

1. v1 lost_segments 负链错侧 → QKI/motif/保守性/class 全部重建；旧字段仅存审计
2. v1 65% 保守率 = bedGraph 中间管线伪影（v2 直算 99.7%，L4 无区分度）
3. L3 饱和（792/824=96%）；class 区分度由 L2 承担
4. GSE143531 不支持候选 PAP 富集先验（Atp2a2 54% 中位）——降级为区域依赖待验证
5. GSE330741 manifest 矛盾未解；对照待 Koester 2026 supplementary
6. PolyASite 旧域名废弃（NXDOMAIN）→ 新址 unibas.ch
7. 13/27 口径=DaPars2 vs Salmon（QAPA 未成）；27 项不含 Atp2a2（已单独补算，五时点同向）
8. FIMO 5.5.9 --thresh 按 p 值截断；tsv 有 footer 行
9. DaPars2 PAS 必在 20%/150nt/5% 修剪窗内；一切 PAS="预测断点"
10. v5.1 曾发生：579→471 优先级 bug、GSE286075 链取列错误重跑、sham1_2 md5 损坏重拉、深度表占位值即改、**day3_rep1 文库异常（insert=150，影响 d3 解读）**——全部在案
11. Aplp1 单转录本→第二方法不可算；总 RNA 塌陷需湿实验控制
12. d1/d21/d60 事件未在 R3 核查窗口（BAM 可随时补切）

## 8. 交付物

见 `REANALYSIS_FINAL_REPORT.md` §九 + 本报告；打包件 `REANALYSIS_RESULTS_20260921.zip`（89.1MB/313 项，含全部表、10 张图、120 切片、GSE143531 原表、PAS atlas）；另有 v5.1 期交付包 `D:\stroke_apa\STROKE_APA_FINAL_DELIVERABLES_20260920.zip`（25,387,318 B）与 IGV 读段核查包 `ATP2A2_IGV_CHECK_PACKAGE_20260920.zip`（29,569,120 B）。

## 9. 校验记录（两层）

### 9.1 本报告（全程报告）的五轮自校验

| 轮 | 方法 | 结果与修正 |
|---|---|---|
| R1 | 报告全部数字/参数逐一对账实际文件（70+ 项） | 首轮 13 处数字缺失（备靶覆盖比、PAS 距离、事件号、背景中位、Sakers 佐证、zip 尺寸等）——全部补入；4 处为检查串写法误报 |
| R2 | 脚本/章节完整性（26 个核心脚本点名 + 章节 1–10） | 首轮 5 个脚本名未点名——补入后全过 |
| R3 | 参数实证：报告声明 vs 实际脚本/冻结文件 grep | 21/21 过（fastp/STAR/-split/0x904/Coverage_threshold=10/ε=0.01/robust+trend/FDR<0.1/|ΔPDUI|≥0.1/4/4/seed42/-Q20-q20/A-run≥6/25% 规则等全部实证） |
| R4 | 红线 grep + 纪律标记 | 首轮 1 处命中为纪律声明引用的禁语字面——改写措辞后零命中 |
| R5 | 新眼通读全文 | 3 处修正：文首节号错引、§9 两层校验混淆、§5.4 审查集构成表述澄清 |

### 9.2 重建报告（REANALYSIS_FINAL_REPORT.md）的五遍审计（09-21 已执行）

| 轮 | 方法 | 结果 |
|---|---|---|
| 1 | 全链数字独立复算（34 检：从 depth/quant.sf/fimo/atlas/交集原始数据重推） | 全过（1 处审计脚本自身路径 bug，修正后全绿） |
| 2 | TODO §16/§18 交付物逐项对账 | 26/26 过 |
| 3 | 跨文书一致 + 设计区间几何（引物区∈shared 等） | 22/22 过（1 处圆整容许度修正） |
| 4 | 纪律红线 grep（9 条禁语）+ 诚实标记 + 旧值勘误语境 | 干净 |
| 5 | 报告逐数字回验（41 项）+ zip 直验 | 全过（2 处检查写法修正） |

两层审计发现的问题全部是**审计/校验脚本自身**的缺陷，无数据与结论错误。

## 10. 核心代码详解（真实片段 + 逐段说明）

> 以下为各阶段承载关键逻辑的代码单元。完整脚本均在 `scripts/`（v5.1 在旧区 scripts/），此处摘录核心段。

### 10.1 比对链（`p2_gse238125_align.sh`，bash）

```bash
fastp -i R1 -I R2 -o trimmed_1 -O trimmed_2       --thread 8 --detect_adapter_for_pe
```
`--detect_adapter_for_pe`：双端测序接头序列未知时由 fastp 从数据自动推断；不开此选项，接头残留会伪增 3' 端覆盖，直接扭曲 APA 覆盖剖面。
```bash
STAR --runThreadN 8 --genomeDir <index> --readFilesCommand zcat      --outSAMtype BAM Unsorted --quantMode GeneCounts      --outTmpDir /home/taylor/STARtmp/${name}_STARtmp
```
`--outTmpDir` 指到 ext4：STAR 临时文件用 FIFO，NTFS（/mnt/d）不能建 FIFO，放在 /mnt 下会直接报错。`Unsorted` 而非 SortedByCoordinate：排序交给后续 `samtools sort -@ 4`，STAR 内部省一次大内存排序。`--quantMode GeneCounts` 顺带产出基因计数，供 QC 用（非主分析）。
```bash
bedtools genomecov -ibam X.sorted.bam -bg -split
samtools view -c -F 0x904 X.sorted.bam
```
`-split`：把 N 操作（剪接）两侧分开计覆盖——不加则内含子被当作覆盖，3'UTR 覆盖剖面系统性失真（试点期 reanalysis_v3 纠错的核心项）。`-F 0x904`：0x4 未比对 + 0x100 次要 + 0x800 补充 = 只数主比对，与 DaPars2 深度归一口径一致。

### 10.2 SAP 层 2 统计（`p2_sap_layer2_limma.R`，R）

```r
pdui_t <- (pdui * (1 - 2*eps) + eps)          # eps = 0.01
logit  <- log(pdui_t / (1 - pdui_t))
fit <- lmFit(logit_mat, design)               # design = ~0 + timepoint
fit <- eBayes(fit, robust = TRUE, trend = TRUE)
```
`pdui*(1-2eps)+eps`：把 [0,1] 映射到 [0.01,0.99]，logit 不发散；不这样做，PDUI=0/1 的事件直接产生 Inf。`~0+timepoint` 无截距参数化：每组系数=该组均值，contrast 直接是组间差。`robust=TRUE`：对少数方差极端的事件做经验贝叶斯权重收缩，防个别噪声事件绑架全局方差估计；`trend=TRUE`：PDUI 估计方差随覆盖下降，需均值-方差趋势。contrast 与联合 BH：
```r
contrasts <- makeContrasts(d1-sham, d3-sham, d7-sham, d21-sham, d60-sham, levels=design)
fit2 <- eBayes(contrasts.fit(fit, contrasts))
# 联合 BH：全部 事件×5 contrast 的 p 值合并一起调 p.adjust(method="BH")
```
事件资格过滤在进模型前完成（4/4 有效 PDUI，NA 不插补）；显著判定 = 联合 padj<0.1 且 |ΔPDUI|≥0.1 且 rep1、rep2 的 ΔPDUI 同号。

### 10.3 链方向区间切分（`01_build_strand_aware_segments.py`）

```python
if strand not in ("+", "-"): raise ValueError(...)
if not (us < pas < ue):      raise ValueError("PAS not strictly inside UTR")
if strand == "-":
    distal = (us, pas); shared = (pas, ue)
else:
    shared = (us, pas); distal = (pas, ue)
if not (0 < distal[1]-distal[0] and 0 < shared[1]-shared[0]): raise ...
if (distal[1]-distal[0]) + (shared[1]-shared[0]) != ue - us:  raise ...
```
四个 raise 就是完成标准的代码化：链域合法、PAS 严格在内、段长为正、和守恒。配置从 yaml 读入而非硬编码——冻结规则变更时只改配置。`--coordinate-config` 存在的意义：让"规则"与"实现"分离，审计时先审配置再审代码。

### 10.4 坐标表构建（`p1_build_event_coordinates.py`）

关键点一：**event_num 的来源**——PDUI 矩阵数据行的 1-based 序号，即 SAP 表整数 event 的出处（Atp2a2=6953 由矩阵第 6954−1 行验证）。关键点二：链方向联接键是**转录本 ID**（genePred 142,604 行 `name -> strand`），匹配率 100%——若用基因符号会因一对多出错。关键点三：proximal_pas_raw 与 proximal_pas_bed 两列并存（P1.1 裁定两者相等），保留字段是为了未来 DaPars2 版本变化时可追溯。

### 10.5 链方向取序（`p2_extract_distal_fastas.py`）

```python
def fetch(f, fai, chrom, start, end):        # .fai 随机访问，不整文件读入
    ln, off, lb, lw = fai[chrom]
    f.seek(off + (start // lb) * lw)          # 行首偏移 = off + 完整行数*行宽(含换行)
    ...
if strand == "-":
    seq = seq.translate(COMP)[::-1]           # 反向互补 = bedtools getfasta -s
```
2.6GB 基因组不能整读，用 .fai 的 (offset, linebases, linewidth) 做逐行 seek。`translate(COMP)[::-1]` 一行完成碱基互补+倒序。头部用简单 ID（E0001/N0001）而非原名加坐标：v1 踩过的坑——bedtools getfasta 生成的 `name::chr:start-end` 头会被 FIMO 当 UCSC 区域解析，导致序列名塌缩。

### 10.6 匹配 null 构建（`p2_build_null_set.py`）

```text
rng = random.Random(42)   # 伪代码（节选）                        # 固定种子：任何人重跑得到同一 null
rng.shuffle(pool)
def pick(c):
    for win, key in ((0.10, |dGC|), (0.25, |dGC|), (1e9, |dLen|)):   # 三档回退
        feas = [p for p in pool if not used and |len差| <= max(win*len, 20)]
        if feas: best = min(feas, key=key); best["_used"]=True; return best
```
设计意图：长度最紧的窗内优先选 GC 最接近的；窗内无解则放宽窗。seed 42 + 贪心确定序 → null 可精确重现（论文审稿可复现性）。结果：824/824，长度差中位 21.5nt、GC 差 0.03pp——null 与候选在两个混杂维度上几乎重合，富集检验的阴性结论因此可信。

### 10.7 保守性逐碱基统计（`p2_conservation.py`）

```python
iv = bw.intervals(chrom, s, e)                 # pyBigWig 取该区间全部连续块
for b0, b1, v in iv:
    ov = min(e, b1) - max(s, b0)               # 块与段的交集
    tot += ov;  pos += ov if v > 0 else 0
L4v2 = pos >= 10
```
逐块求交而非逐碱基循环——大区间也只遍历几十个块。`score>0 碱基数≥10` 与 v1"合并保守元件交叠≥10bp"在碱基集合上恒等（相邻 score>0 块合并与否不改变覆盖碱基），因此与冻结规则兼容，却完全绕开了 v1 的 bedGraph 中间管线（该管线正是 65% 伪影的来源）。

### 10.8 覆盖比与判定（`p3_coverage_stats.py` / `p3_review_table.py`）

```python
vals = [depth[s].get(p, 0) for p in range(a, b)]     # 缺位=0：均匀分母
ratio = mean(distal) / mean(shared)                  # 描述性读段证据
advance = (drop7 >= 0.25) and conc7 and arich=="low" and probe=="ok"
```
判据全部机械化并预先写死在代码里——advance/hold/drop 不由人临场判断。`confidence` 三档（high：降幅≥40% 且 d3+d7 全一致；medium：d7 一致且深度足；low：深度薄）把"方向对但证据薄"与"证据扎实"分开，防低深度事件混入高置信行列。

### 10.9 GSE143531 技术重复合并（`p4_1_gse143531.py`）

```python
groups.setdefault((animal, part), []).append(gsm)
agg[gene][key] = sum(samples.get(gsm, 0) for gsm in groups[key])
```
48 条目按 (动物, 分区) 求和合并为 6 组——**计数数据求和技术重复是合法的**（恢复单动物的分子总数），前提是分析以动物为单位；若做线性建模则应保留重复层级。百分位检验：把每个候选的 mean log2(PAP/Full) 放进全基因背景分布里排位（6,762 个可用基因），避免"负值=不富集"的误读。

### 10.10 审计脚本（`audit_pass1_numbers.py`）

原则：**从原始数据重推，不读自己的汇总表**。例如覆盖比不读 `bam_coverage_review.tsv` 而是从 6 个 depth 文件重新解析计算；sha256 抽验直接重算文件哈希比对 P0 冻结值。任何一处不一致都会让脚本以非零退出——这是"五遍校验"第一遍的机械保障。
