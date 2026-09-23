# Atp2a2 等候选位点读段核查包（IGV package）

> 生成：2026-09-20 ｜ 课题：卒中→星形胶质细胞 APA（pipeline v5.1 干实验已关闭，G1–G5 全过）
> 用途：湿实验聚焦前的最后一道干实验核查——在**原始比对读段**层面确认 Atp2a2（主靶）与 Ndrg2、Agpat3（备靶）的 APA 事件真实可靠。**不需要全量重跑。**

---

## 一、背景（为什么要做这次核查）

现有结果表只能告诉我们"算法报了什么"（PDUI/ΔPDUI/padj），不能单独判断覆盖变化是否由比对伪影、末端外显子混淆或单只样本造成。已知的两个事实使这次核查更必要：

1. 第二工具交叉验证的 **13/27 同向**来自 DaPars2 vs **Salmon 定量推算的长异构体比例**（QAPA 构建未成功，不得写作 "DaPars2/QAPA"）；且那 27 项里**不含 Atp2a2**。
2. Atp2a2 是湿实验拟主靶（见下表），其证据主张只在事件层（显著缩短 + 丢失段 QKI 重叠 + PAP 注释；class = D，丢失段无保守元件，如实登记）。

## 二、三个候选（本次切片对象）

| 基因 | event | 3'UTR 注释区间（DaPars2 bed） | 显著时点 | 最强 ΔPDUI（padj_min） | class | 逐样本方向 |
|---|---|---|---|---|---|---|
| **Atp2a2** | 6953 | chr5:122,453,512-122,457,153 | d3/d7 | −0.125（0.028） | D | d3 −0.11/−0.14；d7 −0.10/−0.11（一致） |
| Ndrg2 | 2072 | chr14:51,905,270-51,906,142 | d1/d7/d21/d60 | −0.280（0.034） | B | 五时点 4 显著，重复全部同向 |
| Agpat3 | 216 | chr10:78,269,177-78,273,689 | d1/d3/d7/d21 | −0.305（0.074） | B | 四时点显著，重复全部同向 |

Atp2a2 的 PDUI 逐样本读数：sham 0.19/0.20 → d3 0.08/0.06 → d7 0.09/0.09 → d21 0.11/0.13 → d60 0.12/0.12（急性窗口缩短后部分回弹）。
完整数值见 `01_computed_results/candidate_focus/`（逐样本 PDUI、逐 contrast 统计、注释矩阵行）。

## 三、包内清单

```
01_computed_results/          已算出的结果（直接引用，无需重算）
  p2_full_pdui_matrix.tsv               9,693 事件 × 12 样本 PDUI 全矩阵
  p2_sap_layer2_events.tsv              SAP 层2 逐事件逐 contrast 统计（limma robust+trend）
  p2_sap_layer2_omnibus.tsv             omnibus 时间效应
  p2_second_tool_top40_signs.tsv        第二工具（Salmon 长异构体比例）top40 逐事件符号
  p2_second_tool_consistency.tsv        全基因组 dLongFrac vs ΔPDUI Spearman（五时点全阳性）
  p3_p4_event_annotation_matrix.tsv     977 显著事件 × L1/L2/L3/L4 + 三轴 + class
  p4_candidate_cards.md                 5 张 B 级候选卡（坐标/motif/保守/QKI）
  p4_top_b_class_candidates.tsv         B 级候选清单
  final_answer_gene_list.tsv            977 基因最终答案表
  candidate_focus/                      ★ 三候选聚焦表（本包核心速查）
  qc/                                   样本级 QC（相关/PCA/分层/合并 QC）
  params/                               比对与 DaPars2 运行脚本（参数记录）
  P2_SAP_frozen.md                      预注册统计计划（冻结于全量统计前）
  data_inventory.tsv                    全部产物的登记总表

02_bam_slices/                读段证据（IGV 直接拖入 .bam 即可，.bai 已附）
  regions.bed / gene_spans.tsv          切片区间定义 = 全基因 ±5kb
  <sample>.cand3.bam + .bai             6 样本：sham1 sham2 day3_rep1 day3_rep2 day7_rep1 day7_rep2
  qc_*_STAR_Log.final.out               STAR 比对 QC（5/6 样本；sham1 日志未存档，数值见 01/qc/）
  qc_*_fastp.log                        fastp QC
  dapars2_full_chr2_job.cfg             全量 DaPars2 实际运行 cfg（per-chr 同构）

03_sample_map_and_refs/       样本对应与参考
  p2_sample_map.tsv                     12 样本 ↔ SRR 运行号
  sequencing_depth_dapars2.tsv          DaPars2 深度表（口径 -F 0x904 主比对）
  gencode_M25_3UTR_for_DaPars2.bed      APA 注释（DaPars2 实际输入，tandem 3'UTR）
  p3_level2_qki_clip_peaks.bed          QKI CLIP peaks（437，可作 IGV track 叠加）
  p4_mpra_counts_long.tsv               GSE330741 MPRA 长表（reporter 对照设计参考；A 级不可用已登记）
  ENVIRONMENT.md                        完整环境与版本清单
```

## 四、在 IGV 里按什么顺序看

1. **Atp2a2（chr5:122,448,513-122,507,225）**：day7/day3 两只 vs sham 两只，看 3' 末端覆盖分布是否整体向预测近端位点移动（event 6953 断点在 3'UTR 注释区间 chr5:122,453,512-122,457,153 内；sham PDUI≈0.19-0.20 → d7≈0.09，即卒中后约 3/4 读段止于近端）。
2. **剪接与比对质量**：末端外显子是否有异常剪接/误比对；DaPars2 事件定义为"最长注释 3'UTR 内断点移位"，基因其他区域的异常与 APA 无关。
3. **单样本驱动检查**：两只重复的 3' 端形态是否同向（SAP 表 d_rep1/d_rep2 列与上面逐样本 PDUI 可直接对照）。
4. **备靶同法**：Ndrg2（chr14:51,900,271-51,919,158）、Agpat3（chr10:78,264,178-78,357,489）。
5. 可选叠加 `p3_level2_qki_clip_peaks.bed`：Atp2a2 预测丢失段的 QKI peak 位于 chr5:122,456,550–650（FDR 0.006）；Ndrg2 在 chr14:51,905,900-51,906,000。

**判定口径（建议）**：三候选全部通过（3' 端移位真实、无伪影、重复同向）→ 定 1 主靶 + 1–2 备靶进入湿实验；若 Atp2a2 读段层面不成立 → 主靶降级为 Ndrg2/Agpat3，如实登记。本核查不产生"卒中后 RNA 去了哪里"的结果——那是湿实验（3'RACE → 异构体定量 → smFISH/分区）的问题。

## 五、版本与可溯源性

- 物种/参考：**GRCm38 (mm10)**，**GENCODE vM25**（GTF = WSL `~/reference_v2/gencode.vM25.annotation.gtf`）
- 比对链：fastp 0.23.4（PE，--detect_adapter_for_pe）→ STAR 2.7.11b（8 线程）→ samtools 1.19.2 sort/index；bedgraph 语义 = bedtools genomecov **-split**
- APA 定量：DaPars2（conda env dapars2 / Py 3.8.20；含 2026-09-18 空覆盖补丁——空覆盖事件返回 NA，统计语义不变）；Coverage_threshold=10
- 统计：SAP v1.1（预注册，commit bd13e58/6579a0b）＝层1 样本级 QC + 层2 limma eBayes(robust=TRUE, trend=TRUE)，联合 BH<0.1、|ΔPDUI|≥0.1、双重复同向
- 原始数据：GSE238125 24 个 FASTQ（24/24 md5 通过）保留于 `D:\stroke_apa_data\fastq\`；12 样本全量 BAM 于 WSL `~/stroke_APA_work/bam/`；bedGraph/depth 冷备于 `D:\stroke_apa_data\ext4_backup_20260918\`
- **如需其余样本或其他基因区域**：用同法 `samtools view -b -L regions.bed` 切片即可，无需重新比对

## 六、可直接转发的一句话

> 原始数据与全量 BAM 均在。本包已含 GSE238125 sham×2、day3×2、day7×2 的 Atp2a2、Ndrg2、Agpat3 全基因±5kb BAM 切片（含索引）、样本对应表、深度表、APA 注释、QKI peak track 与全部参数版本。请 IGV 核查：3' 端覆盖移位是否真实、剪接/比对是否正常、两重复是否同向；暂不做全量重分析。
