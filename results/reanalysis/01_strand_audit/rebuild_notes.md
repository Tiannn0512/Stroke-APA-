# P1 链方向与坐标体系审计 — 完成记录（2026-09-21）

## P1.1 DaPars2 breakpoint 1bp 转换裁定：**无需转换**

证据链（DaPars2_Multi_Sample_Multi_Chr.py，补丁版 src 已迁 `data/tools/DaPars2/src/`）：

1. **Loci = 原始 BED 坐标（0-based half-open）**：`UTR_pos = "%s:%s-%s" % (chr, region_start, region_end)` 在 20% 修剪与 ±1 收缩**之前**构造（Load_Target_Wig_files, ~L325）。
2. **Predicted_Proximal_APA = 分析窗内基因组坐标**：搜索区在 20% 端修剪 + `start+=1 / end-=1` 后的窗内取点（± strand 方向相反，负链降序遍历），`selected_break_point` 直接输出为基因组坐标（De_Novo_3UTR_Coverage_estimation, L255-258, 275-278）。搜索边距 = `search_point_start=150`、`search_point_end=5%·len`（L235-236）。
3. 两者同为 BED 尺度 → **proximal_pas_raw = proximal_pas_bed**，事件表两列均携带。
4. 人工事件校验：Atp2a2 PAS 122,456,498 ∈ 负链搜索区 [122,453,693, 122,457,002] ✓。

**注意事项（写入惯例）**：PAS 永远落在 20%/150nt/5% 修剪窗内，不可能出现在搜索侧两端；坐标分辨率 = 覆盖剖面拟合的 ±1nt 点，一律表述为"预测断点"。

## P1.2/P1.3 执行结果

- 输入表：`input_links/dapars2_event_coordinates.tsv`（9,693 事件；strand 来自 GENCODE vM25 genePred 按转录本 ID 联接；0 重复 / 0 缺链 / 0 解析失败）。
- 单元测试：`scripts/test_01_segments.py` **17/17 通过**（正/负链合成事件、4 种边界拒绝、Atp2a2 锚点精确复现、Atp2a2 新旧 QKI peak 侧别断言、Ndrg2/Agpat3 负链低坐标侧 + 长度和恒等）。
- 全量：9,693/9,693 成功（+ 4,913 / − 4,780），0 invalid；完成标准五项全过（正长度 / 相邻无重叠 / 长度和恒等 / event_id 唯一 / strand ∈ {+,-}）。

## 三候选校正后区间（P2 重算 QKI/motif/保守性的输入）

| 基因 | UTR (BED) | proximal_pas (BED) | distal（负链=低坐标侧） | shared |
|---|---|---|---|---|
| Atp2a2 (6953) | chr5:122453512-122457153 | 122456498 | 122453512-122456498 | 122456498-122457153 |
| Ndrg2 (2072) | chr14:51905270-51906142 | 51905641 | 51905270-51905641 | 51905641-51906142 |
| Agpat3 (216) | chr10:78269177-78273689 | 78271713 | 78269177-78271713 | 78271713-78273689 |

侧别快查（正式重叠待 P2.1 bedtools）：Atp2a2 QKI peak 122,454,200-250 ⊂ distal ✓，旧引用 122,456,550-650 ⊂ shared ✓；**Ndrg2 旧 QKI peak 51,905,900-51,906,000 ⊂ shared**（distal 止于 51,905,641）——与 Proposal §2.2"Ndrg2/Agpat3 旧 QKI 证据校正后不成立"一致。

## P2 重建结果补记（2026-09-21）

- **QKI（P2.1）**：distal 有 QKI peak（FDR<0.05）的事件 = 42/9,693（977 显著集中 9 个：8 个 class B + Qk 自身事件 event 3015 落 none）；三条完成标准全过（Atp2a2 distal 命中 122,454,200-250 FDR 1.04e-05；旧引用 550-650 区域确在 shared；Ndrg2/Agpat3 distal 零命中）。v1 的 11/977 系错侧区间产物。
- **Motif（P2.3）**：FIMO 5.5.9 `--thresh 1e-4` 按 p 值截断；824 候选 distal 中 792（96%）有冻结 8 家族强命中——**L3 接近饱和**。vs 长度/GC 匹配 null（824/824，长度差中位 21.5nt、GC 差 0.03pp）：**16 家族全部无富集**（最小 p=0.054，全部 q≥0.71）。v1 的 motif 阴性结论在正确坐标下成立。
- **保守性（P2.4）**：直接从官方 bw 按 pyBigWig 逐碱基统计（score>0 碱基 ≥10 = 冻结 L4 二值）：9,667/9,693（99.7%）为 yes——**L4 在此阈值下无区分度**；v1 的 65% 判定为 bedGraph 中间管线伪影。均值得分（连续量）更有信息量：候选 distal 中位均值 0.299。
- **Class v2**：B=8 / C=783 / D=1 / none=185（v1 错侧口径：B=5/C=475/D=250/none=247）；346/977 class 变更、L2 变更 8。**B 由 L2 驱动**（L3/L4 饱和）。
- **Top B 候选**：Sirt2(4时点)、Aplp1(3,PAP+)、Cnp(2)、Pea15a(2)、**Atp2a2(2,PAP+，D→B)**、Kazn(1,PAP+)、Plec(1,PAP+)、Fam107a(1)。**Ndrg2/Agpat3 B→C**（旧 QKI 证据消失）——与 Proposal §2.2 预测一致。
- **待获取**：PolyASite mm10 PAS 表（站点网络不可达），PAS 距离排名字段暂 NA，不阻断 Gate 1。
- **PAS 距离（补全）**：旧域名 polyasite.ethz.ch 已废弃（权威 NXDOMAIN），新址 polyasite.unibas.ch；下载 GRCm38.96 atlas（301,006 clusters）。977 显著事件同链最近 PAS 距离：<50bp 302 / <200 386 / <500 197 / ≥500 92；**Atp2a2 预测 PAS 距已知 cluster 仅 13bp（强支持）**，Sirt2 159 / Aplp1 182 / Cnp 203 / Pea15a 130。

## P3 BAM/IGV 事件核查结果（2026-09-21）

- **P3.1**：6/6 BAM quickcheck 通过；idxstats 全 59 染色体在位（BAI 时间戳警告为 9p 拷贝所致，内容已 md5 抽验一致）。
- **P3.3 覆盖比（-Q20 -q20，distal_mean/shared_mean）**：
  - **Atp2a2：sham 0.302 → d3 0.084 → d7 0.085（-71.7%，d3+d7 双重复全一致）= 10 事件中最强读段证据**
  - Agpat3：0.471 → 0.201 → 0.083（-82.3%，全一致）——class C 但读段证据极强
  - Ndrg2：0.542 → (d3 杂) → 0.365（-32.7%）；Sirt2 0.381→0.24（-36.9%）；Aplp1 0.299→0.159（-46.9%，⚠ 总 RNA 同步暴跌 4823→242，比例解释需湿实验控制总 RNA）；Kazn/Plec/Fam107a 方向一致但深度薄（d7 总深 71–204）= 低置信
  - Cnp（-19.3%）、Pea15a（-22.4%）未过 25% 阈 → hold
- **P3.2**：10 基因固定尺度覆盖图（plots/*.png，公共 y 轴、PAS/边界/QKI 标注）+ igv_snapshot_commands.txt（用户 GUI 快照用；IGV 亦可直接加载 data/bam 全量 BAM）
- **富 A 风险**：10/10 low（max A-run ≤5，A% ≤34%）——内部引物风险全部排除
- **PAS 支持**：Ndrg2 8bp（审查集内最近；977 全景另有 0–1bp 事件）、Atp2a2 13bp、Pea15a 130、Kazn 125、Sirt2 159、Aplp1 182、Cnp 203、Plec 231、Fam107a 235、Agpat3 171
- **P3.4 审查表**（event_review_table.tsv）：advance 8（high：Atp2a2、Agpat3；medium：Sirt2、Aplp1、Ndrg2；low 置信：Kazn、Plec、Fam107a）/ hold 2（Cnp、Pea15a）/ drop 0。**P3 完成标准达成：Atp2a2 = advance 且证据链完整（读段+PAS+QKI+探针可行+无伪影）。**
- 主靶建议（供 Gate 1 用户决定）：**Atp2a2 主靶**（class B + QKI + PAS 13bp + 71.7% 一致下降 + PAP+），备靶 **Agpat3**（读段最强、class C 如实）与 **Sirt2 或 Aplp1**（B 级新候选）；Ndrg2 保留观察（PAS 8bp 但 QKI 证据消失、d3 不一致）。
