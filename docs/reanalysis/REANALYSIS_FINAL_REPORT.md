# 重建阶段最终报告（REANALYSIS FINAL REPORT）

> 生成：2026-09-21 ｜ 审计：五遍（数字独立复算 / 交付物对账 / 跨文书一致+几何校验 / 纪律红线 grep / 报告逐数回验）
> 工作区：`D:\stroke_apa_reanalysis\`（git 库）｜ 治理文件：`01_研究Proposal.md` + `02_项目TODO.md`（用户 2026-09-21 定稿）
> 决策记录：Gate 1 靶标由用户 2026-09-21 批准（"按你推荐的拍板，批准执行"）

---

## 一、任务缘起与执行总览

用户方在复盘中发现：**Atp2a2/Ndrg2/Agpat3 均位于负链，v5.1 干实验阶段的 `lost_segments.bed` 把丢失段定义在 [近端位点→高坐标] 一侧——对负链基因这是错侧**。该错误波及 lost 段序列、QKI/RBP 重叠、motif、保守性与候选分级（不波及逐样本 PDUI/ΔPDUI/时间轨迹）。据此用户定稿了三关验证方案（Proposal）与 P0–P10 执行清单（TODO），本轮完成其中的全部计算阶段：

| 阶段 | 内容 | 状态 | 关键产出 |
|---|---|---|---|
| P0 | 输入冻结与清单 | ✅ | sha256 清单 21 项、双 zip 清单+校验、软件版本、provenance |
| P1 | 链方向与坐标审计 | ✅ | 1bp 裁定（无需转换）、9,693 事件区间表、单元测试 17/17 |
| P2 | QKI/motif/保守性/分级重建 | ✅ | 三条完成标准全过；motif 阴性复现；L3/L4 饱和如实登记；class v2 |
| P3 | 六样本 BAM 核查 | ✅ | 覆盖比全表、10 张固定尺度图、富 A 全 low、审查表 8 advance/2 hold |
| P4 | 公共 PAP/MPRA 辅助 | ✅ | GSE143531 百分位分析、GSE330741 矛盾如实登记 |
| **Gate 1** | **定靶** | **✅ 用户批准** | **主靶 Atp2a2；备靶 Agpat3、Aplp1；Sirt2 替补；Ndrg2 观察** |

数据底座：12 样本全量 BAM+BAI、GRCm38 genome.fa、GENCODE vM25 GTF/genePred、PolyASite 2.0 atlas（301,006 clusters）、salmon 12 样本定量、补丁版 DaPars2 源码——**全部在 `D:\stroke_apa_reanalysis\data\`（已自 WSL 迁出并校验）**。原始 FASTQ 已于 09-20 删除（公开数据；重下校验元数据保留于 `D:\stroke_apa_data\`）。

## 二、P1 链方向与坐标审计

**1bp 转换裁定：无需转换。** 源码证据（`DaPars2_Multi_Sample_Multi_Chr.py`）：`Loci` 在 20% 修剪与 ±1 收缩**之前**由原始 BED 坐标构造；`Predicted_Proximal_APA` = 修剪分析窗内的基因组坐标（搜索边距 150nt/5%·len，负链降序遍历但输出仍为基因组坐标）。两者同为 BED 0-based 尺度 → `proximal_pas_raw = proximal_pas_bed`。人工事件校验：Atp2a2 PAS 122,456,498 ∈ 负链搜索区 [122,453,693, 122,457,002] ✓。注意事项：PAS 永远位于修剪窗内（不可能落在搜索侧两端），分辨率 ±1nt，一律表述为"**预测断点**"。

**事件表与区间**：9,693 事件全部拿到链方向（genePred 按转录本 ID 联接；0 重复/0 缺链/0 解析失败）。`01_build_strand_aware_segments.py` 按冻结规则切分 distal/shared；单元测试 **17/17 通过**（含正/负链合成事件、四种边界拒绝、Atp2a2 锚点精确复现、新旧 QKI peak 侧别断言）。全量 **9,693/9,693 成功（+4,913/−4,780），0 invalid**，五项完成标准（正长度/相邻无重叠/长度和恒等/ID 唯一/链域合法）逐行全过。

## 三、P2 注释与候选分级重建

### QKI 重叠（P2.1）——三条完成标准全过
| 标准 | 结果 |
|---|---|
| Atp2a2 distal 命中已知 peak | ✅ chr5:122,454,200–250，FDR 1.04e-05 |
| 旧引用 122,456,550–650 只在 shared | ✅（实为两个 50nt bin，FDR 0.0055/0.0709） |
| Ndrg2/Agpat3 distal 零命中 | ✅ 旧 QKI 证据校正后消失 |

全景：**42/9,693 事件** distal 含 QKI peak（FDR<0.05）；977 显著事件中 **9 个** distal 有 peak——8 个进入 class B，第 9 个为 Qk 自身事件（event 3015，L3v2=no 落 none；Qk 自调控回环观察线索）。v1 的"11/977"系错侧区间产物，作废。

### Motif（P2.3）——阴性结论在正确坐标下复现
候选集规则精确复原：sig & dPDUI≤−0.1 的 (事件,时间点) 对 = **1,541**，与旧 lost_segments.bed 行数完全一致；去重 **824 唯一事件**取链方向 distal 序列。匹配 null：824/824（长度差中位 21.5nt、GC 差 0.03pp，seed 42）。FIMO 5.5.9（`--thresh 1e-4` 为 p 值截断）：候选 61,488 命中全部 p<1e-4，**792/824 有冻结 8 家族强命中——L3 接近饱和**。vs 匹配 null 的 Fisher+BH：**16 家族全部无富集**（最小 p=0.054 TARDBP，全部 q≥0.71）。技术注意：FIMO 5.5.9 的 tsv 末尾有 footer 行需跳过。

### 保守性（P2.4）——含一条 v1 勘误
官方 mm10.60way.phastCons.bw（并行分块重下 4.58GB）+ pyBigWig 逐碱基统计：**L4（score>0 碱基 ≥10）= 9,667/9,693（99.7%）**——该阈值在 UTR 上**无区分度**（候选 distal 中位均值得分 0.299，连续量更有信息量）。**勘误：v1 的 65%（1,004/1,541）判定为 bedGraph 中间管线伪影**，非生物学差异。

### Class v2 与新旧对照
按冻结规则（B=L2&L3&L4；C=L3&L4；D=L3；A 不可达）重判 977 事件：**B=8 / C=783 / D=1 / none=185**（v1 错侧口径 B=5/C=475/D=250/none=247）。346 个事件 class 翻转；L2 状态变更 8 个。由于 L3/L4 饱和，**class 由 L2（distal QKI peak）驱动**。Top-B（class→时点数→|ΔPDUI|→PAP）：Sirt2(4)、Aplp1(3,PAP+)、Cnp(2)、Pea15a(2)、**Atp2a2(2,PAP+)**、Kazn(1,PAP+)、Plec(1,PAP+)、Fam107a(1)；**Atp2a2 D→B**；**Ndrg2/Agpat3 B→C**（程序化确认 Proposal §2.2 预测）。

## 四、P3 六样本 BAM 核查（sham×2 / d3×2 / d7×2）

**P3.1**：6/6 quickcheck 通过，idxstats 全染色体在位。**富 A 内部引物风险：10/10 low**（最大 A-run ≤5、A% ≤34%）——3'RACE 无风险位点。

**P3.3 覆盖比**（distal_mean/shared_mean，MQ20+BQ20；判读依据：缩短事件卒中后 distal 占比应下降）：

| 基因 | sham | d3 | d7 | d7 降幅 | 一致性 | 置信 |
|---|---|---|---|---|---|---|
| **Atp2a2** | 0.302 | 0.084 | 0.085 | **-71.7%** | d3+d7 全一致 | **high** |
| **Agpat3** | 0.471 | 0.201 | 0.083 | **-82.3%** | 全一致 | **high** |
| Aplp1 | 0.299 | 0.291 | 0.159 | -46.9% | d7 一致 | medium（⚠总 RNA 同步暴跌 4823→242，需控制） |
| Sirt2 | 0.381 | 0.453 | 0.240 | -36.9% | d7 一致 | medium |
| Ndrg2 | 0.542 | (杂) | 0.365 | -32.7% | d7 一致 | medium |
| Kazn / Plec / Fam107a | — | — | — | 28–33% | d7 一致 | **low（深度薄，d7 总深 71–204）** |
| Cnp / Pea15a | — | — | — | 19.3% / 22.4% | 一致 | hold（<25% 阈） |

Atp2a2 附加观察：shared 侧组均值深度 sham→d7 反而上升（291→565），**排除"总 RNA 下降"假阳性解释**。已知 d3_rep1 文库异常（insert=150）与 d3 部分杂音相容。

**P3.2/P3.4**：10 基因固定尺度覆盖图（`results/03_bam_review/plots/*.png`，公共 y 轴 + PAS/distal/QKI 标注）+ IGV 快照指令单（GUI 截图可随时补做，不阻断）。事件审查表：**8 advance（high 2）/ 2 hold / 0 drop**。完成标准达成：至少一候选过读段+探针核查 ✓；Atp2a2 的 advance 决定证据链完整 ✓。

## 五、P4 公共数据辅助

### GSE143531（健康海马 PAP-TRAP；3 动物 × Full/PAP × 8 技术重复，合并为 6 组）
- **可测性（Gate 1 加分项判定要求）**：主靶与全部备靶在 PAP 隔室**可检出** ✓
- **富集方向（如实）**：相对全基因组 log2(PAP/Full) 分布（6,762 可用基因，中位 −3.42），候选**均无 PAP 富集**：Atp2a2 **54.0% 百分位（中位）**、Agpat3 51.8%、Aplp1 32.8%、Sirt2 15.0%、Ndrg2 13.9%、Cnp 4.7%、Pea15a 22.0%、Fam107a 6.8%、Plec 95.8%（唯一富集者，但 P3 侧低置信）；Kazn 2/3 动物不可检出
- **结论**：v1 的 PAP 富集先验（源自 GSE74456 皮层 PAPTRAP）**不迁移**至健康海马数据集，降级为"**区域依赖性待验证**"；项目假设（卒中诱导的异构体分布改变）不依赖基线富集，将由脑片空间实验在皮层卒中样本中直接测量

### GSE330741（SN-MPRA）
本地 651 基因/1,631 tile 计数表中 **Glt1/Slc1a2、Sparc、Hsbp1 及全部候选均不存在**（并见 Excel 日期损坏基因名如 "9-Sep"）——与论文描述的 focused tiling 文库矛盾，**manifest 核对仍未解决**（沿袭 v5.1 P4-04 BLOCKED）。reporter 阳性对照需从 Koester 2026（PMID 42094343）正式 supplementary 获取或改用文献元件。**不阻断 Gate 1**（MPRA 不在证据链内）。

## 六、Gate 1 定靶（用户已批准）

| 角色 | 基因 | event | class v2 | 核心证据 |
|---|---|---|---|---|
| **主靶** | **Atp2a2** | 6953 | B | QKI distal peak FDR 1.04e-5 + PAS 距已知 **13bp** + 覆盖比 **-71.7%**（d3/d7 双重复全一致）+ PAP+ + 富 A low + 探针可行 |
| 备靶 1 | Agpat3 | 216 | C | 覆盖比 -82.3% 全一致（4 时点显著；QKI 证据消失如实） |
| 备靶 2 | Aplp1 | 7821 | B | 覆盖比 -46.9% + 自有 distal QKI peak FDR 0.02；⚠总 RNA 塌陷需在湿实验设计中控 |
| 替补 | Sirt2 | — | B | 4 时点；PAP 百分位 15% 如实 |
| 观察 | Ndrg2 | 2072 | C | PAS 8bp 极佳；QKI 消失、d3 不一致 |

必须满足项六项全 ✓；加分项 4/5（第二方法已补算：**Atp2a2 salmon 长异构体占比（重复1/2）sham 0.719/0.541 → d1 0.059/0.105 → d3 0.130/0.112 → d7 0.039/0.039 → d21 0.066/0.059 → d60 0.275/0.127，与 PDUI 轨迹逐时点同向、双重复一致**；Agpat3 同向；Aplp1 注释仅 1 条转录本，转录本级不可算，如实登记）。

**湿实验设计区间**（详见 `results/05_target_nomination/target_nomination_v1.md`，三靶全部给出）：3'RACE 外/内引物区、ddPCR common/distal 扩增子区、smFISH 铺瓦区（BED 0-based，负链已核对）。示例主靶：distal 段 chr5:122,453,512–122,456,498（2,986bp）；shared 段 chr5:122,456,498–122,457,153（655bp）；PAS 122,456,498。

## 七、新旧结论对照（v1 错侧 → v2 校正）

| 项目 | v1（错侧） | v2（校正后） |
|---|---|---|
| Atp2a2 QKI | peak 122,456,550–650"在丢失段" | 该区域在 shared；**真 peak 122,454,200–250 在 distal**（FDR 1.04e-5） |
| Atp2a2 class | D | **B** |
| Ndrg2/Agpat3 QKI | "丢失段重叠" | **消失**（旧 peak 落 shared），class B→C |
| 有 distal QKI 的事件 | 11/977 | **42/9,693（977 内 9 个：8 个 class B + Qk 自身事件落 none）** |
| Motif 富集 | 阴性 | **阴性复现**（16 家族 q≥0.71） |
| 保守率 | 65%（1,004/1,541） | **99.7%**——v1 系 bedGraph 中间管线伪影；L4 阈值无区分度 |
| L3 覆盖 | 730/977 | 792/824 候选（96%）——饱和，无区分度 |
| 候选排序 | Ndrg2≈Agpat3>Cnp≈Pea15a>Fam107a | **Atp2a2 主靶**；Agpat3/Aplp1 备靶；Ndrg2 观察 |

**不变的结论**：逐样本 PDUI/ΔPDUI/时间轨迹不受区间方向影响——Atp2a2 缩短（sham 0.19/0.20 → d3 0.08/0.06 → d7 0.09/0.09）与全局缩短主导叙事维持。

## 八、诚实登记（负面结果、限制与勘误全清单）

1. L3（motif）与 L4（保守性 ≥10bp 二值）在本次数据上**均无区分度**（96% / 99.7% 基线）；class 区分度由 L2 承担
2. **GSE143531 不支持候选的 PAP 富集先验**（Atp2a2 中位）；先验降级为区域依赖待验证
3. **GSE330741 manifest 矛盾未解**；reporter 阳性对照待 Koester 2026 supplementary
4. **Aplp1 第二方法不可算**（单转录本注释）；其总 RNA 塌陷需湿实验控制
5. **Kazn 健康海马 PAP 2/3 动物不可检出**；Plec/Fam107a/Kazn 深度薄判 low 置信
6. **PolyASite 旧域名 polyasite.ethz.ch 已废弃**（权威 NXDOMAIN）；现址 polyasite.unibas.ch；本轮数据取自 GRCm38.96 atlas（301,006 clusters）
7. **13/27 口径**：DaPars2 vs **Salmon** 长异构体比例（QAPA 未构建成功）；27 项不含 Atp2a2（本轮已单独补算并同向）
8. FIMO 5.5.9 `--thresh` 按 p 值截断；fimo.tsv 末尾有 footer 行
9. DaPars2 PAS 永远落在 20%/150nt/5% 修剪窗内（结构性限制）；一切 PAS 均为"预测断点"
10. 本轮 6 样本核查窗口为 sham/d3/d7——d1/d21/d60 事件（如 Fam107a 的 d1 显著）未在 BAM 核查窗口内，需要时可从全量 BAM 同法补切
11. v1 遗留资产中 `lost_segments.bed` 及基于它的 QKI/motif/保守性/class 字段仅作审计对照，**禁止用于选靶**（已复制为 `lost_segments.OLD_WRONGSTRAND.bed` 存证）

## 九、交付物清单

- 定靶：`results/05_target_nomination/target_nomination_v1.md`（含三靶全部设计区间）
- P1：`results/01_strand_audit/`（区间表、distal/shared BED+sorted、rebuild_notes.md）+ `config/coordinate_convention.yaml` + `scripts/01_build_strand_aware_segments.py` + `scripts/test_01_segments.py`
- P2：`results/02_candidate_rebuild/`（QKI 交集×2+逐事件汇总、fimo 双扫描+家族富集+null QC、保守性表、`event_annotation_matrix.v2.tsv`、新旧差异表×3、`candidate_cards_v2.md`、`candidate_ranking_v2.tsv`、`pas_distance_977.tsv`、`salmon_gate1_baseline.tsv`）
- P3：`results/03_bam_review/`（quickcheck/idxstats、regions、6×depth、覆盖比表、富 A 表、10×PNG、IGV 指令单、`event_review_table.tsv`）
- P4：`results/04_public_pap/`（GSE143531 汇总+百分位表）
- P0：`results/00_inventory/`（sha256 清单、zip 清单、软件版本、provenance）
- 审计：`scripts/audit_pass1_numbers.py`（34 检全过）、`scripts/audit_pass23.py`（48 检全过）、红线 grep 记录
- 数据：`data/`（bam 12 样本、references 含 genome.fa/GTF/PAS/phastCons bw/salmon_idx、salmon、dapars2_full、qc、GSE143531）；`.gitignore` 排除大文件

## 十、下一步（湿实验 WP1）

1. 主靶 Atp2a2 3'RACE 引物对（外侧 + 嵌套）按提名文书区间交合作方 Primer3 设计；附 no-RT/no-template 对照与已知 3'端阳性对照
2. ddPCR common/distal 扩增子设计与标准品线性验证
3. 动物实验对齐 GSE238125 模型（永久性远端 MCA 结扎 + 7min 双侧 CCA 阻断），d7 主时点，每组先导 4 只，随机化+盲法按 Proposal §5.1.4/§11
4. Gate 2A 判定：≥2 个有 poly(A) 接合证据的 RNA 3'端 + 独立动物比例同向
5. 停走规则全程生效（Proposal §6 结果解释框架 / TODO Gate 3 表）
