# 一页决策摘要（任务书 §8.6，2026-09-23）

## 决定：**ADVANCE Atp2a2 → 3'RACE**

| 判据（任务书 §9） | 状态 | 证据 |
|---|---|---|
| 表格 join 修复 | ✅ | `event_review_table.v2.tsv`：Atp2a2 PDUI 0.19/0.20、0.08/0.06、0.09/0.09 全部实测（v1 全 NA 的列名 bug 已修）；second_method 接入 salmon 并标注 same-cohort cross-check；probe 字段降为 "region identified; oligo design pending" |
| 局部 BAM/IGV 读段无足以解释结果的明确问题 | ✅ | 六样本：MAPQ 中位 255、低质量比对 <0.7%、soft-clip ≤0.40%、UTR 区剪接读段 0%、错配 NM 中位 0；六样本覆盖均延伸至注释 UTR 末端（122457152）；distal 下降形态与"约 9% 长版本残留"自洽；无重叠基因/异常剪接可解释信号（指标表 `Atp2a2_read_metrics.tsv`、逐样本图 `igv/Atp2a2_*.png`） |
| 序列级引物可特异设计 | ✅ | 外侧/嵌套 GSP、common/distal ddPCR 引物全部产出（Primer3，Tm 59.7–62.0）；**distal 引物零脱靶**（已移位避开 ENSMUST00000177974.7 重叠 UTR）；特异性筛查为转录组精确匹配代理，Primer-BLAST 终检待实验负责人 |

## 同轮修正的辅助结论（不影响上述决定，但影响论文表述）

- **GSE143531 归一化重算**（edgeR TMM，6 文库）：PAP 文库深度约为 Full 的 1/10，
  v1 未归一化百分位被深度差主导；归一化后 Atp2a2 log2FC=-0.273、FDR=0.683、
  百分位 55.4%——"健康海马中不偏向 PAP"的定性结论成立并强化。
  该数据集仅承担 gene-level exploratory 排序，不是 Atp2a2 进 3'RACE 的依据
  （依据 = DaPars2 + 覆盖比 + salmon 三者方向一致 + QKI/PAS 结构支持）。
- **PAP 基因列表来源**已审计：2,729 列表 = 课题组对 GSE74456 的自重分析
  （相对富集口径），Atp2a2 为成员（4.4 倍，padj 8e-7）；Sakers 官方补充表未对表。

## 备靶状态

- **Agpat3**（备靶 1）：未做同深度核查与引物设计；如 Atp2a2 在 3'RACE 失败，
  按相同管线（BAM 核查→assay design）切换，预计 1–2 个工作日。
- **Aplp1**（备靶 2）：distal 区仅 163bp 且总覆盖同步下降，按任务书 §7.4 不作等价替代。
- **Sirt2/Plec/Ndrg2**：维持替补/观察。

## 缺什么（进入下一步前的开放项）

1. Primer-BLAST 网页终检（引物序列已给，5 分钟人工操作）
2. ddPCR 扩增效率/线性范围验证（湿实验侧）
3. 正式样本量计算所需的先导方差（Gate 2A 通过后由首批 4+4 只估计）
4. Primer-BLAST 之外的转运相关 RBP 文献支持（QKI 之外）——非阻断

## 对外表述基线（任务书 §10）

> GSE238125 的星形胶质细胞 RNA-seq 分析将 Atp2a2 event 6953 提名为卒中相关长 3'UTR
> 使用下降的候选。链方向校正后的 BAM 覆盖和同一队列的 Salmon 定量方向一致。该位点及
> 其长、短 RNA 版本仍需在独立动物中通过 3'端测序和异构体检测确认，PAP 中的空间分布
> 尚待原位实验测量。
