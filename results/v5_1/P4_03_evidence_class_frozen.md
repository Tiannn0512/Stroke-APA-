# P4-03 Evidence Class 判据（冻结版 v1.0）

> 冻结时点：2026-09-20，**在计算任何交集/候选判定之前**（P4-01 尚未运行；P4-03 时点规则经 2026-09-19 审计修订并由用户采纳）。
> 冻结后本文件只读；任何修改须显式登记理由与时间。
> 输入工件全部先于本冻结存在：fimo_lost 扫描（09-20）、QKI peaks bed（09-17）、层2 事件表（09-19）。

## 1. 证据层定义（逐事件）

| 层 | 判据（冻结） | 来源 |
|---|---|---|
| L1 文献 | 事件基因 ∈ `results/p3_level1_pairs.tsv` 的 12 对靶基因 | 已有 |
| L2 CLIP | 事件丢失段与 GSE147119 QKI peak 重叠（`results/p3_level2_qki_intersect.tsv`，已算） | 已有 |
| L3 motif | fimo p < **1e-4**（扫描阈 1e-3 的 10×收紧）命中于丢失段，且 **family ∈ 冻结清单**（见 §2） | 已有扫描，事后按冻结参数过滤 |
| L4 保守性 | 丢失段与 mm10 phastCons 多物种**保守元件**重叠 ≥10 bp（二值） | 下载后计算，规则先冻结 |
| L-MPRA | **本周期不可用**（P4-04 BLOCKED）——A 级判定本轮不可达，如实登记 | — |

**冻结清单（L3 family）= ATtRACT 16 family ∩ P1 目录 31 regulator − FUS（L1 阴性）**：
QKI, ELAVL1, PTBP1, MBNL2, NOVA1, NOVA2, TARDBP, HNRNPA2B1（8 个；MBNL1/PTBP2/KHDRBS1/HNRNPA1/ZFP36/PABPN1 不在 concordant 集内，排除）。
> 精确化说明：以 P1 目录的 concordant 集为准做交集，脚本按 gene_name 大小写不敏感匹配（Qk→QKI）。

## 2. Class 判定（可枚举，无权重）

- **A 级**：L-MPRA 阳性 —— 本轮不可达（登记：MPRA 阻断所致，非候选之过）。
- **B 级**：L2 ✓ + L3 ✓ + L4 ✓
- **C 级**：L3 ✓ + L4 ✓（无 L2）
- **D 级**：L3 ✓（无 L2/L4）
- 无 L3 命中 → 无 class（登记为 none）。

## 3. 候选排序规则（冻结）

1. class 升序（B > C > D > none）；
2. 联合 BH padj 升序；
3. |ΔPDUI| 降序；
4. 显著时间点数降序；
→ 取前 **0–5** 为 top 候选（0 = 如实报零）。

## 4. 表述约束（照抄红线）

L1 事件可写"调控"；L2 只写"结合证据一致/预测丢失一段 QKI 结合区"；L3 只写"预测相关"；全文结论上限 = "altered localization potential"。阴性层照登：L3 家族层面富集阴性、MPRA 不可用、top40 单事件抽验功效不足。
