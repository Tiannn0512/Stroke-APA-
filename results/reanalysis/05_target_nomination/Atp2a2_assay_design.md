# Atp2a2 实验设计（任务 D 交付，2026-09-23）

> 生成：`scripts/p4d_atp2a2_assay_design.py`（Primer3 2.3.1 via primer3-py；特异性筛查
> = 引物全长及反向互补序列在 GENCODE vM25 转录组（142,604 条）中的精确匹配）
> 坐标：mm10 BED 0-based half-open；Atp2a2 为负链，所有序列按转录方向给出（基因组切片的反向互补）。

## 1. 结构性发现（设计前必须知道的转录本背景）

GENCODE vM25 注释 Atp2a2 有 **3 条带 3'UTR 的转录本**：

| 转录本 | 3'UTR (BED) | 与事件 6953 的关系 |
|---|---|---|
| ENSMUST00000179939.7 | 122453512–122457153 | **事件载体**（远端 PAS≈122453512，预测近端 PAS 122456498） |
| ENSMUST00000177974.7 | 122453513–122454293 | **设计风险源**：其 3'UTR 压在事件 distal 段前部（QKI peak 122454200–250 亦在其内） |
| ENSMUST00000031423.9 | 122456339–122457153 | 仅 shared 区后段，与 common 检测共存属预期 |

设计含义：distal 检测必须放在 **122454293 以远**（179939 独有区段），否则会交叉检测
177974 的 RNA。首轮设计的 distal 引物（122454237–122454394）确实命中了 177974——
已按此发现移位重设计。

## 2. 引物设计结果（Primer3 输出，完整字段见 Atp2a2_assay_design.tsv）

| 实验 | 引物 | 序列（转录方向） | Tm | 产物 | 脱靶精确命中 |
|---|---|---|---|---|---|
| 3'RACE 外侧 GSP | forward | TTGGCTGGTGAAGGAGGTTTCATA | 62.0 | — | 1（Atp2a2 其他转录本，见 §3） |
| 3'RACE 外侧 GSP（备选） | forward | TATTGGCTGGTGAAGGAGGTTTCA | 62.0 | — | 1 |
| 3'RACE 嵌套 GSP | forward | TGTCTCGCTTGTGCAGAAAACATT | 62.0 | — | 1 |
| 3'RACE 嵌套 GSP（备选） | forward | TCGGCTTCTTCTGTAGGTCTGATG | 61.9 | — | 1 |
| qPCR common 对 | F / R | AGTCCTGCTACCTCAGTCAGTA / CAGACCTACAGAAGAAGCCGAG | 59.7/60.2 | 区间内 | 1（=31423，common 检测设计预期） |
| qPCR distal 对 | F / R | TGACTATGGCACTAACCACCAC / GAAGGCATTGGAGGAAAACACC | 60.0/60.0 | 区间内 | **0** |

外侧 GSP 区 chr5:122456598–122457098；嵌套区 chr5:122456518–122456698；
common 扩增子区 chr5:122456558–122457103；**distal 扩增子区（修订后）chr5:122454500–122456300**
（179939 独有区段，原提名区间 122453562–122454712 的前半与 177974 的 3'UTR 重叠，已移位）。

## 3. 风险与待办

1. **177974 转录本**：若其有表达，"distal 检测"（修订后）不会检出它；它的 3' 端与
   179939 的远端 PAS 同区（~122453512），解释上属"使用远端位点的另一转录本"。对
   event-6953 特异的 APA 测量，排除它是正确选择；同时建议湿实验第一轮用
   3'RACE 直接观察该区域的全部 poly(A) 接合，一次性厘清转录本结构。
2. **smFISH distal 探针铺瓦区**：同理应改用 179939 独有区段 chr5:122454293–122456400
   （原提名 122453552–122454512 与 177974 3'UTR 重叠，已不建议）。
3. **Primer-BLAST 未做**（本地无 BLAST 服务）：特异性筛查为转录组精确匹配代理；
   引物订单前应由实验负责人在 NCBI Primer-BLAST 对 mm10 + RefSeq/GenBank 做最终确认。
4. ddPCR 上机前须完成：标准品/梯度验证扩增效率与线性范围；no-RT/no-template 对照；
   熔解曲线或凝胶确认单产物。
5. 本设计只覆盖 common/distal 两类检测；短版本直接检测（junction 特异）需等 3'RACE
   给出真实接合序列后再补设计（任务书 §7.3 口径）。
