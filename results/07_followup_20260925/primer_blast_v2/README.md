# primer_blast_v2/README.md — NCBI Primer-BLAST 补跑记录（2026-09-25）

数据库：RefSeq Reference / mRNA（GRCm39 现行版本；**辅助性质**——主判定为
mm10/GRCm38 + GENCODE vM25 的本地逐碱基与 ≤1 错配扫描，见 f4 与 assay_design_verified.tsv）。

| job | 模式 | 结果 | 判读 |
|---|---|---|---|
| long_junction_v2_pair.html | pair (F=TGTGGTGTTTTCCTCCAATGCCT / R=CACGCACCCGAACACCCTTATAT) | intended targets 上 **product length = 1692** | 与本地基因组产物一致 ✓ |
| qPCR_S_spliceform_pair.html | pair (J_S / distal_F) | **No target templates** | 预期：J_S 跨 S 独有剪接接头，参考序列中不存在连续拷贝 → 无任何基因组扩增 ✓ |
| qPCR_L_orthogonal_pair.html | pair (J_L / distal_F) | **No target templates** | 同上（L 接头跨接） ✓ |
| outerGSP_v2_single | 单引物 specificity | **NCBI 不再接受脚本化单引物提交**（3 次尝试均退回表单，存证 `*_submit_fail*.html`） | 见下 |
| nestedGSP_v2_single | 同上 | 同上 | 见下 |

## 两个 RACE GSP 的等效核查（05 任务书允许"Primer-BLAST 或等效近似匹配检查"）

1. 本地（主）：outer_v2/nested_v2 在 5 条 Atp2a2 转录本 + Gm30970 + chr5 基因座
   （122.40-122.51Mb）双链精确与 ≤1 错配扫描，均唯一命中目标且 Gm30970 最小错配
   13/12（f4_assay_verification.py 日志）。
2. 残余缺口：chr5 基因座以外的小鼠全基因组近似匹配未做（如需可由实验负责人在下单前
   用浏览器手动跑一次 Primer-BLAST 单引物 specificity——网页交互可用，脚本通道已关）。
3. 06 包声称的"GSP 单引物 Primer-BLAST 完成"实未成功（06/primer_blast/ 只有 3 个
   相同的失败提交页），本轮已如实修正记录。
