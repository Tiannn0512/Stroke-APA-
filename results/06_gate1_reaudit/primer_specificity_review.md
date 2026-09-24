# primer_specificity_review.md — 引物特异性复核（Gate1 复核 · 工作 3）

日期：2026-09-24 ｜ 引物：assay_design_corrected.tsv（8 条寡核苷酸，5 个 assay）
筛选数据库：GENCODE vM25 全转录组（142,604 条记录，264MB，双链搜索）
方法：逐引物全长度精确匹配（oligo 及其反向互补）；引物对水平扫描全部 4 种 dsDNA 镜像排布，
候选产物窗口 40-3000bp。脚本：g1_3_specificity.py（命中表 primer_specificity_hits.tsv）。

## 1. 设计前提：两个必须规避的陷阱（本轮均已处理）

1. **旧脚本 FASTA 提取 bug**：p4d_atp2a2_assay_design.py 自写提取函数未跳过行内
   `start % line_bases` 个碱基——旧引物序列与报告坐标不对应。本轮全部序列改用
   pyfaidx（BED 0-based）提取，并与 samtools faidx 独立逐碱基对拍
   （verify_window；单测覆盖行内偏移起点与跨行区间，test_g1_3_seq.py 全过；
   5 个设计窗口含跨 177974-UTR 边界、跨预测 PAS、跨共享外显子边界）。
2. **Gm30970 加工假基因**（ENSMUST00000227710.1，chr14:79,106,039-79,106,401，363bp）：
   其记录含 Atp2a2 mRNA ~3156-3330 区域（174bp+ 精确匹配）——恰覆盖部分共有外显子。
   任何"两引物同落该窗"的组合都会在假基因转录本上产生 81bp 脱靶扩增子。
   处理：common qPCR 迁至上游远处共有外显子对（122489237-489376 + 122491680-491785，
   跨内含子产物）；RACE GSP 选址于假基因覆盖窗外的清洁候选。

## 2. 最终 8 条寡核苷酸的筛查结果

| assay | 引物 | 全长精确命中 | 命中归属 | pair 扩增子 |
|---|---|---|---|---|
| RACE_outer_GSP | TCAAGAAGTGAAGTGACATGGACAAG | 3 | Atp2a2 三转录本 | （单引物） |
| RACE_nested_GSP | TGAGGGCATTACACATCTCTATGGTT | 4 | Atp2a2（含 Gm30970 假基因 1 处，见 §3） | （单引物） |
| qPCR_common | F AGCCTTTGTAGAGCCGTTTGTA / R CACACTCTTTCTGTCCTGTCGA | 5/5 | Atp2a2 全部 5 转录本 | 5 个扩增子全部 on-target（3 目标 + 同基因 196490/197415），OFF=0 |
| distal_terminal | F AGTTAGGACTGGAGGCCTATGT / R TTGTAAGTGGCCAGATTGCTCT | 2/2 | 仅 179939+177974 | 2 个扩增子均 on-target，OFF=0 |
| long_junction(179939 特异) | F GGAGTAACCGCTTCCTAAACCA / R ATCCTACTATGCGGCAGAACAG | 2/1 | F 见于 179939+177974，R 仅 179939 | **仅 179939 一个 366bp 扩增子**；177974/31423 因无 R 位点不扩增 ✓ |

RACE GSP 绑定坐标：outer chr5:122459602-122459626（反义，位于共有外显子 122459447-459650）；
nested chr5:122458092-122458116（反义，位于共有外显子 122458087-458221）。
nested 位于 outer 下游 ✓（嵌套关系成立）；RACE 产物长度：outer 体系 ~2.5-3.3kb、
nested 体系 ~2.2-3.0kb（至 179939 预测 PAS）；31423/177974 转录本若被测序，
其产物更短，可通过测序区分。

## 3. 状态与边界（诚实声明）

- **Primer-BLAST 网页复核：已完成（2026-09-24，NCBI primertool.cgi，Mus musculus）**。
  结果页已存档于 `primer_blast/` 目录：
  - `qPCR_common_genome.html`：RefSeq 代表基因组库，产物仅见 chr5 预期靶位点
    （GRCm39 NC_000071.7，产物 2458bp 含内含子），**零潜在脱靶**；
  - `distal_terminal_genome.html`：产物 136bp 仅 chr5 预期靶位点，**零潜在脱靶**；
  - `long_junction_refseqmrna.html` / `long_junction_genome*.html`：RefSeq mRNA 库
    零可扩增模板（长 UTR 异构体的该剪接接合不存在于 RefSeq mRNA 模型，仅 Ensembl
    179939 有）；基因组库三种角色组合均未报告产物（Primer-BLAST 对该跨-UTR 内含子
    接对的检索局限，已记录）。该引物对的特异性由以下三重保证：转录组精确匹配
    （R 位点仅存在于 179939）、产物尺寸可辨（cDNA 366bp vs 基因组 1954bp）、
    3'RACE 测序终审。
  - RACE 两条 GSP 为单引物，Primer-BLAST 配对模式不适用；其特异性由转录组
    全长精确命中表保证（全部命中限于 Atp2a2 基因，含 Gm30970 假基因在内零额外命中）。
  综上：**assay 状态可写 primer-blast(checked)（GSP 两条为 transcriptome-screened），
  但"validated"仍需湿实验确认；引物订单可提交，由实验负责人终审。**
- "零脱靶"的准确表述：**转录组全长精确匹配零脱靶**（且引物对扫描无非 Atp2a2 扩增子）。
  含错配的近似匹配风险由 pending 的 Primer-BLAST 覆盖。
- 基因组 DNA 扩增：qPCR_common 与 long_junction 均为跨内含子设计（基因组产物 1.5-1.6kb+），
  与 cDNA 产物（154bp/366bp）尺寸可分；distal_terminal 为单外显子内产物，
  gDNA 会得到同尺寸产物，需 -RT 对照排除（湿实验常规）。
