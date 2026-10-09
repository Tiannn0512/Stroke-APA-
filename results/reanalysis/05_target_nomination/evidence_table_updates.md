# 证据表更新（对应外部 02_evidence_table.md 的差异行，2026-09-23）

> 用法：把下列行并入课题证据表对应位置；原表保留，不覆盖。所有行按本轮预设口径标注证据等级。

| # | 证据项 | 本轮更新 | 证据等级 | 来源文件 |
|---|---|---|---|---|
| G1 | GSE143531 PAP 可测性 | **作废重算**：v1 未归一化百分位被文库深度差（PAP≈Full 的 1/10）主导；edgeR TMM 归一化后 Atp2a2 log2FC=-0.273、FDR=0.683、百分位 55.4% | gene-level exploratory（n=3 文库/组） | results/04_public_pap/edgeR_TMM_PAP_vs_Full.tsv、candidate_gene_summary.v2.tsv、gse143531_reanalysis.md |
| G2 | PAP 富集先验 | 维持"候选在健康海马无 PAP 富集"结论（归一化后成立）；Plec 为两脑区唯一一致富集基因（log2FC +3.86/+4.41） | 同上 | 同上 |
| G3 | PAP 基因列表来源 | 2,729 列表 = 课题组自重分析 GSE74456（相对富集口径）；Atp2a2 成员（4.4 倍，padj 8e-7）；Sakers 官方补充表未对表 | 列表来源审计 | results/04_public_pap/pap_gene_set_provenance.tsv |
| G4 | 事件审查表 | v2 表修复 PDUI 全 NA（列名 bug）与 second_method 未联接；Atp2a2 PDUI/覆盖比/salmon 三源齐备 | 数据修复 | results/03_bam_review/event_review_table.v2.tsv、event_review_table_QC.md |
| G5 | Atp2a2 读段核查 | 六样本读段级指标全净（MAPQ 中位 255、soft-clip≤0.4%、剪接 0%、NM 0）；覆盖均达注释 UTR 端；无伪影解释信号 | 读段级描述证据 | results/03_bam_review/Atp2a2_read_review.tsv/.md、igv/Atp2a2_*.png、Atp2a2_read_metrics.tsv |
| G6 | Atp2a2 引物/检测设计 | 3'RACE 外/内 GSP、common/distal ddPCR 引物序列级设计完成；distal 引物零脱靶（已避开 177974 重叠 UTR）；177974 转录本记为设计风险 | in-silico 设计（待 Primer-BLAST 终检） | results/05_target_nomination/Atp2a2_assay_design.tsv/.md、Atp2a2_specificity_screen.txt |
| G7 | 决定 | **Atp2a2 ADVANCE → 3'RACE**（三判据全过）；Agpat3 备靶待命；空间/机制分支按 Gate 2A/2B/3 顺序解锁 | 决策 | results/05_target_nomination/decision_summary_v1.md |

对外表述基线：见 decision_summary_v1.md 末节（预设 §10 建议句式）。
