# -*- coding: utf-8 -*-
"""Rebuild inventory rows 38-63 with correct UTF-8 text (heredoc mojibake repair),
keeping 579->471 correction applied."""
import io

p = r"D:\stroke_apa\results\data_inventory.tsv"
lines = io.open(p, encoding="utf-8").read().split("\n")
keep = lines[:37]  # rows 1-37 (header + 2026-09-16..09-18), verified clean

rows = [
 '2026-09-19\treviews/Stroke_APA_v5.1_论证逻辑审计报告_20260918.md\t外部审计（微信原件拷贝入库）\t849 行\tMajor revision/有条件 Go；v5.1r 修订依据\t审计',
 '2026-09-19\tresults/P2_SAP_frozen.md\t本次编写（审计 §3/§4/§6 细则）\tv1.0/v1.1（R 例外修订见文内）\tP2 统计分析计划冻结（全量统计前）\tP2-06a',
 '2026-09-19\tresults/p1_08_grid_205v409.tsv\tscripts/p1_08_205v409.py\t5 变体 × 5 列\t205 vs 409 判据分解网格\tP1-08',
 '2026-09-19\tresults/p1_286075_paper_criteria_repro.tsv\tscripts/p1_08_205v409.py (V1)\t10,862 基因，55 sig\t原文判据复现全表（anchor Hspa1a/Hspa1b rank1/2）\tP1-08',
 '2026-09-19\tresults/p1_08_286075_205v409_appendix.md\t本次编写\t4 节\t205 vs 409 正式复现附录\tP1-08',
 '2026-09-20\tresults/p3_m4_set1_paptrap_de.tsv + p3_m4_set1_pap_enriched.txt\tscripts/p3_m4_reference_sets.py (Set1)\t15,550 入模；2,729 PAP 富集\tPAP_localized 参考轴\tM4 参考集1',
 '2026-09-20\tresults/p3_m4_set2_endfoot_de.txt + _background.txt\tscripts/p3_m4_reference_sets.py (Set2)\t409 DE；14,320 背景\tendfoot_stroke_responsive 参考轴\tM4 参考集2',
 '2026-09-20\tresults/p3_m4_set3_stroke_responsive_{cortex,whitematter}.txt + zone_de_*.tsv\tscripts/p3_m4_reference_sets.py (Set3)\t皮层并集 4,773；白质并集 10,434\tstroke_spatial_zone_response 参考轴（10 zone contrast 全存）\tM4 参考集3',
 '2026-09-20\tresults/p1_07_sensitivity_summary.tsv + p1_07_variant_{DE,subpop,regulators}_*.tsv\tscripts/p1_07_doublet_ambient.py\tscrublet 1,624 双联体(2.85%)；ambient 2,849 (5.0%)；B/C 变体 astro 3,046/3,042 vs 原 3,045；DE_sig 1,310/1,311 vs 原 1,280（overlap 96%/94%）；31 regulator 方向翻转 = 0\tdoublet/ambient 敏感性：P1 结论稳健\tP1-07',
 '2026-09-20\tresults/p2_full_pdui_matrix.tsv\tp2_full_merge_qc.py（20/20 chr 合并；depth 22/22 日志核验）\t9,693 事件 × 12 样本 + 17 列\t全量 PDUI 矩阵（主表）\tP2-06',
 '2026-09-20\tresults/p2_full_merge_qc.tsv\t同上\t35 行指标\tsham r=0.9257；d1 4.84:1 缩短；d3 0.73:1（含技术异常，见层1）；d7 7.23:1；d21 6.58:1；d60 4.58:1\tP2-06',
 '2026-09-20\tresults/p2_full_qc_sample_level.tsv + _correlation.tsv + _pca.tsv + _stratified.tsv\tp2_sap_layer1_qc.py\t12 样本 QC + 混杂 + 分层\tday3 双样本 insert=150 异常（d3_rep1 dup 0.215）；d1 缩短随表达量升高而增强（Q4 5.87 vs Q1 2.17）；PC1 36.8%\tSAP 层1',
 '2026-09-20\tresults/p2_sap_layer2_events.tsv + _omnibus.tsv\tp2_sap_layer2_limma.R（limma robust+trend）\t9,693 事件检验\t显著（联合BH<0.1 & |ΔPDUI|≥0.1 & 双重复同向）：d1 292 / d3 159 / d7 435 / d21 496 / d60 392；并集 977（≥2 时间点 471，09-20 审计修正）；omnibus 4,910/5,923\tSAP 层2',
 '2026-09-20\tresults/p2_gate_g2_evidence_extras.tsv + p2_gate_g2_top_events.tsv\tp2_gate_g2_evidence.py\tLOSO/方向拆分/top40\tLOSO 稳定；d3 翻转系 d3_rep1 技术假象（剔除后 2.13:1 缩短）\tGate G2',
 '2026-09-20\tresults/p3_level2_qki_intersect.tsv + _summary.tsv\tp3_02_m4_enrich.py\t11/977 事件丢失段含 QKI peak（10 短缩）；OR=2.71 p=0.0125\tP3-02 L2：QKI 结合区丢失富集显著\tP3-02',
 '2026-09-20\tresults/p3_m4_enrichment.tsv\tp3_02_m4_enrich.py\t4 轴 × 2 基因集 Fisher\t重叠层面不富集（OR 0.96–1.16；仅 WM 边缘 p≈0.045）\tM4 Fisher',
 '2026-09-20\tresults/p3_m4_gsea.tsv\tp3_m4_gsea.R（fgsea 10k perm）\t24 检验 = 4 轴 × 6 排序\tPAP_localized 一致负向富集：d7 padj 1e-4 / d3 0.0018 / omnibus 0.0042 / d1 0.044；WM 同向 4 项过线；endfoot 轴 day3 反向 +1.28 后全阴\tM4 GSEA',
 '2026-09-20\tsalmon 转录组定量 12/12（WSL salmon_idx_vM25 + quant）\tp2_second_tool_chain2.sh\t142,604 转录本索引；12 quant.sf\t第二工具（G2-d 抽验）输入\tP2-G1(d)',
 '2026-09-20\tresults/p2_second_tool_consistency.tsv\tp2_second_tool_consistency.py（salmon TPM 长短异构体占比 vs DaPars2 dPDUI）\t5,499–5,847 基因/时点\t全 5 时点 Spearman 正相关显著（rho 0.041–0.139，p≤0.0023）——双管线方向一致；qapa 自动合并失败（空注释），走注册的自定义兜底\tP2-G1(d)',
 '2026-09-20\tresults/p2_second_tool_top40_signs.tsv\t同上\t27/40 可测；13/27 同号\t单事件一致性 = 功效不足（非矛盾）；单事件真实性归湿实验\tP2-G1(d)',
 '2026-09-20\tresults/p3_level3_family_enrichment.tsv\tp3_fimo_finalize_v4.py（vs UTR 匹配 null，Poisson 密度检验）\t15 family × lost 1.68Mb vs null 1.00Mb\tL3 全阴性：enrichment 0.89–1.07，无 family 过 BH（富集方向）；5 family 轻度缺失（TARDBP 0.92 等）\tP3-03 关闭（负结果照登）',
 '2026-09-20\t[作废] results/p3_level3_family_null_enrichment.tsv（随机基因组 null）\tp3_fimo_finalize.py v1\t设计无效：UTR vs 随机基因组 = 平凡 73×\t负结果/无效设计照登，勿引用\t作废记录',
 '2026-09-20\tresults/p3_level3_motif.tsv\tp3_fimo_finalize.py v3\t6,532 event×family 丢失段命中（描述性）\t逐事件 motif 命中表（fimo p 截断于 1e-3，不作逐事件检验——登记）\t描述层',
 '2026-09-20\tresults/p3_p4_event_annotation_matrix.tsv\tp3_04_matrix.py（L1/L2/L3/L4 逐事件标注，冻结规则）\t977 事件；class B=5 / C=475 / D=250 / none=247；方向 缩短811/延长153/混合13\t逐事件关联矩阵（P3-04 + P4 输入）\tP3-04',
 '2026-09-20\tresults/p4_top_b_class_candidates.tsv\tp3_04_matrix.py\t5 事件\tB 级候选（CLIP+motif+保守）：Ndrg2/Cnp/Fam107a/Pea15a/Agpat3，全短缩、多时点\tP4-08 预表',
 '2026-09-20\tphastCons60way elements（bigWig→bedGraph→score>0 合并 6,788 万元件）\tp3_l4_conservation.sh + p3_l4_merge/overlap.py\t丢失段 CONSERVED 1,004/1,541（65%）\tL4 保守层；随机基因组 null 版本已作废\tL4',
]
keep += rows
keep.append('2026-09-20\tnote: triple-factor (QKI-loss + PAP + shortened) = Atp2a2 (SERCA2) only; see final_qa_pap_triple_factor.tsv\t最终问答提取（Atp2a2 已补入终报告与路线书）\t—\t三要素案例\t补充')
keep.append("")

io.open(p, "w", encoding="utf-8", newline="\n").write("\n".join(keep))
print(f"inventory rebuilt: {len(keep)-1} rows (clean utf-8)")
