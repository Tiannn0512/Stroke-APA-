# P1-01 冻结参数（astrocyte 定义与 QC）— 冻结于 2026-09-16，先于任何 DE 计算

## QC 过滤（每样本内）
- min_genes = 200（细胞表达的基因数下限）
- min_counts = 500（细胞 UMI 总数下限）
- max_mt_pct = 20（线粒体基因 UMI 占比上限）
- doublet：Scrublet 每样本独立运行（sim_doublet_ratio=2, n_neighbors=30, min_counts=3, min_cells=3）；predicted_doublet=True 剔除。若 scrublet 在本环境不可用，降级方案 = 仅 MAD 离群过滤（counts/genes log2 MAD>3 剔除）并在此登记偏差。

## 标准化与聚类
- normalize_total(target_sum=1e4) → log1p → HVG（scanpy seurat flavor，n_top_genes=2000）→ PCA(n=50) → neighbors(n_neighbors=15, n_pcs=30) → leiden(resolution=1.0, flavor=igraph, n_iterations=2)
- UMAP(min_dist=0.3)

## astrocyte 簇判定（簇级，log1p CPM=1e4 尺度）
- 簇平均表达 Aqp4 ≥ 1 且 Slc1a3 ≥ 1 → astrocyte 簇
- 或：Aqp4 ≥ 0.5 且 Gfap ≥ 1 且 Slc1a3 ≥ 0.5 → astrocyte 簇（覆盖 Gfap^high 但 Slc1a3 偏低的反应性簇）
- 排除项：若某簇同时高表达微胶标志（C1qa/Cx3cr1 均值≥1）或内皮标志（Pecam1≥1）则不入选（防误纳）。

## astrocyte 亚群（细胞级）
- reactive_score = score_genes(Gene set: Serpina3n, C3, Gfap, Lgals3, S100a10? 否——冻结为最小集：Serpina3n, C3, Gfap；score 用 scanpy.tl.score_genes，ctrl_size=50, n_bins=25)
- reactive：reactive_score > 0.5；homeostatic：≤ 0.5
- pseudo-bulk 单元 = 样本 × 亚群（homeostatic/reactive/overall 三口径），每单元细胞数 < 20 则该单元不进 DE 并登记。

## DE
- GSE174574：pydeseq2，design = "~ condition"（overall astrocyte 口径，样本级 pseudo-bulk n=3+3）；亚群口径同设计，仅在两条件各 ≥3 个有效单元时跑。
- GSE286075（配对）：pydeseq2，design = "~ mouse + condition"（CTR/Stroke 各 3，同鼠配对）。
- 显著性：padj < 0.05（DE）；报告 log2FC。APA 事件的 FDR<0.1 属 P2，不在本冻结内。

## 冻结纪律
- 本文件在 P1-02/04/05 运行前 commit；此后修改须在 TODO【记录】登记理由。
