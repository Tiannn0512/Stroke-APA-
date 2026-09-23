# event_review_table.v2 QC

- join key: gene_symbol（ranking v2）<- 每子表同名列；PDUI 经 gene_symbol -> PDUI 矩阵 Gene 列（v1 bug：查询了不存在的 {sample}_PDUI 列，现改为矩阵实际列名）
- second_method 来源：results/02_candidate_rebuild/salmon_gate1_baseline.tsv，标注为 same-cohort computational cross-check（同一 GSE238125 队列的计算交叉核对，非独立验证）
- probe_status 统一为 "region identified; oligo design pending"（尚无 oligo 序列，不得写 feasible/validated）
- Atp2a2 验收：PDUI sham 0.19/0.20、day3 0.08/0.06、day7 0.09/0.09（非 NA）
