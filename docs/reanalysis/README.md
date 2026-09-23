# stroke_apa_reanalysis

缺血性卒中—星形胶质细胞 APA—PAP 定位：重建与湿实验验证工作区（2026-09-21 起）。

- 治理文件：`01_研究Proposal.md`（科学方案）+ `02_项目TODO.md`（执行清单 P0–P10）
- 关键背景：v5.1 干实验阶段的 lost-segment 定义存在**负链方向错误**，QKI/motif/保守性/候选分级全部待链方向重建（Proposal §2.2、TODO §1.2）；Atp2a2 人工校正锚点已写入 `config/coordinate_convention.yaml`
- 大文件不入 git（.gitignore）：data/bam（12 样本全量 BAM，自 WSL 迁入）、data/references（GRCm38/fa、GTF、salmon_idx）
- P0 已完成：输入 sha256 冻结、zip 清单、软件版本（results/00_inventory/）
