# envs —— WSL conda 环境导出（2026-09-23）

由 `conda env export -n <env>` 生成，还原：`conda env create -f <name>.yml`。

| yml | 用途 | 关键包 |
|---|---|---|
| `dapars2.yml` | APA 定量 | python 2/3 + DaPars2 依赖（源码本体在 `tools/DaPars2_patched/`） |
| `rlimma.yml` | P2 统计 | R 4.3.3 + limma 3.58.1（edgeR 4.0.16 需另装，见 REGENERATE.md §1） |
| `bioinfo.yml` | 杂项分析 | pysam、pyBigWig、primer3-py 等 |
| `salmon.yml` | 转录本定量 | salmon 1.10.3（注意：salmon 在 envs 顶层，不在 bioinfo 内） |
| `meme.yml` | motif 扫描 | MEME suite / FIMO 5.5.9 |
| `qapa.yml` | QAPA 第二工具（构建未成功，仅存档） | qapa |

注意：
- conda yml 依赖解析在 Linux/WSL 侧有效；Windows 原生不可用。
- scRNA 处理（scanpy/pydeseq2）当年跑在 base 或临时环境，未单独冻结；若需重建按 `results/v5_1` 脚本头部的 import 清单手工安装（scanpy、pydeseq2、scrublet 降级路径=MAD 过滤，见 `docs/v5_1/P1_frozen_params.md`）。
- 软件版本权威记录：`results/reanalysis/00_inventory/software_versions.tsv` 与 `results/v5_1` 冻结文档。
