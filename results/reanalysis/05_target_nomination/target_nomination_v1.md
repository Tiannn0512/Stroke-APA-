# Gate 1 靶标提名 target_nomination_v1

生成：2026-09-21 ｜ 决策人：用户（"按你推荐的拍板，批准执行"）
依据：P1 链方向重建 + P2 注释重建（QKI/motif/保守性/PAS）+ P3 六样本 BAM 核查 + P4 公共数据辅助

## 决定

| 角色 | 基因 | event | class v2 | 关键证据 | 置信 |
|---|---|---|---|---|---|
| 主靶 | Atp2a2 | 6953 | B | day3;day7；覆盖比 d7 降幅 71.7%（双重复一致）；PAS 距已知 13bp；QKI distal chr5:122454200-122454250(fdr=1.04e-05) | high |
| 备靶 1 | Agpat3 | 216 | C | day1;day3;day7;day21；覆盖比 d7 降幅 82.3%（双重复一致）；PAS 距已知 171bp；QKI distal  | high |
| 备靶 2 | Aplp1 | 7821 | B | day1;day7;day21；覆盖比 d7 降幅 46.9%（双重复一致）；PAS 距已知 182bp；QKI distal chr7:30435000-30435050(fdr=2.02e-02) | medium |

替换候选（第 4 序位）：Sirt2（class B、4 时点、PAP 百分位 15%——如实记录其健康海马 PAP 低于中位）。
Ndrg2：保留观察（PAS 8bp 极佳，但校正后 QKI 证据消失、d3 覆盖比不一致、PAP 百分位 13.9%）。
Cnp/Pea15a：hold（覆盖比降幅 <25%）；Kazn/Plec/Fam107a：low 置信（深度不足）；Kazn 另在健康海马 PAP 中 2/3 动物不可检出。

## 各靶标设计区间（mm10，BED 0-based half-open，负链基因转录方向=基因组低→高坐标的反向）

3'RACE 基因特异正向引物区 = 预测近端 PAS 的转录上游（负链=更高坐标）；common 扩增子/探针位于 shared 段；distal 扩增子/探针位于 distal 段且避开 shared 边界 ±100bp。

### Atp2a2（主靶）event 6953 — chr5 -
- 3'UTR 全长：122453512–122457153（DaPars2 BED）
- distal 段（仅长异构体）：**chr5:122453512-122456498**（2986bp）
- shared 段（全部异构体）：**chr5:122456498-122457153**（655bp）
- 预测近端 PAS：122456498（PolyASite 最近 13bp）
- 3'RACE 外侧引物区：chr5:122456598-122457098；嵌套内侧引物区：紧邻 PAS 转录上游 200bp 内
- ddPCR/qPCR common 扩增子区：chr5:122456558-122457103
- ddPCR/qPCR distal 扩增子区：chr5:122453562-122454712
- smFISH 探针铺瓦区：common chr5:122456558-122457153；distal chr5:122453552-122454512

### Agpat3（备靶 1）event 216 — chr10 -
- 3'UTR 全长：78269177–78273689（DaPars2 BED）
- distal 段（仅长异构体）：**chr10:78269177-78271713**（2536bp）
- shared 段（全部异构体）：**chr10:78271713-78273689**（1976bp）
- 预测近端 PAS：78271713（PolyASite 最近 171bp）
- 3'RACE 外侧引物区：chr10:78271813-78272313；嵌套内侧引物区：紧邻 PAS 转录上游 200bp 内
- ddPCR/qPCR common 扩增子区：chr10:78271773-78273639
- ddPCR/qPCR distal 扩增子区：chr10:78269227-78270377
- smFISH 探针铺瓦区：common chr10:78271773-78273689；distal chr10:78269217-78270177

### Aplp1（备靶 2）event 7821 — chr7 -
- 3'UTR 全长：30434981–30435295（DaPars2 BED）
- distal 段（仅长异构体）：**chr7:30434981-30435144**（163bp）
- shared 段（全部异构体）：**chr7:30435144-30435295**（151bp）
- 预测近端 PAS：30435144（PolyASite 最近 182bp）
- 3'RACE 外侧引物区：chr7:30435244-30435245；嵌套内侧引物区：紧邻 PAS 转录上游 200bp 内
- ddPCR/qPCR common 扩增子区：chr7:30435204-30435245
- ddPCR/qPCR distal 扩增子区：chr7:30435031-30435094
- smFISH 探针铺瓦区：common chr7:30435204-30435295；distal chr7:30435021-30435104

## Gate 1 必须满足项核对（主靶 Atp2a2）

| 条件 | 状态 |
|---|---|
| 链方向重建+单元测试通过（P1） | ✓ 17/17 |
| BAM 无明显比对/注释伪影（P3，10 基因富 A 全 low、quickcheck 6/6） | ✓ |
| 两个发现样本方向一致（d3 0.084/0.083；d7 0.085/0.086 覆盖比；PDUI 0.08/0.06、0.09/0.09） | ✓ |
| distal 区可设计特异探针/引物（2,986bp distal，655bp shared，均 >150bp） | ✓ |
| 预测 PAS 附近无严重内部引物风险（max A-run 3、A% 24） | ✓ |
| d7 主时间点同向（覆盖比与 PDUI 同向缩短） | ✓ |

加分项：多时点同向 ✓（d3+d7）；公共 PAS ✓（13bp）；GSE143531 基因级可测 ✓；校正后 distal RBP ✓（QKI FDR 1.04e-5）；第二方法 ✓（补算完成 2026-09-21：salmon TPM 长异构体占比 sham 0.719/0.541 → d3 0.130/0.112 → d7 0.039/0.039 → d60 0.275/0.127，与 DaPars2 PDUI 轨迹逐时点同向、双重复全一致，见 salmon_gate1_baseline.tsv；Aplp1 注释仅 1 条转录本，转录本级比例不可算，如实登记 absent；Agpat3 同法已算）。

## P4 公共数据辅助的如实记录

1. **GSE143531（健康海马 PAP-TRAP，3 动物×2 分区×8 技术重复，已合并）**：全部主/备候选**可检出** ✓（Gate 1 加分项判定要求）。但相对全基因组分布，**候选未见 PAP 富集**：Atp2a2 位于 54% 百分位（中位）、Agpat3 52%、Aplp1 33%、Sirt2 15%、Plec 96%（唯一富集者，但 BAM 侧低置信）。结论：**v1 的 PAP 富集先验不跨区域/协议迁移，降级为"区域依赖性待验证"**；项目假设（卒中诱导的异构体分布改变）不依赖基线富集，将由脑片空间实验在皮层直接测量。
2. **GSE330741（SN-MPRA）**：本地 651 基因/1,631 tile 计数表中 Glt1/Slc1a2、Sparc、Hsbp1 及全部候选**均不存在**（并有 Excel 日期损坏基因名如 "9-Sep"）——与论文描述的 focused tiling 文库矛盾，**manifest 核对仍未解决**（沿袭 v5.1 P4-04 BLOCKED）。阳性对照需从 Koester 2026（PMID 42094343）正式 supplementary 获取或改用文献元件；**不阻断 Gate 1**（MPRA 不在证据链内，仅背景）。

## 下一步（进入湿实验 WP1）

- 主靶 Atp2a2 3'RACE 引物对（外侧+嵌套）按上表区间交湿实验合作方用 Primer3 设计（Tm 58–62、产物 80–250bp、外显子跨剪接不适用——3'UTR 区连续，靠位置特异性）。
- ddPCR common/distal 扩增子同法设计并用标准品验证线性范围后，再订 smFISH 探针组。
- 动物模型尽量对齐 GSE238125（永久性远端 MCA 结扎+7 分钟双侧 CCA 阻断），d7 主时点，每组先导 4 只。
