#!/usr/bin/env python3
"""Gate 1: target_nomination_v1.md generator.

User decision 2026-09-21 (approved): main = Atp2a2; backups = Agpat3 + Aplp1;
Sirt2 recorded as alternate. Pulls exact coordinates from ranking v2 and emits
per-target common/distal regions, 3'RACE primer zones, qPCR/smFISH design intervals.
"""
import csv, os

BASE = r"D:\stroke_apa_reanalysis"
R = os.path.join(BASE, "results", "02_candidate_rebuild")
OUT = os.path.join(BASE, "results", "05_target_nomination")
os.makedirs(OUT, exist_ok=True)

rank = {r["gene_symbol"]: r for r in csv.DictReader(open(os.path.join(R, "candidate_ranking_v2.tsv"), encoding="utf-8"), delimiter="\t")}
review = {r["gene_symbol"]: r for r in csv.DictReader(open(os.path.join(BASE, "results", "03_bam_review", "event_review_table.tsv"), encoding="utf-8"), delimiter="\t")}

def segs(sym):
    r = rank[sym]
    d = list(map(int, r["distal"].split("-")))
    s = list(map(int, r["shared"].split("-")))
    return r, d, s

targets = [("Atp2a2", "主靶"), ("Agpat3", "备靶 1"), ("Aplp1", "备靶 2")]
L = [
    "# Gate 1 靶标提名 target_nomination_v1",
    "",
    "生成：2026-09-21 ｜ 决策人：用户（\"按你推荐的拍板，批准执行\"）",
    "依据：P1 链方向重建 + P2 注释重建（QKI/motif/保守性/PAS）+ P3 六样本 BAM 核查 + P4 公共数据辅助",
    "",
    "## 决定",
    "",
    "| 角色 | 基因 | event | class v2 | 关键证据 | 置信 |",
    "|---|---|---|---|---|---|",
]
for sym, role in targets:
    r = rank[sym]
    rv = review[sym]
    L.append(f"| {role} | {sym} | {r['event_num']} | {r['class_v2']} | {r['contrasts_sig']}；覆盖比 d7 降幅 {rv['d7_drop_pct']}%（双重复一致）；PAS 距已知 {r['pas_dist_bp']}bp；QKI distal {r['L2v2_peaks'][:60]} | {rv['confidence']} |")
L += [
    "",
    "替换候选（第 4 序位）：Sirt2（class B、4 时点、PAP 百分位 15%——如实记录其健康海马 PAP 低于中位）。",
    "Ndrg2：保留观察（PAS 8bp 极佳，但校正后 QKI 证据消失、d3 覆盖比不一致、PAP 百分位 13.9%）。",
    "Cnp/Pea15a：hold（覆盖比降幅 <25%）；Kazn/Plec/Fam107a：low 置信（深度不足）；Kazn 另在健康海马 PAP 中 2/3 动物不可检出。",
    "",
    "## 各靶标设计区间（mm10，BED 0-based half-open，负链基因转录方向=基因组低→高坐标的反向）",
    "",
    "3'RACE 基因特异正向引物区 = 预测近端 PAS 的转录上游（负链=更高坐标）；common 扩增子/探针位于 shared 段；distal 扩增子/探针位于 distal 段且避开 shared 边界 ±100bp。",
    "",
]
for sym, role in targets:
    r, d, s = segs(sym)
    pas = int(r["proximal_pas_bed"])
    if r["strand"] == "-":
        fwd_zone = (pas + 100, min(s[1] - 50, pas + 600))
        common_amp = (pas + 60, s[1] - 50)
        distal_amp = (d[0] + 50, min(d[1] - 50, d[0] + 1200))
        fish_common = (pas + 60, s[1])
        fish_distal = (d[0] + 40, min(d[1] - 40, d[0] + 1000))
    else:
        fwd_zone = (max(s[0] + 50, pas - 600), pas - 100)
        common_amp = (s[0] + 50, pas - 60)
        distal_amp = (max(d[0] + 50, d[1] - 1200), d[1] - 50)
        fish_common = (s[0], pas - 60)
        fish_distal = (max(d[0], d[1] - 1000), d[1] - 40)
    L += [
        f"### {sym}（{role}）event {r['event_num']} — {r['chrom']} {r['strand']}",
        f"- 3'UTR 全长：{d[0] if r['strand']=='-' else s[0]}–{s[1] if r['strand']=='-' else d[1]}（DaPars2 BED）",
        f"- distal 段（仅长异构体）：**{r['chrom']}:{d[0]}-{d[1]}**（{d[1]-d[0]}bp）",
        f"- shared 段（全部异构体）：**{r['chrom']}:{s[0]}-{s[1]}**（{s[1]-s[0]}bp）",
        f"- 预测近端 PAS：{pas}（PolyASite 最近 {r['pas_dist_bp']}bp）",
        f"- 3'RACE 外侧引物区：{r['chrom']}:{fwd_zone[0]}-{fwd_zone[1]}；嵌套内侧引物区：紧邻 PAS 转录上游 200bp 内",
        f"- ddPCR/qPCR common 扩增子区：{r['chrom']}:{common_amp[0]}-{common_amp[1]}",
        f"- ddPCR/qPCR distal 扩增子区：{r['chrom']}:{distal_amp[0]}-{distal_amp[1]}",
        f"- smFISH 探针铺瓦区：common {r['chrom']}:{fish_common[0]}-{fish_common[1]}；distal {r['chrom']}:{fish_distal[0]}-{fish_distal[1]}",
        "",
    ]

L += [
    "## Gate 1 必须满足项核对（主靶 Atp2a2）",
    "",
    "| 条件 | 状态 |",
    "|---|---|",
    "| 链方向重建+单元测试通过（P1） | ✓ 15/15 |",
    "| BAM 无明显比对/注释伪影（P3，10 基因富 A 全 low、quickcheck 6/6） | ✓ |",
    "| 两个发现样本方向一致（d3 0.084/0.083；d7 0.085/0.086 覆盖比；PDUI 0.08/0.06、0.09/0.09） | ✓ |",
    "| distal 区可设计特异探针/引物（2,986bp distal，655bp shared，均 >150bp） | ✓ |",
    "| 预测 PAS 附近无严重内部引物风险（max A-run 3、A% 24） | ✓ |",
    "| d7 主时间点同向（覆盖比与 PDUI 同向缩短） | ✓ |",
    "",
    "加分项：多时点同向 ✓（d3+d7）；公共 PAS ✓（13bp）；GSE143531 基因级可测 ✓；校正后 distal RBP ✓（QKI FDR 1.04e-5）；第二方法直接覆盖 ✗（absent，如实登记）。",
    "",
    "## P4 公共数据辅助的如实记录",
    "",
    "1. **GSE143531（健康海马 PAP-TRAP，3 动物×2 分区×8 技术重复，已合并）**：全部主/备候选**可检出** ✓（Gate 1 加分项判定要求）。但相对全基因组分布，**候选未见 PAP 富集**：Atp2a2 位于 54% 百分位（中位）、Agpat3 52%、Aplp1 33%、Sirt2 15%、Plec 96%（唯一富集者，但 BAM 侧低置信）。结论：**v1 的 PAP 富集先验不跨区域/协议迁移，降级为\"区域依赖性待验证\"**；项目假设（卒中诱导的异构体分布改变）不依赖基线富集，将由脑片空间实验在皮层直接测量。",
    "2. **GSE330741（SN-MPRA）**：本地 651 基因/1,631 tile 计数表中 Glt1/Slc1a2、Sparc、Hsbp1 及全部候选**均不存在**（并有 Excel 日期损坏基因名如 \"9-Sep\"）——与论文描述的 focused tiling 文库矛盾，**manifest 核对仍未解决**（沿袭 v5.1 P4-04 BLOCKED）。阳性对照需从 Koester 2026（PMID 42094343）正式 supplementary 获取或改用文献元件；**不阻断 Gate 1**（MPRA 不在证据链内，仅背景）。",
    "",
    "## 下一步（进入湿实验 WP1）",
    "",
    "- 主靶 Atp2a2 3'RACE 引物对（外侧+嵌套）按上表区间交湿实验合作方用 Primer3 设计（Tm 58–62、产物 80–250bp、外显子跨剪接不适用——3'UTR 区连续，靠位置特异性）。",
    "- ddPCR common/distal 扩增子同法设计并用标准品验证线性范围后，再订 smFISH 探针组。",
    "- 动物模型尽量对齐 GSE238125（永久性远端 MCA 结扎+7 分钟双侧 CCA 阻断），d7 主时点，每组先导 4 只。",
]
open(os.path.join(OUT, "target_nomination_v1.md"), "w", encoding="utf-8").write("\n".join(L) + "\n")
print("written:", os.path.join(OUT, "target_nomination_v1.md"))
