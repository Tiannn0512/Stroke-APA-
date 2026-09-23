# -*- coding: utf-8 -*-
"""Gate G4 candidate cards: per B-class candidate, pull coordinates, QKI peak overlap,
frozen-family motif hits (p<1e-4), conservation bp, per-timepoint dPDUI.
Output: results/p4_candidate_cards.md (also printed)."""
import numpy as np, pandas as pd

OUT = "/mnt/d/stroke_apa/results"
cands = pd.read_csv(f"{OUT}/p4_top_b_class_candidates.tsv", sep="\t")
mat = pd.read_csv(f"{OUT}/p2_full_pdui_matrix.tsv", sep="\t")
mat1 = mat.copy(); mat1.index = mat1.index + 1
ev = pd.read_csv(f"{OUT}/p2_sap_layer2_events.tsv", sep="\t")
qki = pd.read_csv(f"{OUT}/p3_level2_qki_intersect.tsv", sep="\t")
cons = pd.read_csv("/home/taylor/stroke_APA_work/fimo_seq/lost_segments_conservation.tsv",
                   sep="\t", header=None, names=["seq_key", "call", "sum"])
lost_bed = pd.read_csv(f"{OUT}/fimo_seqs/lost_segments.bed", sep="\t", header=None,
                       names=["chr", "start", "end", "seq_key"])
fimo = pd.read_csv("/home/taylor/stroke_APA_work/fimo_seq/fimo_lost/fimo.tsv", sep="\t")
fimo["family"] = fimo["motif_id"].astype(str).str.split("_").str[0]
fimo["seq_key"] = fimo["sequence_name"].astype(str).str.split("_at_").str[0]
FROZEN = {"QKI","ELAVL1","PTBP1","MBNL2","NOVA1","NOVA2","TARDBP","HNRNPA2B1"}
TPS = ["day1","day3","day7","day21","day60"]

def event_coords(eid):
    r = mat1.loc[int(eid)]
    loci = r["Loci"]; prox = float(r["Predicted_Proximal_APA"])
    ch, s, e = loci.split(":")[0], int(loci.split(":")[1].split("-")[0]), int(loci.split(":")[1].split("-")[1])
    return ch, s, e, int(np.floor(prox))

md = "# Gate G4 候选卡（B 级，冻结规则自动选出 5/977）\n\n"
md += "> 每卡四要素：事件（mm10 坐标 + ΔPDUI 轨迹）/ regulator 及证据级别 / class 判定 / 三轴归属。表述按 L2 纪律：\"预测丢失一段 QKI 结合区\"。\n\n"
for _, c in cands.iterrows():
    eid = str(int(c["event"]))
    sym = c["symbol"]
    evr = ev[(ev["event"].astype(int) == int(eid))]
    ch, s, e, prox = event_coords(eid)
    sub = evr[["contrast","dPDUI","padj_joint"]].sort_values("contrast")
    traj = "; ".join(f"{r['contrast']} ΔPDUI={r['dPDUI']:+.3f} (padj {r['padj_joint']:.3f})"
                     for _, r in sub.iterrows())
    # QKI peaks overlapping this event's lost segment
    qp = qki[(qki["event"].astype(str) == eid)]
    qtxt = "; ".join(f"{r['peak_chr']}:{r['peak_start']}-{r['peak_end']}(FDR {r['peak_fdr']:.1e})"
                     for _, r in qp.iterrows()) or "n/a"
    # conservation
    cname = lost_bed[lost_bed["seq_key"].astype(str).str.startswith(sym + "|") &
                     lost_bed["seq_key"].astype(str).str.endswith("|" + eid)]
    cons_bp = ""
    for _, cr in cname.iterrows():
        hit = cons[cons["seq_key"] == cr["seq_key"]]
        if len(hit):
            cons_bp = f"{hit['call'].iloc[0]} (overlap {int(hit['sum'].iloc[0])}bp)"
    # frozen-family motif hits p<1e-4
    fh = fimo[(fimo["seq_key"].astype(str).str.endswith("|" + eid)) &
              (fimo["p-value"] < 1e-4) & (fimo["family"].isin(FROZEN))]
    fams = fh.groupby("family")["p-value"].min().sort_values()
    ftxt = ", ".join(f"{f}({p:.1e})" for f, p in fams.items()) or "n/a"
    # axis flags
    axes = []
    for ax, fp in [("PAP", "p3_m4_set1_pap_enriched.txt"), ("endfoot", "p3_m4_set2_endfoot_de.txt"),
                   ("zone_cortex", "p3_m4_set3_stroke_responsive_cortex.txt"),
                   ("zone_WM", "p3_m4_set3_stroke_responsive_whitematter.txt")]:
        sset = set(x.strip() for x in open(f"{OUT}/{fp}") if x.strip())
        if sym in sset:
            axes.append(ax)
    md += f"## {sym}（event {eid}）\n"
    md += f"- **事件**：{ch}:{s}-{e}，推断近端 PAS = {prox}；丢失段 = {prox}-{e}\n"
    md += f"- **ΔPDUI 轨迹**：{traj}\n"
    md += f"- **L2（QKI）**：丢失段与 QKI peak 重叠：{qtxt}\n"
    md += f"- **L3（冻结家族命中，p<1e-4）**：{ftxt}\n"
    md += f"- **L4（保守性）**：{cons_bp or 'n/a'}\n"
    md += f"- **class**：B（CLIP+motif+保守）｜**三轴归属**：{', '.join(axes) or 'none'}\n\n"

with open(f"{OUT}/p4_candidate_cards.md", "w") as f:
    f.write(md)
print(md[:1500])
print("... [full cards written to p4_candidate_cards.md]")
