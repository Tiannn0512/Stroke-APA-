# -*- coding: utf-8 -*-
"""AUDIT PASS 1: verify every quantitative claim in FINAL_REPORT against artifact files."""
import numpy as np, pandas as pd

OUT = "/mnt/d/stroke_apa/results"
ok, fail = [], []
def chk(name, expected, actual):
    (ok if str(expected) == str(actual) else fail).append(f"{name}: report={expected} actual={actual}")

# --- P2 layer ---
m = pd.read_csv(f"{OUT}/p2_full_pdui_matrix.tsv", sep="\t")
chk("events_total", 9693, len(m))
qc = pd.read_csv(f"{OUT}/p2_full_merge_qc.tsv", sep="\t")
qv = dict(zip(qc.metric, qc.value))
chk("sham_r", 0.9257, round(qv["sham1_sham2_pearson_r"], 4))
for tp, ratio in [("day1",4.84),("day3",0.73),("day7",7.23),("day21",6.58),("day60",4.58)]:
    chk(f"ratio_{tp}", ratio, qv[f"{tp}_ratio"])

ev = pd.read_csv(f"{OUT}/p2_sap_layer2_events.tsv", sep="\t")
sig = ev[ev["sig_joint"] == True]
cnt = sig["contrast"].value_counts()
for tp, n in [("day1",292),("day3",159),("day7",435),("day21",496),("day60",392)]:
    chk(f"sig_{tp}", n, int(cnt[tp]))
chk("union_sig", 977, sig["event"].nunique())
w = sig.pivot_table(index="event", columns="contrast", values="dPDUI", aggfunc="first")
n2 = w.apply(lambda r: sum(1 for v in r if pd.notna(v)) >= 2 and
             ((r.dropna() < 0).all() or (r.dropna() > 0).all()), axis=1).sum()
chk("ge2tp_concordant", 579, int(n2))
om = pd.read_csv(f"{OUT}/p2_sap_layer2_omnibus.tsv", sep="\t")
chk("omnibus", "4910/5923", f"{int((om['padj_omnibus']<0.1).sum())}/{len(om)}")

# --- P3 ---
qsum = pd.read_csv(f"{OUT}/p3_level2_qki_intersect_summary.tsv", sep="\t")
qs = dict(zip(qsum.metric, qsum.value))
chk("L2_events_hit", 11, int(qs["sig_events_overlapping_qki_peak"]))
chk("L2_OR", 2.71, round(float(qs["fisher_OR"]), 2))
chk("L2_p", "1.25e-02", qs["fisher_p"])

fam = pd.read_csv(f"{OUT}/p3_level3_family_enrichment.tsv", sep="\t")
chk("L3_enrich_range", "0.89-1.07", f"{fam['enrichment'].min():.2f}-{fam['enrichment'].max():.2f}")
chk("L3_fams_pass_q", 0, int((fam["q_bh_families"] < 0.1).sum()))

cons = pd.read_csv("/home/taylor/stroke_APA_work/fimo_seq/lost_segments_conservation.tsv",
                   sep="\t", header=None, names=["k","call","sum"])
chk("L4_conserved", "1004/1541", f"{(cons['call']=='CONSERVED').sum()}/{len(cons)}")

gs = pd.read_csv(f"{OUT}/p3_m4_gsea.tsv", sep="\t")
pap = gs[gs["axis"]=="PAP_localized"].set_index("contrast")
for tp, v in [("day7","9.99900009999e-05"),("day3","0.00182386659718603"),("day1","0.0440955904409559")]:
    chk(f"gsea_PAP_{tp}", v, str(pap.loc[tp,"padj"]))
chk("gsea_omnibus_PAP", "0.0041995800419958", str(pap.loc["omnibus_directed","padj"]))
chk("gsea_endfoot_d3_pos", "+1.28", f"{gs[(gs.axis=='endfoot_stroke_responsive')&(gs.contrast=='day3')]['NES'].iloc[0]:+.2f}")

st = pd.read_csv(f"{OUT}/p2_second_tool_consistency.tsv", sep="\t")
chk("salmon_rho_range", "0.041-0.139", f"{st['spearman_rho'].min():.3f}-{st['spearman_rho'].max():.3f}")
chk("salmon_all_sig", 5, int(len(st)))

# --- P4 ---
mm = pd.read_csv(f"{OUT}/p3_p4_event_annotation_matrix.tsv", sep="\t")
chk("matrix_events", 977, len(mm))
chk("class_B5_C475_D250_none247", "5/475/250/247",
    f"{(mm['class']=='B').sum()}/{(mm['class']=='C').sum()}/{(mm['class']=='D').sum()}/{(mm['class']=='none').sum()}")
chk("direction_811_153_13", "811/153/13",
    f"{(mm['direction']=='shortening').sum()}/{(mm['direction']=='lengthening').sum()}/{(mm['direction']=='mixed').sum()}")
for ax, n in [("axis_PAP",161),("axis_endfoot",34),("axis_zone_cortex",241),("axis_zone_WM",429)]:
    chk(f"axis_{ax}", n, int((mm[ax]=="yes").sum()))
chk("Atp2a2_class_D", "D", mm[mm["symbol"]=="Atp2a2"]["class"].iloc[0])

# --- P1 (from committed tables) ---
cat = pd.read_csv(f"{OUT}/p1_regulator_catalog.tsv", sep="\t")
chk("P1_31_regulators", 31, len(cat))
chk("P1_concordant_17", 17, int((cat["direction_consistency"]=="concordant").sum()))

# --- G2(d) ---
chk("top40_concordance", "13/27", "13/27")  # from p2_second_tool_top40_signs.tsv
sg = pd.read_csv(f"{OUT}/p2_second_tool_top40_signs.tsv", sep="\t")
chk("top40_recount", "13/27", f"{(sg['call']=='agree').sum()}/{len(sg)}")

print(f"PASS: {len(ok)}  FAIL: {len(fail)}")
for f in fail:
    print("  FAIL:", f)
for o in ok:
    print("  ok:", o)
