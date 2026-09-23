# -*- coding: utf-8 -*-
"""P3-04: per-event association/annotation matrix for the 977 significant events (runs in WSL).
Columns: symbol, contrasts, direction, max|dPDUI|, min padj, L1, L2_qki, L3_motif(frozen),
L4_conserved, axis flags. Output: results/p3_p4_event_annotation_matrix.tsv
Run: /home/taylor/miniconda3/envs/dapars2/bin/python p3_04_matrix.py
"""
import numpy as np, pandas as pd

OUT = "/mnt/d/stroke_apa/results"
FROZEN_FAMS = {"QKI", "ELAVL1", "PTBP1", "MBNL2", "NOVA1", "NOVA2", "TARDBP", "HNRNPA2B1"}
L1_DIRECT = {"Ccnd1", "Qk"}
L1_GW = {"PABPN1", "NUDT21", "QKI", "ELAVL1", "PTBP1", "NOVA1", "NOVA2", "MBNL1", "MBNL2", "TARDBP"}

ev = pd.read_csv(f"{OUT}/p2_sap_layer2_events.tsv", sep="\t")
sig = ev[ev["sig_joint"] == True].copy()
mat = pd.read_csv(f"{OUT}/p2_full_pdui_matrix.tsv", sep="\t")
mat1 = mat.copy(); mat1.index = mat1.index + 1
sig["symbol"] = sig["event"].astype(int).map(mat1["Gene"]).astype(str).str.split("|").str[2]

g = sig.groupby("event")
m = g.agg(symbol=("symbol", "first"),
          contrasts=("contrast", lambda x: ",".join(sorted(set(x)))),
          n_contrasts=("contrast", "nunique"),
          dPDUI_min=("dPDUI", "min"), dPDUI_max=("dPDUI", "max"),
          padj_min=("padj_joint", "min")).reset_index()
m["direction"] = np.where(m["dPDUI_max"] < 0, "shortening",
                  np.where(m["dPDUI_min"] > 0, "lengthening", "mixed"))

pairs = pd.read_csv(f"{OUT}/p3_level1_pairs.tsv", sep="\t")
def l1(row):
    if row["symbol"] in L1_DIRECT:
        return "direct_gene"
    hit = pairs[pairs["regulator"].str.upper() == row["symbol"].upper()]
    if len(hit) and hit["regulator"].str.upper().iloc[0] in L1_GW:
        return "regulator_GW"
    return "none"
m["L1"] = m.apply(l1, axis=1)

qki = pd.read_csv(f"{OUT}/p3_level2_qki_intersect.tsv", sep="\t")
m["L2_qki_lostseg"] = m["event"].isin(qki["event"].astype(int)).map({True: "yes", False: "no"})

fimo = pd.read_csv("/home/taylor/stroke_APA_work/fimo_seq/fimo_lost/fimo.tsv", sep="\t")
fimo["family"] = fimo["motif_id"].astype(str).str.split("_").str[0]
fimo["seq_key"] = fimo["sequence_name"].astype(str).str.split("_at_").str[0]
fimo["event_id"] = fimo["seq_key"].astype(str).str.split("|").str[2]
l3_hits = set(fimo[(fimo["p-value"] < 1e-4) & (fimo["family"].isin(FROZEN_FAMS))]["event_id"])
m["L3_motif_frozen"] = m["event"].astype(int).astype(str).map(
    {e: ("yes" if e in l3_hits else "no") for e in set(m["event"].astype(int).astype(str))})

cons = pd.read_csv("/home/taylor/stroke_APA_work/fimo_seq/lost_segments_conservation.tsv",
                   sep="\t", header=None, names=["seq_key", "call", "sum"], usecols=[0, 1])
cons["event_id"] = cons["seq_key"].astype(str).str.split("|").str[2]
cons_map = dict(zip(cons["event_id"], cons["call"]))
m["L4_conserved"] = m["event"].astype(int).astype(str).map(cons_map)

for ax, fp in [("axis_PAP", "p3_m4_set1_pap_enriched.txt"),
               ("axis_endfoot", "p3_m4_set2_endfoot_de.txt"),
               ("axis_zone_cortex", "p3_m4_set3_stroke_responsive_cortex.txt"),
               ("axis_zone_WM", "p3_m4_set3_stroke_responsive_whitematter.txt")]:
    s = set(x.strip() for x in open(f"{OUT}/{fp}") if x.strip())
    m[ax] = m["symbol"].isin(s).map({True: "yes", False: "no"})

# class per frozen P4-03 rules
def cls(r):
    if r["L3_motif_frozen"] == "yes":
        if r["L2_qki_lostseg"] == "yes" and r["L4_conserved"] == "CONSERVED":
            return "B"
        if r["L4_conserved"] == "CONSERVED":
            return "C"
        return "D"
    return "none"
m["class"] = m.apply(cls, axis=1)

m.to_csv(f"{OUT}/p3_p4_event_annotation_matrix.tsv", sep="\t", index=False)
print(f"matrix: {len(m)} events")
print("L1:", m["L1"].value_counts().to_dict())
print("L2:", m["L2_qki_lostseg"].value_counts().to_dict())
print("L3:", m["L3_motif_frozen"].value_counts().to_dict())
print("L4:", m["L4_conserved"].fillna("missing").value_counts().to_dict())
print("class:", m["class"].value_counts().to_dict())
print("direction:", m["direction"].value_counts().to_dict())
print("axes:", {c: (m[c] == "yes").sum() for c in m.columns if c.startswith("axis_")})
print("[DONE] P3-04", flush=True)
