# -*- coding: utf-8 -*-
"""Final answer gene lists."""
import pandas as pd
m = pd.read_csv("/mnt/d/stroke_apa/results/p3_p4_event_annotation_matrix.tsv", sep="\t")
g = m.groupby("symbol").agg(n_events=("event", "count"), min_padj=("padj_min", "min"),
                            direction=("direction", "first"), cls=("class", "first"),
                            L2=("L2_qki_lostseg", "first"), PAP=("axis_PAP", "first")).reset_index()
g = g.sort_values("min_padj")
print("unique genes:", len(g))
print("--- top 30 by padj ---")
print(" ".join(g.head(30)["symbol"]))
print("--- QKI-overlap genes (11 events) ---")
print(" ".join(m[m.L2_qki_lostseg == "yes"]["symbol"].unique()))
print("--- B-class 5 ---")
print(" ".join(m[m["class"] == "B"]["symbol"].unique()))
print("--- PAP-axis members among sig genes ---")
print((g["PAP"] == "yes").sum())
g.to_csv("/mnt/d/stroke_apa/results/final_answer_gene_list.tsv", sep="\t", index=False)
print("[saved] results/final_answer_gene_list.tsv")
