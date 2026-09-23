# -*- coding: utf-8 -*-
"""Build FIMO input BEDs v2: SIMPLE ids (no '|' no '::') + mapping TSVs.
Outputs: results/fimo_seqs/{pas_events.bed,lost_segments.bed} + {pas_map.tsv,lost_map.tsv}"""
import numpy as np, pandas as pd

OUT = "/mnt/d/stroke_apa/results"
import os
os.makedirs(f"{OUT}/fimo_seqs", exist_ok=True)
mat = pd.read_csv(f"{OUT}/p2_full_pdui_matrix.tsv", sep="\t")
mat1 = mat.copy(); mat1.index = mat1.index + 1  # layer2 event ids are 1-based
ev = pd.read_csv(f"{OUT}/p2_sap_layer2_events.tsv", sep="\t")
sig = ev[ev["sig_joint"] == True].copy()
sig["Gene"] = sig["event"].astype(int).map(mat1["Gene"])
sig["symbol"] = sig["Gene"].astype(str).str.split("|").str[2]
sig["loci"] = sig["event"].astype(int).map(mat1["Loci"])
sig["prox"] = sig["event"].astype(int).map(mat1["Predicted_Proximal_APA"]).astype(float)
sp = sig["loci"].str.extract(r"(chr[^:]+):(\d+)-(\d+)")
sig["chr"], sig["utr_end"] = sp[0], sp[2].astype(int)
sig = sig.dropna(subset=["prox", "chr"])
sig["prox"] = np.floor(sig["prox"]).astype(int)

# PAS windows (one per event, union events)
pas = sig.drop_duplicates("event").copy().reset_index(drop=True)
pas["start"] = (pas["prox"] - 250).clip(lower=0)
pas["end"] = pas["prox"] + 250
pas["fid"] = ["pas%05d" % i for i in range(len(pas))]
pas[["chr", "start", "end", "fid"]].to_csv(f"{OUT}/fimo_seqs/pas_events.bed",
                                           sep="\t", header=False, index=False)
pas[["fid", "symbol", "contrast", "event", "dPDUI"]].to_csv(f"{OUT}/fimo_seqs/pas_map.tsv",
                                                            sep="\t", index=False)
print(f"[beds v2] pas windows: {len(pas)}")

# lost segments (shortening significant, non-trivial)
lost = sig[sig["dPDUI"] < 0].copy()
lost = lost[lost["utr_end"] > lost["prox"] + 10].reset_index(drop=True)
lost["fid"] = ["lost%05d" % i for i in range(len(lost))]
lost[["chr", "prox", "utr_end", "fid"]].to_csv(f"{OUT}/fimo_seqs/lost_segments.bed",
                                               sep="\t", header=False, index=False)
lost[["fid", "symbol", "contrast", "event", "dPDUI"]].to_csv(f"{OUT}/fimo_seqs/lost_map.tsv",
                                                             sep="\t", index=False)
print(f"[beds v2] lost segments: {len(lost)} ({(lost['utr_end']-lost['prox']).sum()/1e6:.1f} Mb)")
