# -*- coding: utf-8 -*-
"""Fix 579->471 precedence bug in G2 evidence extras (registered correction 2026-09-20)."""
import pandas as pd
OUT = "/mnt/d/stroke_apa/results"
ev = pd.read_csv(f"{OUT}/p2_sap_layer2_events.tsv", sep="\t")
sig = ev[ev["sig_joint"] == True]
w = sig.pivot_table(index="event", columns="contrast", values="dPDUI", aggfunc="first")
n2 = w.apply(lambda r: sum(1 for v in r if pd.notna(v)) >= 2 and
             ((r.dropna() < 0).all() or (r.dropna() > 0).all()), axis=1).sum()
print(f"corrected >=2-timepoint concordant events: {int(n2)} (was 579 via precedence bug)")
