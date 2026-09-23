# -*- coding: utf-8 -*-
"""Final Q&A extraction: PAP-member significant genes + triple-factor check
(localized + shortened + QKI-binding-loss)."""
import pandas as pd

OUT = "/mnt/d/stroke_apa/results"
m = pd.read_csv(f"{OUT}/p3_p4_event_annotation_matrix.tsv", sep="\t")

pap = m[m["axis_PAP"] == "yes"].copy()
print(f"PAP members among 977 sig genes: {len(pap)}")
print("top 20 by padj:")
top = pap.sort_values("padj_min").head(20)
print(" ".join(f"{r.symbol}({r.padj_min:.3f})" for _, r in top.iterrows()))

qki = m[m["L2_qki_lostseg"] == "yes"]
triple = qki[qki["axis_PAP"] == "yes"]
print(f"\ntriple-factor (QKI-loss + PAP-localized): {len(triple)}")
print(triple[["symbol", "contrasts", "direction", "dPDUI_min", "padj_min"]].to_string(index=False))

ef = m[m["axis_endfoot"] == "yes"]
print(f"\nendfoot-response members among sig: {len(ef)}")
print(" ".join(ef.sort_values("padj_min")["symbol"].head(20)))

wm = m[m["axis_zone_WM"] == "yes"]
triple_wm = qki[qki["axis_zone_WM"] == "yes"]
print(f"\nWM members: {len(wm)}; QKI-loss+WM: {len(triple_wm)} -> {' '.join(triple_wm['symbol'])}")
