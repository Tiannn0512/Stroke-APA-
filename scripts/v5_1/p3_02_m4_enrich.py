# -*- coding: utf-8 -*-
"""P3-02 L2 intersect + M4 three-axis Fisher enrichment.
1) QKI CLIP peaks x lost-distal-segment of significant APA events
   (lost segment = inferred proximal APA -> annotated UTR end; shortening = lower distal usage).
2) Union significant events -> gene set -> Fisher exact vs 4 annotation axes
   (PAP_localized / endfoot_stroke_responsive / zone cortex / zone white matter).
Run: Windows python.
"""
import numpy as np, pandas as pd
from scipy.stats import fisher_exact

OUT = r"D:\stroke_apa\results"

# ---------- inputs ----------
mat = pd.read_csv(f"{OUT}\\p2_full_pdui_matrix.tsv", sep="\t")
# R layer2 wrote event ids as 1-BASED row numbers; re-index matrix to 1-based to align
mat1 = mat.copy(); mat1.index = mat1.index + 1
ev = pd.read_csv(f"{OUT}\\p2_sap_layer2_events.tsv", sep="\t")
sig = ev[ev["sig_joint"] == True].copy()
sig["Gene"] = sig["event"].astype(int).map(mat1["Gene"])
sig["symbol"] = sig["Gene"].astype(str).str.split("|").str[2]
sig["loci"] = sig["event"].astype(int).map(mat1["Loci"])
sig["prox"] = sig["event"].astype(int).map(mat1["Predicted_Proximal_APA"]).astype(float)

peaks = pd.read_csv(f"{OUT}\\p3_level2_qki_clip_peaks.bed", sep="\t", header=None,
                    names=["chr", "start", "end", "gene", "fdr", "strand"])

# parse event loci
sp = sig["loci"].str.extract(r"(chr[^:]+):(\d+)-(\d+)")
sig["chr"], sig["utr_start"], sig["utr_end"] = sp[0], sp[1].astype(int), sp[2].astype(int)

# ---------- 1. L2: lost distal segment x QKI peaks ----------
# lost segment for a shortening event = [proximal_APA, utr_end]; for lengthening the same segment gains usage.
seg = sig.dropna(subset=["prox"]).copy()
seg["seg_start"] = np.floor(seg["prox"]).astype(int)
seg["seg_end"] = seg["utr_end"]
rows = []
peaks_by_chr = {c: g for c, g in peaks.groupby("chr")}
for _, r in seg.iterrows():
    pk = peaks_by_chr.get(r["chr"])
    if pk is None:
        continue
    hit = pk[(pk["end"] > r["seg_start"]) & (pk["start"] < r["seg_end"])]
    for _, p in hit.iterrows():
        rows.append(dict(event=int(r["event"]), contrast=r["contrast"], symbol=r["symbol"],
                         dPDUI=r["dPDUI"], padj_joint=r["padj_joint"], chr=r["chr"],
                         seg_start=r["seg_start"], seg_end=r["seg_end"],
                         peak_chr=p["chr"], peak_start=p["start"], peak_end=p["end"],
                         peak_gene=p["gene"], peak_fdr=p["fdr"]))
inter = pd.DataFrame(rows)
inter.to_csv(f"{OUT}\\p3_level2_qki_intersect.tsv", sep="\t", index=False)

n_events_total = seg["event"].nunique()
n_events_hit = inter["event"].nunique()
n_short_hit = inter[inter["dPDUI"] < 0]["event"].nunique()
print(f"[L2] events with lost/gained segment: {n_events_total}; "
      f"events overlapping >=1 QKI peak: {n_events_hit} ({n_events_hit/n_events_total*100:.1f}%); "
      f"shortening-hit: {n_short_hit}", flush=True)
print(f"[L2] peak-side genes hit: {inter['peak_gene'].nunique()}; "
      f"event-side genes: {inter['symbol'].nunique()}", flush=True)

# background rate: all testable events (4/4) same overlap test -> enrichment of QKI overlap among sig
ev_all = ev.copy()
ev_all["Gene"] = ev_all["event"].astype(int).map(mat1["Gene"])
ev_all["symbol"] = ev_all["Gene"].astype(str).str.split("|").str[2]
ev_all["loci"] = ev_all["event"].astype(int).map(mat1["Loci"])
ev_all["prox"] = ev_all["event"].astype(int).map(mat1["Predicted_Proximal_APA"]).astype(float)
sp2 = ev_all["loci"].str.extract(r"(chr[^:]+):(\d+)-(\d+)")
ev_all["chr"], ev_all["utr_end"] = sp2[0], sp2[2].astype(int)
seg_all = ev_all.dropna(subset=["prox"])
hit_all = set()
for c, g in seg_all.groupby("chr"):
    pk = peaks_by_chr.get(c)
    if pk is None:
        continue
    starts = pk["start"].values; ends = pk["end"].values
    for _, r in g.iterrows():
        s0, e0 = int(np.floor(r["prox"])), r["utr_end"]
        if ((ends > s0) & (starts < e0)).any():
            hit_all.add(int(r["event"]))
sig_ids = set(sig["event"].astype(int))
sig_in = len(sig_ids & hit_all); sig_out = len(sig_ids) - sig_in
bck_ids = set(seg_all["event"].astype(int)) - sig_ids
bck_in = len(bck_ids & hit_all); bck_out = len(bck_ids) - bck_in
odds, pval = fisher_exact([[sig_in, sig_out], [bck_in, bck_out]])
print(f"[L2 enrichment] sig hit {sig_in}/{len(sig_ids)} vs background {bck_in}/{len(bck_ids)}; "
      f"OR={odds:.2f} p={pval:.2e}", flush=True)
with open(f"{OUT}\\p3_level2_qki_intersect_summary.tsv", "w") as f:
    f.write("metric\tvalue\n")
    f.write(f"sig_events_with_segment\t{n_events_total}\n")
    f.write(f"sig_events_overlapping_qki_peak\t{n_events_hit}\n")
    f.write(f"shortening_sig_events_overlapping\t{n_short_hit}\n")
    f.write(f"enrichment_sig_hit\t{sig_in}/{len(sig_ids)}\n")
    f.write(f"enrichment_background_hit\t{bck_in}/{len(bck_ids)}\n")
    f.write(f"fisher_OR\t{odds:.3f}\n"); f.write(f"fisher_p\t{pval:.3e}\n")

# ---------- 2. M4 three-axis Fisher ----------
universe = set(mat["Gene"].astype(str).str.split("|").str[2].dropna())
sig_genes = set(sig["symbol"].dropna())
sig_short = set(sig[sig["dPDUI"] < 0]["symbol"].dropna())
axes = {
    "PAP_localized": f"{OUT}\\p3_m4_set1_pap_enriched.txt",
    "endfoot_stroke_responsive": f"{OUT}\\p3_m4_set2_endfoot_de.txt",
    "zone_cortex": f"{OUT}\\p3_m4_set3_stroke_responsive_cortex.txt",
    "zone_whitematter": f"{OUT}\\p3_m4_set3_stroke_responsive_whitematter.txt",
}
rows = []
for name, fp in axes.items():
    axset = set(x.strip() for x in open(fp) if x.strip())
    for tag, q in [("all_sig", sig_genes), ("shortening_sig", sig_short)]:
        a = len(q & axset); b = len(q - axset)
        c = len(universe - q) - (len((universe - q) & axset)); d = len((universe - q) & axset)
        # fix: c = universe&axset - q∩axset
        c = len((universe & axset) - q); d = len((universe - q) - axset)
        odds, pval = fisher_exact([[a, b], [c, d]])
        rows.append(dict(axis=name, gene_set=tag, genes=len(q), overlap=a,
                         odds_ratio=round(odds, 3), p_value=f"{pval:.3e}",
                         pct_in_axis=round(a / len(q) * 100, 1)))
enr = pd.DataFrame(rows)
enr.to_csv(f"{OUT}\\p3_m4_enrichment.tsv", sep="\t", index=False)
print("\n[M4 enrichment]")
print(enr.to_string(index=False), flush=True)
print("[DONE] P3-02 + M4 enrichment", flush=True)
