# -*- coding: utf-8 -*-
"""P3-03 finalize v3: mapping-table joins (simple ids), family rates vs matched null,
manual BH (old scipy), descriptive event x family table.
Outputs: results/p3_level3_motif.tsv + p3_level3_family_enrichment.tsv"""
import re
import numpy as np, pandas as pd

BASE = "/home/taylor/stroke_APA_work/fimo_seq"
OUT = "/mnt/d/stroke_apa/results"

def load_hits(path):
    d = pd.read_csv(path, sep="\t", comment="#", low_memory=False)
    d = d[d["motif_id"].notna()].copy()
    d["family"] = d["motif_id"].astype(str).str.split("_").str[0]
    return d

def bh(p):
    p = np.asarray(p, float); n = len(p)
    o = np.argsort(p)
    q = p[o] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[o] = np.clip(q, 0, 1)
    return out

pas = load_hits(f"{BASE}/fimo_pas/fimo.tsv")
lost = load_hits(f"{BASE}/fimo_lost/fimo.tsv")
null = load_hits(f"{BASE}/fimo_null/fimo.tsv")
print(f"hits: pas={len(pas)} lost={len(lost)} null={len(null)}", flush=True)

# sequence ids: 'symbol|contrast|eventid_at_chrN:start-end' -> key = part before '_at_'
def key_of(d):
    return d["sequence_name"].astype(str).str.split("_at_").str[0]

pas["seq_key"] = key_of(pas)
lost["seq_key"] = key_of(lost)
null["seq_key"] = key_of(null)

lost_bed = pd.read_csv(f"{OUT}/fimo_seqs/lost_segments.bed", sep="\t", header=None,
                       names=["chr", "start", "end", "seq_id"])
pas_bed = pd.read_csv(f"{OUT}/fimo_seqs/pas_events.bed", sep="\t", header=None,
                      names=["chr", "start", "end", "seq_id"])

lost_f = lost.groupby(["seq_key", "family"], as_index=False)["p-value"].min()
lost_f = lost_f.merge(lost_bed, left_on="seq_key", right_on="seq_id", how="left")
lost_f["symbol"] = lost_f["seq_id"].astype(str).str.split("|").str[0]
lost_f["contrast"] = lost_f["seq_id"].astype(str).str.split("|").str[1]
lost_f["event"] = lost_f["seq_id"].astype(str).str.split("|").str[2]
pas_f = pas.groupby(["seq_key", "family"], as_index=False)["p-value"].min()
pas_f = pas_f.merge(pas_bed, left_on="seq_key", right_on="seq_id", how="left")
pas_f["symbol"] = pas_f["seq_id"].astype(str).str.split("|").str[0]
pas_f["event"] = pas_f["seq_id"].astype(str).str.split("|").str[2]

summary = lost_f[["seq_id", "symbol", "contrast", "event", "family", "p-value"]].copy()
summary = summary.rename(columns={"p-value": "lost_best_p"})
pas_min = pas_f.groupby(["event", "family"], as_index=False)["p-value"].min().rename(
    columns={"p-value": "pas_window_best_p"})
summary = summary.merge(pas_min, on=["event", "family"], how="left")
summary["interpretation"] = "lost_segment_of_shortened_event"
summary.to_csv(f"{OUT}/p3_level3_motif.tsv", sep="\t", index=False)
n_ev_hit = summary["event"].nunique()
n_matched = summary["event"].notna().sum()
print(f"[level3] event-family lost rows: {len(summary)}; matched: {n_matched}; "
      f"events with >=1 lost hit: {n_ev_hit}", flush=True)
print(summary.groupby("family")["event"].nunique().sort_values(ascending=False).head(8).to_string(),
      flush=True)

# ---- family-level rate test vs matched null ----
n_lost = lost_bed["seq_id"].nunique()
n_null = len(pd.read_csv(f"{OUT}/fimo_seqs/null_segments.bed", sep="\t", header=None))
try:
    from scipy.stats import binomtest
    binom_pfun = lambda k, n, p: binomtest(k, n, p).pvalue
except ImportError:
    from scipy.stats import binom_test
    binom_pfun = lambda k, n, p: binom_test(k, n=n, p=p)
rows = []
fams = sorted(null["family"].unique())
for fam_ in fams:
    hit_lost = lost[lost["family"] == fam_]["sequence_name"].nunique()
    hit_null = null[null["family"] == fam_]["sequence_name"].nunique()
    p_null = hit_null / max(n_null, 1)
    rate_lost = hit_lost / max(n_lost, 1)
    pv = binom_pfun(hit_lost, n_lost, p_null)
    rows.append(dict(family=fam_, lost_hit_regions=hit_lost, null_hit_regions=hit_null,
                     rate_lost=round(rate_lost, 4), rate_null=round(p_null, 4),
                     enrichment=round(rate_lost / max(p_null, 1e-9), 2), binom_p=f"{pv:.2e}"))
enr = pd.DataFrame(rows)
enr["q_bh_families"] = bh([float(x) for x in enr["binom_p"]])
enr = enr.sort_values("enrichment", ascending=False)
enr.to_csv(f"{OUT}/p3_level3_family_enrichment.tsv", sep="\t", index=False)
print("\n[family enrichment vs matched null]")
print(enr.to_string(index=False), flush=True)
print("[DONE] P3-03 finalize v3", flush=True)
