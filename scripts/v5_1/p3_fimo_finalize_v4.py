# -*- coding: utf-8 -*-
"""P3-03 finalize v4 (final design):
Family-level motif enrichment in lost segments vs UTR-derived matched null
(non-significant events' distal segments, length-matched).
Statistic: hits per Mb (density) per family; Poisson test with exposure = bases.
Outputs: results/p3_level3_family_enrichment.tsv (primary) + p3_level3_motif.tsv (descriptive, kept from v3).
"""
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

lost = load_hits(f"{BASE}/fimo_lost/fimo.tsv")
null = load_hits(f"{BASE}/fimo_null_utr/fimo.tsv")

def bases(bed):
    b = pd.read_csv(bed, sep="\t", header=None, names=["chr", "s", "e", "name"])
    return float((b["e"] - b["s"]).sum())
bases_lost = bases(f"{OUT}/fimo_seqs/lost_segments.bed")
bases_null = bases(f"{OUT}/fimo_seqs/null_utr.bed")
print(f"hits: lost={len(lost)} null_utr={len(null)}; bases: lost={bases_lost/1e6:.2f}Mb "
      f"null_utr={bases_null/1e6:.2f}Mb", flush=True)

try:
    from scipy.stats import binomtest
    binom_pfun = lambda k, n, p: binomtest(k, n, p).pvalue
except ImportError:
    from scipy.stats import binom_test
    binom_pfun = lambda k, n, p: binom_test(k, n=n, p=p)
from scipy.special import gammaincc
def poisson_sf(k, lam):  # P(X >= k) under Poisson(lam)
    return float(gammaincc(k, lam))

n_lost_regions = len(pd.read_csv(f"{OUT}/fimo_seqs/lost_segments.bed", sep="\t", header=None))
n_null_regions = len(pd.read_csv(f"{OUT}/fimo_seqs/null_utr.bed", sep="\t", header=None))
rows = []
fams = sorted(null["family"].unique())
for fam_ in fams:
    hl = len(lost[lost["family"] == fam_])
    hn = len(null[null["family"] == fam_])
    dens_l = hl / (bases_lost / 1e6)
    dens_n = hn / (bases_null / 1e6)
    lam = dens_n * (bases_lost / 1e6)
    pv = poisson_sf(hl, lam)
    rows.append(dict(family=fam_, hits_lost=hl, hits_null_utr=hn,
                     density_lost_perMb=round(dens_l, 1), density_null_perMb=round(dens_n, 1),
                     enrichment=round(dens_l / max(dens_n, 1e-9), 2), poisson_p=f"{pv:.2e}"))
enr = pd.DataFrame(rows)
enr["q_bh_families"] = bh([float(x) for x in enr["poisson_p"]])
enr = enr.sort_values("poisson_p")
enr.to_csv(f"{OUT}/p3_level3_family_enrichment.tsv", sep="\t", index=False)
print(enr.to_string(index=False), flush=True)
sig = enr[enr["q_bh_families"] < 0.1]
print(f"\nfamilies with density enrichment q<0.1: {len(sig)} -> {', '.join(sig['family'])}", flush=True)
print("[DONE] finalize v4", flush=True)
