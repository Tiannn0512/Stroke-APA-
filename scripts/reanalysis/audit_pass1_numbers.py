#!/usr/bin/env python3
"""AUDIT PASS 1 — independent recomputation of every load-bearing number.

Recomputes from primary sources (depth tsvs, quant.sf, fimo.tsv, atlas.gz,
intersections, conservation tsv, sap/pdui inputs) and compares against the
values quoted in the nomination doc / rebuild notes / review table.
Prints PASS/FAIL per check. Exit 1 on any FAIL.
"""
import bisect, csv, gzip, math, os, sys

BASE = r"D:\stroke_apa_reanalysis"
R = os.path.join(BASE, "results", "02_candidate_rebuild")
S1 = os.path.join(BASE, "results", "01_strand_audit")
BR = os.path.join(BASE, "results", "03_bam_review")
IL = os.path.join(BASE, "input_links")
FAILS = []

def check(label, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {label}" + (f"  | {detail}" if detail and not ok else ""))
    if not ok:
        FAILS.append(label)

def rd(p):
    return list(csv.DictReader(open(p, encoding="utf-8"), delimiter="\t"))

# ---------- 1. segments / strand ----------
ev = rd(os.path.join(IL, "dapars2_event_coordinates.tsv"))
seg = rd(os.path.join(S1, "strand_aware_event_segments.tsv"))
check("events = 9693", len(ev) == 9693, str(len(ev)))
plus = sum(1 for r in ev if r["strand"] == "+")
check("strand +4913/-4780", plus == 4913 and len(ev) - plus == 4780)
ok = True
for r in seg:
    d0, d1 = int(r["distal_start"]), int(r["distal_end"])
    s0, s1_ = int(r["shared_start"]), int(r["shared_end"])
    us, ue = int(r["utr_start"]), int(r["utr_end"])
    if not (0 < d1 - d0 and 0 < s1_ - s0): ok = False; break
    if min(d1, s1_) - max(d0, s0) > 0: ok = False; break
    if (d1 - d0) + (s1_ - s0) != ue - us: ok = False; break
check("segments: positive/adjacent/sum for all 9693", ok and len(seg) == 9693)
a = next(r for r in seg if r["gene_symbol"] == "Atp2a2")
check("Atp2a2 anchors distal 122453512-122456498 / shared 122456498-122457153",
      (a["distal_start"], a["distal_end"], a["shared_start"], a["shared_end"]) ==
      ("122453512", "122456498", "122456498", "122457153"))

# ---------- 2. coverage ratios recomputed from depth ----------
S6 = ["sham1", "sham2", "day3_rep1", "day3_rep2", "day7_rep1", "day7_rep2"]
regions = {}
for line in open(os.path.join(BR, "regions.bed"), encoding="utf-8"):
    p = line.rstrip("\n").split("\t")
    regions[p[3].split("|")[0]] = (p[0], int(p[1]) + 5000, int(p[2]) - 5000)
dep = {g: {s: {} for s in S6} for g in regions}
for s in S6:
    for line in open(os.path.join(BR, f"depth_{s}.tsv"), encoding="utf-8"):
        c, pos, v = line.rstrip("\n").split("\t")
        pos = int(pos)
        for g, (cc, us, ue) in regions.items():
            if c == cc and us <= pos < ue:
                dep[g][s][pos] = int(v)
rank = {r["gene_symbol"]: r for r in rd(os.path.join(R, "candidate_ranking_v2.tsv"))}
def mdepth(g, s, side):
    d = list(map(int, rank[g][side].split("-")))
    vals = [dep[g][s].get(p, 0) for p in range(d[0], d[1])]
    return sum(vals) / len(vals)
exp = {"Atp2a2": [0.334, 0.270, 0.086, 0.083, 0.086, 0.085],
       "Agpat3": [0.555, 0.387, 0.195, 0.206, 0.095, 0.072],
       "Ndrg2":  [0.544, 0.541, 0.769, 0.430, 0.378, 0.352]}
for g, expv in exp.items():
    got = [round(mdepth(g, s, "distal") / mdepth(g, s, "shared"), 3) for s in S6]
    check(f"coverage ratios {g} recomputed", got == expv, f"got {got} exp {expv}")
sham_m = (mdepth("Atp2a2", "sham1", "distal") / mdepth("Atp2a2", "sham1", "shared") +
          mdepth("Atp2a2", "sham2", "distal") / mdepth("Atp2a2", "sham2", "shared")) / 2
d7_m = (mdepth("Atp2a2", "day7_rep1", "distal") / mdepth("Atp2a2", "day7_rep1", "shared") +
        mdepth("Atp2a2", "day7_rep2", "distal") / mdepth("Atp2a2", "day7_rep2", "shared")) / 2
drop = round(100 * (sham_m - d7_m) / sham_m, 1)
check("Atp2a2 d7 drop 71.7%", abs(drop - 71.7) < 0.15, str(drop))

# ---------- 3. salmon baseline recomputed ----------
txs_atp = {"ENSMUST00000031423.9", "ENSMUST00000177974.7", "ENSMUST00000179939.7"}
txs_agp = {"ENSMUST00000001240.11", "ENSMUST00000105389.7", "ENSMUST00000166360.7",
           "ENSMUST00000105390.7", "ENSMUST00000105388.7", "ENSMUST00000105387.7",
           "ENSMUST00000219932.1"}
exp_sal = {"Atp2a2": [0.719, 0.541, 0.13, 0.112, 0.039, 0.039, 0.066, 0.059, 0.275, 0.127],
           "Agpat3": [0.096, 0.127, 0.067, 0.042, 0.041, 0.034, 0.045, 0.04, 0.044, 0.041]}
S10 = ["sham1", "sham2", "day3_rep1", "day3_rep2", "day7_rep1", "day7_rep2",
       "day21_rep1", "day21_rep2", "day60_rep1", "day60_rep2"]
sal_rows = rd(os.path.join(R, "salmon_gate1_baseline.tsv"))
for g, txs in (("Atp2a2", txs_atp), ("Agpat3", txs_agp)):
    got = []
    for s in S10:
        with open(os.path.join(BASE, "data", "salmon", s, "quant.sf"), encoding="utf-8") as f:
            hdr = f.readline().split("\t"); ti, pi = hdr.index("Name"), hdr.index("TPM")
            tot = lon = 0.0
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) < len(hdr): continue
                if p[ti] in txs:
                    v = float(p[pi]); tot += v
                    if p[ti] == ("ENSMUST00000179939.7" if g == "Atp2a2" else "ENSMUST00000001240.11"):
                        lon = v
        got.append(round(lon / tot, 3))
    check(f"salmon long_frac {g} recomputed", got == exp_sal[g], f"got {got}")

# ---------- 4. QKI ----------
di = {}
for line in open(os.path.join(R, "qki_distal_intersections.tsv"), encoding="utf-8"):
    p = line.rstrip("\n").split("\t")
    di.setdefault(p[3], []).append((p[6], int(p[7]), int(p[8]), float(p[10])))
a_hits = di.get(next(r["event_id"] for r in seg if r["gene_symbol"] == "Atp2a2"), [])
check("Atp2a2 distal peak 122454200-122454250 FDR 1.0416e-05",
      any(h[1] == 122454200 and h[2] == 122454250 and abs(h[3] - 1.04156605028874e-05) < 1e-12 for h in a_hits), str(a_hits))
n = next(r["event_id"] for r in seg if r["gene_symbol"] == "Ndrg2")
gg = next(r["event_id"] for r in seg if r["gene_symbol"] == "Agpat3")
check("Ndrg2/Agpat3 distal QKI = none", not di.get(n) and not di.get(gg))
summ = rd(os.path.join(R, "qki_per_event_summary.tsv"))
n05 = sum(1 for r in summ if r["L2v2_distal_qki"] == "yes")
check("events with distal QKI FDR<0.05 = 42", n05 == 42, str(n05))

# ---------- 5. PAS distances recomputed from atlas ----------
idx = {}
with gzip.open(os.path.join(BASE, "data", "references", "PAS", "atlas.clusters.2.0.GRCm38.96.bed.gz"), "rt", encoding="utf-8") as f:
    for line in f:
        p = line.split("\t")
        if len(p) < 6: continue
        c = p[0] if p[0].startswith("chr") else "chr" + p[0]
        idx.setdefault((c, p[5]), []).append(int(p[1]))
for k in idx: idx[k].sort()
def nearest(chrom, strand, pos):
    arr = idx.get((chrom, strand))
    if not arr: return None
    i = bisect.bisect_left(arr, pos)
    best = None
    for j in (i - 1, i):
        if 0 <= j < len(arr):
            d = abs(arr[j] - pos)
            best = d if best is None else min(best, d)
    return best
exp_pas = {"Atp2a2": 13, "Ndrg2": 8, "Agpat3": 171, "Sirt2": 159, "Aplp1": 182,
           "Cnp": 203, "Pea15a": 130, "Kazn": 125, "Plec": 231, "Fam107a": 235}
for g, expd in exp_pas.items():
    r = rank[g]
    got = nearest(r["chrom"], r["strand"], int(r["proximal_pas_bed"]))
    check(f"PAS distance {g} = {expd}", got == expd, str(got))

# ---------- 6. motif recount ----------
hits = 0; segs8 = set(); F8 = {"QKI","ELAVL1","PTBP1","MBNL2","NOVA1","NOVA2","TARDBP","HNRNPA2B1"}
with open(os.path.join(R, "fimo_distal", "fimo.tsv"), encoding="utf-8") as f:
    hdr = f.readline().lstrip("#").strip().split("\t")
    ci = {c: i for i, c in enumerate(hdr)}
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < len(hdr): continue
        hits += 1
        if float(p[ci["p-value"]]) < 1e-4 and p[ci["motif_id"]].split("_")[0].upper() in F8:
            segs8.add(p[ci["sequence_name"]])
check("fimo distal hits = 61488", hits == 61488, str(hits))
check("frozen8-hit segments = 792/824", len(segs8) == 792, str(len(segs8)))
enr = rd(os.path.join(R, "motif_family_enrichment.tsv"))
minp = min(float(r["fisher_p"]) for r in enr)
minq = min(float(r["fisher_q_bh"]) for r in enr)
check("enrichment: min p=0.0537, all q>=0.71, 16 families",
      len(enr) == 16 and abs(minp - 0.0537) < 5e-4 and minq >= 0.71, f"minp={minp} minq={minq}")

# ---------- 7. conservation ----------
cons = rd(os.path.join(R, "distal_conservation.tsv"))
ny = sum(1 for r in cons if r["L4v2_conserved"] == "yes")
check("L4v2 conserved 9667/9693 (99.7%)", ny == 9667 and len(cons) == 9693, str(ny))
ca = next(r for r in cons if r["gene_symbol"] == "Atp2a2")
check("Atp2a2 distal pos_bases recorded", int(ca["distal_pos_bases"]) > 0, ca["distal_pos_bases"])

# ---------- 8. class / diffs ----------
m2 = rd(os.path.join(R, "event_annotation_matrix.v2.tsv"))
cc = {"B": 0, "C": 0, "D": 0, "none": 0}
for r in m2: cc[r["class_v2"]] += 1
check("class v2 B=8 C=783 D=1 none=185", cc == {"B": 8, "C": 783, "D": 1, "none": 185}, str(cc))
old = {r["event"]: r for r in rd(os.path.join(IL, "p3_p4_event_annotation_matrix.tsv"))}
dcls = sum(1 for r in m2 if r["event_num"] in old and old[r["event_num"]]["class"] != r["class_v2"])
check("class flips = 346", dcls == 346, str(dcls))

# ---------- 9. percentiles (recompute background) ----------
gsm_meta = {}
with gzip.open(os.path.join(BASE, "data", "gse143531_series_matrix.txt.gz"), "rt", encoding="utf-8") as f:
    geos = parts = animals = None
    for line in f:
        if line.startswith("!Sample_geo_accession"): geos = [x.strip('"') for x in line.strip().split("\t")[1:]]
        elif line.startswith("!Sample_characteristics_ch1"):
            vals = [x.strip('"') for x in line.strip().split("\t")[1:]]
            if vals and vals[0].startswith("cell part:"): parts = [v.split(":",1)[1].strip() for v in vals]
            elif vals and vals[0].startswith("animal id:"): animals = [v.split(":",1)[1].strip().replace("Animal ","") for v in vals]
for gsm, part, an in zip(geos, parts, animals): gsm_meta[gsm] = (an, part)
D14 = os.path.join(BASE, "data", "gse143531")
mat = {}
for fn in sorted(os.listdir(D14)):
    if not fn.endswith(".tsv.gz"): continue
    gsm = fn.split("_")[0]
    with gzip.open(os.path.join(D14, fn), "rt", encoding="utf-8") as f:
        rdr = csv.reader(f, delimiter="\t"); next(rdr)
        for row in rdr:
            if len(row) >= 2: mat.setdefault(row[0], {})[gsm] = int(row[1])
groups = {}
for gsm, (an, part) in gsm_meta.items(): groups.setdefault((an, part), []).append(gsm)
bg = []
for gid, samples in mat.items():
    rs = []
    for an in ("C1", "C2", "C3"):
        fu = sum(samples.get(g, 0) for g in groups[(an, "Full Astrocytes")])
        pa = sum(samples.get(g, 0) for g in groups[(an, "PAP")])
        if fu >= 10 and pa >= 1: rs.append(math.log2(pa / fu))
    if len(rs) == 3: bg.append(sum(rs) / 3)
bg.sort()
sym2gid = {}
for r in ev: sym2gid.setdefault(r["gene_symbol"], r["event_id"].split("|")[1])
def pct(sym):
    gid = sym2gid[sym]
    rs = []
    for an in ("C1", "C2", "C3"):
        fu = sum(mat[gid].get(g, 0) for g in groups[(an, "Full Astrocytes")])
        pa = sum(mat[gid].get(g, 0) for g in groups[(an, "PAP")])
        if fu >= 10 and pa >= 1: rs.append(math.log2(pa / fu))
    m = sum(rs) / 3
    return round(100.0 * bisect.bisect_left(bg, m) / len(bg), 1)
check("percentiles Atp2a2 54.0 / Agpat3 51.8 / Aplp1 32.8 / Sirt2 15.0 / Plec 95.8",
      pct("Atp2a2") == 54.0 and pct("Agpat3") == 51.8 and pct("Aplp1") == 32.8
      and pct("Sirt2") == 15.0 and pct("Plec") == 95.8,
      f"{pct('Atp2a2')},{pct('Agpat3')},{pct('Aplp1')},{pct('Sirt2')},{pct('Plec')}")

# ---------- 10. review table + P0 hash spot ----------
rv = rd(os.path.join(BR, "event_review_table.tsv"))
adv = sum(1 for r in rv if r["decision"] == "advance")
hi = sum(1 for r in rv if r["confidence"] == "high")
hold = sum(1 for r in rv if r["decision"].startswith("hold"))
check("review: 8 advance (2 high) / 2 hold / 0 drop", adv == 8 and hi == 2 and hold == 2 and len(rv) == 10,
      f"adv={adv} hi={hi} hold={hold}")
inv = rd(os.path.join(BASE, "results", "00_inventory", "input_files_sha256.tsv"))
import hashlib
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()
spot = [r for r in inv if r["file"] in ("input_links/p2_full_pdui_matrix.tsv", "01_研究Proposal.md")]
okh = all(sha(os.path.join(BASE, r["file"])) == r["sha256"] for r in spot)
check("P0 sha256 spot-verify (2 files)", okh and len(inv) >= 21, f"{len(inv)} rows")
zp = r"D:\stroke_apa\ATP2A2_IGV_CHECK_PACKAGE_20260920.zip"
check("IGV zip size 29569120", os.path.getsize(zp) == 29569120 if os.path.exists(zp) else False)

print()
print("AUDIT PASS 1:", "ALL PASS" if not FAILS else f"{len(FAILS)} FAILURES -> {FAILS}")
sys.exit(1 if FAILS else 0)
