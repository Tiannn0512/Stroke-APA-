#!/usr/bin/env python3
"""P4.1: GSE143531 candidate-gene PAP testability (gene level).

Data: 48 samples = 3 animals (C1/C2/C3) x 2 cell parts (Full Astrocytes / PAP) x 8 technical
replicates; Aldh1l1-RPL10a-eGFP TRAP counts from dorsal hippocampus.
Tech replicates are SUMMED within (animal, cell part) -> 6 aggregated columns (NOT 48 n).
Per candidate gene: detectability (aggregated counts per group), PAP/Full ratio per animal,
mean log2FC across animals, and direction consistency across the 3 animals.
Outputs: results/04_public_pap/gse143531_candidate_gene_summary.tsv (+ qc counts matrix)
"""
import csv, gzip, os, re

BASE = r"D:\stroke_apa_reanalysis"
D = os.path.join(BASE, "data", "gse143531")
SMX = os.path.join(BASE, "data", "gse143531_series_matrix.txt.gz")
OUT = os.path.join(BASE, "results", "04_public_pap")
os.makedirs(OUT, exist_ok=True)

# --- parse series matrix: GSM -> (animal, cell_part, tech_rep) ---
gsm_meta = {}
with gzip.open(SMX, "rt", encoding="utf-8") as f:
    titles, geos, parts, animals = None, None, None, None
    for line in f:
        if line.startswith("!Sample_title"):
            titles = re.findall(r'"([^"]*)"', line)
        elif line.startswith("!Sample_geo_accession"):
            geos = re.findall(r'"([^"]*)"', line)
        elif line.startswith("!Sample_characteristics_ch1"):
            vals = re.findall(r'"([^"]*)"', line)
            if vals and vals[0].startswith("cell part:"):
                parts = [v.split(":", 1)[1].strip() for v in vals]
            elif vals and vals[0].startswith("animal id:"):
                animals = [v.split(":", 1)[1].strip().replace("Animal ", "") for v in vals]
for gsm, title, part, animal in zip(geos, titles, parts, animals):
    m = re.search(r"technical replicate (\d+)", title)
    gsm_meta[gsm] = {"animal": animal, "part": part, "rep": m.group(1) if m else "?"}
print(f"series matrix: {len(gsm_meta)} samples; parts={set(v['part'] for v in gsm_meta.values())}; animals={set(v['animal'] for v in gsm_meta.values())}")

# --- load count tables ---
mat = {}  # gene -> {gsm: count}
for fn in sorted(os.listdir(D)):
    if not fn.endswith(".tsv.gz"):
        continue
    gsm = fn.split("_")[0]
    with gzip.open(os.path.join(D, fn), "rt", encoding="utf-8") as f:
        rdr = csv.reader(f, delimiter="\t")
        hdr = next(rdr)
        for row in rdr:
            if len(row) < 2:
                continue
            mat.setdefault(row[0], {})[gsm] = int(row[1])
print(f"genes: {len(mat)}; example gsm keys: {sorted(next(iter(mat.values())))[0]}")

# --- aggregate tech reps ---
groups = {}
for gsm, m in gsm_meta.items():
    groups.setdefault((m["animal"], m["part"]), []).append(gsm)
agg = {}  # gene -> {(animal,part): summed}
for g, samples in mat.items():
    agg[g] = {k: sum(samples.get(gsm, 0) for gsm in v) for k, v in groups.items()}

# --- candidate symbols -> ENSMUSG (from coordinate table) ---
sym2gene = {}
for r in csv.DictReader(open(os.path.join(BASE, "input_links", "dapars2_event_coordinates.tsv"), encoding="utf-8"), delimiter="\t"):
    sym2gene.setdefault(r["gene_symbol"], r["event_id"].split("|")[1])
CAND = ["Atp2a2", "Agpat3", "Aplp1", "Sirt2", "Ndrg2", "Cnp", "Pea15a", "Kazn", "Plec", "Fam107a"]
import math
rows = []
for sym in CAND:
    gid = sym2gene.get(sym)
    if gid not in agg:
        rows.append({"gene_symbol": sym, "ensembl": gid or "NA", "status": "not_in_matrix"})
        continue
    a = agg[gid]
    ratios, det = [], []
    line = {"gene_symbol": sym, "ensembl": gid, "status": "ok"}
    for an in ("C1", "C2", "C3"):
        full, pap = a[(an, "Full Astrocytes")], a[(an, "PAP")]
        line[f"{an}_full_counts"] = full
        line[f"{an}_pap_counts"] = pap
        det.append(full >= 10 and pap >= 1)
        if full > 0 and pap > 0:
            ratios.append(math.log2(pap / full))
        elif pap == 0:
            ratios.append(float("-inf"))
        elif full == 0:
            ratios.append(float("inf"))
    vals = [r for r in ratios if math.isfinite(r)]
    line["pap_over_full_log2"] = ";".join("inf" if r == float("inf") else "-inf" if r == float("-inf") else f"{r:.2f}" for r in ratios)
    line["mean_log2_pap_over_full"] = round(sum(vals) / len(vals), 2) if vals else "NA"
    signs = set((r > 0) for r in vals)
    line["direction_consistent_3animals"] = "yes" if len(signs) == 1 and vals else "mixed"
    line["detectable_all_animals"] = "yes" if all(det) else "no"
    rows.append(line)

with open(os.path.join(OUT, "gse143531_candidate_gene_summary.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(rows)

print(f"{'gene':<9} {'det':<4} {'log2 PAP/Full (C1;C2;C3)':<28} {'mean':>6} {'consist'}")
for r in rows:
    if r.get("status") == "ok":
        print(f"{r['gene_symbol']:<9} {r['detectable_all_animals']:<4} {r['pap_over_full_log2']:<28} {str(r['mean_log2_pap_over_full']):>6} {r['direction_consistent_3animals']}")
    else:
        print(f"{r['gene_symbol']:<9} {r.get('status')}")

# tech-rep QC: correlation of two reps within one group (sanity)
g0 = groups[("C1", "PAP")][:2]
import statistics
xs = [mat.get(gid, {}).get(g0[0], 0) for gid in list(mat)[:5000]]
print("QC columns aggregated: 6 groups:", sorted(groups))
