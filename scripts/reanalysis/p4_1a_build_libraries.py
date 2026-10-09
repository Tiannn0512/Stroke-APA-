#!/usr/bin/env python3
"""Task A step 1-2  : GSE143531 library metadata,
technical-rep merging, QC, and a 6-library count matrix.

Units (per task brief + Mazaré 2020 methods): each biological LIBRARY = 8 technical
sequencing files summed by gene ID. PAP libraries reportedly pool 2 mice each;
Full libraries 1 mouse. C1/C2/C3 GEO labels alone do NOT prove PAP/Full pairing ->
analysis unit = library, primary model unpaired (decided in R script).

Outputs (results/04_public_pap/):
  library_metadata.tsv / counts_by_biological_library.tsv / library_qc.tsv
  counts_per_file_raw.tsv
"""
import csv, gzip, os, re

BASE = r"D:\stroke_apa_reanalysis"
D = os.path.join(BASE, "data", "gse143531")
SMX = os.path.join(BASE, "data", "gse143531_series_matrix.txt.gz")
OUT = os.path.join(BASE, "results", "04_public_pap")
os.makedirs(OUT, exist_ok=True)

geos = titles = parts = animals = None
with gzip.open(SMX, "rt", encoding="utf-8") as f:
    for line in f:
        if line.startswith("!Sample_geo_accession"): geos = re.findall(r'"([^"]*)"', line)
        elif line.startswith("!Sample_title"): titles = re.findall(r'"([^"]*)"', line)
        elif line.startswith("!Sample_characteristics_ch1"):
            vals = re.findall(r'"([^"]*)"', line)
            if vals and vals[0].startswith("cell part:"): parts = [v.split(":",1)[1].strip() for v in vals]
            elif vals and vals[0].startswith("animal id:"): animals = [v.split(":",1)[1].strip().replace("Animal ","") for v in vals]

meta_rows = []
for gsm, title, part, animal in zip(geos, titles, parts, animals):
    m = re.match(r"(\d{4})_(\d+)([a-z]+):\s*(.*?)\s+technical replicate (\d+)$", title)
    lib = m.group(1) + m.group(2)   # e.g. 2017_545 -> 2017545（与文件名前缀一致）
    suffix = m.group(3)
    rep = int(m.group(5))
    meta_rows.append({
        "GSM": gsm, "file_prefix": f"{gsm}_expression_output_expression_{lib}{suffix}",
        "library_id": lib, "cell_part": part, "animal_label": animal, "technical_rep": rep,
        "pool_size": {"PAP": 2, "Full Astrocytes": 1}[part],
        "biological_unit": f"{lib}_{part.replace(' ','')}",
        "paired_status": "unpaired (label match not sufficient; pool composition differs by fraction)",
    })
with open(os.path.join(OUT, "library_metadata.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(meta_rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(meta_rows)

mat = {}; sizes = []
for fn in sorted(os.listdir(D)):
    if not fn.endswith(".tsv.gz"): continue
    gsm = fn.split("_")[0]
    tot = 0; ngenes = 0
    with gzip.open(os.path.join(D, fn), "rt", encoding="utf-8") as f:
        rdr = csv.reader(f, delimiter="\t"); next(rdr)
        for row in rdr:
            if len(row) < 2: continue
            c = int(row[1]); tot += c; ngenes += 1
            key = fn.replace(".tsv.gz", "")
            mat.setdefault(row[0], {})[key] = c
    sizes.append({"file": fn, "gsm": gsm, "total_counts": tot, "n_genes_rows": ngenes})
with open(os.path.join(OUT, "counts_per_file_raw.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["file","gsm","total_counts","n_genes_rows"], delimiter="\t")
    w.writeheader(); w.writerows(sizes)

prefix2meta = {r["file_prefix"]: r for r in meta_rows}
gsm2lib = {r["GSM"]: r["library_id"] for r in meta_rows}
libs = {}   # library_id -> gene -> summed counts（8 个技术文件按 gene ID 求和）
for g, per in mat.items():
    for key, c in per.items():
        lib = gsm2lib[key.split("_")[0]]
        libs.setdefault(lib, {})[g] = libs.get(lib, {}).get(g, 0) + c
lib_ids = sorted(libs)
with open(os.path.join(OUT, "counts_by_biological_library.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["gene_id"] + lib_ids)
    for g in sorted(mat):
        w.writerow([g] + [libs[k].get(g, 0) for k in lib_ids])

lib2meta = {r["library_id"]: r for r in meta_rows}
qc = []
for k in lib_ids:
    m = lib2meta[k]
    tot = sum(libs[k].values())
    nz = sum(1 for v in libs[k].values() if v >= 1)
    qc.append({"library_id": k, "biological_unit": m["biological_unit"], "cell_part": m["cell_part"],
               "animal_label": m["animal_label"], "n_tech_files": 8, "total_counts": tot,
               "genes_detected_ge1": nz, "paired_status": m["paired_status"]})
with open(os.path.join(OUT, "library_qc.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(qc[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(qc)

print("libraries:")
for r in qc:
    print(f"  {r['library_id']}  {r['cell_part']:<16} animal={r['animal_label']:<3} counts={r['total_counts']:>12,}  genes>=1: {r['genes_detected_ge1']:,}")
fulls = [s["total_counts"] for s in sizes if prefix2meta[s["file"].replace(".tsv.gz","")]["cell_part"] == "Full Astrocytes"]
paps  = [s["total_counts"] for s in sizes if prefix2meta[s["file"].replace(".tsv.gz","")]["cell_part"] == "PAP"]
print(f"per-file: Full {min(fulls):,}-{max(fulls):,}; PAP {min(paps):,}-{max(paps):,}")
