#!/usr/bin/env python3
# P6 package: extract candidate-focus rows (Atp2a2/Ndrg2/Agpat3) from the three main tables.
import csv, os, sys

BASE = r"D:\stroke_apa"
OUT = os.path.join(BASE, "igv_package", "01_computed_results", "candidate_focus")
os.makedirs(OUT, exist_ok=True)
SYMBOLS = {"Atp2a2", "Ndrg2", "Agpat3"}

def read_tsv(path):
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def write_tsv(path, rows, fieldnames):
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

# 1) annotation matrix rows for the three candidates
matrix = read_tsv(os.path.join(BASE, "results", "p3_p4_event_annotation_matrix.tsv"))
hit = [r for r in matrix if r["symbol"] in SYMBOLS]
write_tsv(os.path.join(OUT, "candidates_event_annotation.tsv"), hit, list(hit[0].keys()))
event_ids = {r["event"] for r in hit}
print(f"matrix rows: {len(hit)}; events: {sorted(event_ids)}")

# 2) SAP per-contrast stats for those events
sap = read_tsv(os.path.join(BASE, "results", "p2_sap_layer2_events.tsv"))
sap_hit = [r for r in sap if r["event"] in event_ids]
write_tsv(os.path.join(OUT, "candidates_sap_stats.tsv"), sap_hit, list(sap[0].keys()))
print(f"sap rows: {len(sap_hit)} (contrasts x events)")

# 3) per-sample PDUI for the three genes (all 12 samples)
pdui_path = os.path.join(BASE, "results", "p2_full_pdui_matrix.tsv")
with open(pdui_path, encoding="utf-8") as f:
    rdr = csv.reader(f, delimiter="\t")
    header = next(rdr)
    keep = [header] + [row for row in rdr if row and row[0].split("|")[-1] in SYMBOLS]
with open(os.path.join(OUT, "candidates_pdui_per_sample.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerows(keep)
print(f"pdui rows: {len(keep)-1}")
