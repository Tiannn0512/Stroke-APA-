#!/usr/bin/env python3
"""P1.1/P1.2 input builder: dapars2_event_coordinates.tsv
Joins the PDUI matrix (event = DaPars2 name `tx|gene|symbol`) with strand from
GENCODE vM25 genePred (keyed by transcript id), and applies the P1.1 coordinate
determination: proximal_pas_raw = Predicted_Proximal_APA as output;
proximal_pas_bed = same value (no 1bp conversion; both on BED 0-based scale).
Output columns follow TODO P1.2: event_id, chrom, utr_start, utr_end,
proximal_pas_raw, proximal_pas_bed, strand, gene_symbol.
"""
import csv, os, sys

BASE = r"D:\stroke_apa_reanalysis"
PDUI = os.path.join(BASE, "input_links", "p2_full_pdui_matrix.tsv")
GENEPRED = os.path.join(BASE, "data", "references", "gencode.vM25.annotation.gtf.genePred")
OUT = os.path.join(BASE, "input_links", "dapars2_event_coordinates.tsv")

# 1) transcript -> strand
tx_strand = {}
with open(GENEPRED, encoding="utf-8") as f:
    for line in f:
        if not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        tx_strand[fields[0]] = fields[2]
print(f"genePred transcripts: {len(tx_strand)}")

# 2) events from PDUI matrix
rows, seen, dups, missing_strand, bad = [], {}, 0, [], []
event_num = 0  # 1-based data-row index; this is the integer `event` used by sap tables / old matrix
with open(PDUI, encoding="utf-8") as f:
    rdr = csv.DictReader(f, delimiter="\t")
    for r in rdr:
        event_num += 1
        name = r["Gene"]
        loci = r["Loci"]
        pas = r["Predicted_Proximal_APA"].strip()
        try:
            chrom, span = loci.split(":")
            us, ue = (int(x) for x in span.split("-"))
            pas_i = int(float(pas))
        except Exception as e:
            bad.append((name, loci, pas, str(e)))
            continue
        strand = tx_strand.get(name.split("|")[0])
        if strand is None:
            missing_strand.append(name)
            continue
        if name in seen:
            dups += 1
            name = f"{name}#dup{dups}"
        seen[name] = True
        rows.append({
            "event_id": name, "event_num": event_num, "chrom": chrom, "utr_start": us, "utr_end": ue,
            "proximal_pas_raw": pas_i, "proximal_pas_bed": pas_i,
            "strand": strand, "gene_symbol": r.get("symbol", name.split("|")[-1]),
        })

with open(OUT, "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(rows)

symbols = {r["gene_symbol"] for r in rows}
for g in ("Atp2a2", "Ndrg2", "Agpat3"):
    for r in rows:
        if r["gene_symbol"] == g:
            print(f"{g}: {r['chrom']}:{r['utr_start']}-{r['utr_end']} strand {r['strand']} pas {r['proximal_pas_raw']}")
print(f"events written: {len(rows)}; dups renamed: {dups}; missing strand: {len(missing_strand)}; bad parse: {len(bad)}")
if bad[:5]:
    print("bad examples:", bad[:5])
