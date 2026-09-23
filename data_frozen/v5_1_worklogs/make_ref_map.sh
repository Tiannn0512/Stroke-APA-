#!/bin/bash
# Build Ensembl(gene_id, version-stripped) -> symbol map from GENCODE vM25 gene rows
GTF=/home/taylor/stroke_APA_pilot/reference/genome/gencode.vM25.annotation.gtf
OUT=/mnt/d/stroke_apa_data/ref_map_vM25.tsv
awk -F'\t' '$3=="gene"' "$GTF" | sed -E 's/.*gene_id "([^"]+)".*gene_name "([^"]+)".*/\1\t\2/' | sed -E 's/\.([0-9]+)\t/\t/' > "$OUT"
wc -l "$OUT"
head -3 "$OUT"
