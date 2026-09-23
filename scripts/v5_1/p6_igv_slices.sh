#!/bin/bash
# P6 IGV read-evidence package: slice candidate-gene loci from 6 BAMs
# (sham x2, day3 x2, day7 x2) for Atp2a2 / Ndrg2 / Agpat3, whole gene +-5kb.
set -euo pipefail
GTF=/home/taylor/reference_v2/gencode.vM25.annotation.gtf
W=/home/taylor/stroke_APA_work
OUT=$W/igv_slices
GENES="Atp2a2 Ndrg2 Agpat3"
SAMPLES="sham1 sham2 day3_rep1 day3_rep2 day7_rep1 day7_rep2"
FLANK=5000

mkdir -p "$OUT"
printf 'gene\tchr\tstart\tend\tregion_start\tregion_end\n' > "$OUT/gene_spans.tsv"
: > "$OUT/regions.bed"

for g in $GENES; do
  line=$(awk -F'\t' -v G="gene_name \"$g\"" '$3=="gene" && index($9,G)>0 {print; exit}' "$GTF")
  [ -n "$line" ] || { echo "MISSING gene $g" >&2; exit 1; }
  chr=$(printf '%s' "$line" | cut -f1)
  s=$(printf '%s' "$line" | cut -f4)
  e=$(printf '%s' "$line" | cut -f5)
  rs=$(( s - FLANK )); [ "$rs" -lt 1 ] && rs=1
  re=$(( e + FLANK ))
  printf '%s\t%s\t%s\t%s\n' "$chr" "$rs" "$re" "$g" >> "$OUT/regions.bed"
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$g" "$chr" "$s" "$e" "$rs" "$re" >> "$OUT/gene_spans.tsv"
done

echo "=== regions ==="
cat "$OUT/regions.bed"

for s in $SAMPLES; do
  echo "[slice] $s"
  samtools view -b -L "$OUT/regions.bed" "$W/bam/${s}.sorted.bam" > "$OUT/${s}.cand3.bam"
  samtools index "$OUT/${s}.cand3.bam"
  # per-sample align QC (STAR) + fastp log if present
  cp "$W/bam/${s}_/Log.final.out" "$OUT/qc_${s}_STAR_Log.final.out" 2>/dev/null || true
  cp "$W/bam/${s}_/fastp.log" "$OUT/qc_${s}_fastp.log" 2>/dev/null || true
  cp "$W/bam/${s}_/fastp.json" "$OUT/qc_${s}_fastp.json" 2>/dev/null || true
done

# DaPars2 full-run config actually used
cp /mnt/d/stroke_apa_data/dapars2/dapars2_gse238125.cfg "$OUT/dapars2_gse238125.cfg" 2>/dev/null || true

echo "=== done ==="
ls -la "$OUT"
