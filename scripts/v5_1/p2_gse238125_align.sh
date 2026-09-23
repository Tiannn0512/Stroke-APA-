#!/bin/bash
# P2-05 GSE238125 alignment + coverage pipeline (WSL side)
# Semantics frozen 2026-09-16, reusing ~/stroke_APA_pilot assets:
#   STAR 2.7.11b-avx2, genome dir = pilot STAR_index, runThreadN 8, BAM Unsorted, GeneCounts
#   Coverage bedgraph = bedtools genomecov -bg -split  (CORRECTED semantics per reanalysis_v3 TODO-03)
#   Depth = primary mapped alignments (samtools view -c -F 0x904)
set -euo pipefail

SAMPLE_MAP="/mnt/d/stroke_apa/scripts/p2_sample_map.tsv"
FASTQ_DIR="/mnt/d/stroke_apa_data/fastq"
# 2026-09-17: all work dirs moved to ext4 (9p/NTFS was throttling fastp/STAR/genomecov;
# cf. sham1: mapping 27min + bedgraph ~60min via /mnt/d). Products archived to D: at the end.
WORK=/home/taylor/stroke_APA_work
BAM_DIR=/home/taylor/stroke_APA_work/bam
COV_DIR=/home/taylor/stroke_APA_work/bedgraph
QC_DIR=/home/taylor/stroke_APA_work/qc
STAR_IDX=/home/taylor/reference_v2/STAR_index

mkdir -p "$WORK" "$BAM_DIR" "$COV_DIR" "$QC_DIR"

# depth accumulator
DEPTH_TSV=/home/taylor/stroke_APA_work/sequencing_depth.tsv
echo -e "sample\tdepth" > "$DEPTH_TSV.new"

# STAR tmp must live on a Linux filesystem: FIFOs cannot be created on /mnt/d (NTFS drvfs)
STAR_TMP_BASE=/home/taylor/STARtmp
mkdir -p "$STAR_TMP_BASE"

# iterate map (skip header)
tail -n +2 "$SAMPLE_MAP" | while IFS=$'\t' read -r name srr; do
  r1="$FASTQ_DIR/${name}_1.fastq.gz"
  r2="$FASTQ_DIR/${name}_2.fastq.gz"
  if [ ! -s "$r1" ] || [ ! -s "$r2" ]; then
    echo "[SKIP] $name ($srr): fastq missing"
    continue
  fi
  # already aligned: only refresh the depth row from existing BAM (stage2 rerun safety)
  if [ -s "$COV_DIR/${name}.bedgraph" ] && [ -s "$BAM_DIR/${name}.sorted.bam" ]; then
    depth=$(samtools view -c -F 0x904 "$BAM_DIR/${name}.sorted.bam")
    echo -e "${name}\t${depth}" >> "$DEPTH_TSV.new"
    echo "[SKIP] $name ($srr): already aligned, depth recomputed"
    continue
  fi
  echo "[STEP] $name ($srr): fastp trim"
  fastp -i "$r1" -I "$r2" \
        -o "$WORK/${name}_1.trimmed.fastq.gz" -O "$WORK/${name}_2.trimmed.fastq.gz" \
        --html "$QC_DIR/${name}_fastp.html" --json "$QC_DIR/${name}_fastp.json" \
        --thread 8 --detect_adapter_for_pe 2> "$QC_DIR/${name}_fastp.log"

  echo "[STEP] $name: STAR align"
  mkdir -p "$BAM_DIR/${name}_"
  rm -rf "$BAM_DIR/${name}_/_STARtmp" "$STAR_TMP_BASE/${name}_STARtmp"
  STAR --runThreadN 8 --genomeDir "$STAR_IDX" \
       --readFilesIn "$WORK/${name}_1.trimmed.fastq.gz" "$WORK/${name}_2.trimmed.fastq.gz" \
       --readFilesCommand zcat \
       --outFileNamePrefix "$BAM_DIR/${name}_/" \
       --outTmpDir "$STAR_TMP_BASE/${name}_STARtmp" \
       --outSAMtype BAM Unsorted --quantMode GeneCounts

  echo "[STEP] $name: sort + index"
  samtools sort -@ 4 -o "$BAM_DIR/${name}.sorted.bam" "$BAM_DIR/${name}_/Aligned.out.bam"
  samtools index "$BAM_DIR/${name}.sorted.bam"
  rm -f "$BAM_DIR/${name}_/Aligned.out.bam"

  echo "[STEP] $name: bedgraph (-split) + depth"
  bedtools genomecov -ibam "$BAM_DIR/${name}.sorted.bam" -bg -split > "$COV_DIR/${name}.bedgraph"
  depth=$(samtools view -c -F 0x904 "$BAM_DIR/${name}.sorted.bam")
  echo -e "${name}\t${depth}" >> "$DEPTH_TSV.new"

  # cleanup trimmed intermediates
  rm -f "$WORK/${name}_1.trimmed.fastq.gz" "$WORK/${name}_2.trimmed.fastq.gz"
  echo "[DONE] $name"
done

sort -k1,1 "$DEPTH_TSV.new" > "$DEPTH_TSV"
rm -f "$DEPTH_TSV.new"
echo "ALL SAMPLES PROCESSED. Depth table:"
cat "$DEPTH_TSV"
