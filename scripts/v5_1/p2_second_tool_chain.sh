#!/bin/bash
# Overnight chain (2026-09-19 night): second-tool consistency assets for Gate G2(d).
# 1) PolyASite mm10 atlas download (independent PAS annotation for top-event spot check)
# 2) salmon transcriptome index (gffread -> transcripts.fa)
# 3) salmon quant x12 on original FASTQs (D:), serialized
# 4) QAPA annotation + qapa quant merge
# Log: /mnt/d/stroke_apa/results/p2_second_tool_chain.log  (also nohup-safe)
set -uo pipefail
LOG=/mnt/d/stroke_apa/results/p2_second_tool_chain.log
exec >> "$LOG" 2>&1
echo "[chain] start $(date)"
source /home/taylor/miniconda3/etc/profile.d/conda.sh

# --- 1. PolyASite ---
mkdir -p /home/taylor/reference_v2/PAS
cd /home/taylor/reference_v2/PAS
if [ ! -s mm10_atlas.clusters.bed ]; then
  curl -sL -A "Mozilla/5.0" --max-time 300 -o mm10_atlas.clusters.bed.gz \
    "https://polyasite.unibas.ch/download/atlas/2.0/GRCm38.atlas.clusters.bed.gz" \
  && gunzip -f mm10_atlas.clusters.bed.gz \
  && echo "[pas] OK $(wc -l < mm10_atlas.clusters.bed) clusters" \
  || echo "[pas] DOWNLOAD FAILED (register as pending)"
else
  echo "[pas] exists"
fi

# --- 2. gffread + transcripts.fa + salmon index ---
conda activate qapa
command -v gffread >/dev/null || conda install -y -c bioconda gffread
cd /home/taylor/reference_v2
if [ ! -s gencode.vM25.transcripts.fa ]; then
  gffread -w gencode.vM25.transcripts.fa -g genome.fa gencode.vM25.annotation.gtf \
    && echo "[txfa] OK $(grep -c '^>' gencode.vM25.transcripts.fa) transcripts"
fi
if [ ! -d salmon_idx_vM25 ]; then
  salmon index -t gencode.vM25.transcripts.fa -i salmon_idx_vM25 -k 31 -p 8 \
    && echo "[salmon_idx] OK"
fi

# --- 3. qapa annotation ---
cd /home/taylor/reference_v2
if [ ! -s qapa_annotation_all.bed ]; then
  qapa annotation --gtf gencode.vM25.annotation.gtf -o qapa_annotation
  if ls qapa_annotation_all.bed >/dev/null 2>&1; then echo "[qapa_anno] OK"; else ls qapa_annotation*; fi
fi

# --- 4. salmon quant x12 ---
W=/home/taylor/stroke_APA_work
mkdir -p "$W/salmon"
SAMPLES="sham1 sham2 day1_rep1 day1_rep2 day3_rep1 day3_rep2 day7_rep1 day7_rep2 day21_rep1 day21_rep2 day60_rep1 day60_rep2"
for s in $SAMPLES; do
  if [ -s "$W/salmon/$s/quant.sf" ]; then echo "[quant] $s exists"; continue; fi
  salmon quant -i /home/taylor/reference_v2/salmon_idx_vM25 -l A -p 8 --validateMappings \
    -1 "/mnt/d/stroke_apa_data/fastq/${s}_1.fastq.gz" -2 "/mnt/d/stroke_apa_data/fastq/${s}_2.fastq.gz" \
    -o "$W/salmon/$s" && echo "[quant] $s OK $(date '+%T')"
done

# --- 5. qapa quant merge ---
cd "$W/salmon"
qapa quant --salmon -p 8 /home/taylor/reference_v2/qapa_annotation_all.bed */quant.sf \
  && echo "[qapa_quant] OK"
echo "[chain] END $(date)"
