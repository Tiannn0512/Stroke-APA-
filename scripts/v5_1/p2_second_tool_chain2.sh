#!/bin/bash
# Overnight chain v2 (2026-09-19): QAPA second-tool assets, correct CLI.
# salmon lives in envs/salmon (full path); qapa in envs/qapa.
# Steps: gtfToGenePred -> qapa build -> qapa fasta -> salmon index -> 12x quant -> qapa quant.
# PolyASite retry with Referer (else pending). Log: results/p2_second_tool_chain2.log
set -uo pipefail
LOG=/mnt/d/stroke_apa/results/p2_second_tool_chain2.log
exec >> "$LOG" 2>&1
echo "[chain2] start $(date)"
QAPA=/home/taylor/miniconda3/envs/qapa/bin/qapa
SALMON=/home/taylor/miniconda3/envs/salmon/bin/salmon
R=/home/taylor/reference_v2
W=/home/taylor/stroke_APA_work
source /home/taylor/miniconda3/etc/profile.d/conda.sh
conda activate qapa 2>/dev/null || true
command -v gtfToGenePred >/dev/null || conda install -y -c bioconda ucsc-gtfToGenePred

cd "$R"
# --- PolyASite retry (non-fatal) ---
mkdir -p PAS && cd PAS
if [ ! -s mm10_atlas.clusters.bed ] && [ ! -s mm10_atlas.clusters.bed.gz ]; then
  curl -sL --max-time 300 -A "Mozilla/5.0 (X11; Linux x86_64)" \
    -e "https://polyasite.unibas.ch/" \
    -o mm10_atlas.clusters.bed.gz \
    "https://polyasite.unibas.ch/download/atlas/2.0/GRCm38.atlas.clusters.bed.gz" \
  && file mm10_atlas.clusters.bed.gz | grep -q gzip \
  && gunzip -f mm10_atlas.clusters.bed.gz && echo "[pas] OK $(wc -l < mm10_atlas.clusters.bed)" \
  || echo "[pas] FAILED again (pending)"
else
  echo "[pas] present"
fi
cd "$R"

# --- genePred + qapa build + fasta ---
if [ ! -s gencode.vM25.annotation.gtf.genePred ]; then
  gtfToGenePred gencode.vM25.annotation.gtf gencode.vM25.annotation.gtf.genePred && echo "[genepred] OK"
fi
if [ ! -s qapa_vm25.bed ]; then
  $QAPA build --db qapa_vm25 -p PAS/mm10_atlas.clusters.bed gencode.vM25.annotation.gtf.genePred > qapa_vm25.bed \
    && echo "[qapa_build] OK $(wc -l < qapa_vm25.bed) UTRs" \
    || { echo "[qapa_build] retry without PolyASite"; $QAPA build --db qapa_vm25 gencode.vM25.annotation.gtf.genePred > qapa_vm25.bed && echo "[qapa_build_nopas] OK"; }
fi
if [ ! -s qapa_vm25.fa ]; then
  $QAPA fasta -f genome.fa qapa_vm25.bed qapa_vm25.fa && echo "[qapa_fasta] OK $(grep -c '^>' qapa_vm25.fa)"
fi

# --- salmon index ---
if [ ! -d salmon_idx_vM25 ]; then
  $SALMON index -t gencode.vM25.transcripts.fa -i salmon_idx_vM25 -k 31 -p 8 && echo "[salmon_idx] OK"
fi

# --- 12x salmon quant (original FASTQs on /mnt/d) ---
mkdir -p "$W/salmon"
SAMPLES="sham1 sham2 day1_rep1 day1_rep2 day3_rep1 day3_rep2 day7_rep1 day7_rep2 day21_rep1 day21_rep2 day60_rep1 day60_rep2"
for s in $SAMPLES; do
  if [ -s "$W/salmon/$s/quant.sf" ]; then echo "[quant] $s exists"; continue; fi
  $SALMON quant -i "$R/salmon_idx_vM25" -l A -p 8 --validateMappings \
    -1 "/mnt/d/stroke_apa_data/fastq/${s}_1.fastq.gz" -2 "/mnt/d/stroke_apa_data/fastq/${s}_2.fastq.gz" \
    -o "$W/salmon/$s" > "$W/salmon/$s.log" 2>&1 \
    && echo "[quant] $s OK $(date '+%T')" || echo "[quant] $s FAIL (see log)"
done

# --- qapa quant merge ---
cd "$W/salmon"
$QAPA quant -F salmon --db "$R/qapa_vm25.fa" */quant.sf > qapa_usage.tsv 2> qapa_quant.err \
  && echo "[qapa_quant] OK $(wc -l < qapa_usage.tsv) rows" \
  || { echo "[qapa_quant] FAIL"; tail -5 qapa_quant.err; }
echo "[chain2] END $(date)"
