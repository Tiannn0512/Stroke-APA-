#!/bin/bash
# DaPars2 slim + per-chromosome parallel runner (user-approved 2026-09-18):
#   split slim bed by chromosome -> 2-way parallel DaPars2 jobs -> merge results
# Memory-safe: 2 concurrent python jobs (~1-2GB each), runs AFTER alignment finishes.
set -uo pipefail
LOG=/mnt/d/stroke_apa_data/dapars2_slim_run.log
R=/home/taylor/reference_v2
W=/home/taylor/stroke_APA_work
OUT=/mnt/d/stroke_apa_data/dapars2
DP=/home/taylor/DaPars2/src
PY=/home/taylor/miniconda3/envs/dapars2/bin/python
SAMPLES="sham1,sham2,day1_rep1,day1_rep2"
TAG=slim4

{
echo "[slimrun] start $(date '+%F %T')"

# 1) wait for alignment to finish (no STAR/fastp/align processes)
while true; do
  BUSY=$(ps aux | grep -E 'STAR-avx2|fastp|gse238125_align' | grep -v grep | wc -l)
  if [ "$BUSY" = "0" ]; then break; fi
  sleep 120
done
echo "[slimrun] alignment done -> DaPars2 slim starts $(date '+%F %T')"

# 2) split slim bed by chromosome
mkdir -p "$R/slim_split" "$OUT/chunks"
awk -F'\t' '{print > "'"$R"'/slim_split/"$1".bed"}' "$R/gencode_M25_3UTR_slim.bed"
CHRS=$(ls "$R/slim_split" | sed 's/.bed//' | tr '\n' ' ')
echo "[slimrun] chromosomes with UTRs: $CHRS"

# 3) per-chromosome job
run_chr () {
  local chr=$1
  local od="$OUT/chunks/$chr"
  mkdir -p "$od"
  cd "$od"
  echo -e "$chr" > "$od/chromosomes_one.txt"
  printf 'Aligned_Wig_files=%s/%s.bedgraph,%s/%s.bedgraph,%s/%s.bedgraph,%s/%s.bedgraph\nOutput_directory=%s\nAnnotated_3UTR=%s/slim_split/%s.bed\nOutput_result_file=dapars2_%s_%s\nsequencing_depth_file=%s/sequencing_depth_dapars2.tsv\nNum_Threads=2\nCoverage_threshold=10\n' \
    "$W" sham1 "$W" sham2 "$W" day1_rep1 "$W" day1_rep2 \
    "$od" "$R" "$chr" "$TAG" "$chr" "$W" > "$od/job.cfg"
  "$PY" "$DP/DaPars2_Multi_Sample_Multi_Chr.py" "$od/job.cfg" "$od/chromosomes_one.txt" \
    > "$od/run.log" 2>&1
  echo "[slimrun] $chr exit=$? $(date '+%T')"
}
export -f run_chr
export R W OUT DP PY TAG

# 4) 2-way parallel over chromosomes
printf '%s\n' $CHRS | xargs -P 2 -I{} bash -c 'run_chr {}'
echo "[slimrun] all chromosome jobs done $(date '+%F %T')"

# 5) merge PDUI tables
MERGED="$OUT/dapars2_${TAG}_merged_3UTR_and_PDUI.txt"
head -1 "$OUT/chunks/chr1/dapars2_${TAG}_chr1_3UTR_and_PDUI.txt" > "$MERGED" 2>/dev/null
for f in "$OUT"/chunks/*/dapars2_${TAG}_*_3UTR_and_PDUI.txt; do
  tail -n +2 "$f" >> "$MERGED"
done
echo "[slimrun] merged -> $MERGED ($(wc -l < "$MERGED") lines)"
echo "[slimrun] END $(date '+%F %T')"
} >> "$LOG" 2>&1
