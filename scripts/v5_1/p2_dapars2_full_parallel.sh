#!/bin/bash
# P2-06 full 12-sample DaPars2, per-chromosome independent instance.
# Usage (WSL): bash p2_dapars2_full_parallel.sh chr2 chr22 ...
# Each chr gets its own output dir + cfg (audit §16.1; parallelism = concurrent invocations).
# Known traps fixed: headerless depth tsv / relative Output_result_file / ext4 output /
# Num_Threads=1 / Coverage_threshold=10 (frozen, do not touch) / slim bed 3-segment names.
set -uo pipefail
R=/home/taylor/reference_v2
W=/home/taylor/stroke_APA_work
DP=/home/taylor/DaPars2/src
PY=/home/taylor/miniconda3/envs/dapars2/bin/python
BASE=$W/dapars2_full
mkdir -p "$BASE"

SAMPLES="day1_rep1 day1_rep2 day21_rep1 day21_rep2 day3_rep1 day3_rep2 day60_rep1 day60_rep2 day7_rep1 day7_rep2 sham1 sham2"
BG=""
for s in $SAMPLES; do
  [ -s "$W/bedgraph/$s.bedgraph" ] || { echo "[FATAL] missing bedgraph: $s"; exit 1; }
  BG="$BG$W/bedgraph/$s.bedgraph,"
done
BG=${BG%,}

for CHR in "$@"; do
  OUT="$BASE/$CHR"
  mkdir -p "$OUT"
  LOG="$OUT/run.log"
  {
    echo "[full:$CHR] start $(date '+%F %T')"
    echo "$CHR" > "$OUT/chr.txt"
    printf 'Aligned_Wig_files=%s\nOutput_directory=%s\nAnnotated_3UTR=%s/gencode_M25_3UTR_slim.bed\nOutput_result_file=dapars2_full_%s\nsequencing_depth_file=%s/sequencing_depth_dapars2.tsv\nNum_Threads=1\nCoverage_threshold=10\n' \
      "$BG" "$OUT" "$R" "$CHR" "$W" > "$OUT/job.cfg"
    ( cd "$OUT" && "$PY" "$DP/DaPars2_Multi_Sample_Multi_Chr.py" job.cfg "$OUT/chr.txt" )
    rc=$?
    echo "[full:$CHR] exit=$rc $(date '+%F %T')"
    ls -la "$OUT" | head -8
  } >> "$LOG" 2>&1
done
