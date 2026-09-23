#!/bin/bash
# DaPars2 single-process run on slim gene-level bed (4-sample pilot).
# Avoids the multiprocessing Manager dict instability seen in per-chr chunked runs.
set -uo pipefail
LOG=/mnt/d/stroke_apa_data/dapars2_single_run.log
R=/home/taylor/reference_v2
W=/home/taylor/stroke_APA_work
OUT=/mnt/d/stroke_apa_data/dapars2_single
DP=/home/taylor/DaPars2/src
PY=/home/taylor/miniconda3/envs/dapars2/bin/python
{
echo "[single] start $(date '+%F %T')"
mkdir -p "$OUT"
cd "$OUT"
printf 'Aligned_Wig_files=%s/bedgraph/sham1.bedgraph,%s/bedgraph/sham2.bedgraph,%s/bedgraph/day1_rep1.bedgraph,%s/bedgraph/day1_rep2.bedgraph\nOutput_directory=%s\nAnnotated_3UTR=%s/gencode_M25_3UTR_slim.bed\nOutput_result_file=dapars2_slim4\nsequencing_depth_file=%s/sequencing_depth_dapars2.tsv\nNum_Threads=1\nCoverage_threshold=10\n' \
  "$W" "$W" "$W" "$W" "$OUT" "$R" "$W" > job.cfg
"$PY" "$DP/DaPars2_Multi_Sample_Multi_Chr.py" job.cfg /home/taylor/reference_v2/chromosomes.txt
echo "[single] exit=$? $(date '+%F %T')"
ls -la "$OUT" | head -10
echo "[single] END $(date '+%F %T')"
} >> "$LOG" 2>&1
