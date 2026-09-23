#!/bin/bash
# Morning serial v2 (all-WSL paths): wait batch1 -> realign day7_rep2 -> depth tables -> slim DaPars2 -> QC
set -uo pipefail
LOG=/mnt/d/stroke_apa_data/morning_serial.log
PYW=/mnt/d/stroke_APA_winpy/python.exe
{
echo "[serial2] start $(date '+%F %T')"

# 1) wait for running alignment batch to end
while true; do
  BUSY=$(ps aux | grep -E 'STAR-avx2|fastp|gse238125_align' | grep -v grep | wc -l)
  if [ "$BUSY" = "0" ]; then break; fi
  sleep 120
done
echo "[serial2] batch1 ended $(date '+%F %T')"

# 2) realign day7_rep2 (others SKIP)
wsl_ex=0
bash /mnt/d/stroke_apa/scripts/p2_gse238125_align.sh >> /mnt/d/stroke_apa_data/align_batch3_day7rep2.log 2>&1 || wsl_ex=$?
echo "[serial2] day7_rep2 align exit=$wsl_ex $(date '+%F %T')"
ls /home/taylor/stroke_APA_work/bedgraph/ | wc -l

# 3) depth tables: rebuild .new -> final (guard against header duplicates)
awk -F'\t' 'NR>1 && $1!="sample" {print $1"\t"$2}' /home/taylor/stroke_APA_work/sequencing_depth.tsv.new | sort -k1,1 > /home/taylor/stroke_APA_work/sequencing_depth_dapars2.tsv
{ echo -e "sample\tdepth"; cat /home/taylor/stroke_APA_work/sequencing_depth_dapars2.tsv; } > /home/taylor/stroke_APA_work/sequencing_depth.tsv
wc -l /home/taylor/stroke_APA_work/sequencing_depth_dapars2.tsv
echo "[serial2] depth refreshed $(date '+%F %T')"

# 4) slim DaPars2 (4-sample pilot, 2-way per-chr), script already has /mnt/d LOG
bash /mnt/d/stroke_apa_data/dapars2_slim_run_serial.sh
echo "[serial2] slim exit=$? $(date '+%F %T')"

# 5) PDUI QC (WSL qapa python, has pydeseq2 not needed here; just pandas)
exec /home/taylor/miniconda3/envs/qapa/bin/python /mnt/d/stroke_apa/scripts/p2_pilot_qc_wsl.py
echo "[serial2] qc exit=$? end $(date '+%F %T')"
} >> "$LOG" 2>&1
