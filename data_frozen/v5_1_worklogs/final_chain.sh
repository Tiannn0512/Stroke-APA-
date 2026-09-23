#!/bin/bash
# Final pilot chain (post-rebuild): fix depth table -> day3_rep1 realign -> DaPars2 subset -> PDUI QC
set -uo pipefail
LOG=/d/stroke_apa_data/final_chain.log
PY="/d/Auto Claw/AutoClaw/resources/python/python.exe"
{
echo "[final] start $(date '+%F %T')"

# 1) regenerate depth table from .new (sort was skipped when align chain died)
wsl.exe bash -c "sort -k1,1 /home/taylor/stroke_APA_work/sequencing_depth.tsv.new > /home/taylor/stroke_APA_work/sequencing_depth.tsv && cat /home/taylor/stroke_APA_work/sequencing_depth.tsv"
echo "[final] depth table restored $(date '+%F %T')"

# 2) align day3_rep1 (skips completed 4 samples)
wsl.exe bash -c "bash /mnt/d/stroke_apa/scripts/p2_gse238125_align.sh"
echo "[final] align exit=$? $(date '+%F %T')"
wsl.exe bash -c "ls -la /home/taylor/stroke_APA_work/bedgraph/; cat /home/taylor/stroke_APA_work/sequencing_depth.tsv"

# 3) DaPars2 pilot subset (4 runs)
wsl.exe bash -c "bash /mnt/d/stroke_apa/scripts/p2_dapars2_run.sh sham1,sham2,day1_rep1,day1_rep2"
echo "[final] dapars2 exit=$? $(date '+%F %T')"
wsl.exe bash -c "ls -la /mnt/d/stroke_apa_data/dapars2/ | head -12"

# 4) PDUI QC
"$PY" /d/stroke_apa/scripts/p2_pilot_qc.py
echo "[final] qc exit=$? end $(date '+%F %T')"
} >> "$LOG" 2>&1
