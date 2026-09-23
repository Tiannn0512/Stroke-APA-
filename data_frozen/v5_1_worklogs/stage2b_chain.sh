#!/bin/bash
# Stage-2b: restart alignment after accidental VM kill (13:40), serialized after m4 build.
#   wait m4 build python exit -> align (auto-skips sham1) -> pilot DaPars2 -> PDUI QC
set -uo pipefail
LOG=/d/stroke_apa_data/stage2b_chain.log
PY="/d/Auto Claw/AutoClaw/resources/python/python.exe"
{
echo "[stage2b] start $(date '+%F %T')"

# 1) wait for stage3's m4 build (p3_m4_reference_sets.py) to finish - avoid RAM contention
while wsl.exe bash -c "ps aux | grep 'p3_m4_reference_sets' | grep -v grep >/dev/null 2>&1"; do
  sleep 90
done
echo "[stage2b] m4 build gone from process table -> align $(date '+%F %T')"

# 2) alignment: skips sham1 (already done), aligns sham2/day1_rep1/day1_rep2/day3_rep1
wsl.exe bash -c "bash /mnt/d/stroke_apa/scripts/p2_gse238125_align.sh"
echo "[stage2b] align exit=$? $(date '+%F %T')"
wsl.exe bash -c "ls -la /home/taylor/stroke_APA_work/bedgraph/; cat /home/taylor/stroke_APA_work/sequencing_depth.tsv"

# 3) pilot DaPars2 (4 runs)
wsl.exe bash -c "bash /mnt/d/stroke_apa/scripts/p2_dapars2_run.sh sham1,sham2,day1_rep1,day1_rep2"
echo "[stage2b] dapars2 exit=$? $(date '+%F %T')"

# 4) PDUI QC
"$PY" /d/stroke_apa/scripts/p2_pilot_qc.py
echo "[stage2b] qc exit=$? end $(date '+%F %T')"
} >> "$LOG" 2>&1
