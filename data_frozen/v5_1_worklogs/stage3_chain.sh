#!/bin/bash
# Stage-3 chain (serialized after alignment to avoid RAM contention with STAR):
#   wait align pipeline end -> M4 reference sets build -> ensure pilot DaPars2 + PDUI QC
# Align end detected by: no STAR/fastp process AND last expected bedgraph (day1_rep2) exists.
set -uo pipefail
LOG=/d/stroke_apa_data/stage3_chain.log
PY="/d/Auto Claw/AutoClaw/resources/python/python.exe"
{
echo "[stage3] start $(date '+%F %T')"

# 1) wait for alignment pipeline to end (60-min stall guard)
while true; do
  BUSY=$(wsl.exe bash -c "ps aux | grep -E 'STAR-avx2|fastp|gse238125_align' | grep -v grep | wc -l" 2>/dev/null | tr -d '\r')
  BUSY=${BUSY:-1}
  if [ "$BUSY" = "0" ]; then break; fi
  sleep 120
done
echo "[stage3] alignment ended $(date '+%F %T')"
ls -la /home/taylor/stroke_APA_work/bedgraph/ 2>/dev/null || wsl.exe bash -c "ls -la /home/taylor/stroke_APA_work/bedgraph/"

# 2) M4 reference sets (retry after earlier OOM alongside STAR)
"$PY" /d/stroke_apa/scripts/p3_m4_reference_sets.py >> /d/stroke_apa/results/p3_m4_build.log 2>&1
echo "[stage3] m4 build exit=$? $(date '+%F %T')"

# 3) ensure pilot DaPars2 ran (stage2 wrapper may or may not be alive; avoid double-run)
if ls /mnt/d/stroke_apa_data/dapars2/dapars2_subset*PDUI* >/dev/null 2>&1 || \
   wsl.exe bash -c "ls /mnt/d/stroke_apa_data/dapars2/*PDUI* >/dev/null 2>&1"; then
  echo "[stage3] pilot DaPars2 output already present (stage2 wrapper did it)"
else
  echo "[stage3] running pilot DaPars2 subset $(date '+%F %T')"
  wsl.exe bash -c "bash /mnt/d/stroke_apa/scripts/p2_dapars2_run.sh sham1,sham2,day1_rep1,day1_rep2"
  echo "[stage3] dapars2 exit=$? $(date '+%F %T')"
fi

# 4) PDUI QC (descriptive)
"$PY" /d/stroke_apa/scripts/p2_pilot_qc.py
echo "[stage3] qc exit=$? end $(date '+%F %T')"
} >> "$LOG" 2>&1
