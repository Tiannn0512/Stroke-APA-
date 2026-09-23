#!/bin/bash
# Watcher: wait for pilot downloader to finish (single ENA pipe -> serialize downloads),
# then launch the 8-run downloader. 15-min stall abort for safety.
set -uo pipefail
LOG=/d/stroke_apa_data/dl_8runs.log
{
echo "[watcher] start $(date '+%F %T')"
last=0; stall=0
while ! grep -q "ALL PILOT FILES DONE" /d/stroke_apa_data/dl_pilot_2.log 2>/dev/null; do
  cur=$(du -sb /d/stroke_apa_data/fastq/*.parts 2>/dev/null | awk '{s+=$1} END{print s+0}')
  if [ "$cur" = "$last" ]; then stall=$((stall+1)); else stall=0; fi
  last=$cur
  if [ "$stall" -ge 15 ]; then
    echo "[watcher] ABORT: pilot stalled 15 min — launching 8runs anyway (resume-safe) $(date '+%F %T')"
    break
  fi
  sleep 60
done
echo "[watcher] pilot done or timed out -> launching 8runs $(date '+%F %T')"
"/d/Auto Claw/AutoClaw/resources/python/python.exe" /d/stroke_apa/scripts/p2_dl_8runs.py
echo "[watcher] 8runs exit=$? $(date '+%F %T')"
} >> "$LOG" 2>&1
