#!/bin/bash
# Watcher: when DaPars2 writes PDUI output -> run pilot QC immediately
LOG=/d/stroke_apa_data/dapars2_watch.log
PY="/d/Auto Claw/AutoClaw/resources/python/python.exe"
{
echo "[dwatch] start $(date '+%F %T')"
while ! ls /mnt/d/stroke_apa_data/dapars2/*_PDUI* >/dev/null 2>&1; do
  sleep 60
done
echo "[dwatch] PDUI files appeared $(date '+%F %T')"
"$PY" /d/stroke_apa/scripts/p2_pilot_qc.py
echo "[dwatch] qc exit=$? end $(date '+%F %T')"
} >> "$LOG" 2>&1
