#!/bin/bash
# Stage-4: after stage2b finishes (its qc line), retry M4 build (RAM now free), then stop.
set -uo pipefail
LOG=/d/stroke_apa_data/stage4_chain.log
PY="/d/Auto Claw/AutoClaw/resources/python/python.exe"
{
echo "[stage4] start $(date '+%F %T')"
while ! grep -q "stage2b] qc exit" /d/stroke_apa_data/stage2b_chain.log 2>/dev/null; do
  sleep 180
done
echo "[stage4] stage2b done -> retry m4 build $(date '+%F %T')"
"$PY" /d/stroke_apa/scripts/p3_m4_reference_sets.py >> /d/stroke_apa/results/p3_m4_build.log 2>&1
echo "[stage4] m4 build exit=$? $(date '+%F %T')"
} >> "$LOG" 2>&1
