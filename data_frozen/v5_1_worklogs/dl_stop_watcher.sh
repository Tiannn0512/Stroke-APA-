#!/bin/bash
# Watcher: once day3_rep2_2 is DONE, kill the 8runs downloader (user wants to try faster
# alternative download routes; resume is safe if restarted later).
LOG=/d/stroke_apa_data/dl_8runs.log
{
echo "[stopwatcher] start $(date '+%F %T')"
while ! grep -q "^\[DONE\] day3_rep2_2" "$LOG" 2>/dev/null; do
  sleep 30
done
echo "[stopwatcher] day3_rep2_2 done -> killing 8runs python $(date '+%F %T')"
PID=$(powershell -NoProfile -Command "Get-CimInstance Win32_Process -Filter \"Name='python.exe'\" | Where-Object {\$_.CommandLine -like '*p2_dl_8runs*'} | Select-Object -ExpandProperty ProcessId" 2>/dev/null | tr -d '\r')
if [ -n "$PID" ]; then
  taskkill //PID $PID //F
  echo "[stopwatcher] killed PID=$PID"
else
  echo "[stopwatcher] no 8runs python found (already exited?)"
fi
echo "[stopwatcher] end"
} >> /d/stroke_apa_data/dl_stop.log 2>&1
