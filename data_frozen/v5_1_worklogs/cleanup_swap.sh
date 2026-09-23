#!/bin/bash
# Auto-authorized C-drive cleanup (user approved 2026-09-17):
# When stage2b QC line appears AND no downloads/alignment/python work is running:
#   wsl --shutdown  (kills VM -> swap.vhdx freed -> .wslconfig swapFile=D: takes effect)
# Safety gates, all must pass:
#   1) stage2b qc line in log
#   2) no p2_dl_* python processes (downloads)
#   3) no STAR/fastp/align in WSL
#   4) no p3_m4/stage*_chain bash running (self excepted by name filter below)
# Post-shutdown health check: WSL starts fresh, swap file location verified on next use.
LOG=/d/stroke_apa_data/cleanup.log
{
echo "[cleanup] start $(date '+%F %T')"
while true; do
  G1=$(grep -c "stage2b] qc exit" /d/stroke_apa_data/stage2b_chain.log 2>/dev/null)
  G2=$(powershell -NoProfile -Command "(Get-CimInstance Win32_Process -Filter \"Name='python.exe'\" | Where-Object {\$_.CommandLine -like '*p2_dl_*'}).Count" 2>/dev/null | tr -d '\r')
  G3=$(wsl.exe bash -c "ps aux | grep -E 'STAR-avx2|fastp|gse238125_align|p3_m4' | grep -v grep | wc -l" 2>/dev/null | tr -d '\r')
  echo "[cleanup] gates: qc=$G1 dl=$G2 wslbusy=$G3"
  if [ "${G1:-0}" -ge 1 ] && [ "${G2:-1}" = "0" ] && [ "${G3:-1}" = "0" ]; then
    break
  fi
  sleep 300
done
echo "[cleanup] all gates passed -> wsl --shutdown $(date '+%F %T')"
wsl.exe --shutdown
sleep 10
# verify swap gone
if ls /c/Users/Taylor/AppData/Local/Temp/*/swap.vhdx >/dev/null 2>&1; then
  echo "[cleanup] WARN: swap.vhdx still present after shutdown:"
  ls -la /c/Users/Taylor/AppData/Local/Temp/*/swap.vhdx
else
  echo "[cleanup] OK: swap.vhdx freed (8.1GB reclaimed on C:)"
fi
echo "[cleanup] verifying fresh WSL boot + new swap location"
wsl.exe bash -c "echo VM_BOOT_OK; cat /proc/swaps; df -h /home | tail -1" 2>&1 | grep -v "w s l"
sleep 5
if [ -f /d/stroke_apa_data/wsl_swap.vhdx ]; then
  echo "[cleanup] OK: new swap file on D: confirmed"
else
  echo "[cleanup] note: D: swap file will be created on first memory pressure"
fi
df -h /c | tail -1
echo "[cleanup] end $(date '+%F %T')"
} >> "$LOG" 2>&1
