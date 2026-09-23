#!/bin/bash
# Night chain: wait sham1_2 redo -> md5 gate -> alignment closed loop (P2-05c/d)
# Any gate failure stops the chain without touching downstream.
set -uo pipefail
LOG=/d/stroke_apa_data/night_chain.log
OUT=/d/stroke_apa_data/fastq/sham1_2.fastq.gz
PARTS=/d/stroke_apa_data/fastq/sham1_2.fastq.gz.redo.parts
M1_EXP=34ced9d056b611e1e117cd456125e92c
M2_EXP=85d69dbe99f36c75b73475ef24699787

{
echo "[chain] start $(date '+%F %T')"

# 1) wait for redo completion (file assembled) with 10-min stall detector
last=0; stall=0
while [ ! -f "$OUT" ] || [ -d "$PARTS" ]; do
  cur=$(du -sb "$PARTS" 2>/dev/null | cut -f1); cur=${cur:-0}
  if [ "$cur" = "$last" ]; then stall=$((stall+1)); else stall=0; fi
  last=$cur
  if [ "$stall" -ge 10 ]; then
    echo "[chain] ABORT: redo stalled 10 min at $cur bytes $(date '+%F %T')"
    exit 2
  fi
  sleep 60
done
echo "[chain] redo finished $(date '+%F %T')"

# 2) md5 gate
M1=$(md5sum /d/stroke_apa_data/fastq/sham1_1.fastq.gz | cut -d' ' -f1)
M2=$(md5sum "$OUT" | cut -d' ' -f1)
echo "[chain] sham1_1 md5=$M1 (want $M1_EXP)"
echo "[chain] sham1_2 md5=$M2 (want $M2_EXP)"
if [ "$M2" != "$M2_EXP" ] || [ "$M1" != "$M1_EXP" ]; then
  echo "[chain] GATE FAIL: md5 mismatch — NOT aligning"
  exit 3
fi
echo "[chain] GATE PASS -> alignment $(date '+%F %T')"

# 3) closed loop: fastp -> STAR -> bedgraph -> depth (processes any complete sample)
# bash -c wrap: prevents MSYS path conversion mangling /mnt/d/... args passed to wsl.exe
wsl.exe bash -c "bash /mnt/d/stroke_apa/scripts/p2_gse238125_align.sh"
echo "[chain] align exit=$? $(date '+%F %T')"
echo "--- bedgraph ---"; ls -la /d/stroke_apa_data/bedgraph/ 2>/dev/null
echo "--- depth ---"; cat /d/stroke_apa_data/sequencing_depth.tsv 2>/dev/null
echo "[chain] end $(date '+%F %T')"
} >> "$LOG" 2>&1
