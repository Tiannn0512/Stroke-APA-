#!/bin/bash
# Stage-2 chain: night_chain end -> pilot ALL DONE -> md5 gate -> align 4 runs -> DaPars2 subset -> PDUI QC
# Sequential by design (only one STAR / one DaPars2 at a time; RAM 23GB constraint).
set -uo pipefail
LOG=/d/stroke_apa_data/stage2_chain.log
PY="/d/Auto Claw/AutoClaw/resources/python/python.exe"
{
echo "[stage2] start $(date '+%F %T')"

# 1) wait for night chain to end — only for an "end" line appended AFTER this run starts
# (old log already contains "[chain] end" from previous failed runs)
NL0=$(wc -l < /d/stroke_apa_data/night_chain.log 2>/dev/null || echo 0)
while ! tail -n +$((NL0 + 1)) /d/stroke_apa_data/night_chain.log 2>/dev/null | grep -q "\[chain\] end"; do
  sleep 120
done
echo "[stage2] night chain ended $(date '+%F %T')"

# 2) wait for pilot ALL DONE (15-min stall detector)
last=0; stall=0
while ! grep -q "ALL PILOT FILES DONE" /d/stroke_apa_data/dl_pilot_2.log 2>/dev/null; do
  cur=$(du -sb /d/stroke_apa_data/fastq/*.parts 2>/dev/null | awk '{s+=$1} END{print s+0}')
  if [ "$cur" = "$last" ]; then stall=$((stall+1)); else stall=0; fi
  last=$cur
  if [ "$stall" -ge 15 ]; then echo "[stage2] ABORT: pilot stalled 15 min at $cur bytes"; exit 2; fi
  sleep 60
done
echo "[stage2] pilot download complete $(date '+%F %T')"

# 3) md5 gate on the 8 pilot files
RES=$("$PY" /d/stroke_apa/scripts/p2_md5_check.py 2>&1 | grep -E "sham|day1")
echo "$RES"
if echo "$RES" | grep -q "MISMATCH"; then echo "[stage2] GATE FAIL: md5 mismatch"; exit 3; fi
NOK=$(echo "$RES" | grep -c "\[OK\]")
echo "[stage2] md5 OK=$NOK (need 8)"
if [ "$NOK" -lt 8 ]; then echo "[stage2] GATE FAIL: files missing"; exit 3; fi
echo "[stage2] GATE PASS -> align $(date '+%F %T')"

# 4) align all complete pilot samples (regenerates consistent 4-row depth table)
# bash -c wrap: prevents MSYS path conversion mangling /mnt/d/... args passed to wsl.exe
wsl.exe bash -c "bash /mnt/d/stroke_apa/scripts/p2_gse238125_align.sh"
echo "[stage2] align exit=$? $(date '+%F %T')"
echo "--- depth ---"; cat /d/stroke_apa_data/sequencing_depth.tsv

# 5) DaPars2 pilot subset (4 runs)
wsl.exe bash -c "bash /mnt/d/stroke_apa/scripts/p2_dapars2_run.sh sham1,sham2,day1_rep1,day1_rep2"
echo "[stage2] dapars2 exit=$? $(date '+%F %T')"
ls -la /mnt/d/stroke_apa_data/dapars2/ 2>/dev/null | head -20

# 6) PDUI QC (descriptive)
"$PY" /d/stroke_apa/scripts/p2_pilot_qc.py
echo "[stage2] qc exit=$? end $(date '+%F %T')"
} >> "$LOG" 2>&1
