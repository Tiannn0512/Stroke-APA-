#!/bin/bash
# Serial chain: wait align_batch2 -> realign day7_rep2 (now complete) -> slim DaPars2 (2-way)
set -uo pipefail
LOG=/d/stroke_apa_data/morning_serial.log
{
echo "[serial] start $(date '+%F %T')"

# 1) wait for the running alignment batch to end (process-based)
while true; do
  BUSY=$(ps aux | grep -E 'STAR-avx2|fastp|gse238125_align' | grep -v grep | wc -l)
  if [ "$BUSY" = "0" ]; then break; fi
  sleep 120
done
echo "[serial] batch1 alignment ended $(date '+%F %T')"

# 2) realign (only day7_rep2 will actually run; others SKIP)
wsl.exe bash -c "bash /mnt/d/stroke_apa/scripts/p2_gse238125_align.sh" >> /d/stroke_apa_data/align_batch3_day7rep2.log 2>&1
echo "[serial] day7_rep2 align exit=$? $(date '+%F %T')"
wsl.exe bash -c "ls /home/taylor/stroke_APA_work/bedgraph/ | wc -l; cat /home/taylor/stroke_APA_work/sequencing_depth.tsv"

# 3) refresh depth tables (align run rewrote .new with all 12 rows)
wsl.exe bash -c "tail -n +2 /home/taylor/stroke_APA_work/sequencing_depth.tsv.new | grep -v '^sample' | sort -k1,1 > /home/taylor/stroke_APA_work/sequencing_depth.tsv; cp /home/taylor/stroke_APA_work/sequencing_depth.tsv /home/taylor/stroke_APA_work/sequencing_depth_dapars2.tsv; wc -l /home/taylor/stroke_APA_work/sequencing_depth_dapars2.tsv"
echo "[serial] depth tables refreshed $(date '+%F %T')"

# 4) slim DaPars2 (2-way per-chromosome) on 12 samples? NO - pilot first: 4 samples per pre-registration
#    (pilot = sham1,sham2,day1_rep1,day1_rep2; full 12-sample run is a separate step after QC)
sed 's/^SAMPLES=.*/SAMPLES="sham1,sham2,day1_rep1,day1_rep2"/' /mnt/d/stroke_apa_data/dapars2_slim_run.sh > /mnt/d/stroke_apa_data/dapars2_slim_run_serial.sh
sed -i 's|^LOG=/mnt/d/stroke_apa_data/dapars2_slim_run.log|LOG=/mnt/d/stroke_apa_data/dapars2_slim_run_serial.log|' /mnt/d/stroke_apa_data/dapars2_slim_run_serial.sh
# strip the alignment-wait block (already serialized here): replace wait loop with no-op
python3 - <<'PYEOF'
import re
p = "/mnt/d/stroke_apa_data/dapars2_slim_run_serial.sh"
t = open(p).read()
t = t.replace("""# 1) wait for alignment to finish (no STAR/fastp/align processes)
while true; do
  BUSY=$(ps aux | grep -E 'STAR-avx2|fastp|gse238125_align' | grep -v grep | wc -l)
  if [ "$BUSY" = "0" ]; then break; fi
  sleep 120
done""", "# wait block removed (serialized by caller)")
open(p, "w").write(t)
print("slim serial script ready")
PYEOF
bash /mnt/d/stroke_apa_data/dapars2_slim_run_serial.sh
echo "[serial] slim dapars2 exit=$? $(date '+%F %T')"

# 5) PDUI QC
"/d/Auto Claw/AutoClaw/resources/python/python.exe" /d/stroke_apa/scripts/p2_pilot_qc.py >> /d/stroke_apa/results/p2_pilot_qc.log 2>&1
echo "[serial] qc done; end $(date '+%F %T')"
} >> "$LOG" 2>&1
