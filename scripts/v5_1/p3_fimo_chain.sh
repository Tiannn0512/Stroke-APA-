#!/bin/bash
# P3-03 FIMO chain (approved 2026-09-20: 1e-3 + BH, audit calibration protocol).
# Step A: sequence extraction (proximal PAS +/-250 of union-977 events; lost distal segments
#         of shortening-significant events) via bedtools getfasta from genome.fa.
# Step B: positive-control calibration on independent QKI CLIP peak sequences (smoke set).
# Step C: FIMO scan with data/attract/regulators.motifs.meme (325 motifs / 16 regulators).
# Log: results/p3_fimo_chain.log ; outputs: results/fimo_*
set -uo pipefail
LOG=/mnt/d/stroke_apa/results/p3_fimo_chain.log
exec >> "$LOG" 2>&1
echo "[fimo] start $(date)"
source /home/taylor/miniconda3/etc/profile.d/conda.sh
conda activate meme 2>/dev/null || { echo "[fimo] meme env activate FAILED"; exit 1; }
R=/home/taylor/reference_v2
W=/home/taylor/stroke_APA_work
OUTD=/mnt/d/stroke_apa/results
WORK=$W/fimo_seq
mkdir -p "$WORK"
command -v bedtools >/dev/null || conda install -y -c bioconda bedtools

# ---- Step A: BEDs from event tables (python side prepped) ----
# /mnt/d/stroke_apa/results/fimo_seqs/{pas_events.bed,lost_segments.bed} must exist (built by
# scripts/p3_fimo_makebeds.py); if missing, build here:
if [ ! -s "$OUTD/fimo_seqs/pas_events.bed" ]; then
  echo "[fimo] building BEDs via p3_fimo_makebeds.py"
  /home/taylor/miniconda3/envs/dapars2/bin/python /mnt/d/stroke_apa/scripts/p3_fimo_makebeds.py || exit 1
fi
cd "$WORK"
bedtools getfasta -fi "$R/genome.fa" -bed "$OUTD/fimo_seqs/pas_events.bed" -name -fo pas_events.raw.fa
bedtools getfasta -fi "$R/genome.fa" -bed "$OUTD/fimo_seqs/lost_segments.bed" -name -fo lost_segments.raw.fa
# strip getfasta's '::chr:start-end' suffix -- FIMO parses it as UCSC coords and mangles names
sed 's/::/_at_/' pas_events.raw.fa > pas_events.fa
sed 's/::/_at_/' lost_segments.raw.fa > lost_segments.fa
rm -f pas_events.raw.fa lost_segments.raw.fa
echo "[fimo] seqs: pas=$(grep -c '^>' pas_events.fa) lost=$(grep -c '^>' lost_segments.fa)"
ls -la pas_events.fa lost_segments.fa | awk '{print $5, $NF}'

# ---- Step B: positive-control calibration (freeze before candidate scan) ----
CAL=$OUTD/fimo_calibration
mkdir -p "$CAL" && cd "$CAL"
fimo --thresh 1e-3 --norc --oc cal_control \
  /mnt/d/stroke_apa/data/attract/regulators.motifs.meme \
  /mnt/d/stroke_apa/data/attract/smoke_peaks.fa
HITS=$(awk 'NR>1' cal_control/fimo.tsv 2>/dev/null | wc -l)
NSEQ=$(grep -c '^>' /mnt/d/stroke_apa/data/attract/smoke_peaks.fa)
echo "[fimo] CALIBRATION (frozen): $HITS motif-hit lines across $NSEQ positive-control QKI peak sequences"
if [ "$HITS" -lt 10 ]; then echo "[fimo] calibration too weak - check before candidates"; fi

# ---- Step C: candidate scans ----
cd "$WORK"
mkdir -p fimo_pas fimo_lost
fimo --thresh 1e-3 --norc --max-stored-scores 1000000 --oc fimo_pas \
  /mnt/d/stroke_apa/data/attract/regulators.motifs.meme pas_events.fa \
  && echo "[fimo] pas scan OK"
fimo --thresh 1e-3 --norc --max-stored-scores 5000000 --oc fimo_lost \
  /mnt/d/stroke_apa/data/attract/regulators.motifs.meme lost_segments.fa \
  && echo "[fimo] lost scan OK"
echo "[fimo] END $(date)"
