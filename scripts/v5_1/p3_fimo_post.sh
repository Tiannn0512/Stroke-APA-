#!/bin/bash
# P3-03 post-processing: null BED -> fimo on null -> family merge + BH finalize.
set -uo pipefail
LOG=/mnt/d/stroke_apa/results/p3_fimo_post.log
exec >> "$LOG" 2>&1
echo "[post] start $(date)"
source /home/taylor/miniconda3/etc/profile.d/conda.sh
conda activate meme 2>/dev/null || conda activate bioinfo 2>/dev/null
W=/home/taylor/stroke_APA_work
R=/home/taylor/reference_v2
OUTD=/mnt/d/stroke_apa/results
cd "$W/fimo_seq"

# 0) chrom sizes (from samtools faidx index)
if [ ! -s "$R/chrom_sizes.txt" ]; then
  awk 'BEGIN{OFS="\t"}{print $1, $2}' "$R/genome.fa.fai" > "$R/chrom_sizes.txt"
fi

# 1) null BED (length + chromosome-matched, seed 42)
if [ ! -s "$OUTD/fimo_seqs/null_segments.bed" ]; then
  /home/taylor/miniconda3/envs/dapars2/bin/python /mnt/d/stroke_apa/scripts/p3_fimo_nullbed.py
fi

# 2) fimo on null sequences
if [ ! -s fimo_null/fimo.tsv ]; then
  mkdir -p fimo_null
  bedtools getfasta -fi "$R/genome.fa" -bed "$OUTD/fimo_seqs/null_segments.bed" -name -fo null_segments.fa
  fimo --thresh 1e-3 --norc --max-stored-scores 5000000 --oc fimo_null \
    /mnt/d/stroke_apa/data/attract/regulators.motifs.meme null_segments.fa \
    && echo "[post] null scan OK"
fi

# 3) finalize: family merge + BH + lost/gained/retained + null enrichment
/home/taylor/miniconda3/envs/dapars2/bin/python /mnt/d/stroke_apa/scripts/p3_fimo_finalize.py
echo "[post] END $(date)"
