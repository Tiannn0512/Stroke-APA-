#!/bin/bash
# null UTR scan (v2) -- script file to avoid quoting issues
source /home/taylor/miniconda3/etc/profile.d/conda.sh
conda activate meme || { echo "ACTIVATE_FAIL"; exit 1; }
cd /home/taylor/stroke_APA_work/fimo_seq
mkdir -p fimo_null_utr
bedtools getfasta -fi /home/taylor/reference_v2/genome.fa \
  -bed /mnt/d/stroke_apa/results/fimo_seqs/null_utr.bed -name -fo null_utr.raw.fa \
  && sed 's/::/_at_/' null_utr.raw.fa > null_utr.fa \
  && fimo --thresh 1e-3 --norc --max-stored-scores 5000000 --oc fimo_null_utr \
     /mnt/d/stroke_apa/data/attract/regulators.motifs.meme null_utr.fa \
  && echo "NULL_SCAN_OK"
