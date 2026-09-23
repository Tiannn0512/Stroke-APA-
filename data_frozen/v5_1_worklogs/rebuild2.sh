#!/bin/bash
# Manual rebuild v2 (linear, absolute paths, no cd):
#   unpack -> STAR index -> DaPars2 3'UTR bed -> chromosomes.txt
set -uo pipefail
LOG=/mnt/d/stroke_apa_data/rebuild2.log
R=/home/taylor/reference_v2
D=/mnt/d/stroke_apa_data/reference_rebuild
PYD=/home/taylor/miniconda3/envs/dapars2/bin/python
{
echo "[rebuild2] start $(date '+%F %T')"
cd "$R" || { echo "[rebuild2] FATAL: no $R"; exit 2; }

echo "[rebuild2] unpacking..."
gunzip -c "$D/gencode.vM25.annotation.gtf.gz" > "$R/gencode.vM25.annotation.gtf"
gunzip -c "$D/GRCm38.primary_assembly.genome.fa.gz" > "$R/genome.fa"
echo "[rebuild2] unpacked: gtf=$(stat -c %s "$R/gencode.vM25.annotation.gtf") fa=$(stat -c %s "$R/genome.fa") $(date '+%F %T')"

echo "[rebuild2] STAR genomeGenerate (sparseD2, ~40-90 min expected)..."
STAR --runMode genomeGenerate --runThreadN 10 \
     --genomeDir "$R/STAR_index" \
     --genomeFastaFiles "$R/genome.fa" \
     --sjdbGTFfile "$R/gencode.vM25.annotation.gtf" \
     --sjdbOverhang 99 \
     --genomeSAsparseD 2 \
     --limitGenomeGenerateRAM 20000000000 \
     --outFileNamePrefix "$R/STARidx_log_"
echo "[rebuild2] STAR index exit=$? $(date '+%F %T')"

# DaPars2 tandem 3'UTR bed: try official generator, else hand-build from GTF
"$PYD" /home/taylor/DaPars2/src/Generate_Annotation.py \
      "$R/gencode.vM25.annotation.gtf" "$R/gencode_M25_3UTR_for_DaPars2.bed" mm10 \
      > "$R/gen_anno.log" 2>&1 \
  && echo "[rebuild2] DaPars2 bed generated" \
  || echo "[rebuild2] Generate_Annotation failed (see gen_anno.log) - hand-build fallback needed"
ls -la "$R/gencode_M25_3UTR_for_DaPars2.bed" 2>/dev/null || true

printf "chr1\nchr2\nchr3\nchr4\nchr5\nchr6\nchr7\nchr8\nchr9\nchr10\nchr11\nchr12\nchr13\nchr14\nchr15\nchr16\nchr17\nchr18\nchr19\nchrX\nchrY\n" > "$R/chromosomes.txt"
echo "[rebuild2] chromosomes.txt done"
echo "[rebuild2] END $(date '+%F %T')"
} >> "$LOG" 2>&1
