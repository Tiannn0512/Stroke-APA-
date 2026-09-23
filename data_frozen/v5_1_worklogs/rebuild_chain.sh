#!/bin/bash
# Reference rebuild chain (after accidental ~/stroke_APA_pilot deletion 2026-09-17):
#   wait for GENCODE download -> build new STAR index in ext4 (RAM-adapted)
#   -> generate DaPars2 3'UTR bed + chromosomes.txt -> rerun day3_rep1 alignment
# New layout: /home/taylor/reference_v2 (NOT under any "pilot" name)
set -uo pipefail
LOG=/d/stroke_apa_data/rebuild.log
{
echo "[rebuild] start $(date '+%F %T')"
cd /home/taylor/reference_v2

# 1) wait for downloads (gtf + fa present and >1GB)
while true; do
  G=$(stat -c %s /mnt/d/stroke_apa_data/reference_rebuild/gencode.vM25.annotation.gtf.gz 2>/dev/null || echo 0)
  F=$(stat -c %s /mnt/d/stroke_apa_data/reference_rebuild/GRCm38.primary_assembly.genome.fa.gz 2>/dev/null || echo 0)
  echo "[rebuild] sizes: gtf=$G fa=$F $(date '+%T')"
  if [ "$G" -gt 10000000 ] && [ "$F" -gt 1000000000 ]; then break; fi
  sleep 60
done
echo "[rebuild] downloads complete $(date '+%F %T')"

# 2) unpack
gunzip -c /mnt/d/stroke_apa_data/reference_rebuild/gencode.vM25.annotation.gtf.gz > gencode.vM25.annotation.gtf
gunzip -c /mnt/d/stroke_apa_data/reference_rebuild/GRCm38.primary_assembly.genome.fa.gz > genome.fa
echo "[rebuild] unpacked $(date '+%F %T')"

# 3) STAR index (24GB RAM constraint -> genomeSAsparseD 2, per pilot lesson)
mkdir -p STAR_index
/home/taylor/miniconda3/envs/bioinfo/bin/python - <<'PYEOF' 2>/dev/null || true
PYEOF
STAR --runMode genomeGenerate --runThreadN 10 \
     --genomeDir /home/taylor/reference_v2/STAR_index \
     --genomeFastaFiles /home/taylor/reference_v2/genome.fa \
     --sjdbGTFfile /home/taylor/reference_v2/gencode.vM25.annotation.gtf \
     --sjdbOverhang 99 \
     --genomeSAsparseD 2 \
     --limitGenomeGenerateRAM 20000000000
echo "[rebuild] STAR index exit=$? $(date '+%F %T')"

# 4) DaPars2 annotation: tandem 3'UTR bed from GTF (same semantics as pilot gencode_M25_3UTR_for_DaPars2.bed)
/home/taylor/miniconda3/envs/dapars2/bin/python /home/taylor/DaPars2/src/Generate_Annotation.py \
  gencode.vM25.annotation.gtf gencode_M25_3UTR_for_DaPars2.bed mm10 2>/dev/null \
  || echo "[rebuild] Generate_Annotation.py failed - will hand-build tandem 3UTR bed"
ls -la gencode_M25_3UTR_for_DaPars2.bed 2>/dev/null || true

# 5) chromosomes.txt (mm10 primary)
printf "chr1\nchr2\nchr3\nchr4\nchr5\nchr6\nchr7\nchr8\nchr9\nchr10\nchr11\nchr12\nchr13\nchr14\nchr15\nchr16\nchr17\nchr18\nchr19\nchrX\nchrY\n" > chromosomes.txt
echo "[rebuild] chromosomes.txt written"

echo "[rebuild] rebuild done $(date '+%F %T')"
} >> "$LOG" 2>&1
