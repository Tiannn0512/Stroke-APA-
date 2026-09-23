#!/bin/bash
# L4 conservation assets: bigWig + bigWigToBedGraph -> elements = score>0 regions (official phastCons elements)
set -uo pipefail
LOG=/mnt/d/stroke_apa/results/p3_l4_conservation.log
exec >> "$LOG" 2>&1
echo "[l4] start $(date)"
D=/home/taylor/reference_v2/phastcons60way
mkdir -p "$D" && cd "$D"
EXE=https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64
[ -s bigWigToBedGraph ] || { curl -sLO "$EXE/bigWigToBedGraph" && chmod +x bigWigToBedGraph; }
[ -s bigWigAverageOverBed ] || { curl -sLO "$EXE/bigWigAverageOverBed" && chmod +x bigWigAverageOverBed; }
ls -la bigWigToBedGraph bigWigAverageOverBed
# bigWig downloaded via Windows proxy to /mnt/d (2026-09-20; WSL direct was 0.15MB/s)
[ -s /mnt/d/stroke_apa_data/mm10.60way.phastCons.bw ] || { echo "[l4] bigWig missing on /mnt/d"; exit 1; }
BW=/mnt/d/stroke_apa_data/mm10.60way.phastCons.bw
ls -la "$BW"
./bigWigToBedGraph "$BW" all_scores.bedgraph
# elements = score>0 (phastCons is 0 outside predicted elements)
awk 'BEGIN{OFS="\t"}$4>0{print $1,$2,$3}' all_scores.bedgraph > phastCons60way_elements.bed
echo "[l4] elements: $(wc -l < phastCons60way_elements.bed) regions"
# lost segment overlap >=10bp
bedtools intersect -a /mnt/d/stroke_apa/results/fimo_seqs/lost_segments.bed \
  -b phastCons60way_elements.bed -wao \
  | awk 'BEGIN{OFS="\t"}{key=$4; ov=$NF; if(NF>10) for(i=12;i<=NF-3;i++){} print}' > /dev/null 2>&1 || true
bedtools intersect -a /mnt/d/stroke_apa/results/fimo_seqs/lost_segments.bed \
  -b phastCons60way_elements.bed -wao \
  | awk 'BEGIN{OFS="\t"}{k=$4; sum[k]+=$NF; if($NF>=10) hit[k]=1} END{for(k in sum) print k, (k in hit)?"CONSERVED":"not", sum[k]}' \
  > lost_segments_conservation.tsv
echo "[l4] scored: $(wc -l < lost_segments_conservation.tsv) segments; conserved: $(awk '$2==\"CONSERVED\"' lost_segments_conservation.tsv | wc -l)"
echo "[l4] END $(date)"
