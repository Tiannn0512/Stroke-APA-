#!/bin/bash
# P2-06 DaPars2 runner for GSE238125 (12 samples: d1/d3/d7/d21/d60 each vs Sham, n=2/点)
# Reuses pilot annotation: gencode_M25_3UTR_for_DaPars2.bed (tandem 3'UTR annotation)
# Coverage semantics = -split bedgraph (matches depth tsv from align pipeline)
# Optional arg $1 = comma-separated sample subset (e.g. "sham1,sham2,day1_rep1,day1_rep2")
# for the pilot 4-run subset cfg; default = all 12.
set -euo pipefail

DAPARS2_SRC=/home/taylor/DaPars2/src
ANNO=/home/taylor/reference_v2/gencode_M25_3UTR_for_DaPars2.bed
COV_DIR=/home/taylor/stroke_APA_work/bedgraph
OUT_BASE=/mnt/d/stroke_apa_data/dapars2
DEPTH_TSV=/home/taylor/stroke_APA_work/sequencing_depth_dapars2.tsv
CHROM=/home/taylor/reference_v2/chromosomes.txt
PY=/home/taylor/miniconda3/envs/dapars2/bin/python

ALL="day1_rep1 day1_rep2 day21_rep1 day21_rep2 day3_rep1 day3_rep2 day60_rep1 day60_rep2 day7_rep1 day7_rep2 sham1 sham2"
SUBSET="${1:-}"
if [ -n "$SUBSET" ]; then
  SAMPLES=$(echo "$SUBSET" | tr ',' ' ')
  TAG="subset$(echo "$SUBSET" | tr ',' '-')"
else
  SAMPLES="$ALL"
  TAG="gse238125"
fi

mkdir -p "$OUT_BASE"

# build bedgraph list (alphabetical sample order as in map)
BG_LIST=""
for name in $SAMPLES; do
  [ -s "$COV_DIR/${name}.bedgraph" ] || { echo "missing bedgraph: $name"; exit 1; }
  BG_LIST="${BG_LIST}${COV_DIR}/${name}.bedgraph,"
done
BG_LIST=${BG_LIST%,}

CFG="$OUT_BASE/dapars2_${TAG}.cfg"
cat > "$CFG" <<EOF
Aligned_Wig_files=$BG_LIST
Output_directory=$OUT_BASE
Annotated_3UTR=$ANNO
Output_result_file=dapars2_${TAG}
sequencing_depth_file=$DEPTH_TSV
Num_Threads=8
Coverage_threshold=10
EOF

echo "[RUN] DaPars2_Multi_Sample_Multi_Chr.py ($TAG)"
cd "$OUT_BASE"
"$PY" "$DAPARS2_SRC/DaPars2_Multi_Sample_Multi_Chr.py" "$CFG" "$CHROM"

echo "[DONE] results in $OUT_BASE"
ls -la "$OUT_BASE" | head
