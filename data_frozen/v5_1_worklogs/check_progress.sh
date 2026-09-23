#!/bin/bash
for d in /mnt/d/stroke_apa_data/dapars2_single_chr*; do
  n=$(cat "$d"/dapars2_slim4_result_temp.chr*.txt 2>/dev/null | wc -l)
  echo "$(basename $d): $n rows"
done
echo "proc: $(pgrep -c -f DaPars2)"
