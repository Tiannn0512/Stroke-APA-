#!/bin/bash
for d in /home/taylor/stroke_APA_work/dapars2_full/chr*/; do
  c=$(basename "$d")
  e=$(grep -c "exit=0" "$d/run.log" 2>/dev/null)
  a=$(grep "assigned events:" "$d/run.log" 2>/dev/null | grep -o "[0-9]*")
  echo "$c done=$e assigned=${a:-0}"
done | sort -V
echo "--- running workers ---"
ps aux | grep "[D]aPars2_Multi" | grep "R+" | sed "s/.*dapars2_full\///;s/\/chr.txt.*//" | sort | uniq -c
free -g | sed -n 2p
