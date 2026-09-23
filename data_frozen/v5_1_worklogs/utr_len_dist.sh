#!/bin/bash
awk -F'\t' '{print $3-$2}' /home/taylor/reference_v2/gencode_M25_3UTR_for_DaPars2.bed | sort -n | awk '
{a[NR]=$1}
END {
  print "n=" NR
  print "min=" a[1]
  print "p25=" a[int(NR*0.25)]
  print "median=" a[int(NR*0.5)]
  print "p75=" a[int(NR*0.75)]
  print "p90=" a[int(NR*0.9)]
  print "max=" a[NR]
  c100=0; c200=0; c500=0; c1000=0
  for(i=1;i<=NR;i++){ if(a[i]>=100)c100++; if(a[i]>=200)c200++; if(a[i]>=500)c500++; if(a[i]>=1000)c1000++ }
  print "ge100=" c100
  print "ge200=" c200
  print "ge500=" c500
  print "ge1000=" c1000
}'
