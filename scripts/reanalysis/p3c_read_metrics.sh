#!/bin/bash
# Task C: per-sample read-level metrics for the Atp2a2 review region (UTR + 1kb flanks)
set -uo pipefail
REG="chr5:122452512-122458153"
OUT=results/03_bam_review
echo -e "sample\treads\tmapq_median\tmapq_lt20_pct\tsoftclip_base_pct\tspliced_read_pct\tnm_median" > $OUT/Atp2a2_read_metrics.tsv
for s in sham1 sham2 day3_rep1 day3_rep2 day7_rep1 day7_rep2; do
  samtools view data/bam/${s}.sorted.bam $REG | awk -v s="$s" '
  {
    n++
    mapq=$5; mq[mq_n++] = mapq; if (mapq<20) low++
    cigar=$6
    sc=0; body=0; spliced=0
    tmp=cigar
    while (match(tmp, /[0-9]+[MIDNSHP=X]/)) {
      op=substr(tmp, RSTART+RLENGTH-1, 1); len=substr(tmp, RSTART, RLENGTH-1)+0
      if (op=="S"||op=="H") sc+=len
      if (op=="M"||op="="||op=="X") body+=len
      if (op=="N") spliced=1
      tmp=substr(tmp, RSTART+RLENGTH)
    }
    softbases+=sc; totalbases+=body
    if (spliced) sp++
    for (i=12; i<=NF; i++) if ($i ~ /^NM:i:/) { nm=substr($i,6); nms[nm_n++]=nm; break }
  }
  END {
    asort(mq); med = mq[int(mq_n/2)]
    asort(nms); nmed = nms[int(nm_n/2)]
    printf "%s\t%d\t%d\t%.2f\t%.2f\t%.2f\t%d\n", s, n, med, 100*low/n, 100*softbases/totalbases, 100*sp/n, nmed
  }' >> $OUT/Atp2a2_read_metrics.tsv
done
echo DONE
