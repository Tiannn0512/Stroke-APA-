"""Replicate DaPars2 load_wig logic for chr1 on the real bedgraph to find why events miss coverage."""
import sys

wig = "/home/taylor/stroke_APA_work/bedgraph/sham1.bedgraph"
CHR = sys.argv[1] if len(sys.argv) > 1 else "chr1"
d = {}
broke_at = None
with open(wig) as fin:
    for line in fin:
        if line[0] == "#" or line[0] == "t":
            continue
        f = line.split("\t")
        chrom = f[0]
        if chrom == CHR:
            if chrom not in d:
                d[chrom] = [[0], [0]]
            s, e = int(f[1]), int(f[2])
            if s > d[chrom][0][-1]:
                d[chrom][0].append(s)
                d[chrom][1].append(0)
            d[chrom][0].append(e)
            d[chrom][1].append(int(float(f[-1])))
        else:
            if len(d) > 0:
                broke_at = chrom
                break
print("dict chroms:", list(d.keys()))
print("broke at:", broke_at)
print("chr1 collected points:", len(d.get(CHR, [[], []])[0]) if CHR in d else 0)
# Xkr4 UTR region probe: chr1:3214481-3216024
if CHR in d:
    starts, covs = d[CHR]
    import bisect
    li = bisect.bisect(starts, 3214481)
    ri = bisect.bisect(starts, 3216024)
    print(f"Xkr4 region idx {li}:{ri} cov points {ri-li+1}")
