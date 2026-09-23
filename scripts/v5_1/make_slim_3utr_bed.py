"""Slim DaPars2 annotation: gene-level longest-3'UTR representative.
Rationale: DaPars2 tandem model fits ONE proximal/distal PAS pair per gene region;
keeping every transcript's UTR duplicates genes many times (52,333 rows -> ~1 per gene).
Keeps the longest UTR per gene symbol key (gene_id|symbol field), preserves format.
Output: /home/taylor/reference_v2/gencode_M25_3UTR_slim.bed (via WSL ext4)
"""
from collections import defaultdict

src = "/home/taylor/reference_v2/gencode_M25_3UTR_for_DaPars2.bed"
dst = "/home/taylor/reference_v2/gencode_M25_3UTR_slim.bed"

best = {}  # gene_key -> (len, line)
n = 0
for line in open(src):
    f = line.rstrip("\n").split("\t")
    if len(f) < 6:
        continue
    n += 1
    gene = f[3].split("|", 1)[1] if "|" in f[3] else f[3]
    L = int(f[2]) - int(f[1])
    if gene not in best or L > best[gene][0]:
        best[gene] = (L, line)

with open(dst, "w") as fo:
    for gene, (L, line) in best.items():
        fo.write(line)
print(f"input rows: {n} | slim rows (gene-level): {len(best)} -> {dst}")
