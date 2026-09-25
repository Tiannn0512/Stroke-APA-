#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""反向定位 S penult 末端的真实基因组坐标。"""
from pyfaidx import Fasta

BASE = "/mnt/d/stroke_apa_reaudit"
GENOME = Fasta(f"{BASE}/reference/chr5.fa")["chr5"]

SEQS, cur = {}, None
with open(f"{BASE}/reference/gencode.vM25.transcripts.fa") as fh:
    for line in fh:
        if line.startswith(">"):
            cur = line[1:].split("|")[0].strip(); SEQS[cur] = []
        else:
            SEQS[cur].append(line.strip())
SEQS = {k: "".join(v).upper() for k, v in SEQS.items()}
S_ = "ENSMUST00000177974.7"
s = SEQS[S_]

def rc(x):
    return x[::-1].translate(str.maketrans("ACGT", "TGCA"))

# S 的 penult = [3069, 3193)
penult = s[3069:3193]
print("penult len:", len(penult))
tail20 = penult[-20:]
print("S penult tail20:", tail20)
probe = rc(tail20)
g = str(GENOME[122450000:122470000]).upper()
i = g.find(probe)
print("rc(tail20) in genome win at:", 122450000 + i if i >= 0 else -1)
if i >= 0:
    print("=> tail20 对应基因组:", 122450000 + i, "..", 122450000 + i + 19,
          "(mRNA 首碱基在高坐标端)")
# 对比不同切片
print("rc(G[122457290:122457303]):", rc(str(GENOME[122457290:122457303]).upper()))
print("rc(G[122457291:122457304]):", rc(str(GENOME[122457291:122457304]).upper()))
print("rc(G[122457303:122457316]):", rc(str(GENOME[122457303:122457316]).upper()))
print("S[3180:3193]             :", s[3180:3193])
# 同样反向核对 terminal 首 13
head13 = s[3193:3206]
print("S terminal head13:", head13, " genomic probe idx:",
      122450000 + g.find(rc(head13)) if g.find(rc(head13)) >= 0 else -1)
