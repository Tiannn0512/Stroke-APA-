#!/usr/bin/env python3
"""P2.3a: build length+GC matched null set from the null pool (seed 42).

For each candidate segment: choose an unused pool segment with minimal |GC diff|
among those with length within +-10% (fallback +-25%, then nearest length).
Deterministic: pool shuffled with random.Random(42) before greedy assignment.
Writes null_matched.fa (simple N#-style ids kept from pool) + null_matched_mapping.tsv.
"""
import csv, os, random

BASE = r"D:\stroke_apa_reanalysis"
OUT = os.path.join(BASE, "results", "02_candidate_rebuild")

cand = list(csv.DictReader(open(os.path.join(OUT, "distal_candidates_mapping.tsv"), encoding="utf-8"), delimiter="\t"))
pool = list(csv.DictReader(open(os.path.join(OUT, "distal_nullpool_mapping.tsv"), encoding="utf-8"), delimiter="\t"))

rng = random.Random(42)
rng.shuffle(pool)

def pick(c):
    cl, cg = int(c["len"]), float(c["gc_pct"])
    for win, key in ((0.10, lambda p: abs(float(p["gc_pct"]) - cg)),
                     (0.25, lambda p: abs(float(p["gc_pct"]) - cg)),
                     (1e9,  lambda p: abs(int(p["len"]) - cl))):
        feas = [p for p in pool if not p.get("_used") and abs(int(p["len"]) - cl) <= max(win * cl, 20)]
        if feas:
            best = min(feas, key=key)
            best["_used"] = True
            return best
    return None

matched = []
for c in cand:
    n = pick(c)
    if n is None:
        print("NO MATCH for", c["seq_id"])
        continue
    matched.append((c, n))

# write matched null fasta in candidate order
pool_fa = {}
sid = None
with open(os.path.join(OUT, "distal_nullpool.fa"), encoding="utf-8") as f:
    for line in f:
        if line.startswith(">"):
            sid = line[1:].strip()
            pool_fa[sid] = []
        else:
            pool_fa[sid].append(line.strip())

with open(os.path.join(OUT, "null_matched.fa"), "w", encoding="utf-8", newline="") as f:
    for c, n in matched:
        f.write(f">{c['seq_id']}_null\n")
        s = "".join(pool_fa[n["seq_id"]])
        for j in range(0, len(s), 60):
            f.write(s[j:j+60] + "\n")

with open(os.path.join(OUT, "null_matched_mapping.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["cand_seq_id", "cand_len", "cand_gc", "null_seq_id", "null_event_num", "null_event_id",
                "null_len", "null_gc", "len_diff", "gc_diff"])
    for c, n in matched:
        w.writerow([c["seq_id"], c["len"], c["gc_pct"], n["seq_id"], n["event_num"], n["event_id"],
                    n["len"], n["gc_pct"], abs(int(c["len"]) - int(n["len"])),
                    abs(float(c["gc_pct"]) - float(n["gc_pct"]))])

import statistics
dl = [abs(int(c["len"]) - int(n["len"])) for c, n in matched]
dg = [abs(float(c["gc_pct"]) - float(n["gc_pct"])) for c, n in matched]
print(f"matched: {len(matched)}/{len(cand)}; median len diff {statistics.median(dl)}; median GC diff {statistics.median(dg):.2f}pp")
