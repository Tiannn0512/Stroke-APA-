#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""f5: 对 v2 新 oligo 补跑 NCBI Primer-BLAST（GRCm39 RefSeq, 辅助性质;
mm10/vM25 本地核查为主——本地 ≤1 错配扫描见 f4）。"""
import re
import subprocess
import time

BASE = "https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi"
OUT = "/mnt/d/stroke_apa_reaudit/results/07_followup_20260925/primer_blast_v2"
import os
os.makedirs(OUT, exist_ok=True)

OUTER_V2 = "GGTTTCCTGAGGCTTTAATTCCTGTC"
NESTED_V2 = "CGGTCCAAGAGTCTCCTTCTACCA"
JUNC_F = "TGTGGTGTTTTCCTCCAATGCCT"
JUNC_R = "CACGCACCCGAACACCCTTATAT"
J_S = "TGGAACAACCCGCAATACTGGAGT"
J_L = "CCATCAACTAACCAATACTGGAGT"
DISTAL_F = "AGTTAGGACTGGAGGCCTATGT"

JOBS = [
    ("outerGSP_v2_single", dict(LEFT=OUTER_V2, RIGHT="", DB="refseq_mrna")),
    ("nestedGSP_v2_single", dict(LEFT=NESTED_V2, RIGHT="", DB="refseq_mrna")),
    ("long_junction_v2_pair", dict(LEFT=JUNC_F, RIGHT=JUNC_R,
                                   DB="refseq_representative_genomes")),
    ("qPCR_S_spliceform_pair", dict(LEFT=J_S, RIGHT=DISTAL_F,
                                    DB="refseq_representative_genomes")),
    ("qPCR_L_orthogonal_pair", dict(LEFT=J_L, RIGHT=DISTAL_F,
                                    DB="refseq_representative_genomes")),
]

def submit(primers):
    fields = [("CMD", "request"), ("ORGANISM", "10090"),
              ("PRIMER_LEFT_INPUT", primers["LEFT"]),
              ("SEARCH_SPECIFIC_PRIMER", "on"),
              ("PRIMER_SPECIFICITY_DATABASE", primers["DB"]),
              ("TOTAL_PRIMER_SPECIFICITY_MISMATCH", "4"),
              ("PRIMER_PRODUCT_MIN", "50"), ("PRIMER_PRODUCT_MAX", "3000")]
    if primers["RIGHT"]:
        fields.append(("PRIMER_RIGHT_INPUT", primers["RIGHT"]))
    cmd = ["curl", "-s", "--max-time", "150", BASE, "-o", "/tmp/pb2_submit.html",
           "-w", "%{http_code}"]
    for k, v in fields:
        cmd += ["-F", f"{k}={v}"]
    code = subprocess.run(cmd, capture_output=True, text=True).stdout.strip()
    h = open("/tmp/pb2_submit.html", encoding="utf-8", errors="replace").read()
    m = re.search(r'job_key=([A-Za-z0-9\-]+)', h)
    return (m.group(1) if m else None), code

def fetch(job_key, out):
    url = f"{BASE}?job_key={job_key}"
    for _ in range(40):
        time.sleep(15)
        code = subprocess.run(["curl", "-s", "--max-time", "90", url, "-o", out,
                               "-w", "%{http_code}"],
                              capture_output=True, text=True).stdout.strip()
        if code != "200":
            continue
        h = open(out, encoding="utf-8", errors="replace").read()
        if "meta http-equiv" in h and "Refresh" in h:
            continue
        return True
    return False

def summarize(path):
    h = open(path, encoding="utf-8", errors="replace").read()
    t = re.sub(r"<script.*?</script>", "", h, flags=re.S)
    t = re.sub(r"<[^>]+>", "\n", t)
    lines = [l.strip() for l in t.split("\n") if l.strip()]
    keys = ("No target templates", "Primer Pair", "Primer pair", "product",
            "Specificity of primers", "target templates were found")
    seen, res = set(), []
    for l in lines:
        if any(k.lower() in l.lower() for k in keys) and l not in seen:
            seen.add(l)
            res.append(l)
    return res[:16]

log = []
for name, p in JOBS:
    jk, code = submit(p)
    print(f"{name}: submit {code} job_key={jk}", flush=True)
    line = f"{name}\tsubmit={code}\tjob_key={jk}"
    if not jk:
        log.append(line + "\tFAILED")
        continue
    out_html = f"{OUT}/{name}.html"
    ok = fetch(jk, out_html)
    lines = summarize(out_html)
    print(f"  done={ok}")
    for l in lines:
        print("   |", l)
    log.append(line + f"\tdone={ok}\t" + " ;; ".join(lines)[:600])

with open(f"{OUT}/primer_blast_v2_summary.tsv", "w") as f:
    f.write("job\tsubmit\tdone\tsummary\n")
    for l in log:
        f.write(l.replace("\t", "\t", 2) + "\n")
print("ALL DONE", flush=True)
