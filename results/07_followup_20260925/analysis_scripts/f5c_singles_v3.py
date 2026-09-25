#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""单引物 Primer-BLAST 第三次尝试：RIGHT 字段存在但留空 + 补充 primer3 常用字段。"""
import re
import subprocess
import time
import os

BASE = "https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi"
OUT = "/mnt/d/stroke_apa_reaudit/results/07_followup_20260925/primer_blast_v2"
os.makedirs(OUT, exist_ok=True)

JOBS = [
    ("outerGSP_v2_single", "GGTTTCCTGAGGCTTTAATTCCTGTC"),
    ("nestedGSP_v2_single", "CGGTCCAAGAGTCTCCTTCTACCA"),
]

def submit(seq, tag):
    fields = [("CMD", "request"), ("ORGANISM", "10090"),
              ("PRIMER_LEFT_INPUT", seq), ("PRIMER_RIGHT_INPUT", ""),
              ("SEARCH_SPECIFIC_PRIMER", "on"),
              ("PRIMER_SPECIFICITY_DATABASE", "refseq_mrna"),
              ("TOTAL_PRIMER_SPECIFICITY_MISMATCH", "4"),
              ("PRIMER_PRODUCT_MIN", "50"), ("PRIMER_PRODUCT_MAX", "3000"),
              ("PRIMER_MIN_TM", "55"), ("PRIMER_OPT_TM", "60"),
              ("PRIMER_MAX_TM", "70"), ("PRIMER_PAIRS_NUMBER", "20")]
    cmd = ["curl", "-s", "--max-time", "150", BASE,
           "-o", f"/tmp/pb3_{tag}.html", "-w", "%{http_code}"]
    for k, v in fields:
        cmd += ["-F", f"{k}={v}"]
    code = subprocess.run(cmd, capture_output=True, text=True).stdout.strip()
    h = open(f"/tmp/pb3_{tag}.html", encoding="utf-8", errors="replace").read()
    m = re.search(r'job_key=([A-Za-z0-9\-]+)', h)
    return code, (m.group(1) if m else None), h

def fetch(job_key, out):
    url = f"{BASE}?job_key={job_key}"
    for _ in range(30):
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
    keys = ("No target templates", "Specificity of primers", "product length",
            "templates were found", "Primer Bound To")
    seen, res = set(), []
    for l in lines:
        if any(k.lower() in l.lower() for k in keys) and l not in seen:
            seen.add(l)
            res.append(l)
    return res[:14]

for name, seq in JOBS:
    code, jk, h = submit(seq, name)
    print(f"{name}: submit={code} job_key={jk}", flush=True)
    if jk:
        out_html = f"{OUT}/{name}.html"
        ok = fetch(jk, out_html)
        for l in summarize(out_html):
            print("   |", l)
    else:
        with open(f"{OUT}/{name}_submit_fail_v3.html", "w",
                  encoding="utf-8", errors="replace") as f:
            f.write(h)
        print("   -> 仍被退回表单，已存证")
print("DONE", flush=True)
