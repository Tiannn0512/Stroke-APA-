#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""补跑两个 GSP 单引物 Primer-BLAST；失败时落盘提交页以便诊断。"""
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
              ("PRIMER_LEFT_INPUT", seq),
              ("SEARCH_SPECIFIC_PRIMER", "on"),
              ("PRIMER_SPECIFICITY_DATABASE", "refseq_mrna"),
              ("TOTAL_PRIMER_SPECIFICITY_MISMATCH", "4"),
              ("PRIMER_PRODUCT_MIN", "50"), ("PRIMER_PRODUCT_MAX", "3000")]
    cmd = ["curl", "-s", "--max-time", "150", BASE,
           "-o", f"/tmp/pb2r_{tag}.html", "-w", "%{http_code}"]
    for k, v in fields:
        cmd += ["-F", f"{k}={v}"]
    code = subprocess.run(cmd, capture_output=True, text=True).stdout.strip()
    h = open(f"/tmp/pb2r_{tag}.html", encoding="utf-8", errors="replace").read()
    m = re.search(r'job_key=([A-Za-z0-9\-]+)', h)
    return code, (m.group(1) if m else None), h

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
    keys = ("No target templates", "Specificity of primers", "product length",
            "target templates were found", "primer pair", "Primer Bound To")
    seen, res = set(), []
    for l in lines:
        if any(k.lower() in l.lower() for k in keys) and l not in seen:
            seen.add(l)
            res.append(l)
    return res[:14]

for name, seq in JOBS:
    ok_all = False
    for attempt in range(3):
        code, jk, h = submit(seq, name)
        print(f"{name}: attempt{attempt} submit={code} job_key={jk}", flush=True)
        if jk:
            out_html = f"{OUT}/{name}.html"
            ok = fetch(jk, out_html)
            for l in summarize(out_html):
                print("   |", l)
            ok_all = ok
            break
        else:
            with open(f"{OUT}/{name}_submit_fail.html", "w",
                      encoding="utf-8", errors="replace") as f:
                f.write(h)
            time.sleep(30)
    print(f"{name} final={ok_all}", flush=True)
print("DONE", flush=True)
