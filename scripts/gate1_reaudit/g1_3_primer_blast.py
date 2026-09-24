#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Primer-BLAST 全自动：提交 → job_key 轮询 → 保存全部结果页"""
import re
import subprocess
import time

BASE = "https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi"
OUT = "/mnt/d/stroke_apa_reaudit/results/06_gate1_reaudit/primer_blast"
R = "/mnt/d/stroke_apa_reaudit/results/06_gate1_reaudit"

JOBS = [
    ("long_junction_refseqmrna", None),   # 已有结果文件
    ("long_junction_genome",
     dict(LEFT="GGAGTAACCGCTTCCTAAACCA", RIGHT="ATCCTACTATGCGGCAGAACAG",
          DB="refseq_representative_genomes")),
    ("qPCR_common_genome",
     dict(LEFT="AGCCTTTGTAGAGCCGTTTGTA", RIGHT="CACACTCTTTCTGTCCTGTCGA",
          DB="refseq_representative_genomes")),
    ("distal_terminal_genome",
     dict(LEFT="AGTTAGGACTGGAGGCCTATGT", RIGHT="TTGTAAGTGGCCAGATTGCTCT",
          DB="refseq_representative_genomes")),
    ("RACE_outer_GSP_mrna",
     dict(LEFT="GACAGGAATTAAAGCCTCAGGAAACC", RIGHT="", DB="refseq_mrna")),
    ("RACE_nested_GSP_mrna",
     dict(LEFT="TGAGGGCATTACACATCTCTATGGTT", RIGHT="", DB="refseq_mrna")),
]


def post(d, out):
    cmd = ["curl", "-s", "--max-time", "150", BASE, "-o", out, "-w", "%{http_code}"]
    code = subprocess.run(cmd, capture_output=True, text=True).stdout.strip()
    return code


def get(url, out):
    code = subprocess.run(["curl", "-s", "--max-time", "90", url, "-o", out, "-w", "%{http_code}"],
                          capture_output=True, text=True).stdout.strip()
    return code


def submit(primers, db):
    f = primers["LEFT"]
    r = primers.get("RIGHT", "")
    fields = [("CMD", "request"), ("ORGANISM", "10090"),
              ("PRIMER_LEFT_INPUT", f), ("SEARCH_SPECIFIC_PRIMER", "on"),
              ("PRIMER_SPECIFICITY_DATABASE", db),
              ("TOTAL_PRIMER_SPECIFICITY_MISMATCH", "4"),
              ("submit", "Get Primer Pairs")]
    if r:
        fields.insert(3, ("PRIMER_RIGHT_INPUT", r))
    cmd = ["curl", "-s", "--max-time", "150", BASE, "-o", "/tmp/pb_submit.html", "-w", "%{http_code}"]
    for k, v in fields:
        cmd += ["-F", f"{k}={v}"]
    code = subprocess.run(cmd, capture_output=True, text=True).stdout.strip()
    h = open("/tmp/pb_submit.html", encoding="utf-8", errors="replace").read()
    m = re.search(r'job_key=([A-Za-z0-9]+)', h)
    return (m.group(1) if m else None), code


def fetch(job_key, out):
    url = f"{BASE}?ctg_time=0&job_key={job_key}"
    for i in range(15):
        code = get(url, out)
        h = open(out, encoding="utf-8", errors="replace").read()
        if "job_key=" in h and "Status" in h and "Submitted" in h or "Waiting" in h:
            time.sleep(20)
            continue
        if "Primer-BLAST Results" in h or "No target templates" in h or "Primer pairs" in h:
            return True
        time.sleep(15)
    return False


def summarize(path):
    h = open(path, encoding="utf-8", errors="replace").read()
    t = re.sub(r"<script.*?</script>", "", h, flags=re.S)
    t = re.sub(r"<[^>]+>", "\n", t)
    lines = [l.strip() for l in t.split("\n") if l.strip()]
    keys = ("No target templates", "Primer Pair", "Primer pair", "product", "Specificity of primers",
            "target templates were found")
    out = []
    for i, l in enumerate(lines):
        if any(k.lower() in l.lower() for k in keys):
            out.append(l)
    # 去重保序
    seen, res = set(), []
    for l in out:
        if l not in seen:
            seen.add(l)
            res.append(l)
    return res[:20]


for name, primers in JOBS:
    out_html = f"{OUT}/{name}.html"
    if name == "long_junction_refseqmrna":
        subprocess.run(["cp", f"{R}/pb_junction_result.html", out_html])
        print(name, "-> 已有（No target templates in refseq_mrna，即 RefSeq mRNA 库零可扩增模板）")
        continue
    jk, code = submit(primers["LEFT"] if False else primers, primers["DB"])
    print(f"{name}: 提交 {code} job_key={jk}")
    if not jk:
        continue
    ok = fetch(jk, out_html)
    lines = summarize(out_html)
    print(f"  完成={ok}")
    for l in lines:
        print("   |", l)
print("全部完成；结果页保存在", OUT)
