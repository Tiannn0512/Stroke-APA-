# -*- coding: utf-8 -*-
"""P2-05 preflight: ENA filereport for PRJNA997998 (GSE238125), GEO suppl listing, ENA speed test."""
import os, re, sys, time, json, urllib.request

BASE = r"D:\stroke_apa\data\GSE238125"
os.makedirs(BASE, exist_ok=True)
UA = {"User-Agent": "Mozilla/5.0 (OpenClaw preflight)"}

API = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA997998"
       "&result=read_run&fields=run_accession,sample_accession,experiment_title,"
       "library_name,read_count,fastq_ftp,fastq_bytes,fastq_md5&format=tsv&limit=0")

def fetch(url, headers=None, timeout=90):
    req = urllib.request.Request(url, headers=headers or UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()

# --- 1. ENA filereport ---
try:
    t0 = time.time()
    data = fetch(API).decode("utf-8", errors="replace")
    print(f"[ENA API] OK, {len(data)} bytes in {time.time()-t0:.1f}s")
    out = os.path.join(BASE, "ena_filereport_PRJNA997998.tsv")
    with open(out, "w", encoding="utf-8", newline="") as f:
        f.write(data)
    print(f"[SAVED] {out}")
except Exception as e:
    print(f"[ENA API] FAIL: {type(e).__name__}: {e}")
    sys.exit(1)

lines = [l for l in data.strip().split("\n") if l.strip()]
hdr = lines[0].split("\t")
rows = [dict(zip(hdr, l.split("\t"))) for l in lines[1:]]
print(f"[RUNS] {len(rows)}")

total = 0
for r in rows:
    fb = r.get("fastq_bytes", "")
    nbytes = sum(int(x) for x in fb.split(";") if x.strip().isdigit()) if fb else 0
    total += nbytes
    ftp = r.get("fastq_ftp", "")
    first = ftp.split(";")[0] if ftp else "(no fastq)"
    print(f"  {r.get('run_accession','?')}\t{r.get('library_name','?')}\t{nbytes/1e9:6.2f} GB\t{first}")
print(f"[TOTAL] {total/1e9:.2f} GB across {len(rows)} runs")

# --- 2. GEO supplementary listing ---
try:
    html = fetch("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE238nnn/GSE238125/suppl/", timeout=60).decode(errors="replace")
    files = [h for h in re.findall(r'href="([^"?][^"]*)"', html) if h not in ("..",) and not h.endswith("/")]
    print(f"[GEO SUPPL] {len(files)} files: {files[:20]}")
except Exception as e:
    print(f"[GEO SUPPL] FAIL: {type(e).__name__}: {e}")

# --- 3. speed test: first 4 MB of run #1 fastq_1 ---
if rows and rows[0].get("fastq_ftp"):
    url = "https://" + rows[0]["fastq_ftp"].split(";")[0]
    try:
        req = urllib.request.Request(url, headers={**UA, "Range": "bytes=0-4194303"})
        t0 = time.time(); n = 0
        with urllib.request.urlopen(req, timeout=120) as r:
            while True:
                chunk = r.read(1 << 20)
                if not chunk:
                    break
                n += len(chunk)
        dt = max(time.time() - t0, 1e-6)
        print(f"[SPEED] {n/1e6:.1f} MB in {dt:.1f}s = {n/1e6/dt:.2f} MB/s  <- {url[:90]}")
    except Exception as e:
        print(f"[SPEED] FAIL: {type(e).__name__}: {e}")
else:
    print("[SPEED] skipped: no fastq_ftp in report")
