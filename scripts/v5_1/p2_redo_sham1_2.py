# -*- coding: utf-8 -*-
"""One-off re-download of corrupt sham1_2 (SRR25403070_2), 8 workers, independent paths.
Run alongside p2_dl_pilot.py (16 workers) — combined 24 connections, acceptable.
Evidence: assembled 2026-09-17 00:00 size-exact but md5 e1badc6b... != expected 85d69dbe...
"""
import os, time, urllib.request
from concurrent.futures import ThreadPoolExecutor

OUT = r"D:\stroke_apa_data\fastq\sham1_2.fastq.gz"
CHUNK = 8 * 1024 * 1024
WORKERS = 8
REQ_TIMEOUT = 600
URL = "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/070/SRR25403070/SRR25403070_2.fastq.gz"
UA = {"User-Agent": "Mozilla/5.0"}

def head_size(url):
    for att in range(6):
        try:
            req = urllib.request.Request(url, method="HEAD", headers=UA)
            with urllib.request.urlopen(req, timeout=60) as r:
                return int(r.headers["Content-Length"])
        except Exception:
            time.sleep(3 * (att + 1))
    raise RuntimeError("HEAD fail")

def fetch_chunk(url, path, start, end):
    exp = end - start + 1
    if os.path.exists(path) and os.path.getsize(path) == exp:
        return "skip"
    tmp = path + ".tmp"
    for att in range(10):
        try:
            req = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end}", **UA})
            with urllib.request.urlopen(req, timeout=REQ_TIMEOUT) as r, open(tmp, "wb") as f:
                while True:
                    b = r.read(1 << 16)
                    if not b:
                        break
                    f.write(b)
            if os.path.getsize(tmp) == exp:
                os.replace(tmp, path)
                return "ok"
            os.remove(tmp)
        except Exception:
            if os.path.exists(tmp):
                os.remove(tmp)
            time.sleep(min(2 * (att + 1), 30))
    return "FAIL"

if os.path.exists(OUT):
    raise SystemExit("[ABORT] output already exists — delete corrupt file first")
size = head_size(URL)
ranges = [(s, min(s + CHUNK, size) - 1) for s in range(0, size, CHUNK)]
pdir = OUT + ".redo.parts"
os.makedirs(pdir, exist_ok=True)
print(f"[FILE] {os.path.basename(OUT)} {size:,}B in {len(ranges)} chunks", flush=True)
t0 = time.time()
fails = []
with ThreadPoolExecutor(max_workers=WORKERS) as ex:
    futs = {ex.submit(fetch_chunk, URL, os.path.join(pdir, f"x.p{s}"), s, e): s for s, e in ranges}
    for i, fut in enumerate(list(futs)):
        if fut.result() == "FAIL":
            fails.append(futs[fut])
if fails:
    print(f"[FAIL] {len(fails)} chunks: {fails[:5]}", flush=True)
    raise SystemExit(2)
tmp_out = OUT + ".tmp"
with open(tmp_out, "wb") as fo:
    for s, e in ranges:
        with open(os.path.join(pdir, f"x.p{s}"), "rb") as fi:
            while True:
                b = fi.read(1 << 22)
                if not b:
                    break
                fo.write(b)
assert os.path.getsize(tmp_out) == size, "assembled size mismatch"
os.replace(tmp_out, OUT)
for f in os.listdir(pdir):
    os.remove(os.path.join(pdir, f))
os.rmdir(pdir)
print(f"[DONE] {os.path.basename(OUT)} re-downloaded clean in {time.time()-t0:.0f}s", flush=True)
