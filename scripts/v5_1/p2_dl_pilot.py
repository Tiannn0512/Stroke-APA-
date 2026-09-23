# -*- coding: utf-8 -*-
"""Resilient parallel downloader for GSE238125 pilot (4 runs x R1/R2, ENA FASTQ).
Resume-safe: 8MB chunks stored as .parts/x.p{start}; reuse legacy parts; assemble on completion.
Bandwidth concentrated on one file at a time; workers=24.
"""
import os, sys, time, urllib.request
from concurrent.futures import ThreadPoolExecutor

OUTDIR = r"D:\stroke_apa_data\fastq"
CHUNK = 8 * 1024 * 1024
WORKERS = 16
REQ_TIMEOUT = 600
BASE = "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/{sub}/{run}/{run}_{r}.fastq.gz"

# (run, sample_name) -- pilot: Sham rep1/2 + day1 rep1/2
PILOT = [("SRR25403070", "sham1"), ("SRR25403071", "sham2"),
         ("SRR25403072", "day1_rep1"), ("SRR25403073", "day1_rep2")]

UA = {"User-Agent": "Mozilla/5.0"}

def head_size(url):
    for att in range(6):
        try:
            req = urllib.request.Request(url, method="HEAD", headers=UA)
            with urllib.request.urlopen(req, timeout=60) as r:
                return int(r.headers["Content-Length"])
        except Exception:
            time.sleep(3 * (att + 1))
    raise RuntimeError(f"HEAD fail {url}")

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
            try:
                if os.path.exists(tmp):
                    os.remove(tmp)
            except OSError:
                pass
            time.sleep(min(2 * (att + 1), 30))
    return "FAIL"

total_done_bytes = 0

def run_one(url, out):
    global total_done_bytes
    if os.path.exists(out):
        print(f"[SKIP] {os.path.basename(out)} already complete", flush=True)
        return True
    size = head_size(url)
    ranges = [(s, min(s + CHUNK, size) - 1) for s in range(0, size, CHUNK)]
    pdir = out + ".parts"
    os.makedirs(pdir, exist_ok=True)
    print(f"[FILE] {os.path.basename(out)} {size/1e9:.2f}GB in {len(ranges)} chunks", flush=True)
    t0 = time.time()
    fails = []
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = {ex.submit(fetch_chunk, url, os.path.join(pdir, f"x.p{s}"), s, e): s for s, e in ranges}
        n = 0
        for fut in futs:
            pass
        for i, fut in enumerate(list(futs)):
            r = fut.result()
            n += 1
            if r == "FAIL":
                fails.append(futs[fut])
            if n % 15 == 0:
                got = sum(os.path.getsize(os.path.join(pdir, f)) for f in os.listdir(pdir) if f.endswith(".p") or ".p" in f)
                el = time.time() - t0
                print(f"  {n}/{len(ranges)} chunks, {got/1e6:.0f}MB, {got/1e6/max(el,1):.2f}MB/s", flush=True)
    if fails:
        print(f"[FAIL] {len(fails)} chunks failed: {fails[:5]}", flush=True)
        return False
    # assemble
    tmp_out = out + ".tmp"
    with open(tmp_out, "wb") as fo:
        for s, e in ranges:
            with open(os.path.join(pdir, f"x.p{s}"), "rb") as fi:
                while True:
                    b = fi.read(1 << 22)
                    if not b:
                        break
                    fo.write(b)
    if os.path.getsize(tmp_out) != size:
        print("[ERR] assembled size mismatch", flush=True)
        os.remove(tmp_out)
        return False
    os.replace(tmp_out, out)
    for f in os.listdir(pdir):
        os.remove(os.path.join(pdir, f))
    os.rmdir(pdir)
    print(f"[DONE] {os.path.basename(out)} {os.path.getsize(out)/1e9:.2f}GB in {time.time()-t0:.0f}s", flush=True)
    return True

for run, name in PILOT:
    for r in ("1", "2"):
        url = BASE.format(sub=run[-3:], run=run, r=r)
        out = os.path.join(OUTDIR, f"{name}_{r}.fastq.gz")
        ok = run_one(url, out)
        if not ok:
            print(f"[ABORT] {out} incomplete; rerun script to resume", flush=True)
            sys.exit(2)
print("ALL PILOT FILES DONE", flush=True)
