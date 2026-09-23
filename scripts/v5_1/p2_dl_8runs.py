# -*- coding: utf-8 -*-
"""P2-05h: download remaining 8 GSE238125 runs (day3/7/21/60, ~32GB) from ENA direct.
Chunk-resume design identical to p2_dl_pilot.py. Launched by dl_8runs_watcher.sh AFTER the
pilot downloader finishes (single pipe ~0.8-1.2MB/s; serial > parallel for ENA).
Route decision 2026-09-17 09:20: ENA direct 1.22MB/s > proxy socks5h://10808 (47KB/s)
> AWS S3 (12-82KB/s) > Azure (dead 409) > GCS (404) > NGDC (no SRR mirror)."""
import os, sys, time, urllib.request
from concurrent.futures import ThreadPoolExecutor

OUTDIR = r"D:\stroke_apa_data\fastq"
CHUNK = 8 * 1024 * 1024
WORKERS = 16
REQ_TIMEOUT = 600
BASE = "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/{sub}/{run}/{run}_{r}.fastq.gz"

RUNS = [("SRR25403062", "day3_rep1"), ("SRR25403063", "day3_rep2"),
        ("SRR25403064", "day7_rep1"), ("SRR25403065", "day7_rep2"),
        ("SRR25403066", "day21_rep1"), ("SRR25403067", "day21_rep2"),
        ("SRR25403068", "day60_rep1"), ("SRR25403069", "day60_rep2")]

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
            if os.path.exists(tmp):
                os.remove(tmp)
            time.sleep(min(2 * (att + 1), 30))
    return "FAIL"

def run_one(url, out):
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
        for i, fut in enumerate(list(futs)):
            if fut.result() == "FAIL":
                fails.append(futs[fut])
            if (i + 1) % 30 == 0:
                got = sum(os.path.getsize(os.path.join(pdir, f)) for f in os.listdir(pdir))
                el = time.time() - t0
                print(f"  {i+1}/{len(ranges)} chunks, {got/1e6:.0f}MB, {got/1e6/max(el,1):.2f}MB/s", flush=True)
    if fails:
        print(f"[FAIL] {len(fails)} chunks failed: {fails[:5]}", flush=True)
        return False
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

for run, name in RUNS:
    for r in ("1", "2"):
        url = BASE.format(sub=run[-3:], run=run, r=r)
        out = os.path.join(OUTDIR, f"{name}_{r}.fastq.gz")
        if not run_one(url, out):
            print(f"[ABORT] {out} incomplete; rerun script to resume", flush=True)
            sys.exit(2)
print("ALL 8RUNS FILES DONE", flush=True)
