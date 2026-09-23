"""Parallel ranged downloader for ENA FASTQ (HTTPS). Resumable, chunk 8MB, 8 workers."""
import os, sys, time, urllib.request
from concurrent.futures import ThreadPoolExecutor

RUNS = {
    "SRR25403070": "sham1", "SRR25403071": "sham2",
    "SRR25403072": "day1_rep1", "SRR25403073": "day1_rep2",
}
BASE = "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/{sub}/{run}/{run}_{r}.fastq.gz"
OUTDIR = r"D:\stroke_apa_data\fastq"
CHUNK = 8 * 1024 * 1024
WORKERS = 8

def head_size(url):
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return int(r.headers["Content-Length"])

def fetch_range(url, path, start, end):
    tmp = path + f".p{start}"
    exp = end - start + 1
    if os.path.exists(tmp) and os.path.getsize(tmp) == exp:
        return None
    for att in range(6):
        try:
            req = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end}", "User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=300) as r, open(tmp, "wb") as f:
                while True:
                    b = r.read(1 << 16)
                    if not b:
                        break
                    f.write(b)
            if os.path.getsize(tmp) == exp:
                return tmp
            os.remove(tmp)
        except Exception:
            time.sleep(2 * (att + 1))
    raise RuntimeError(f"chunk fail {start}")

for run, name in RUNS.items():
    for r in ("1", "2"):
        url = BASE.format(sub=run[-3:], run=run, r=r)
        out = os.path.join(OUTDIR, f"{name}_{r}.fastq.gz")
        if os.path.exists(out):
            print("exists, skip", out, flush=True)
            continue
        size = head_size(url)
        print(f"{run}_{r}: {size/1e9:.2f} GB", flush=True)
        ranges = [(s, min(s + CHUNK, size) - 1) for s in range(0, size, CHUNK)]
        parts_dir = out + ".parts"
        os.makedirs(parts_dir, exist_ok=True)
        t0 = time.time()
        with ThreadPoolExecutor(max_workers=WORKERS) as ex:
            futs = [ex.submit(fetch_range, url, os.path.join(parts_dir, "x"), s, e) for s, e in ranges]
            for i, f in enumerate(futs):
                f.result()
                if (i + 1) % 20 == 0:
                    print(f"  {i+1}/{len(ranges)} chunks, {time.time()-t0:.0f}s", flush=True)
        with open(out + ".tmp", "wb") as fout:
            for s, e in ranges:
                with open(os.path.join(parts_dir, f"x.p{s}"), "rb") as fi:
                    fout.write(fi.read())
        assert os.path.getsize(out + ".tmp") == size
        os.replace(out + ".tmp", out)
        for f in os.listdir(parts_dir):
            os.remove(os.path.join(parts_dir, f))
        os.rmdir(parts_dir)
        print(f"DONE {out} {os.path.getsize(out)/1e9:.2f}GB in {time.time()-t0:.0f}s", flush=True)
print("ALL DONE")
