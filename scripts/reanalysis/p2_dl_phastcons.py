#!/usr/bin/env python3
"""Parallel Range downloader for mm10.60way.phastCons.bw (UCSC). 24 threads, 8MB chunks, retries."""
import os, sys, threading, urllib.request, queue, time

URL = "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/phastCons60way/mm10.60way.phastCons.bw"
OUT = r"D:\stroke_apa_reanalysis\data\references\mm10.60way.phastCons.bw"
CHUNK = 8 * 1024 * 1024
THREADS = 24

def get_size():
    req = urllib.request.Request(URL, method="HEAD")
    with urllib.request.urlopen(req, timeout=60) as r:
        return int(r.headers["Content-Length"])

def fetch(part_path, start, end, retries=6):
    for attempt in range(retries):
        try:
            have = os.path.getsize(part_path) if os.path.exists(part_path) else 0
            if have >= end - start + 1:
                return
            req = urllib.request.Request(URL, headers={"Range": f"bytes={start+have}-{end}"})
            with urllib.request.urlopen(req, timeout=120) as r, open(part_path, "ab") as f:
                while True:
                    b = r.read(1 << 18)
                    if not b:
                        break
                    f.write(b)
            if os.path.getsize(part_path) >= end - start + 1:
                return
        except Exception as e:
            time.sleep(3 * (attempt + 1))
    raise RuntimeError(f"chunk {start}-{end} failed after {retries} retries")

def main():
    total = get_size()
    n_chunks = (total + CHUNK - 1) // CHUNK
    print(f"total {total} bytes, {n_chunks} chunks", flush=True)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)

    def part_size(i):
        p = f"{OUT}.part{i:04d}"
        return os.path.getsize(p) if os.path.exists(p) else 0

    def complete(i):
        start = i * CHUNK
        return part_size(i) >= min(start + CHUNK, total) - start

    t0 = time.time()
    for round_ in range(4):
        missing = [i for i in range(n_chunks) if not complete(i)]
        print(f"round {round_}: {len(missing)} chunks to fetch", flush=True)
        if not missing:
            break
        q = queue.Queue()
        for i in missing:
            q.put(i)
        done, lock = [0], threading.Lock()
        def worker():
            while True:
                try:
                    i = q.get_nowait()
                except queue.Empty:
                    return
                start = i * CHUNK
                end = min(start + CHUNK, total) - 1
                try:
                    fetch(f"{OUT}.part{i:04d}", start, end)
                except Exception as e:
                    print(f"chunk {i} failed: {e}", flush=True)
                with lock:
                    done[0] += 1
                    if done[0] % 10 == 0 or done[0] == len(missing):
                        print(f"  round {round_}: {done[0]}/{len(missing)}", flush=True)
        threads = [threading.Thread(target=worker) for _ in range(THREADS)]
        for t in threads: t.start()
        for t in threads: t.join()

    missing = [i for i in range(n_chunks) if not complete(i)]
    if missing:
        raise RuntimeError(f"chunks still incomplete after rounds: {missing[:20]}")

    with open(OUT, "wb") as out:
        for i in range(n_chunks):
            with open(f"{OUT}.part{i:04d}", "rb") as p:
                while True:
                    b = p.read(1 << 22)
                    if not b:
                        break
                    out.write(b)
    for i in range(n_chunks):
        os.remove(f"{OUT}.part{i:04d}")
    assert os.path.getsize(OUT) == total, "final size mismatch"
    print(f"DONE {OUT} {os.path.getsize(OUT)} bytes in {time.time()-t0:.0f}s", flush=True)

if __name__ == "__main__":
    main()
