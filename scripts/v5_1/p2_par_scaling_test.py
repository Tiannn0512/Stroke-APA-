# -*- coding: utf-8 -*-
"""60s parallel-scaling test on ENA fastq URL: 32 workers x 2MB ranges, report aggregate MB/s."""
import os, time, urllib.request, threading

URL = "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/070/SRR25403070/SRR25403070_2.fastq.gz"
OUT = r"D:\stroke_apa_data\fastq\_par_test.bin"
NWORK, CHUNK, WINDOW = 32, 2 * 1024 * 1024, 60

lock = threading.Lock()
done_bytes = 0
errors = 0
t_end = time.time() + WINDOW

def worker(wid):
    global done_bytes, errors
    off = wid * CHUNK
    while time.time() < t_end:
        try:
            req = urllib.request.Request(URL, headers={"Range": f"bytes={off}-{off+CHUNK-1}", "User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=30) as r:
                while time.time() < t_end:
                    b = r.read(1 << 16)
                    if not b:
                        break
                    with lock:
                        done_bytes += len(b)
            off += NWORK * CHUNK
        except Exception:
            with lock:
                errors += 1
            time.sleep(1)

threads = [threading.Thread(target=worker, args=(i,), daemon=True) for i in range(NWORK)]
t0 = time.time()
for t in threads: t.start()
for t in threads: t.join(timeout=WINDOW + 15)
dt = max(time.time() - t0, 1e-6)
print(f"aggregate: {done_bytes/1e6:.1f} MB in {dt:.0f}s = {done_bytes/1e6/dt:.2f} MB/s with {NWORK} conns, errors={errors}")
try:
    os.remove(OUT)
except OSError:
    pass
