"""Parallel ranged downloader for GSE174574_RAW.tar (NCBI supports Range)."""
import os
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor

URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174574/suppl/GSE174574_RAW.tar"
OUT = r"D:\stroke_apa\data\GSE174574\GSE174574_RAW.tar"
CHUNK = 8 * 1024 * 1024          # 8 MB per part
WORKERS = 8
TOTAL = 263_086_080              # from filelist.txt

parts_dir = OUT + ".parts"
os.makedirs(parts_dir, exist_ok=True)

ranges = []
for start in range(0, TOTAL, CHUNK):
    end = min(start + CHUNK, TOTAL) - 1
    ranges.append((start, end))

def fetch(rng):
    start, end = rng
    path = os.path.join(parts_dir, f"{start:012d}.part")
    if os.path.exists(path) and os.path.getsize(path) == end - start + 1:
        return f"skip {start}"
    req = urllib.request.Request(URL, headers={"Range": f"bytes={start}-{end}", "User-Agent": "Mozilla/5.0"})
    last_err = None
    for attempt in range(6):
        try:
            t0 = time.time()
            with urllib.request.urlopen(req, timeout=300) as r, open(path + ".tmp", "wb") as f:
                while True:
                    b = r.read(1 << 16)
                    if not b:
                        break
                    f.write(b)
            got = os.path.getsize(path + ".tmp")
            if got != end - start + 1:
                os.remove(path + ".tmp")
                last_err = f"SHORT {start}: {got}/{end-start+1}"
                time.sleep(2 * (attempt + 1))
                continue
            os.replace(path + ".tmp", path)
            dt = time.time() - t0
            return f"ok {start} ({got/1e6:.1f}MB in {dt:.1f}s = {got/1e6/dt:.2f}MB/s)"
        except Exception as e:
            last_err = f"retry {start} attempt {attempt+1}: {type(e).__name__}"
            time.sleep(3 * (attempt + 1))
    return f"FAIL {start}: {last_err}"

t0 = time.time()
done = 0
with ThreadPoolExecutor(max_workers=WORKERS) as ex:
    for msg in ex.map(fetch, ranges):
        done += 1
        if done % 5 == 0 or not msg.startswith(("ok", "skip")):
            print(f"[{done}/{len(ranges)}] {msg}", flush=True)

# verify all parts
missing = [r for r in ranges if not (os.path.exists(os.path.join(parts_dir, f"{r[0]:012d}.part"))
           and os.path.getsize(os.path.join(parts_dir, f"{r[0]:012d}.part")) == r[1] - r[0] + 1)]
print("missing parts:", len(missing), flush=True)
if missing:
    raise SystemExit(f"retry needed for {len(missing)} parts")

# assemble
with open(OUT + ".new", "wb") as out:
    for start, end in ranges:
        with open(os.path.join(parts_dir, f"{start:012d}.part"), "rb") as f:
            out.write(f.read())
size = os.path.getsize(OUT + ".new")
print(f"assembled {size} bytes (expect {TOTAL})", flush=True)
if size == TOTAL:
    if os.path.exists(OUT):
        os.remove(OUT)
    os.rename(OUT + ".new", OUT)
    print(f"DONE: {OUT}", flush=True)
else:
    raise SystemExit("size mismatch")
print(f"total wall time {time.time()-t0:.0f}s", flush=True)
