"""Resume/retry failed chunks for sham1_1, then assemble. Reuses logic from p2_download_fastq."""
import os, time, urllib.request
from concurrent.futures import ThreadPoolExecutor

URL = "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/070/SRR25403070/SRR25403070_1.fastq.gz"
PARTS = r"D:\stroke_apa_data\fastq\sham1_1.fastq.gz.parts"
OUT = r"D:\stroke_apa_data\fastq\sham1_1.fastq.gz"
CHUNK = 8 * 1024 * 1024
WORKERS = 8

size = 2090117617  # from ENA filereport
ranges = [(s, min(s + CHUNK, size) - 1) for s in range(0, size, CHUNK)]

def fetch_range(rng):
    start, end = rng
    tmp = os.path.join(PARTS, f"x.p{start}")
    exp = end - start + 1
    if os.path.exists(tmp) and os.path.getsize(tmp) == exp:
        return None
    for att in range(8):
        try:
            req = urllib.request.Request(URL, headers={"Range": f"bytes={start}-{end}", "User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=300) as r, open(tmp, "wb") as f:
                while True:
                    b = r.read(1 << 16)
                    if not b:
                        break
                    f.write(b)
            if os.path.getsize(tmp) == exp:
                return start
            os.remove(tmp)
        except Exception:
            time.sleep(3 * (att + 1))
    raise RuntimeError(f"chunk fail {start}")

todo = [r for r in ranges if not (os.path.exists(os.path.join(PARTS, f"x.p{r[0]}")) and os.path.getsize(os.path.join(PARTS, f"x.p{r[0]}")) == r[1] - r[0] + 1)]
print("chunks to fetch:", len(todo), flush=True)
with ThreadPoolExecutor(max_workers=WORKERS) as ex:
    for res in ex.map(fetch_range, todo):
        pass
# verify
missing = [r for r in ranges if not (os.path.exists(os.path.join(PARTS, f"x.p{r[0]}")) and os.path.getsize(os.path.join(PARTS, f"x.p{r[0]}")) == r[1] - r[0] + 1)]
print("still missing:", len(missing), flush=True)
if missing:
    raise SystemExit(1)
with open(OUT + ".tmp", "wb") as fout:
    for s, e in ranges:
        with open(os.path.join(PARTS, f"x.p{s}"), "rb") as fi:
            fout.write(fi.read())
assert os.path.getsize(OUT + ".tmp") == size
os.replace(OUT + ".tmp", OUT)
for f in os.listdir(PARTS):
    os.remove(os.path.join(PARTS, f))
os.rmdir(PARTS)
print("ASSEMBLED:", OUT, os.path.getsize(OUT), flush=True)
