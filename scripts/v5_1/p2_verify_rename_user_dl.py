# -*- coding: utf-8 -*-
"""Verify and rename user-downloaded SRR-named FASTQ files.
Map: 064=day7_rep1, 065=day7_rep2, 066=day21_rep1, 067=day21_rep2, 068=day60_rep1, 069=day60_rep2
md5 reference: D:\stroke_apa_data\fastq_md5_all.tsv (rows: SRR  url1;url2  md5_1;md5_2)
"""
import os, hashlib, shutil

D = r"D:\stroke_apa_data\fastq"
MAP = {"SRR25403064": "day7_rep1", "SRR25403065": "day7_rep2",
       "SRR25403066": "day21_rep1", "SRR25403067": "day21_rep2",
       "SRR25403068": "day60_rep1", "SRR25403069": "day60_rep2"}

md5_ref = {}
with open(r"D:\stroke_apa_data\fastq_md5_all.tsv", encoding="utf-8") as f:
    next(f)
    for ln in f:
        parts = ln.rstrip("\n").split("\t")
        if len(parts) >= 3 and parts[0]:
            md5_ref[parts[0]] = parts[2].split(";")

def md5sum(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()

candidates = [fn for fn in os.listdir(D) if fn.startswith("SRR2540306")]
print("candidate files:", candidates)
for fn in sorted(candidates):
    src = os.path.join(D, fn)
    base = fn.replace(" (1)", "").replace(" (2)", "")
    srr = base.split("_")[0]
    r = "1" if "_1.fastq.gz" in base else "2"
    name = MAP.get(srr)
    if not name:
        print(f"[SKIP] {fn}: unknown SRR")
        continue
    exp = md5_ref.get(srr, ["?", "?"])[int(r) - 1]
    got = md5sum(src)
    if got == exp:
        dst = os.path.join(D, f"{name}_{r}.fastq.gz")
        if os.path.exists(dst):
            print(f"[DUP] {fn} md5 OK but {name}_{r} already exists -> delete duplicate")
            os.remove(src)
        else:
            os.rename(src, dst)
            print(f"[OK+RENAMED] {fn} -> {name}_{r}.fastq.gz")
    else:
        print(f"[BAD] {fn}: md5 mismatch (got {got[:8]}, want {exp[:8]}) -> quarantined")
        os.rename(src, src + ".badmd5")

print("\n=== final inventory ===")
for r_ in sorted(MAP.values()):
    for r in ("1", "2"):
        p = os.path.join(D, f"{r_}_{r}.fastq.gz")
        print(("OK  " if os.path.exists(p) else "MISS"), f"{r_}_{r}")
