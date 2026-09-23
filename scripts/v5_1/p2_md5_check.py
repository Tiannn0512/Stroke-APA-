# -*- coding: utf-8 -*-
"""MD5 verification of downloaded GSE238125 FASTQ files against ENA fastq_md5."""
import hashlib, os, sys, csv

MD5_TSV = r"D:\stroke_apa_data\fastq_md5_all.tsv"
FASTQ_DIR = r"D:\stroke_apa_data\fastq"
MAP_TSV = r"D:\stroke_apa\scripts\p2_sample_map.tsv"

# name -> srr
name2srr = {}
with open(MAP_TSV, encoding="utf-8") as f:
    next(f)
    for line in f:
        if line.strip():
            name, srr = line.rstrip("\n").split("\t")
            name2srr[name] = srr

# srr -> [md5_r1, md5_r2]
srr2md5 = {}
with open(MD5_TSV, encoding="utf-8") as f:
    next(f)
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 3 and parts[0]:
            srr2md5[parts[0]] = parts[2].split(";")

def md5sum(path, buf=1 << 22):
    h = hashlib.md5()
    with open(path, "rb") as f:
        while True:
            b = f.read(buf)
            if not b:
                break
            h.update(b)
    return h.hexdigest()

n_ok = n_bad = n_missing = 0
for name, srr in sorted(name2srr.items()):
    for r, md5_exp in zip(("1", "2"), srr2md5.get(srr, ["?", "?"])):
        path = os.path.join(FASTQ_DIR, f"{name}_{r}.fastq.gz")
        if not os.path.exists(path):
            print(f"[MISSING] {name}_{r}")
            n_missing += 1
            continue
        got = md5sum(path)
        status = "OK" if got == md5_exp else "MISMATCH"
        if status == "OK":
            n_ok += 1
        else:
            n_bad += 1
        print(f"[{status}] {name}_{r} {os.path.getsize(path):,}B md5={got}")
print(f"\nsummary: ok={n_ok} bad={n_bad} missing={n_missing}")
sys.exit(1 if (n_bad or n_missing) else 0)
