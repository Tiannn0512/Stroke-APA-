#!/usr/bin/env python3
"""P0.1 input freeze: sha256 + size inventory of input_links, zip listings, package checksums."""
import hashlib, zipfile, os, csv, sys

BASE = r"D:\stroke_apa_reanalysis"
OUT = os.path.join(BASE, "results", "00_inventory")
os.makedirs(OUT, exist_ok=True)

def sha256(path, buf=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(buf)
            if not b:
                break
            h.update(b)
    return h.hexdigest()

rows = []
link_dir = os.path.join(BASE, "input_links")
for name in sorted(os.listdir(link_dir)):
    p = os.path.join(link_dir, name)
    if os.path.isfile(p):
        rows.append({"file": f"input_links/{name}", "bytes": os.path.getsize(p), "sha256": sha256(p)})

# governing docs
for name in ["01_研究Proposal.md", "02_项目TODO.md"]:
    p = os.path.join(BASE, name)
    rows.append({"file": name, "bytes": os.path.getsize(p), "sha256": sha256(p)})

with open(os.path.join(OUT, "input_files_sha256.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["file", "bytes", "sha256"], delimiter="\t")
    w.writeheader()
    w.writerows(rows)

# zip listings + checksums
zipdir = r"D:\stroke_apa"
for zn, listing in [("ATP2A2_IGV_CHECK_PACKAGE_20260920.zip", "inventory_igv_zip.txt"),
                    ("STROKE_APA_FINAL_DELIVERABLES_20260920.zip", "inventory_final_zip.txt")]:
    zp = os.path.join(zipdir, zn)
    with zipfile.ZipFile(zp) as z:
        names = z.namelist()
    with open(os.path.join(OUT, listing), "w", encoding="utf-8") as f:
        f.write("\n".join(names) + "\n")
    rows.append({"file": zn, "bytes": os.path.getsize(zp), "sha256": sha256(zp)})

with open(os.path.join(OUT, "package_sha256.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["file", "bytes", "sha256"], delimiter="\t")
    w.writeheader()
    w.writerows(rows[-2:])

print(f"input files hashed: {len(rows)-2}; zips listed+hashed: 2")
