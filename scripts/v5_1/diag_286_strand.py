"""Diagnostic: GSE286075 ReadsPerGene strand-column distribution (which column did P1-05 actually need?)."""
import gzip, os

d = r"D:\stroke_apa\data\GSE286075\extract"
for fn in sorted(os.listdir(d)):
    if not fn.endswith(".tab.gz"):
        continue
    with gzip.open(os.path.join(d, fn), "rt") as f:
        rows = [ln.rstrip("\n").split("\t") for ln in f]
    print("==", fn)
    for r in rows[:4]:
        print("   summary:", r[:2])
    genes = [r for r in rows if not r[0].startswith("N_")]
    c2 = sum(int(r[1]) for r in genes)  # col2 unstranded
    c3 = sum(int(r[2]) for r in genes)  # col3 1st-read-strand
    c4 = sum(int(r[3]) for r in genes)  # col4 2nd-read-strand
    print(f"   sums: col2(unstranded)={c2:,}  col3(1st)={c3:,}  col4(2nd)={c4:,}  col3+col4={c3+c4:,}")
