#!/usr/bin/env python3
"""Stream-merge the 20.9GB per-base element records into true elements (memory-safe)."""
c = s = e = None
n_in = n_out = 0
with open("/home/taylor/reference_v2/phastcons60way/elements_merged.bed", "w") as out:
    for ln in open("/home/taylor/reference_v2/phastcons60way/phastCons60way_elements.bed"):
        ch, a, b = ln.split("\t")[:3]
        a, b = int(a), int(b)
        n_in += 1
        if ch == c and a <= e:
            if b > e:
                e = b
        else:
            if c is not None:
                out.write(f"{c}\t{s}\t{e}\n")
                n_out += 1
            c, s, e = ch, a, b
    if c is not None:
        out.write(f"{c}\t{s}\t{e}\n")
        n_out += 1
print(f"in={n_in} merged={n_out}")
