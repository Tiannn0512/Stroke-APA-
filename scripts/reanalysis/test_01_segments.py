#!/usr/bin/env python3
"""P1.3 unit tests for 01_build_strand_aware_segments.py.

Test set (TODO P1.3): one synthetic plus-strand event, one synthetic minus-strand
event, plus the real Atp2a2 / Ndrg2 / Agpat3 rows from the generated coordinate
table. Atp2a2 must reproduce the frozen anchors exactly.

Exit code 0 = all tests pass; 1 = any failure.
"""
import csv, importlib.util, os, sys

BASE = r"D:\stroke_apa_reanalysis"
spec = importlib.util.spec_from_file_location(
    "builder", os.path.join(BASE, "scripts", "01_build_strand_aware_segments.py"))
builder = importlib.util.module_from_spec(spec)
spec.loader.exec_module(builder)

rules = builder.load_rules(os.path.join(BASE, "config", "coordinate_convention.yaml"))
failures = []

def check(label, cond, detail=""):
    print(f"  [{'PASS' if cond else 'FAIL'}] {label}" + (f"  {detail}" if detail and not cond else ""))
    if not cond:
        failures.append(label)

print("== synthetic plus-strand event ==")
ev = {"strand": "+", "utr_start": "1000", "utr_end": "2000", "proximal_pas_bed": "1500"}
seg = builder.build_segments(ev, rules)
check("shared = [1000,1500)", seg["shared"] == (1000, 1500), str(seg["shared"]))
check("distal = [1500,2000)", seg["distal"] == (1500, 2000), str(seg["distal"]))

print("== synthetic minus-strand event ==")
ev = {"strand": "-", "utr_start": "1000", "utr_end": "2000", "proximal_pas_bed": "1500"}
seg = builder.build_segments(ev, rules)
check("distal = [1000,1500)", seg["distal"] == (1000, 1500), str(seg["distal"]))
check("shared = [1500,2000)", seg["shared"] == (1500, 2000), str(seg["shared"]))

print("== edge cases (must be rejected) ==")
for label, ev in [
    ("pas == utr_start (minus)", {"strand": "-", "utr_start": "1000", "utr_end": "2000", "proximal_pas_bed": "1000"}),
    ("pas == utr_end (plus)", {"strand": "+", "utr_start": "1000", "utr_end": "2000", "proximal_pas_bed": "2000"}),
    ("pas outside UTR", {"strand": "+", "utr_start": "1000", "utr_end": "2000", "proximal_pas_bed": "3000"}),
    ("bad strand", {"strand": "*", "utr_start": "1000", "utr_end": "2000", "proximal_pas_bed": "1500"}),
]:
    try:
        builder.build_segments(ev, rules)
        check(f"reject: {label}", False)
    except ValueError:
        check(f"reject: {label}", True)

print("== real events from coordinate table ==")
real = {}
with open(os.path.join(BASE, "input_links", "dapars2_event_coordinates.tsv"), encoding="utf-8") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["gene_symbol"] in ("Atp2a2", "Ndrg2", "Agpat3"):
            real[r["gene_symbol"]] = r
check("all three candidates present", set(real) == {"Atp2a2", "Ndrg2", "Agpat3"}, str(list(real)))

a = real["Atp2a2"]
seg = builder.build_segments(a, rules)
check("Atp2a2 anchors: distal = 122453512-122456498",
      seg["distal"] == (122453512, 122456498), str(seg["distal"]))
check("Atp2a2 anchors: shared = 122456498-122457153",
      seg["shared"] == (122456498, 122457153), str(seg["shared"]))
check("Atp2a2 unit-test QKI peak inside distal (122454200-122454250)",
      seg["distal"][0] <= 122454200 and 122454250 <= seg["distal"][1])
check("Atp2a2 old-cited peak on shared side (122456550-650)",
      seg["shared"][0] <= 122456550 and 122456650 <= seg["shared"][1])

for g in ("Ndrg2", "Agpat3"):
    r = real[g]
    seg = builder.build_segments(r, rules)
    check(f"{g}: strand minus, distal on low-coordinate side",
          seg["distal"] == (int(r["utr_start"]), int(r["proximal_pas_bed"])))
    check(f"{g}: lengths sum to UTR length",
          (seg["distal"][1] - seg["distal"][0]) + (seg["shared"][1] - seg["shared"][0])
          == int(r["utr_end"]) - int(r["utr_start"]))

print()
if failures:
    print(f"FAILED: {len(failures)} -> {failures}")
    sys.exit(1)
print("ALL TESTS PASS")
