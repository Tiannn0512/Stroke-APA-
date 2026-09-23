#!/bin/bash
# Defensive patch: skip events with empty coverage (DaPars2 upstream assumes non-empty;
# zero-coverage UTRs in low-depth samples crash the division). Semantics unchanged:
# empty events become NA later anyway. Patched 2026-09-18, recorded in TODO P2-06.
F=/home/taylor/DaPars2/src/DaPars2_Multi_Sample_Multi_Chr.py
cp "$F" "$F.bak_20260918"
python3 - "$F" <<'PYEOF'
import sys
p = sys.argv[1]
t = open(p).read()
old = """    for i in range(num_samples):
        curr_Region_Coverage_raw = All_Samples_curr_3UTR_coverages[i]
        curr_Region_Coverage = curr_Region_Coverage_raw/weight_for_second_coverage[i]"""
new = """    for i in range(num_samples):
        curr_Region_Coverage_raw = All_Samples_curr_3UTR_coverages[i]
        if len(curr_Region_Coverage_raw) == 0:
            return 'Na', 'Na', 'Na'
        curr_Region_Coverage = curr_Region_Coverage_raw/weight_for_second_coverage[i]"""
assert old in t, "pattern not found"
t = t.replace(old, new, 1)
open(p, "w").write(t)
print("patched OK")
PYEOF
