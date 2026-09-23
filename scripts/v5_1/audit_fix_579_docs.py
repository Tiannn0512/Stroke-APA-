# -*- coding: utf-8 -*-
"""Registered correction: >=2-timepoint concordant events 579 -> 471.
Root cause: operator-precedence bug in p2_gate_g2_evidence.py
(`A and B or C` counted single-contrast all-positive events).
"""
import io

FIX_NOTE = "（2026-09-20 审计修正：原 579 系 evidence 脚本 and/or 优先级 bug 虚增，正确值 471）"

# 1. TODO P2 完成记录
p = r"D:\stroke_apa\TODO_pipeline_v5.1.md"
t = io.open(p, encoding="utf-8").read()
before = t
t = t.replace("≥2 时间点 579", f"≥2 时间点 471{FIX_NOTE}")
if t != before:
    io.open(p, "w", encoding="utf-8", newline="").write(t)
    print("TODO fixed")

# 2. memory/2026-09-20.md
p = r"D:\stroke_apa\memory\2026-09-20.md"
t = io.open(p, encoding="utf-8").read()
before = t
t = t.replace("≥2 时间点 579", f"≥2 时间点 471{FIX_NOTE}")
if t != before:
    io.open(p, "w", encoding="utf-8", newline="").write(t)
    print("memory fixed")

# 3. inventory
p = r"D:\stroke_apa\results\data_inventory.tsv"
t = io.open(p, encoding="utf-8").read()
before = t
t = t.replace("≥2 时间点 579", f"≥2 时间点 471{FIX_NOTE}")
if t != before:
    io.open(p, "w", encoding="utf-8", newline="").write(t)
    print("inventory fixed")

# 4. G2 evidence extras: append correction (do not edit historical numbers in place)
p = r"D:\stroke_apa\results\p2_gate_g2_evidence_extras.tsv"
with open(p, "a", encoding="utf-8") as f:
    f.write(f"# CORRECTION 2026-09-20: union_events_significant_at_>=2_timepoints was 579 via "
            f"and/or precedence bug; correct value 471 (see audit_fix_579.py)\n")
print("extras annotated")

# 5. FINAL_REPORT: insert corrected count in §二 row 2
p = r"D:\stroke_apa\STROKE_APA_FINAL_REPORT.md"
t = io.open(p, encoding="utf-8").read()
old = "**977 显著**（联合 BH<0.1 & \\|ΔPDUI\\|≥0.1 & 双重复同向）；d1–d60 全程缩短主导"
new = "**977 显著**（联合 BH<0.1 & \\|ΔPDUI\\|≥0.1 & 双重复同向；**471 个在 ≥2 时间点同向显著**）；d1–d60 全程缩短主导"
if old in t:
    t = t.replace(old, new)
    io.open(p, "w", encoding="utf-8", newline="").write(t)
    print("report updated")
else:
    print("report pattern not found - check manually")
