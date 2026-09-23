# -*- coding: utf-8 -*-
"""Repair GBK-contaminated files (bash heredoc appends), normalize to UTF-8,
apply 579->471 correction, scan repo text files for encoding issues."""
import glob, io, os

def smart_decode(path):
    raw = open(path, "rb").read()
    try:
        return raw.decode("utf-8"), False
    except UnicodeDecodeError as e:
        part1 = raw[:e.start].decode("utf-8", errors="ignore")
        part2 = raw[e.start:].decode("gbk", errors="replace")
        fixed = part1 + part2
        open(path, "wb").write(fixed.encode("utf-8"))
        return fixed, True

# 1. data_inventory.tsv
p = r"D:\stroke_apa\results\data_inventory.tsv"
t, was = smart_decode(p)
t2 = t.replace("≥2 时间点 579", "≥2 时间点 471（09-20 审计修正：原 579 系脚本 and/or 优先级 bug 虚增）")
io.open(p, "w", encoding="utf-8", newline="").write(t2)
print(f"inventory: encoding_repaired={was}; 579_fixed={t != t2}")

# 2. G2 evidence extras annotation
p = r"D:\stroke_apa\results\p2_gate_g2_evidence_extras.tsv"
t, was = smart_decode(p)
if "CORRECTION" not in t:
    t += ("\n# CORRECTION 2026-09-20: union_events_significant_at_>=2_timepoints was 579 via "
          "and/or precedence bug; correct value 471\n")
    io.open(p, "w", encoding="utf-8", newline="").write(t)
print(f"extras: encoding_repaired={was}; annotated=True")

# 3. scan all text files for encoding issues
bad = []
for pat in ["*.md", "*.tsv", "*.py", "*.sh", "*.html"]:
    for f in glob.glob(rf"D:\stroke_apa\{pat}") + glob.glob(rf"D:\stroke_apa\results\{pat}") + \
             glob.glob(rf"D:\stroke_apa\memory\{pat}") + glob.glob(rf"D:\stroke_apa\scripts\{pat}") + \
             glob.glob(rf"D:\stroke_apa\reviews\{pat}"):
        raw = open(f, "rb").read()
        try:
            raw.decode("utf-8")
        except UnicodeDecodeError:
            bad.append(f)
            fixed, _ = smart_decode(f)
print(f"encoding scan: {len(bad)} repaired -> {bad}")

# 4. FINAL_REPORT: add 471 to section 2
p = r"D:\stroke_apa\STROKE_APA_FINAL_REPORT.md"
t = io.open(p, encoding="utf-8").read()
needle = "双重复同向）；d1–d60 全程缩短主导"
add = "双重复同向；**471 个在 ≥2 时间点同向显著**）；d1–d60 全程缩短主导"
if "471 个" not in t and needle in t:
    t = t.replace(needle, add)
    io.open(p, "w", encoding="utf-8", newline="").write(t)
    print("report: 471 inserted")
else:
    print("report: 471 already present or pattern missing ->", "471 个" in t)
