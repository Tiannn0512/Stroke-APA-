# -*- coding: utf-8 -*-
"""Repair mixed-encoding memory file, then apply the 579->471 correction."""
import io

p = r"D:\stroke_apa\memory\2026-09-20.md"
raw = open(p, "rb").read()
try:
    txt = raw.decode("utf-8")
    print("file is clean utf-8")
except UnicodeDecodeError as e:
    idx = e.start
    print(f"bad byte at {idx}; part1 utf-8 ({idx}B) + part2 gbk")
    part1 = raw[:idx].decode("utf-8")
    part2 = raw[idx:].decode("gbk", errors="replace")
    txt = part1 + part2
    open(p, "wb").write(txt.encode("utf-8"))
    print("re-encoded to utf-8")

n_579 = txt.count("579")
txt = txt.replace("≥2 时间点 579", "≥2 时间点 471（2026-09-20 审计修正：原 579 系 evidence 脚本 and/or 优先级 bug 虚增，正确值 471）")
open(p, "wb").write(txt.encode("utf-8"))
print(f"579 occurrences before fix: {n_579}; fixed")
