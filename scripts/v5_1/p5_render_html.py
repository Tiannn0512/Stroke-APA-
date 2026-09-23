# -*- coding: utf-8 -*-
"""Render STROKE_APA_FINAL_REPORT.md -> STROKE_APA_FINAL_REPORT.html (tables+formatting)."""
import io, re, html as H

src = io.open("STROKE_APA_FINAL_REPORT.md", encoding="utf-8").read()
lines = src.split("\n")
out, in_table = [], False
for ln in lines:
    e = H.escape(ln)
    e = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", e)
    if ln.startswith("|"):
        cells = [c.strip() for c in e.strip("|").split("|")]
        if all(re.fullmatch(r"[-: ]+", c) for c in cells):
            continue
        if not in_table:
            out.append('<table border=1 cellspacing=0 cellpadding=4>')
            in_table = True
            out.append("<tr>" + "".join(f"<th>{c}</th>" for c in cells) + "</tr>")
        else:
            out.append("<tr>" + "".join(f"<td>{c}</td>" for c in cells) + "</tr>")
        continue
    if in_table:
        out.append("</table>")
        in_table = False
    if ln.startswith("# "):
        out.append(f"<h1>{e[2:]}</h1>")
    elif ln.startswith("## "):
        out.append(f"<h2>{e[3:]}</h2>")
    elif ln.startswith("> "):
        out.append(f"<blockquote>{e[2:]}</blockquote>")
    elif ln.startswith("- "):
        out.append(f"<li>{e[2:]}</li>")
    elif re.fullmatch(r"\d+\. .*", ln):
        out.append(f"<li>{e[3:]}</li>")
    elif ln.strip() == "---":
        out.append("<hr>")
    else:
        out.append(f"<p>{e}</p>")
doc = ('<!DOCTYPE html><html><head><meta charset="utf-8"><title>Stroke APA Final Report v5.1r</title>'
       '<style>body{font-family:Segoe UI,Arial,sans-serif;max-width:1000px;margin:2em auto;line-height:1.6;padding:0 1em}'
       'table{border-collapse:collapse;font-size:.95em}th{background:#eef}td,th{border:1px solid #aaa;padding:4px 8px}'
       'h1{border-bottom:3px solid #2a5d9f}h2{color:#2a5d9f;border-bottom:1px solid #ccc}'
       'blockquote{background:#f4f7fb;padding:.6em 1em;border-left:4px solid #2a5d9f}</style>'
       '</head><body>' + "\n".join(out) + "</body></html>")
io.open("STROKE_APA_FINAL_REPORT.html", "w", encoding="utf-8").write(doc)
print("HTML regenerated OK, bytes:", len(doc))
