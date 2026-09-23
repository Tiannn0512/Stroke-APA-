#!/usr/bin/env python3
# Extract paragraph text from a docx (zip -> word/document.xml), preserving order.
import re, sys, zipfile

path = sys.argv[1]
z = zipfile.ZipFile(path)
xml = z.read("word/document.xml").decode("utf-8")

# split into paragraphs; keep table cell boundaries readable
paras = re.split(r"</w:p>", xml)
for p in paras:
    texts = re.findall(r"<w:t[^>]*>([^<]*)</w:t>", p)
    line = "".join(texts).strip()
    if line:
        print(line)
