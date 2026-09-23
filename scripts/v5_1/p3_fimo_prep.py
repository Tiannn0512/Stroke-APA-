"""P3-03 prep: convert ATtRACT pwm.txt (">Matrix_id width" + probability rows)
-> canonical MEME minimal format, keeping only our regulator panel motifs.
Outputs: data/attract/regulators.motifs.meme (FIMO-ready single file).
"""
import pandas as pd

PWM_TXT = r"D:\stroke_apa\data\attract\pwm.txt"
MAP_TSV = r"D:\stroke_apa\results\p3_level3_attract_regulator_motifs.tsv"
OUT = r"D:\stroke_apa\data\attract\regulators.motifs.meme"

keep = pd.read_csv(MAP_TSV, sep="\t")
id2gene = dict(zip(keep["Matrix_id"].astype(str), keep["Gene_name"]))
print("regulator Matrix_ids to keep:", len(id2gene))

out_blocks, kept = [], []
cur_id, cur_rows = None, []
def flush():
    global cur_id, cur_rows
    if cur_id and cur_rows and cur_id in id2gene:
        gene = id2gene[cur_id]
        cons = "".join("ACGT"[max(range(4), key=lambda i: r[i])] for r in cur_rows)
        out_blocks.append(
            f"MOTIF {gene}_{cur_id}\n"
            f"letter-probability matrix: alength= 4 w= {len(cur_rows)} nsites= 20 E= 0\n" +
            "\n".join("  ".join(f"{x:.4f}" for x in r) for r in cur_rows) +
            f"\nURL attract.cnic.es; consensus={cons}")
        kept.append((gene, cur_id, len(cur_rows), cons))
    cur_id, cur_rows = None, []

for ln in open(PWM_TXT, encoding="utf-8"):
    ln = ln.rstrip("\n")
    if ln.startswith(">"):
        flush()
        head = ln[1:].split()
        cur_id = head[0]
    elif ln.strip() and cur_id is not None:
        vals = [float(x) for x in ln.split()]
        if len(vals) == 4:
            cur_rows.append(vals)
flush()

with open(OUT, "w", newline="\n") as f:
    f.write("MEME version 4\nALPHABET= ACGT\nstrands: + -\n"
            "Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n")
    f.write("\n\n".join(out_blocks) + "\n")

print(f"kept {len(out_blocks)} motifs / {len(set(g for g,_,_,_ in kept))} regulators -> {OUT}")
for g, mid, w, cons in kept:
    if g in ("QKI", "PABPN1", "NUDT21"):
        print(f"  {g}_{mid} w={w} consensus={cons}")
