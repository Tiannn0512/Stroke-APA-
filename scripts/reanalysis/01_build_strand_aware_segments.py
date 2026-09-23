#!/usr/bin/env python3
"""01_build_strand_aware_segments.py  (TODO P1.2)

Strand-aware split of each DaPars2 3'UTR event into shared (proximal-PAS-constitutive)
and distal-specific segments, per the frozen convention in config/coordinate_convention.yaml:

    plus_distal_rule:  [proximal_pas, utr_end)
    minus_distal_rule: [utr_start, proximal_pas)
    shared = the remaining (termination-codon-proximal) side.

Coordinates: BED 0-based half-open. proximal_pas_raw and proximal_pas_bed are both
carried through (P1.1 determination: identical; DaPars2 outputs BED-scale genomic
coordinates, no 1bp conversion).

Outputs (results/01_strand_audit/):
    strand_aware_event_segments.tsv
    distal_segments.strand_aware.bed
    shared_segments.strand_aware.bed
    invalid_event_coordinates.tsv
"""
import argparse, csv, os, sys

def load_rules(config_path):
    rules = {}
    with open(config_path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if ":" in line and not line.startswith("#"):
                key, _, val = line.partition(":")
                val = val.strip().strip('"').strip("'")
                if val.startswith("[") or val in ("+", "-") or key in (
                    "genome", "annotation", "bed_coordinates", "event_coordinate_source"):
                    rules[key.strip()] = val
    return rules

def build_segments(event, rules):
    """Return dict of segments or raise ValueError with reason."""
    strand = event["strand"]
    us, ue = int(event["utr_start"]), int(event["utr_end"])
    pas = int(event["proximal_pas_bed"])
    if strand not in ("+", "-"):
        raise ValueError(f"strand not +/-: {strand!r}")
    if not (us < pas < ue):
        raise ValueError(f"proximal PAS not strictly inside UTR: {us} <= {pas} <= {ue} fails")
    if rules.get("minus_distal_rule") == "[utr_start, proximal_pas)":
        if strand == "-":
            distal = (us, pas)
            shared = (pas, ue)
        else:
            shared = (us, pas)
            distal = (pas, ue)
    else:
        raise ValueError("config minus_distal_rule not the frozen convention")
    if not (0 < distal[1] - distal[0] and 0 < shared[1] - shared[0]):
        raise ValueError("non-positive segment length")
    if (distal[1] - distal[0]) + (shared[1] - shared[0]) != ue - us:
        raise ValueError("segment lengths do not sum to UTR length")
    lo, hi = sorted((distal, shared))
    if lo[1] > hi[0]:  # overlap check (segments are adjacent by construction)
        raise ValueError("segments overlap")
    return {"shared": shared, "distal": distal}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", required=True)
    ap.add_argument("--coordinate-config", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    rules = load_rules(args.coordinate_config)
    out = args.out_dir
    os.makedirs(out, exist_ok=True)

    seg_rows, invalid = [], []
    distal_bed, shared_bed = [], []
    seen = {}
    dups = 0
    with open(args.events, encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            eid = r["event_id"]
            if eid in seen:
                dups += 1
                invalid.append({**r, "reason": f"duplicate event_id (first at row {seen[eid]})"})
                continue
            seen[eid] = len(seen) + 1
            try:
                seg = build_segments(r, rules)
            except ValueError as e:
                invalid.append({**r, "reason": str(e)})
                continue
            name = eid  # tx|gene|symbol already unique; symbol included
            distal_bed.append((r["chrom"], seg["distal"][0], seg["distal"][1], name, 0, r["strand"]))
            shared_bed.append((r["chrom"], seg["shared"][0], seg["shared"][1], name, 0, r["strand"]))
            seg_rows.append({
                "event_id": eid, "event_num": r.get("event_num", ""), "chrom": r["chrom"], "strand": r["strand"],
                "gene_symbol": r["gene_symbol"],
                "utr_start": r["utr_start"], "utr_end": r["utr_end"],
                "proximal_pas_raw": r["proximal_pas_raw"], "proximal_pas_bed": r["proximal_pas_bed"],
                "shared_start": seg["shared"][0], "shared_end": seg["shared"][1],
                "distal_start": seg["distal"][0], "distal_end": seg["distal"][1],
                "shared_len": seg["shared"][1] - seg["shared"][0],
                "distal_len": seg["distal"][1] - seg["distal"][0],
            })

    fields = list(seg_rows[0].keys()) if seg_rows else ["event_id"]
    with open(os.path.join(out, "strand_aware_event_segments.tsv"), "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader(); w.writerows(seg_rows)

    for base, rows in (("distal_segments.strand_aware.bed", distal_bed),
                       ("shared_segments.strand_aware.bed", shared_bed)):
        with open(os.path.join(out, base), "w", encoding="utf-8", newline="") as f:
            for r in rows:
                f.write("\t".join(str(x) for x in r) + "\n")

    with open(os.path.join(out, "invalid_event_coordinates.tsv"), "w", encoding="utf-8", newline="") as f:
        if invalid:
            w = csv.DictWriter(f, fieldnames=list(invalid[0].keys()), delimiter="\t")
            w.writeheader(); w.writerows(invalid)
        else:
            f.write("reason\tempty - no invalid events\n")

    print(f"segments written : {len(seg_rows)}")
    print(f"invalid events   : {len(invalid)} (see invalid_event_coordinates.tsv)")
    print(f"duplicate ids    : {dups}")
    plus = sum(1 for r in seg_rows if r["strand"] == "+")
    print(f"strand counts    : + {plus} / - {len(seg_rows)-plus}")

if __name__ == "__main__":
    main()
