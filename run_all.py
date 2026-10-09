#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Stroke-APA pipeline stage list.

    python run_all.py           # list the stages
    python run_all.py --run     # with run notes
"""

from __future__ import annotations

import argparse
import sys

STAGES = [
    ("P0", "sample audit & QC",
     "scripts/v5_1/", "configs/frozen_P1_params.md",
     "results/v5_1/p0_cell_audit_raw.tsv, p1_qc_summary.tsv"),
    ("P1", "clustering, astrocyte calling, 3-layer DE, regulator catalog",
     "scripts/v5_1/", "configs/frozen_P1_params.md",
     "results/v5_1/p1_*.tsv, p1_astrocyte_clusters/"),
    ("P2", "DaPars2 APA discovery (SAP layer 2)",
     "scripts/v5_1/ + tools/DaPars2_patched", "configs/frozen_P2_SAP.md",
     "results/v5_1/p2_gate_g2_top_events.tsv, data/p2_sap_layer2_events.tsv"),
    ("P3", "regulator evidence: QKI eCLIP / motif / conservation",
     "scripts/v5_1/", "configs/frozen_evidence_class.md",
     "results/v5_1/p3_level*"),
    ("P4", "localization layer & candidate cards",
     "scripts/v5_1/", "configs/frozen_evidence_class.md",
     "results/v5_1/p4_*.tsv, p4_candidate_cards.md"),
    ("R1", "strand audit of all 9,693 event intervals",
     "scripts/reanalysis/", "configs/coordinate_convention.yaml",
     "results/reanalysis/01_strand_audit/"),
    ("R2", "candidate rebuild (QKI / motif / conservation / ranking v2)",
     "scripts/reanalysis/", "configs/coordinate_convention.yaml",
     "results/reanalysis/02_candidate_rebuild/"),
    ("R3", "BAM re-review (coverage ratios, read-level checks)",
     "scripts/reanalysis/", "configs/coordinate_convention.yaml",
     "results/reanalysis/03_bam_review/"),
    ("R4", "public PAP re-checks (GSE143531 edgeR TMM recompute)",
     "scripts/reanalysis/", "configs/frozen_evidence_class.md",
     "results/reanalysis/04_public_pap/"),
    ("R5", "target nomination & first assay design",
     "scripts/reanalysis/", "configs/coordinate_convention.yaml",
     "results/reanalysis/05_target_nomination/"),
    ("G1", "re-audit of the lead event: parsing fixes, structures, IGV evidence",
     "scripts/gate1_reaudit/", "configs/coordinate_convention.yaml",
     "results/06_gate1_reaudit/"),
    ("F1", "follow-up: PAS evidence, assay verification, cohort scan",
     "results/07_followup_20260925/analysis_scripts/",
     "configs/coordinate_convention.yaml",
     "results/07_followup_20260925/"),
]

WSL_STAGES = {"P0", "P2", "G1", "F1"}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true", help="also print run notes")
    args = ap.parse_args()

    print(f"Stroke-APA stages ({len(STAGES)})")
    print("=" * 78)
    for sid, name, script, frozen, outputs in STAGES:
        env = " [WSL]" if sid in WSL_STAGES else ""
        print(f"[{sid:>2}]{env} {name}")
        print(f"     scripts : {script}")
        print(f"     params  : {frozen}")
        print(f"     outputs : {outputs}")
    if args.run:
        print("-" * 78)
        print("Notes:")
        print("  * Heavy stages ran under WSL with the conda envs in envs/.") 
        print("  * Raw data is not included: see data/DATA_MANIFEST.md and")
        print("    docs/REGENERATE.md for sources and re-download commands.")
        print("  * Diff any re-run output against results/ before replacing.")
        print("  * Tests: python -m pytest tests/ -q")
    return 0


if __name__ == "__main__":
    sys.exit(main())
