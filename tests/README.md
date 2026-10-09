# tests — regression tests

These tests encode the two parsing bugs found by the 2026-09-24 audit; they
must stay green for any re-run of the read-level analyses.

| File | Covers |
|---|---|
| `test_g1_1.py` | CIGAR parsing (`cigar_op_bases`): splice-N / soft-clip / hard-clip / deletion accounting; junction extraction; **part C** is a real-BAM regression (needs a local BAM slice; skipped if absent) |
| `test_g1_3_seq.py` | FASTA extraction (`verify_window`): synthetic fixed-width FASTA with line offsets and cross-line windows, plus real-genome spot checks (needs `chr5.fa`; paths at the top of the file) |

Run:

```bash
python -m pytest tests/ -q          # or run each file directly: python tests/test_g1_1.py
```

The modules under test (`g1_1_read_metrics.py`, `g1_3_seq_utils.py`) are
copied here from the Gate-1 audit workspace so the tests are self-contained;
the originals live in `results/06_gate1_reaudit/` provenance and in the
release zip `v2026-09-24-gate1-reaudit`.
