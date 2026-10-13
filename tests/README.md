# Tests

Tests for the two read-level parsers. They were written after the 2026-09-24
audit found bugs in CIGAR parsing and FASTA extraction, so they are kept to
catch regressions on re-runs.

| File | What it tests |
|---|---|
| `test_g1_1.py` | CIGAR parsing: splice N, soft/hard clip, deletion. Part C checks junction counts against a real BAM slice (skipped if no local BAM) |
| `test_g1_3_seq.py` | FASTA extraction: fixed-width FASTA with line offsets and cross-line windows, plus real-genome spot checks (needs `chr5.fa`, path at the top of the file) |

```bash
python -m pytest tests/ -q
```

The modules under test (`g1_1_read_metrics.py`, `g1_3_seq_utils.py`) are
copied here so the tests run standalone.
