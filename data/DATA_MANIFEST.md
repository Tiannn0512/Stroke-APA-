# DATA_MANIFEST — sources and checksums

No raw data is redistributed in this repository. This manifest lists every
external input, where to get it, and (for local references) the SHA-256 of
the copies used in the analyses. Regeneration order and commands:
`docs/REGENERATE.md`.

## Discovery cohort

| Item | Value |
|---|---|
| GEO accession | **GSE238125** (mouse cortex FACS-sorted astrocytes, MCAO; sham/d1/d3/d7/d21/d60, n=2 per time point; this project analyzes sham/d3/d7 with the shared sham) |
| FASTQ | via ENA (see `docs/REGENERATE.md` §4 for the exact re-download commands) |
| Sample map | `data/input_links/p2_sample_map.tsv` |
| Alignment | STAR 2.7.11b on mm10/GRCm38 + GENCODE vM25 (transcriptome-command documented in the reproduction manual) |

## Public re-check cohorts (gene-level)

GSE74456 (Sakers 2017 PAPTRAP), GSE143531 (Mazaré 2020), GSE225110
(Aldh1l1-RiboTag stroke), GSE330741 (AAV9-MPRA, 2026). Download URLs,
per-file excerpts, and verdicts: `results/07_followup_20260925/notes/gse_recheck_log.md`
(133 MB download cache not committed; re-downloadable from the URLs therein).

## Reference files (local copies used; SHA-256)

| File | SHA-256 |
|---|---|
| chr5.fa (UCSC mm10) | `1b970a617aec9b406ac823fad5a44538d1a9180e56d62ee121bf084ecc75ad5b` |
| gencode.vM25.annotation.gtf.gz | `e6c9704cb5832b8a1c71431ecac9ea30057377c72cdc9a12dc5d0405ab066d5a` |
| gencode.vM25.transcripts.fa | `145475bc09af4192cac6aa6757dcbcda24c657fba815b380f17b1b72780e0f16` |

## Ancillary resources

| Resource | Version / date | URL |
|---|---|---|
| PolyASite 2.0 atlas | GRCm38.96, downloaded 2026-09-25 | https://polyasite.unibas.ch/download/atlas/2.0/GRCm38.96/atlas.clusters.2.0.GRCm38.96.bed.gz |
| GENCODE vM25 | release M25 (GRCm38) | https://www.gencodegenes.org/mouse/ |
| QKI eCLIP (ENCODE) | per `data/input_links/p3_level2_qki_clip_peaks.bed` provenance | https://www.encodeproject.org/ |
| ATtRACT | backup 26.08.2020 (in release zips) | https://attract.cnr.it/ |

## Frozen input anchors (in `data/input_links/`)

DaPars2 event coordinates, PDUI matrix, SAP layer-2 events, second-tool
consistency, QKI peaks (BED), motif family mapping, frozen RBP motifs (MEME),
3′UTR BED for DaPars2, sample map, MPRA counts, PAP gene-set — these are the
hash-frozen inputs that the downstream stages consumed; see
`results/reanalysis/00_inventory/` for the full SHA-256 inventory.
