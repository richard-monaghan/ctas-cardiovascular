# Colocalisation

## 1. Re-summary onto the final loci (runs from packaged inputs)
`01_coloc_resummary.py` reads the cell-level coloc.abf results (`inputs/colocalisation/V2_COLOC_LONG_RESULTS.tsv`: 104 loci × 7 tissues × 16 traits, signed 7 Aug 2026) and the final locus table. It writes the final cell, locus, tier, tissue and trait summaries behind Supplementary Table 9, Fig. 4 and Supplementary Fig. 4. All 12 analytic tables are byte-identical to the originals (`tests/reproduce_packaged_outputs.sh`).

## 2. coloc.abf engine (needs public external data)
`coloc_engine/v2_coloc_run.R` (with `v2_coloc_config.R`) runs coloc.abf (coloc 5.2.3, `beta_EUR` mode) for each canonical locus's statistical representative:
- 7 GTEx v8 cis-eQTL datasets from the eQTL Catalogue (study QTS000015: QTD000256 heart LV, QTD000251 heart AA, QTD000131 aorta, QTD000136 coronary, QTD000356 blood, QTD000266 liver, QTD000261 kidney cortex);
- all 16 Pan-UKBB traits;
- hg19 → hg38 via UCSC liftOver;
- MHC loci skipped;
- at least 50 matched SNPs per test.

Set these environment variables:

| Variable | Content |
|---|---|
| `CTAS_EQTL_SLICES_V1` | folder with `<gene>__<tissue>.tsv` eQTL slices for the 48 targets marked `cached=True` in `inputs/colocalisation/v2_target_coords.tsv` |
| `CTAS_EQTL_SLICES_V2` | the same for the other 56 targets |
| `CTAS_PANUKBB_BGZ_DIR` | local Pan-UKBB `<phenocode>-both_sexes.tsv.bgz` + `.tbi` for the 12 non-remote traits (the 4 others are read directly from the Pan-UKBB S3 bucket) |
| `CTAS_CHAIN_HG19_TO_HG38`, `CTAS_LIFTOVER_BIN` | UCSC chain file and liftOver binary |
| `CTAS_COLOC_OUT` | output folder (checkpoints per locus) |

**Slices.** A slice is `tabix https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/<QTD>/<QTD>.all.tsv.gz <chr>:<wstart>-<wend>`, using the window in `v2_target_coords.tsv` (gene body ± 1 Mb, GRCh38). The exact fetch used in the 16 Nov 2026 re-execution is recorded with URLs, access times and SHA-256 values; see the release notes `COLOC_REFETCH_REPORT.md`.

**Verified (16 Nov 2026):** all 104 per-locus checkpoints reproduced byte-for-byte (337 refetched + 391 original slices; checksums in `inputs/colocalisation/EQTL_SLICE_MANIFEST.tsv`; fetch script `fetch_eqtl_slices.py`).
