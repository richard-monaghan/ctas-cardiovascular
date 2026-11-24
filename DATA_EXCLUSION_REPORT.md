# Data exclusion report: CTAS v1.0.0 public repository (prepared 16 Nov 2026)

Evidence comes from `packaging/04_exclusion_scan.sh`, run on `repo/` after assembly. The raw output is in `packaging/exclusion_scan_output.txt`. A copy of this report is shipped as `repo/DATA_EXCLUSION_REPORT.md`.

## Scan results

The repository holds 145 files and 21.4 MB; the largest file is 3.2 MB (`results/core/final_C_gene_results_seed42.tsv`).

*Update, 17 Nov 2026 (v1.0.0 release candidate):* the rescan of the completed repository finds 267 files, 50.5 MB, largest file 4.1 MB (`inputs/literature/hgnc_complete_set_20260919.txt.gz`, a documented exception), with E1–E9 and E11–E12 unchanged. E10 and the FinnGen wording are updated below.

| ID | Excluded class | Result |
|---|---|---|
| E1 | Raw GWAS, genotype or index files and MAGMA intermediates (`.gz .bgz .tbi .bed .bim .fam .bam .vcf .magma .annot .genes.out .genes.raw .parquet .h5`) | **0 files** |
| E2 | Saved R objects, including the six 10⁹ null runs (`.rds .RData .rda`) | **0 files**. The runs are regenerable (README; the packaged engine is proven bit-identical at 10⁷). `.gitignore` blocks `work/` and `*.rds`. |
| E3 | Working caches and OS junk (`.DS_Store`, `__pycache__`, `.Rhistory`, `work/`, `cache/`) | **0** |
| E4 | Files > 10 MB | **0** |
| E5 | Individual-level identifier columns (`eid`, `IID`, `FID`) in tsv/csv/json | **0** |
| E6 | Superseded earlier-release result wording outside `provenance/` (the earlier release's headline gene and locus counts and tier names) | **0** |
| E7 | Superseded earlier-release result files (101-locus tables, `frozen_results`, 244/101 tables) | **0** |
| E8 | Absolute home-directory paths in text files | **0** (17 Nov 2026). The 2 Explorer build scripts that held them (`build_explorer_data.py`, `update_method_page.py`) now take explicit arguments or environment variables; the served site is unchanged. `provenance/` paths are masked. |
| E9 | Absolute paths inside the XLSX workbooks | **0** |
| E10 | E-mail addresses | **2**, both intentional: the corresponding contact address of R.M.M. in `README.md` and `CITATION.cff` (added 17 Nov 2026 at the author's instruction). No other addresses. |
| E11 | Credential keywords | 3 files (two Explorer data files and this report), all biology text in the Explorer (words such as secretory and secretion); this report also matches because it names them. **No credentials.** |
| E12 | Per-variant tables (> 100,000 rows) | **0**. The largest tables are gene-level (19,080 rows). |

## What *is* included and why it is allowed
- **Gene-level MAGMA P matrix (19,080 × 16).** A derived summary computed in-house from public Pan-UK Biobank summary statistics. It is not per-variant or individual-level data.
- **FinnGen-derived gene-level P values** (ST11 `gene_ledger`, 19,080 genes × 10 endpoints) and FinnGen support flags in the Explorer. These are derived from FinnGen R12 **public** summary statistics. No FinnGen variant-level rows are included. FinnGen's official acknowledgement sentence is used in `README.md` and `.zenodo.json`, with the Kurki *et al.* 2023 citation; FinnGen publishes no licence or redistribution clause for the public summary statistics (verified 16 Nov 2026).
- **Locked locus map.** The locus map keeps identifiers only. The earlier-release score columns and the superseded V2-lock scores (Σ_A, 10⁶) were removed (`packaging/02_derive_locus_map.R`). Legacy `ST2-n` identifiers are labels, not earlier-release results.
- **`provenance/shared_functions/00_utility_functions_REFERENCE_COPY.R`.** Code from the earlier release (it includes an ACAT helper), kept because the V2 sensitivity engines source `fit_gaussian_copula()` from it. It contains no results.

## Explicitly excluded from the source workspace
- **Raw downloads and reference panels.** Pan-UKBB raw downloads and MAGMA-formatted files, `magma_cvd_results/` per-trait outputs and annotations, FinnGen R12 summary statistics and lifted files, eQTL/GTEx data, the 1000G reference and the MAGMA binary.
- **Saved null objects.** `large_null_rerun_20261005/work/*.rds` (~1.1 GB including the 10⁷/10⁸ runs).
- **Earlier-release deposit.** `zenodo_deposit/` in full: the superseded copula/ACAT method, its 244/101 results, tables and figures.
- **Superseded intermediates.** The pre-final V2 lock tables (`gene_stringency_scores.csv`; the q and tier columns of `v2_gene_to_locus.csv` and `v2_locus_master.csv`), the superseded FinnGen endpoint-version assets and the GRCh38-mismatch FinnGen outputs.
- **Manuscript and editorial files.** Manuscript DOCX/PDF, PPTX decks, review and QA folders, CONTINUITY notes and memory files.
- **Workspace housekeeping.** Protected-hash ledgers containing absolute workspace paths.

Excluded scripts and the reason for each are in `SCRIPT_PROVENANCE.tsv` (decision `EXCLUDED`).
