# provenance/: canonical generators that need inputs not packaged here

Each script is the canonical (last effective) version of the generator for a released output. The code is verbatim except that machine-specific absolute paths are masked. Every masked copy was checked to differ from its original **only** on lines containing a placeholder:

| Placeholder | Meaning |
|---|---|
| `<CTAS_WORKSPACE>` | the original analysis workspace root |
| `<CTAS_V2_PAPER_WORKSPACE>` | the original manuscript-asset workspace |
| `<REFERENCE_DIR>` | directory holding `reference/g1000_eur.*`, `reference/NCBI37.3.gene.loc` and `reference/hg38ToHg19.over.chain` |
| `<HOME>` | the original home directory |

In shell scripts the same placeholders appear as `${CTAS_WORKSPACE:?…}` (and so on), so that a script run without setting them stops with a message rather than misreading `<…>` as a redirection.

The "Re-executed" column shows what the authors re-ran on 16 Nov 2026 against the original workspace, with outputs written to a scratch directory and compared with the frozen originals.

| Folder | Generates | Re-executed (16 Nov 2026) | Inputs not packaged |
|---|---|---|---|
| `stage0_magma/` | 16 per-trait MAGMA gene tables (Pan-UKBB) | No. The MAGMA log flags match `run_magma_cvd.sh` (12 traits) and `process_4_existing_cvd_FIXED.sh` (4 traits). | Pan-UKBB downloads, 1000G EUR, NCBI37.3, MAGMA |
| `stage0_effective_n_magma/` | effective-N gene matrix used by the effective-N sensitivity analysis | No | as above |
| `stage1_locus_map/` | locked locus map | Runnable cleaned builder in `analysis/locus_map/` (byte-identical from minimal packaged inputs); originals archived here | — |
| `stage3_colocalisation/` | coloc.abf engine and re-summary (Supplementary Table 9, Fig. 4, Supplementary Fig. 4) | The re-summary and the engine are now runnable copies in `analysis/colocalisation/`. The engine was re-executed on 16 Nov 2026 with the surviving Aug 2026 slices plus refetched slices; see the release notes `COLOC_REFETCH_REPORT.md`. | Pan-UKBB `.bgz`, eQTL Catalogue slices |
| `stage3_finngen/` | FinnGen R12 GRCh38→GRCh37 conversion and MAGMA for 10 endpoints (`10`–`14`), analysis (`20`), comparison/ledger (`21`), six-family measure (`100`), Supplementary Table 11 (`133`) | **The corrected ten-endpoint run was made and validated on 4 Nov 2026** (pre-registered; 99.99% BIM agreement; 120/120 sentinels; QA passed). Its records are in `run_records_20261104/`, and its 180-file hash ledger was re-verified on 16 Nov 2026 (180/180). That run was not repeated during the 16 Nov packaging; everything downstream of its MAGMA outputs runs inside this repository (`analysis/finngen/`) with byte-identical results. | FinnGen R12 summary statistics (for steps 10–14 only) |
| `stage3_sensitivity/` | leave-one-trait-out, effective-N, obstructive panel | Runnable copies in `analysis/sensitivity/` (byte-identical from packaged inputs; obstructive also needs the public NCBI37.3 file) | NCBI37.3 (obstructive only) |
| `stage3_literature/` | four-comparator support | Runnable copies in `analysis/literature/` (byte-identical from the four publications' files; 58/91 reproduced) | comparator publication tables (not redistributed; see `inputs/literature/`) |
| `stage3_architecture/` | archived originals of `analysis/architecture/` and the leave-one-trait-out export | Yes | (runnable copy in `analysis/architecture/`) |
| `shared_functions/` | `fit_gaussian_copula()`, used by the sensitivity engines | n/a | n/a |
| `stage4_tables/` | Supplementary Tables 1–12, Tables 1–3 | **Yes**: all 12 supplementary workbooks are cell-identical to their staged builds | stage-3 intermediates |
| `stage4_figures/`, `stage4_figures_families/` | figure panels and their source data | **Yes**: 27/27 panel source tables byte-identical | stage-3 intermediates |
| `stage5_terminology_passes/` | label-only edits (I80/I82, I73) applied to the staged tables, figures and Explorer | No (presentation only). The staged → released differences equal the documented edits. | — |

The runnable analyses are in `analysis/` (see the top-level `README.md`).
