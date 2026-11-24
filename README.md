# CTAS: Cross-Trait Aggregation Score for cardiovascular disease

Code, derived results and interactive Explorer for:

> **Dependence-aware aggregation maps shared genetic architecture across cardiovascular traits.**
> Richard M. Monaghan, Simon G. Williams, Hui Guo, Bernard D. Keavney. *bioRxiv preprint forthcoming.*

**Version 1.0.0**, released 2026-11-24. Archived at Zenodo: [10.5281/zenodo.24499456](https://doi.org/10.5281/zenodo.24499456).

- **Repository:** https://github.com/richard-monaghan/ctas-cardiovascular
- **Explorer:** https://richard-monaghan.github.io/ctas-cardiovascular/
- **Contact:** Richard M. Monaghan (corresponding author; ORCID: https://orcid.org/0000-0002-3323-1337). Questions about the repository may be submitted through GitHub Issues: https://github.com/richard-monaghan/ctas-cardiovascular/issues

## Scientific purpose

Shared genetic effects across related cardiovascular diseases can reveal biology that crosses diagnostic boundaries. However, correlation between traits (shared participants, shared risk factors, overlapping definitions) makes naive cross-trait recurrence hard to interpret, because the same evidence is counted repeatedly.

CTAS asks one calibrated question per gene: does it show more coordinated positive evidence across a correlated trait panel than the panel's dependence alone would produce? The method:

1. **Estimate trait dependence.** It is estimated from gene-level MAGMA P values for 16 Pan-UK Biobank cardiovascular traits, as a transformed rank-correlation matrix R = 2 sin(πρ_S/6) (minimum eigenvalue 0.50).
2. **Decorrelate.** Each gene's 16-trait z-score profile is decorrelated by ZCA-cor whitening.
3. **Score.** CTAS = M × log₂(1 + PR) is computed from the positive whitened evidence (magnitude M, participation ratio PR).
4. **Calibrate.** The score is calibrated against 10⁹ null draws whitened and scored identically, and controlled genome-wide by Benjamini–Hochberg FDR.

**Primary result:**
- **Genes:** 223 / 99 / 38 at q < 0.05 / 0.01 / 0.001.
- **Loci:** 91 / 46 / 21 at the same thresholds.
- **Stability:** membership is identical across null seeds 42, 1 and 2.

## What this repository reproduces (demonstrated, not assumed)

| Level | Analyses | How to check |
|---|---|---|
| **Reproduced exactly from the packaged inputs** | The central CTAS analysis: Σ matrices, observed scores, 10⁹-draw nulls, empirical P, BH q, the 223/99/38 genes and the 91/46/21 loci with their representatives, seed stability and the raw-Spearman sensitivity. **Locus map:** rebuilt from minimal discovery-stage inputs. **Behaviour simulations** (Supplementary Table 4). **Family architecture** (Supplementary Table 3 data). **Sensitivity analyses:** leave-one-trait-out (Supplementary Table 7) and effective-N (Supplementary Table 6), from packaged historical intermediates and a regenerated 10⁶ null. **FinnGen external support** downstream of the preserved 4 Nov 2026 MAGMA run (Supplementary Table 11, Fig. 5b, Table 3 column). **Colocalisation re-summary** (Supplementary Table 9, Fig. 4). **Supplementary Tables 1, 2 and 5.** **Supplementary Table 12 counts.** | `Rscript tests/smoke_test.R`; `tests/reproduce_packaged_outputs.sh` (~1 min; 62 intermediates byte-identical); `./run_core_pipeline.sh` (~10 min, 10⁹ nulls) |
| **Reproduced exactly with public external inputs** | **Obstructive-panel sensitivity** (Supplementary Table 8): needs the public NCBI37.3 gene-location file, which must match the SHA-256 in `inputs/EXTERNAL_INPUT_CHECKSUMS.tsv` (19,427 lines); any other copy is rejected with a FAIL. **Comparator/literature overlap** (Supplementary Table 10, Fig. 5a; 58/91): needs the four publications' supplementary files, which are not redistributed; see `inputs/literature/`. **Table 1** case counts and heritabilities: 112/112 fields match the Pan-UKBB phenotype manifest. | `CTAS_NCBI37_GENE_LOC=… tests/reproduce_packaged_outputs.sh`; `analysis/literature/run_literature.sh <sources>`; `analysis/table1/derive_table1_from_panukbb.py <manifest>` |
| **Colocalisation coloc.abf** | `analysis/colocalisation/coloc_engine/` re-executes the original engine. **Reproduced exactly with public external inputs:** the 16 Nov 2026 re-execution (refetched + original eQTL Catalogue slices; Pan-UKBB summary statistics) gave 104/104 per-locus results byte-identical to the original run. | `analysis/colocalisation/README_COLOC.md` |
| **Preserved and validated; not re-executed here** | Pan-UKBB MAGMA (stage 0a; gene-matrix assembly verified byte-identical). **FinnGen GRCh38→GRCh37 conversion and the ten MAGMA analyses:** run and validated on 4 Nov 2026; records and hashes in `provenance/stage3_finngen/run_records_20261104/`; outputs packaged. | `provenance/` |
| **Presentation only** | Terminology passes, PowerPoint assembly and rendering, Table 2 (conceptual), curated text, Explorer data build | — |

## Repository layout

| Path | Contents |
|---|---|
| `config/ctas_paths.R` | The only path definitions; everything is relative to the repository root |
| `R/` | Canonical CTAS functions (`common.R`) and the chunked null engine (`engine.R`) |
| `analysis/core/` | Core pipeline, steps 00–03 |
| `analysis/simulations/` | Behaviour simulations (Supplementary Table 4) |
| `analysis/architecture/` | Family-share architecture with the final whitener |
| `analysis/tables/` | Builder for Supplementary Tables 1, 2 and 5 from the frozen core results |
| `analysis/locus_map/`, `analysis/sensitivity/`, `analysis/finngen/`, `analysis/colocalisation/`, `analysis/literature/`, `analysis/table1/` | Runnable copies of the downstream generators (see the table above) |
| `tests/` | Smoke test, packaged-output reproduction, freeze and workbook comparisons |
| `inputs/` | Derived inputs: gene × trait matrix, frozen Σ_A, locus map. See `inputs/INPUT_MANIFEST.md`. |
| `results/core/` | Frozen primary results |
| `results/supplementary_tables/`, `results/main_tables/` | Supplementary Tables 1–12 and the Table 3 workbook, as released |
| `results/figure_source_data/` | The exact data plotted in every figure panel |
| `explorer/site/` | Static interactive Explorer, build `4042a25e88` |
| `provenance/` | All other generators, by stage (see `provenance/README.md`) |
| `env/` | Recorded R/Python environments and external tool versions |

## Input requirements

The runnable analyses need only `inputs/` and `results/core/`. Trait columns use internal keys (for example `PAD` = ICD-10 I73, `DVT` = I80, `VTE` = I82); the full mapping is in `inputs/INPUT_MANIFEST.md`.

**Software:**
- **R:** R ≥ 4.5 with `Matrix`, `digest` and `cluster` (recorded: R 4.5.1, Matrix 1.7-4, digest 0.6.37, cluster 2.1.8.2); the coloc engine also needs `data.table` 1.17.8 and `coloc` 5.2.3, plus `tabix`.
- **Python 3:** pandas, numpy and openpyxl (recorded: Python 3.13.2).
- **System:** a Unix-like OS (`mclapply`), with `shasum` or `sha256sum`.

## Analysis order

| Step | Entry point | Output | Reproduced from packaged inputs? |
|---|---|---|---|
| 0a | `provenance/stage0_magma/` | 16 per-trait MAGMA gene tables | No: needs Pan-UKBB downloads, the 1000G EUR reference, NCBI37.3 and MAGMA v1.10 |
| 0b | `analysis/core/00_assemble_gene_matrix.R <magma_dir>` | gene × trait matrix | Needs the step-0a outputs. Verified byte-identical to `inputs/`. |
| 1 | `analysis/core/01_sigma_matrices.R` | Σ_A check, Σ_C, eigenvalues | **Yes** |
| 2 | `analysis/core/02_run_large_null.R 1e9` | 10⁹-draw nulls | **Yes**: bit-identical to the original at 10⁷; ~10 min for 10⁹ |
| 3 | `analysis/core/03_freeze_primary.R` | `results/core/` tables | **Yes**: byte-identical from the original 10⁹ runs |
| 4a | `analysis/simulations/run_behaviour_simulations.R` | Supplementary Table 4 data | **Yes**: byte-identical |
| 4b | `analysis/architecture/20_architecture.R` | architecture source data | **Yes**: 8/8 primary outputs byte-identical. The comparison with superseded summaries is optional. |
| 4c | `analysis/locus_map/build_locus_map.R` | locus map from discovery-stage inputs | **Yes** (byte-identical) |
| 4d | `analysis/sensitivity/` (00 → run_all_loo → run_v2_effective_n → obstructive) | Supplementary Tables 6–8 source data | **Yes** (obstructive also needs NCBI37.3) |
| 4e | `analysis/finngen/` (20 → 21_ledger → 100 → st11/133) | Supplementary Table 11, Fig. 5b | **Yes**, from the preserved FinnGen MAGMA outputs |
| 4f | `analysis/colocalisation/` (coloc_engine → 01_coloc_resummary) | Supplementary Table 9, Fig. 4 | re-summary **yes**; engine needs public eQTL/GWAS data |
| 4g | `analysis/literature/run_literature.sh <sources>` | Supplementary Table 10, Fig. 5a | needs the publications' files |
| 4h | `analysis/table1/derive_table1_from_panukbb.py <manifest>` | Table 1 check | needs the Pan-UKBB manifest |
| 5 | `analysis/tables/41_st_primary.py`, `provenance/stage4_*`, `provenance/stage5_terminology_passes/` | tables and figures | Supplementary Tables 1, 2 and 5 **yes** (cell-identical); the rest need stage-4c intermediates. The final PDFs were rendered in PowerPoint. |

## Expected final counts and checksums

| Quantity | Expected |
|---|---|
| Genes tested | 19,080 |
| Genes at q < 0.05 / 0.01 / 0.001 | 223 / 99 / 38 |
| Loci at q < 0.05 / 0.01 / 0.001 | 91 / 46 / 21 |
| Σ_C minimum eigenvalue | 0.5022736914630489 (0.50) |
| Genes with a tail-extrapolated P (seed 42) | 1 (LPA) |
| `results/core/final_C_gene_results_seed42.tsv` SHA-256 | `12c7a0e0d8cb9dc747d598296e416337137eb9a13d9e16e319ffd5567af78446` |
| `results/core/final_C_locus_results_seed42.tsv` SHA-256 | `5dc73ff679b29198a3f1d609cc0302398bfa78eda1f559fa0bcee750a5e97802` |
| `inputs/UKBB_16trait_gene_matrix.tsv` SHA-256 | `dfbbcc240350c656f4917adfaf23256892425706eb27bfdea83bb595e44b87b6` |

`CHECKSUMS.sha256` lists every file.

## Regenerating the 10⁹-draw null

The six saved null objects (~87 MB each) are not distributed. `02_run_large_null.R` regenerates them using L'Ecuyer-CMRG streams (seeds 42, 1 and 2; one stream per 2×10⁶-draw chunk), so the results do not depend on the number of cores. Set `CTAS_CORES` to control parallelism.

Bit-identity was demonstrated on the release machine (R 4.5.1, macOS arm64) at M = 10⁷. Other platforms or BLAS libraries may differ in the last floating-point digits. The gene and locus sets are identical across three seeds at 10⁹.

## Explorer

`explorer/site/` is a self-contained static site:

```bash
cd explorer/site && python3 -m http.server 8000
```

**Deployment:** GitHub Pages at https://richard-monaghan.github.io/ctas-cardiovascular/. The workflow `.github/workflows/pages.yml` publishes `explorer/site/` unchanged (build `4042a25e88`); it does not rebuild or rewrite the site. All internal links and data paths are relative; the site was tested locally under the `/ctas-cardiovascular/` path prefix (17 Nov 2026). `explorer/site/scripts/` holds the original data-build scripts. They were not re-executed for this release. Their hard-coded workspace paths were replaced by explicit arguments: `update_method_page.py [method.html]`, and `CTAS_WORKSPACE` / `CTAS_V2_PAPER_WORKSPACE` for `build_explorer_data.py`. The original workspaces are not part of this repository. The served site files are unchanged.

## Data availability and licences

- Pan-UK Biobank summary statistics: https://pan.ukbb.broadinstitute.org
- FinnGen R12 public summary statistics: https://www.finngen.fi/en/access_results. "We want to acknowledge the participants and investigators of the FinnGen study". Cite Kurki *et al.*, *Nature* 613, 508–518 (2023); see `inputs/finngen_r12/README_FINNGEN.md`.
- Comparator studies (Qiao 2025, Liu 2023, Foley 2021, PLEIO 2021): source tables are not redistributed; see `inputs/literature/COMPARATOR_SOURCE_MANIFEST.tsv`.
- HGNC gene nomenclature (CC0): https://www.genenames.org
- No individual-level data are included.

### Licences

| Material | Licence | File |
|---|---|---|
| Software and source code | MIT, Copyright (c) 2026 The University of Manchester | `LICENSE` |
| Project-generated derived results and documentation | CC BY 4.0 | `LICENSE-CC-BY-4.0.txt` |
| Third-party inputs and data derived from them | Their original terms; not relicensed | see `LICENSE_SCOPE.md` |

`LICENSE_SCOPE.md` explains which material falls under each licence and lists the upstream terms that continue to apply.

## Citation

Please cite the article (bioRxiv preprint forthcoming) and this archive (Zenodo, [10.5281/zenodo.24499456](https://doi.org/10.5281/zenodo.24499456)). See `CITATION.cff`.

## Funding

This work was supported by British Heart Foundation Programme Grant RG/F/21/110050. B.D.K. was supported by a British Heart Foundation personal chair. H.G. was supported by the Engineering and Physical Sciences Research Council (EPSRC) through the Causality in Healthcare AI Hub (CHAI), grant number EP/Y028856/1.

## Acknowledgements

We thank Korbinian Strimmer for advice on the whitening approach. We want to acknowledge the participants and investigators of the FinnGen study. We acknowledge the Pan-UK Biobank (Pan-UKB) project and UK Biobank. The Pan-UKB summary statistics were generated using the UK Biobank Resource (project ID 31063).
