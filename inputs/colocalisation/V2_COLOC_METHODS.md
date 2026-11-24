# CTAS V2 colocalisation — methods

**Date:** 2026-08-07. Authority scripts: `_source_assets/scripts/v2_coloc_config.R`
(single method source), `v2_coloc_run.R` (engine). This supersedes the V1 coloc
implementation (`eqtl_vignette/run_coloc_all_101.R`) for the V2 manuscript.

## Design
Definitive design: **104 canonical loci (V2L001–V2L104) × 1 statistical
representative gene per locus × 7 GTEx cis-eQTL datasets × 16 Pan-UKBB CVD traits**,
minus documented non-evaluable combinations and the MHC rule. Loci and target genes
are the frozen `panel_v2_locus_lock` objects; no locus is re-derived and no
105-locus object is used.

## eQTL data
Seven valid GTEx cis-eQTL datasets from the eQTL Catalogue (GRCh38):
heart left ventricle (QTD000256), heart atrial appendage (QTD000251), artery aorta
(QTD000131), artery coronary (QTD000136), whole blood (QTD000356), liver
(QTD000266), kidney cortex (QTD000261). Sample sizes N = 382, 372, 387, 213, 670,
208, 73 respectively.

The V1 "eqtlgen" slot (QTS000026/QTD000439 = Schmiedel_2018 Tfh_memory, N≈73,
mislabelled as whole blood) is **excluded everywhere** — inputs, denominators,
tables and figures. GTEx whole blood is the sole whole-blood dataset.

cis-window = representative gene body ± 1 Mb (hg38, Ensembl GRCh38 symbol lookup).
Slices for the 48 unchanged targets reuse the cached V1 slices unchanged; the 56
new/changed targets were fetched fresh to `eqtl_vignette/gtex_slices_v2/` by tabix
range query. eQTL rows are filtered to the target ENSG with `type == "SNP"` and
`maf > 0.01`.

## GWAS data — one consistent EUR effect-size mode
All 16 traits use a single, consistent European-ancestry effect-size
implementation from the Pan-UKBB flat files (GRCh37): `beta_EUR`, `se_EUR`,
`neglog10_pval_EUR` (QC), `low_confidence_EUR` (excluded), with EUR case/control
sample sizes from the effective-N authority and an EUR **weighted** allele
frequency `(N_cases·af_cases_EUR + N_controls·af_controls_EUR)/N`. This replaces
the two V1 inconsistencies: (i) `beta_meta`/`af_cases_meta` for 12 traits and
p-value-only input for CAD_I25/AF_I48/HF_I50/Hypertension_401, now all beta_EUR;
(ii) the `af_cases_meta` shortcut, now the documented EUR weighted frequency.
Column positions are resolved per file from each header by name (layouts differ,
e.g. 28 vs 52 columns) — never by fixed index. 12 traits are read from local
`.bgz`; 4 (icd10-I25, icd10-I48, icd10-I50, phecode-401) by remote tabix from the
Pan-UKBB S3 flat files.

## Harmonisation and colocalisation
GWAS variants are lifted hg19 → hg38 (`liftOver` + `hg19ToHg38.over.chain`);
pre/post-liftover and unmapped counts are retained per cell. GWAS and eQTL are
harmonised on `chr:pos:ref:alt` with allele match, strand-flip (β sign and MAF
inverted) and complement handling; the lowest-p duplicate per SNP is kept. Tests
with < 50 matched SNPs are not run.

Colocalisation uses `coloc.abf` (coloc 5.2.3): GWAS `type = "cc"` with case
fraction `s`; eQTL `type = "quant"`; priors `p1 = p2 = 1e-4`, `p12 = 1e-5`.
Posterior probabilities PP0–PP4 are recorded per locus × tissue × trait.
**Supported** colocalisation = PP4 ≥ 0.5; **strong** = PP4 ≥ 0.8; a **distinct-variant
(PP3-only)** locus has no PP4 ≥ 0.5 but at least one PP3 ≥ 0.5.

## MHC
Genes whose hg38 body lies within chr6:28,500,000–33,500,000 are flagged; standard
single-causal-variant coloc is not defensible there, so they are recorded
`not_assessed_MHC` and excluded from the evaluable denominator only (never called
"no colocalisation"). Of the two candidates, only **V2L036/ATP6V1G2** (body 31.54 Mb)
is inside; **V2L037/C6orf106** (body 34.59 Mb) is outside and assessed normally.

## Locus classification
Each locus is assigned the best evidence across all evaluated cells:
`strong_PP4` (max PP4 ≥ 0.8), `moderate_PP4` (0.5 ≤ max PP4 < 0.8), `PP3_only`
(no PP4 ≥ 0.5 but max PP3 ≥ 0.5), `neither` (assessed, nothing reached threshold),
`not_assessed` (MHC or no evaluable cell). No V1 tiers are used in the V2 summary.

## Reproducibility / validation
The rebuilt engine reproduces the stored V1 numbers exactly in a legacy-input mode
(PP4 correlation 1.00000, max |ΔPP4| = 0.0000 over 848 subset cells, both beta and
p-value paths). See `V2_COLOC_VALIDATION.md`. Per-locus checkpoints make the run
fully resumable.
