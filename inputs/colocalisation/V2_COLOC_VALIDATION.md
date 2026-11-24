# CTAS V2 colocalisation — validation record

**Date:** 2026-08-07. **Status:** validation bridge PASSED; full run in progress.
Scripts: `_source_assets/scripts/v2_coloc_{config,run,validation_bridge}.R`.
Bridge outputs: `_source_assets/coloc_review/validation/`.

## 1. Engine faithfulness — legacy reproduction of stored V1 (Step A)

To prove the rebuilt liftover / allele-harmonisation / `coloc.abf` plumbing is
faithful before trusting any corrected number, the V1 method was re-implemented
exactly — 12 traits via Pan-UKBB `beta_meta`/`se_meta`/`af_cases_meta`, and the 4
traits CAD_I25 / AF_I48 / HF_I50 / Hypertension_401 via the MAGMA p-value path —
and run on a deterministic 8-locus subset of the 48 unchanged targets, then
compared cell-by-cell to the stored V1 results (`results_all_101/coloc_pp4_all_101.tsv`).

**Deterministic subset (spans every required category):**

| Locus | Gene | V1 category | rationale |
|-------|------|-------------|-----------|
| V2L034 | PHACTR1 | strong_PP4 | canonical strong (CAD/aorta) |
| V2L017 | CEP68 | strong_PP4 | non-coronary strong (AF) |
| V2L004 | CELSR2 | strong_PP4 | lipid locus |
| V2L033 | FBN2 | strong_PP4 | non-coronary strong (HF) |
| V2L069 | VPS11 | PP3_only | distinct-variant example |
| V2L040 | CENPW | PP3_only | distinct-variant example |
| V2L095 | LDLR | negative | lipid, negative |
| V2L055 | ZFPM2 | negative | non-coronary negative (DVT) |

**Result — exact reproduction:**
- 848 locus × tissue × trait cells compared.
- PP4 Pearson correlation legacy vs stored V1 = **1.00000**.
- max |ΔPP4| = **0.0000**; mean = 0.0000; median = 0.0000.
- 848 / 848 cells (100.0%) within |ΔPP4| < 0.01.
- Exact in **both** GWAS paths: `bgz_beta` (636 cells) and `magma_pval` (212 cells),
  each max |ΔPP4| = 0.0000.

The rebuilt engine reproduces V1 bit-for-bit. Any difference in the corrected run
is therefore a genuine consequence of the design correction, not an implementation change.

## 2. Effect of the correction — legacy vs corrected beta_EUR (Step B)

Same 8 loci, corrected mode = one consistent EUR effect-size implementation for all
16 traits (`beta_EUR`/`se_EUR`, EUR weighted allele frequency, EUR case/control N),
Schmiedel Tfh_memory excluded.

- 848 cells in both; PP4 correlation legacy vs corrected = **0.9975**.
- mean |ΔPP4| = 0.0035; max |ΔPP4| = 0.157.
- Class transitions (strong ≥0.8 / moderate 0.5–0.8 / below): only **3 / 848** cells
  change class — 2 below→moderate (gained), 1 strong→moderate (FBN2·heart_AA·AF_I48,
  0.888→0.799). 815 below→below, 25 strong→strong, 5 moderate→moderate.
- The largest single shifts (|ΔPP4| ≈ 0.06–0.16) are concentrated on the four traits
  that V1 had assessed by p-values (CAD/AF/HF/HTN), now assessed by proper effect
  sizes — small and interpretable, in both directions.

**Interpretation.** Switching to the consistent beta_EUR design leaves colocalisation
conclusions essentially unchanged on the validated subset, while removing the two V1
inconsistencies (mixed p-value/beta inputs; Tfh mislabel). Strong CAD colocalisations
(PHACTR1, CELSR2) and the strong AF colocalisation at CEP68 (0.997) are preserved; the
PP3-only loci (VPS11, CENPW) and negatives (LDLR, ZFPM2) remain sub-threshold.

## 3. Non-evaluable handling
- MHC gene-body check (hg38): only **V2L036/ATP6V1G2** (body 31.54 Mb) lies inside
  chr6:28.5–33.5 Mb → `not_assessed_MHC`. **V2L037/C6orf106** (body 34.59 Mb) is outside
  → assessed normally.
- `no_ensg_in_window` (target gene absent from a tissue's cis eQTLs) and
  `too_few_matched_snps` (<50 harmonised SNPs) are recorded as non-evaluable for that
  cell — never counted as "no colocalisation".

## 5. Defect found and fixed during the run (transparency)
While running the full 104 loci, three cached loci (HOXD13, TBX20, APOA5) crashed with
"subscript out of bounds". Root cause: the eQTL slices were pre-allocated as a length-7
list and populated with `eqtl[[ti]] <- s$dt`; when a tissue slice was absent (`s$dt`
NULL), `[[<-`NULL **deletes** the list element in R, shrinking the list and misaligning
tissue indices for any locus with an absent slice. Fix: only assign when non-NULL, leaving
the pre-allocated NULL in place (`v2_coloc_run.R`, slice-loading loop).

This defect could only ever crash (the trait loop always reaches index 7, which is
out of bounds on a shrunken list), so **no mislabeled checkpoint was ever persisted** —
affected loci wrote nothing and were re-run after the fix. All 65 checkpoints written
before the fix were independently verified to have 112 rows and 7 distinct tissues, and
all 104 final checkpoints are clean. The validation-bridge subset (8 all-7-slice loci)
was unaffected, so §1–§2 stand unchanged.

## 4. Verdict
**PASS.** Legacy reproduction is exact; the corrected design is stable and its
departures are confined to the intended corrections. Cleared to run the definitive
104-locus × 7-dataset × 16-trait beta_EUR analysis.
