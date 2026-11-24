# CTAS V2 colocalisation — sign-off package (R4)

**Date:** 2026-08-07. **Status:** COMPLETE — awaiting Richard's sign-off. F4/SF4 NOT built.
Deliverables: `Supplementary Tables/ST10_V2_colocalisation.xlsx` (11 sheets) + TSV mirrors
(`_source_assets/coloc_review/ST10_tsv_mirrors/`); locus abstraction + long results +
stratified summaries (`_source_assets/coloc_review/`); methods / validation records.

## What was run
Definitive design executed exactly: **104 canonical loci × 1 representative gene × 7 GTEx
cis-eQTL datasets × 16 Pan-UKBB CVD traits = 11,648 locus×tissue×trait cells**, one
consistent EUR beta/SE mode, Schmiedel Tfh_memory excluded, `coloc.abf`
(p1=p2=1e-4, p12=1e-5), hg19→hg38 liftOver, chr:pos:ref:alt harmonisation. Full method in
`V2_COLOC_METHODS.md`.

## Headline results
- **Locus classes (104):** 19 strong_PP4 (max PP4 ≥ 0.8), 14 moderate_PP4 (0.5–0.8),
  28 PP3_only (distinct-variant), 38 neither, 5 not_assessed.
- **Supported (PP4 ≥ 0.5): 33 / 104 loci** (99 evaluable; 5 non-evaluable).
- **Supported cells:** 164 (PP4 ≥ 0.5), of which 92 strong (PP4 ≥ 0.8).
- **Inherited vs new:** inherited loci 28/67 supported (17 strong, 11 moderate);
  new loci 5/37 supported (2 strong, 3 moderate) — colocalisation concentrates at
  inherited loci, as expected.
- **Tissues (supported cells):** aorta 41 (26 strong), heart atrial appendage 36,
  whole blood 26, coronary 21, heart LV 20, liver 17, kidney cortex fewest.
- **Traits (supported cells):** Hypertension 39 (27 strong), CAD 32, Angina 24, AF 21,
  MI 16, then the rest.
- **19 strong loci** include PHACTR1 (CAD, aorta, PP4≈1.00), CELSR2 (CAD, liver),
  TCF21 (MI, aorta), MRAS (CAD, coronary), ARHGAP42 (CAD, coronary), CEP68 (AF),
  FBN2 (HF), LPL (Angina), ORMDL3 (AF), FES, GTF2I, C6orf106 (HTN), plus two new-V2
  loci (SMC6, UPF3A) reaching strong.

## Non-evaluable (5, each with a stated reason — never "no colocalisation")
- **V2L036 / ATP6V1G2** — hg38 gene body (31.54 Mb) inside MHC chr6:28.5–33.5 Mb →
  `not_assessed_MHC`. (V2L037 / C6orf106, body 34.59 Mb, is outside → assessed, strong HTN hit.)
- **V2L019 HOXD13, V2L035 HDGFL1, V2L064 OR8H1, V2L065 KCNK4** — target gene has no
  cis-eQTL (`type==SNP`, maf>0.01) in any of the 7 datasets → not eQTL-evaluable.
These 5 are excluded from the evaluable denominator only.

## Validation (see `V2_COLOC_VALIDATION.md`)
- **Engine faithfulness:** rebuilt engine reproduces stored V1 exactly in legacy-input
  mode — PP4 correlation **1.00000**, max |ΔPP4| **0.0000** over 848 subset cells (both
  the beta and p-value paths).
- **Effect of the correction:** legacy→beta_EUR correlation 0.9975; only 3/848 cells
  change class; **at the locus level, V2 support agrees with V1 for all 48 unchanged
  targets (48/48, no gains, no losses)**. The corrections remove the V1 inconsistencies
  (mixed p-value/beta inputs; Tfh mislabel) without overturning any V1 conclusion.
- A list-indexing defect found mid-run (assigning NULL to a list element deleted it,
  misaligning tissues for loci with an absent slice) was fixed; the error-or-nothing
  property guaranteed no mislabeled checkpoint was ever persisted, and all 104 final
  checkpoints are clean (112 rows / 7 tissues). See `V2_COLOC_VALIDATION.md §5`.

## Provenance / integrity
- Canonical lock untouched (`v2_locus_master.csv` md5 unchanged).
- V1 coloc files (`run_coloc_all_101.R`, `results_all_101/`) untouched; V2 slices in a
  parallel `gtex_slices_v2/`; cached V1 slices reused read-only.
- Fetch QC: 392/392 new slices, 0 transport failures (`logs_v2/v2_slice_fetch_qc.tsv`).

## Requested sign-off
Please confirm before F4/SF4 are designed:
1. Locus-classification scheme (strong/moderate/PP3-only/neither/not_assessed) and the
   PP4≥0.5 / ≥0.8 thresholds.
2. Handling of the 5 non-evaluable loci (MHC + 4 no-eQTL).
3. The descriptive (non-inferential) framing of the V1↔V2 comparison.
Panel options for F4/SF4 are drafted in `R4_F4_SF4_PANEL_OPTIONS.md` (not built).
