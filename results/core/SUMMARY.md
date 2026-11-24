# Task 1: primary estimator-C freeze (21 Oct 2026)

**Status: PASS. The one initial FAIL (minimum eigenvalue) is RESOLVED: Richard confirmed on 21 Oct 2026 that 0.50 is the paper-facing minimum eigenvalue for C (the expected 0.52 was Sigma_A's value).** The original FAIL row is kept in `validation_checks.tsv` for audit, with a PASS row `MIN_EIGEN_C_PAPER_FACING_0.50_RESOLVED` beside it. Nothing was forced. Results derive from the saved 10^9 runs (none regenerated); all 187 protected files are hash-identical before and after.

## Primary result (estimator C, M = 10^9, seed 42, strict `<`)
- Genes: **223 / 99 / 38** at q<0.05 / 0.01 / 0.001 (19,080 genes, unique). No gene sits exactly on a threshold.
- Loci: **91 / 46 / 21**; all map to the locked V2L universe (0 outside; 104 locked). Locus tiers: 45 q<0.05, 25 q<0.01, 21 q<0.001 (most-stringent member).
- Representative = highest-CTAS significant member; no exact ties; it coincides with the lowest-q gene in all 91 loci.

## Seed stability (C)
Gene and locus **membership sets** (not counts) are identical for seeds 42, 1, 2 at every tier: symmetric difference 0 everywhere (`seed_stability_C.tsv`; `seed_membership_differences_C.tsv` records none). A is also stable across seeds (symmetric difference 0).

## Eigenvalue / floor
- Min eigenvalue of Sigma_C = **0.5022736914630489** (paper-facing **0.50**); max 1.8934. The earlier expectation of ~0.52 does not apply to C (it is Sigma_A's value); resolved as above.
- Source of the mismatch: 0.52 is the minimum eigenvalue of the **raw-Spearman Sigma_A (0.5206325)**, not of Sigma_C (verified in `eigen_floor_check.tsv`). The 0.52 in CONTINUITY (Korbinian reply, 06_NEXT_ACTIONS) and the Bernard SOURCE_AND_QA note appears to quote A's value. If C is primary, the paper should say 0.50.
- Floor never activates (0 eigenvalues below 1e-10). Floor vs no-floor whiteners: max |diff| = 0 (CTAS values identical too). Whitener symmetric, W Sigma W = I to 1.1e-14. Sigma_A file matches recomputation to 5e-16.

## Raw-Spearman sensitivity (A, sensitivity only)
Genome-wide Spearman 0.99996171, Pearson 0.99996083. Genes A/C: 240/223, 107/99, 45/38 (C is a strict subset of A at every tier; C-only = 0); gene Jaccard 0.929 / 0.925 / 0.844; loci A/C 98/91, 50/46, 26/21; locus Jaccard 0.929 / 0.920 / 0.808. Estimator choice has little effect on genome-wide ranking; C is the prespecified primary analysis.

## Core retention
16/16 former q<0.001 loci and 32/32 genes remain at q<0.001 for seeds 42, 1 and 2 (`former_core_retention.tsv`).

## Notes / deviations
- Validation row `MIN_EIGEN_POSITIVE_APPROX_0.52` retained as FAIL for audit (observed 0.5023, expected 0.52); resolved by the user-approved 0.50 row; `INFO_MIN_EIGEN_SIGMA_A_REFERENCE` documents the likely source. No tolerance was widened to hide it.
- Whitener tolerance 1e-10 is justified because the floor cannot act unless an eigenvalue < 1e-10; observed difference is exactly 0.
- `zca_cor_whitener_C.csv` is the historical (floored) whitener; identical to the no-floor one.
- Extra files beyond the list: `session_info.txt`. `MANIFEST.tsv` cannot hold its own checksum.
- `common.R` and every locked source were not edited.

## Next gate
Gate resolved. Update the 0.52 mentions downstream (Korbinian reply, manuscript Methods, CONTINUITY) to 0.50 when those tasks run. Review this summary before FinnGen and colocalisation tasks begin.
