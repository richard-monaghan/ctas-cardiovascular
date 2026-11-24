# RERUN_PRESPECIFICATION — corrected FinnGen external-support rerun (frozen 2026-11-04, before any corrected outcome is read)

Authority for approval: Richard's instruction in chat, 2026-11-04 ("Approve the corrected FinnGen rerun ... with the safeguards below"). Motivation: `target_pipeline/apoe_finngen_forensic_01/` (GRCh38 SNP positions annotated against GRCh37 gene locations in the accepted run).
Nothing outside this folder is written. The original FinnGen files, outputs, the APOE audit, final CTAS results, manuscript, figures, tables and Explorer are read-only. No old gene-level MAGMA result is reused as an input.

## 1. Fixed inputs
- Endpoints (FinnGen R12; label = file code): AF=I9_AF, Angina=I9_ANGINA, CAD=I9_CORATHER, Cardiomyopathy=I9_CARDMYO, HF=I9_HEARTFAIL, Hypertension=I9_HYPTENS, MI=I9_MI_STRICT, Stroke=I9_STR, TIA=I9_TIA, VTE=I9_VTE.
- Raw files: `<CTAS_WORKSPACE>/CTAS_FinnGen/data/finngen_gwas/finngen_R12_<code>.gz` (GRCh38; columns `#chrom pos ref alt rsids nearest_genes pval mlogp beta sebeta af_alt af_alt_cases af_alt_controls`).
- Gene locations: `<REFERENCE_DIR>/reference/NCBI37.3.gene.loc` (GRCh37). LD reference: `<REFERENCE_DIR>/reference/g1000_eur` (1000G Phase 3 EUR, GRCh37; 22,665,064 SNPs, 503 individuals) with `g1000_eur.synonyms` auto-detected.
- Chain: `<REFERENCE_DIR>/reference/hg38ToHg19.over.chain`.
- MAGMA: `<CTAS_WORKSPACE>/MAGMA_Mac/magma`, "MAGMA v1.10 (custom)" (x86_64, run under Rosetta), default SNPwise-mean model.
- CTAS sets: `outputs/final_c_primary_freeze_20261021/final_C_gene_results_seed42.tsv`, `final_C_locus_results_seed42.tsv` (223/99/38 genes; 91/46/21 loci), unchanged.
- Frozen parameters carried over from the accepted run (from its logs): `--annotate window=0,0`; `--pval ... N=500348` for all ten endpoints; support = gene-level MAGMA P<0.01 in >=3 of 10 endpoints; genes absent from FinnGen = not evaluable (never zero).

## 2. Only intended change
SNP coordinates are lifted GRCh38 -> GRCh37 before annotation. Everything else is kept. Unavoidable/minor differences are logged in DEVIATIONS.tsv: (a) the MAGMA SNP-location/P input is written without a header line and read with `use=1,4` (the old header line only caused a harmless "CHR not recognised" skip); (b) tab-explicit parsing replaces whitespace splitting, so rows with an empty rsids field no longer shift columns.

## 3. Conversion rules
- Parse strictly by tab; every data row must have exactly 13 fields, header must equal the expected 13 names; any violation = STOP.
- Keep rows with non-empty rsids (the previous filter `rsid != "." and != ""`), pval numeric with 0<P<=1, P floored at 1e-300 (as before). Multi-rsID strings (`rs1,rs2`) are kept verbatim as the SNP ID (as before; MAGMA cannot match them to the reference, so they are inert in the gene test but appear in annotation).
- Lift each retained row with the chain: same-chromosome chains only; where blocks overlap choose the highest-score chain; unmapped, different-chromosome or indel-boundary positions are dropped and counted (failed lift). rsID, REF, ALT, P, original GRCh38 position and new GRCh37 position are preserved in the converted layer (`input_lifted/<EP>.lifted.tsv.gz`). MAGMA input per endpoint: `SNP CHR BP37 P`.
- Chromosome codes are kept as in the raw file (23 = X).

## 4. Gates (any failure = STOP; MAGMA gene analysis is not run; report instead)
G1 (per endpoint) Coordinate validation vs the exact 1000G BIM: among rsIDs that occur exactly once in the converted layer and exactly once in the BIM, share with identical (chromosome, BP) must be >= 95.0%. (Threshold allows for dbSNP rsID merges/recycling between the BIM build and R12; the unlifted agreement is reported alongside.) Also report successful/failed lift, duplicated rsIDs, multimapped, reference-allele mismatches (REF/ALT vs BIM alleles among position-matched rsIDs).
G2 Sentinels, all ten endpoints: rs429358 -> chr19:45,411,941 and rs7412 -> chr19:45,412,079 in the converted layer and in the MAGMA input; both listed under APOE in each `.genes.annot`. Additional sentinels (expected GRCh37 position read from the BIM; expected gene = gene.loc interval containing it): rs7903146 (TCF7L2), rs10455872 (LPA), rs12740374 (CELSR2/SORT1), rs11591147 (PCSK9), rs1333049 (9p21), rs3184504 (SH2B3), rs6025 (F5), rs1801133 (MTHFR), rs9939609 (FTO), rs699 (AGT). A sentinel with no row in an endpoint or no BIM entry is reported "not_testable" and does not by itself fail; any testable sentinel with the wrong position or gene fails G2. rs429358/rs7412 are mandatory.
G3 Genome-wide assignment check, per endpoint: 100,000 random (seed 20261104) SNP->gene assignments from the annotation file whose SNP has a BIM position: >=95% must lie inside the assigned gene's GRCh37 interval.
G4 (after MAGMA, per endpoint) logs complete with "End time"; genes analysed >=18,000; APOE NSNPS>=10.

## 5. Quality control reported (not gates)
Gene count, NSNPS distribution (median, quartiles), missingness vs the old run, duplicates. Canonical positive-control check (fixed now): in CAD or MI, at least 3 of {APOE, LPA, CDKN2B, CELSR2, LDLR} have gene P<1e-5; failure is reported prominently, not tuned.

## 6. Analysis (unchanged definitions; code adapted from `outputs/final_ctas_finngen_20261021/code/01_finngen.py`)
Per-gene ledger for 19,080 CTAS-table genes; support = n endpoints with P<0.01; nested sets q<0.05/q<0.01/q<0.001 (223/99/38), background q>=0.05; k>=1,2,3; one-sided Fisher vs background, BH across 9 tests; disjoint bands. Added (new statistics, fixed now): Wilson 95% CI for proportions; Woolf 95% CI for odds ratios (Haldane 0.5 if a cell is zero). Locus-representative companion for 91/46/21; strict-core 21 (Table 3 FinnGen column = n endpoints P<0.01, min P and its endpoint per representative as in `strict_core_21_locus_finngen.tsv`). Anchors reported explicitly: APOE, LPA, CELSR2, CDKN2B, LDLR, plus all 21 core representatives.
Old values (from `outputs/final_ctas_finngen_20261021/`) are shown beside corrected values and labelled `INVALID_UNCORRECTED_BUILD`.

## 7. Interpretation rules (fixed now)
- Enrichment "retained" if BH q<0.05 for the >=3-endpoint test in all three nested sets; "disappears" if none; otherwise "partly retained".
- "Increases/decreases": compare corrected vs old proportion and OR point estimates; flag as material when the corrected >=3-endpoint proportion in any CTAS set differs from the old by >=5 percentage points, or when the q<0.01 -> q<0.001 plateau (|difference| < 2 percentage points in the old run) changes by >=5 percentage points, or when any strict-core representative changes between "supported >=1 endpoint" and not, for >=3 loci.
- No threshold, window, N, endpoint, gene set or definition will change after corrected results are seen. Any necessary deviation goes in DEVIATIONS.tsv with its reason.
