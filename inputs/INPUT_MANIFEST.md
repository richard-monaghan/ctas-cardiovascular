# Input manifest

## A. Packaged (derived, redistributable)

| File | What it is | Provenance | SHA-256 |
|---|---|---|---|
| `UKBB_16trait_gene_matrix.tsv` | 19,080 genes × 16 traits. One-sided MAGMA v1.10 gene-level P (column `P`, not re-derived from Z). `NA` = gene absent from that trait's MAGMA output (3,144 cells). | Byte copy of the frozen input recorded in the primary freeze (`input_checksums.tsv`). Re-assembled byte-for-byte from the 16 `*.genes.out` files by `analysis/core/00_assemble_gene_matrix.R` (release QA). | `dfbbcc240350c656f4917adfaf23256892425706eb27bfdea83bb595e44b87b6` |
| `sigma_A_raw_spearman.csv` | Σ_A, the 16 × 16 Spearman correlation of normal scores (pairwise complete) | Byte copy of the frozen cache used by the primary freeze. Recomputed from the matrix by `sigma_A()` to within 5×10⁻¹⁶. | `a5bd0f610cc36c098d85965f46d7452a650b5ff103655e1db46c69399c963799` |
| `locus_map/gene_neighbourhood_map.csv` | 1,145 genes → stable neighbourhood id (`ST2-n` or `NEW-n`), with GRCh37 coordinates (NCBI37.3) | Column subset (`gene, CHR, START, END, locus_id`) of the locked `gene_locus_assignment.csv` (SHA `b3cc0621…`). The superseded score columns were dropped; values are unchanged. | see `CHECKSUMS.sha256` |
| `locus_map/stable_to_canonical_locus.csv` | 104 stable-id → canonical `V2Lnnn` pairs | Unique pairs from the locked `v2_gene_to_locus.csv` (SHA `bc97aee0…`) | see `CHECKSUMS.sha256` |
| `locus_map/canonical_locus_master.csv` | 104 canonical loci with coordinates and legacy identifiers | Column subset of the locked `v2_locus_master.csv` (SHA `f8a47252…`). The superseded tiers and score columns from the pre-final lock were dropped. | see `CHECKSUMS.sha256` |

## A2. Inputs added to complete reproduction (16 Nov 2026)

| Path | What it is | Label |
|---|---|---|
| `locus_map/discovery_stage/anchor_regions_earlier_release_101.tsv` | 101 anchor regions (identifier, chromosome, start, end) from an earlier release of this analysis | Frozen discovery-stage input (coordinates only) |
| `locus_map/discovery_stage/gene_universe_1145.tsv` | the 1,145 genes covered by the locus map, with the pre-final q<0.05 flag that selected the 104 canonical loci | Frozen discovery-stage input (identifiers and flag only) |
| `locus_map/discovery_stage/gene_coordinates_NCBI37.3_subset.tsv` | NCBI37.3 (GRCh37) coordinates of those genes (MAGMA gene-location file) | Public reference subset |
| `historical_prefinal/prefinal_*.csv` (4 files) | Pre-final gene table, gene→locus assignment, gene→locus and locus master, reduced to the columns read by the sensitivity and comparator engines | **Historical / frozen intermediates** (pre-final Σ_A, 10⁶-draw analysis of Aug 2026). Used only by the sensitivity and comparator analyses, as in the paper. They are not the final results. |
| `effective_n/gene_matrix_effective_n.tsv`, `effective_n/effective_n_table.csv` | Gene × trait MAGMA P matrix rerun with trait-specific effective N, and the N table | Frozen derived input (MAGMA rerun on public Pan-UKBB data; `provenance/stage0_effective_n_magma/`) |
| `finngen_r12/` | 10 FinnGen R12 MAGMA gene-level outputs and two corrected-only summaries | Frozen derived input; see `finngen_r12/README_FINNGEN.md` (acknowledgement and citation) |
| `colocalisation/` | Cell-level coloc.abf results, locus abstraction, target coordinates/windows, manifest, the original method, sign-off and validation notes, and `EQTL_SLICE_MANIFEST.tsv` (URL, dataset, window, access date and SHA-256 for all 728 eQTL slices) | Frozen derived input (7 Aug 2026 run); see `analysis/colocalisation/README_COLOC.md` |
| `literature/hgnc_complete_set_20260919.txt.gz` | HGNC complete set as downloaded on 19 Sep 2026 (latest record 2026-09-18) | Public reference, CC0 (genenames.org); packaged because no archive holds this exact snapshot |
| `literature/foley_snp_hg19.tsv` | GRCh37 positions of the 43 Foley 2021 candidate SNPs (1000 Genomes EUR reference) | Derived; regenerable with `analysis/literature/00_transform_sources.py` |
| `literature/COMPARATOR_SOURCE_MANIFEST.tsv`, `literature/FOLEY_TRANSCRIPTION.md` | Exact source tables, selection rules, checksums and retrieval/transcription instructions for the four comparator studies | Documentation. **The comparator source tables are not redistributed.** |

The 10⁶-draw historical reference null used by the effective-N analysis is **not** packaged. `analysis/sensitivity/00_reference_null_1e6.R` regenerates it byte-identically from `sigma_A_raw_spearman.csv` and the gene matrix (seed 42).

**Locus map version 2026-08-05 (locked).** This is a versioned derived input.

How it was built:
1. **Gene coordinates.** Each gene is placed at its NCBI37.3 (GRCh37) coordinates.
2. **Anchored loci.** Genes within 1 Mb of one of 101 anchor regions are assigned to that anchor's locus. The anchors are the locus set of an earlier release of this analysis (identifier `ST2-n`).
3. **New loci.** Remaining genes are merged in sequence into new loci (`NEW-n`) whenever the next gene lies within 1 Mb of the current locus end.
4. **Gene universe.** The map covers 1,145 genes that were significant in, or candidates for, the earlier analyses.
5. **Canonical loci.** The 104 loci containing at least one gene at q < 0.05 in the pre-final analysis were given canonical identifiers `V2L001`–`V2L104`, after a deterministic sort by chromosome, start, end and representative gene.

The final 91 loci all lie within these 104 (verified, 0 outside).

`analysis/locus_map/build_locus_map.R` rebuilds all three files byte-identically from the minimal discovery-stage inputs above, and the 91/46/21 loci follow from the rebuilt map (`tests/reproduce_packaged_outputs.sh`). The original builders are archived in `provenance/stage1_locus_map/`. The column `legacy_v1_locus_id` in the frozen tables holds the earlier release's anchor identifier.

## B. Trait keys

| Internal key | Paper label | Pan-UKBB phenotype |
|---|---|---|
| AF_I48 | AF, atrial fibrillation | icd10-I48 (both sexes) |
| Angina | Ang | icd10-I20 |
| Aortic_Aneurysm | AoAn | icd10-I71 |
| Aortic_Valve | AoV | icd10-I35 |
| Atherosclerosis | Ath | icd10-I70 |
| CAD_I25 | CAD | icd10-I25 |
| Cardiomyopathy | CMP | icd10-I42 |
| **DVT** | **I80**, phlebitis/thrombophlebitis (complete I80 category) | icd10-I80 |
| HF_I50 | HF | icd10-I50 |
| Hypertension_401 | HTN | phecode-401 |
| Ischemic_Stroke | IS | icd10-I63 |
| MI | MI | icd10-I21 |
| Mitral_Valve | MV | icd10-I34 |
| **PAD** | **I73**, other peripheral vascular disease | icd10-I73 |
| TIA | TIA | icd10-G45 |
| **VTE** | **I82**, venous embolism/thrombosis (complete I82 category) | icd10-I82 |

The keys DVT, PAD and VTE are historical identifiers. They are kept in code and data so that the files reproduce exactly. They are **not** clinical labels; use the paper labels.

## C. Obtainable from the original source, not redistributed

| Input | Needed for | Where to obtain | Why not included |
|---|---|---|---|
| Pan-UKBB per-phenotype summary statistics (16 files, GRCh37) | Step 0a (MAGMA) | `https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/` (see `provenance/stage0_magma/download_cvd_panukbb_final.sh`) | Raw GWAS files: bulky and available from the source (release decision) |
| 1000 Genomes phase 3 EUR PLINK reference (`g1000_eur`) | MAGMA LD | MAGMA website (CTG, VU Amsterdam) | Third-party reference |
| NCBI37.3 gene locations (`NCBI37.3.gene.loc`) | MAGMA annotation | MAGMA website | Third-party reference. Required SHA-256 `c628e49f6674492c61140d3df8b850cd1bf963f87df505c175d6ada743b5d845` (19,427 lines; `inputs/EXTERNAL_INPUT_CHECKSUMS.tsv`); pass it via `CTAS_NCBI37_GENE_LOC`. |
| MAGMA v1.10 binary | Step 0a; FinnGen gene analysis | https://ctg.cncr.nl/software/magma | Third-party licence |
| FinnGen R12 summary statistics (10 endpoints) | FinnGen external-support stage | https://www.finngen.fi/en/access_results | Raw GWAS files (release decision); subject to FinnGen terms |
| UCSC `hg38ToHg19.over.chain` | FinnGen GRCh38 → GRCh37 conversion | UCSC Genome Browser downloads | Third-party |
| GTEx v8 cis-eQTL summary data (via the eQTL Catalogue) and the per-cell coloc.abf results | Colocalisation stage | eQTL Catalogue (EMBL-EBI) | Bulky third-party data; upstream coloc run scripts not packaged (see release notes) |
| Comparator study tables (Qiao 2025, Liu 2023, Foley 2021, PLEIO) | Literature-overlap stage | The original publications' supplements | Third-party copyright |

## D. Not redistributed and not needed
- Saved 10⁹-draw null objects (`run_{A,C}_seed{42,1,2}_M1e+09.rds`, ~87 MB each). Regenerate them with step 2.
- Individual-level data of any kind: none was used.
