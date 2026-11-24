# Licence scope

This repository uses two licences for material created by the CTAS project. Third-party material keeps its own terms.

| Material | Licence |
|---|---|
| **Software and source code**: R, Python and shell code in `R/`, `analysis/`, `tests/`, `provenance/`, `config/`, `env/`, `run_core_pipeline.sh` and the Explorer's code (`explorer/site/assets/`, `explorer/site/scripts/`, `explorer/site/*.html`) | MIT. Copyright (c) 2026 The University of Manchester. See `LICENSE`. |
| **Project-generated derived results and documentation**: CTAS scores, P and q values, gene and locus sets, supplementary and main tables, figure source data, Explorer data files (`explorer/site/data/`), README and other documentation written for this project | CC BY 4.0. See `LICENSE-CC-BY-4.0.txt`. Please attribute by citing the article and the Zenodo archive (`CITATION.cff`). |
| **Third-party inputs, and the portions of derived files that reproduce them** | Their original terms, listed below. This repository does not relicense them and claims no ownership of them. |

## Third-party material and upstream terms

| Source | What is included here | Upstream terms |
|---|---|---|
| Pan-UK Biobank summary statistics (https://pan.ukbb.broadinstitute.org) | Gene-level MAGMA results derived from them (`inputs/UKBB_16trait_gene_matrix.tsv` and downstream results). No variant-level summary statistics. | CC BY 4.0. Pan-UKB summary statistics were generated using the UK Biobank Resource (project ID 31063). |
| FinnGen release 12 public summary statistics (https://www.finngen.fi/en/access_results) | Gene-level MAGMA results and summaries derived from them (`inputs/finngen_r12/`, Supplementary Table 11, Explorer support flags). No variant-level data. | FinnGen's terms of use. Required acknowledgement: "We want to acknowledge the participants and investigators of the FinnGen study". Cite Kurki *et al.*, *Nature* 613, 508–518 (2023). See `inputs/finngen_r12/README_FINNGEN.md`. |
| eQTL Catalogue / GTEx v8 cis-eQTL summary statistics (https://www.ebi.ac.uk/eqtl/) | Colocalisation results derived from them and slice manifests. No eQTL summary statistics. | eQTL Catalogue: CC BY 4.0. GTEx Portal terms apply to GTEx data. |
| HGNC gene nomenclature (https://www.genenames.org) | `inputs/literature/hgnc_complete_set_20260919.txt.gz` (unmodified snapshot) | CC0 |
| NCBI37.3 gene locations, 1000 Genomes phase 3 reference, MAGMA v1.10 | Not redistributed. Users obtain them from the original sources (`inputs/INPUT_MANIFEST.md`, `inputs/EXTERNAL_INPUT_CHECKSUMS.tsv`). | Original providers' terms |
| Comparator publications (Foley 2021, Liu 2023, Qiao 2025, PLEIO 2021) | Not redistributed. Only selection rules, checksums and retrieval instructions are included (`inputs/literature/`). | Publishers' terms |

If a file combines project-generated results with third-party-derived values, the CC BY 4.0 licence covers only the project's contribution. The third-party-derived values remain subject to the upstream terms above.
