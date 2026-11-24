# Foley et al. 2021: transcription of the 43 coronary heart disease regions

**Source:** Foley CN *et al.*, *Nature Communications* 2021, doi:10.1038/s41467-020-20885-8, main-text **Tables 1, 2 and 3** (the 43 CHD colocalisation regions). The article provides no machine-readable supplement for these regions, so they were transcribed by hand (May 2026). The transcription is not redistributed here, because the publisher licence could not be verified automatically (16 Nov 2026).

**To recreate it**, write a tab-separated file `foley_2021/foley_2021_chd43_regions.tsv` (UTF-8, Unix line endings) with this header and one row per region, in the order the regions appear in Tables 1, 2 and 3:

```
chr	locus	traits	candidate_snp	nearest_gene	known_chd_locus	ppfc	ppe	source_table
```

| Column | Content |
|---|---|
| chr | chromosome |
| locus | the locus/gene label as printed (for example `ABCG8, ABCG5`) |
| traits | the colocalising traits as printed (for example `CHD, LDL`) |
| candidate_snp | the candidate SNP rsID |
| nearest_gene | the nearest gene |
| known_chd_locus | `Yes`/`No` |
| ppfc | posterior probability of full colocalisation |
| ppe | posterior probability explained by the candidate SNP |
| source_table | `Table1`, `Table2` or `Table3` |

The file used for the paper has 43 rows (18 from Table1, 20 from Table2, 5 from Table3) and SHA-256 `d4666db11f3a957e346c53c707cbbf811414cedf763e8d896bf53ea323ef8cf6`. A matching checksum confirms an exact transcription.

The analysis uses only `candidate_snp` (GRCh37 position from the 1000 Genomes EUR reference, packaged as `foley_snp_hg19.tsv`) and the gene labels (secondary gene-symbol call). It does not use the posterior probabilities.
