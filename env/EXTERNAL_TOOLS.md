# External tools (versions as recorded)

| Tool | Version | Used for | Record |
|---|---|---|---|
| MAGMA | v1.10 ("custom" macOS build) | gene-level P (Pan-UKBB, FinnGen) | MAGMA log headers (`Welcome to MAGMA v1.10 (custom)`) |
| 1000 Genomes EUR reference | phase 3, 503 individuals | MAGMA LD reference | MAGMA logs |
| NCBI gene locations | NCBI37.3 | MAGMA annotation, locus coordinates | MAGMA logs |
| UCSC chain | hg38ToHg19 | FinnGen GRCh38→GRCh37 conversion (custom Python mapping; liftOver binary not used) | manuscript Methods; provenance/stage3_finngen/10_lift.py |
| UCSC liftOver | command-line binary, hg19ToHg38 chain | colocalisation coordinate conversion | manuscript Methods (version not recorded) |
| coloc | 5.2.3 | coloc.abf | manuscript Methods |
| Microsoft PowerPoint | macOS, AppleScript PDF export | final table/figure PDF rendering | provenance/stage4_tables/pptx_to_pdf.sh |
