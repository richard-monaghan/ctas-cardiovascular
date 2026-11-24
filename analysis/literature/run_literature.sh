#!/usr/bin/env bash
# Rebuilds the comparator (literature-overlap) analysis from the four publications' supplementary files.
# Usage: analysis/literature/run_literature.sh <sources_dir> [path/to/g1000_eur.bim]
#   <sources_dir> must contain (see inputs/literature/COMPARATOR_SOURCE_MANIFEST.tsv for checksums and retrieval):
#     qiao_2025/41467_2025_62419_MOESM19_ESM.xlsx   liu_2023_ten_cvd/mmc2.xlsx   liu_2023_ten_cvd/mmc3.xlsx
#     pleio_2021/mmc2.xlsx                           foley_2021/foley_2021_chd43_regions.tsv (transcribed; see FOLEY_TRANSCRIPTION.md)
#   The optional BIM (1000 Genomes EUR, GRCh37) regenerates the Foley SNP coordinates; otherwise the packaged file is used.
set -euo pipefail
cd "$(dirname "$0")/../.."; export CTAS_REPO="$PWD"; SRC="$(cd "$1" && pwd)"
W="$PWD/work/literature"; rm -rf "$W"; mkdir -p "$W/sources/processed" "$W/sources/qiao_2025" "$W/sources/liu_2023_ten_cvd" "$W/pleio/data_raw" "$W/data_raw"
export CTAS_LIT_WORK="$W" CTAS_LIT_SOURCES="$W/sources" CTAS_LIT_PLEIO="$W/pleio"
cp "$SRC/qiao_2025/41467_2025_62419_MOESM19_ESM.xlsx" "$W/sources/qiao_2025/"; cp "$SRC/liu_2023_ten_cvd/mmc2.xlsx" "$SRC/liu_2023_ten_cvd/mmc3.xlsx" "$W/sources/liu_2023_ten_cvd/"
cp "$SRC/foley_2021/foley_2021_chd43_regions.tsv" "$W/sources/processed/"; cp "$SRC/pleio_2021/mmc2.xlsx" "$W/pleio/data_raw/mmc2.xlsx"
if [ $# -ge 2 ]; then python3 analysis/literature/00_transform_sources.py "$SRC" "$W/sources" "$2"; mv "$W/sources/foley_snp_hg19.tsv" "$W/data_raw/"
else python3 analysis/literature/00_transform_sources.py "$SRC" "$W/sources"; cp inputs/literature/foley_snp_hg19.tsv "$W/data_raw/"; fi
gunzip -c inputs/literature/hgnc_complete_set_20260919.txt.gz > "$W/data_raw/hgnc_complete_set.txt"
python3 analysis/literature/01_extract_pleio.py
python3 analysis/literature/01_build_clean_inputs.py
python3 analysis/literature/02_match_and_results.py
python3 analysis/literature/03_pleio_permutation.py
python3 analysis/literature/04_four_study_union_final91.py
