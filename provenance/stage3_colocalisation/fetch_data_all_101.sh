#!/bin/bash
# fetch_data_all_101.sh
# Pull ±1 Mb cis-window eQTL slices for all 101 CTAS lead genes × 8 tissues.
#
# This is the unified fetch script for the full 101-locus coloc run.
# It supersedes the two-stage overlap fetch (fetch_data_all_overlap.sh +
# fetch_data_liver_kidney.sh).  All 8 tissues are fetched in a single pass.
# Existing slices in gtex_slices/ are reused automatically (skip-if-present
# logic preserved from the overlap scripts).
#
# eQTLGen slot note (identical to overlap pipeline):
#   eQTLGen whole blood is not available in the eQTL Catalogue tabix FTP.
#   The slot labeled "eqtlgen" uses QTS000026/QTD000439 (Schmiedel_2018
#   Tfh_memory, N≈73) for naming consistency with the overlap run.
#   See run_coloc_all_101.R header for the scientific caveat.
#
# Outputs:
#   coords/all_101_gene_coords.tsv     — gene, ensg, chr, start, end, wstart, wend
#   logs/coord_lookup_failures_101.log — genes whose Ensembl lookup failed
#   logs/tabix_errors_all_101.log
#   logs/eqtlgen_errors_all_101.log
#   results_all_101/all_101_download_qc.tsv
#
# Run from any directory:
#   bash ${CTAS_WORKSPACE:?set CTAS_WORKSPACE}/eqtl_vignette/fetch_data_all_101.sh

set -u

WORKDIR=${CTAS_WORKSPACE:?set CTAS_WORKSPACE}/eqtl_vignette
cd "$WORKDIR"
mkdir -p gtex_slices coords logs metadata results_all_101

# ── Sanity check ─────────────────────────────────────────────────────────────
if ! command -v tabix >/dev/null 2>&1; then
  echo "ERROR: tabix not found. Run: brew install htslib"; exit 1
fi
echo "tabix: $(tabix --version 2>&1 | head -1)"
echo

# ── Config ────────────────────────────────────────────────────────────────────
PAD=1000000
MANIFEST=metadata/all_101_lead_genes_manifest.tsv
COORD_FILE=coords/all_101_gene_coords.tsv

if [ ! -f "$MANIFEST" ]; then
  echo "ERROR: $MANIFEST not found."; exit 1
fi

# 6 GTEx tissues (QTS000015 study ID)
TISSUES_LABEL=(heart_LV heart_AA aorta coronary blood liver kidney_cortex)
TISSUES_QTD=(QTD000256 QTD000251 QTD000131 QTD000136 QTD000356 QTD000266 QTD000261)

# Fixed eQTLGen-slot URL: QTS000026/QTD000439 = Schmiedel_2018 Tfh_memory
EQTLGEN_URL="https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000026/QTD000439/QTD000439.all.tsv.gz"

# ── Step 1: Ensembl hg38 coordinate lookup ─────────────────────────────────
echo "=== [1/3] Ensembl hg38 coordinate lookup for 101 CTAS lead genes ==="
echo -e "gene\tensg\tchr\tstart\tend\twstart\twend" > "$COORD_FILE"

> logs/coord_lookup_failures_101.log

NGENES=0
while IFS=$'\t' read -r LOCUS_ID GENE TIER SCORE CHR LSTART LEND GENES LOO; do
  RESP=$(curl -sS --max-time 30 \
    -H "Accept: application/json" \
    "https://rest.ensembl.org/lookup/symbol/homo_sapiens/${GENE}?expand=0" 2>/dev/null)

  PARSED=$(echo "$RESP" | python3 -c "
import sys, json
try:
    d = json.load(sys.stdin)
    gid  = d.get('id','')
    seq  = d.get('seq_region_name','')
    st   = d.get('start','')
    en   = d.get('end','')
    print(gid, seq, st, en, sep='\t')
except:
    print('', '', '', '', sep='\t')
" 2>/dev/null)

  ENSG=$(echo "$PARSED" | cut -f1)
  GCHR=$(echo "$PARSED" | cut -f2)
  GSTART=$(echo "$PARSED" | cut -f3)
  GEND=$(echo "$PARSED" | cut -f4)

  if [ -z "$ENSG" ] || [ -z "$GCHR" ] || [ -z "$GSTART" ] || [ -z "$GEND" ]; then
    echo "  ${GENE}: FAILED (Ensembl lookup returned no data)" | tee -a logs/coord_lookup_failures_101.log
    continue
  fi

  # Skip non-autosomal
  if ! [[ "$GCHR" =~ ^[0-9]+$ ]]; then
    echo "  ${GENE}: SKIPPED (non-autosomal chr=${GCHR})" | tee -a logs/coord_lookup_failures_101.log
    continue
  fi

  WSTART=$((GSTART - PAD))
  [ $WSTART -lt 1 ] && WSTART=1
  WEND=$((GEND + PAD))

  echo "  ${GENE}: chr${GCHR}:${WSTART}-${WEND}  ${ENSG}"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$GENE" "$ENSG" "$GCHR" "$GSTART" "$GEND" "$WSTART" "$WEND" >> "$COORD_FILE"
  NGENES=$((NGENES + 1))

done < <(tail -n +2 "$MANIFEST")

echo
echo "Coords written: ${NGENES}/101 genes"
NFAIL=$(wc -l < logs/coord_lookup_failures_101.log 2>/dev/null || echo 0)
[ "$NFAIL" -gt 0 ] && echo "Failed lookups: ${NFAIL} (see logs/coord_lookup_failures_101.log)"
echo

# ── Step 2: GTEx tabix pulls (7 tissues × 101 genes) ───────────────────────
echo "=== [2/3] GTEx tabix pulls (7 tissues × up to 101 genes) ==="
echo "    Skips files already present in gtex_slices/"
echo
> logs/tabix_errors_all_101.log

NI=0
for L in "${TISSUES_LABEL[@]}"; do
  QTD="${TISSUES_QTD[$NI]}"
  # liver and kidney_cortex are in QTS000015 (GTEx v8) — same study as the others
  URL="https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/${QTD}/${QTD}.all.tsv.gz"
  echo "--- Tissue: ${L}  (${QTD}) ---"

  while IFS=$'\t' read -r GENE ENSG GCHR GSTART GEND WSTART WEND; do
    OUT="gtex_slices/${GENE}__${L}.tsv"
    if [ -s "$OUT" ]; then
      echo "  ${GENE}: already present ($(wc -l < "$OUT") rows) — skip"
      continue
    fi
    REGION="${GCHR}:${WSTART}-${WEND}"
    printf "  %-14s  %s ... " "${GENE}" "${REGION}"
    if tabix "$URL" "$REGION" > "$OUT" 2>>logs/tabix_errors_all_101.log; then
      ROWS=$(wc -l < "$OUT")
      SIZE=$(du -h "$OUT" | cut -f1)
      echo "OK  (${ROWS} rows, ${SIZE})"
    else
      echo "FAILED"
      rm -f "$OUT"
      echo "${GENE}__${L}: tabix FAILED" >> logs/tabix_errors_all_101.log
    fi
    sleep 2
  done < <(tail -n +2 "$COORD_FILE")

  NI=$((NI + 1))
  echo
done

# ── Step 3: eQTLGen-slot pulls (QTD000439 Schmiedel_2018 Tfh_memory) ────────
echo "=== [3/3] eQTLGen-slot pulls (QTD000439 / Schmiedel_2018 Tfh_memory) ==="
echo "    URL: ${EQTLGEN_URL}"
echo "    NOTE: This is NOT eQTLGen blood. See script header for explanation."
echo
> logs/eqtlgen_errors_all_101.log

while IFS=$'\t' read -r GENE ENSG GCHR GSTART GEND WSTART WEND; do
  OUT="gtex_slices/${GENE}__eqtlgen.tsv"
  if [ -s "$OUT" ]; then
    echo "  ${GENE}: already present ($(wc -l < "$OUT") rows) — skip"
    continue
  fi
  REGION="${GCHR}:${WSTART}-${WEND}"
  printf "  %-14s  %s ... " "${GENE}" "${REGION}"
  if tabix "$EQTLGEN_URL" "$REGION" > "$OUT" 2>>logs/eqtlgen_errors_all_101.log; then
    ROWS=$(wc -l < "$OUT")
    SIZE=$(du -h "$OUT" | cut -f1)
    echo "OK  (${ROWS} rows, ${SIZE})"
  else
    echo "FAILED"
    rm -f "$OUT"
    echo "${GENE}__eqtlgen: tabix FAILED" >> logs/eqtlgen_errors_all_101.log
  fi
  sleep 2
done < <(tail -n +2 "$COORD_FILE")

# ── Download QC summary ───────────────────────────────────────────────────────
echo
echo "================================================="
echo "Download QC Summary"
echo "================================================="

QC_FILE="results_all_101/all_101_download_qc.tsv"
ALL_TISSUES=("${TISSUES_LABEL[@]}" eqtlgen)
printf "gene\ttissue\tstatus\tnrows\n" > "$QC_FILE"

while IFS=$'\t' read -r GENE ENSG GCHR GSTART GEND WSTART WEND; do
  for T in "${ALL_TISSUES[@]}"; do
    F="gtex_slices/${GENE}__${T}.tsv"
    if [ -s "$F" ]; then
      ROWS=$(wc -l < "$F")
      STATUS="ok"
    else
      ROWS=0
      STATUS="missing"
    fi
    printf "%s\t%s\t%s\t%s\n" "$GENE" "$T" "$STATUS" "$ROWS" >> "$QC_FILE"
  done
done < <(tail -n +2 "$COORD_FILE")

echo "QC table: $QC_FILE"
echo

OK_CT=$(grep -c "	ok	" "$QC_FILE" || echo 0)
MISS_CT=$(grep -c "	missing	" "$QC_FILE" || echo 0)
TOTAL=$((OK_CT + MISS_CT))
echo "Slices OK:      $OK_CT / $TOTAL"
echo "Slices missing: $MISS_CT / $TOTAL"
echo

if [ -s logs/tabix_errors_all_101.log ]; then
  echo "GTEx tabix errors (tail):"
  tail -10 logs/tabix_errors_all_101.log
fi
if [ -s logs/eqtlgen_errors_all_101.log ]; then
  echo "eQTLGen errors (tail):"
  tail -10 logs/eqtlgen_errors_all_101.log
fi
if [ -s logs/coord_lookup_failures_101.log ]; then
  echo "Coord lookup failures:"
  cat logs/coord_lookup_failures_101.log
fi

echo
echo "=== DONE ==="
echo "Next step: Rscript ${CTAS_WORKSPACE:?set CTAS_WORKSPACE}/eqtl_vignette/run_coloc_all_101.R"
