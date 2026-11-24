#!/bin/bash

################################################################################
# Process 4 existing CVD traits with harmonized 19K gene annotation
# FIXED: Uses correct MAGMA syntax (use=1,4 instead of ncol=4)
################################################################################

set -e

echo "═══════════════════════════════════════════════════════════════"
echo "Processing 4 existing CVD traits with 19K gene annotation"
echo "═══════════════════════════════════════════════════════════════"
echo ""

# Directories
GWAS_DIR="${CTAS_WORKSPACE:?set CTAS_WORKSPACE}/gwas_data_existing4"
FORMATTED_DIR="${CTAS_WORKSPACE:?set CTAS_WORKSPACE}/gwas_magma_formatted_existing4"
OUTPUT_DIR="${CTAS_WORKSPACE:?set CTAS_WORKSPACE}/magma_cvd_results"

# Reference files
GENE_LOC="${REFERENCE_DIR:?set REFERENCE_DIR}/reference/NCBI37.3.gene.loc"
BFILE="${REFERENCE_DIR:?set REFERENCE_DIR}/reference/g1000_eur"
MAGMA="${CTAS_WORKSPACE:?set CTAS_WORKSPACE}/MAGMA_Mac/magma"

# Create directories
mkdir -p "$FORMATTED_DIR"

# Define traits (name:filename_base:sample_size)
declare -a TRAITS=(
  "CAD_I25:icd10-I25-both_sexes:500000"
  "AF_I48:icd10-I48-both_sexes:500000"
  "HF_I50:icd10-I50-both_sexes:500000"
  "Hypertension_401:phecode-401-both_sexes:500000"
)

################################################################################
# Process each trait
################################################################################

for TRAIT_INFO in "${TRAITS[@]}"; do
  NAME=$(echo "$TRAIT_INFO" | cut -d: -f1)
  TRAIT=$(echo "$TRAIT_INFO" | cut -d: -f2)
  N=$(echo "$TRAIT_INFO" | cut -d: -f3)

  echo "──────────────────────────────────────────────────────────────────"
  echo "Processing: $NAME"
  echo "──────────────────────────────────────────────────────────────────"

  # Skip if already completed
  if [ -f "${OUTPUT_DIR}/${NAME}.genes.out" ]; then
    GENES=$(tail -n +2 "${OUTPUT_DIR}/${NAME}.genes.out" | wc -l | tr -d ' ')
    echo "✓ Already completed ($GENES genes)"
    echo ""
    continue
  fi

  # Check for input file (.tsv.gz or .tsv.bgz)
  if [ -f "${GWAS_DIR}/${TRAIT}.tsv.gz" ]; then
    GWAS_FILE="${GWAS_DIR}/${TRAIT}.tsv.gz"
    DECOMPRESS="gunzip -c"
  elif [ -f "${GWAS_DIR}/${TRAIT}.tsv.bgz" ]; then
    GWAS_FILE="${GWAS_DIR}/${TRAIT}.tsv.bgz"
    DECOMPRESS="gunzip -c"
  else
    echo "ERROR: GWAS file not found for $TRAIT"
    echo "  Looked for: ${GWAS_DIR}/${TRAIT}.tsv.{gz,bgz}"
    continue
  fi

  FORMATTED_FILE="${FORMATTED_DIR}/${NAME}.magma"

  # Step 1: Format GWAS file for MAGMA
  if [ ! -f "${FORMATTED_FILE}" ]; then
    echo "Step 1: Formatting GWAS file..."
    echo "  Input: $(basename $GWAS_FILE) ($(du -h $GWAS_FILE | cut -f1))"

    $DECOMPRESS "$GWAS_FILE" | awk 'BEGIN {FS="\t"; OFS="\t"}
    NR==1 {
      # Find column indices
      for (i=1; i<=NF; i++) {
        if ($i == "chr") chr_col = i
        if ($i == "pos") pos_col = i
        if ($i == "ref") ref_col = i
        if ($i == "alt") alt_col = i
        if ($i == "neglog10_pval_meta_hq") pval_col = i
        if ($i == "neglog10_pval_meta" && pval_col == "") pval_col = i
      }
      print "SNP", "CHR", "POS", "PVAL"
      next
    }
    {
      chr = $chr_col
      pos = $pos_col
      ref = $ref_col
      alt = $alt_col
      neglog10_pval = $pval_col

      # Skip if missing critical fields or pval
      if (chr == "" || pos == "" || ref == "" || alt == "" || neglog10_pval == "NA") next

      # Convert neglog10 to pval
      pval = 10^(-neglog10_pval)

      # Skip invalid pvals
      if (pval <= 0 || pval > 1) next

      # Create SNP ID
      snp_id = chr ":" pos ":" ref ":" alt

      print snp_id, chr, pos, pval
    }' > "${FORMATTED_FILE}.tmp"

    # Check if we got any data
    LINE_COUNT=$(wc -l < "${FORMATTED_FILE}.tmp" | tr -d ' ')
    if [ "$LINE_COUNT" -lt 100 ]; then
      echo "ERROR: Only $LINE_COUNT lines in formatted file!"
      rm "${FORMATTED_FILE}.tmp"
      exit 1
    fi

    mv "${FORMATTED_FILE}.tmp" "$FORMATTED_FILE"
    echo "  ✓ Formatted $(($LINE_COUNT - 1)) SNPs"
  else
    LINE_COUNT=$(wc -l < "$FORMATTED_FILE" | tr -d ' ')
    echo "Step 1: Using existing formatted file ($(($LINE_COUNT - 1)) SNPs)"
  fi

  # Step 2: Convert SNP IDs to rsIDs
  RSID_FILE="${OUTPUT_DIR}/${NAME}_rsids.magma"
  if [ ! -f "$RSID_FILE" ]; then
    echo "Step 2: Converting SNP IDs to rsIDs..."

    bash ${CTAS_WORKSPACE:?set CTAS_WORKSPACE}/scripts/create_gwas_with_rsids.sh \
      "$FORMATTED_FILE" \
      "${BFILE}.bim" \
      "$RSID_FILE" \
      > "${OUTPUT_DIR}/${NAME}_rsid_conversion.log" 2>&1

    if [ ! -f "$RSID_FILE" ]; then
      echo "ERROR: rsID conversion failed! Check ${OUTPUT_DIR}/${NAME}_rsid_conversion.log"
      exit 1
    fi

    RSID_COUNT=$(wc -l < "$RSID_FILE" | tr -d ' ')
    echo "  ✓ Converted to rsIDs ($RSID_COUNT SNPs)"
  else
    RSID_COUNT=$(wc -l < "$RSID_FILE" | tr -d ' ')
    echo "Step 2: Using existing rsID file ($RSID_COUNT SNPs)"
  fi

  # Step 3: Annotate SNPs to genes (19K genes)
  if [ ! -f "${OUTPUT_DIR}/${NAME}_annot.genes.annot" ]; then
    echo "Step 3: Annotating SNPs to genes (19K annotation)..."

    "$MAGMA" \
      --annotate \
      --snp-loc "$RSID_FILE" \
      --gene-loc "$GENE_LOC" \
      --out "${OUTPUT_DIR}/${NAME}_annot" \
      > "${OUTPUT_DIR}/${NAME}_annot.log" 2>&1

    if [ ! -f "${OUTPUT_DIR}/${NAME}_annot.genes.annot" ]; then
      echo "ERROR: Annotation failed! Check ${OUTPUT_DIR}/${NAME}_annot.log"
      exit 1
    fi

    ANNOT_GENES=$(tail -n +2 "${OUTPUT_DIR}/${NAME}_annot.log" | grep "gene definitions read" | awk '{print $1}')
    echo "  ✓ Annotated to $ANNOT_GENES genes"
  else
    echo "Step 3: Using existing annotation"
  fi

  # Step 4: Run gene-based analysis (FIXED: use=1,4 instead of ncol=4)
  echo "Step 4: Running gene-based analysis (N=$N)..."

  "$MAGMA" \
    --bfile "$BFILE" \
    --pval "$RSID_FILE" use=1,4 N=$N \
    --gene-annot "${OUTPUT_DIR}/${NAME}_annot.genes.annot" \
    --out "${OUTPUT_DIR}/${NAME}" \
    > "${OUTPUT_DIR}/${NAME}.log" 2>&1

  if [ -f "${OUTPUT_DIR}/${NAME}.genes.out" ]; then
    GENES=$(tail -n +2 "${OUTPUT_DIR}/${NAME}.genes.out" | wc -l | tr -d ' ')
    MIN_P=$(awk 'NR>1 {if (min=="" || $9<min) min=$9} END {printf "%.2e", min}' "${OUTPUT_DIR}/${NAME}.genes.out")
    SIG_GENES=$(awk 'NR>1 && $9 < 2.5e-6' "${OUTPUT_DIR}/${NAME}.genes.out" | wc -l | tr -d ' ')

    echo "  ✓ Complete"
    echo "    - Genes analyzed: $GENES"
    echo "    - Min P-value: $MIN_P"
    echo "    - Significant genes (p<2.5e-6): $SIG_GENES"
  else
    echo "ERROR: Gene analysis failed! Check ${OUTPUT_DIR}/${NAME}.log"
    cat "${OUTPUT_DIR}/${NAME}.log"
    exit 1
  fi

  echo ""
done

################################################################################
# QC Summary
################################################################################

echo "═══════════════════════════════════════════════════════════════"
echo "QC SUMMARY - 4 Existing CVD Traits (19K gene annotation)"
echo "═══════════════════════════════════════════════════════════════"
echo ""

printf "%-20s %10s %15s %15s\n" "Trait" "Genes" "Min_P" "Sig_Genes"
echo "────────────────────────────────────────────────────────────────"

for TRAIT_INFO in "${TRAITS[@]}"; do
  NAME=$(echo "$TRAIT_INFO" | cut -d: -f1)
  OUTFILE="${OUTPUT_DIR}/${NAME}.genes.out"

  if [ -f "$OUTFILE" ]; then
    GENES=$(tail -n +2 "$OUTFILE" | wc -l | tr -d ' ')
    MIN_P=$(awk 'NR>1 {if (min=="" || $9<min) min=$9} END {printf "%.2e", min}' "$OUTFILE")
    SIG_GENES=$(awk 'NR>1 && $9 < 2.5e-6' "$OUTFILE" | wc -l | tr -d ' ')

    printf "%-20s %10s %15s %15s\n" "$NAME" "$GENES" "$MIN_P" "$SIG_GENES"
  else
    printf "%-20s %10s %15s %15s\n" "$NAME" "MISSING" "-" "-"
  fi
done

echo ""
echo "✓ Processing complete"
echo "  Output directory: $OUTPUT_DIR"
echo "  Results: CAD_I25.genes.out, AF_I48.genes.out, HF_I50.genes.out, Hypertension_401.genes.out"
echo "═══════════════════════════════════════════════════════════════"
