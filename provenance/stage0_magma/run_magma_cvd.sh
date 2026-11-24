#!/bin/bash
################################################################################
# RUN MAGMA GENE-BASED ANALYSIS ON CVD TRAITS
#
# Processes 12 CVD disease traits for MAIN CVD-CTAS analysis
# (excludes I61 hemorrhagic stroke and I26 PE per categorization plan)
#
# Date: February 17, 2026
################################################################################

set -e

echo "═══════════════════════════════════════════════════════════════"
echo "MAGMA GENE-BASED ANALYSIS - CVD TRAITS"
echo "═══════════════════════════════════════════════════════════════"
echo ""

# Paths
BASE_DIR="${CTAS_WORKSPACE:?set CTAS_WORKSPACE}"
GWAS_DIR="${BASE_DIR}/gwas_magma_formatted"
OUTPUT_DIR="${BASE_DIR}/magma_cvd_results"
MAGMA="${BASE_DIR}/MAGMA_Mac/magma"
REF_DATA="${REFERENCE_DIR:?set REFERENCE_DIR}/reference/g1000_eur"
GENE_LOC="${REFERENCE_DIR:?set REFERENCE_DIR}/reference/NCBI37.3.gene.loc"

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}/logs"

# Verify files exist
if [ ! -f "${MAGMA}" ]; then
  echo "ERROR: MAGMA not found: ${MAGMA}"
  exit 1
fi

if [ ! -f "${GENE_LOC}" ]; then
  echo "ERROR: Gene location file not found: ${GENE_LOC}"
  exit 1
fi

if [ ! -f "${REF_DATA}.bed" ]; then
  echo "ERROR: LD reference not found: ${REF_DATA}.bed"
  exit 1
fi

echo "MAGMA executable: ${MAGMA}"
echo "Gene locations:   ${GENE_LOC}"
echo "LD reference:     ${REF_DATA}"
echo "Input directory:  ${GWAS_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# Define traits for MAIN CVD-CTAS (12 traits)
# Excludes I61 (hemorrhagic stroke) and I26 (PE)
TRAITS=(
  "icd10-I21-both_sexes:MI"
  "icd10-I20-both_sexes:Angina"
  "icd10-I63-both_sexes:Ischemic_Stroke"
  "icd10-G45-both_sexes:TIA"
  "icd10-I82-both_sexes:VTE"
  "icd10-I42-both_sexes:Cardiomyopathy"
  "icd10-I35-both_sexes:Aortic_Valve"
  "icd10-I34-both_sexes:Mitral_Valve"
  "icd10-I73-both_sexes:PAD"
  "icd10-I71-both_sexes:Aortic_Aneurysm"
  "icd10-I70-both_sexes:Atherosclerosis"
  "icd10-I80-both_sexes:DVT"
)

TOTAL=${#TRAITS[@]}
CURRENT=0
SUCCESSFUL=0
FAILED=0

echo "Processing ${TOTAL} CVD traits..."
echo ""

for TRAIT_INFO in "${TRAITS[@]}"; do
  CURRENT=$((CURRENT + 1))

  # Parse trait info
  TRAIT_FILE=$(echo ${TRAIT_INFO} | cut -d':' -f1)
  TRAIT_NAME=$(echo ${TRAIT_INFO} | cut -d':' -f2)

  INPUT_FILE="${GWAS_DIR}/${TRAIT_FILE}.magma.gz"
  OUTPUT_PREFIX="${OUTPUT_DIR}/${TRAIT_NAME}"

  echo "═══════════════════════════════════════════════════════════════"
  echo "[${CURRENT}/${TOTAL}] ${TRAIT_NAME}"
  echo "═══════════════════════════════════════════════════════════════"
  echo ""

  # Check if input exists
  if [ ! -f "${INPUT_FILE}" ]; then
    echo "  ✗ Input file not found: ${INPUT_FILE}"
    echo ""
    FAILED=$((FAILED + 1))
    continue
  fi

  # Skip if output already exists
  if [ -f "${OUTPUT_PREFIX}.genes.out" ]; then
    GENES=$(tail -n +2 "${OUTPUT_PREFIX}.genes.out" | wc -l | tr -d ' ')
    if [ ${GENES} -gt 10000 ]; then
      echo "  ✓ Already completed (${GENES} genes)"
      echo ""
      SUCCESSFUL=$((SUCCESSFUL + 1))
      continue
    else
      echo "  ⚠ Output exists but incomplete, re-running..."
      rm -f "${OUTPUT_PREFIX}".*
    fi
  fi

  # Uncompress input for MAGMA
  echo "  Uncompressing input..."
  INPUT_UNCOMPRESSED="${OUTPUT_DIR}/${TRAIT_NAME}_input.magma"
  gunzip -c "${INPUT_FILE}" > "${INPUT_UNCOMPRESSED}"

  if [ ! -f "${INPUT_UNCOMPRESSED}" ]; then
    echo "  ✗ Failed to uncompress input"
    echo ""
    FAILED=$((FAILED + 1))
    continue
  fi

  VARIANTS=$(tail -n +2 "${INPUT_UNCOMPRESSED}" | wc -l | tr -d ' ')
  echo "  Input: ${VARIANTS} variants"
  echo ""

  ################################################################################
  # STEP 0.5: Convert chr:pos:ref:alt IDs to rsIDs
  ################################################################################

  echo "  Step 0.5: Converting SNP IDs to rsIDs..."
  INPUT_RSIDS="${OUTPUT_DIR}/${TRAIT_NAME}_rsids.magma"

  bash "${BASE_DIR}/scripts/create_gwas_with_rsids.sh" \
    "${INPUT_UNCOMPRESSED}" \
    "${REF_DATA}.bim" \
    "${INPUT_RSIDS}" \
    > "${OUTPUT_DIR}/logs/${TRAIT_NAME}_rsid_conversion.log" 2>&1

  if [ ! -f "${INPUT_RSIDS}" ]; then
    echo "  ✗ rsID conversion failed"
    echo ""
    rm -f "${INPUT_UNCOMPRESSED}"
    FAILED=$((FAILED + 1))
    continue
  fi

  RSID_VARIANTS=$(wc -l < "${INPUT_RSIDS}" | tr -d ' ')
  echo "  ✓ Converted ${RSID_VARIANTS} variants to rsIDs"
  echo ""

  ################################################################################
  # STEP 1: Annotate SNPs to genes
  ################################################################################

  echo "  Step 1: Annotating SNPs to genes..."

  "${MAGMA}" \
    --annotate window=0,0 \
    --snp-loc "${INPUT_RSIDS}" \
    --gene-loc "${GENE_LOC}" \
    --out "${OUTPUT_PREFIX}_annot" \
    > "${OUTPUT_DIR}/logs/${TRAIT_NAME}_annot.log" 2>&1

  if [ $? -ne 0 ]; then
    echo "  ✗ Annotation failed (check logs/${TRAIT_NAME}_annot.log)"
    echo ""
    rm -f "${INPUT_UNCOMPRESSED}"
    FAILED=$((FAILED + 1))
    continue
  fi

  if [ ! -f "${OUTPUT_PREFIX}_annot.genes.annot" ]; then
    echo "  ✗ Annotation file not created"
    echo ""
    rm -f "${INPUT_UNCOMPRESSED}"
    FAILED=$((FAILED + 1))
    continue
  fi

  ANNOTATED_GENES=$(tail -n +2 "${OUTPUT_PREFIX}_annot.genes.annot" | wc -l | tr -d ' ')
  echo "  ✓ Annotated ${ANNOTATED_GENES} genes"
  echo ""

  ################################################################################
  # STEP 2: Run gene-based analysis
  ################################################################################

  echo "  Step 2: Running gene-based analysis (with LD correction)..."

  # Sample size for Pan-UKBB (conservative estimate)
  SAMPLE_SIZE=500000

  "${MAGMA}" \
    --bfile "${REF_DATA}" \
    --pval "${INPUT_RSIDS}" use=1,4 N=${SAMPLE_SIZE} \
    --gene-annot "${OUTPUT_PREFIX}_annot.genes.annot" \
    --out "${OUTPUT_PREFIX}" \
    > "${OUTPUT_DIR}/logs/${TRAIT_NAME}.log" 2>&1

  if [ $? -ne 0 ]; then
    echo "  ✗ MAGMA analysis failed (check logs/${TRAIT_NAME}.log)"
    echo ""
    rm -f "${INPUT_UNCOMPRESSED}"
    FAILED=$((FAILED + 1))
    continue
  fi

  # Cleanup temporary files
  rm -f "${INPUT_UNCOMPRESSED}"
  rm -f "${INPUT_RSIDS}"

  # Verify output
  if [ -f "${OUTPUT_PREFIX}.genes.out" ]; then
    TOTAL_GENES=$(tail -n +2 "${OUTPUT_PREFIX}.genes.out" | wc -l | tr -d ' ')
    SIG_GENES=$(awk '$9 < 2.5e-6' "${OUTPUT_PREFIX}.genes.out" | wc -l)

    echo "  ✓ Analysis complete"
    echo "    Total genes: ${TOTAL_GENES}"
    echo "    Significant (p<2.5e-6): ${SIG_GENES}"
    echo ""

    SUCCESSFUL=$((SUCCESSFUL + 1))
  else
    echo "  ✗ Output file not created"
    echo ""
    FAILED=$((FAILED + 1))
  fi

done

echo "═══════════════════════════════════════════════════════════════"
echo "SUMMARY"
echo "═══════════════════════════════════════════════════════════════"
echo ""
echo "Total traits:  ${TOTAL}"
echo "Successful:    ${SUCCESSFUL}"
echo "Failed:        ${FAILED}"
echo ""

if [ ${SUCCESSFUL} -eq ${TOTAL} ]; then
  echo "✓ All ${TOTAL} traits processed successfully!"
elif [ ${SUCCESSFUL} -gt 0 ]; then
  echo "⚠ ${SUCCESSFUL}/${TOTAL} traits completed (${FAILED} failed)"
else
  echo "✗ All analyses failed"
  exit 1
fi

echo ""
echo "Output directory: ${OUTPUT_DIR}"
echo ""
echo "Output files:"
ls -lh "${OUTPUT_DIR}"/*.genes.out 2>/dev/null | awk '{print "  " $9 " (" $5 ")"}'

echo ""
echo "═══════════════════════════════════════════════════════════════"
echo "NEXT STEPS"
echo "═══════════════════════════════════════════════════════════════"
echo ""
echo "1. Combine with existing 4 CVD trait MAGMA results"
echo "   - CAD (I25), AF (I48), HF (I50), Hypertension (401)"
echo ""
echo "2. Create gene × trait matrix for 16 traits"
echo ""
echo "3. Calculate CVD-CTAS scores"
echo ""
echo "Done!"
