#!/bin/bash
# Download CVD GWAS from Pan-UK Biobank
# Using pattern discovered from MI download
# Date: 2026-02-14

# Create download directory
DOWNLOAD_DIR="${CTAS_WORKSPACE:?set CTAS_WORKSPACE}/gwas_panukbb_cvd"
mkdir -p ${DOWNLOAD_DIR}
cd ${DOWNLOAD_DIR}

# Base URL from the MI example (corrected)
BASE_URL="https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files"

echo "========================================="
echo "DOWNLOADING CVD GWAS FROM PAN-UKBB"
echo "========================================="
echo ""
echo "Download directory: ${DOWNLOAD_DIR}"
echo ""

# Function to download
download_trait() {
  local code=$1
  local description=$2
  local priority=$3

  local filename="icd10-${code}-both_sexes.tsv.bgz"

  echo "---"
  echo "[${priority}] Downloading: ${description}"
  echo "Code: ${code}"
  echo "File: ${filename}"

  if [ -f "${filename}" ]; then
    SIZE=$(stat -f%z "${filename}" 2>/dev/null || stat -c%s "${filename}" 2>/dev/null)
    if [ ${SIZE} -gt 100000 ]; then
      echo "  ✓ Already exists ($(numfmt --to=iec ${SIZE} 2>/dev/null || echo ${SIZE} bytes))"
      return 0
    else
      echo "  ⚠ Exists but too small, re-downloading..."
      rm -f "${filename}"
    fi
  fi

  curl -f -L -C - -o "${filename}" "${BASE_URL}/${filename}"

  if [ $? -eq 0 ]; then
    SIZE=$(stat -f%z "${filename}" 2>/dev/null || stat -c%s "${filename}" 2>/dev/null)
    if [ ${SIZE} -gt 100000 ]; then
      echo "  ✓ Downloaded successfully ($(numfmt --to=iec ${SIZE} 2>/dev/null || echo ${SIZE} bytes))"
    else
      echo "  ✗ Download failed or file too small (${SIZE} bytes)"
      rm -f "${filename}"
      return 1
    fi
  else
    echo "  ✗ Download failed (404 or network error)"
    rm -f "${filename}"
    return 1
  fi
}

# =========================================
# HIGH PRIORITY - Ischemic Heart Disease & Stroke (7)
# =========================================

echo ""
echo "========================================="
echo "HIGH PRIORITY TRAITS"
echo "========================================="

download_trait "I21" "Myocardial Infarction" "HIGH"
download_trait "I20" "Angina Pectoris" "HIGH"
download_trait "I63" "Ischemic Stroke" "HIGH"
download_trait "I61" "Hemorrhagic Stroke" "HIGH"
download_trait "I64" "Stroke (All Types)" "HIGH"
download_trait "G45" "Transient Ischemic Attack (TIA)" "HIGH"

# Try both I82 (VTE) and I80/I26 separately
download_trait "I82" "Venous Thromboembolism" "HIGH"

# =========================================
# MEDIUM PRIORITY - Structural & Vascular (8)
# =========================================

echo ""
echo "========================================="
echo "MEDIUM PRIORITY TRAITS"
echo "========================================="

download_trait "I42" "Cardiomyopathy" "MEDIUM"
download_trait "I35" "Aortic Valve Disease" "MEDIUM"
download_trait "I73" "Peripheral Artery Disease" "MEDIUM"
download_trait "I71" "Aortic Aneurysm" "MEDIUM"
download_trait "I70" "Atherosclerosis" "MEDIUM"
download_trait "I34" "Mitral Valve Disease" "MEDIUM"
download_trait "I80" "Deep Vein Thrombosis" "MEDIUM"
download_trait "I26" "Pulmonary Embolism" "MEDIUM"

# =========================================
# Summary
# =========================================

echo ""
echo "========================================="
echo "DOWNLOAD SUMMARY"
echo "========================================="

DOWNLOADED=$(ls -1 icd10-*.tsv.bgz 2>/dev/null | wc -l | tr -d ' ')
echo "Files successfully downloaded: ${DOWNLOADED}"
echo ""

if [ ${DOWNLOADED} -gt 0 ]; then
  echo "✓ Downloaded ${DOWNLOADED} GWAS files"
  echo ""
  echo "Downloaded files with sizes:"
  ls -lh icd10-*.tsv.bgz 2>/dev/null | awk '{print $9, "\t", $5}'
else
  echo "⚠ No files downloaded successfully"
fi

echo ""
echo "========================================="
echo "NEXT STEPS"
echo "========================================="
echo ""

if [ ${DOWNLOADED} -ge 7 ]; then
  echo "✓ Excellent! Got ${DOWNLOADED} new CVD traits"
  echo ""
  echo "Next steps:"
  echo "1. Format for MAGMA: Rscript scripts/format_panukbb_for_magma.R"
  echo "2. Run MAGMA: bash scripts/run_magma_panukbb_cvd.sh"
  echo "3. Calculate CVD-CTAS with $(( 13 + DOWNLOADED )) total traits"
elif [ ${DOWNLOADED} -gt 0 ]; then
  echo "⚠ Got ${DOWNLOADED} files - may want to try manual download for missing"
else
  echo "✗ Download failed - check network or try manual download"
  echo "   URL pattern: ${BASE_URL}/icd10-{CODE}-both_sexes.tsv.bgz"
fi

echo ""
echo "Done!"
