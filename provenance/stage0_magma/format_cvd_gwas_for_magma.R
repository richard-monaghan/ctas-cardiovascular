#!/usr/bin/env Rscript
# Format Downloaded CVD GWAS for MAGMA Input
# Date: 2026-02-14

library(data.table)
library(dplyr)

# Directories
input_dir <- "<CTAS_WORKSPACE>/gwas_downloads_cvd"
output_dir <- "<CTAS_WORKSPACE>/gwas_magma_cvd"

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("========================================\n")
cat("FORMATTING CVD GWAS FOR MAGMA\n")
cat("========================================\n\n")

cat("Input directory:", input_dir, "\n")
cat("Output directory:", output_dir, "\n\n")

# Get all GWAS files
gwas_files <- list.files(input_dir, pattern = "\\.tsv\\.bgz$", full.names = TRUE)

if (length(gwas_files) == 0) {
  stop("ERROR: No GWAS files found in ", input_dir, "\n",
       "Please run download_cvd_gwas.sh first!")
}

cat("Found", length(gwas_files), "GWAS files to format\n\n")

# Track formatting success
success_count <- 0
fail_count <- 0
failed_files <- character()

# Process each GWAS file
for (i in seq_along(gwas_files)) {
  file <- gwas_files[i]
  trait_name <- basename(file)

  cat("---\n")
  cat(sprintf("[%d/%d] Processing: %s\n", i, length(gwas_files), trait_name))

  tryCatch({
    # Read GWAS file
    # Neale Lab format has columns: variant, minor_allele, minor_AF, low_confidence_variant,
    #                                n_complete_samples, AC, ytx, beta, se, tstat, pval
    cat("  Reading file...\n")

    # First, check file structure
    header <- fread(file, nrows = 1)
    cat("  Columns found:", paste(names(header), collapse = ", "), "\n")

    # Read full file with required columns
    if ("variant" %in% names(header) && "pval" %in% names(header)) {
      # Read variant and pval columns
      gwas <- fread(file, select = c("variant", "pval"))

      # Parse variant (format: chr:pos:ref:alt)
      cat("  Parsing variant IDs...\n")
      variant_parts <- tstrsplit(gwas$variant, ":", fixed = TRUE)

      gwas$CHR <- as.integer(variant_parts[[1]])
      gwas$BP <- as.integer(variant_parts[[2]])
      gwas$SNP <- gwas$variant  # Keep original variant ID

      # Select and rename for MAGMA
      gwas_magma <- gwas %>%
        select(SNP, CHR, BP, P = pval) %>%
        filter(!is.na(P), P > 0, !is.na(CHR), !is.na(BP))

      cat("  SNPs after filtering:", nrow(gwas_magma), "\n")

      if (nrow(gwas_magma) < 100000) {
        warning("  WARNING: Less than 100k SNPs retained - check filtering!")
      }

      # Create output filename
      outfile <- file.path(output_dir,
                          paste0(gsub("\\.tsv\\.bgz$", "", basename(file)), "_magma.txt"))

      # Write MAGMA input
      fwrite(gwas_magma, outfile, sep = "\t", quote = FALSE)

      cat("  ✓ Written:", nrow(gwas_magma), "SNPs to", basename(outfile), "\n")
      success_count <- success_count + 1

    } else {
      stop("Required columns 'variant' and 'pval' not found")
    }

  }, error = function(e) {
    cat("  ✗ ERROR:", conditionMessage(e), "\n")
    fail_count <- fail_count + 1
    failed_files <<- c(failed_files, basename(file))
  })

  cat("\n")
}

# Summary
cat("========================================\n")
cat("FORMATTING SUMMARY\n")
cat("========================================\n\n")

cat("Total files processed:", length(gwas_files), "\n")
cat("Successfully formatted:", success_count, "\n")
cat("Failed:", fail_count, "\n\n")

if (fail_count > 0) {
  cat("Failed files:\n")
  for (f in failed_files) {
    cat("  -", f, "\n")
  }
  cat("\n")
}

# List output files
output_files <- list.files(output_dir, pattern = "_magma\\.txt$")
cat("Output files created:", length(output_files), "\n")
cat("Location:", output_dir, "\n\n")

if (success_count == 19) {
  cat("✓ All 19 new CVD traits successfully formatted!\n\n")
  cat("Next step:\n")
  cat("  Run MAGMA analysis: bash scripts/run_magma_cvd_all.sh\n\n")
} else if (success_count > 0) {
  cat("⚠ Partial success - check errors above\n\n")
} else {
  cat("✗ No files formatted successfully - check input directory\n\n")
}

cat("========================================\n")
cat("Done!\n")
