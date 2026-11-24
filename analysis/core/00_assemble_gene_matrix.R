#!/usr/bin/env Rscript
################################################################################
# 00_assemble_gene_matrix.R -- Stage 0b: assemble the 19,080 x 16 gene-level P matrix
# from the 16 per-trait MAGMA v1.10 *.genes.out files (column P, used as-is; no
# re-derivation from ZSTAT). Written for the public release (no original assembly
# script survives); validated in QA to reproduce inputs/UKBB_16trait_gene_matrix.tsv
# byte-for-byte from the original MAGMA outputs.
# Usage: Rscript analysis/core/00_assemble_gene_matrix.R <dir_with_genes.out> [out.tsv]
################################################################################
args <- commandArgs(TRUE)
if (length(args) < 1) stop("usage: 00_assemble_gene_matrix.R <magma_dir> [out.tsv]")
mdir <- args[1]; out <- if (length(args) > 1) args[2] else "work/UKBB_16trait_gene_matrix.rebuilt.tsv"
TRAITS <- c("AF_I48", "Angina", "Aortic_Aneurysm", "Aortic_Valve", "Atherosclerosis", "CAD_I25", "Cardiomyopathy", "DVT",
            "HF_I50", "Hypertension_401", "Ischemic_Stroke", "MI", "Mitral_Valve", "PAD", "TIA", "VTE")   # internal keys; see inputs/INPUT_MANIFEST.md
tabs <- lapply(TRAITS, function(t) { x <- read.table(file.path(mdir, paste0(t, ".genes.out")), header = TRUE, colClasses = "character")
  setNames(x$P, x$GENE) })
genes <- sort(unique(unlist(lapply(tabs, names))))
M <- sapply(tabs, function(p) unname(p[genes])); colnames(M) <- TRAITS
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
write.table(data.frame(GENE = genes, M, check.names = FALSE), out, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat(sprintf("%d genes x %d traits -> %s\n", length(genes), length(TRAITS), out))
