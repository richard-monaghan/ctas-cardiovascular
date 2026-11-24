#!/usr/bin/env Rscript
################################################################################
# 01_sigma_matrices.R -- Stage 1: trait-dependence matrices (new wrapper for release).
# Recomputes Sigma_A (raw Spearman of normal scores; canonical sigma_A() in R/common.R)
# from the gene matrix, checks it against the frozen inputs/sigma_A_raw_spearman.csv,
# then forms the primary transformed matrix Sigma_C = 2 sin(pi * Sigma_A / 6) and its
# eigenvalues. The frozen Sigma_A file (not the recomputation) is what later stages read,
# so results are bit-reproducible. Writes work/sigma/.
################################################################################
source(file.path(Sys.getenv("CTAS_REPO", "."), "config", "ctas_paths.R")); source(file.path(ROOT, "R", "common.R"))
od <- file.path(WORK_DIR, "sigma"); dir.create(od, showWarnings = FALSE, recursive = TRUE)
P <- read_matrix(); SA_re <- sigma_A(P)
SA <- as.matrix(read.csv(file.path(INPUTS, "sigma_A_raw_spearman.csv"), row.names = 1))
d <- max(abs(SA_re - SA)); SC <- sigma_C(SA); dimnames(SC) <- dimnames(SA)
ev <- eigen(SC, symmetric = TRUE, only.values = TRUE)$values
write.csv(SA_re, file.path(od, "sigma_A_recomputed.csv")); write.csv(SC, file.path(od, "sigma_C.csv"))
cat(sprintf("genes %d | max|Sigma_A recomputed - frozen| = %.3g | min eigenvalue Sigma_C = %.17g (%.2f) | off-diagonal range %.3f to %.3f\n",
  nrow(P), d, min(ev), min(ev), min(SC[upper.tri(SC)]), max(SC[upper.tri(SC)])))
if (d > 1e-8) { cat("FAIL: recomputed Sigma_A differs from the frozen input\n"); quit(status = 1) }
cat("PASS\n")
