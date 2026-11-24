#!/usr/bin/env Rscript
################################################################################
# panel_independence_stringency: magnitude-aware independence-weighted score.
#
# Builds on panel_independence_score/'s feasibility-check result (GO on
# calibratability) by fixing Variant B's discrimination weakness: instead of
# a pure "positive participation ratio" (which gives even null genes a
# substantial baseline, since ~half of 16 whitened dimensions are positive by
# chance), this combines a MAGNITUDE term with a SPREAD term, in the same
# functional PATTERN as the project's own established CTAS score
# (CTAS = -log10(ACAT_p) * log2(1+N_sig), i.e. magnitude x log-breadth) -
# reused here, not invented: Score_C = sqrt(total whitened excess) *
# log2(1 + effective_independent_axes).
#
# STILL a single full-16-trait statistic per gene - NO panel search anywhere.
#
# Read-only w.r.t. everything outside panel_independence_stringency/. Writes
# only here.
################################################################################

BREADTH_THRESHOLD <- 0.01  # unused for Score_C itself (continuous, no threshold),
                            # kept only for compatibility with the discrete Score_A
                            # baseline reused from panel_independence_score/
M_NULL <- 1e6
DEFAULT_SEED <- 42

script_dir <- "<CTAS_WORKSPACE>/panel_independence_stringency"
prior_dir <- "<CTAS_WORKSPACE>/panel_independence_score"
screen_dir <- "<CTAS_WORKSPACE>/panel_screen"

Sigma <- as.matrix(read.csv(file.path(screen_dir, "cache", "sigma_matrix.csv"), row.names = 1))
sigma_checksum <- digest::digest(Sigma)
trait_names <- rownames(Sigma)
n_traits <- length(trait_names)
cat(sprintf("Sigma loaded from panel_screen/cache/sigma_matrix.csv, checksum=%s\n", sigma_checksum))
cat("Never refit in this analysis - reused directly from the canonical source.\n\n")

eig <- eigen(Sigma, symmetric = TRUE)
Sigma_inv_sqrt <- eig$vectors %*% diag(1 / sqrt(pmax(eig$values, 1e-10))) %*% t(eig$vectors)

# ---- Magnitude-aware score: Score_C = sqrt(sum(w)) * log2(1 + PR(w)) -------
# w_k = max(0, z_white_k)^2 (only the positive/"excess evidence" direction),
# PR(w) = (sum w)^2 / sum(w^2) in [0,16] (0 or 1 if <=1 nonzero component),
# magnitude = sqrt(sum w) (a z-like scale of total independent-axis excess).
# This directly mirrors CTAS's own (magnitude x log-breadth) structure, reused
# rather than invented, substituting ACAT-significance -> whitened magnitude
# and discrete N_sig -> continuous participation ratio.
score_C <- function(z_vec) {
  z_white <- as.vector(Sigma_inv_sqrt %*% z_vec)
  w <- pmax(0, z_white)^2
  s <- sum(w)
  if (s == 0) return(0)
  pr <- s^2 / sum(w^2)
  sqrt(s) * log2(1 + pr)
}

################################################################################
# REAL GENES
################################################################################
cat("Scoring real genes (single full-16-trait statistic per gene, no search)...\n")
pval_matrix <- read.delim("<CTAS_WORKSPACE>/analytic_null_validation/inputs/UKBB_16trait_gene_matrix.tsv",
                           check.names = FALSE)
gene_names <- pval_matrix$GENE
pval_data <- as.matrix(pval_matrix[, -1])
stopifnot(identical(colnames(pval_data), trait_names))
n_genes <- nrow(pval_data)

pval_clipped <- pmin(pmax(pval_data, 1e-300), 1 - 1e-10)
z_matrix <- qnorm(pval_clipped, lower.tail = FALSE)  # numerically safe (see panel_independence_score note)
z_matrix[is.na(z_matrix)] <- 0

Score_C_real <- sapply(seq_len(n_genes), function(i) score_C(z_matrix[i, ]))
cat(sprintf("  done. Score_C range: [%.3f, %.3f]\n", min(Score_C_real), max(Score_C_real)))

################################################################################
# NULL GENES
################################################################################
cat(sprintf("\nDrawing M_NULL=%.0e null genes ~ N(0,Sigma) (seed=%d, matches panel_screen exactly)...\n",
            M_NULL, DEFAULT_SEED))
L <- chol(Sigma)
set.seed(DEFAULT_SEED)
draws <- matrix(rnorm(M_NULL * n_traits), nrow = M_NULL, ncol = n_traits) %*% L
colnames(draws) <- trait_names

Score_C_null <- apply(draws, 1, score_C)
cat(sprintf("  done. Null Score_C range: [%.3f, %.3f]\n", min(Score_C_null), max(Score_C_null)))

saveRDS(list(Score_C_real = Score_C_real, Score_C_null = Score_C_null, gene_names = gene_names,
             sigma_checksum = sigma_checksum, trait_names = trait_names, M_NULL = M_NULL,
             pval_data = pval_data), file.path(script_dir, "magnitude_scores_cache.rds"))
cat("\nCached magnitude_scores_cache.rds\n")
