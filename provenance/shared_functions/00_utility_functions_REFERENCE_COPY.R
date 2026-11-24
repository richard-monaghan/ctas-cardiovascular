#!/usr/bin/env Rscript
################################################################################
# Utility Functions for CVD-CTAS Robustness Checks
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(copula)
  library(Matrix)
})

################################################################################
# ACAT Function
################################################################################

acat <- function(pvalues, weights = NULL) {
  pvalues <- pvalues[!is.na(pvalues)]
  if (length(pvalues) == 0) return(NA)

  pvalues[pvalues == 0] <- 1e-300
  pvalues[pvalues == 1] <- 1 - 1e-10

  if (is.null(weights)) {
    weights <- rep(1/length(pvalues), length(pvalues))
  } else {
    weights <- weights / sum(weights)
  }

  cauchy_stats <- tan((0.5 - pvalues) * pi)
  T_acat <- sum(weights * cauchy_stats)
  p_acat <- 0.5 - atan(T_acat) / pi

  return(p_acat)
}

################################################################################
# Copula Fitting Function
################################################################################

fit_gaussian_copula <- function(pval_matrix, verbose = TRUE) {
  # Clean extreme values
  pval_clean <- pmax(pval_matrix, 1e-300)
  pval_clean <- pmin(pval_clean, 1 - 1e-10)

  # Convert to pseudo-observations
  if (verbose) cat("  Converting to pseudo-observations...\n")
  pseudo_obs <- apply(pval_clean, 2, function(x) {
    rank(x, ties.method = "average", na.last = "keep") / (sum(!is.na(x)) + 1)
  })

  # Estimate correlation
  if (verbose) cat("  Estimating correlation matrix...\n")
  cor_matrix <- cor(qnorm(pseudo_obs), use = "pairwise.complete.obs", method = "spearman")

  # Ensure positive definite
  if (min(eigen(cor_matrix)$values) < 0.01) {
    if (verbose) cat("  Making correlation matrix positive definite...\n")
    cor_matrix <- as.matrix(nearPD(cor_matrix, corr = TRUE)$mat)
  }

  if (verbose) {
    cat(sprintf("  Correlation range: %.3f to %.3f\n",
                min(cor_matrix[lower.tri(cor_matrix)]),
                max(cor_matrix[lower.tri(cor_matrix)])))
  }

  return(list(
    cor_matrix = cor_matrix,
    empirical_quantiles = lapply(1:ncol(pval_clean), function(i) {
      sort(pval_clean[!is.na(pval_clean[, i]), i])
    })
  ))
}

################################################################################
# Transform Uniform to Empirical
################################################################################

transform_to_empirical <- function(uniform_vec, empirical_quantiles, trait_idx) {
  quantile(empirical_quantiles[[trait_idx]], probs = uniform_vec, type = 1)
}

################################################################################
# Run Copula Simulations
################################################################################

run_copula_simulations <- function(copula_obj, empirical_quantiles,
                                   n_genes, n_traits, n_sims = 1000,
                                   verbose = TRUE) {

  null_breadth_001 <- matrix(0, nrow = n_genes, ncol = n_sims)
  null_breadth_0001 <- matrix(0, nrow = n_genes, ncol = n_sims)

  for (i in 1:n_sims) {
    if (verbose && i %% 200 == 0) cat(sprintf("  Simulation %d/%d\n", i, n_sims))

    sim_uniforms <- rCopula(n_genes, copula_obj)
    sim_pvals <- matrix(0, nrow = n_genes, ncol = n_traits)

    for (j in 1:n_traits) {
      sim_pvals[, j] <- transform_to_empirical(sim_uniforms[, j],
                                                empirical_quantiles, j)
    }

    null_breadth_001[, i] <- rowSums(sim_pvals < 0.01)
    null_breadth_0001[, i] <- rowSums(sim_pvals < 0.001)
  }

  return(list(
    null_001 = null_breadth_001,
    null_0001 = null_breadth_0001
  ))
}

################################################################################
# Calculate Enrichment
################################################################################

calculate_enrichment <- function(obs_breadth, null_matrix, max_N = 6) {
  results <- data.frame()

  for (N in 1:max_N) {
    observed <- sum(obs_breadth >= N)
    expected <- mean(colSums(null_matrix >= N))
    expected_sd <- sd(colSums(null_matrix >= N))
    fold_enrichment <- observed / expected
    z_score <- (observed - expected) / expected_sd

    results <- rbind(results, data.frame(
      N = N,
      observed = observed,
      expected = expected,
      expected_sd = expected_sd,
      fold_enrichment = fold_enrichment,
      z_score = z_score
    ))
  }

  return(results)
}

################################################################################
# Pretty Print Results
################################################################################

print_enrichment_table <- function(results, title = "Enrichment Results") {
  cat("\n")
  cat(strrep("=", 80), "\n")
  cat(title, "\n")
  cat(strrep("=", 80), "\n\n")

  for (i in 1:nrow(results)) {
    cat(sprintf("N≥%d: obs=%d, exp=%.1f±%.1f, fold=%.2fx, Z=%.1f\n",
                results$N[i], results$observed[i], results$expected[i],
                results$expected_sd[i], results$fold_enrichment[i],
                results$z_score[i]))
  }
  cat("\n")
}

################################################################################
# Standard Enrichment Plot
################################################################################

plot_enrichment_comparison <- function(results_list, labels, title = "Enrichment Comparison") {

  plot_data <- bind_rows(lapply(seq_along(results_list), function(i) {
    results_list[[i]] %>% mutate(analysis = labels[i])
  }))

  p <- ggplot(plot_data, aes(x = N, y = fold_enrichment, color = analysis)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
    labs(x = "Minimum number of traits (N≥)",
         y = "Fold enrichment",
         title = title,
         color = "Analysis") +
    theme_minimal() +
    theme(legend.position = "bottom")

  return(p)
}

cat("Utility functions loaded successfully.\n")
