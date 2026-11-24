#!/usr/bin/env Rscript
################################################################################
# panel_gene_cloud_xaxis_search: build the master table, construct the 5
# candidate x-axes, quantify each, and produce comparison figures.
################################################################################

suppressPackageStartupMessages(library(ggplot2))
have_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
script_dir <- "<CTAS_WORKSPACE>/panel_gene_cloud_xaxis_search"

master <- read.csv("<CTAS_WORKSPACE>/panel_gene_space_figures/gene_space_master_table.csv")
weighted <- read.csv(file.path(script_dir, "gene_panel_weighted_breadth.csv"))

m <- merge(master, weighted, by = "gene")
stopifnot(nrow(m) == nrow(master))  # every FE>=3 union gene must have a weighted-breadth row

# ---- candidate x-axes -------------------------------------------------------
m$raw_panel_count <- m$n_panels_fold3
m$log_raw_panel_count <- log10(m$n_panels_fold3 + 1)
# panel_size_weighted_breadth, effective_dimension_weighted_breadth,
# resolved_signal_mass already present from the join (computed in
# compute_weighted_breadth.R, reusing panel_screen/panel_pi_depth.csv's own
# size/effective_dimension/pi_depth per eligible panel - not recomputed here)

write.csv(m, file.path(script_dir, "gene_cloud_xaxis_master_table.csv"), row.names = FALSE)
cat(sprintf("Master table: %d genes x %d columns\n", nrow(m), ncol(m)))

candidates <- c("raw_panel_count", "log_raw_panel_count", "panel_size_weighted_breadth",
                "effective_dimension_weighted_breadth", "resolved_signal_mass")
candidate_labels <- c(raw_panel_count = "Raw panel count (n_panels_fe3)",
                      log_raw_panel_count = "log10(n_panels_fe3 + 1)",
                      panel_size_weighted_breadth = "Panel-size-weighted breadth",
                      effective_dimension_weighted_breadth = "Effective-dimension-weighted breadth",
                      resolved_signal_mass = "Resolved signal mass (sum -log10(panel pi_depth))")

################################################################################
# QUANTIFICATION
################################################################################
quantify <- function(x, y) {
  n <- length(x)
  rho <- cor(x, y, method = "spearman")
  n_unique <- length(unique(x))
  tab <- sort(table(x), decreasing = TRUE)
  prop_top10 <- sum(head(tab, 10)) / n
  discreteness_score <- 1 - (n_unique / n)  # 0 = fully continuous/unique, 1 = fully collapsed
  data.frame(candidate = NA, spearman_rho = rho, n_unique_x = n_unique, n_genes = n,
             prop_in_top10_positions = prop_top10, discreteness_score = discreteness_score)
}

results <- do.call(rbind, lapply(candidates, function(cand) {
  r <- quantify(m[[cand]], m$magnitude_score)
  r$candidate <- cand
  r
}))

# ---- judgement: fixed absolute cutoffs did not discriminate here (all 5
# candidates have very similar Spearman rho, 0.65-0.66 - none is remotely
# "too circular", and none is anywhere near "too discrete" in absolute terms
# either) - so the judgement is made RELATIVE to raw_panel_count, which is
# the explicit baseline this task is trying to improve on. A candidate is
# "good" only if it reduces discreteness meaningfully (>15% relative
# reduction) while keeping rho within a small tolerance of the baseline
# (not meaningfully more circular); "too circular with y" is reserved for a
# rho increase large enough to suggest redundancy with Score_C; "too
# discrete" for anything WORSE than the raw baseline; anything else is
# "acceptable" (a modest, non-decisive improvement or a wash). --------------
baseline_discreteness <- results$discreteness_score[results$candidate == "raw_panel_count"]
baseline_rho <- results$spearman_rho[results$candidate == "raw_panel_count"]
results$judgement <- sapply(seq_len(nrow(results)), function(i) {
  if (results$candidate[i] %in% c("raw_panel_count", "log_raw_panel_count")) return("baseline (raw breadth)")
  rel_discreteness_change <- (results$discreteness_score[i] - baseline_discreteness) / baseline_discreteness
  rho_gap <- results$spearman_rho[i] - baseline_rho
  if (rho_gap > 0.15) return("too circular with y")
  if (rel_discreteness_change > 0.1) return("too discrete (worse than raw breadth)")
  if (rel_discreteness_change < -0.15) return("good (meaningfully less discrete)")
  return("acceptable (modest improvement)")
})
results <- results[, c("candidate", "spearman_rho", "n_unique_x", "n_genes",
                        "prop_in_top10_positions", "discreteness_score", "judgement")]
write.csv(results, file.path(script_dir, "xaxis_comparison.csv"), row.names = FALSE)
cat("\n=== xaxis_comparison.csv ===\n")
print(results, digits = 3)

################################################################################
# FIGURES - same visual style for every candidate
################################################################################
label_genes <- union(m$gene[order(-m$magnitude_score)][1:15], intersect(c("ENPEP", "UPF3A"), m$gene))
cat(sprintf("\nLabelled genes: %s\n", paste(label_genes, collapse = ", ")))

base_theme <- theme_bw(base_size = 11) +
  theme(legend.position = "right", panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))

save_both <- function(p, name, w = 9, h = 7) {
  ggsave(file.path(script_dir, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(script_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 220)
}

make_plot <- function(xcol, xlab, use_log_scale = FALSE) {
  df <- m
  df$is_v2 <- df$survives_at_q05
  df$is_244 <- df$in_16trait_ctas
  p <- ggplot(df, aes(x = .data[[xcol]], y = magnitude_score)) +
    geom_point(data = df[!df$is_v2, ], color = "grey75", alpha = 0.45, size = 1.5) +
    geom_point(data = df[df$is_v2, ], color = "#B2182B", alpha = 0.65, size = 1.8) +
    geom_point(data = df[df$is_v2 & df$is_244, ], shape = 1, color = "black", alpha = 0.5, size = 2.4, stroke = 0.5)
  if (use_log_scale) p <- p + scale_x_log10()
  p <- p + labs(title = sprintf("Score_C vs. %s", xlab), x = xlab, y = "Score_C (V2 magnitude score)",
                subtitle = "Grey = FE≥3 union; red = V2 q<0.05; black outline = original 244") +
    base_theme
  lab_df <- df[df$gene %in% label_genes, ]
  if (have_ggrepel) {
    p <- p + geom_point(data = lab_df, aes(x = .data[[xcol]], y = magnitude_score), color = "black",
                         shape = 21, size = 2.4, stroke = 0.9) +
      ggrepel::geom_text_repel(data = lab_df, aes(x = .data[[xcol]], y = magnitude_score, label = gene),
                                color = "black", size = 3, fontface = "bold", max.overlaps = 30,
                                segment.size = 0.3, min.segment.length = 0)
  } else {
    p <- p + geom_point(data = lab_df, aes(x = .data[[xcol]], y = magnitude_score), color = "black",
                         shape = 21, size = 2.4, stroke = 0.9) +
      geom_text(data = lab_df, aes(x = .data[[xcol]], y = magnitude_score, label = gene), color = "black",
                size = 3, fontface = "bold", vjust = -0.8)
  }
  p
}

p1 <- make_plot("raw_panel_count", "Raw panel count (n_panels_fe3)", use_log_scale = TRUE)
save_both(p1, "xaxis_raw_panel_count")

p2 <- make_plot("log_raw_panel_count", expression(log[10]*"(n"[panels]*" at FE">=3*" + 1)"))
save_both(p2, "xaxis_log_raw_panel_count")

p3 <- make_plot("panel_size_weighted_breadth", "Panel-size-weighted breadth", use_log_scale = TRUE)
save_both(p3, "xaxis_panel_size_weighted_breadth")

p4 <- make_plot("effective_dimension_weighted_breadth", "Effective-dimension-weighted breadth", use_log_scale = TRUE)
save_both(p4, "xaxis_effective_dimension_weighted_breadth")

p5 <- make_plot("resolved_signal_mass", "Resolved signal mass (sum of -log10(panel pi_depth))", use_log_scale = TRUE)
save_both(p5, "xaxis_resolved_signal_mass")

cat("\nAll 5 candidate figures written (.png + .pdf).\n")
