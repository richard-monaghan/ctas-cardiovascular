#!/usr/bin/env Rscript
################################################################################
# panel_independence_stringency: empirical p/q for Score_C, stringency sweep,
# locus collapsing, FinnGen/literature overlap, coronary discount, anchor
# retention, and 10-trait nesting.
################################################################################

script_dir <- "<CTAS_WORKSPACE>/panel_independence_stringency"
r <- readRDS(file.path(script_dir, "magnitude_scores_cache.rds"))
M_NULL <- r$M_NULL

# ---- empirical p + exponential-tail extrapolation for the handful above the
# direct-measurement floor (same method as panel_independence_score/) -------
floor_p <- 1 / (M_NULL + 1)
null_sorted <- sort(r$Score_C_null)
rank_ge <- length(null_sorted) - findInterval(r$Score_C_real, null_sorted)
null_p <- (1 + rank_ge) / (1 + M_NULL)

at_floor <- null_p <= floor_p
cat(sprintf("Genes at the direct-measurement floor: %d of %d\n", sum(at_floor), length(null_p)))

band <- quantile(r$Score_C_null, c(0.99, 0.99999))
xg <- seq(band[1], band[2], length.out = 30)
Pg <- sapply(xg, function(x) mean(r$Score_C_null >= x))
keep <- Pg > 0
fit_df <- data.frame(xg = xg[keep], Pg = Pg[keep])
fit_exp <- lm(log(Pg) ~ xg, data = fit_df)

extrap_decades <- rep(0, length(null_p))
if (any(at_floor)) {
  null_p[at_floor] <- pmin(floor_p, exp(predict(fit_exp, newdata = data.frame(xg = r$Score_C_real[at_floor]))))
  extrap_decades[at_floor] <- pmax(0, log10(floor_p) - log10(null_p[at_floor]))
}
q_BH <- p.adjust(null_p, method = "BH")

################################################################################
# Gene metadata: is_core, coronary_composition, locus, FinnGen, literature
################################################################################
ctas16 <- read.delim("<CTAS_WORKSPACE>/manuscript/tables/supplementary_table_2_244_genes.tsv")
ctas16_genes <- ctas16$Gene
is_core <- r$gene_names %in% ctas16_genes

family_map <- read.delim("<CTAS_WORKSPACE>/manuscript/trait_family_coverage/tables/T0_trait_family_mapping.tsv")
coronary_traits <- family_map$trait[family_map$family == "Coronary"]
coronary_idx <- match(coronary_traits, r$trait_names)
crossed_bool <- (r$pval_data < 0.01) & !is.na(r$pval_data)
coronary_composition <- sapply(seq_len(nrow(crossed_bool)), function(i) {
  crossed <- which(crossed_bool[i, ])
  if (length(crossed) == 0) return("none_crossed")
  if (all(crossed %in% coronary_idx)) return("coronary_only")
  if (any(crossed %in% coronary_idx)) return("coronary_plus_other")
  return("no_coronary")
})

# ---- locus mapping (from ST2 = 101 loci; KDF1/ZPR1 flagged as a pre-existing
# discrepancy in the source table - see README) ------------------------------
loci_tab <- read.delim("<CTAS_WORKSPACE>/manuscript/tables/supplementary_table_1_101_loci.tsv")
gene_to_locus <- do.call(rbind, lapply(seq_len(nrow(loci_tab)), function(i) {
  genes <- trimws(strsplit(loci_tab$Genes[i], ",")[[1]])
  data.frame(gene = genes, locus_id = loci_tab$Locus_ID[i], stringsAsFactors = FALSE)
}))
locus_lookup <- setNames(gene_to_locus$locus_id, gene_to_locus$gene)
gene_locus <- locus_lookup[r$gene_names]
missing_locus_but_core <- is_core & is.na(gene_locus)
if (any(missing_locus_but_core)) {
  cat(sprintf("[FLAG FOR HUMAN carried through] %d core genes missing from ST2's gene->locus map (pre-existing source discrepancy): %s\n",
              sum(missing_locus_but_core), paste(r$gene_names[missing_locus_but_core], collapse = ", ")))
}

# ---- FinnGen per-gene replication (reusing the exact established method from
# scripts/generate_finngen_count_summaries.R, recomputed here read-only from
# the same raw FinnGen MAGMA outputs) ---------------------------------------
cat("\nComputing per-gene FinnGen replication (reusing the established method)...\n")
finngen_traits <- c("AF", "Angina", "CAD", "Cardiomyopathy", "HF", "Hypertension", "MI", "Stroke", "TIA", "VTE")
finngen_dir <- "<CTAS_WORKSPACE>/CTAS_FinnGen/results/magma_results"
finngen_list <- list()
for (t in finngen_traits) {
  f <- file.path(finngen_dir, paste0(t, ".genes.out"))
  if (file.exists(f)) finngen_list[[t]] <- read.table(f, header = TRUE, stringsAsFactors = FALSE)
}
finngen_all <- do.call(rbind, finngen_list)
finngen_sig_count <- table(finngen_all$GENE[finngen_all$P < 0.01])
n_finngen_traits <- as.integer(finngen_sig_count[r$gene_names])
n_finngen_traits[is.na(n_finngen_traits)] <- 0
names(n_finngen_traits) <- r$gene_names
cat(sprintf("  FinnGen universe covers %d of %d genes\n", sum(r$gene_names %in% finngen_all$GENE), length(r$gene_names)))

# ---- literature overlap (locus-level, ST11) --------------------------------
lit_tab <- read.delim("<CTAS_WORKSPACE>/manuscript/tables/supplementary_table_10_literature_overlap.tsv")
lit_novel_lookup <- setNames(lit_tab$In_any_prior_study, lit_tab$Locus_ID)

################################################################################
# Assemble gene table
################################################################################
gene_table <- data.frame(
  gene = r$gene_names,
  locus = gene_locus,
  is_core = is_core,
  coronary_composition = coronary_composition,
  crossed_traits = rowSums(crossed_bool),
  magnitude_score = r$Score_C_real,
  null_p = null_p,
  q_BH = q_BH,
  extrapolation_decades_below_band = extrap_decades,
  finngen_n_traits = as.integer(n_finngen_traits[r$gene_names]),
  finngen_support = as.integer(n_finngen_traits[r$gene_names]) >= 1,
  literature_support = lit_novel_lookup[as.character(gene_locus)] == "yes",
  novel_flag = !is.na(gene_locus) & (lit_novel_lookup[as.character(gene_locus)] == "no")
)
gene_table$survives_at_q05 <- gene_table$q_BH < 0.05
gene_table$survives_at_q01 <- gene_table$q_BH < 0.01
gene_table$survives_at_q001 <- gene_table$q_BH < 0.001
gene_table <- gene_table[order(gene_table$q_BH), ]
write.csv(gene_table, file.path(script_dir, "gene_stringency_scores.csv"), row.names = FALSE)

cat(sprintf("\nGenes passing q<0.05: %d (core=%d) | q<0.01: %d (core=%d) | q<0.001: %d (core=%d)\n",
            sum(gene_table$survives_at_q05), sum(gene_table$survives_at_q05 & gene_table$is_core),
            sum(gene_table$survives_at_q01), sum(gene_table$survives_at_q01 & gene_table$is_core),
            sum(gene_table$survives_at_q001), sum(gene_table$survives_at_q001 & gene_table$is_core)))

################################################################################
# DISTILLATION CURVE
################################################################################
q_levels <- c(0.05, 0.01, 0.001)
curve_rows <- list()
for (q in q_levels) {
  surv <- gene_table[gene_table$q_BH < q, ]
  loci_surv <- unique(na.omit(surv$locus))
  curve_rows[[length(curve_rows) + 1]] <- data.frame(
    q_level = q, n_genes = nrow(surv), n_loci = length(loci_surv),
    genes = paste(surv$gene, collapse = ";"),
    loci = paste(sort(loci_surv), collapse = ";")
  )
}
distillation <- do.call(rbind, curve_rows)
write.csv(distillation, file.path(script_dir, "distillation_curve.csv"), row.names = FALSE)
cat("\nDistillation curve:\n")
print(distillation[, c("q_level", "n_genes", "n_loci")], row.names = FALSE)

################################################################################
# IRREDUCIBLE CORE (strictest stringency, q<0.001)
################################################################################
core_strict <- gene_table[gene_table$q_BH < 0.001, ]
core_strict <- core_strict[order(core_strict$q_BH), ]
anchors <- c("MARK3", "LPA", "ATP1B1", "CFDP1")
core_strict$is_anchor <- core_strict$gene %in% anchors
write.csv(core_strict, file.path(script_dir, "irreducible_core.csv"), row.names = FALSE)

cat(sprintf("\nIrreducible core (q<0.001): n=%d genes, %d loci\n", nrow(core_strict), length(unique(na.omit(core_strict$locus)))))
cat("Anchors present:", paste(anchors[anchors %in% core_strict$gene], collapse = ", "), "\n")
cat("Anchors MISSING:", paste(setdiff(anchors, core_strict$gene), collapse = ", "), "(none if blank)\n")

################################################################################
# 10-TRAIT NESTING (ST14)
################################################################################
st14 <- read.delim("<CTAS_WORKSPACE>/manuscript/Supplementary_Table_14_panel_reduction_sensitivity.tsv")
st14_retained_loci <- st14$locus_id[st14$in_10trait_CADsingleton_architecture == "Y"]
core_loci <- unique(na.omit(core_strict$locus))
nested_in_10trait <- core_loci %in% st14_retained_loci
cat(sprintf("\n10-trait nesting: %d of %d irreducible-core loci also retained in the 10-trait-reduced panel\n",
            sum(nested_in_10trait), length(core_loci)))
if (any(!nested_in_10trait)) {
  cat("Core loci NOT retained in 10-trait panel (nesting break - reported, not discarded):",
      paste(core_loci[!nested_in_10trait], collapse = ", "), "\n")
}

################################################################################
# CORONARY DISCOUNT within irreducible core
################################################################################
cat("\nCoronary composition within irreducible core:\n")
print(table(core_strict$coronary_composition))

################################################################################
# EXTERNAL VALIDATION: FinnGen + literature, core vs dropped-from-g_min
################################################################################
cat(sprintf("\nIrreducible core external support: FinnGen=%d/%d, literature=%d/%d\n",
            sum(core_strict$finngen_support), nrow(core_strict),
            sum(core_strict$literature_support, na.rm = TRUE), nrow(core_strict)))

gmin_genes <- read.csv("<CTAS_WORKSPACE>/panel_precision_gate/gene_precision_gate.csv")
gmin_pass05 <- gmin_genes$gene[gmin_genes$pass_q05]
dropped_from_gmin <- setdiff(gmin_pass05, core_strict$gene)
dropped_tab <- gene_table[gene_table$gene %in% dropped_from_gmin, ]
cat(sprintf("\nGenes in old g_min q<0.05 list (n=%d) NOT in the new irreducible core (n=%d dropped):\n",
            length(gmin_pass05), length(dropped_from_gmin)))
cat(sprintf("  Of the dropped genes: FinnGen-supported=%d/%d, literature-supported=%d/%d\n",
            sum(dropped_tab$finngen_support, na.rm = TRUE), nrow(dropped_tab),
            sum(dropped_tab$literature_support, na.rm = TRUE), nrow(dropped_tab)))

saveRDS(list(gene_table = gene_table, distillation = distillation, core_strict = core_strict,
             nested_in_10trait = nested_in_10trait, core_loci = core_loci,
             st14_retained_loci = st14_retained_loci, dropped_from_gmin = dropped_from_gmin,
             dropped_tab = dropped_tab, extrap_decades = extrap_decades),
        file.path(script_dir, "stringency_summary.rds"))
cat("\nSaved all outputs.\n")
