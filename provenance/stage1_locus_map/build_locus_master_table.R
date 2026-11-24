#!/usr/bin/env Rscript
################################################################################
# panel_locus_cloud: locus-level rebuild of the CTAS V2 gene cloud.
#
# Locus assignment reuses the EXACT method already established in
# panel_noncoronary_characterisation/characterise_branch.R (genomic overlap
# with the original 101 CTAS loci via NCBI37.3.gene.loc coordinates and a
# 1Mb merge rule identical to scripts/collapse_to_loci.R), extended here to
# ALL 1,143 FE>=3 union genes (that prior analysis only covered the 88
# non-244 branch genes) plus the 2 V2-significant genes outside the FE union
# (ENPEP, UPF3A). Gene-level panel-support metrics (panel_size_weighted_
# breadth, effective_dimension_weighted_breadth, resolved_signal_mass,
# magnitude_score/Score_C, q_BH) are REUSED directly from
# panel_gene_cloud_xaxis_search/ - not recomputed.
#
# Read-only w.r.t. everything outside panel_locus_cloud/. Writes only here.
################################################################################

script_dir <- "<CTAS_WORKSPACE>/panel_locus_cloud"
MERGE_DISTANCE <- 1e6  # matches scripts/collapse_to_loci.R exactly, reused not reinvented

################################################################################
# STEP 0: assemble the gene-level table (1,143 FE-union genes + ENPEP/UPF3A)
################################################################################
gene_tab <- read.csv("<CTAS_WORKSPACE>/panel_gene_cloud_xaxis_search/gene_cloud_xaxis_master_table.csv")
gs <- read.csv("<CTAS_WORKSPACE>/panel_independence_stringency/gene_stringency_scores.csv")

extra_genes <- c("ENPEP", "UPF3A")
stopifnot(all(extra_genes %in% gs$gene), !any(extra_genes %in% gene_tab$gene))
extra_rows <- gs[gs$gene %in% extra_genes, c("gene", "magnitude_score", "q_BH", "survives_at_q05",
                                              "survives_at_q01", "survives_at_q001", "coronary_composition")]
extra_rows$in_16trait_ctas <- FALSE
extra_rows$n_panels_fold3 <- 0
extra_rows$raw_panel_count <- 0
extra_rows$log_raw_panel_count <- 0
extra_rows$panel_size_weighted_breadth <- 0
extra_rows$effective_dimension_weighted_breadth <- 0
extra_rows$resolved_signal_mass <- 0

common_cols <- c("gene", "in_16trait_ctas", "magnitude_score", "q_BH", "survives_at_q05",
                  "survives_at_q01", "survives_at_q001", "coronary_composition",
                  "n_panels_fold3", "raw_panel_count", "log_raw_panel_count",
                  "panel_size_weighted_breadth", "effective_dimension_weighted_breadth", "resolved_signal_mass")
gene_all <- rbind(gene_tab[, common_cols], extra_rows[, common_cols])
cat(sprintf("Gene-level table: %d genes (%d FE-union + %d V2-significant-outside-union)\n",
            nrow(gene_all), nrow(gene_tab), length(extra_genes)))

################################################################################
# STEP 1: locus assignment for ALL genes (reusing the exact established method)
################################################################################
gene_loc_path <- "<REFERENCE_DIR>/reference/NCBI37.3.gene.loc"
if (!file.exists(gene_loc_path)) stop("[FLAG FOR HUMAN] Gene-location reference file not found: ", gene_loc_path)
gene_loc <- read.table(gene_loc_path, header = FALSE, stringsAsFactors = FALSE,
                        col.names = c("GENE", "CHR", "START", "END", "STRAND", "GENE_NAME"))

loci_101 <- read.delim("<CTAS_WORKSPACE>/manuscript/tables/supplementary_table_1_101_loci.tsv")

gene_all_loc <- merge(gene_all, gene_loc[, c("GENE", "CHR", "START", "END")], by.x = "gene", by.y = "GENE", all.x = TRUE)
missing_coord <- gene_all_loc$gene[is.na(gene_all_loc$CHR)]
if (length(missing_coord) > 0) {
  cat(sprintf("[FLAG FOR HUMAN] %d genes have no coordinates in the reference file and cannot be locus-mapped: %s\n",
              length(missing_coord), paste(missing_coord, collapse = ", ")))
}
gene_all_loc <- gene_all_loc[!is.na(gene_all_loc$CHR), ]

find_overlapping_locus <- function(chr, start, end) {
  cand <- loci_101[loci_101$Chr == chr, ]
  if (nrow(cand) == 0) return(NA_integer_)
  hit <- cand[start <= (cand$End_pos + MERGE_DISTANCE) & end >= (cand$Start_pos - MERGE_DISTANCE), ]
  if (nrow(hit) == 0) return(NA_integer_)
  hit$Locus_ID[which.min(abs(hit$Start_pos - start))]
}
gene_all_loc$overlapping_locus_id <- mapply(find_overlapping_locus, gene_all_loc$CHR, gene_all_loc$START, gene_all_loc$END)
gene_all_loc$locus_status <- ifelse(is.na(gene_all_loc$overlapping_locus_id), "new", "overlap")

new_genes <- gene_all_loc[gene_all_loc$locus_status == "new", ]
new_genes <- new_genes[order(new_genes$CHR, new_genes$START), ]
if (nrow(new_genes) > 0) {
  new_locus_id <- integer(nrow(new_genes)); new_locus_id[1] <- 1
  last_chr <- new_genes$CHR[1]; last_end <- new_genes$END[1]; cur <- 1
  if (nrow(new_genes) > 1) {
    for (i in 2:nrow(new_genes)) {
      if (new_genes$CHR[i] == last_chr && (new_genes$START[i] - last_end) <= MERGE_DISTANCE) {
        new_locus_id[i] <- cur; last_end <- max(last_end, new_genes$END[i])
      } else {
        cur <- cur + 1; new_locus_id[i] <- cur; last_chr <- new_genes$CHR[i]; last_end <- new_genes$END[i]
      }
    }
  }
  new_genes$new_locus_num <- new_locus_id
}
gene_all_loc$locus_id <- ifelse(gene_all_loc$locus_status == "overlap",
                                 paste0("ST2-", gene_all_loc$overlapping_locus_id),
                                 paste0("NEW-", new_genes$new_locus_num[match(gene_all_loc$gene, new_genes$gene)]))

cat(sprintf("Locus assignment: %d genes overlap an existing 101-locus, %d form %d new loci\n",
            sum(gene_all_loc$locus_status == "overlap"), nrow(new_genes),
            if (nrow(new_genes) > 0) max(new_genes$new_locus_num) else 0))

# ---- cross-validate against the existing (smaller) locus mappings, where
# they overlap, as a sanity check that this re-derivation is consistent -----
existing_map <- read.csv("<CTAS_WORKSPACE>/panel_architecture_map/gene_map_clusters.csv")
common_genes <- intersect(gene_all_loc$gene, existing_map$gene)
existing_locus <- setNames(existing_map$locus, existing_map$gene)[common_genes]
new_locus_for_common <- setNames(gene_all_loc$locus_id, gene_all_loc$gene)[common_genes]
match_rate <- mean(existing_locus == new_locus_for_common, na.rm = TRUE)
cat(sprintf("Cross-validation vs. panel_architecture_map/gene_map_clusters.csv locus assignment (n=%d genes in common): %.1f%% match\n",
            length(common_genes), 100 * match_rate))
mismatch <- common_genes[!is.na(existing_locus) & !is.na(new_locus_for_common) & existing_locus != new_locus_for_common]
if (length(mismatch) > 0) cat(sprintf("  Mismatched genes (%d): %s\n", length(mismatch), paste(mismatch, collapse = ", ")))

write.csv(gene_all_loc, file.path(script_dir, "gene_locus_assignment.csv"), row.names = FALSE)

################################################################################
# STEP 2: LOCUS-LEVEL COLLAPSE
#
# Rule (explicit, documented in README): locus Score_C = MAX gene Score_C in
# the locus; locus q_BH = MIN gene q_BH; all panel-support x-axis candidates
# = MAX across genes in the locus (not SUM), to avoid conflating locus GENE
# COUNT with panel-space support - a large multi-gene locus should not look
# more "supported" merely because it contains more genes; MAX isolates the
# locus's single most-supported representative, mirroring the pre-existing
# "lead gene" convention already used throughout this project's own 101-locus
# definitions (each locus is represented by its lead/most-significant gene).
################################################################################
agg_cols <- c("magnitude_score", "raw_panel_count", "log_raw_panel_count",
              "panel_size_weighted_breadth", "effective_dimension_weighted_breadth", "resolved_signal_mass")
locus_agg <- do.call(rbind, lapply(split(gene_all_loc, gene_all_loc$locus_id), function(d) {
  best_gene <- d$gene[which.max(d$magnitude_score)]
  data.frame(locus_id = d$locus_id[1],
             representative_gene = best_gene,
             n_genes_in_locus = nrow(d),
             locus_Score_C = max(d$magnitude_score),
             locus_q_BH = min(d$q_BH),
             in_original_101 = grepl("^ST2-", d$locus_id[1]),
             contains_v2_hit = any(d$survives_at_q05),
             contains_original_244 = any(d$in_16trait_ctas),
             raw_panel_count = max(d$raw_panel_count),
             log_raw_panel_count = max(d$log_raw_panel_count),
             panel_size_weighted_breadth = max(d$panel_size_weighted_breadth),
             effective_dimension_weighted_breadth = max(d$effective_dimension_weighted_breadth),
             resolved_signal_mass = max(d$resolved_signal_mass),
             genes_in_locus = paste(d$gene, collapse = ";"))
}))
rownames(locus_agg) <- NULL
locus_agg <- locus_agg[order(-locus_agg$locus_Score_C), ]
write.csv(locus_agg, file.path(script_dir, "locus_master_table.csv"), row.names = FALSE)

cat(sprintf("\nLocus-level master table: %d loci (%d overlap the original 101, %d new)\n",
            nrow(locus_agg), sum(locus_agg$in_original_101), sum(!locus_agg$in_original_101)))
cat(sprintf("Loci containing >=1 V2 q<0.05 gene: %d\n", sum(locus_agg$contains_v2_hit)))
cat(sprintf("Loci containing >=1 original-244 gene: %d\n", sum(locus_agg$contains_original_244)))
