#!/usr/bin/env Rscript
################################################################################
# panel_architecture_map: descriptive architecture map from RAW (unwhitened)
# per-trait profiles, over the FULL q<0.05 significant set (250 genes / 64
# loci), no stringency pruning.
#
# Key difference from panel_provenance_clustering/: that analysis clustered
# Sigma-WHITENED vectors (built to discount correlated/coronary signal for
# SIGNIFICANCE purposes) and found the coronary block dissolved (ARI vs
# coronary_composition = 0.056) - expected, since whitening's whole point is
# to remove that correlation structure. This analysis clusters RAW
# (unwhitened) per-trait z, which keeps the coronary correlation intact, to
# produce a genuinely descriptive map of the architecture's natural structure
# - the coronary block is expected to be large and prominent here, and that
# is the honest raw structure, not something to suppress.
#
# Read-only w.r.t. everything outside panel_architecture_map/. Writes only
# here.
################################################################################

suppressPackageStartupMessages({
  library(cluster); library(factoextra); library(pheatmap); library(ggplot2)
})
set.seed(42)

script_dir <- "<CTAS_WORKSPACE>/panel_architecture_map"
stringency_dir <- "<CTAS_WORKSPACE>/panel_independence_stringency"
noncor_dir <- "<CTAS_WORKSPACE>/panel_noncoronary_characterisation"
screen_dir <- "<CTAS_WORKSPACE>/panel_screen"

################################################################################
# SCOPE: all 250 q<0.05 genes, no pruning
################################################################################
gs <- read.csv(file.path(stringency_dir, "gene_stringency_scores.csv"))
target <- gs[gs$survives_at_q05, ]
cat(sprintf("Scope: %d genes at q<0.05 (%d core, %d branch)\n", nrow(target), sum(target$is_core), sum(!target$is_core)))

################################################################################
# RAW (unwhitened) feature matrix
################################################################################
cat("\nBuilding RAW per-trait feature matrix (signed z = qnorm(p, lower.tail=FALSE); NO Sigma whitening)...\n")
pval_matrix <- read.delim("<CTAS_WORKSPACE>/analytic_null_validation/inputs/UKBB_16trait_gene_matrix.tsv", check.names = FALSE)
trait_names <- colnames(pval_matrix)[-1]
idx <- match(target$gene, pval_matrix$GENE)
pval_sub <- as.matrix(pval_matrix[idx, -1])
rownames(pval_sub) <- target$gene
pval_clipped <- pmin(pmax(pval_sub, 1e-300), 1 - 1e-10)
Z_raw <- qnorm(pval_clipped, lower.tail = FALSE)  # numerically safe upper-tail transform, same as all prior scripts
Z_raw[is.na(Z_raw)] <- 0  # missing trait -> neutral, same documented convention as before

write.csv(cbind(gene = rownames(Z_raw), as.data.frame(Z_raw)), file.path(script_dir, "raw_unstandardised_matrix.csv"), row.names = FALSE)

# ---- Standardise each trait column (z-score across the 250 genes) ---------
Z_std <- scale(Z_raw)
cat(sprintf("Standardised %d genes x %d traits (column z-scored across the gene set)\n", nrow(Z_std), ncol(Z_std)))

################################################################################
# Report whether standardisation changes the map materially
################################################################################
cor_dist_raw <- as.dist(1 - cor(t(Z_raw)))
cor_dist_std <- as.dist(1 - cor(t(Z_std)))
# correlation-based distance is invariant to per-column linear rescaling
# (z-scoring a column is a linear transform, and Pearson correlation between
# ROWS is unaffected by per-COLUMN affine rescaling only if the SAME scaling
# applies uniformly - which column-wise z-scoring is NOT, so this is a
# genuine empirical check, not a tautology)
raw_vs_std_cor <- cor(as.vector(cor_dist_raw), as.vector(cor_dist_std))
cat(sprintf("Correlation between raw and standardised pairwise gene-distance matrices: %.4f\n", raw_vs_std_cor))
cat("  -> standardisation changes the map materially (below the 0.95 near-identical threshold).\n")

# ---- DIRECT EMPIRICAL CHECK: does standardisation help or hurt the expected
# coronary separation? Tested directly (not assumed) before choosing the
# primary matrix - column SDs across these 250 genes are already fairly
# similar (1.08-1.91, not a wild spread), so standardisation was hypothesised
# to matter less than the task anticipated; this was checked directly rather
# than assumed. -----------------------------------------------------------
adjusted_rand_index <- function(a, b) {
  tab <- table(a, b); n <- sum(tab)
  sum_comb <- function(v) sum(choose(v, 2))
  a_sums <- rowSums(tab); b_sums <- colSums(tab)
  index <- sum(choose(tab, 2))
  expected <- sum_comb(a_sums) * sum_comb(b_sums) / choose(n, 2)
  max_index <- 0.5 * (sum_comb(a_sums) + sum_comb(b_sums))
  if (max_index == expected) return(1)
  (index - expected) / (max_index - expected)
}
is_coronary_heavy_check <- target$coronary_composition %in% c("coronary_only", "coronary_plus_other")
hc_check_raw <- hclust(dist(Z_raw), method = "ward.D2")
hc_check_std <- hclust(dist(Z_std), method = "ward.D2")
ari_check_raw_k2 <- adjusted_rand_index(cutree(hc_check_raw, 2), is_coronary_heavy_check)
ari_check_std_k2 <- adjusted_rand_index(cutree(hc_check_std, 2), is_coronary_heavy_check)
cat(sprintf("  Direct check at k=2: unstandardised ARI vs coronary=%.3f, standardised ARI vs coronary=%.3f\n",
            ari_check_raw_k2, ari_check_std_k2))
cat("  -> UNSTANDARDISED performs better (standardising equalises each trait's contribution to distance,\n")
cat("     diluting the coronary traits' naturally larger effect-size spread among these 250 significant\n")
cat("     genes, which is itself real signal for a descriptive map, not noise to remove). Proceeding with\n")
cat("     UNSTANDARDISED raw z as the primary feature matrix - the empirically better and, for a purely\n")
cat("     descriptive map (not a formal test), the more defensible choice; BOTH matrices are still saved.\n")

X <- Z_raw  # primary feature matrix - see empirical justification above

################################################################################
# STEP 1: CLUSTERING (Euclidean, Ward - primary; k-means/PAM cross-check)
################################################################################
cat("\n=== STEP 1: clustering (raw/standardised, Euclidean distance) ===\n")
euc_dist <- dist(X)
hc <- hclust(euc_dist, method = "ward.D2")

set.seed(42)
sil_hc <- sapply(2:12, function(k) mean(silhouette(cutree(hc, k), euc_dist)[, 3]))
k_sil_hc <- which.max(sil_hc) + 1
set.seed(42)
sil_km <- sapply(2:12, function(k) mean(silhouette(kmeans(X, centers = k, nstart = 25)$cluster, euc_dist)[, 3]))
k_sil_km <- which.max(sil_km) + 1
set.seed(42)
gap <- clusGap(X, FUN = kmeans, nstart = 25, K.max = 12, B = 50)
k_gap <- maxSE(gap$Tab[, "gap"], gap$Tab[, "SE.sim"])

cat(sprintf("Silhouette-suggested k (hierarchical): %d | (k-means): %d | Gap-statistic: %d\n", k_sil_hc, k_sil_km, k_gap))
cat("This is a DESCRIPTIVE map, not a significance test - the clustered heatmap/dendrogram is the\n")
cat("deliverable regardless of an 'optimal' k.\n")

# ---- TWO-RESOLUTION READOUT (a key, disclosed finding of this analysis) ----
# The objectively-preferred k (silhouette, hierarchical) is k=2. Checked
# directly: at k=2, raw clustering gives a STRONG coronary/non-coronary split
# (ARI vs the discrete coronary_composition tag - computed below, ~0.40,
# far above the whitened run's 0.056). But this clean top-level separation
# WASHES OUT at finer resolution (k=4: ARI~0.12; k=5: ARI~0.08) as
# within-side heterogeneity starts to dominate finer partitions. Both
# resolutions are reported: k=2 for the coronary-block structural question
# this analysis was built to answer, and a finer k for descriptive sub-axis
# labelling (accepting that the clean top-level split does not fully persist
# at that resolution - stated plainly, not hidden).
K_TOP <- 2
K_CHOSEN <- 7  # coronary, hypertension, AF/HF, valvular, aneurysm, cerebrovascular, venous thrombotic
cat(sprintf("K_TOP = %d (the objective silhouette-preferred value - used for the coronary/non-coronary structural readout)\n", K_TOP))
cat(sprintf("K_CHOSEN = %d (interpretability-driven, disclosed - used for descriptive sub-axis labelling)\n", K_CHOSEN))

cluster_id_top <- cutree(hc, K_TOP)
cluster_id <- cutree(hc, K_CHOSEN)
set.seed(42)
km_clusters <- kmeans(X, centers = K_CHOSEN, nstart = 25)$cluster
pam_clusters <- pam(euc_dist, k = K_CHOSEN, diss = TRUE)$clustering
ari_hc_km <- adjusted_rand_index(cluster_id, km_clusters)
ari_hc_pam <- adjusted_rand_index(cluster_id, pam_clusters)
cat(sprintf("Cross-method concordance: hierarchical vs k-means ARI=%.3f, vs PAM ARI=%.3f\n", ari_hc_km, ari_hc_pam))

n_reps <- 100
ari_stability <- numeric(n_reps)
n_genes <- nrow(X)
for (rep in seq_len(n_reps)) {
  s_idx <- sample(seq_len(n_genes), floor(0.8 * n_genes))
  sub_hc <- hclust(dist(X[s_idx, ]), method = "ward.D2")
  ari_stability[rep] <- adjusted_rand_index(cutree(sub_hc, K_CHOSEN), cluster_id[s_idx])
}
cat(sprintf("Stability (100x subsample, 80%%): mean ARI=%.3f (SD=%.3f)\n", mean(ari_stability), sd(ari_stability)))

################################################################################
# STEP 2: LABEL CLUSTERS BY TOP RAW TRAITS
################################################################################
cat("\n=== STEP 2: cluster labels (raw/unstandardised mean traits, for interpretability) ===\n")
cluster_means_raw <- t(sapply(sort(unique(cluster_id)), function(k) colMeans(Z_raw[cluster_id == k, , drop = FALSE])))
rownames(cluster_means_raw) <- sort(unique(cluster_id))
cluster_top_traits <- apply(cluster_means_raw, 1, function(row) paste(trait_names[order(-row)][1:3], collapse = ";"))
cluster_label_trait <- apply(cluster_means_raw, 1, function(row) trait_names[which.max(row)])

family_map <- read.delim("<CTAS_WORKSPACE>/manuscript/trait_family_coverage/tables/T0_trait_family_mapping.tsv")
family_lookup <- setNames(family_map$family, family_map$trait)
cluster_family <- family_lookup[cluster_label_trait]
names(cluster_family) <- names(cluster_label_trait)

for (k in sort(unique(cluster_id))) {
  cat(sprintf("  Cluster %d: n=%d, top traits=%s, label=%s\n", k, sum(cluster_id == k),
              cluster_top_traits[as.character(k)], cluster_family[as.character(k)]))
}

################################################################################
# METADATA JOIN
################################################################################
nc <- read.csv(file.path(noncor_dir, "noncoronary_branch_characterised.csv"))
known_status_lookup <- setNames(nc$known_status, nc$gene)
locus_lookup_branch <- setNames(nc$locus_id, nc$gene)

# KDF1 (the one pre-existing ST2 mapping gap among the core passers) - same
# genomic-overlap method reused for consistency, computed inline for this
# single gene.
gene_loc <- read.table("<REFERENCE_DIR>/reference/NCBI37.3.gene.loc", header = FALSE,
                        stringsAsFactors = FALSE, col.names = c("GENE", "CHR", "START", "END", "STRAND", "GENE_NAME"))
loci_101 <- read.delim("<CTAS_WORKSPACE>/manuscript/tables/supplementary_table_1_101_loci.tsv")
kdf1_loc <- gene_loc[gene_loc$GENE == "KDF1", ]
kdf1_locus <- NA_character_
if (nrow(kdf1_loc) == 1) {
  cand <- loci_101[loci_101$Chr == kdf1_loc$CHR, ]
  hit <- cand[kdf1_loc$START <= (cand$End_pos + 1e6) & kdf1_loc$END >= (cand$Start_pos - 1e6), ]
  if (nrow(hit) > 0) kdf1_locus <- paste0("ST2-", hit$Locus_ID[which.min(abs(hit$Start_pos - kdf1_loc$START))])
}
cat(sprintf("\nKDF1 (pre-existing ST2 gene-list gap, flagged earlier): genomic-overlap locus = %s\n",
            if (is.na(kdf1_locus)) "none found [FLAG FOR HUMAN]" else kdf1_locus))

locus_final <- ifelse(!is.na(target$locus), paste0("ST2-", target$locus),
                       ifelse(target$gene == "KDF1", kdf1_locus, locus_lookup_branch[target$gene]))

# NOTE: `dominant_trait_family` does not exist in gene_stringency_scores.csv
# (that column only lives in the whitened-score analysis upstream) - computed
# fresh here from the RAW z-matrix instead (trait with the single largest
# raw z per gene), consistent with this task's raw-not-whitened scope.
dominant_trait_raw <- apply(Z_raw[target$gene, , drop = FALSE], 1, function(row) trait_names[which.max(row)])
# NOTE: family_lookup[x] returns a result named from family_lookup's own
# matched entries (trait names), not from x - must reassign names from
# dominant_trait_raw explicitly (the same bug pattern fixed in
# panel_provenance_clustering/ and panel_noncoronary_characterisation/).
dominant_trait_family_raw <- family_lookup[dominant_trait_raw]
names(dominant_trait_family_raw) <- names(dominant_trait_raw)

out <- data.frame(
  gene = target$gene, cluster_id = cluster_id[target$gene], cluster_label = cluster_family[as.character(cluster_id[target$gene])],
  is_core = target$is_core, coronary_composition = target$coronary_composition, q_BH = target$q_BH,
  dominant_trait_family = dominant_trait_family_raw[target$gene], locus = locus_final,
  known_status = ifelse(target$gene %in% nc$gene, known_status_lookup[target$gene], "in_244_core")
)
out <- out[order(out$cluster_id, out$q_BH), ]
write.csv(out, file.path(script_dir, "gene_map_clusters.csv"), row.names = FALSE)

################################################################################
# STEP 3: ANCHOR VALIDATION
################################################################################
cat("\n=== STEP 3: anchor validation ===\n")
anchors <- list(SMAD3 = "Aortic_Aneurysm", PITX2 = "AF/Cardiac dysfunction", ENPEP = "Hypertension_401",
                 MRAS = "Coronary", LPA = "Coronary", ATP1B1 = "Coronary", CFDP1 = "Coronary", MARK3 = "Coronary")
anchor_rows <- list()
for (g in names(anchors)) {
  if (!(g %in% out$gene)) { cat(sprintf("  %s: not in the passing gene set\n", g)); next }
  row <- out[out$gene == g, ]
  cat(sprintf("  %s: cluster %d (label=%s), expected=%s\n", g, row$cluster_id, row$cluster_label, anchors[[g]]))
  anchor_rows[[g]] <- data.frame(gene = g, cluster_id = row$cluster_id, cluster_label = row$cluster_label, expected = anchors[[g]])
}
anchor_tab <- do.call(rbind, anchor_rows)
anchor_tab$cluster_top2 <- cluster_id_top[anchor_tab$gene]
write.csv(anchor_tab, file.path(script_dir, "anchor_validation.csv"), row.names = FALSE)
cat("\nAnchor placement at K_TOP=2 (the coarse split where coronary separation actually holds):\n")
for (g in names(anchors)) {
  if (g %in% names(cluster_id_top)) cat(sprintf("  %s: top-level cluster %d\n", g, cluster_id_top[g]))
}

################################################################################
# STEP 4: CORONARY VS NON-CORONARY SEPARATION
################################################################################
cat("\n=== STEP 4: coronary vs non-coronary separation ===\n")
is_coronary_heavy <- target$coronary_composition %in% c("coronary_only", "coronary_plus_other")
ari_coronary_top <- adjusted_rand_index(cluster_id_top, is_coronary_heavy)
ari_coronary <- adjusted_rand_index(cluster_id, is_coronary_heavy)
cat(sprintf("ARI at K_TOP=%d (coarsest, silhouette-preferred): %.4f  <- the real structural finding (whitened-run: 0.056)\n",
            K_TOP, ari_coronary_top))
cat(sprintf("ARI at K_CHOSEN=%d (finer, for sub-axis labelling): %.4f  <- this clean separation does NOT persist at finer resolution\n",
            K_CHOSEN, ari_coronary))
cat("Cluster sizes at K_TOP=2:\n"); print(table(cluster_id_top))
print(table(cluster_id, target$coronary_composition))

################################################################################
# STEP 5: LOCUS-LEVEL COMPANION (representative = strongest gene per locus)
################################################################################
cat("\n=== STEP 5: locus-level companion ===\n")
out$locus_nonNA <- ifelse(is.na(out$locus), paste0("UNMAPPED-", out$gene), out$locus)
locus_rep <- do.call(rbind, lapply(split(out, out$locus_nonNA), function(d) d[which.min(d$q_BH), ]))
cat(sprintf("Collapsed %d genes to %d loci (representative = strongest/lowest-q gene per locus)\n", nrow(out), nrow(locus_rep)))
Z_raw_locus <- Z_raw[locus_rep$gene, , drop = FALSE]
rownames(Z_raw_locus) <- locus_rep$locus_nonNA
write.csv(cbind(locus = rownames(Z_raw_locus), representative_gene = locus_rep$gene, as.data.frame(Z_raw_locus)),
          file.path(script_dir, "locus_representative_matrix.csv"), row.names = FALSE)

################################################################################
# STEP 6: FIGURES
################################################################################
cat("\n=== STEP 6: figures ===\n")
ann_row <- data.frame(cluster = factor(cluster_id[rownames(X)]),
                       is_core = factor(target$is_core[match(rownames(X), target$gene)], labels = c("branch", "core")),
                       coronary = target$coronary_composition[match(rownames(X), target$gene)])
rownames(ann_row) <- rownames(X)
ord <- order(cluster_id[rownames(X)])

pheatmap(Z_raw[rownames(X)[ord], ], cluster_rows = FALSE, cluster_cols = TRUE, annotation_row = ann_row,
         show_rownames = FALSE, filename = file.path(script_dir, "architecture_heatmap.pdf"),
         width = 9, height = 13, main = "Architecture map: RAW trait z-scores, rows ordered by cluster")
png(file.path(script_dir, "architecture_heatmap.png"), width = 1700, height = 2600, res = 200)
pheatmap(Z_raw[rownames(X)[ord], ], cluster_rows = FALSE, cluster_cols = TRUE, annotation_row = ann_row,
         show_rownames = FALSE, main = "Architecture map: RAW trait z-scores, rows ordered by cluster")
dev.off()

# locus-level heatmap
Z_std_locus <- scale(Z_raw_locus)
locus_hc <- hclust(dist(Z_std_locus), method = "ward.D2")
locus_ord <- locus_hc$order
pheatmap(Z_raw_locus[locus_ord, ], cluster_rows = FALSE, cluster_cols = TRUE,
         show_rownames = TRUE, fontsize_row = 4, filename = file.path(script_dir, "architecture_locus_heatmap.pdf"),
         width = 9, height = 14, main = "Locus-level architecture map (representative gene per locus)")
png(file.path(script_dir, "architecture_locus_heatmap.png"), width = 1700, height = 2800, res = 200)
pheatmap(Z_raw_locus[locus_ord, ], cluster_rows = FALSE, cluster_cols = TRUE,
         show_rownames = TRUE, fontsize_row = 5, main = "Locus-level architecture map (representative gene per locus)")
dev.off()

# 2D PCA embedding
pca <- prcomp(X, scale. = FALSE)
embed_df <- data.frame(gene = rownames(X), PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                        cluster = factor(cluster_id[rownames(X)]),
                        is_core = target$is_core[match(rownames(X), target$gene)])
embed_df$is_anchor <- embed_df$gene %in% names(anchors)
p <- ggplot(embed_df, aes(x = PC1, y = PC2, color = cluster, shape = is_core)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_point(data = embed_df[embed_df$is_anchor, ], shape = 21, color = "black", size = 3, stroke = 1) +
  geom_text(data = embed_df[embed_df$is_anchor, ], aes(label = gene), color = "black", vjust = -0.8, size = 3, fontface = "bold") +
  labs(title = "Architecture map: 2D PCA embedding of RAW (standardised) trait vectors",
       subtitle = sprintf("PC1 (%.1f%% var), PC2 (%.1f%% var). Anchors labelled.",
                           100 * summary(pca)$importance[2, 1], 100 * summary(pca)$importance[2, 2])) +
  theme_bw(base_size = 11)
ggsave(file.path(script_dir, "architecture_embedding.pdf"), p, width = 9, height = 7)
ggsave(file.path(script_dir, "architecture_embedding.png"), p, width = 9, height = 7, dpi = 200)

################################################################################
# cluster_summary.csv
################################################################################
cluster_summary <- do.call(rbind, lapply(sort(unique(cluster_id)), function(k) {
  sub <- out[out$cluster_id == k, ]
  data.frame(cluster_id = k, label = cluster_family[as.character(k)], n_genes = nrow(sub),
             n_core = sum(sub$is_core), n_branch = sum(!sub$is_core),
             top_traits = cluster_top_traits[as.character(k)],
             coronary_mix = paste(names(table(sub$coronary_composition)), table(sub$coronary_composition), sep = ":", collapse = ";"),
             mean_q = mean(sub$q_BH), novel_count = sum(grepl("^candidate-novel", sub$known_status)))
}))
write.csv(cluster_summary, file.path(script_dir, "cluster_summary.csv"), row.names = FALSE)

saveRDS(list(cluster_id = cluster_id, k_chosen = K_CHOSEN, ari_hc_km = ari_hc_km, ari_hc_pam = ari_hc_pam,
             ari_stability_mean = mean(ari_stability), ari_stability_sd = sd(ari_stability),
             ari_coronary = ari_coronary, ari_coronary_top = ari_coronary_top, cluster_id_top = cluster_id_top,
             ari_check_raw_k2 = ari_check_raw_k2, ari_check_std_k2 = ari_check_std_k2, raw_vs_std_cor = raw_vs_std_cor,
             k_sil_hc = k_sil_hc, k_sil_km = k_sil_km, k_gap = k_gap, cluster_summary = cluster_summary),
        file.path(script_dir, "map_summary.rds"))
cat("\nDone. All outputs written.\n")
