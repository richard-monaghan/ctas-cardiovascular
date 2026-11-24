#!/usr/bin/env Rscript
################################################################################
# Deterministic recomputation of the descriptive locus architecture for the
# final 91 CTAS loci with the accepted final ZCA-cor whitener (Task 1). Same
# formulas as _source_assets/scripts/R3_build_architecture.R (family metric,
# RMS sensitivity, raw and Hellinger ordination, clustering diagnostics, seed
# 20260805); only the whitener and the locus set change. Writes source_data/ only.
################################################################################
suppressPackageStartupMessages(library(cluster))
set.seed(20260805)
CT  <- "<CTAS_WORKSPACE>"
T1  <- file.path(CT, "outputs/final_c_primary_freeze_20261021")
OUT <- file.path(CT, "outputs/final_ctas_tables_20261021/source_data")
R3  <- "<CTAS_V2_PAPER_WORKSPACE>/_source_assets/R3_review"

pm <- read.delim(file.path(CT, "analytic_null_validation/inputs/UKBB_16trait_gene_matrix.tsv"), check.names = FALSE)
P <- as.matrix(pm[, -1]); rownames(P) <- pm$GENE
trait_cols <- c("AF_I48","Angina","Aortic_Aneurysm","Aortic_Valve","Atherosclerosis","CAD_I25","Cardiomyopathy","DVT",
                "HF_I50","Hypertension_401","Ischemic_Stroke","MI","Mitral_Valve","PAD","TIA","VTE")
stopifnot(identical(colnames(P), trait_cols))
Z <- qnorm(pmin(pmax(P, 1e-300), 1 - 1e-10), lower.tail = FALSE); Z[is.na(Z)] <- 0      # canonical z, missing -> 0
Wc <- as.matrix(read.csv(file.path(T1, "zca_cor_whitener_C.csv"), row.names = 1, check.names = FALSE))[trait_cols, trait_cols]
gene <- read.delim(file.path(T1, "final_C_gene_results_seed42.tsv"), stringsAsFactors = FALSE)
loc  <- read.delim(file.path(T1, "final_C_locus_results_seed42.tsv"), stringsAsFactors = FALSE)
stopifnot(nrow(loc) == 91)
rownames(gene) <- gene$gene

# ---- whitened vectors for all 223 significant genes; M/PR must reproduce Task 1
sig <- gene$gene[gene$q_BH < 0.05]
ZW <- Z[sig, , drop = FALSE] %*% Wc; colnames(ZW) <- trait_cols
pos <- pmax(ZW, 0); s <- rowSums(pos^2)
M <- sqrt(s); PR <- ifelse(s > 0, s^2 / rowSums(pos^4), 0)
dM <- max(abs(M - gene[sig, "magnitude_M"])); dPR <- max(abs(PR - gene[sig, "participation_ratio_PR"]))
dCT <- max(abs(M * log2(1 + PR) - gene[sig, "CTAS"]))
cat(sprintf("[M/PR check, 223 genes] max|dM|=%.2e max|dPR|=%.2e max|dCTAS|=%.2e\n", dM, dPR, dCT))
stopifnot(dM < 1e-10, dPR < 1e-10, dCT < 1e-10)

families <- list(coronary = c("CAD_I25","MI","Angina"), pressure = "Hypertension_401", arrhythmia = "AF_I48",
  heart_failure = c("HF_I50","Cardiomyopathy"), valve = c("Aortic_Valve","Mitral_Valve"),
  aortic_vascular = c("Aortic_Aneurysm","PAD","Atherosclerosis"), cerebrovascular = c("Ischemic_Stroke","TIA"),
  thrombotic = c("DVT","VTE"))
fam_names <- names(families)
reps <- loc$best_gene_by_final_C_CTAS
signed <- ZW[reps, , drop = FALSE]; posv <- pmax(signed, 0)
fam_agg <- function(r, fun) sapply(families, function(tr) fun(r[tr]))
A_mean <- t(apply(posv, 1, function(r) fam_agg(r, function(x) sum(x) / length(x))))
A_rms  <- t(apply(posv, 1, function(r) fam_agg(r, function(x) sqrt(mean(x^2)))))
shares <- function(A) A / rowSums(A)
S_mean <- shares(A_mean); S_rms <- shares(A_rms); stopifnot(all(abs(rowSums(S_mean) - 1) < 1e-9))
Neff <- 1 / rowSums(S_mean^2); Neff_rms <- 1 / rowSums(S_rms^2)
maxshare <- apply(S_mean, 1, max); largest <- fam_names[apply(S_mean, 1, which.max)]
second <- apply(S_mean, 1, function(x) sort(x, decreasing = TRUE)[2])
lead_label <- ifelse(second >= 0.9 * maxshare, paste(largest, fam_names[apply(S_mean, 1, function(x) order(-x)[2])], sep = " / "), largest)
colnames(A_mean) <- paste0("A_", fam_names); colnames(S_mean) <- paste0("S_", fam_names)

# ---- ranked locus object (representative CTAS desc; ties by empirical p then symbol)
ctas <- loc$best_gene_CTAS; ord <- order(-ctas, loc$best_gene_empirical_p, reps)
tier <- ifelse(loc$n_q001_genes > 0, "q<0.001", ifelse(loc$n_q01_genes > 0, "q<0.01", "q<0.05"))
L <- data.frame(canonical_locus_id = loc$canonical_locus_id, ctas_rank = NA, chromosome = loc$chromosome, locus_start = loc$locus_start, locus_end = loc$locus_end,
  representative_gene = reps, representative_CTAS = ctas, magnitude_M = M[reps], participation_ratio_PR = PR[reps],
  representative_empirical_p = loc$best_gene_empirical_p, representative_q_BH = loc$best_gene_q_BH, locus_q_tier = tier,
  in_q05 = TRUE, in_q01 = loc$n_q01_genes > 0, in_q001 = loc$n_q001_genes > 0, n_q05_genes = loc$n_q05_genes, all_q05_genes = loc$all_q05_genes,
  Neff_family_breadth = Neff, max_family_share = maxshare, largest_share_family = largest, leading_family_display = lead_label,
  A_mean, S_mean, stringsAsFactors = FALSE)
L <- L[ord, ]; L$ctas_rank <- seq_len(nrow(L))
S_o <- S_mean[ord, , drop = FALSE]; Srms_o <- S_rms[ord, , drop = FALSE]; Neff_o <- Neff[ord]; Neffrms_o <- Neff_rms[ord]
prof <- data.frame(canonical_locus_id = L$canonical_locus_id, representative_gene = L$representative_gene,
  setNames(as.data.frame(signed[ord, ]), paste0("w_", trait_cols)), setNames(as.data.frame(posv[ord, ]), paste0("v_", trait_cols)), check.names = FALSE)

# ---- ordination and clustering diagnostics (verbatim recipe; raw shares + Hellinger)
ari <- function(a, b) { t <- table(a, b); n <- sum(t); ai <- rowSums(t); bj <- colSums(t)
  idx <- sum(choose(t, 2)); e <- sum(choose(ai, 2)) * sum(choose(bj, 2)) / choose(n, 2); m <- (sum(choose(ai, 2)) + sum(choose(bj, 2))) / 2; (idx - e) / (m - e) }
qlev <- ifelse(L$in_q001, 3, ifelse(L$in_q01, 2, 1))
ordination <- function(Mx, seed = 20260805) {
  set.seed(seed); pca <- prcomp(Mx, center = TRUE, scale. = FALSE); ve <- pca$sdev^2 / sum(pca$sdev^2)
  D <- dist(Mx); hc <- hclust(D, method = "ward.D2"); rows <- list()
  for (k in 2:10) { ch <- cutree(hc, k); ck <- kmeans(Mx, centers = k, nstart = 50, iter.max = 100)$cluster; cp <- pam(D, k = k, diss = TRUE)$clustering
    rows[[length(rows) + 1]] <- data.frame(k = k, sil_hclust = mean(silhouette(ch, D)[, 3]), sil_kmeans = mean(silhouette(ck, D)[, 3]), sil_pam = mean(silhouette(cp, D)[, 3]),
      ARI_hc_km = ari(ch, ck), ARI_hc_pam = ari(ch, cp), ARI_km_pam = ari(ck, cp)) }
  sil <- do.call(rbind, rows)
  sil$boot_stability_pam <- sapply(2:10, function(k) { full <- pam(D, k = k, diss = TRUE)$clustering
    mean(replicate(200, { idx <- sort(sample(nrow(Mx), round(0.8 * nrow(Mx)))); ari(full[idx], pam(dist(Mx[idx, ]), k = k, diss = TRUE)$clustering) })) })
  drivers <- data.frame(variable = c("CTAS", "Neff_family_breadth", paste0("S_", fam_names), "q_level_ordinal"),
    cor_PC1 = c(cor(L$representative_CTAS, pca$x[, 1]), cor(Neff_o, pca$x[, 1]), sapply(1:8, function(j) cor(S_o[, j], pca$x[, 1])), cor(qlev, pca$x[, 1])),
    cor_PC2 = c(cor(L$representative_CTAS, pca$x[, 2]), cor(Neff_o, pca$x[, 2]), sapply(1:8, function(j) cor(S_o[, j], pca$x[, 2])), cor(qlev, pca$x[, 2])))
  ms <- max(sil$sil_pam); st <- sil$boot_stability_pam[which.max(sil$sil_pam)]
  list(pca = pca, ve = ve, sil = sil, drivers = drivers, max_sil = ms, stab = st, verdict = if (ms < 0.5 || st < 0.6) "OVERLAPPING_CONTINUUM" else "POSSIBLE_DISCRETE_STRUCTURE_FLAG")
}
raw <- ordination(S_o); hel <- ordination(sqrt(S_o))
dist_cor_rh <- cor(as.vector(dist(S_o)), as.vector(dist(sqrt(S_o))))
proc <- function(X, Y) { Xc <- scale(X, scale = FALSE); Yc <- scale(Y, scale = FALSE); s <- svd(t(Xc) %*% Yc); sum(s$d) / sqrt(sum(Xc^2) * sum(Yc^2)) }
pca_df <- data.frame(canonical_locus_id = L$canonical_locus_id, ctas_rank = L$ctas_rank, representative_gene = L$representative_gene, CTAS = L$representative_CTAS,
  Neff_family_breadth = Neff_o, locus_q_tier = L$locus_q_tier, setNames(as.data.frame(hel$pca$x[, 1:8]), paste0("hellinger_PC", 1:8)),
  setNames(as.data.frame(raw$pca$x[, 1:8]), paste0("raw_PC", 1:8)))
tag <- function(o, t) rbind(
  data.frame(metric = paste0("PCA_var_explained_", t), key = paste0("PC", 1:8), value = o$ve[1:8]),
  data.frame(metric = paste0("silhouette_pam_", t), key = paste0("k=", o$sil$k), value = o$sil$sil_pam),
  data.frame(metric = paste0("silhouette_hclust_", t), key = paste0("k=", o$sil$k), value = o$sil$sil_hclust),
  data.frame(metric = paste0("silhouette_kmeans_", t), key = paste0("k=", o$sil$k), value = o$sil$sil_kmeans),
  data.frame(metric = paste0("ARI_hclust_kmeans_", t), key = paste0("k=", o$sil$k), value = o$sil$ARI_hc_km),
  data.frame(metric = paste0("ARI_hclust_pam_", t), key = paste0("k=", o$sil$k), value = o$sil$ARI_hc_pam),
  data.frame(metric = paste0("ARI_kmeans_pam_", t), key = paste0("k=", o$sil$k), value = o$sil$ARI_km_pam),
  data.frame(metric = paste0("boot_stability_pam_", t), key = paste0("k=", o$sil$k), value = o$sil$boot_stability_pam),
  data.frame(metric = paste0("driver_cor_PC1_", t), key = o$drivers$variable, value = o$drivers$cor_PC1),
  data.frame(metric = paste0("driver_cor_PC2_", t), key = o$drivers$variable, value = o$drivers$cor_PC2))
struct <- rbind(tag(hel, "hellinger_final"), tag(raw, "raw_sensitivity"),
  data.frame(metric = "raw_vs_hellinger", key = c("pairwise_distance_correlation", "procrustes_concordance", "PC1_abscor", "PC2_abscor", "verdict_agreement"),
    value = c(dist_cor_rh, proc(raw$pca$x[, 1:2], hel$pca$x[, 1:2]), abs(cor(raw$pca$x[, 1], hel$pca$x[, 1])), abs(cor(raw$pca$x[, 2], hel$pca$x[, 2])), as.numeric(raw$verdict == hel$verdict))),
  data.frame(metric = "verdict_hellinger_final", key = hel$verdict, value = NA), data.frame(metric = "verdict_raw_sensitivity", key = raw$verdict, value = NA))
prof_cor <- sapply(seq_len(nrow(S_o)), function(i) if (sd(S_o[i, ]) == 0 || sd(Srms_o[i, ]) == 0) NA else cor(S_o[i, ], Srms_o[i, ]))
pca_rms <- prcomp(Srms_o, center = TRUE, scale. = FALSE)
sens <- data.frame(check = c("median_per_locus_profile_correlation", "min_per_locus_profile_correlation", "effective_breadth_correlation", "pairwise_distance_correlation", "PCA_PC1_concordance_abscor", "PCA_PC2_concordance_abscor"),
  value = c(median(prof_cor, na.rm = TRUE), min(prof_cor, na.rm = TRUE), cor(Neff_o, Neffrms_o), cor(as.vector(dist(S_o)), as.vector(dist(Srms_o))), abs(cor(raw$pca$x[, 1], pca_rms$x[, 1])), abs(cor(raw$pca$x[, 2], pca_rms$x[, 2]))))
sens$verdict <- ""; sens <- rbind(sens, data.frame(check = "verdict", value = NA, verdict = if (median(prof_cor, na.rm = TRUE) > 0.9 && sens$value[4] > 0.9 && sens$value[3] > 0.9) "STABLE" else "MATERIAL_INSTABILITY_FLAG"))

# ---- raw-P descriptive family coverage (method-independent; P<0.01)
pr <- P[L$representative_gene, , drop = FALSE]
famflag <- sapply(families, function(tr) apply(pr[, tr, drop = FALSE] < 0.01, 1, function(x) any(x, na.rm = TRUE)))
nfam <- rowSums(famflag)
cov <- data.frame(canonical_locus_id = L$canonical_locus_id, representative_gene = L$representative_gene, in_q01 = L$in_q01, in_q001 = L$in_q001,
  setNames(as.data.frame(pr), paste0("p_", trait_cols)), setNames(as.data.frame(famflag + 0), paste0("fam_", fam_names)), n_families_p01 = nfam,
  family_count_category = ifelse(nfam >= 4, "4+", as.character(nfam)), check.names = FALSE)
cats <- c("0", "1", "2", "3", "4+")
covsum <- do.call(rbind, lapply(list(`q<0.05` = rep(TRUE, nrow(L)), `q<0.01` = L$in_q01, `q<0.001` = L$in_q001), function(ix) {
  tb <- table(factor(cov$family_count_category[ix], levels = cats)); data.frame(n_loci = sum(ix), family_count = cats, n = as.integer(tb), pct = 100 * as.integer(tb) / sum(ix)) }))
covsum$locus_set <- rep(c("q<0.05", "q<0.01", "q<0.001"), each = 5)

# ---- comparison with the earlier (raw-Spearman-whitened) derived summaries
old <- read.csv(file.path(R3, "R3_locus_ranked_architecture.csv"), stringsAsFactors = FALSE)
oldcov <- read.csv(file.path(R3, "R3_descriptive_family_coverage.csv"), stringsAsFactors = FALSE)
oi <- match(L$canonical_locus_id, old$canonical_locus_id)
stopifnot(!anyNA(oi), all(old$statistical_representative_gene[oi] == L$representative_gene))
cmp <- data.frame(canonical_locus_id = L$canonical_locus_id, representative_gene = L$representative_gene,
  representative_matches_earlier = old$statistical_representative_gene[oi] == L$representative_gene,
  CTAS_final = L$representative_CTAS, CTAS_earlier = old$CTAS[oi], Neff_final = L$Neff_family_breadth, Neff_earlier = old$Neff_family_breadth[oi],
  delta_Neff = L$Neff_family_breadth - old$Neff_family_breadth[oi], leading_final = L$largest_share_family, leading_earlier = old$largest_share_family[oi],
  leading_family_changed = L$largest_share_family != old$largest_share_family[oi],
  n_families_p01_final = nfam, n_families_p01_earlier = oldcov$n_families[match(L$canonical_locus_id, oldcov$canonical_locus_id)])
cmp$raw_coverage_identical <- cmp$n_families_p01_final == cmp$n_families_p01_earlier
cat(sprintf("[compare] leading-family changes %d/91; max|dNeff| %.4f; median|dNeff| %.4f; raw-coverage identical %d/91\n",
  sum(cmp$leading_family_changed), max(abs(cmp$delta_Neff)), median(abs(cmp$delta_Neff)), sum(cmp$raw_coverage_identical)))
cat(sprintf("[ordination] Hellinger best sil %.3f stab %.3f -> %s | raw -> %s\n", hel$max_sil, hel$stab, hel$verdict, raw$verdict))

wr <- function(d, f) write.table(d, file.path(OUT, f), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
wr(L, "architecture_locus_ranked.tsv"); wr(prof, "architecture_trait_profiles.tsv"); wr(pca_df, "architecture_ordination.tsv")
wr(struct, "architecture_structure_diagnostics.tsv"); wr(sens, "architecture_family_metric_sensitivity.tsv")
wr(cov, "architecture_raw_family_coverage.tsv"); wr(covsum, "architecture_raw_family_coverage_summary.tsv"); wr(cmp, "architecture_vs_earlier_summaries.tsv")
wr(data.frame(check = c("max_abs_dM_223", "max_abs_dPR_223", "max_abs_dCTAS_223"), value = c(dM, dPR, dCT)), "architecture_MPR_reproduction.tsv")
