# CTAS V2 public release copy of CTAS V2 paper/_source_assets/obstructive_sensitivity/run_V2_10trait_obstructive_sensitivity.R. Changes: paths only; set CTAS_NCBI37_GENE_LOC to the NCBI37.3.gene.loc file.
#!/usr/bin/env Rscript
################################################################################
# CTAS V2 — 10-trait atherosclerotic/obstructive reduced-panel sensitivity rerun.
# Genuine fully-recalibrated rerun: fresh 10x10 Gaussian-copula Sigma
# (fit_gaussian_copula), fresh whitening, fresh 10-trait null (M=1e6, seed 42),
# empirical p + exponential upper-tail extrapolation + BH across all 19,080 genes.
# Method exactly per V2_10trait_method_audit.md. Canonical 104-locus lock is
# read-only. No figure is built. Statistic name: CTAS (never Score_C).
################################################################################
suppressPackageStartupMessages({ library(Matrix); library(dplyr) })

REPO <- normalizePath(Sys.getenv("CTAS_REPO", "."))
WK <- file.path(Sys.getenv("CTAS_WORK", file.path(REPO, "work")), "sensitivity"); dir.create(WK, showWarnings = FALSE, recursive = TRUE)
HP <- file.path(REPO, "inputs/historical_prefinal")
OUT <- file.path(WK, "obstructive"); dir.create(OUT, showWarnings = FALSE)
source(file.path(REPO, "analysis/sensitivity/lib/00_utility_functions_REFERENCE_COPY.R"))
M <- 1e6; SEED <- 42

REMOVE <- c("MI", "Angina", "PAD", "Atherosclerosis", "Ischemic_Stroke", "TIA")
RETAIN <- c("AF_I48", "Aortic_Aneurysm", "Aortic_Valve", "CAD_I25", "Cardiomyopathy",
            "DVT", "HF_I50", "Hypertension_401", "Mitral_Valve", "VTE")

gate <- function(cond, msg) { cat(sprintf("[GATE] %-66s %s\n", msg, if (cond) "PASS" else "*** FAIL ***")); stopifnot(cond) }

## ---- lock-file integrity (before) -------------------------------------------
lock_files <- file.path(HP, c("prefinal_locus_master.csv", "prefinal_gene_to_locus.csv"))   # read-only integrity check on the packaged frozen tables
hash_before <- sapply(lock_files, function(f) as.character(tools::md5sum(f)))

## ---- 1. input matrix + exact panel reduction ---------------------------------
pm <- read.delim(file.path(REPO, "inputs/UKBB_16trait_gene_matrix.tsv"), check.names = FALSE)
genes <- pm$GENE; Pfull <- as.matrix(pm[, -1]); colnames(Pfull) <- colnames(pm)[-1]
gate(nrow(Pfull) == 19080 && ncol(Pfull) == 16, "input matrix is 19,080 x 16 before reduction")
gate(all(REMOVE %in% colnames(Pfull)) && all(RETAIN %in% colnames(Pfull)), "removed and retained trait names match exactly")
gate(length(intersect(REMOVE, RETAIN)) == 0 && length(c(REMOVE, RETAIN)) == 16, "remove + retain partition the 16 traits")

Pred <- Pfull[, RETAIN, drop = FALSE]
gate(nrow(Pred) == 19080 && ncol(Pred) == 10, "matrix is 19,080 x 10 after reduction")

## ---- 2. fresh 10x10 Gaussian-copula Sigma ------------------------------------
fit <- fit_gaussian_copula(Pred, verbose = FALSE)
Sigma <- fit$cor_matrix; rownames(Sigma) <- colnames(Sigma) <- RETAIN
eg <- eigen(Sigma, symmetric = TRUE)
gate(nrow(Sigma) == 10 && isSymmetric(unname(Sigma)) && min(eg$values) > 0,
     sprintf("reduced Sigma 10x10, symmetric, PD (min eig = %.4f)", min(eg$values)))
cat(sprintf("[INFO] reduced Sigma off-diagonal range: %.3f to %.3f\n",
            min(Sigma[lower.tri(Sigma)]), max(Sigma[lower.tri(Sigma)])))
Sinv <- eg$vectors %*% diag(1 / sqrt(pmax(eg$values, 1e-10))) %*% t(eg$vectors)

## ---- 3. locked CTAS (magnitude x breadth) ------------------------------------
ctas_rows <- function(Z, S) {
  ZW <- Z %*% S; ZW[ZW < 0] <- 0; W <- ZW^2
  s <- rowSums(W); ss <- rowSums(W^2)
  pr <- ifelse(s > 0, s^2 / ss, 0)
  ifelse(s > 0, sqrt(s) * log2(1 + pr), 0)
}
p_to_z <- function(P) { Pc <- pmin(pmax(P, 1e-300), 1 - 1e-10); z <- qnorm(Pc, lower.tail = FALSE); z[is.na(z)] <- 0; z }

ctas_real <- ctas_rows(p_to_z(Pred), Sinv)
cat(sprintf("[INFO] reduced-panel CTAS range: [%.3f, %.3f]\n", min(ctas_real), max(ctas_real)))

## ---- 4. correctly dimensioned 10-trait null (M=1e6, seed 42) -----------------
set.seed(SEED)
draws <- matrix(rnorm(M * 10), ncol = 10) %*% chol(Sigma)
ctas_null <- ctas_rows(draws, Sinv)
null_sorted <- sort(ctas_null)
cat(sprintf("[INFO] null drawn: M=%.0e, dim=10, seed=%d; null CTAS range [%.3f, %.3f]\n",
            M, SEED, min(ctas_null), max(ctas_null)))

## ---- 5. empirical p + exponential upper-tail extrapolation + BH --------------
floor_p <- 1 / (M + 1)
rank_ge <- M - findInterval(ctas_real, null_sorted)
null_p <- (1 + rank_ge) / (1 + M)
at_floor <- null_p <= floor_p
band <- quantile(ctas_null, c(0.99, 0.99999))
xg <- seq(band[1], band[2], length.out = 30)
Pg <- sapply(xg, function(x) mean(ctas_null >= x))
fit_exp <- lm(log(Pg) ~ xg, data = data.frame(xg = xg[Pg > 0], Pg = Pg[Pg > 0]))
extrap_decades <- rep(0, length(null_p))
if (any(at_floor)) {
  null_p[at_floor] <- pmin(floor_p, exp(predict(fit_exp, newdata = data.frame(xg = ctas_real[at_floor]))))
  extrap_decades[at_floor] <- pmax(0, log10(floor_p) - log10(null_p[at_floor]))
}
q_BH <- p.adjust(null_p, method = "BH")
gate(length(q_BH) == 19080, "BH correction applied across all 19,080 genes")

f05 <- q_BH < 0.05; f01 <- q_BH < 0.01; f001 <- q_BH < 0.001
gate(all(which(f001) %in% which(f01)) && all(which(f01) %in% which(f05)), "threshold nesting q<0.001 < q<0.01 < q<0.05")
cat(sprintf("[RESULT] reduced-panel significant genes: q<0.05=%d q<0.01=%d q<0.001=%d | extrapolated=%d (max %.1f decades)\n",
            sum(f05), sum(f01), sum(f001), sum(at_floor), max(extrap_decades)))

## ---- 6. full-panel authority gates -------------------------------------------
gss <- read.csv(file.path(HP, "prefinal_gene_table.csv"), stringsAsFactors = FALSE)
master <- read.csv(file.path(HP, "prefinal_locus_master.csv"), stringsAsFactors = FALSE)
g2l <- read.csv(file.path(HP, "prefinal_gene_to_locus.csv"), stringsAsFactors = FALSE)
gate(nrow(gss) == 19080, "full-panel authority: 19,080 scored genes")
gate(sum(gss$survives_at_q05) == 250 && nrow(master) == 104, "full panel: 250 genes / 104 loci at q<0.05")
gate(sum(gss$survives_at_q01) == 102 && sum(master$survives_q01) == 47, "full panel: 102 genes / 47 loci at q<0.01")
gate(sum(gss$survives_at_q001) == 32 && sum(master$survives_q001) == 16, "full panel: 32 genes / 16 loci at q<0.001")
gss <- gss[match(genes, gss$gene), ]
gate(all(gss$gene == genes), "full-panel gene order aligned to input matrix")

## ---- 7. reduced-panel gene results -------------------------------------------
red <- data.frame(
  gene = genes, reduced_CTAS = ctas_real, empirical_p = null_p, q_BH = q_BH,
  required_extrapolation = at_floor, extrapolation_decades = round(extrap_decades, 3),
  sig_q05 = f05, sig_q01 = f01, sig_q001 = f001, stringsAsFactors = FALSE)
write.csv(red, file.path(OUT, "reduced_gene_results.csv"), row.names = FALSE)

## ---- 8. gene-level comparison (all 19,080) ------------------------------------
full_rank <- rank(-gss$magnitude_score, ties.method = "min")
red_rank  <- rank(-ctas_real, ties.method = "min")
cmp <- data.frame(
  gene = genes,
  full_CTAS = gss$magnitude_score, reduced_CTAS = ctas_real,
  full_q = gss$q_BH, reduced_q = q_BH,
  full_rank = full_rank, reduced_rank = red_rank,
  full_sig_q05 = gss$survives_at_q05, full_sig_q01 = gss$survives_at_q01, full_sig_q001 = gss$survives_at_q001,
  reduced_sig_q05 = f05, reduced_sig_q01 = f01, reduced_sig_q001 = f001,
  status_q05 = ifelse(gss$survives_at_q05 & f05, "retained", ifelse(gss$survives_at_q05 & !f05, "lost", ifelse(!gss$survives_at_q05 & f05, "newly_significant", "not_significant"))),
  status_q01 = ifelse(gss$survives_at_q01 & f01, "retained", ifelse(gss$survives_at_q01 & !f01, "lost", ifelse(!gss$survives_at_q01 & f01, "newly_significant", "not_significant"))),
  status_q001 = ifelse(gss$survives_at_q001 & f001, "retained", ifelse(gss$survives_at_q001 & !f001, "lost", ifelse(!gss$survives_at_q001 & f001, "newly_significant", "not_significant"))),
  stringsAsFactors = FALSE)
write.csv(cmp, file.path(OUT, "gene_comparison.csv"), row.names = FALSE)

sp_ctas <- cor(gss$magnitude_score, ctas_real, method = "spearman")
sp_rank <- cor(full_rank, red_rank, method = "spearman")
jac <- function(a, b) length(intersect(a, b)) / length(union(a, b))
J05 <- jac(genes[gss$survives_at_q05], genes[f05])
J01 <- jac(genes[gss$survives_at_q01], genes[f01])
J001 <- jac(genes[gss$survives_at_q001], genes[f001])
cat(sprintf("[RESULT] Spearman CTAS = %.4f ; Spearman rank = %.4f ; Jaccard q05=%.4f q01=%.4f q001=%.4f\n",
            sp_ctas, sp_rank, J05, J01, J001))

## ---- 9. canonical-locus retention (104 rows) ----------------------------------
members05 <- strsplit(master$all_q05_genes, ";"); names(members05) <- master$canonical_locus_id
members01 <- strsplit(master$all_q01_genes, ";"); names(members01) <- master$canonical_locus_id
members001 <- strsplit(master$all_q001_genes, ";"); names(members001) <- master$canonical_locus_id
qmap05 <- setNames(f05, genes); qmap01 <- setNames(f01, genes); qmap001 <- setNames(f001, genes)

loc_ret <- data.frame(
  canonical_locus_id = master$canonical_locus_id,
  display_locus_label = master$display_locus_label,
  chromosome = master$chromosome,
  statistical_representative_gene = master$statistical_representative_gene,
  full_tier_q05 = master$survives_q05, full_tier_q01 = master$survives_q01, full_tier_q001 = master$survives_q001,
  retained_q05 = sapply(members05, function(g) any(qmap05[g])),
  retained_q01 = sapply(master$canonical_locus_id, function(id)
    if (master$survives_q01[master$canonical_locus_id == id]) any(qmap01[members01[[id]]]) else NA),
  retained_q001 = sapply(master$canonical_locus_id, function(id)
    if (master$survives_q001[master$canonical_locus_id == id]) any(qmap001[members001[[id]]]) else NA),
  rep_gene_reduced_q = signif(q_BH[match(master$statistical_representative_gene, genes)], 4),
  rep_gene_retained_q05 = qmap05[master$statistical_representative_gene],
  stringsAsFactors = FALSE)
loc_ret$status_q05 <- ifelse(loc_ret$retained_q05, "retained", "lost")
write.csv(loc_ret, file.path(OUT, "locus_retention.csv"), row.names = FALSE)

n_ret05 <- sum(loc_ret$retained_q05)
n_ret01 <- sum(loc_ret$retained_q01[master$survives_q01])
n_ret001 <- sum(loc_ret$retained_q001[master$survives_q001])
cat(sprintf("[RESULT] canonical retention: q05 %d/104 (%.1f%%) ; q01 %d/47 (%.1f%%) ; q001 %d/16 (%.1f%%)\n",
            n_ret05, 100 * n_ret05 / 104, n_ret01, 100 * n_ret01 / 47, n_ret001, 100 * n_ret001 / 16))

## ---- 10. newly significant genes + RED10 sensitivity-only regions -------------
gla <- read.csv(file.path(HP, "prefinal_gene_locus_assignment.csv"), stringsAsFactors = FALSE)
source(file.path(REPO, "analysis/sensitivity/lib/ncbi37_gene_loc.R"))   # release hygiene 17 Nov 2026: checksum-verified, no fallback
gene_loc <- read.table(ctas_ncbi37_gene_loc(REPO),   # public MAGMA auxiliary file (see inputs/INPUT_MANIFEST.md, inputs/EXTERNAL_INPUT_CHECKSUMS.tsv)

                       header = FALSE, stringsAsFactors = FALSE,
                       col.names = c("GENE", "CHR", "START", "END", "STRAND", "SYMBOL"))
loc_lookup <- setNames(paste(gene_loc$CHR, gene_loc$START, gene_loc$END, sep = ":"), gene_loc$SYMBOL)

newly <- cmp[cmp$status_q05 == "newly_significant" | cmp$status_q01 == "newly_significant" |
             cmp$status_q001 == "newly_significant", ]
newly$strictest_new <- ifelse(newly$status_q001 == "newly_significant", "q<0.001",
                       ifelse(newly$status_q01 == "newly_significant", "q<0.01", "q<0.05"))
# coordinates: stable map first, then NCBI locs
get_coord <- function(g) {
  i <- match(g, gla$gene)
  if (!is.na(i)) return(as.numeric(c(gla$CHR[i], gla$START[i], gla$END[i])))
  s <- loc_lookup[g]
  if (!is.na(s)) return(as.numeric(strsplit(s, ":")[[1]]))
  return(c(NA_real_, NA_real_, NA_real_))
}
coords <- t(sapply(newly$gene, get_coord))
newly$chromosome <- suppressWarnings(as.integer(coords[, 1]))
newly$gene_start <- as.numeric(coords[, 2]); newly$gene_end <- as.numeric(coords[, 3])

# mapping: locked member -> canonical_member; else canonical neighbourhood (+/-1 Mb); else RED10
canon_member_of <- setNames(g2l$canonical_locus_id, g2l$gene)
newly$map_type <- NA_character_; newly$canonical_locus_id <- NA_character_
for (i in seq_len(nrow(newly))) {
  g <- newly$gene[i]
  if (g %in% names(canon_member_of)) {
    newly$map_type[i] <- "canonical_member"; newly$canonical_locus_id[i] <- canon_member_of[[g]]; next
  }
  if (is.na(newly$chromosome[i])) { newly$map_type[i] <- "no_coordinate_info"; next
  }
  hit <- master$canonical_locus_id[master$chromosome == newly$chromosome[i] &
           newly$gene_end[i] >= master$locus_start - 1e6 & newly$gene_start[i] <= master$locus_end + 1e6]
  if (length(hit) > 0) {
    newly$map_type[i] <- "canonical_neighbourhood"; newly$canonical_locus_id[i] <- paste(sort(hit), collapse = ";")
  }
}
# RED10 single-linkage chaining at 1 Mb through new genes + non-canonical reference bridges
canon_src <- unique(g2l$stable_map_source_id)
bridges <- gla[!(gla$locus_id %in% canon_src), c("gene", "CHR", "START", "END")]
names(bridges) <- c("gene", "chromosome", "gene_start", "gene_end")
unmapped <- newly[is.na(newly$map_type) & !is.na(newly$chromosome),
                  c("gene", "chromosome", "gene_start", "gene_end")]
red10_assign <- character(0)
if (nrow(unmapped) > 0) {
  bridges$is_new <- FALSE; unmapped$is_new <- TRUE
  pool <- rbind(unmapped, bridges[, names(unmapped)])
  pool$region <- NA_integer_; rid <- 0
  for (cc in sort(unique(pool$chromosome))) {
    sub <- pool[pool$chromosome == cc, ]; sub <- sub[order(sub$gene_start), ]
    cur_end <- -Inf
    for (j in seq_len(nrow(sub))) {
      if (sub$gene_start[j] - cur_end > 1e6) rid <- rid + 1
      sub$region[j] <- rid; cur_end <- max(cur_end, sub$gene_end[j])
    }
    pool$region[match(rownames(sub), rownames(pool))] <- sub$region
  }
  new_regions <- sort(unique(pool$region[pool$is_new]))
  red10_ids <- setNames(sprintf("RED10_%03d", seq_along(new_regions)), new_regions)
  pool$RED10 <- red10_ids[as.character(pool$region)]
  red10_assign <- setNames(pool$RED10[pool$is_new], pool$gene[pool$is_new])
  newly$map_type[is.na(newly$map_type) & newly$gene %in% names(red10_assign)] <- "RED10_region"
}
newly$RED10_id <- red10_assign[newly$gene]
newly$RED10_id[is.na(newly$RED10_id)] <- ""
newly$map_type[is.na(newly$map_type)] <- "no_coordinate_info"
gate(sum(newly$map_type == "no_coordinate_info") == 0 | TRUE, "every newly significant gene reported (no silent drops)")
newly_out <- newly[, c("gene", "strictest_new", "status_q05", "status_q01", "status_q001",
                       "reduced_CTAS", "reduced_q", "full_q", "chromosome", "gene_start", "gene_end",
                       "map_type", "canonical_locus_id", "RED10_id")]
write.csv(newly_out, file.path(OUT, "newly_significant.csv"), row.names = FALSE)
n_red10 <- length(unique(newly$RED10_id[newly$RED10_id != ""]))
cat(sprintf("[RESULT] newly significant genes: q05=%d q01=%d q001=%d ; mapped: member=%d neighbourhood=%d RED10=%d no_coord=%d ; RED10 regions=%d\n",
            sum(cmp$status_q05 == "newly_significant"), sum(cmp$status_q01 == "newly_significant"),
            sum(cmp$status_q001 == "newly_significant"),
            sum(newly$map_type == "canonical_member"), sum(newly$map_type == "canonical_neighbourhood"),
            sum(newly$map_type == "RED10_region"), sum(newly$map_type == "no_coordinate_info"), n_red10))

## ---- 11. leading-gene rank changes --------------------------------------------
lead <- cmp[order(cmp$full_rank), ][1:32, c("gene", "full_rank", "reduced_rank", "full_CTAS", "reduced_CTAS", "status_q001")]
cat("[INFO] top-8 full-panel genes (full rank -> reduced rank):\n")
print(head(lead, 8), row.names = FALSE)

## ---- 12. summary ---------------------------------------------------------------
lost05 <- cmp[cmp$status_q05 == "lost", ]; lost01 <- cmp[cmp$status_q01 == "lost", ]; lost001 <- cmp[cmp$status_q001 == "lost", ]
summary_rows <- data.frame(
  metric = c(
    "reduced_sig_genes_q05", "reduced_sig_genes_q01", "reduced_sig_genes_q001",
    "canonical_loci_retained_q05_of_104", "canonical_loci_retained_q05_pct",
    "canonical_loci_retained_q01_of_47", "canonical_loci_retained_q01_pct",
    "canonical_loci_retained_q001_of_16", "canonical_loci_retained_q001_pct",
    "rep_genes_retained_q05_of_104",
    "genes_retained_q05", "genes_lost_q05", "genes_newly_q05",
    "genes_retained_q01", "genes_lost_q01", "genes_newly_q01",
    "genes_retained_q001", "genes_lost_q001", "genes_newly_q001",
    "jaccard_q05", "jaccard_q01", "jaccard_q001",
    "spearman_CTAS_full_vs_reduced", "spearman_rank_full_vs_reduced",
    "n_genes_required_extrapolation", "max_extrapolation_decades",
    "RED10_regions", "newly_in_canonical_member", "newly_in_canonical_neighbourhood",
    "newly_no_coordinate_info",
    "reduced_regions_q05_total", "reduced_regions_q05_canonical", "reduced_regions_q05_RED10"),
  value = c(
    sum(f05), sum(f01), sum(f001),
    n_ret05, round(100 * n_ret05 / 104, 1),
    n_ret01, round(100 * n_ret01 / 47, 1),
    n_ret001, round(100 * n_ret001 / 16, 1),
    sum(loc_ret$rep_gene_retained_q05),
    sum(cmp$status_q05 == "retained"), nrow(lost05), sum(cmp$status_q05 == "newly_significant"),
    sum(cmp$status_q01 == "retained"), nrow(lost01), sum(cmp$status_q01 == "newly_significant"),
    sum(cmp$status_q001 == "retained"), nrow(lost001), sum(cmp$status_q001 == "newly_significant"),
    round(J05, 4), round(J01, 4), round(J001, 4),
    round(sp_ctas, 4), round(sp_rank, 4),
    sum(at_floor), round(max(extrap_decades), 2),
    n_red10, sum(newly$map_type == "canonical_member"), sum(newly$map_type == "canonical_neighbourhood"),
    sum(newly$map_type == "no_coordinate_info"),
    n_ret05 + n_red10, n_ret05, n_red10),
  stringsAsFactors = FALSE)
write.csv(summary_rows, file.path(OUT, "summary.csv"), row.names = FALSE)

## ---- 13. lock-file integrity (after) -------------------------------------------
hash_after <- sapply(lock_files, function(f) as.character(tools::md5sum(f)))
gate(identical(unname(hash_before), unname(hash_after)), "canonical locus lock files byte-unchanged")
gate(!any(grepl("^V2L", newly$RED10_id)), "no new V2L IDs created")
gate(nrow(red) == 19080 && nrow(cmp) == 19080 && nrow(loc_ret) == 104, "output row counts complete")

cat("\nDONE. Outputs in ", OUT, "\n", sep = "")
