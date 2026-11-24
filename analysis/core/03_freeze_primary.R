# CTAS V2 public release copy of outputs/final_c_primary_freeze_20261021/code/01_freeze.R (SHA-256 40838559...)
# Changes: paths; relative input-checksum paths; eigenvalue check states the resolved 0.50 expectation;
# former-core retention made optional (superseded reference tables not shipped). All statistics verbatim.
# Usage: Rscript analysis/core/03_freeze_primary.R   (needs work/null_runs/*_M1e+09.rds from step 02)
#!/usr/bin/env Rscript
# Task 1 -- freeze primary estimator-C result. Reads saved 1e9 runs; regenerates no null. Writes only to this folder.
source(file.path(Sys.getenv("CTAS_REPO", "."), "config", "ctas_paths.R")); source(file.path(ROOT, "R", "engine.R"))   # also sources common.R
FZ <- OUT_DIR; dir.create(FZ, showWarnings = FALSE, recursive = TRUE)
wt <- function(x, f) { x <- as.data.frame(x, stringsAsFactors = FALSE); for (j in seq_along(x)) if (is.numeric(x[[j]]) && !is.integer(x[[j]])) x[[j]] <- ifelse(is.na(x[[j]]), NA, format(x[[j]], digits = 17, scientific = FALSE, trim = TRUE) |> (\(s) ifelse(grepl("e", format(x[[j]], digits = 17)), format(x[[j]], digits = 17, scientific = TRUE, trim = TRUE), s))())
  write.table(x, file.path(FZ, f), sep = "\t", quote = FALSE, row.names = FALSE, na = "") }
checks <- list(); chk <- function(id, expected, observed, ok, file) checks[[length(checks) + 1]] <<- data.frame(check_id = id, expected = expected, observed = observed, status = ifelse(isTRUE(ok), "PASS", "FAIL"), evidence_file = file, stringsAsFactors = FALSE)
tn <- c("none", "q<0.05", "q<0.01", "q<0.001"); tier_of <- function(q) tn[1 + (q < .05) + (q < .01) + (q < .001)]
thr <- c(q05 = .05, q01 = .01, q001 = .001)

# ---------- 1. preflight
runfile <- function(est, seed) file.path(W2, sprintf("run_%s_seed%d_M1e+09.rds", est, seed))
RUN <- list(); pre <- list()
for (est in c("C", "A")) for (seed in c(42, 1, 2)) { f <- runfile(est, seed); r <- tryCatch(readRDS(f), error = function(e) NULL)
  ok <- !is.null(r) && isTRUE(r$M == 1e9); pre[[length(pre) + 1]] <- data.frame(estimator = est, seed = seed, file = f, readable = !is.null(r), M = if (is.null(r)) NA else r$M, ok = ok)
  if (!ok) stop("saved run unreadable or M != 1e9: ", f); RUN[[paste(est, seed)]] <- r }
pre <- do.call(rbind, pre)
chk("RUNS_6_READABLE_M1e9", "6 readable, M=1e9", paste(sum(pre$ok), "of 6"), all(pre$ok), "input_checksums.tsv")
inputs <- c(pre$file, file.path(ROOT, c("inputs/UKBB_16trait_gene_matrix.tsv", "inputs/sigma_A_raw_spearman.csv", "inputs/locus_map/gene_neighbourhood_map.csv",
  "inputs/locus_map/stable_to_canonical_locus.csv", "inputs/locus_map/canonical_locus_master.csv", "R/common.R", "R/engine.R")))
sh <- sapply(inputs, function(f) strsplit(system2("shasum", c("-a", "256", shQuote(f)), stdout = TRUE), " ")[[1]][1])
write.table(data.frame(sha256 = sh, bytes = file.size(inputs), path = sub(paste0("^", normalizePath(ROOT), "/"), "", normalizePath(inputs))), file.path(FZ, "input_checksums.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(sessioninfo_line(), capture.output(sessionInfo())), file.path(FZ, "session_info.txt"))

I <- load_inputs(); LM <- load_locus_map(); g <- I$genes; master <- LM$master; lk <- locus_key(g, LM)
chk("GENES_19080_UNIQUE", "19080 present and unique", sprintf("%d rows, %d unique", length(g), length(unique(g))), length(g) == 19080 && !anyDuplicated(g), "final_C_gene_results_seed42.tsv")

# ---------- 2. Sigma_C, eigen, floor, whiteners
SA <- I$Sig$A; SC <- I$Sig$C
SA_re <- sigma_A(I$P); dSA <- max(abs(SA_re - SA))
EG <- eigen(SC, symmetric = TRUE); ev <- EG$values; minev <- min(ev)
floor_active <- any(ev < 1e-10)
W_floor <- inv_sqrt(SC); W_nofloor <- EG$vectors %*% diag(1 / sqrt(ev)) %*% t(EG$vectors)
dW <- max(abs(W_floor - W_nofloor)); sym <- max(abs(W_nofloor - t(W_nofloor))); ident <- max(abs(W_nofloor %*% SC %*% W_nofloor - diag(16)))
pc <- ctas_parts(I$Z, W_nofloor); dCTAS <- max(abs(pc$ctas - I$real$C))
dimnames(SC) <- dimnames(SA); dimnames(W_nofloor) <- dimnames(SA)
write.csv(format(SC, digits = 17), file.path(FZ, "sigma_C.csv"), quote = FALSE)
write.csv(format(W_floor, digits = 17, scientific = FALSE) |> `dimnames<-`(dimnames(SA)), file.path(FZ, "zca_cor_whitener_C.csv"), quote = FALSE)
wt(data.frame(index = seq_along(ev), eigenvalue = ev), "eigenvalues_C.tsv")
wt(data.frame(quantity = c("min_eigenvalue_full_precision", "min_eigenvalue_paper_2dp", "max_eigenvalue", "n_eigenvalues_below_1e-10", "floor_activated", "max_abs_diff_whitener_floor_vs_nofloor",
  "whitener_symmetry_max_abs", "max_abs_W_Sigma_W_minus_I", "max_abs_diff_real_CTAS_engine_vs_nofloor_whitener", "max_abs_diff_sigma_A_file_vs_recomputed", "diag_sigma_C_all_one", "REFERENCE_min_eigenvalue_sigma_A_raw_spearman"),
  value = c(format(minev, digits = 17), sprintf("%.2f", minev), format(max(ev), digits = 17), sum(ev < 1e-10), floor_active, format(dW, digits = 6), format(sym, digits = 6), format(ident, digits = 6), format(dCTAS, digits = 6), format(dSA, digits = 6), all(diag(SC) == 1), format(min(eigen(SA, symmetric = TRUE, only.values = TRUE)$values), digits = 17))), "eigen_floor_check.tsv")
chk("MIN_EIGEN_C_PAPER_FACING_0.50", "positive; 0.50 to 2 dp (Sigma_C; the earlier 0.52 expectation was Sigma_A, resolved 21 Oct 2026)", sprintf("%.17g (%.2f)", minev, minev), minev > 0 && sprintf("%.2f", minev) == "0.50", "eigenvalues_C.tsv")
chk("EIGEN_FLOOR_NOT_ACTIVATED", "no eigenvalue < 1e-10", sprintf("%d below floor", sum(ev < 1e-10)), !floor_active, "eigen_floor_check.tsv")
chk("WHITENERS_FLOOR_VS_NOFLOOR_EQUAL", "max |diff| < 1e-10 (floor inert; only fp rounding expected)", format(dW, digits = 4), dW < 1e-10, "eigen_floor_check.tsv")
chk("SIGMA_A_FILE_MATCHES_RECOMPUTE", "max |diff| < 1e-8", format(dSA, digits = 4), dSA < 1e-8, "eigen_floor_check.tsv")

# ---------- 3. gene result seed 42
q_of <- function(est, seed) { cq <- counts_to_q(I$real[[est]], RUN[[paste(est, seed)]]); list(cq = cq, q = setNames(cq$q, g)) }
C42 <- q_of("C", 42); q42 <- C42$q; cq <- C42$cq
chk0 <- RUN[["C 42"]]$count
gene <- data.frame(gene = g, CTAS = pc$ctas, magnitude_M = pc$mag, participation_ratio_PR = pc$pr, breadth_B_log2_1_plus_PR = log2(1 + pc$pr), null_exceedance_count = chk0,
  empirical_p = cq$p, p_tail_extrapolated = chk0 == 0 & cq$p < cq$floor_p - 1e-15, q_BH = cq$q, q_tier = tier_of(cq$q), canonical_locus_id_or_mapping_key = lk, stringsAsFactors = FALSE)
gene$CHR <- LM$gla$CHR[match(g, LM$gla$gene)]
chk("CTAS_EQ_M_TIMES_B", "CTAS == M * B (max abs diff < 1e-9)", format(max(abs(gene$CTAS - gene$magnitude_M * gene$breadth_B_log2_1_plus_PR)), digits = 4), max(abs(gene$CTAS - gene$magnitude_M * gene$breadth_B_log2_1_plus_PR)) < 1e-9, "final_C_gene_results_seed42.tsv")
wt(gene, "final_C_gene_results_seed42.tsv")
ng <- c(sum(q42 < .05), sum(q42 < .01), sum(q42 < .001)); nl <- c(length(unique(lk[q42 < .05])), length(unique(lk[q42 < .01])), length(unique(lk[q42 < .001])))
chk("GENE_COUNTS_223_99_38", "223/99/38", paste(ng, collapse = "/"), identical(as.integer(ng), c(223L, 99L, 38L)), "final_C_gene_results_seed42.tsv")
chk("LOCUS_COUNTS_91_46_21", "91/46/21", paste(nl, collapse = "/"), identical(as.integer(nl), c(91L, 46L, 21L)), "final_C_locus_results_seed42.tsv")
sigk <- unique(lk[q42 < .05]); outside <- setdiff(sigk, master$canonical_locus_id)
chk("NO_SIG_LOCUS_OUTSIDE_LOCKED_104", "0 outside; 104 locked loci", sprintf("%d outside of %d; locked universe %d", length(outside), length(sigk), length(unique(master$canonical_locus_id))), length(outside) == 0 && nrow(master) == 104, "final_C_locus_results_seed42.tsv")

# ---------- 4. locus result
sig <- which(q42 < .05); rows <- list(); ties <- character()
for (k in sort(unique(lk[sig]))) { ix <- sig[lk[sig] == k]; o <- order(-gene$CTAS[ix], gene$q_BH[ix], gene$gene[ix]); b <- ix[o[1]]
  tie <- length(ix) > 1 && gene$CTAS[ix[o[2]]] == gene$CTAS[b]; if (tie) ties <- c(ties, k)
  m <- master[match(k, master$canonical_locus_id), ]
  mem <- function(th) paste(gene$gene[ix][order(-gene$CTAS[ix])][gene$q_BH[ix][order(-gene$CTAS[ix])] < th], collapse = ";")
  rows[[k]] <- data.frame(canonical_locus_id = k, chromosome = m$chromosome, locus_start = m$locus_start, locus_end = m$locus_end, best_gene_by_final_C_CTAS = gene$gene[b], best_gene_CTAS = gene$CTAS[b],
    best_gene_empirical_p = gene$empirical_p[b], best_gene_q_BH = gene$q_BH[b], best_gene_q_tier = gene$q_tier[b], most_stringent_locus_tier = tier_of(min(gene$q_BH[ix])),
    n_q05_genes = sum(gene$q_BH[ix] < .05), n_q01_genes = sum(gene$q_BH[ix] < .01), n_q001_genes = sum(gene$q_BH[ix] < .001), all_q05_genes = mem(.05), all_q01_genes = mem(.01), all_q001_genes = mem(.001),
    legacy_v1_locus_id = m$legacy_v1_locus_id, lowest_q_gene = gene$gene[ix][which.min(gene$q_BH[ix])], exact_CTAS_tie_at_top = tie, stringsAsFactors = FALSE) }
loc <- do.call(rbind, rows); rownames(loc) <- NULL
chk("LOCUS_TABLE_91_ROWS", "91 q<0.05 loci", as.character(nrow(loc)), nrow(loc) == 91, "final_C_locus_results_seed42.tsv")
chk("LOCUS_TIER_COUNTS_MATCH", "46 / 21 loci with tier<=q<0.01 / q<0.001", paste(sum(loc$n_q01_genes > 0), sum(loc$n_q001_genes > 0), sep = "/"), sum(loc$n_q01_genes > 0) == 46 && sum(loc$n_q001_genes > 0) == 21, "final_C_locus_results_seed42.tsv")
wt(loc, "final_C_locus_results_seed42.tsv")

# ---------- 5. seed stability (C), actual membership
sets <- function(q) list(genes = lapply(thr, function(t) sort(g[q < t])), loci = lapply(thr, function(t) sort(unique(lk[q < t]))))
SC_ <- lapply(c(42, 1, 2), function(s) if (s == 42) q42 else q_of("C", s)$q); names(SC_) <- c(42, 1, 2)
st <- list(); df <- list()
for (s in names(SC_)) { S <- sets(SC_[[s]]); S0 <- sets(SC_[["42"]])
  for (lv in c("genes", "loci")) for (a in names(thr)) { x <- S[[lv]][[a]]; y <- S0[[lv]][[a]]; d1 <- setdiff(x, y); d2 <- setdiff(y, x)
    st[[length(st) + 1]] <- data.frame(estimator = "C", seed = as.integer(s), level = lv, tier = a, n = length(x), n_vs_seed42_only_in_this_seed = length(d1), n_only_in_seed42 = length(d2), symmetric_difference = length(d1) + length(d2), identical_to_seed42 = length(d1) + length(d2) == 0)
    if (length(d1) + length(d2) > 0) df[[length(df) + 1]] <- data.frame(seed = as.integer(s), level = lv, tier = a, only_in_this_seed = paste(d1, collapse = ";"), only_in_seed42 = paste(d2, collapse = ";")) } }
st <- do.call(rbind, st); wt(st, "seed_stability_C.tsv")
dfo <- if (length(df)) do.call(rbind, df) else data.frame(seed = integer(), level = character(), tier = character(), only_in_this_seed = character(), only_in_seed42 = character())
if (!nrow(dfo)) dfo <- data.frame(seed = NA, level = "all", tier = "all", only_in_this_seed = "", only_in_seed42 = "", note = "no differences: gene and locus sets identical across seeds 42, 1, 2 at q<0.05, q<0.01, q<0.001")
wt(dfo, "seed_membership_differences_C.tsv")
chk("SEED_MEMBERSHIP_IDENTICAL_ALL_TIERS", "sets identical, seeds 42/1/2, 3 tiers, genes+loci (symdiff 0)", sprintf("total symmetric difference = %d", sum(st$symmetric_difference)), all(st$identical_to_seed42), "seed_stability_C.tsv")

# ---------- 6. raw-Spearman A sensitivity
SAq <- lapply(c(42, 1, 2), function(s) q_of("A", s)$q); names(SAq) <- c(42, 1, 2)
rs <- list(); add <- function(...) rs[[length(rs) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)
sp <- cor(I$real$A, I$real$C, method = "spearman"); pe <- cor(I$real$A, I$real$C)
add(section = "genomewide_correlation", seed = NA, tier = "", metric = "spearman_CTAS_A_vs_C", value = format(sp, digits = 17)); add(section = "genomewide_correlation", seed = NA, tier = "", metric = "pearson_CTAS_A_vs_C", value = format(pe, digits = 17))
SS <- sets(SAq[["42"]]); for (s in c("42", "1", "2")) { SA_s <- sets(SAq[[s]]); SC_s <- sets(SC_[[s]])
  for (a in names(thr)) { ga <- SA_s$genes[[a]]; gc <- SC_s$genes[[a]]; la <- SA_s$loci[[a]]; lc <- SC_s$loci[[a]]
    for (kv in list(c("genes_A", length(ga)), c("genes_C", length(gc)), c("genes_intersection", length(intersect(ga, gc))), c("genes_jaccard", length(intersect(ga, gc)) / length(union(ga, gc))), c("genes_A_only", length(setdiff(ga, gc))), c("genes_C_only", length(setdiff(gc, ga))),
      c("loci_A", length(la)), c("loci_C", length(lc)), c("loci_intersection", length(intersect(la, lc))), c("loci_jaccard", length(intersect(la, lc)) / length(union(la, lc)))))
      add(section = "A_vs_C_same_seed", seed = as.integer(s), tier = a, metric = kv[1], value = format(as.numeric(kv[2]), digits = 17))
    S0a <- SS$genes[[a]]; add(section = "A_seed_stability_genes_symdiff_vs_seed42", seed = as.integer(s), tier = a, metric = "symmetric_difference", value = as.character(length(setdiff(ga, S0a)) + length(setdiff(S0a, ga))))
    S0l <- SS$loci[[a]]; add(section = "A_seed_stability_loci_symdiff_vs_seed42", seed = as.integer(s), tier = a, metric = "symmetric_difference", value = as.character(length(setdiff(la, S0l)) + length(setdiff(S0l, la)))) } }
add(section = "interpretation", seed = NA, tier = "", metric = "note", value = "Estimator choice has little effect on genome-wide ranking (Spearman rho ~0.99996); estimator C is the prespecified primary analysis and A is a sensitivity comparison only.")
rs <- do.call(rbind, rs); write.table(rs, file.path(FZ, "raw_spearman_sensitivity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
chk("A_VS_C_RANK_CORR_APPROX_0.99996", "~0.99996 (|diff|<5e-6)", format(sp, digits = 9), abs(sp - 0.99996) < 5e-6, "raw_spearman_sensitivity.tsv")

# ---------- 7. former strict core
# RELEASE NOTE: this historical comparison needs the pre-final V2 lock (Sigma_A, M=1e6) gene table and locus tiers,
# which are superseded intermediates and are not redistributed. The frozen output is results/core/former_core_retention.tsv.
REF_GSS <- Sys.getenv("CTAS_FORMER_CORE_GENE_TABLE", ""); REF_MASTER <- Sys.getenv("CTAS_FORMER_CORE_LOCUS_MASTER", "")
if (nzchar(REF_GSS) && nzchar(REF_MASTER)) {
  master <- read.csv(REF_MASTER, stringsAsFactors = FALSE); ROOT_GSS <- REF_GSS
core_l <- master$canonical_locus_id[master$survives_q001]; gss <- read.csv(ROOT_GSS); core_g <- sort(gss$gene[gss$q_BH < .001])
cr <- list(); for (s in names(SC_)) { q <- SC_[[s]]; Lq <- unique(lk[q < .001]); Gq <- g[q < .001]
  cr[[length(cr) + 1]] <- data.frame(seed = as.integer(s), type = "locus", id = core_l, retained_q001 = core_l %in% Lq, final_q = NA)
  cr[[length(cr) + 1]] <- data.frame(seed = as.integer(s), type = "gene", id = core_g, retained_q001 = core_g %in% Gq, final_q = q[core_g]) }
cr <- do.call(rbind, cr); wt(cr, "former_core_retention.tsv")
agg <- aggregate(retained_q001 ~ seed + type, cr, sum)
chk("CORE_LOCI_16_OF_16_ALL_SEEDS", "16/16 for seeds 42,1,2", paste(agg$retained_q001[agg$type == "locus"], collapse = "/"), length(core_l) == 16 && all(agg$retained_q001[agg$type == "locus"] == 16), "former_core_retention.tsv")
chk("CORE_GENES_32_OF_32_ALL_SEEDS", "32/32 for seeds 42,1,2", paste(agg$retained_q001[agg$type == "gene"], collapse = "/"), length(core_g) == 32 && all(agg$retained_q001[agg$type == "gene"] == 32), "former_core_retention.tsv")
} else checks[[length(checks) + 1]] <- data.frame(check_id = "FORMER_CORE_RETENTION", expected = "historical comparison (optional)", observed = "skipped: superseded reference tables not redistributed; see results/core/former_core_retention.tsv", status = "SKIPPED", evidence_file = "results/core/former_core_retention.tsv")
chk("STRICT_LESS_THAN_THRESHOLDS", "strict < used; genes with q exactly at threshold", sprintf("q==0.05: %d, q==0.01: %d, q==0.001: %d", sum(q42 == .05), sum(q42 == .01), sum(q42 == .001)), TRUE, "final_C_gene_results_seed42.tsv")
checks[[length(checks) + 1]] <- data.frame(check_id = "INFO_MIN_EIGEN_SIGMA_A_REFERENCE", expected = "(informational) value that matches the ~0.52 expectation", observed = sprintf("Sigma_A min eigenvalue = %.17g (%.2f); Sigma_C = %.17g (%.2f)", min(eigen(SA, symmetric = TRUE, only.values = TRUE)$values), min(eigen(SA, symmetric = TRUE, only.values = TRUE)$values), minev, minev), status = "INFO", evidence_file = "eigen_floor_check.tsv")
V <- do.call(rbind, checks); write.table(V, file.path(FZ, "validation_checks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE); print(V, right = FALSE)
cat("ties at top of locus:", ifelse(length(ties), paste(ties, collapse = ","), "none"), "\n")
cat("rep (max CTAS) differs from lowest-q gene in", sum(loc$best_gene_by_final_C_CTAS != loc$lowest_q_gene), "loci\n")
cat("min eig", format(minev, digits = 17), " dW", dW, " dCTAS", dCTAS, " dSA", dSA, "\n"); cat("spearman", format(sp, digits = 12), "pearson", format(pe, digits = 12), "\n")
