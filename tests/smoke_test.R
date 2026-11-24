#!/usr/bin/env Rscript
################################################################################
# smoke_test.R -- fast reproduction check (no 10^9 simulation; ~30 s).
# Run from the repository root:  Rscript tests/smoke_test.R
# Verifies, from the packaged materials only:
#   S01 every file listed in CHECKSUMS.sha256 is present and unchanged
#   S02 Sigma_A recomputed from the gene matrix equals the frozen input; Sigma_C min eigenvalue 0.50
#   S03 observed CTAS, magnitude M and PR recomputed with the canonical code equal results/core
#   S04 empirical P = (1 + null exceedances)/(1 + 10^9) for every non-extrapolated gene
#   S05 Benjamini-Hochberg q recomputed from the packaged P equals q_BH; 223/99/38 genes
#   S06 locus keys recomputed from the packaged locus map equal results/core; 91/46/21 loci
#   S07 locus representatives (highest-CTAS q<0.05 member) equal the published 91-locus table
#   S08 seed stability: identical membership for seeds 42/1/2 (published table)
#   S09 Explorer build tag 4042a25e88 and 19,080-gene / 91-locus counts in its data
# Exit status 0 only if every check passes.
################################################################################
source(file.path(Sys.getenv("CTAS_REPO", "."), "config", "ctas_paths.R")); source(file.path(ROOT, "R", "common.R"))
res <- list(); chk <- function(id, ok, detail) { res[[length(res) + 1]] <<- data.frame(check = id, status = if (isTRUE(ok)) "PASS" else "FAIL", detail = detail); cat(sprintf("%-4s %-5s %s\n", id, if (isTRUE(ok)) "PASS" else "FAIL", detail)) }
sha256 <- function(f) { cmd <- if (nzchar(Sys.which("sha256sum"))) c("sha256sum", "") else c("shasum", "-a 256")
  sub(" .*", "", system(paste(cmd[1], cmd[2], shQuote(f)), intern = TRUE)) }

# S01 checksums
ck <- read.table(file.path(ROOT, "CHECKSUMS.sha256"), header = FALSE, col.names = c("sha", "path"), stringsAsFactors = FALSE, sep = "", quote = "", comment.char = "")
ck$path <- sub("^\\*", "", ck$path); have <- file.exists(file.path(ROOT, ck$path))
obs <- ifelse(have, vapply(file.path(ROOT, ck$path), sha256, ""), NA)
chk("S01", all(have) && all(obs == ck$sha), sprintf("%d files listed; %d missing; %d changed", nrow(ck), sum(!have), sum(have & obs != ck$sha)))

# S02 Sigma
P <- read_matrix(); g <- rownames(P)
SA <- as.matrix(read.csv(file.path(INPUTS, "sigma_A_raw_spearman.csv"), row.names = 1)); SC <- sigma_C(SA); dimnames(SC) <- dimnames(SA)
dSA <- max(abs(sigma_A(P) - SA)); mev <- min(eigen(SC, symmetric = TRUE, only.values = TRUE)$values)
chk("S02", dSA < 1e-8 && sprintf("%.2f", mev) == "0.50", sprintf("max|Sigma_A recomputed - frozen| = %.2g; min eigenvalue Sigma_C = %.6f", dSA, mev))

# S03 observed statistic
G <- read.delim(file.path(ROOT, "results/core/final_C_gene_results_seed42.tsv"), stringsAsFactors = FALSE)
stopifnot(identical(G$gene, g))
pc <- ctas_parts(na0(p_to_z_na(P)), inv_sqrt(SC))
d3 <- c(max(abs(pc$ctas - G$CTAS)), max(abs(pc$mag - G$magnitude_M)), max(abs(pc$pr - G$participation_ratio_PR)))
chk("S03", all(d3 < 1e-9), sprintf("19,080 genes; max abs diff CTAS %.2g, M %.2g, PR %.2g", d3[1], d3[2], d3[3]))

# S04 empirical P from exceedance counts (M = 10^9)
M <- 1e9; ext <- as.logical(G$p_tail_extrapolated); p_re <- (1 + G$null_exceedance_count) / (1 + M)
rel <- max(abs(p_re[!ext] - G$empirical_p[!ext]) / G$empirical_p[!ext])
chk("S04", rel < 1e-12 && sum(ext) == 1, sprintf("max rel diff %.2g over %d genes; %d gene(s) with tail-extrapolated P (%s)", rel, sum(!ext), sum(ext), paste(G$gene[ext], collapse = ",")))

# S05 BH q and gene counts
q_re <- p.adjust(ifelse(ext, G$empirical_p, p_re), "BH"); dq <- max(abs(q_re - G$q_BH) / G$q_BH)
ng <- sapply(c(.05, .01, .001), function(a) sum(q_re < a))
chk("S05", dq < 1e-10 && identical(as.integer(ng), c(223L, 99L, 38L)), sprintf("max rel diff q %.2g; genes q<0.05/0.01/0.001 = %s", dq, paste(ng, collapse = "/")))

# S06 locus keys and counts
LM <- load_locus_map(); lk <- locus_key(g, LM)
nl <- sapply(c(.05, .01, .001), function(a) length(unique(lk[q_re < a])))
chk("S06", identical(lk, G$canonical_locus_id_or_mapping_key) && identical(as.integer(nl), c(91L, 46L, 21L)) && all(grepl("^V2L", unique(lk[q_re < .05]))),
  sprintf("keys identical for all genes: %s; loci q<0.05/0.01/0.001 = %s", identical(lk, G$canonical_locus_id_or_mapping_key), paste(nl, collapse = "/")))

# S07 representatives
L <- read.delim(file.path(ROOT, "results/core/final_C_locus_results_seed42.tsv"), stringsAsFactors = FALSE)
sig <- which(q_re < .05); rep <- tapply(sig, lk[sig], function(ix) g[ix][order(-G$CTAS[ix], G$q_BH[ix], g[ix])][1]); rep <- setNames(as.vector(rep), names(rep))
chk("S07", nrow(L) == 91 && identical(unname(rep[L$canonical_locus_id]), L$best_gene_by_final_C_CTAS) && sum(L$n_q01_genes > 0) == 46 && sum(L$n_q001_genes > 0) == 21,
  sprintf("%d loci; representatives match: %s", nrow(L), identical(unname(rep[L$canonical_locus_id]), L$best_gene_by_final_C_CTAS)))

# S08 seed stability (published)
S <- read.delim(file.path(ROOT, "results/core/seed_stability_C.tsv"))
chk("S08", all(S$identical_to_seed42) && nrow(S) == 18, sprintf("%d seed x level x tier rows; total symmetric difference %d", nrow(S), sum(S$symmetric_difference)))

# S09 Explorer
ex <- file.path(ROOT, "explorer/site"); idx <- readLines(file.path(ex, "index.html"), warn = FALSE)
bt <- trimws(readLines(file.path(ex, "data/build_version.txt"), warn = FALSE)[1])
loci_json <- paste(readLines(file.path(ex, "data/loci.json"), warn = FALSE), collapse = "")
nv2l <- length(unique(regmatches(loci_json, gregexpr("V2L[0-9]{3}", loci_json))[[1]]))
chk("S09", grepl("4042a25e88", bt) && any(grepl("v=4042a25e88", idx)) && nv2l == 91, sprintf("build_version.txt = %s; index.html cache tag present: %s; distinct V2L ids in loci.json = %d", bt, any(grepl("v=4042a25e88", idx)), nv2l))

R <- do.call(rbind, res); dir.create(WORK_DIR, showWarnings = FALSE)
write.table(R, file.path(WORK_DIR, "smoke_test_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\nSMOKE TEST: %s (%d/%d PASS)\n", if (all(R$status == "PASS")) "PASS" else "FAIL", sum(R$status == "PASS"), nrow(R)))
quit(status = if (all(R$status == "PASS")) 0 else 1)
