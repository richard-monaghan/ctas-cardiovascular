# CTAS V2 public release copy of panel_v2_effective_n_sensitivity/run_v2_effective_n.R. Changes: paths only (needs the reference null from 00_reference_null_1e6.R).
#!/usr/bin/env Rscript
################################################################################
# V2 CTAS effective-N sensitivity — direct bridge through the LOCKED V2 pipeline.
#
# Verification analysis only. The ONLY upstream change is the real gene x trait
# p-value matrix (flat-N -> trait-specific effective-N). Every V2 parameter, the
# locked Sigma, the locked 1e6-draw null (seed 42) and the tail/BH procedures are
# reused verbatim. See PIPELINE_PROVENANCE.md.
#
# Writes ONLY into panel_v2_effective_n_sensitivity/. Reads everything else
# read-only. Does not touch any locked object.
################################################################################
suppressPackageStartupMessages({ library(digest) })
options(stringsAsFactors = FALSE)

REPO <- normalizePath(Sys.getenv("CTAS_REPO", "."))
WK <- file.path(Sys.getenv("CTAS_WORK", file.path(REPO, "work")), "sensitivity"); dir.create(WK, showWarnings = FALSE, recursive = TRUE)
HP <- file.path(REPO, "inputs/historical_prefinal")
OUT  <- file.path(WK, "effective_n")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

FLAT_TSV  <- file.path(REPO, "inputs/UKBB_16trait_gene_matrix.tsv")
EFF_TSV   <- file.path(REPO, "inputs/effective_n/gene_matrix_effective_n.tsv")
SIGMA_CSV <- file.path(REPO, "inputs/sigma_A_raw_spearman.csv")
CACHE_RDS <- file.path(WK, "magnitude_scores_cache.rds")
SEALED_GSS<- file.path(HP, "prefinal_gene_table.csv")
GLA_CSV   <- file.path(HP, "prefinal_gene_locus_assignment.csv")
MASTER_CSV<- file.path(HP, "prefinal_locus_master.csv")
G2L_CSV   <- file.path(HP, "prefinal_gene_to_locus.csv")
ENTAB_CSV <- file.path(REPO, "inputs/effective_n/effective_n_table.csv")

msg <- function(...) cat(sprintf(...), "\n")
GATES <- list()
setgate <- function(name, pass, detail="") { GATES[[name]] <<- list(pass=pass, detail=detail)
  msg("[GATE] %-42s %s  %s", name, if (pass) "PASS" else "*** FAIL ***", detail) }

################################################################################
# 0. Load locked objects
################################################################################
Sigma <- as.matrix(read.csv(SIGMA_CSV, row.names = 1))
sigma_obj_checksum <- digest(Sigma)
trait_names <- rownames(Sigma)
n_traits <- length(trait_names)
msg("Sigma loaded: %dx%d  object-checksum=%s", nrow(Sigma), ncol(Sigma), sigma_obj_checksum)

eig <- eigen(Sigma, symmetric = TRUE)
Sigma_inv_sqrt <- eig$vectors %*% diag(1 / sqrt(pmax(eig$values, 1e-10))) %*% t(eig$vectors)

r <- readRDS(CACHE_RDS)                 # locked cache: Score_C_null (1e6, seed42), gene_names, pval_data (flat-N)
M_NULL <- r$M_NULL
Score_C_null <- r$Score_C_null
gene_names <- r$gene_names              # canonical (flat-N) gene order
stopifnot(identical(r$trait_names, trait_names))
msg("Locked null reused: M_NULL=%.0e  |Score_C_null|=%d  cache sigma_checksum=%s (match=%s)",
    M_NULL, length(Score_C_null), r$sigma_checksum, identical(r$sigma_checksum, sigma_obj_checksum))

################################################################################
# Frozen pipeline functions (verbatim from the locked scripts)
################################################################################
p_to_z <- function(pval_data) {
  pc <- pmin(pmax(pval_data, 1e-300), 1 - 1e-10)
  z <- qnorm(pc, lower.tail = FALSE)
  z[is.na(z)] <- 0
  z
}
score_C <- function(z_vec) {
  z_white <- as.vector(Sigma_inv_sqrt %*% z_vec)
  w <- pmax(0, z_white)^2
  s <- sum(w)
  if (s == 0) return(0)
  pr <- s^2 / sum(w^2)
  sqrt(s) * log2(1 + pr)
}
score_all <- function(pval_data) {
  z <- p_to_z(pval_data)
  vapply(seq_len(nrow(z)), function(i) score_C(z[i, ]), numeric(1))
}
# empirical p + exponential tail extrapolation + BH  (verbatim from stringency_sweep.R)
emp_p_q <- function(Score_C_real, Score_C_null, M_NULL) {
  floor_p <- 1 / (M_NULL + 1)
  null_sorted <- sort(Score_C_null)
  rank_ge <- length(null_sorted) - findInterval(Score_C_real, null_sorted)
  null_p <- (1 + rank_ge) / (1 + M_NULL)
  at_floor <- null_p <= floor_p
  band <- quantile(Score_C_null, c(0.99, 0.99999))
  xg <- seq(band[1], band[2], length.out = 30)
  Pg <- sapply(xg, function(x) mean(Score_C_null >= x))
  keep <- Pg > 0
  fit_df <- data.frame(xg = xg[keep], Pg = Pg[keep])
  fit_exp <- lm(log(Pg) ~ xg, data = fit_df)
  extrap_decades <- rep(0, length(null_p))
  if (any(at_floor)) {
    null_p[at_floor] <- pmin(floor_p, exp(predict(fit_exp, newdata = data.frame(xg = Score_C_real[at_floor]))))
    extrap_decades[at_floor] <- pmax(0, log10(floor_p) - log10(null_p[at_floor]))
  }
  q_BH <- p.adjust(null_p, method = "BH")
  list(null_p = null_p, q_BH = q_BH, extrap = extrap_decades, at_floor = at_floor)
}

################################################################################
# 1. Load and align both matrices to the canonical (flat-N) gene order
################################################################################
flat_df <- read.delim(FLAT_TSV, check.names = FALSE)
eff_df  <- read.delim(EFF_TSV,  check.names = FALSE)
stopifnot(identical(colnames(flat_df), colnames(eff_df)))
stopifnot(identical(colnames(flat_df)[-1], trait_names))

flat_mat <- as.matrix(flat_df[, -1]); rownames(flat_mat) <- flat_df$GENE
eff_mat  <- as.matrix(eff_df[, -1]);  rownames(eff_mat)  <- eff_df$GENE
# align both to the locked gene order (gene_names) and trait order
stopifnot(setequal(rownames(flat_mat), gene_names), setequal(rownames(eff_mat), gene_names))
flat_al <- flat_mat[gene_names, trait_names]
eff_al  <- eff_mat[gene_names,  trait_names]
# baseline sanity: the aligned flat matrix must equal the locked cache pval_data
stopifnot(identical(dim(flat_al), dim(r$pval_data)))
flat_vs_cache_max <- max(abs(flat_al - r$pval_data), na.rm = TRUE)
flat_vs_cache_na  <- sum(is.na(flat_al) != is.na(r$pval_data))
setgate("flat_matrix_matches_cache", flat_vs_cache_max == 0 && flat_vs_cache_na == 0,
        sprintf("max|dP|=%.3g NA-mismatch=%d", flat_vs_cache_max, flat_vs_cache_na))

################################################################################
# STEP 1 — INPUT COMPATIBILITY GATE
################################################################################
en <- read.csv(ENTAB_CSV)
# canonical Table 1 mapping (full name / abbreviation / phenotype)
tinfo <- data.frame(
  trait = c("CAD_I25","MI","Angina","AF_I48","HF_I50","Cardiomyopathy","Hypertension_401",
            "Atherosclerosis","PAD","Ischemic_Stroke","TIA","DVT","VTE","Aortic_Aneurysm",
            "Aortic_Valve","Mitral_Valve"),
  full  = c("Coronary artery disease","Myocardial infarction","Angina pectoris","Atrial fibrillation",
            "Heart failure","Cardiomyopathy","Hypertension","Atherosclerosis","Peripheral artery disease",
            "Ischaemic stroke","Transient ischaemic attack","Deep vein thrombosis","Venous thromboembolism",
            "Aortic aneurysm","Aortic valve disease","Mitral valve disease"),
  abbr  = c("CAD","MI","Ang","AF","HF","CMP","HTN","Ath","PAD","IS","TIA","DVT","VTE","AoAn","AoV","MV"),
  pheno = c("ICD-10 I25","ICD-10 I21","ICD-10 I20","ICD-10 I48","ICD-10 I50","ICD-10 I42","Phecode 401",
            "ICD-10 I70","ICD-10 I73","ICD-10 I63","ICD-10 G45","ICD-10 I80","ICD-10 I82","ICD-10 I71",
            "ICD-10 I35","ICD-10 I34"))

per_trait <- do.call(rbind, lapply(trait_names, function(tn) {
  f <- flat_al[, tn]; e <- eff_al[, tn]
  comparable <- !is.na(f) & !is.na(e)
  d <- abs(f - e)
  enr <- en[en$trait == tn, ]
  ti  <- tinfo[tinfo$trait == tn, ]
  cross <- function(thr) sum((f[comparable] < thr) != (e[comparable] < thr))
  data.frame(
    trait = tn, abbreviation = ti$abbr, full_name = ti$full, pan_ukbb_phenotype = ti$pheno,
    N_cases = enr$N_cases_meta, N_controls = enr$N_controls_meta,
    N_eff = enr$N_eff_meta, flat_N = enr$flat_N, N_eff_over_flat_N = enr$ratio_meta_over_flat,
    n_comparable_genes = sum(comparable),
    n_missing_flat = sum(is.na(f)), n_missing_eff = sum(is.na(e)),
    median_abs_delta_p = median(d[comparable]),
    max_abs_delta_p = max(d[comparable]),
    n_cross_p05 = cross(0.05), n_cross_p01 = cross(0.01), n_cross_p001 = cross(0.001))
}))
write.csv(per_trait, file.path(OUT, "per_trait_input_comparison.csv"), row.names = FALSE)

# overall input summary
na_pattern_identical <- all(is.na(flat_al) == is.na(eff_al))
comparable_cells <- sum(!is.na(flat_al) & !is.na(eff_al))
all_delta <- abs(flat_al - eff_al)
input_summary <- data.frame(
  metric = c("n_genes_flat","n_genes_eff","gene_sets_identical","gene_order_identical_after_align",
             "n_traits","trait_order_identical","duplicate_genes_flat","duplicate_genes_eff",
             "na_cells_flat","na_cells_eff","na_pattern_identical","comparable_cells",
             "p_in_unit_interval_flat","p_in_unit_interval_eff",
             "median_abs_delta_p_all","max_abs_delta_p_all","n_changed_cells",
             "n_cross_p05_total","n_cross_p01_total","n_cross_p001_total"),
  value = c(nrow(flat_al), nrow(eff_al), setequal(rownames(flat_mat), rownames(eff_mat)),
            TRUE, n_traits, identical(colnames(flat_al), colnames(eff_al)),
            sum(duplicated(rownames(flat_mat))), sum(duplicated(rownames(eff_mat))),
            sum(is.na(flat_al)), sum(is.na(eff_al)), na_pattern_identical, comparable_cells,
            all(flat_al >= 0 & flat_al <= 1, na.rm = TRUE), all(eff_al >= 0 & eff_al <= 1, na.rm = TRUE),
            median(all_delta, na.rm = TRUE), max(all_delta, na.rm = TRUE),
            sum(all_delta > 0, na.rm = TRUE),
            sum(per_trait$n_cross_p05), sum(per_trait$n_cross_p01), sum(per_trait$n_cross_p001)))
write.csv(input_summary, file.path(OUT, "input_comparison_summary.csv"), row.names = FALSE)

# changed / non-comparable cells only (gene x trait long form)
chg_idx <- which((all_delta > 0) | (is.na(flat_al) != is.na(eff_al)), arr.ind = TRUE)
if (nrow(chg_idx) > 0) {
  disc <- data.frame(
    gene = gene_names[chg_idx[, 1]], trait = trait_names[chg_idx[, 2]],
    p_flat = flat_al[chg_idx], p_eff = eff_al[chg_idx])
  disc$abs_delta_p <- abs(disc$p_flat - disc$p_eff)
  disc$crosses_p05  <- (disc$p_flat < 0.05)  != (disc$p_eff < 0.05)
  disc$crosses_p01  <- (disc$p_flat < 0.01)  != (disc$p_eff < 0.01)
  disc$crosses_p001 <- (disc$p_flat < 0.001) != (disc$p_eff < 0.001)
  disc <- disc[order(-disc$abs_delta_p), ]
} else {
  disc <- data.frame(gene=character(), trait=character(), p_flat=numeric(), p_eff=numeric(),
                     abs_delta_p=numeric(), crosses_p05=logical(), crosses_p01=logical(), crosses_p001=logical())
}
write.csv(disc, file.path(OUT, "gene_input_discrepancies.csv"), row.names = FALSE)

input_ok <- setequal(rownames(flat_mat), rownames(eff_mat)) && identical(colnames(flat_al), colnames(eff_al)) &&
  na_pattern_identical && sum(duplicated(rownames(eff_mat))) == 0 &&
  all(eff_al >= 0 & eff_al <= 1, na.rm = TRUE) &&
  sum(per_trait$n_cross_p05) == 0 && sum(per_trait$n_cross_p01) == 0 && sum(per_trait$n_cross_p001) == 0
setgate("input_compatibility", input_ok,
        sprintf("comparable_cells=%d changed_cells=%d max|dP|=%.3g crossings(05/01/001)=%d/%d/%d",
                comparable_cells, sum(all_delta > 0, na.rm=TRUE), max(all_delta, na.rm=TRUE),
                sum(per_trait$n_cross_p05), sum(per_trait$n_cross_p01), sum(per_trait$n_cross_p001)))

################################################################################
# STEP 2 — BASELINE REPRODUCTION GATE (flat-N through the frozen pipeline)
################################################################################
CTAS_flat <- score_all(flat_al)
# sanity vs the cached Score_C_real
setgate("ctas_recompute_matches_cache", max(abs(CTAS_flat - r$Score_C_real)) == 0,
        sprintf("max|dCTAS|=%.3g", max(abs(CTAS_flat - r$Score_C_real))))
bp <- emp_p_q(CTAS_flat, Score_C_null, M_NULL)
base <- data.frame(gene = gene_names, CTAS = CTAS_flat, null_p = bp$null_p, q_BH = bp$q_BH,
                   s05 = bp$q_BH < 0.05, s01 = bp$q_BH < 0.01, s001 = bp$q_BH < 0.001)
sealed <- read.csv(SEALED_GSS)
m <- merge(base, sealed[, c("gene","magnitude_score","null_p","q_BH",
                            "survives_at_q05","survives_at_q01","survives_at_q001")],
           by = "gene", suffixes = c("", "_sealed"))
d_ctas <- max(abs(m$CTAS - m$magnitude_score))
d_p    <- max(abs(m$null_p - m$null_p_sealed))
d_q    <- max(abs(m$q_BH - m$q_BH_sealed))
mem_ok <- all(m$s05 == m$survives_at_q05) && all(m$s01 == m$survives_at_q01) && all(m$s001 == m$survives_at_q001)
base_counts <- c(q05 = sum(base$s05), q01 = sum(base$s01), q001 = sum(base$s001))
seal_counts <- c(q05 = sum(sealed$survives_at_q05), q01 = sum(sealed$survives_at_q01), q001 = sum(sealed$survives_at_q001))
write.csv(m[order(m$q_BH), ], file.path(OUT, "baseline_reproduction_gene_comparison.csv"), row.names = FALSE)

setgate("baseline_gene_counts", all(base_counts == c(250,102,32)) && all(seal_counts == c(250,102,32)),
        sprintf("recomputed=%d/%d/%d sealed=%d/%d/%d", base_counts[1],base_counts[2],base_counts[3],
                seal_counts[1],seal_counts[2],seal_counts[3]))
setgate("baseline_values_match_sealed", d_ctas < 1e-9 && d_p < 1e-9 && d_q < 1e-9 && mem_ok,
        sprintf("max|dCTAS|=%.3g max|dp|=%.3g max|dq|=%.3g membership_ok=%s", d_ctas, d_p, d_q, mem_ok))

################################################################################
# Locus machinery (locked map, reused; representative rule verbatim)
################################################################################
gla <- read.csv(GLA_CSV)
gla$CHR <- suppressWarnings(as.integer(gla$CHR))
gla <- gla[!is.na(gla$CHR) & gla$CHR %in% 1:22, ]
master <- read.csv(MASTER_CSV)
g2l <- read.csv(G2L_CSV)
# stable-map locus_id -> canonical V2L id (from the locked gene->locus file)
stable2canon <- unique(g2l[, c("stable_map_source_id","canonical_locus_id")])
canon_stable_ids <- unique(g2l$stable_map_source_id)

# generic: given a per-gene table (gene, CTAS, empirical_p, q_BH, s05/s01/s001),
# reproduce canonical loci per threshold + representative per locus via the LOCKED map.
locus_from_genes <- function(gt) {
  sig <- gt[gt$s05, ]
  # attach stable-map locus_id (locked) and canonical id
  sig$stable_id <- gla$locus_id[match(sig$gene, gla$gene)]
  sig$canonical_locus_id <- stable2canon$canonical_locus_id[match(sig$stable_id, stable2canon$stable_map_source_id)]
  unassignable <- sig$gene[is.na(sig$stable_id)]                 # not in the locked map at all
  new_stable   <- sig$gene[!is.na(sig$stable_id) & is.na(sig$canonical_locus_id)]  # stable id w/o baseline canonical
  # per-locus survival + representative (rule: min emp p, max CTAS, min q, gene)
  loci <- do.call(rbind, lapply(split(sig, sig$canonical_locus_id), function(s) {
    if (all(is.na(s$canonical_locus_id))) return(NULL)
    o <- order(s$empirical_p, -s$CTAS, s$q_BH, s$gene)
    s <- s[o, ]
    data.frame(canonical_locus_id = s$canonical_locus_id[1],
               survives_q05 = TRUE,
               survives_q01 = any(s$s01),
               survives_q001 = any(s$s001),
               n_q05_genes = sum(s$s05), n_q01_genes = sum(s$s01), n_q001_genes = sum(s$s001),
               representative = s$gene[1],
               representative_empirical_p = s$empirical_p[1],
               representative_CTAS = s$CTAS[1])
  }))
  list(loci = loci[!is.na(loci$canonical_locus_id), ], unassignable = unassignable, new_stable = new_stable)
}

# baseline gene table for locus reproduction (uses recomputed baseline scores)
base_gt <- data.frame(gene = base$gene, CTAS = base$CTAS, empirical_p = base$null_p, q_BH = base$q_BH,
                      s05 = base$s05, s01 = base$s01, s001 = base$s001)
base_loci <- locus_from_genes(base_gt)
bl <- base_loci$loci
base_locus_counts <- c(q05 = sum(bl$survives_q05), q01 = sum(bl$survives_q01), q001 = sum(bl$survives_q001))
# compare against the locked master
master_counts <- c(q05 = sum(master$survives_q05), q01 = sum(master$survives_q01), q001 = sum(master$survives_q001))
# representative agreement vs locked master
bl$locked_rep <- master$statistical_representative_gene[match(bl$canonical_locus_id, master$canonical_locus_id)]
rep_match <- sum(bl$representative == bl$locked_rep)
write.csv(data.frame(threshold=c("q<0.05","q<0.01","q<0.001"),
                     recomputed_baseline=as.integer(base_locus_counts),
                     locked_master=as.integer(master_counts)),
          file.path(OUT, "baseline_reproduction_locus_comparison.csv"), row.names = FALSE)
setgate("baseline_locus_counts", all(base_locus_counts == c(104,47,16)) && all(master_counts == c(104,47,16)) &&
          length(base_loci$unassignable)==0 && rep_match == nrow(bl),
        sprintf("recomputed=%d/%d/%d master=%d/%d/%d reps_match=%d/%d unassignable=%d",
                base_locus_counts[1],base_locus_counts[2],base_locus_counts[3],
                master_counts[1],master_counts[2],master_counts[3], rep_match, nrow(bl),
                length(base_loci$unassignable)))

baseline_pass <- all(sapply(GATES[c("input_compatibility","baseline_gene_counts",
                                     "baseline_values_match_sealed","baseline_locus_counts")], function(g) g$pass))

################################################################################
# STEP 3 — EFFECTIVE-N SCORING (only the input matrix changes)
################################################################################
CTAS_eff <- score_all(eff_al)
ep <- emp_p_q(CTAS_eff, Score_C_null, M_NULL)
eff <- data.frame(gene = gene_names, CTAS = CTAS_eff, null_p = ep$null_p, q_BH = ep$q_BH,
                  extrap = ep$extrap, s05 = ep$q_BH < 0.05, s01 = ep$q_BH < 0.01, s001 = ep$q_BH < 0.001)
eff_counts <- c(q05 = sum(eff$s05), q01 = sum(eff$s01), q001 = sum(eff$s001))

################################################################################
# STEP 4 — GENE-LEVEL COMPARISON
################################################################################
g <- data.frame(gene = gene_names,
                CTAS_flat = base$CTAS, CTAS_eff = eff$CTAS,
                emp_p_flat = base$null_p, emp_p_eff = eff$null_p,
                q_flat = base$q_BH, q_eff = eff$q_BH)
g$rank_flat <- rank(-g$CTAS_flat, ties.method = "min")
g$rank_eff  <- rank(-g$CTAS_eff,  ties.method = "min")
g$delta_CTAS <- g$CTAS_eff - g$CTAS_flat
g$delta_rank <- g$rank_eff - g$rank_flat
g$q05_flat <- base$s05; g$q05_eff <- eff$s05
g$q01_flat <- base$s01; g$q01_eff <- eff$s01
g$q001_flat <- base$s001; g$q001_eff <- eff$s001
g <- g[order(g$q_flat), ]
write.csv(g, file.path(OUT, "v2_gene_comparison_all.csv"), row.names = FALSE)

memb <- function(fl, ef) list(retained = sum(fl & ef), lost = sum(fl & !ef), gained = sum(!fl & ef),
                              jaccard = sum(fl & ef) / sum(fl | ef),
                              flat_n = sum(fl), eff_n = sum(ef))
m05 <- memb(g$q05_flat, g$q05_eff); m01 <- memb(g$q01_flat, g$q01_eff); m001 <- memb(g$q001_flat, g$q001_eff)
changed <- g[(g$q05_flat != g$q05_eff) | (g$q01_flat != g$q01_eff) | (g$q001_flat != g$q001_eff), ]
write.csv(changed, file.path(OUT, "v2_gene_membership_changes.csv"), row.names = FALSE)

pear <- cor(g$CTAS_flat, g$CTAS_eff, method = "pearson")
spear <- cor(g$CTAS_flat, g$CTAS_eff, method = "spearman")
rank_spear <- cor(g$rank_flat, g$rank_eff, method = "spearman")
nlp_pear <- cor(-log10(g$emp_p_flat), -log10(g$emp_p_eff), method = "pearson")
nlp_spear <- cor(-log10(g$emp_p_flat), -log10(g$emp_p_eff), method = "spearman")
gene_summary <- data.frame(
  metric = c("n_genes","CTAS_pearson","CTAS_spearman","rank_spearman","neglog10p_pearson","neglog10p_spearman",
             "max_abs_delta_CTAS","median_abs_delta_CTAS","max_abs_delta_q","median_abs_delta_q",
             "max_abs_rank_movement","n_genes_rank_changed",
             "q05_flat","q05_eff","q05_retained","q05_lost","q05_gained","q05_jaccard",
             "q01_flat","q01_eff","q01_retained","q01_lost","q01_gained","q01_jaccard",
             "q001_flat","q001_eff","q001_retained","q001_lost","q001_gained","q001_jaccard"),
  value = c(nrow(g), pear, spear, rank_spear, nlp_pear, nlp_spear,
            max(abs(g$delta_CTAS)), median(abs(g$delta_CTAS)),
            max(abs(g$q_eff - g$q_flat)), median(abs(g$q_eff - g$q_flat)),
            max(abs(g$delta_rank)), sum(g$delta_rank != 0),
            m05$flat_n, m05$eff_n, m05$retained, m05$lost, m05$gained, m05$jaccard,
            m01$flat_n, m01$eff_n, m01$retained, m01$lost, m01$gained, m01$jaccard,
            m001$flat_n, m001$eff_n, m001$retained, m001$lost, m001$gained, m001$jaccard))
write.csv(gene_summary, file.path(OUT, "v2_gene_summary.csv"), row.names = FALSE)

top20 <- g[order(-abs(g$delta_rank), g$q_flat), ][1:20,
          c("gene","CTAS_flat","CTAS_eff","delta_CTAS","rank_flat","rank_eff","delta_rank","q_flat","q_eff")]
write.csv(top20, file.path(OUT, "v2_rank_movement_top20.csv"), row.names = FALSE)

################################################################################
# STEP 5 — LOCUS-LEVEL COMPARISON (locked map + representative rule)
################################################################################
eff_gt <- data.frame(gene = eff$gene, CTAS = eff$CTAS, empirical_p = eff$null_p, q_BH = eff$q_BH,
                     s05 = eff$s05, s01 = eff$s01, s001 = eff$s001)
eff_loci_res <- locus_from_genes(eff_gt)
el <- eff_loci_res$loci

# locked canonical locus sets per threshold (from master)
lock_q05 <- master$canonical_locus_id[master$survives_q05]
lock_q01 <- master$canonical_locus_id[master$survives_q01]
lock_q001<- master$canonical_locus_id[master$survives_q001]
eff_q05  <- el$canonical_locus_id[el$survives_q05]
eff_q01  <- el$canonical_locus_id[el$survives_q01]
eff_q001 <- el$canonical_locus_id[el$survives_q001]

locus_row <- function(lock, effs, thr) {
  data.frame(threshold = thr, locked_count = length(lock), effective_n_count = length(effs),
             retained = length(intersect(lock, effs)), lost = length(setdiff(lock, effs)),
             gained = length(setdiff(effs, lock)),
             jaccard = length(intersect(lock, effs)) / length(union(lock, effs)))
}
locus_summary <- rbind(locus_row(lock_q05, eff_q05, "q<0.05"),
                       locus_row(lock_q01, eff_q01, "q<0.01"),
                       locus_row(lock_q001, eff_q001, "q<0.001"))
# representative changes
el$locked_rep <- master$statistical_representative_gene[match(el$canonical_locus_id, master$canonical_locus_id)]
rep_changes <- el[el$representative != el$locked_rep & !is.na(el$locked_rep),
                  c("canonical_locus_id","locked_rep","representative","representative_empirical_p","representative_CTAS")]
names(rep_changes) <- c("canonical_locus_id","locked_representative","effective_n_representative",
                        "eff_representative_empirical_p","eff_representative_CTAS")
locus_summary$n_representative_changes <- c(nrow(rep_changes), NA, NA)
write.csv(locus_summary, file.path(OUT, "v2_locus_summary.csv"), row.names = FALSE)
write.csv(rep_changes, file.path(OUT, "representative_gene_changes.csv"), row.names = FALSE)

# full per-locus comparison table (locked vs effective-N)
comp <- merge(master[, c("canonical_locus_id","chromosome","statistical_representative_gene",
                         "survives_q05","survives_q01","survives_q001",
                         "n_q05_genes","n_q01_genes","n_q001_genes")],
              el[, c("canonical_locus_id","survives_q05","survives_q01","survives_q001",
                     "n_q05_genes","n_q01_genes","n_q001_genes","representative")],
              by = "canonical_locus_id", all = TRUE, suffixes = c("_locked","_eff"))
write.csv(comp[order(comp$canonical_locus_id), ], file.path(OUT, "v2_locus_comparison_all.csv"), row.names = FALSE)

locus_changes <- comp[ (comp$survives_q05_locked %in% TRUE) != (comp$survives_q05_eff %in% TRUE) |
                       (comp$survives_q01_locked %in% TRUE) != (comp$survives_q01_eff %in% TRUE) |
                       (comp$survives_q001_locked %in% TRUE) != (comp$survives_q001_eff %in% TRUE) |
                       (comp$statistical_representative_gene != comp$representative), ]
write.csv(locus_changes, file.path(OUT, "v2_locus_membership_changes.csv"), row.names = FALSE)

################################################################################
# SECONDARY — re-estimate Sigma from the effective-N matrix (clearly labelled)
################################################################################
fit_sigma <- function(pval_matrix) {          # verbatim logic of fit_gaussian_copula()
  pc <- pmin(pmax(pval_matrix, 1e-300), 1 - 1e-10)
  pseudo <- apply(pc, 2, function(x) rank(x, ties.method = "average", na.last = "keep") / (sum(!is.na(x)) + 1))
  cm <- cor(qnorm(pseudo), use = "pairwise.complete.obs", method = "spearman")
  if (min(eigen(cm)$values) < 0.01) cm <- as.matrix(Matrix::nearPD(cm, corr = TRUE)$mat)
  cm
}
Sigma_flat_refit <- fit_sigma(flat_al)        # should ~equal the locked Sigma (same procedure, same input)
Sigma_eff_refit  <- fit_sigma(eff_al)
rownames(Sigma_eff_refit) <- colnames(Sigma_eff_refit) <- trait_names
sig_refit_vs_lock <- max(abs(Sigma_flat_refit - Sigma))    # provenance check on the estimator
sig_max_abs   <- max(abs(Sigma_eff_refit - Sigma))
sig_frob      <- sqrt(sum((Sigma_eff_refit - Sigma)^2))
offdiag <- upper.tri(Sigma)
sig_offdiag_cor <- cor(Sigma[offdiag], Sigma_eff_refit[offdiag])
# whitening change: recompute CTAS with the refit-eff Sigma on the eff matrix
eig2 <- eigen(Sigma_eff_refit, symmetric = TRUE)
Sinv2 <- eig2$vectors %*% diag(1/sqrt(pmax(eig2$values,1e-10))) %*% t(eig2$vectors)
score_C2 <- function(z_vec){ zw <- as.vector(Sinv2 %*% z_vec); w <- pmax(0,zw)^2; s <- sum(w); if(s==0) return(0); sqrt(s)*log2(1+s^2/sum(w^2)) }
z_eff <- p_to_z(eff_al)
CTAS_eff_refitS <- vapply(seq_len(nrow(z_eff)), function(i) score_C2(z_eff[i,]), numeric(1))
sec <- emp_p_q(CTAS_eff_refitS, Score_C_null, M_NULL)   # note: null still locked-Sigma; labelled approximate
sec_counts <- c(q05=sum(sec$q_BH<0.05), q01=sum(sec$q_BH<0.01), q001=sum(sec$q_BH<0.001))
sigma_secondary <- data.frame(
  metric = c("estimator_refit_flat_vs_locked_max_abs","sigma_eff_vs_locked_max_abs","sigma_eff_vs_locked_frobenius",
             "offdiag_correlation","CTAS_pearson_lockedS_vs_refitS_on_eff",
             "secondary_q05_genes","secondary_q01_genes","secondary_q001_genes"),
  value = c(sig_refit_vs_lock, sig_max_abs, sig_frob, sig_offdiag_cor,
            cor(CTAS_eff, CTAS_eff_refitS), sec_counts[1], sec_counts[2], sec_counts[3]))
write.csv(sigma_secondary, file.path(OUT, "sigma_secondary_sensitivity.csv"), row.names = FALSE)
write.csv(as.data.frame(Sigma_eff_refit), file.path(OUT, "sigma_effective_n_refit.csv"))

################################################################################
# GATE 3 — only-input-changed (structural)
################################################################################
setgate("only_input_changed", identical(r$sigma_checksum, sigma_obj_checksum) && flat_vs_cache_max == 0,
        "locked Sigma + locked null reused; only the p-matrix differs")

################################################################################
# Results summary object (for the report + workbook)
################################################################################
res <- list(
  date = as.character(Sys.Date()),
  checksums = list(
    sigma_csv = tools::md5sum(SIGMA_CSV)[[1]], sigma_obj = sigma_obj_checksum,
    flat_tsv = tools::md5sum(FLAT_TSV)[[1]], eff_tsv = tools::md5sum(EFF_TSV)[[1]],
    cache_rds = tools::md5sum(CACHE_RDS)[[1]], sealed_gss = tools::md5sum(SEALED_GSS)[[1]],
    master = tools::md5sum(MASTER_CSV)[[1]], g2l = tools::md5sum(G2L_CSV)[[1]],
    gla = tools::md5sum(GLA_CSV)[[1]]),
  seed = 42, M_NULL = M_NULL,
  eff_n_range = range(en$N_eff_meta),
  max_gene_p_movement = max(all_delta, na.rm = TRUE),
  input_summary = input_summary, per_trait = per_trait,
  base_counts = base_counts, eff_counts = eff_counts,
  base_locus_counts = base_locus_counts,
  eff_locus_counts = c(q05=length(eff_q05), q01=length(eff_q01), q001=length(eff_q001)),
  gene_summary = gene_summary, locus_summary = locus_summary,
  rep_changes = rep_changes, membership_changes = changed,
  sigma_secondary = sigma_secondary, gates = GATES,
  cor = list(pearson=pear, spearman=spear, rank_spearman=rank_spear,
             nlp_pearson=nlp_pear, nlp_spearman=nlp_spear),
  m05=m05, m01=m01, m001=m001,
  n_changed_cells = sum(all_delta > 0, na.rm=TRUE), comparable_cells = comparable_cells)
saveRDS(res, file.path(OUT, "results_summary.rds"))

################################################################################
# Final gate print
################################################################################
msg("\n================= GATE SUMMARY =================")
for (nm in names(GATES)) msg("  %-42s %s", nm, if (GATES[[nm]]$pass) "PASS" else "*** FAIL ***")
overall <- all(sapply(GATES, function(x) x$pass))
msg("OVERALL: %s", if (overall) "PASS" else "FAIL")
msg("\nEffective-N gene counts (q05/q01/q001): %d / %d / %d  (locked 250/102/32)",
    eff_counts[1], eff_counts[2], eff_counts[3])
msg("Effective-N locus counts:               %d / %d / %d  (locked 104/47/16)",
    length(eff_q05), length(eff_q01), length(eff_q001))
msg("CTAS pearson=%.10f spearman=%.10f rank_spearman=%.10f", pear, spear, rank_spear)
msg("max|dCTAS|=%.3g  max|d q|=%.3g  max rank move=%d  genes rank-changed=%d",
    max(abs(g$delta_CTAS)), max(abs(g$q_eff-g$q_flat)), max(abs(g$delta_rank)), sum(g$delta_rank!=0))
msg("q05 ret/lost/gain=%d/%d/%d  q01=%d/%d/%d  q001=%d/%d/%d",
    m05$retained,m05$lost,m05$gained, m01$retained,m01$lost,m01$gained, m001$retained,m001$lost,m001$gained)
msg("locus rep changes=%d ; Sigma eff-refit max|d|=%.4g frob=%.4g offdiag_cor=%.5f",
    nrow(rep_changes), sig_max_abs, sig_frob, sig_offdiag_cor)
saveRDS(overall, file.path(OUT, "overall_pass.rds"))
