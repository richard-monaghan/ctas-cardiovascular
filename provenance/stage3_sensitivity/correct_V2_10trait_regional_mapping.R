#!/usr/bin/env Rscript
################################################################################
# CORRECTION — regional mapping + locus retention for the completed V2 10-trait
# obstructive sensitivity. Reuses the cached, scientifically valid 19,080-gene
# statistics (reduced_gene_results.csv / gene_comparison.csv). NO CTAS/null/
# whitening/p-value/BH rerun. Fixes:
#   1) NCBI37.3.gene.loc column order: SYMBOL is column 1 (was misread as col 6).
#   2) Dual retention metrics: PRIMARY regional_locus_retention (any reduced-panel
#      significant gene mapping into the frozen canonical stable neighbourhood) and
#      SECONDARY locked_member_retention (original strict member-based diagnostic).
#   3) Cross-threshold core context (16- and 47-locus cores at looser thresholds).
# Canonical locus lock remains read-only (md5-verified).
################################################################################
suppressPackageStartupMessages({ library(dplyr) })

CTAS <- "<CTAS_WORKSPACE>"
OUT <- "<CTAS_V2_PAPER_WORKSPACE>/_source_assets/obstructive_sensitivity"
MERGE_DISTANCE <- 1e6

gate <- function(cond, msg) { cat(sprintf("[GATE] %-68s %s\n", msg, if (cond) "PASS" else "*** FAIL ***")); stopifnot(cond) }

## ---- cached gene statistics (unchanged; validated) ----------------------------
red <- read.csv(file.path(OUT, "reduced_gene_results.csv"), stringsAsFactors = FALSE)
cmp <- read.csv(file.path(OUT, "gene_comparison.csv"), stringsAsFactors = FALSE)
gate(nrow(red) == 19080 && nrow(cmp) == 19080, "cached 19,080-gene results reused (no rerun)")
gate(sum(red$sig_q05) == 274 && sum(red$sig_q01) == 77 && sum(red$sig_q001) == 0,
     "cached 274/77/0 gene counts unchanged")
genes <- red$gene
min_q <- min(red$q_BH)
cat(sprintf("[INFO] minimum reduced-panel q = %.5f (zero q<0.001 genes is real)\n", min_q))

## ---- lock integrity (before) ---------------------------------------------------
lock_files <- file.path(CTAS, "panel_v2_locus_lock",
                        c("v2_locus_master.csv", "v2_gene_to_locus.csv", "V2_LOCUS_LOCK.md", "build_v2_locus_lock.R"))
hash_before <- as.character(tools::md5sum(lock_files))

master <- read.csv(file.path(CTAS, "panel_v2_locus_lock", "v2_locus_master.csv"), stringsAsFactors = FALSE)
g2l <- read.csv(file.path(CTAS, "panel_v2_locus_lock", "v2_gene_to_locus.csv"), stringsAsFactors = FALSE)
gla <- read.csv(file.path(CTAS, "panel_locus_cloud", "gene_locus_assignment.csv"), stringsAsFactors = FALSE)
gate(nrow(master) == 104, "canonical source is exactly 104 loci")

## ---- corrected coordinate reference --------------------------------------------
gene_loc <- read.table("<REFERENCE_DIR>/reference/NCBI37.3.gene.loc",
                       header = FALSE, stringsAsFactors = FALSE,
                       col.names = c("SYMBOL", "CHR", "START", "END", "STRAND", "ENTREZ_ID"))
loc_idx <- match(genes, gene_loc$SYMBOL)
loc_df <- data.frame(gene = genes,
                     chr = suppressWarnings(as.integer(gene_loc$CHR[loc_idx])),
                     start = gene_loc$START[loc_idx], end = gene_loc$END[loc_idx],
                     stringsAsFactors = FALSE)

## ---- stable-source -> canonical map --------------------------------------------
src2canon <- setNames(g2l$canonical_locus_id, g2l$stable_map_source_id)  # may repeat per source; collapse
src2canon <- tapply(g2l$canonical_locus_id, g2l$stable_map_source_id, function(x) x[1])
canon_member <- setNames(g2l$canonical_locus_id, g2l$gene)

## ---- per-gene regional disposition ---------------------------------------------
# priority: (1) locked member; (2) stable-map assignment to canonical source;
# (3) corrected coords within canonical neighbourhood (+/-1 Mb); (4) RED10 candidate.
gla_locus <- setNames(gla$locus_id, gla$gene)
sig <- red[red$sig_q05, ]
dispo <- data.frame(gene = sig$gene, reduced_q = sig$q_BH, sig_q01 = sig$sig_q01,
                    map_type = NA_character_, canonical_locus_id = NA_character_,
                    stringsAsFactors = FALSE)
for (i in seq_len(nrow(dispo))) {
  g <- dispo$gene[i]
  if (g %in% names(canon_member)) {
    dispo$map_type[i] <- "locked_member"; dispo$canonical_locus_id[i] <- canon_member[[g]]; next
  }
  st <- gla_locus[g]
  if (!is.na(st) && st %in% names(src2canon)) {
    dispo$map_type[i] <- "stable_map_assignment"; dispo$canonical_locus_id[i] <- src2canon[[st]]; next
  }
  j <- match(g, loc_df$gene)
  if (!is.na(loc_df$chr[j])) {
    hit <- master$canonical_locus_id[master$chromosome == loc_df$chr[j] &
             loc_df$end[j] >= master$locus_start - MERGE_DISTANCE &
             loc_df$start[j] <= master$locus_end + MERGE_DISTANCE]
    if (length(hit) > 0) {
      dispo$map_type[i] <- "canonical_neighbourhood"; dispo$canonical_locus_id[i] <- paste(sort(hit), collapse = ";")
      next
    }
  }
  dispo$map_type[i] <- "RED10_candidate"
}
gate(sum(dispo$map_type == "stable_map_assignment") + sum(dispo$map_type == "locked_member") > 0,
     "stable source IDs used before proximity fallback")

## ---- 38 previously coordinate-less genes now resolved --------------------------
prev_noc <- c("ADH1A","AHSA2","BLK","BTBD16","C3orf18","CDC16","CLCNKA","CSF1R","DBNL","DEDD2",
              "FASTKD1","FNDC1","HEMK1","LRP1B","MFHAS1")
prev_noc_all <- read.csv(file.path(OUT, "superseded_pre_mapping_correction_20260806", "newly_significant.csv"),
                         stringsAsFactors = FALSE)
prev_noc_all <- prev_noc_all$gene[prev_noc_all$map_type == "no_coordinate_info"]
resolved <- prev_noc_all[prev_noc_all %in% gene_loc$SYMBOL]
gate(length(prev_noc_all) == 38 && length(resolved) == 38,
     sprintf("all 38 previously unresolved genes receive corrected coordinates (%d/38)", length(resolved)))
gate(all(loc_df$chr[match(prev_noc_all, loc_df$gene)] %in% 1:22), "corrected coordinates are autosomal and sane")

## ---- RED10 regions (rebuilt, not copied) ---------------------------------------
cand <- dispo[dispo$map_type == "RED10_candidate", ]
cand$chr <- loc_df$chr[match(cand$gene, loc_df$gene)]
cand$start <- loc_df$start[match(cand$gene, loc_df$gene)]
cand$end <- loc_df$end[match(cand$gene, loc_df$gene)]
gate(!any(is.na(cand$chr)), "no coordinate-less gene counted as a region")
canon_src <- names(src2canon)
bridges <- gla[!(gla$locus_id %in% canon_src), c("gene", "CHR", "START", "END")]
names(bridges) <- c("gene", "chr", "start", "end")
bridges$chr <- suppressWarnings(as.integer(bridges$chr))
bridges <- bridges[!is.na(bridges$chr) & bridges$chr %in% 1:22, ]
red10_assign <- character(0)
if (nrow(cand) > 0) {
  cand$is_new <- TRUE; bridges$is_new <- FALSE
  pool <- rbind(cand[, c("gene", "chr", "start", "end", "is_new")], bridges)
  pool <- pool[!is.na(pool$chr), ]
  pool$region <- NA_integer_; rid <- 0
  for (cc in sort(unique(pool$chr))) {
    sub <- pool[pool$chr == cc, ]; sub <- sub[order(sub$start), ]
    cur_end <- -Inf
    for (jj in seq_len(nrow(sub))) {
      if (sub$start[jj] - cur_end > MERGE_DISTANCE) rid <- rid + 1
      sub$region[jj] <- rid; cur_end <- max(cur_end, sub$end[jj])
    }
    pool$region[match(rownames(sub), rownames(pool))] <- sub$region
  }
  new_regions <- sort(unique(pool$region[pool$is_new]))
  red10_ids <- setNames(sprintf("RED10_%03d", seq_along(new_regions)), new_regions)
  red10_assign <- setNames(red10_ids[as.character(pool$region)], pool$gene)
  red10_assign <- red10_assign[names(red10_assign) %in% cand$gene]
}
dispo$RED10_id <- red10_assign[dispo$gene]
dispo$RED10_id[is.na(dispo$RED10_id)] <- ""
dispo$map_type[dispo$map_type == "RED10_candidate"] <- "RED10_region"
n_red10 <- length(unique(dispo$RED10_id[dispo$RED10_id != ""]))

## ---- dual locus-retention metrics ----------------------------------------------
members05 <- strsplit(master$all_q05_genes, ";"); names(members05) <- master$canonical_locus_id
members01 <- strsplit(master$all_q01_genes, ";"); names(members01) <- master$canonical_locus_id
members001 <- strsplit(master$all_q001_genes, ";"); names(members001) <- master$canonical_locus_id
qm05 <- setNames(red$sig_q05, genes); qm01 <- setNames(red$sig_q01, genes); qm001 <- setNames(red$sig_q001, genes)

# regional retention: canonical locus receives >=1 mapped significant gene at threshold
regional_at <- function(sig_flag) {
  ids <- dispo$canonical_locus_id[sig_flag[dispo$gene] & dispo$canonical_locus_id != "" & !is.na(dispo$canonical_locus_id)]
  unique(unlist(strsplit(paste(ids, collapse = ";"), ";")))
}
reg05 <- regional_at(qm05); reg01 <- regional_at(qm01); reg001 <- regional_at(qm001)
gate(length(reg001) == 0, "regional retention at q<0.001 is empty (no q<0.001 genes)")

loc_ret <- data.frame(
  canonical_locus_id = master$canonical_locus_id,
  display_locus_label = master$display_locus_label,
  chromosome = master$chromosome,
  statistical_representative_gene = master$statistical_representative_gene,
  full_tier_q01 = master$survives_q01, full_tier_q001 = master$survives_q001,
  # PRIMARY: regional retention (stable-neighbourhood representation)
  regional_retained_q05 = master$canonical_locus_id %in% reg05,
  regional_retained_q01 = ifelse(master$survives_q01, master$canonical_locus_id %in% reg01, NA),
  regional_retained_q001 = ifelse(master$survives_q001, master$canonical_locus_id %in% reg001, NA),
  # SECONDARY: strict locked-member retention
  locked_member_retained_q05 = sapply(members05, function(g) any(qm05[g])),
  locked_member_retained_q01 = sapply(master$canonical_locus_id, function(id)
    if (master$survives_q01[master$canonical_locus_id == id]) any(qm01[members01[[id]]]) else NA),
  locked_member_retained_q001 = sapply(master$canonical_locus_id, function(id)
    if (master$survives_q001[master$canonical_locus_id == id]) any(qm001[members001[[id]]]) else NA),
  rep_gene_reduced_q = signif(red$q_BH[match(master$statistical_representative_gene, genes)], 4),
  rep_gene_retained_q05 = qm05[master$statistical_representative_gene],
  stringsAsFactors = FALSE)
loc_ret$regional_status_q05 <- ifelse(loc_ret$regional_retained_q05, "retained", "lost")
write.csv(loc_ret, file.path(OUT, "locus_retention.csv"), row.names = FALSE)

RR05 <- sum(loc_ret$regional_retained_q05)
RR01 <- sum(loc_ret$regional_retained_q01[master$survives_q01])
RR001 <- sum(loc_ret$regional_retained_q001[master$survives_q001])
LM05 <- sum(loc_ret$locked_member_retained_q05)
LM01 <- sum(loc_ret$locked_member_retained_q01[master$survives_q01])
LM001 <- sum(loc_ret$locked_member_retained_q001[master$survives_q001])
cat(sprintf("[RESULT] PRIMARY regional retention: q05 %d/104 (%.1f%%) ; q01 %d/47 (%.1f%%) ; q001 %d/16\n",
            RR05, 100 * RR05 / 104, RR01, 100 * RR01 / 47, RR001))
cat(sprintf("[RESULT] SECONDARY locked-member retention: q05 %d/104 (%.1f%%) ; q01 %d/47 (%.1f%%) ; q001 %d/16\n",
            LM05, 100 * LM05 / 104, LM01, 100 * LM01 / 47, LM001))

## ---- cross-threshold core context -----------------------------------------------
core16 <- master$canonical_locus_id[master$survives_q001]
core47 <- master$canonical_locus_id[master$survives_q01]
core16_at_q01 <- sum(core16 %in% reg01); core16_at_q05 <- sum(core16 %in% reg05)
core47_at_q05 <- sum(core47 %in% reg05)
cat(sprintf("[RESULT] cross-threshold: 16-core represented at reduced q01 = %d/16, at reduced q05 = %d/16 ; 47-core at reduced q05 = %d/47\n",
            core16_at_q01, core16_at_q05, core47_at_q05))

## ---- newly significant + disposition for every reduced q05 gene -----------------
newly <- cmp[cmp$status_q05 == "newly_significant" | cmp$status_q01 == "newly_significant" |
             cmp$status_q001 == "newly_significant", ]
newly$strictest_new <- ifelse(newly$status_q001 == "newly_significant", "q<0.001",
                       ifelse(newly$status_q01 == "newly_significant", "q<0.01", "q<0.05"))
di <- match(newly$gene, dispo$gene)
newly$map_type <- dispo$map_type[di]
newly$canonical_locus_id <- dispo$canonical_locus_id[di]
newly$canonical_locus_id[is.na(newly$canonical_locus_id)] <- ""
newly$RED10_id <- dispo$RED10_id[di]
lj <- match(newly$gene, loc_df$gene)
newly$chromosome <- loc_df$chr[lj]; newly$gene_start <- loc_df$start[lj]; newly$gene_end <- loc_df$end[lj]
gate(!any(is.na(newly$map_type)), "every newly significant gene has a documented regional disposition")
gate(!any(is.na(newly$chromosome)), "every newly significant gene has corrected coordinates")
newly_out <- newly[, c("gene", "strictest_new", "status_q05", "status_q01", "status_q001",
                       "reduced_CTAS", "reduced_q", "full_q", "chromosome", "gene_start", "gene_end",
                       "map_type", "canonical_locus_id", "RED10_id")]
write.csv(newly_out, file.path(OUT, "newly_significant.csv"), row.names = FALSE)
gate(nrow(dispo) == 274 && all(dispo$map_type %in% c("locked_member", "stable_map_assignment",
      "canonical_neighbourhood", "RED10_region")), "every reduced q<0.05 gene has a documented regional disposition")

mt <- table(dispo$map_type)
cat(sprintf("[RESULT] disposition of 274 reduced q05 genes: %s\n",
            paste(sprintf("%s=%d", names(mt), as.integer(mt)), collapse = ", ")))
cat(sprintf("[RESULT] RED10 regions (rebuilt): %d covering %d genes\n",
            n_red10, sum(dispo$map_type == "RED10_region")))
nmt <- table(newly$map_type)
cat(sprintf("[RESULT] newly significant (127) disposition: %s\n",
            paste(sprintf("%s=%d", names(nmt), as.integer(nmt)), collapse = ", ")))

## ---- summary --------------------------------------------------------------------
J05 <- 0.3899; J01 <- 0.3358  # unchanged cached values, re-derived below for the record
jac <- function(a, b) length(intersect(a, b)) / length(union(a, b))
J05c <- jac(genes[cmp$full_sig_q05], genes[red$sig_q05])
J01c <- jac(genes[cmp$full_sig_q01], genes[red$sig_q01])
sp <- cor(cmp$full_CTAS, cmp$reduced_CTAS, method = "spearman")
gate(abs(J05c - J05) < 1e-3 && abs(J01c - J01) < 1e-3 && abs(sp - 0.8325) < 1e-3,
     "cached Spearman and Jaccard values unchanged")

summary_rows <- data.frame(
  metric = c(
    "reduced_sig_genes_q05", "reduced_sig_genes_q01", "reduced_sig_genes_q001",
    "min_reduced_q",
    "regional_loci_retained_q05_of_104", "regional_loci_retained_q05_pct",
    "regional_loci_retained_q01_of_47", "regional_loci_retained_q01_pct",
    "regional_loci_retained_q001_of_16",
    "locked_member_retained_q05_of_104", "locked_member_retained_q05_pct",
    "locked_member_retained_q01_of_47", "locked_member_retained_q01_pct",
    "locked_member_retained_q001_of_16",
    "core16_represented_at_reduced_q01", "core16_represented_at_reduced_q05",
    "core47_represented_at_reduced_q05",
    "rep_genes_retained_q05_of_104",
    "genes_retained_q05", "genes_lost_q05", "genes_newly_q05",
    "genes_retained_q01", "genes_lost_q01", "genes_newly_q01",
    "jaccard_q05", "jaccard_q01",
    "spearman_CTAS_full_vs_reduced",
    "disposition_locked_member", "disposition_stable_map_assignment",
    "disposition_canonical_neighbourhood", "disposition_RED10_region",
    "newly_in_canonical_regions", "newly_in_RED10_regions",
    "RED10_regions",
    "reduced_regions_q05_total", "reduced_regions_q05_canonical", "reduced_regions_q05_RED10"),
  value = c(
    274, 77, 0,
    round(min_q, 5),
    RR05, round(100 * RR05 / 104, 1),
    RR01, round(100 * RR01 / 47, 1),
    RR001,
    LM05, round(100 * LM05 / 104, 1),
    LM01, round(100 * LM01 / 47, 1),
    LM001,
    core16_at_q01, core16_at_q05, core47_at_q05,
    sum(loc_ret$rep_gene_retained_q05),
    sum(cmp$status_q05 == "retained"), sum(cmp$status_q05 == "lost"), sum(cmp$status_q05 == "newly_significant"),
    sum(cmp$status_q01 == "retained"), sum(cmp$status_q01 == "lost"), sum(cmp$status_q01 == "newly_significant"),
    round(J05c, 4), round(J01c, 4),
    round(sp, 4),
    as.integer(mt[["locked_member"]]), as.integer(mt[["stable_map_assignment"]]),
    as.integer(mt[["canonical_neighbourhood"]]), as.integer(mt[["RED10_region"]]),
    sum(newly$map_type != "RED10_region"), sum(newly$map_type == "RED10_region"),
    n_red10,
    RR05 + n_red10, RR05, n_red10),
  stringsAsFactors = FALSE)
write.csv(summary_rows, file.path(OUT, "summary.csv"), row.names = FALSE)

## ---- lock integrity (after) ------------------------------------------------------
hash_after <- as.character(tools::md5sum(lock_files))
gate(identical(hash_before, hash_after), "canonical locus lock files byte-unchanged")
gate(!any(grepl("^V2L", dispo$RED10_id)), "no new V2L IDs created")

cat("\nCORRECTION DONE. Rebuilt: locus_retention.csv, newly_significant.csv, summary.csv\n")
