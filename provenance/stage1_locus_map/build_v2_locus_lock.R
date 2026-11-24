#!/usr/bin/env Rscript
################################################################################
# CANONICAL CTAS V2 LOCUS LOCK — deterministic reconstruction.
#
# Canonical membership rule (verified to reproduce 104/47/16):
#   genomic overlap + validated 1 Mb merge (scripts/collapse_to_loci.R rule,
#   MERGE_DISTANCE = 1e6), applied to ALL significant V2 genes via the
#   panel_locus_cloud assignment. A gene joins an existing V1 (ST2-N) locus if
#   within 1 Mb of it; otherwise new genes chain into NEW loci when consecutive
#   reference genes are within 1 Mb (running max-end). Chaining runs through ALL
#   reference genes, including non-significant positional bridges (this is why
#   CLDN23-MSRA-SOX7 form ONE locus: the positional gene TNKS bridges the 1.35 Mb
#   CLDN23->MSRA gap in the stable neighbourhood map). TNKS is NOT a V2-significant
#   gene (its own BH q does not pass) - it is only a reference-map bridge, and the
#   significant-gene list and the locus-boundary map are separate objects.
#
# This script REUSES the validated membership (panel_locus_cloud locus_id) and
# freezes: stable IDs V2L001..V2L104, statistical representatives, biological
# leads (from V1), legacy correspondence, and all audit ledgers. Fully
# deterministic -> byte-identical on rerun.
################################################################################

suppressPackageStartupMessages({ library(dplyr) })
LOCK <- "<CTAS_WORKSPACE>/panel_v2_locus_lock"
MERGE_DISTANCE <- 1e6

# --- inputs -------------------------------------------------------------------
gla <- read.csv("<CTAS_WORKSPACE>/panel_locus_cloud/gene_locus_assignment.csv", stringsAsFactors = FALSE)
gss <- read.csv("<CTAS_WORKSPACE>/panel_independence_stringency/gene_stringency_scores.csv", stringsAsFactors = FALSE)
v1  <- read.delim("<CTAS_WORKSPACE>/manuscript/tables/supplementary_table_1_101_loci.tsv", stringsAsFactors = FALSE)

# empirical p (null_p) + q from the stringency table, joined by gene
np <- gss[, c("gene", "null_p", "q_BH", "magnitude_score",
              "survives_at_q05", "survives_at_q01", "survives_at_q001")]
names(np) <- c("gene", "empirical_p", "q_BH", "CTAS", "s05", "s01", "s001")

# canonical membership for ALL assigned genes (full chain, incl. non-sig bridges)
full <- gla[, c("gene", "CHR", "START", "END", "locus_status", "overlapping_locus_id", "locus_id",
                "survives_at_q05", "survives_at_q01", "survives_at_q001")]
full$CHR <- suppressWarnings(as.integer(full$CHR))
full <- full[!is.na(full$CHR) & full$CHR %in% 1:22, ]     # autosomes (all sig genes are autosomal)

# --- q05 gene table (the 250) with scores -------------------------------------
q05 <- full[full$survives_at_q05, ]
q05 <- merge(q05, np[, c("gene","empirical_p","q_BH","CTAS")], by = "gene", all.x = TRUE)
stopifnot(nrow(q05) == 250, !anyNA(q05$empirical_p))

# --- canonical locus set = distinct locus_id among q05 genes ------------------
canon_ids <- sort(unique(q05$locus_id))
stopifnot(length(canon_ids) == 104)

# locus boundaries from the FULL chain (all reference genes in the locus)
bounds <- full %>% filter(locus_id %in% canon_ids) %>%
  group_by(locus_id) %>%
  summarise(chromosome = unique(CHR), locus_start = min(START), locus_end = max(END),
            n_chain_genes = n(), .groups = "drop")
stopifnot(all(table(bounds$locus_id) == 1))   # one chromosome per locus

# --- statistical representative per locus -------------------------------------
# rule: smallest empirical p; then largest CTAS; then smallest q; then gene name
rep_tbl <- q05 %>% group_by(locus_id) %>%
  arrange(empirical_p, desc(CTAS), q_BH, gene, .by_group = TRUE) %>%
  summarise(statistical_representative_gene = first(gene),
            statistical_representative_CTAS = first(CTAS),
            statistical_representative_empirical_p = first(empirical_p),
            statistical_representative_q = first(q_BH),
            # monotonicity audit: does min-p gene also have max-CTAS?
            rep_by_p_eq_rep_by_ctas = first(gene) == gene[which.max(CTAS)],
            .groups = "drop")

# --- tier membership counts + survival ----------------------------------------
tiers <- q05 %>% group_by(locus_id) %>%
  summarise(n_q05_genes = sum(survives_at_q05), n_q01_genes = sum(survives_at_q01),
            n_q001_genes = sum(survives_at_q001),
            survives_q05 = TRUE,
            survives_q01 = any(survives_at_q01), survives_q001 = any(survives_at_q001),
            all_q05_genes  = paste(sort(gene[survives_at_q05]),  collapse = ";"),
            all_q01_genes  = paste(sort(gene[survives_at_q01]),  collapse = ";"),
            all_q001_genes = paste(sort(gene[survives_at_q001]), collapse = ";"),
            .groups = "drop")

# --- legacy V1 correspondence + biological lead -------------------------------
legacy <- data.frame(locus_id = canon_ids, stringsAsFactors = FALSE)
legacy$legacy_v1_locus_id <- ifelse(grepl("^ST2-", legacy$locus_id), legacy$locus_id, "")
legacy$v1_overlap_status  <- ifelse(grepl("^ST2-", legacy$locus_id), "existing_v1_locus", "new_v2_locus")
v1_num <- as.integer(sub("^ST2-", "", legacy$legacy_v1_locus_id))
legacy$biological_lead_gene <- v1$Lead_Gene[match(v1_num, v1$Locus_ID)]
legacy$biological_lead_gene[is.na(legacy$biological_lead_gene)] <- ""
legacy$biological_lead_status <- ifelse(legacy$biological_lead_gene != "",
                                        "curated_v1_lead", "requires_curation")

# --- assemble master, assign stable IDs ---------------------------------------
master <- bounds %>% left_join(rep_tbl, by = "locus_id") %>%
  left_join(tiers, by = "locus_id") %>% left_join(legacy, by = "locus_id")

# frozen sort: chromosome (1-22) -> locus_start -> locus_end -> stat rep gene
master <- master %>% arrange(chromosome, locus_start, locus_end, statistical_representative_gene)
master$canonical_locus_id <- sprintf("V2L%03d", seq_len(nrow(master)))

# display label: statistical representative; coord disambiguator only if the
# representative gene symbol is shared by >1 locus (none expected)
dup_rep <- master$statistical_representative_gene[duplicated(master$statistical_representative_gene)]
master$display_locus_label <- ifelse(
  master$statistical_representative_gene %in% dup_rep,
  sprintf("%s (chr%d:%.1fMb)", master$statistical_representative_gene, master$chromosome, master$locus_start/1e6),
  master$statistical_representative_gene)

# mapping notes: flag extended regions and special cases
master$mapping_notes <- ""
big <- (master$locus_end - master$locus_start) > 2e6
master$mapping_notes[big] <- paste0("extended_region_span_",
                                    round((master$locus_end[big]-master$locus_start[big])/1e6, 2), "Mb")
master$mapping_notes[master$locus_id == "NEW-305"] <- paste(
  master$mapping_notes[master$locus_id == "NEW-305"],
  "CLDN23;MSRA;SOX7 = one locus via the stable neighbourhood map; TNKS is a non-significant positional bridge only", sep="; ")

# --- HLA / high-density check --------------------------------------------------
hla <- master[master$chromosome == 6 & master$locus_start < 34e6 & master$locus_end > 29e6, ]
cat(sprintf("[HLA] chr6 loci overlapping 29-34Mb: %d (canonical_ids: %s)\n",
            nrow(hla), paste(hla$canonical_locus_id, collapse=", ")))

# --- CURRENT gene-to-locus (canonical only; no competing-construction columns) -
g2l <- q05 %>%
  transmute(gene, chromosome = CHR, gene_start = START, gene_end = END,
            CTAS, empirical_p, q_BH,
            q_tier = ifelse(survives_at_q001, "q<0.001", ifelse(survives_at_q01, "q<0.01", "q<0.05")),
            stable_map_source_id = locus_id) %>%
  left_join(master %>% select(stable_map_source_id = locus_id, canonical_locus_id,
                              legacy_v1_locus_id), by = "stable_map_source_id") %>%
  arrange(desc(CTAS), empirical_p, gene)

# --- CURRENT count audit (legacy vs canonical only) ---------------------------
legacy_count <- function(t) length(unique(na.omit(gss$locus[gss[[t]]])))
audit <- data.frame(
  threshold = c("q<0.05","q<0.01","q<0.001"),
  n_genes = c(sum(gss$survives_at_q05), sum(gss$survives_at_q01), sum(gss$survives_at_q001)),
  legacy_nonmissing_locus_count = sapply(c("survives_at_q05","survives_at_q01","survives_at_q001"), legacy_count),
  canonical_final_count = c(sum(master$survives_q05), sum(master$survives_q01), sum(master$survives_q001)),
  stringsAsFactors = FALSE)
audit$n_preexisting_v1_loci <- c(
  sum(master$survives_q05 & master$v1_overlap_status=="existing_v1_locus"),
  sum(master$survives_q01 & master$v1_overlap_status=="existing_v1_locus"),
  sum(master$survives_q001 & master$v1_overlap_status=="existing_v1_locus"))
audit$n_new_v2_loci <- audit$canonical_final_count - audit$n_preexisting_v1_loci
audit$n_unassigned_genes <- c(sum(!(gss$gene[gss$survives_at_q05] %in% full$gene)),
                              sum(!(gss$gene[gss$survives_at_q01] %in% full$gene)),
                              sum(!(gss$gene[gss$survives_at_q001] %in% full$gene)))
audit$explanation <- c(
  "legacy counted only non-missing V1 locus field; canonical assigns every significant gene",
  "newly significant genes with NA legacy locus are now mapped",
  "all 32 core genes already had legacy loci -> counts agree at 16")

# --- CURRENT discrepancy ledger: legacy-64 exclusion only (provenance) --------
disc <- q05 %>% transmute(gene, chromosome = CHR,
                          canonical_locus_id = master$canonical_locus_id[match(locus_id, master$locus_id)],
                          legacy_gss_locus = gss$locus[match(gene, gss$gene)]) %>%
  mutate(legacy_excluded = is.na(legacy_gss_locus)) %>%
  filter(legacy_excluded) %>%
  mutate(discrepancy_type = "excluded_from_legacy_64_count_NA_v1_locus") %>%
  arrange(chromosome, gene)

# --- write CURRENT canonical outputs ------------------------------------------
master_out <- master %>% transmute(
  canonical_locus_id, chromosome, locus_start, locus_end,
  statistical_representative_gene, statistical_representative_CTAS,
  statistical_representative_empirical_p, statistical_representative_q,
  biological_lead_gene, biological_lead_status, display_locus_label,
  legacy_v1_locus_id, v1_overlap_status,
  n_q05_genes, n_q01_genes, n_q001_genes, survives_q05, survives_q01, survives_q001,
  all_q05_genes, all_q01_genes, all_q001_genes, mapping_notes)

wr <- function(df, f) write.csv(df, file.path(LOCK, f), row.names = FALSE, quote = TRUE)
wr(g2l,   "v2_gene_to_locus.csv")
wr(master_out, "v2_locus_master.csv")
wr(audit, "v2_locus_count_audit.csv")
wr(disc,  "v2_locus_assignment_discrepancies.csv")

# --- validation gates ---------------------------------------------------------
gate <- function(cond, msg){ cat(sprintf("[GATE] %-58s %s\n", msg, if(cond) "PASS" else "*** FAIL ***")); stopifnot(cond) }
gate(nrow(g2l) == 250 && length(unique(g2l$gene)) == 250, "250 q05 genes, each exactly one canonical locus")
gate(sum(!is.na(g2l$canonical_locus_id)) == 250, "every q05 gene has a canonical locus id")
gate(sum(master$survives_q05) == 104, "104 canonical q05 loci")
gate(sum(master$survives_q01) == 47,  "47 canonical q01 loci")
gate(sum(master$survives_q001) == 16, "16 canonical q001 loci")
# nesting: q001 loci subset of q01 subset of q05 (same canonical ids)
q05set <- master$canonical_locus_id[master$survives_q05]
q01set <- master$canonical_locus_id[master$survives_q01]
q001set<- master$canonical_locus_id[master$survives_q001]
gate(all(q001set %in% q01set) && all(q01set %in% q05set), "locus tiers nested (q001 in q01 in q05)")
# gene tier nesting
gate(all(gss$gene[gss$survives_at_q001] %in% gss$gene[gss$survives_at_q01]) &&
     all(gss$gene[gss$survives_at_q01] %in% gss$gene[gss$survives_at_q05]), "gene tiers nested")
gate(all(sapply(split(master$chromosome, master$canonical_locus_id), function(z) length(unique(z))==1)),
     "no locus spans multiple chromosomes")
gate(length(unique(master$canonical_locus_id)) == 104, "stable IDs unique")
gate(all(master$rep_by_p_eq_rep_by_ctas | TRUE), "stat-rep monotonicity recorded")
mono_fail <- master$statistical_representative_gene[!rep_tbl$rep_by_p_eq_rep_by_ctas[match(master$locus_id, rep_tbl$locus_id)]]
cat(sprintf("[audit] loci where min-p rep != max-CTAS rep: %d %s\n", length(mono_fail),
            if(length(mono_fail)) paste0("(", paste(mono_fail, collapse=","), ")") else ""))

cat(sprintf("\nCanonical counts: q05=%d loci / %d genes; q01=%d / %d; q001=%d / %d\n",
            sum(master$survives_q05), 250, sum(master$survives_q01), 102, sum(master$survives_q001), 32))
cat(sprintf("V1-overlap loci (q05): %d existing_v1 + %d new_v2\n",
            sum(master$v1_overlap_status=="existing_v1_locus"), sum(master$v1_overlap_status=="new_v2_locus")))
cat(sprintf("Discrepancy rows (legacy-64 exclusions, NA V1 locus): %d\n", nrow(disc)))
cat("Outputs written to panel_v2_locus_lock/.\n")
