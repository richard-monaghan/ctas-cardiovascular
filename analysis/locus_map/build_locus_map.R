#!/usr/bin/env Rscript
################################################################################
# build_locus_map.R -- rebuilds the locked locus map (version 2026-08-05) from minimal discovery-stage inputs.
#
# Cleaned combination of the two original builders (archived verbatim in provenance/stage1_locus_map/):
#   panel_locus_cloud/build_locus_master_table.R  (STEP 1: neighbourhood assignment; logic verbatim)
#   panel_v2_locus_lock/build_v2_locus_lock.R     (canonical loci, boundaries, frozen sort, IDs, mapping notes)
# Removed: columns and steps that do not affect the map (scores, panel-support metrics, representatives,
# tier counts, audit tables, cross-validation printouts). The original final sort used the representative gene
# only as a 4th tie-break; with no (chromosome, start, end) ties among the 104 loci (asserted below) it cannot
# change the order, so pre-final scores are not needed.
#
# Inputs (inputs/locus_map/discovery_stage/, frozen discovery-stage inputs):
#   anchor_regions_earlier_release_101.tsv  101 anchor regions (identifier, chromosome, start, end) from an earlier
#                                           release of this analysis
#   gene_universe_1145.tsv                  genes covered by the map, with the pre-final q<0.05 flag that defined
#                                           which neighbourhoods became the 104 canonical loci
#   gene_coordinates_NCBI37.3_subset.tsv    NCBI37.3 (GRCh37) coordinates of those genes (MAGMA gene-location file)
# Outputs (work/locus_map/): the three files in inputs/locus_map/, byte-identical (checked by tests/reproduce_packaged_outputs.sh)
################################################################################
source(file.path(Sys.getenv("CTAS_REPO", "."), "config", "ctas_paths.R"))
D <- file.path(LOCUSMAP, "discovery_stage"); O <- file.path(WORK_DIR, "locus_map"); dir.create(O, showWarnings = FALSE, recursive = TRUE)
MERGE_DISTANCE <- 1e6

gene_all <- read.delim(file.path(D, "gene_universe_1145.tsv"), stringsAsFactors = FALSE)[, c("gene", "prefinal_q05")]
gene_loc <- read.delim(file.path(D, "gene_coordinates_NCBI37.3_subset.tsv"), stringsAsFactors = FALSE,
                       colClasses = c(GENE = "character", CHR = "character", START = "integer", END = "integer"))
loci_101 <- read.delim(file.path(D, "anchor_regions_earlier_release_101.tsv"))

# ---- STEP 1 (verbatim logic): anchor overlap within 1 Mb, else sequential 1 Mb merge of the remaining genes
gene_all_loc <- merge(gene_all, gene_loc[, c("GENE", "CHR", "START", "END")], by.x = "gene", by.y = "GENE", all.x = TRUE)
stopifnot(!anyNA(gene_all_loc$CHR))
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
new_genes <- new_genes[order(new_genes$CHR, new_genes$START), ]          # CHR is character, as in the original
new_locus_id <- integer(nrow(new_genes)); new_locus_id[1] <- 1
last_chr <- new_genes$CHR[1]; last_end <- new_genes$END[1]; cur <- 1
for (i in 2:nrow(new_genes)) {
  if (new_genes$CHR[i] == last_chr && (new_genes$START[i] - last_end) <= MERGE_DISTANCE) {
    new_locus_id[i] <- cur; last_end <- max(last_end, new_genes$END[i])
  } else { cur <- cur + 1; new_locus_id[i] <- cur; last_chr <- new_genes$CHR[i]; last_end <- new_genes$END[i] }
}
new_genes$new_locus_num <- new_locus_id
gene_all_loc$locus_id <- ifelse(gene_all_loc$locus_status == "overlap", paste0("ST2-", gene_all_loc$overlapping_locus_id),
                                paste0("NEW-", new_genes$new_locus_num[match(gene_all_loc$gene, new_genes$gene)]))

# ---- canonical loci (lock logic): neighbourhoods containing >=1 pre-final q<0.05 gene, autosomes only
full <- gene_all_loc[, c("gene", "CHR", "START", "END", "locus_id", "prefinal_q05")]
full$CHR <- suppressWarnings(as.integer(full$CHR)); full <- full[!is.na(full$CHR) & full$CHR %in% 1:22, ]
canon_ids <- sort(unique(full$locus_id[full$prefinal_q05])); stopifnot(length(canon_ids) == 104)
fc <- full[full$locus_id %in% canon_ids, ]
bounds <- do.call(rbind, lapply(split(fc, fc$locus_id), function(d) { stopifnot(length(unique(d$CHR)) == 1)
  data.frame(locus_id = d$locus_id[1], chromosome = d$CHR[1], locus_start = min(d$START), locus_end = max(d$END), stringsAsFactors = FALSE) }))
stopifnot(!anyDuplicated(paste(bounds$chromosome, bounds$locus_start, bounds$locus_end)))   # 4th sort key (representative) cannot matter
master <- bounds[order(bounds$chromosome, bounds$locus_start, bounds$locus_end), ]
master$canonical_locus_id <- sprintf("V2L%03d", seq_len(nrow(master)))
master$legacy_v1_locus_id <- ifelse(grepl("^ST2-", master$locus_id), master$locus_id, "")
master$v1_overlap_status <- ifelse(grepl("^ST2-", master$locus_id), "existing_v1_locus", "new_v2_locus")
master$mapping_notes <- ""
big <- (master$locus_end - master$locus_start) > 2e6
master$mapping_notes[big] <- paste0("extended_region_span_", round((master$locus_end[big] - master$locus_start[big]) / 1e6, 2), "Mb")
master$mapping_notes[master$locus_id == "NEW-305"] <- paste(master$mapping_notes[master$locus_id == "NEW-305"],
  "CLDN23;MSRA;SOX7 = one locus via the stable neighbourhood map; TNKS is a non-significant positional bridge only", sep = "; ")

# ---- write in the packaged format (all fields quoted text, as in inputs/locus_map/)
chr <- function(d) { d[] <- lapply(d, as.character); d }
write.csv(chr(gene_all_loc[, c("gene", "CHR", "START", "END", "locus_id")]), file.path(O, "gene_neighbourhood_map.csv"), row.names = FALSE)
s2c <- data.frame(stable_map_source_id = master$locus_id, canonical_locus_id = master$canonical_locus_id, stringsAsFactors = FALSE)
s2c <- s2c[order(s2c$canonical_locus_id, s2c$stable_map_source_id), ]
write.csv(chr(s2c), file.path(O, "stable_to_canonical_locus.csv"), row.names = FALSE)
write.csv(chr(master[, c("canonical_locus_id", "chromosome", "locus_start", "locus_end", "legacy_v1_locus_id", "v1_overlap_status", "mapping_notes")]),
          file.path(O, "canonical_locus_master.csv"), row.names = FALSE)
cat(sprintf("%d genes; %d neighbourhoods; %d canonical loci -> %s\n", nrow(gene_all_loc), length(unique(gene_all_loc$locus_id)), nrow(master), O))
