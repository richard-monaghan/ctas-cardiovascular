#!/usr/bin/env Rscript
################################################################################
# CTAS V2 colocalisation — Phase 5 locus-level abstraction & summaries.
# Reads all per-locus checkpoints, emits the long results table, the 104-row
# locus abstraction, and stratified summaries (tissue/trait/inherited-vs-new/
# q-tier) plus a descriptive V1 comparison. No V1 tiers.
################################################################################
suppressPackageStartupMessages({ library(data.table) })
source("<CTAS_V2_PAPER_WORKSPACE>/_source_assets/scripts/v2_coloc_config.R")
OUT <- V2$COLOC_DIR
mani <- fread(V2$MANIFEST); coords <- fread(V2$COORDS)

ckpts <- list.files(V2$CKPT_DIR, pattern="^V2L\\d+\\.tsv$", full.names=TRUE)
stopifnot(length(ckpts) > 0)
long <- rbindlist(lapply(ckpts, fread), fill=TRUE)
setorder(long, canonical_locus_id, tissue, -PP4)
fwrite(long, file.path(OUT, "V2_COLOC_LONG_RESULTS.tsv"), sep="\t")

SUP <- 0.5; STRONG <- 0.8
ass <- long[disposition=="assessed"]
locabs <- long[, {
  a <- .SD[disposition=="assessed"]
  n_ass <- nrow(a)
  mx4 <- if (n_ass) max(a$PP4) else NA_real_
  mx3 <- if (n_ass) max(a$PP3) else NA_real_
  wi  <- if (n_ass) which.max(a$PP4) else NA_integer_
  is_mhc_any <- any(disposition=="not_assessed_MHC")
  cls <- if (is_mhc_any) "not_assessed" else if (n_ass==0) "not_assessed"
         else if (mx4>=STRONG) "strong_PP4" else if (mx4>=SUP) "moderate_PP4"
         else if (mx3>=SUP) "PP3_only" else "neither"
  .(gene=gene[1], ensg=ensg[1],
    n_cells=.N, n_assessed=n_ass,
    n_PP4_sup=if(n_ass) sum(a$PP4>=SUP) else 0L,
    n_PP4_strong=if(n_ass) sum(a$PP4>=STRONG) else 0L,
    n_PP3_sup=if(n_ass) sum(a$PP3>=SUP) else 0L,
    max_PP4=round(ifelse(is.na(mx4),NA,mx4),4), max_PP3=round(ifelse(is.na(mx3),NA,mx3),4),
    best_tissue=if(n_ass) a$tissue[wi] else NA_character_,
    best_trait =if(n_ass) a$trait[wi]  else NA_character_,
    coloc_class=cls, mhc=is_mhc_any)
}, by=canonical_locus_id]

locabs <- merge(mani[, .(canonical_locus_id, chromosome, target_category, v1_overlap_status,
                         legacy_v1_locus_id, v1_lead_gene, survives_q05, survives_q01, survives_q001)],
                locabs, by="canonical_locus_id")
locabs[, inherited := target_category != "new_v2_locus"]
setorder(locabs, canonical_locus_id)
fwrite(locabs, file.path(OUT, "V2_COLOC_LOCUS_ABSTRACTION.tsv"), sep="\t")

## ---- stratified summaries ---------------------------------------------------
by_class <- locabs[, .N, by=coloc_class][order(-N)]
fwrite(by_class, file.path(OUT, "V2_COLOC_SUMMARY_by_class.tsv"), sep="\t")

by_tissue <- ass[, .(n_assessed=.N, n_sup=sum(PP4>=SUP), n_strong=sum(PP4>=STRONG),
                     max_PP4=round(max(PP4),3)), by=tissue][order(-n_sup)]
fwrite(by_tissue, file.path(OUT, "V2_COLOC_SUMMARY_by_tissue.tsv"), sep="\t")

by_trait <- ass[, .(n_assessed=.N, n_sup=sum(PP4>=SUP), n_strong=sum(PP4>=STRONG),
                    max_PP4=round(max(PP4),3)), by=trait][order(-n_sup)]
fwrite(by_trait, file.path(OUT, "V2_COLOC_SUMMARY_by_trait.tsv"), sep="\t")

by_inherit <- locabs[, .(n_loci=.N, strong=sum(coloc_class=="strong_PP4"),
                         moderate=sum(coloc_class=="moderate_PP4"),
                         pp3_only=sum(coloc_class=="PP3_only"),
                         neither=sum(coloc_class=="neither"),
                         not_assessed=sum(coloc_class=="not_assessed")),
                     by=.(inherited)]
fwrite(by_inherit, file.path(OUT, "V2_COLOC_SUMMARY_by_inherited.tsv"), sep="\t")

qtier <- rbindlist(list(
  locabs[survives_q05==TRUE,  .(tier="q<0.05", n=.N, has_sup=sum(coloc_class %in% c("strong_PP4","moderate_PP4")))],
  locabs[survives_q01==TRUE,  .(tier="q<0.01", n=.N, has_sup=sum(coloc_class %in% c("strong_PP4","moderate_PP4")))],
  locabs[survives_q001==TRUE, .(tier="q<0.001",n=.N, has_sup=sum(coloc_class %in% c("strong_PP4","moderate_PP4")))]))
fwrite(qtier, file.path(OUT, "V2_COLOC_SUMMARY_by_qtier.tsv"), sep="\t")

## ---- descriptive V1 comparison (48 unchanged / 19 changed / 37 new) --------
v1res <- fread(file.path(V2$VIGNETTE, "results_all_101", "coloc_pp4_all_101.tsv"))
v1_per_gene <- v1res[, .(v1_max_PP4=max(PP4), v1_sup=any(PP4>=SUP), v1_strong=any(PP4>=STRONG)), by=gene]
cmp <- merge(locabs, v1_per_gene, by.x="gene", by.y="gene", all.x=TRUE)
cmp_out <- cmp[, .(canonical_locus_id, gene, target_category, coloc_class,
                   max_PP4, v1_max_PP4=round(v1_max_PP4,4),
                   v2_sup = coloc_class %in% c("strong_PP4","moderate_PP4"), v1_sup)]
setorder(cmp_out, canonical_locus_id)
fwrite(cmp_out, file.path(OUT, "V2_COLOC_V1_DESCRIPTIVE_COMPARISON.tsv"), sep="\t")

cat("=== Phase 5 abstraction complete ===\n")
cat(sprintf("long results: %d rows ; loci: %d\n", nrow(long), nrow(locabs)))
cat("\nclass distribution:\n"); print(by_class)
cat(sprintf("\nsupported (strong+moderate) loci: %d / %d\n",
            sum(locabs$coloc_class %in% c("strong_PP4","moderate_PP4")), nrow(locabs)))
cat("\nby inherited vs new:\n"); print(by_inherit)
cat("\ntop tissues:\n"); print(head(by_tissue))
cat("\ntop traits:\n"); print(head(by_trait))
