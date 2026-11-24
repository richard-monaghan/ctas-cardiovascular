#!/usr/bin/env Rscript
# Export the signed leave-one-trait-out per-gene q matrix (16 fully recalibrated 15-trait runs; raw-Spearman working matrices; 10^6-draw nulls) to TSV. Read-only.
x <- readRDS("<CTAS_WORKSPACE>/panel_v2_leave_one_out/loo_results.rds")
O <- "<CTAS_WORKSPACE>/outputs/final_ctas_tables_20261021/work/"
q <- as.data.frame(x$Qmat); q <- cbind(gene = rownames(x$Qmat), q)
write.table(q, paste0(O, "loto_q_by_omitted_trait.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cm <- as.data.frame(x$Cmat); write.table(cbind(gene = rownames(x$Cmat), cm), paste0(O, "loto_ctas_by_omitted_trait.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(gene = names(x$full$ctas), CTAS_full16_sens = x$full$ctas, q_full16_sens = x$full$q), paste0(O, "loto_full16_reference.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(dim(x$Qmat), "eig_min_full", x$full$eig_min, "\n")
