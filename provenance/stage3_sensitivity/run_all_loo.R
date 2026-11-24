#!/usr/bin/env Rscript
################################################################################
# Run all 16 leave-one-trait-out reruns; compute locus/representative retention,
# genome-wide stability, stable-16 subset and newly-significant diagnostics
# against the LOCKED 104-locus canonical object. Writes CSVs for the ST6 builder.
################################################################################
source("<CTAS_WORKSPACE>/panel_v2_leave_one_out/loo_engine.R")
OUT <- "<CTAS_V2_PAPER_WORKSPACE>/_source_assets/R3_review"  # figure-ready area (R4-adjacent)
LOO <- "<CTAS_WORKSPACE>/panel_v2_leave_one_out"

# ---- baseline: full result + locked loci -------------------------------------
full <- if (file.exists(file.path(LOO,"full16_reference.rds"))) readRDS(file.path(LOO,"full16_reference.rds")) else run_ctas(Pfull)
names(full$q) <- names(full$p) <- names(full$ctas) <- genes
gss <- read.csv("<CTAS_WORKSPACE>/panel_independence_stringency/gene_stringency_scores.csv", stringsAsFactors=FALSE)
master <- read.csv("<CTAS_WORKSPACE>/panel_v2_locus_lock/v2_locus_master.csv", stringsAsFactors=FALSE)
stopifnot(nrow(master)==104)
members <- strsplit(master$all_q05_genes, ";")           # locked member genes per locus
names(members) <- master$canonical_locus_id
canon_q05 <- gss$gene[gss$survives_at_q05]; stopifnot(length(canon_q05)==250)

# ---- run 16 LOO --------------------------------------------------------------
Qmat <- matrix(NA_real_, nrow=length(genes), ncol=16, dimnames=list(genes, trait_names))
Cmat <- Qmat
loo_meta <- data.frame()
for (j in seq_along(trait_names)) {
  tt <- trait_names[j]
  r <- run_ctas(Pfull[, -j, drop=FALSE])
  Qmat[,j] <- r$q; Cmat[,j] <- r$ctas
  loo_meta <- rbind(loo_meta, data.frame(omitted_trait=tt, n_traits=ncol(Pfull)-1,
    sigma_symmetric=isSymmetric(unname(r$sigma)), sigma_min_eig=round(r$eig_min,4),
    cor_min=round(r$cor_rng[1],3), cor_max=round(r$cor_rng[2],3),
    n_q05=sum(r$q<0.05), n_q01=sum(r$q<0.01), n_q001=sum(r$q<0.001)))
  cat(sprintf("[LOO -%s] q05=%d q01=%d q001=%d\n", tt, sum(r$q<0.05), sum(r$q<0.01), sum(r$q<0.001)))
}

# ---- locus retention (primary + tiers + representative) ----------------------
mem_min_q <- function(gset, j) min(Qmat[gset, j])            # best (smallest) member q in run j
ret05 <- sapply(seq_along(trait_names), function(j) sapply(members, function(g) mem_min_q(g,j) < 0.05))
ret01 <- sapply(seq_along(trait_names), function(j) sapply(members, function(g) mem_min_q(g,j) < 0.01))
ret001<- sapply(seq_along(trait_names), function(j) sapply(members, function(g) mem_min_q(g,j) < 0.001))
repret<- sapply(seq_along(trait_names), function(j) Qmat[master$statistical_representative_gene, j] < 0.05)
colnames(ret05) <- colnames(ret01) <- colnames(ret001) <- colnames(repret) <- trait_names
rownames(repret) <- master$canonical_locus_id

n_ret <- rowSums(ret05)                                     # retained runs / 16 (primary)
locus_tab <- data.frame(
  canonical_locus_id=master$canonical_locus_id, chromosome=master$chromosome,
  statistical_representative_gene=master$statistical_representative_gene,
  n_member_genes=master$n_q05_genes,
  full_CTAS_rep=round(master$statistical_representative_CTAS,3),
  full_q_rep=master$statistical_representative_q,
  n_runs_retained=n_ret, retention_proportion=round(n_ret/16,3),
  rep_n_runs_retained=rowSums(repret),
  stringsAsFactors=FALSE)
locus_tab <- cbind(locus_tab, as.data.frame(ret05+0))       # 16 retention flags (0/1)
locus_tab <- locus_tab[order(-locus_tab$n_runs_retained, order(match(locus_tab$canonical_locus_id, master$canonical_locus_id))),]

# ---- per-omitted-trait genome-wide stability + retention ---------------------
full05 <- genes[full$q<0.05]
omit_sum <- do.call(rbind, lapply(seq_along(trait_names), function(j){
  loo05 <- genes[Qmat[,j]<0.05]
  jac <- length(intersect(full05,loo05))/length(union(full05,loo05))
  newly <- setdiff(loo05, canon_q05)
  data.frame(omitted_trait=trait_names[j],
    n_loci_retained_q05=sum(ret05[,j]), locus_retention_pct=round(100*mean(ret05[,j]),1),
    n_rep_retained=sum(repret[,j]), n_loci_retained_q01=sum(ret01[,j]), n_loci_retained_q001=sum(ret001[,j]),
    gene_q05_jaccard=round(jac,4), n_loo_q05_genes=length(loo05),
    ctas_spearman_full_vs_loo=round(cor(full$ctas, Cmat[,j], method="spearman"),4),
    n_newly_significant_genes=length(newly))
}))
omit_sum <- omit_sum[order(omit_sum$n_loci_retained_q05),]

# ---- stable 16/16 + newly significant diagnostic ----------------------------
stable16 <- locus_tab[locus_tab$n_runs_retained==16, c("canonical_locus_id","chromosome","statistical_representative_gene","n_member_genes","full_CTAS_rep")]
newly_diag <- do.call(rbind, lapply(seq_along(trait_names), function(j){
  loo05 <- genes[Qmat[,j]<0.05]; newly <- setdiff(loo05, canon_q05)
  if(length(newly)==0) return(NULL)
  data.frame(omitted_trait=trait_names[j], gene=newly, loo_CTAS=round(Cmat[newly,j],3),
             loo_q=signif(Qmat[newly,j],3), row.names=NULL)
}))

wr <- function(d,f) write.csv(d, file.path(OUT,f), row.names=FALSE)
wr(locus_tab, "ST6_locus_retention_104.csv")
wr(omit_sum,  "ST6_omitted_trait_summary.csv")
wr(loo_meta,  "ST6_genomewide_stability.csv")
wr(stable16,  "ST6_stable_16_of_16.csv")
wr(if(is.null(newly_diag)) data.frame(omitted_trait=character(),gene=character()) else newly_diag, "ST6_newly_significant_diagnostic.csv")
# representative-level detail (rep CTAS/q per run)
repdet <- data.frame(canonical_locus_id=master$canonical_locus_id, representative_gene=master$statistical_representative_gene,
                     full_CTAS=round(master$statistical_representative_CTAS,3), full_q=master$statistical_representative_q)
for(j in seq_along(trait_names)){ repdet[[paste0("q_omit_",trait_names[j])]] <- signif(Qmat[master$statistical_representative_gene,j],3) }
wr(repdet, "ST6_representative_retention.csv")
saveRDS(list(Qmat=Qmat, Cmat=Cmat, ret05=ret05, full=full), file.path(LOO,"loo_results.rds"))

cat(sprintf("\nSUMMARY: mean loci retained/run = %.1f/104 (%.1f%%); stable 16/16 = %d loci\n",
    mean(omit_sum$n_loci_retained_q05), 100*mean(omit_sum$n_loci_retained_q05)/104, nrow(stable16)))
cat(sprintf("largest change (fewest retained): omit %s -> %d loci; smallest change: omit %s -> %d loci\n",
    omit_sum$omitted_trait[1], omit_sum$n_loci_retained_q05[1],
    omit_sum$omitted_trait[nrow(omit_sum)], omit_sum$n_loci_retained_q05[nrow(omit_sum)]))
cat(sprintf("CTAS spearman range: %.3f-%.3f ; gene q05 Jaccard range: %.3f-%.3f\n",
    min(omit_sum$ctas_spearman_full_vs_loo), max(omit_sum$ctas_spearman_full_vs_loo),
    min(omit_sum$gene_q05_jaccard), max(omit_sum$gene_q05_jaccard)))
