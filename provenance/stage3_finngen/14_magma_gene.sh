#!/bin/bash
# MAGMA gene analysis (SNPwise-mean, default), frozen parameters. Usage: 14_magma_gene.sh EP
cd ${CTAS_WORKSPACE:?set CTAS_WORKSPACE}/target_pipeline/finngen_corrected_rerun_01
EP=$1
${CTAS_WORKSPACE:?set CTAS_WORKSPACE}/MAGMA_Mac/magma --bfile ${REFERENCE_DIR:?set REFERENCE_DIR}/reference/g1000_eur --pval input_lifted/$EP.magma37.txt use=1,4 N=500348 --gene-annot annot/${EP}_annot.genes.annot --out magma_out/$EP > logs/magma_$EP.stdout 2>&1
