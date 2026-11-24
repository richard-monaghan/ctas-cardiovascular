#!/bin/bash
# MAGMA annotation (window 0,0) from corrected GRCh37 SNP locations. Usage: 12_annotate.sh EP
cd ${CTAS_WORKSPACE:?set CTAS_WORKSPACE}/target_pipeline/finngen_corrected_rerun_01
EP=$1
${CTAS_WORKSPACE:?set CTAS_WORKSPACE}/MAGMA_Mac/magma --annotate window=0,0 --snp-loc input_lifted/$EP.magma37.txt --gene-loc ${REFERENCE_DIR:?set REFERENCE_DIR}/reference/NCBI37.3.gene.loc --out annot/${EP}_annot > logs/annot_$EP.stdout 2>&1
