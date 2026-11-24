#!/usr/bin/env python3
# CTAS V2 public release: extract of target_pipeline/finngen_corrected_rerun_01/code/21_compare.py (lines 3, 7 and 72-74,
# "ledger deliverable"). Re-reads the corrected gene ledger with the same reader and writes
# CORRECTED_FINNGEN_GENE_TRAIT_LEDGER.tsv, which 100_families.py consumes. The two comparison columns
# old_n_endpoints_p01_INVALID / old_supported_ge3_INVALID came from the superseded GRCh38-mismatch run (not redistributed)
# and are omitted; 100_families.py does not read them. The read/write round trip is kept so float formatting is identical.
import os, pandas as pd
R = os.path.abspath(os.environ.get('CTAS_REPO', '.')) + '/'; A = R + 'work/finngen/analysis/'
rd = lambda p: pd.read_csv(p, sep='\t', keep_default_na=False, na_values=[''], dtype={'canonical_locus_id_or_mapping_key': str})
L = rd(A + 'ctas_gene_finngen_ledger.tsv')
L.to_csv(A + 'CORRECTED_FINNGEN_GENE_TRAIT_LEDGER.tsv', sep='\t', index=False)
print('ledger rows', len(L))
