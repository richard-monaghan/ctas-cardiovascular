#!/usr/bin/env python3
"""Literature stage 0 (written for the public release, P6; validated byte-identical against the 30 Apr 2026 extracts).
Transforms the comparator publications' supplementary tables (downloaded by the user; see inputs/literature/
COMPARATOR_SOURCE_MANIFEST.tsv) into the processed tables read by 01_build_clean_inputs.py:
  Qiao 2025  Supplementary Data 17, section a) only (38 CVD-only MTAG loci)      -> processed/qiao_2025_sd17a_loci.tsv
  Liu 2023   Table S2 pleiotropic gene records (264) -> 180 unique genes          -> processed/liu_2023_pleiotropic_genes.tsv
  Liu 2023   Table S3 pleiotropic SNP records (690, all rows)                     -> processed/liu_2023_pleiotropic_snps.tsv
  Foley 2021 candidate SNPs -> GRCh37 positions from the 1000 Genomes EUR BIM (optional; packaged result otherwise)
Usage: python3 analysis/literature/00_transform_sources.py <sources_dir> <out_dir> [g1000_eur.bim]"""
import sys, os, collections, openpyxl, csv
src, out = sys.argv[1], sys.argv[2]; os.makedirs(os.path.join(out, 'processed'), exist_ok=True)
fmt = lambda v: '' if v is None else str(v)
# Qiao: rows of section a) under the header used by the original extract (later labels are positional; 01 reads chr/start/stop/locus/SNP)
ws = openpyxl.load_workbook(os.path.join(src, 'qiao_2025/41467_2025_62419_MOESM19_ESM.xlsx'), read_only=True, data_only=True)['Supplementary Data 17']
rows = list(ws.iter_rows(values_only=True)); i = [k for k, r in enumerate(rows) if r[0] == 'a)'][0] + 3; sec = []
while rows[i][0] is not None: sec.append(rows[i]); i += 1
hdr = 'locus_id\tchr\tstart\tstop\tintersections\tposterior_prob\tregional_prob\tcandidate_snp\tposterior_explained_by_snp\thyprcoloc_traits\tregion_annotation\tvariant_annotation\tm_value_MTAG_AF\tm_value_MTAG_CAD\tm_value_MTAG_VTE\tm_value_MTAG_Stroke\tm_value_MTAG_HF\tm_value_MTAG_PAD\tcol19\tcol20\n'
open(os.path.join(out, 'processed/qiao_2025_sd17a_loci.tsv'), 'w').write(hdr + ''.join('\t'.join(fmt(c) for c in r[:20]) + '\n' for r in sec))
# Liu genes: union of the two traits of every pairwise record; trait pairs = distinct 'Trait' labels; sorted by (-n_traits, gene)
ws = openpyxl.load_workbook(os.path.join(src, 'liu_2023_ten_cvd/mmc2.xlsx'), read_only=True, data_only=True)['Table S2']
g = collections.OrderedDict()
for t1, t2, tp, gene, *_ in [r for r in list(ws.iter_rows(values_only=True))[2:] if r[3]]:
    d = g.setdefault(gene, {'traits': set(), 'pairs': []}); d['traits'].update([t1, t2])
    if tp not in d['pairs']: d['pairs'].append(tp)
open(os.path.join(out, 'processed/liu_2023_pleiotropic_genes.tsv'), 'w').write('gene\tn_traits\ttraits\tn_trait_pairs\ttrait_pairs\n' + ''.join(
    f"{k}\t{len(v['traits'])}\t{'|'.join(sorted(v['traits']))}\t{len(v['pairs'])}\t{'|'.join(sorted(v['pairs']))}\n" for k, v in sorted(g.items(), key=lambda x: (-len(x[1]['traits']), x[0]))))
# Liu SNPs: every record, selected columns, sheet order
ws = openpyxl.load_workbook(os.path.join(src, 'liu_2023_ten_cvd/mmc3.xlsx'), read_only=True, data_only=True)['Table S3']
open(os.path.join(out, 'processed/liu_2023_pleiotropic_snps.tsv'), 'w').write('trait_1\ttrait_2\tsnp\tchr\tpos\tinc_allele\tdec_allele\n' + ''.join(
    '\t'.join(fmt(r[i]) for i in (0, 1, 3, 4, 5, 6, 7)) + '\n' for r in list(ws.iter_rows(values_only=True))[2:] if r[3]))
# Foley candidate SNP positions (GRCh37) from the 1000 Genomes EUR BIM, first occurrence, sorted by rsID
if len(sys.argv) > 3:
    reg = list(csv.DictReader(open(os.path.join(src, 'foley_2021/foley_2021_chd43_regions.tsv')), delimiter='\t'))
    rs = sorted(set(r['candidate_snp'] for r in reg)); pos = {}
    for l in open(sys.argv[3]):
        p = l.split('\t')
        if p[1] in rs and p[1] not in pos: pos[p[1]] = (p[0], p[3])
    open(os.path.join(out, 'foley_snp_hg19.tsv'), 'w').write('rsID\tchr_hg19\tpos_hg19\n' + ''.join(f'{r}\t{pos[r][0]}\t{pos[r][1]}\n' for r in rs))
print('processed tables written to', out)
