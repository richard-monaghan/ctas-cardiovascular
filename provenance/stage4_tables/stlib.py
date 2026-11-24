"""Shared helpers for the supplementary-table builders: workbook writing (README first, frozen header,
filters, widths, full-precision values), number-ledger registration and common inputs."""
import os, csv, math
import numpy as np, pandas as pd
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

CT = '<CTAS_WORKSPACE>/'; OUT = CT + 'outputs/final_ctas_tables_20261021/'
T1 = CT + 'outputs/final_c_primary_freeze_20261021/'; T2 = CT + 'outputs/final_ctas_finngen_20261021/'; T3 = CT + 'outputs/final_ctas_colocalisation_20261021/'
PAPER = '<CTAS_V2_PAPER_WORKSPACE>/'; SRC = PAPER + '_source_assets/'
MATRIX = CT + 'analytic_null_validation/inputs/UKBB_16trait_gene_matrix.tsv'
HF = Font(bold=True, color='FFFFFF'); HFILL = PatternFill('solid', fgColor='0072B2')
SENS_NOTE = ('This is a precomputed sensitivity analysis. Its runs used their established raw-Spearman working trait correlation matrices and '
             '10^6-draw sensitivity nulls; they were not rerun under the final primary analysis (transformed correlation matrix, ZCA-cor whitening, '
             '10^9-draw null). Retention and membership summaries are recomputed here against the final 223/99/38 CTAS genes and 91/46/21 CTAS loci.')

def gene_table():
    g = pd.read_csv(T1 + 'final_C_gene_results_seed42.tsv', sep='\t', dtype={'canonical_locus_id_or_mapping_key': str})
    for k, t in [('in_q05', .05), ('in_q01', .01), ('in_q001', .001)]: g[k] = g.q_BH < t
    return g

def locus_table():
    l = pd.read_csv(T1 + 'final_C_locus_results_seed42.tsv', sep='\t')
    l['in_q05'] = l.n_q05_genes > 0; l['in_q01'] = l.n_q01_genes > 0; l['in_q001'] = l.n_q001_genes > 0
    return l

def final_members(g):
    """locus -> list of final q<0.05 member genes"""
    s = g[g.in_q05]
    return {k: list(v) for k, v in s.groupby('canonical_locus_id_or_mapping_key').gene}

class Book:
    def __init__(self, asset):
        self.asset = asset; self.wb = openpyxl.Workbook(); self.wb.remove(self.wb.active); self.ledger = []
    def readme(self, title, pairs):
        ws = self.wb.create_sheet('README', 0); ws.cell(1, 1, title).font = Font(bold=True, size=12)
        for i, (k, v) in enumerate(pairs, 3):
            ws.cell(i, 1, k).font = Font(bold=True); c = ws.cell(i, 2, v); c.alignment = Alignment(wrap_text=True, vertical='top')
        ws.column_dimensions['A'].width = 28; ws.column_dimensions['B'].width = 130
    def sheet(self, name, df, description):
        assert name not in self.wb.sheetnames and len(name) <= 31, name
        ws = self.wb.create_sheet(name); ws.cell(1, 1, description).alignment = Alignment(wrap_text=False)
        ws.cell(1, 1).font = Font(italic=True, color='444444')
        cols = list(df.columns)
        for j, h in enumerate(cols, 1):
            c = ws.cell(2, j, str(h)); c.font = HF; c.fill = HFILL; c.alignment = Alignment(wrap_text=True, vertical='top')
        for i, row in enumerate(df.itertuples(index=False), 3):
            for j, v in enumerate(row, 1):
                if v is None or (isinstance(v, float) and math.isnan(v)): continue
                if isinstance(v, (np.bool_, bool)): v = 'TRUE' if v else 'FALSE'
                elif isinstance(v, np.integer): v = int(v)
                elif isinstance(v, np.floating): v = float(v)
                ws.cell(i, j, v)
        ws.freeze_panes = 'B3' if len(cols) > 1 else 'A3'
        if len(df): ws.auto_filter.ref = f'A2:{get_column_letter(len(cols))}{len(df) + 2}'
        for j, h in enumerate(cols, 1):
            sample = [str(x) for x in df.iloc[:200, j - 1].tolist()] if len(df) else []
            w = max([len(str(h))] + [len(s) for s in sample]) + 2
            ws.column_dimensions[get_column_letter(j)].width = max(9, min(55, w))
        return ws
    def num(self, location, metric, displayed, full=None, denominator='', unit='', source=''):
        self.ledger.append(dict(asset=self.asset, location=location, metric=metric, displayed_value=displayed,
                                full_precision_value='' if full is None else full, denominator=denominator, unit=unit, authoritative_source=source))
    def save(self, fname):
        self.wb.save(OUT + 'Supplementary Tables/' + fname)
        pd.DataFrame(self.ledger).to_csv(OUT + f'work/ledger_{self.asset}.tsv', sep='\t', index=False)
        print('saved', fname, self.wb.sheetnames)

def headline_pairs():
    return [('Headline (final CTAS)', '19,080 tested genes; 223 / 99 / 38 CTAS genes and 91 / 46 / 21 CTAS loci at q<0.05 / q<0.01 / q<0.001 (nested sets). '
             'Primary analysis: transformed trait correlation matrix Sigma = 2 sin(pi rho_S / 6), ZCA-cor whitening, empirical P from 10^9 null draws, genome-wide Benjamini-Hochberg q.')]
