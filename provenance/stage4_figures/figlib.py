"""Shared figure helpers: house style (Arial, CTAS palette), panel export (vector PDF + 600-dpi PNG),
panel source-data export and number-ledger registration. Every displayed number is registered with its source."""
import os, json, math
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CT = '<CTAS_WORKSPACE>/'; OUT = CT + 'outputs/ctas_final_unified_20261106/work/figout/'
T1 = CT + 'outputs/final_c_primary_freeze_20261021/'; T2 = CT + 'outputs/final_ctas_finngen_20261021/'
T3 = CT + 'outputs/final_ctas_colocalisation_20261021/'; T4 = CT + 'outputs/final_ctas_tables_20261021/'
SIGNED_ST = '<CTAS_V2_PAPER_WORKSPACE>/Supplementary Tables/'
BLUE, SKY, DARK, GREY, LGREY, VERM, INK, AXIS = '#0072B2', '#56B4E9', '#084B7F', '#999999', '#D9D9D9', '#D55E00', '#000000', '#666666'
GREEN, ORANGE, YELLOW, PURPLE = '#009E73', '#E69F00', '#F0E442', '#CC79A7'
plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'axes.titlesize': 9, 'axes.labelsize': 9, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
    'legend.fontsize': 7.5, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': AXIS, 'axes.linewidth': 0.6,
    'xtick.color': INK, 'ytick.color': INK, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
TRAIT_ABBR = {'AF_I48': 'AF', 'Angina': 'Ang', 'Aortic_Aneurysm': 'AoAn', 'Aortic_Valve': 'AoV', 'Atherosclerosis': 'Ath', 'CAD_I25': 'CAD', 'Cardiomyopathy': 'CMP', 'DVT': 'I80',
    'HF_I50': 'HF', 'Hypertension_401': 'HTN', 'Ischemic_Stroke': 'IS', 'MI': 'MI', 'Mitral_Valve': 'MV', 'PAD': 'I73', 'TIA': 'TIA', 'VTE': 'I82'}
FAMILIES = [('coronary', ['CAD_I25', 'MI', 'Angina']), ('pressure', ['Hypertension_401']), ('arrhythmia', ['AF_I48']), ('heart_failure', ['HF_I50', 'Cardiomyopathy']),
    ('valve', ['Aortic_Valve', 'Mitral_Valve']), ('aortic_vascular', ['Aortic_Aneurysm', 'PAD', 'Atherosclerosis']), ('cerebrovascular', ['Ischemic_Stroke', 'TIA']), ('thrombotic', ['DVT', 'VTE'])]
FAM_LABEL = {'coronary': 'Coronary', 'pressure': 'Pressure', 'arrhythmia': 'Arrhythmia', 'heart_failure': 'Heart failure', 'valve': 'Valve', 'aortic_vascular': 'Aortic/vascular',
    'cerebrovascular': 'Cerebrovascular', 'thrombotic': 'Thrombotic'}
TISSUE = {'aorta': 'Aorta', 'heart_AA': 'Heart AA', 'blood': 'Whole blood', 'coronary': 'Coronary', 'heart_LV': 'Heart LV', 'liver': 'Liver', 'kidney_cortex': 'Kidney cortex'}

LEDGER = []; SOURCES = []
def num(fig, panel, quantity, displayed, source, calc, observed, expected=None):
    ok = (expected is None) or (str(observed) == str(expected)) or (isinstance(observed, (int, float)) and isinstance(expected, (int, float)) and math.isclose(observed, expected, rel_tol=1e-12))
    LEDGER.append(dict(figure=fig, panel=panel, quantity=quantity, displayed=displayed, source_file=source, source_columns_or_calculation=calc,
                       expected='' if expected is None else expected, observed=observed, status='PASS' if ok else 'FAIL'))
def src(fig, panel, quantity, source, calc, rebuild):
    SOURCES.append(dict(figure=fig, panel=panel, plotted_quantity=quantity, source_file=source, source_columns_or_calculation=calc, rebuild_status=rebuild))
def save(fig_obj, name, data=None):
    fig_obj.savefig(OUT + f'panels/{name}.png', dpi=600, facecolor='white'); fig_obj.savefig(OUT + f'panels/{name}.pdf', facecolor='white'); plt.close(fig_obj)
    if data is not None: data.to_csv(OUT + f'source_data/panel_{name}.tsv', sep='\t', index=False)
def dump(tag):
    pd.DataFrame(LEDGER).to_csv(OUT + f'work/ledger_{tag}.tsv', sep='\t', index=False); pd.DataFrame(SOURCES).to_csv(OUT + f'work/sources_{tag}.tsv', sep='\t', index=False)
def pct(a, b): return 100 * a / b
