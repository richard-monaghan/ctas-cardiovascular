import openpyxl, csv, re
N = '<CTAS_WORKSPACE>/outputs/ctas_final_unified_20261106/supp_tables/Supplementary Tables/'
NOTE = ('Trait-key note (6 Nov 2026, I73 label)',
        'The internal trait key PAD used in this workbook (column headers, row labels and file names, where present) is retained unchanged so that joins and code are unaffected. '
        'It denotes the complete three-character Pan-UK Biobank ICD-10 category I73, other peripheral vascular diseases (Raynaud syndrome, thromboangiitis obliterans [Buerger], other specified and unspecified peripheral vascular disease); '
        'it is not a peripheral artery disease-specific phenotype and is labelled I73 in paper-facing text.')
log = []
def readme(ws):
    r = ws.max_row + 1
    ws.cell(row=r, column=1).value = NOTE[0]; ws.cell(row=r, column=2).value = NOTE[1]
    for col in (1, 2):
        src = ws.cell(row=r - 1, column=col); dst = ws.cell(row=r, column=col)
        if src.has_style: dst._style = src._style
def run(fn, edits=None):
    wb = openpyxl.load_workbook(N + fn)
    if edits:
        for ws in wb:
            for row in ws.iter_rows():
                for c in row:
                    if isinstance(c.value, str):
                        n = edits(ws.title, c)
                        if n is not None and n != c.value:
                            log.append((fn, ws.title, c.coordinate, c.value, n)); c.value = n
    readme(wb['README']); wb.save(N + fn)
def st6(t, c):
    if t == 'trait_effective_N':
        if c.column == 1 and c.value == 'Peripheral artery disease': return 'Other peripheral vascular disease'
        if c.column == 2 and c.value == 'PAD': return 'I73'
def st8(t, c):
    if t == 'README' and 'removed: MI, angina, PAD,' in c.value:
        return c.value.replace('removed: MI, angina, PAD,', 'removed: MI, angina, I73,')
run('ST6_CTAS_effective_N_sensitivity.xlsx', st6)
run('ST8_CTAS_obstructive_cluster_sensitivity.xlsx', st8)
for fn in ['ST3_CTAS_locus_architecture.xlsx', 'ST7_CTAS_leave_one_trait_out.xlsx', 'ST9_CTAS_colocalisation.xlsx', 'ST12_CTAS_trait_missingness.xlsx']:
    run(fn)
csv.writer(open(N + '../SUPP_EDIT_LOG.tsv', 'w'), delimiter='\t', lineterminator='\n').writerows([('file', 'sheet', 'cell', 'old', 'new')] + [tuple(map(str, l)) for l in log])
print(len(log))
for l in log: print(l[0], l[2], '|', l[3][:60], '->', l[4][:60])
