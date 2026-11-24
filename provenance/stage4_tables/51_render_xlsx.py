#!/usr/bin/env python3
"""Review renders of each staged workbook for visual QA: README in full, every sheet with <= 60 rows in full, and the
header, first rows and last rows of long ledgers. Also performs structural checks (duplicate sheet names, hidden sheets,
formulas, error values, blank required identifiers, duplicated primary keys). Writes work/xlsx_review/<ST>.png and work/xlsx_struct.tsv."""
import os, glob, re
import openpyxl, pandas as pd
from PIL import Image, ImageDraw, ImageFont
O = '<CTAS_WORKSPACE>/outputs/final_ctas_tables_20261021/'; os.makedirs(O + 'work/xlsx_review', exist_ok=True)
font = ImageFont.truetype('/System/Library/Fonts/Menlo.ttc', 13); bold = ImageFont.truetype('/System/Library/Fonts/Menlo.ttc', 14)
KEYS = {'gene_architecture': 'gene', 'gene_null_results': 'gene', 'gene_membership': 'gene', 'locus_summary': 'canonical_locus_id', 'locus_colocalisation': 'canonical_locus_id',
        'four_study_matrix': 'canonical_locus_id', 'gene_ledger': 'gene', 'locus_retention_91': 'canonical_locus_id', 'reduced_gene_results': 'gene', 'per_gene_availability': 'gene', 'strict_core_21': 'canonical_locus_id'}
struct = []
for f in sorted(glob.glob(O + 'Supplementary Tables/*.xlsx')):
    name = os.path.basename(f).split('_')[0]; wb = openpyxl.load_workbook(f)
    lines = [(f'===== {os.path.basename(f)}  sheets: {", ".join(wb.sheetnames)}', True)]
    dup = len(wb.sheetnames) != len(set(wb.sheetnames)); hidden = [ws.title for ws in wb.worksheets if ws.sheet_state != 'visible']; formulas = errs = 0; blank_ids = dupkeys = 0
    for ws in wb.worksheets:
        rows = list(ws.iter_rows(values_only=True))
        for r in rows:
            for v in r:
                if isinstance(v, str) and v.startswith('='): formulas += 1
                if isinstance(v, str) and re.match(r'^#(REF|VALUE|N/A|DIV/0|NAME)', v): errs += 1
        if ws.title in KEYS and ws.title != 'README':
            hdr = list(rows[1]); k = hdr.index(KEYS[ws.title]) if KEYS[ws.title] in hdr else None
            if k is not None:
                ids = [r[k] for r in rows[2:]]
                blank_ids += sum(1 for x in ids if x in (None, '')); 
                if ws.title not in ('gene_membership',) and not (name == 'ST7' and False): dupkeys += len(ids) - len(set(ids))
        lines.append((f'--- [{ws.title}] {ws.max_row - (0 if ws.title == "README" else 2)} data rows x {ws.max_column} cols', True))
        show = rows if (ws.title == 'README' or len(rows) <= 62) else rows[:7] + [('...',)] + rows[-3:]
        for r in show:
            if ws.title == 'README':
                import textwrap
                k, v = (list(r) + [None, None])[:2]; txt = f'{k or ""} | {v or ""}' if v else str(k or '')
                for j, w in enumerate(textwrap.wrap(txt, 200) or ['']): lines.append(((' ' * 4 if j else '') + w, False))
                continue
            cells = ['' if v is None else (f'{v:.6g}' if isinstance(v, float) else str(v)) for v in r][:12]
            lines.append((' | '.join(c[:44] for c in cells)[:330], False))
    struct.append(dict(workbook=os.path.basename(f), sheets=len(wb.sheetnames), duplicate_sheet_names=dup, hidden_sheets=';'.join(hidden), formula_cells=formulas, error_values=errs,
                       blank_required_ids=blank_ids, duplicated_primary_keys=dupkeys, readme_first=wb.sheetnames[0] == 'README'))
    H = 20 * len(lines) + 20; img = Image.new('RGB', (2600, H), 'white'); d = ImageDraw.Draw(img); y = 10
    for t, b_ in lines: d.text((10, y), t, fill=(0, 60, 140) if b_ else (0, 0, 0), font=bold if b_ else font); y += 20
    img.save(O + f'work/xlsx_review/{name}.png')
pd.DataFrame(struct).to_csv(O + 'work/xlsx_struct.tsv', sep='\t', index=False); print(pd.DataFrame(struct).to_string())
