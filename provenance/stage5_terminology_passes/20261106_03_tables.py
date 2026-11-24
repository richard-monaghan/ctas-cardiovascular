import re, csv, copy, openpyxl
from pptx import Presentation
from pptx.oxml.ns import qn
from lxml import etree
N = '<CTAS_WORKSPACE>/outputs/ctas_final_unified_20261106/'
log = []

# ---------------- Table 1
p = Presentation(N + 'tables/Table1_061126_final.pptx')
n1 = 0
for sh in p.slides[0].shapes:
    if sh.has_table:
        t = sh.table
        for r in range(len(t.rows)):
            for para in t.cell(r, 0).text_frame.paragraphs:
                for run in para.runs:
                    if run.text == 'Peripheral artery disease (PAD)':
                        run.text = 'Other peripheral'
                        br = etree.SubElement(run._r.getparent(), qn('a:br')); run._r.addnext(br)
                        r2 = copy.deepcopy(run._r); br.addnext(r2)
                        r2.find(qn('a:t')).text = 'vascular disease (I73)'
                        log.append(('T1', 'pptx-cell', 'row %d col 0' % r, 'Peripheral artery disease (PAD)', 'Other peripheral vascular disease (I73) [two lines]'))
                        n1 += 1
assert n1 == 1
p.save(N + 'tables/Table1_061126_final.pptx')

# ---------------- Table 3 pptx
OLDF = 'HTN, hypertension; I80,'
p = Presentation(N + 'tables/Table3_061126_final.pptx')
def fix(t):
    o = t
    t = t.replace('ten prespecified FinnGen R12 endpoints', 'ten FinnGen R12 endpoints')
    t = t.replace('PAD, peripheral artery disease; ', '')  # footnote entry removed here, I73 inserted below
    if 'PAD, peripheral artery disease' in o and 'Ath, atherosclerosis; PAD' in o:  # notes order variant
        t = o.replace('ten prespecified FinnGen R12 endpoints', 'ten FinnGen R12 endpoints').replace('PAD, peripheral artery disease', 'I73, other peripheral vascular disease')
    elif 'PAD, peripheral artery disease' in o:
        t = t.replace(OLDF, 'HTN, hypertension; I73, other peripheral vascular disease; I80,')
        assert 'I73, other' in t
    t = re.sub(r'\bPAD\b', 'I73', t)
    return t
for si, s in enumerate(p.slides):
    for sh in s.shapes:
        if sh.has_table:
            for r in range(len(sh.table.rows)):
                for c in range(len(sh.table.columns)):
                    for para in sh.table.cell(r, c).text_frame.paragraphs:
                        for run in para.runs:
                            n = fix(run.text)
                            if n != run.text:
                                log.append(('T3', 'pptx-cell', 's%d r%d c%d' % (si + 1, r, c), run.text, n)); run.text = n
        elif sh.has_text_frame:
            for para in sh.text_frame.paragraphs:
                for run in para.runs:
                    n = fix(run.text)
                    if n != run.text:
                        log.append(('T3', 'pptx-text', 's%d %s' % (si + 1, sh.name), run.text[-260:], n[-260:])); run.text = n
    if s.has_notes_slide:
        for para in s.notes_slide.notes_text_frame.paragraphs:
            for run in para.runs:
                n = fix(run.text)
                if n != run.text:
                    log.append(('T3', 'pptx-notes', 's%d' % (si + 1), run.text[-260:], n[-260:])); run.text = n
p.save(N + 'tables/Table3_061126_final.pptx')

# ---------------- Table 3 xlsx
wb = openpyxl.load_workbook(N + 'tables/Table3_061126_final.xlsx')
for ws in wb:
    for row in ws.iter_rows():
        for c in row:
            if isinstance(c.value, str):
                o = c.value; n = o
                if re.search(r'\bPAD\b', n) or 'prespecified FinnGen' in n:
                    if 'PAD, peripheral artery disease' in n:
                        n = n.replace('PAD, peripheral artery disease; ', '').replace('HTN, hypertension; I80,', 'HTN, hypertension; I73, other peripheral vascular disease; I80,')
                        n = n.replace('PAD, peripheral artery disease', 'I73, other peripheral vascular disease') if 'I73, other' not in n else n
                    n = n.replace('ten prespecified FinnGen R12 endpoints', 'ten FinnGen R12 endpoints')
                    n = re.sub(r'\bPAD\b', 'I73', n)
                if n != o:
                    log.append(('T3', 'xlsx', ws.title + '!' + c.coordinate, o[-200:], n[-200:])); c.value = n
rd = wb['README']; r = rd.max_row + 1
rd.cell(row=r, column=1).value = 'Label and wording update (6 Nov 2026)'
rd.cell(row=r, column=2).value = ('The Pan-UKBB trait with internal key PAD is ICD-10 I73, other peripheral vascular diseases (Raynaud syndrome, thromboangiitis obliterans, other specified and unspecified peripheral vascular disease); it is labelled I73 in this table. '
                                  'The FinnGen columns refer to ten related FinnGen R12 endpoints grouped into six clinical families; this is related-endpoint external support, not exact replication.')
wb.save(N + 'tables/Table3_061126_final.xlsx')
csv.writer(open(N + 'tables/TABLE_EDIT_LOG.tsv', 'w'), delimiter='\t', lineterminator='\n').writerows([('table', 'type', 'location', 'old', 'new')] + [tuple(map(str, l)) for l in log])
print(len(log), [l[1] for l in log].count('xlsx'))
