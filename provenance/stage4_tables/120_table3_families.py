import openpyxl,copy,pandas as pd
from pptx import Presentation
W='<CTAS_WORKSPACE>/outputs/finngen_corrected_assets_20261104/'; B=W+'backup_endpoint_version_20261104/tables/'; OUTN='Table3_finngen_corrected_041126'
LR=pd.read_csv(W+'families/FAMILY_LOCUS_SUPPORT.tsv',sep='\t',keep_default_na=False).set_index('representative_gene')
def cell(g):
    r=LR.loc[g]; nf=int(float(r.n_families_p01)); ne=int(float(r.n_endpoints_p01)); fams=r.families.replace(';',', ')
    return f'{nf}/6 ({ne}/10)',fams
prs=Presentation(B+OUTN+'.pptx'); ex=None
for s in prs.slides:
    for sh in s.shapes:
        if sh.has_table:
            for r in list(sh.table.rows)[1:]:
                if len(r.cells[5].text_frame.paragraphs)==2: ex=copy.deepcopy(r.cells[5].text_frame.paragraphs[1]._p); break
        if ex is not None: break
chg=[]
for s in prs.slides:
    for sh in s.shapes:
        if sh.has_table:
            hdr=sh.table.rows[0].cells[5]; hdr.text_frame.paragraphs[0].runs[0].text='FinnGen families / 6 (endpoints / 10)'
            for r in list(sh.table.rows)[1:]:
                g=r.cells[0].text.strip(); k,e=cell(g); tf=r.cells[5].text_frame; old=r.cells[5].text.replace('\n',' | ')
                tf.paragraphs[0].runs[0].text=k
                if e:
                    if len(tf.paragraphs)==1: tf._txBody.append(copy.deepcopy(ex))
                    tf.paragraphs[1].runs[0].text=e
                    for extra in tf.paragraphs[1].runs[1:]: extra._r.getparent().remove(extra._r)
                else:
                    for p in tf.paragraphs[1:]: p._p.getparent().remove(p._p)
                chg.append((g,old,r.cells[5].text.replace('\n',' | ')))
        elif sh.has_text_frame and 'FinnGen:' in sh.text_frame.text:
            done=False
            for p in sh.text_frame.paragraphs:
                for ru in p.runs:
                    t=ru.text; a=t.find('FinnGen: number of ten prespecified'); 
                    if a>=0:
                        b=t.find('Prior multi-trait:')
                        new=('FinnGen: number of six clinical families (coronary: Ang, CAD, MI; arrhythmia: AF; myocardial: CMP, HF; pressure: HTN; cerebrovascular: Stroke, TIA; venous: VTE) with gene-level P<0.01 in at least one member endpoint, supported families named, followed in brackets by the number of the ten prespecified FinnGen R12 endpoints with P<0.01; endpoints were grouped so that correlated endpoints are not counted as independent support. All 21 representatives were evaluable and P values are from the corrected GRCh37 analysis (4 November 2026). CTAS was not rerun in FinnGen, so this is independent-cohort external support. ')
                        ru.text=(t[:a]+new+t[b:]).replace('this is convergence, not replication','this is convergence, not independent confirmation'); done=True
            assert done,'footnote not found'
for _s in prs.slides:
    if _s.has_notes_slide:
        _t=_s.notes_slide.notes_text_frame.text; _a=_t.find('- FinnGen (P<0.01 / 10):'); _b=_t.find('- Prior multi-trait:')
        assert _a>=0 and _b>_a
        _t=_t[:_a]+'- FinnGen (families / 6 (endpoints / 10)): the number of six clinical families (coronary: Ang, CAD, MI; arrhythmia: AF; myocardial: CMP, HF; pressure: HTN; cerebrovascular: Stroke, TIA; venous: VTE) with gene-level P<0.01 in at least one member endpoint, supported families named, with the number of the ten prespecified FinnGen R12 endpoints at P<0.01 in brackets. Endpoints are grouped so that correlated endpoints are not counted as independent support. All 21 representatives are evaluable. CTAS was not rerun in FinnGen, so this is independent-cohort external support.\n'+_t[_b:]
        _t=_t.replace('not independent replication','not independent confirmation')
        _s.notes_slide.notes_text_frame.text=_t
prs.save(W+'tables/'+OUTN+'.pptx')
wb=openpyxl.load_workbook(B+OUTN+'.xlsx'); ws=wb['Table 3']; hdr=[c.value for c in ws[3]]; ki=hdr.index('FinnGen P<0.01 (/10)')
ws.cell(row=3,column=ki+1).value='FinnGen families / 6 (endpoints / 10)'
for row in ws.iter_rows(min_row=4):
    g=row[0].value
    if g in LR.index: k,e=cell(g); row[ki].value=k+(': '+e if e else '')
ev=wb['Evidence data']; h=[c.value for c in ev[3]]; c1=len(h)+1
ev.cell(row=3,column=c1).value='fg_n_families_of_6'; ev.cell(row=3,column=c1+1).value='fg_families'
for src_c,dst in ((ev.cell(row=3,column=h.index('fg_n')+1),ev.cell(row=3,column=c1)),(ev.cell(row=3,column=h.index('fg_n')+1),ev.cell(row=3,column=c1+1))):
    dst.font=copy.copy(src_c.font); dst.fill=copy.copy(src_c.fill); dst.alignment=copy.copy(src_c.alignment)
for row in ev.iter_rows(min_row=4):
    g=row[h.index('gene')].value
    if g in LR.index: ev.cell(row=row[0].row,column=c1).value=int(float(LR.at[g,'n_families_p01'])); ev.cell(row=row[0].row,column=c1+1).value=LR.at[g,'families']
n3=int((LR.loc[[x for x in LR.index if x in set(r[0].value for r in ws.iter_rows(min_row=4))]].n_families_p01.astype(float)>=3).sum())
for row in wb['README'].iter_rows():
    if row[0].value=='FinnGen': row[1].value='All 21 representatives evaluable (corrected GRCh37 rerun). Cells show the number of six clinical families (and of ten endpoints) with gene-level P<0.01, families named. Independent-cohort external support, not a re-estimation of CTAS.'
wb['README'].append(['Version note','FinnGen column rebuilt on 2026-11-04 on the clinical-family measure (families: coronary = Angina, CAD, MI; arrhythmia = AF; myocardial = Cardiomyopathy, HF; pressure = Hypertension; cerebrovascular = Stroke, TIA; venous = VTE). Endpoint counts are in brackets and in Evidence data fg_n. The earlier endpoint-only version is superseded.'])
for row in wb['Validation'].iter_rows():
    v=row[0].value
    if isinstance(v,str) and v.startswith('FinnGen 21 evaluable'): row[0].value='FinnGen 21 evaluable (family and endpoint counts from the corrected rerun ledger)'
    if v=='FinnGen fields equal the corrected FinnGen rerun ledger': row[0].value='FinnGen family and endpoint fields equal the corrected rerun ledger'
for _ws in wb:
    for _row in _ws.iter_rows():
        for _c in _row:
            if isinstance(_c.value,str) and 'replication' in _c.value: _c.value=_c.value.replace('Convergence, not replication.','Convergence, not independent confirmation.').replace('not replication','not independent confirmation')
wb.save(W+'tables/'+OUTN+'.xlsx')
for x in chg: print(x)
