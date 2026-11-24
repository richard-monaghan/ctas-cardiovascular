import re,copy,shutil,openpyxl
from pptx import Presentation
S='<CTAS_WORKSPACE>/outputs/table3_final_20261104_families/Table3_final_041126_families'
D='<CTAS_WORKSPACE>/outputs/ctas_label_correction_20261104/table3/Table3_final_041126_families_labels'
ABB_OLD='DVT, deep vein thrombosis; VTE, venous thromboembolism;'
ABB_NEW='I80, phlebitis/thrombophlebitis (complete ICD-10 I80 category); I82, venous embolism/thrombosis (complete ICD-10 I82 category);'
NOTE=' In the FinnGen column, VTE denotes the separately defined FinnGen venous thromboembolism endpoint (I9_VTE), which is distinct from I82.'
def cellfix(t): return re.sub(r'\bVTE\b','I82',re.sub(r'\bDVT\b','I80',t))
log=[]
p=Presentation(S+'.pptx')
for si,s in enumerate(p.slides):
    for sh in s.shapes:
        if sh.has_table:
            t=sh.table
            for r in range(len(t.rows)):
                for c in range(len(t.columns)):
                    for para in t.cell(r,c).text_frame.paragraphs:
                        for run in para.runs:
                            n=cellfix(run.text)
                            if n!=run.text: log.append(('pptx-cell',si,r,c,run.text,n));run.text=n
        elif sh.has_text_frame:
            for para in sh.text_frame.paragraphs:
                for run in para.runs:
                    FO='DVT, deep vein thrombosis; HF, heart failure; HTN, hypertension; IS, ischaemic stroke; MI, myocardial infarction; MV, mitral valve disease; PAD, peripheral artery disease; TIA, transient ischaemic attack; VTE, venous thromboembolism.'
                    FN='HF, heart failure; HTN, hypertension; I80, phlebitis/thrombophlebitis (complete ICD-10 I80 category); I82, venous embolism/thrombosis (complete ICD-10 I82 category); IS, ischaemic stroke; MI, myocardial infarction; MV, mitral valve disease; PAD, peripheral artery disease; TIA, transient ischaemic attack.'
                    if FO in run.text:
                        o=run.text;n=o.replace(FO,FN)
                        n=n.replace('venous: VTE)','venous: VTE, the FinnGen I9_VTE endpoint)') 
                        n=n.rstrip()+NOTE
                        log.append(('pptx-footnote',si,'','',o[-220:],n[-420:]));run.text=n
    if s.has_notes_slide:
        tf=s.notes_slide.notes_text_frame
        for para in tf.paragraphs:
            for run in para.runs:
                if 'DVT, deep vein thrombosis' in run.text:
                    o=run.text;n=o.replace(ABB_OLD,ABB_NEW).rstrip()+NOTE;log.append(('pptx-notes',si,'','',o[-200:],n[-300:]));run.text=n
p.save(D+'.pptx')
wb=openpyxl.load_workbook(S+'.xlsx')
for ws in wb:
    for row in ws.iter_rows():
        for c in row:
            if isinstance(c.value,str):
                n=c.value
                if ws.title in('Table 3','Evidence data') and not re.search(r'(Angina|CAD|AF|HF|MI|TIA)=',n):  # exclude FinnGen endpoint P strings
                    n=cellfix(n)
                if n!=c.value: log.append(('xlsx',ws.title,c.coordinate,'',c.value,n));c.value=n
rd=wb['README'];r=rd.max_row+1
rd.cell(row=r,column=1).value='Label correction (4 Nov 2026)'
rd.cell(row=r,column=2).value='Pan-UKBB traits formerly shown as DVT and VTE are labelled by their ICD-10 categories I80 (phlebitis and thrombophlebitis) and I82 (other venous embolism and thrombosis). FinnGen VTE (I9_VTE) is a separately defined FinnGen endpoint and is unchanged. Label-only change; no value was recalculated.'
wb.save(D+'.xlsx')
import csv
csv.writer(open('<CTAS_WORKSPACE>/outputs/ctas_label_correction_20261104/table3/TABLE3_EDIT_LOG.tsv','w'),delimiter='\t',lineterminator='\n').writerows([('type','a','b','c','old','new')]+[tuple(map(str,l)) for l in log])
print(len(log),[l[0] for l in log].count('xlsx'))
