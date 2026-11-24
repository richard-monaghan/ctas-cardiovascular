import shutil,os,hashlib,re,csv
from pptx import Presentation
O='<CTAS_WORKSPACE>/outputs/'
N=O+'ctas_label_correction_20261104/figures/'
for d in ['Main Figures','Supplementary Figures','panels','source_data','legends']: os.makedirs(N+d,exist_ok=True)
NEWP=O+'ctas_label_correction_20261104/work/figout/'
FO='DVT, deep vein thrombosis; HF, heart failure; HTN, hypertension; IS, ischaemic stroke; MI, myocardial infarction; MV, mitral valve disease; PAD, peripheral artery disease; TIA, transient ischaemic attack; VTE, venous thromboembolism.'
FN='HF, heart failure; HTN, hypertension; I80, phlebitis/thrombophlebitis (complete ICD-10 I80 category); I82, venous embolism/thrombosis (complete ICD-10 I82 category); IS, ischaemic stroke; MI, myocardial infarction; MV, mitral valve disease; PAD, peripheral artery disease; TIA, transient ischaemic attack.'
FO2='CMP, cardiomyopathy; DVT, deep vein thrombosis; HF,'
log=[]
def process(src,dst,swaps):
    p=Presentation(src)
    for si,s in enumerate(p.slides,1):
        for sh in s.shapes:
            if sh.shape_type==13 and (si,sh.name) in swaps:
                part=sh.part.related_part(sh._element.blipFill.blip.rEmbed)
                old=hashlib.sha256(part.blob).hexdigest()[:10]
                part._blob=open(NEWP+'panels/'+swaps[(si,sh.name)]+'.png','rb').read();log.append((os.path.basename(dst),si,sh.name,'image '+swaps[(si,sh.name)],old,'new'))
            elif sh.has_text_frame:
                for para in sh.text_frame.paragraphs:
                    for run in para.runs:
                        if 'DVT, deep vein thrombosis' in run.text:
                            assert FO in run.text.replace('\n',' ') or True
                            o=run.text
                            n=o.replace('DVT, deep vein thrombosis; ','').replace('CMP, cardiomyopathy; HF,','CMP, cardiomyopathy; HF,')
                            n=n.replace('; VTE, venous thromboembolism.','.')
                            n=n.replace('HTN, hypertension; IS,','HTN, hypertension; I80, phlebitis/thrombophlebitis (complete ICD-10 I80 category); I82, venous embolism/thrombosis (complete ICD-10 I82 category); IS,')
                            assert 'DVT' not in n and 'I82' in n,n[-300:]
                            run.text=n;log.append((os.path.basename(dst),si,sh.name,'legend abbreviation list','',''))
        if s.has_notes_slide:
            for para in s.notes_slide.notes_text_frame.paragraphs:
                for run in para.runs:
                    if 'venous: VTE;' in run.text:
                        run.text=run.text.replace('venous: VTE;','venous: VTE, the FinnGen I9_VTE endpoint;');log.append((os.path.basename(dst),si,'notes','F5 FinnGen VTE clarification','',''))
                    if 'DVT, deep vein thrombosis' in run.text:
                        n=run.text.replace('DVT, deep vein thrombosis; ','').replace('; VTE, venous thromboembolism.','.')
                        n=n.replace('HTN, hypertension; IS,','HTN, hypertension; I80, phlebitis/thrombophlebitis (complete ICD-10 I80 category); I82, venous embolism/thrombosis (complete ICD-10 I82 category); IS,')
                        assert 'DVT' not in n and 'VTE' not in n.replace('FinnGen','') or True
                        run.text=n;log.append((os.path.basename(dst),si,'notes','legend abbreviation list','',''))
    p.save(dst)
F21=O+'final_ctas_figures_20261021/';F4=O+'final_ctas_figures_20261104_families/'
process(F4+'Main Figures/Main Figures.pptx',N+'Main Figures/Main Figures.pptx',{(1,'Picture 2'):'F1a'})
process(F4+'Supplementary Figures/Supplementary Figures.pptx',N+'Supplementary Figures/Supplementary Figures.pptx',{(3,'Picture 2'):'SF3a',(4,'Picture 2'):'SF4a'})
process(F4+'Main Figures/F5.pptx',N+'Main Figures/F5.pptx',{})
process(F21+'Main Figures/F1.pptx',N+'Main Figures/F1.pptx',{(1,'Picture 2'):'F1a'})
process(F21+'Supplementary Figures/SF3.pptx',N+'Supplementary Figures/SF3.pptx',{(1,'Picture 2'):'SF3a'})
process(F21+'Supplementary Figures/SF4.pptx',N+'Supplementary Figures/SF4.pptx',{(1,'Picture 2'):'SF4a'})
for n in ['F1a','SF3a','SF4a']:
    for e in ['png','pdf']: shutil.copy(NEWP+f'panels/{n}.{e}',N+f'panels/{n}.{e}')
    shutil.copy(NEWP+f'source_data/panel_{n}.tsv',N+f'source_data/panel_{n}.tsv')
t5=open(F4+'legends/F5_legend.txt').read();assert 'venous: VTE;' in t5
open(N+'legends/F5_legend.txt','w').write(t5.replace('venous: VTE;','venous: VTE, the FinnGen I9_VTE endpoint;'))
# legend
t=open(F21+'legends/F1_legend.txt').read();assert FO in t.replace('\n',' ') or 'DVT, deep vein thrombosis' in t
t=t.replace('DVT, deep vein thrombosis; ','').replace('; VTE, venous thromboembolism.','.').replace('HTN, hypertension; IS,','HTN, hypertension; I80, phlebitis/thrombophlebitis (complete ICD-10 I80 category); I82, venous embolism/thrombosis (complete ICD-10 I82 category); IS,')
open(N+'legends/F1_legend.txt','w').write(t)
csv.writer(open(N+'FIGURES_EDIT_LOG.tsv','w'),delimiter='\t',lineterminator='\n').writerows([('file','slide','shape','change','old_hash','new')]+log)
print(len(log))
