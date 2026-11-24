from pptx import Presentation
from pptx.util import Emu
from pptx.oxml.ns import qn
from lxml import etree
SRC='<CTAS_WORKSPACE>/outputs/final_ctas_tables_20261021/Main Tables/Table1.pptx'
DST='<CTAS_WORKSPACE>/outputs/ctas_label_correction_20261104/table1/Table1_labels_041126.pptx'
p=Presentation(SRC);s=p.slides[0]
for sh in s.shapes:
    if sh.has_table:
        t=sh.table
        for r in range(len(t.rows)):
            c=t.cell(r,0)
            for para in c.text_frame.paragraphs:
                for run in list(para.runs):
                    pair={'Deep vein thrombosis (DVT)':('Phlebitis/','thrombophlebitis (I80)'),'Venous thromboembolism (VTE)':('Venous embolism/','thrombosis (I82)')}.get(run.text)
                    if pair:
                        import copy
                        run.text=pair[0]
                        br=etree.SubElement(run._r.getparent(),qn('a:br'));run._r.addnext(br)
                        r2=copy.deepcopy(run._r);br.addnext(r2)
                        r2.find(qn('a:t')).text=pair[1]
    elif sh.has_text_frame:
        run=sh.text_frame.paragraphs[0].runs[0]
        assert run.text.endswith('summary statistic used.'),run.text[-40:]
        run.text=run.text+' I80 and I82 denote the complete three-character Pan-UK Biobank ICD-10 categories (phlebitis and thrombophlebitis; other venous embolism and thrombosis), including all subcodes, and are distinct from the broader FinnGen venous thromboembolism endpoint (I9_VTE) used for external support.'
        sh.height=Emu(1000000)
p.save(DST);print('saved')
