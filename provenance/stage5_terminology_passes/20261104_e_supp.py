import openpyxl,os,csv,re
O='<CTAS_WORKSPACE>/outputs/'
N=O+'ctas_label_correction_20261104/supp_tables/Supplementary Tables/';os.makedirs(N,exist_ok=True)
S21=O+'final_ctas_tables_20261021/Supplementary Tables/';S4=O+'final_ctas_tables_20261104_families/Supplementary Tables/'
KEYNOTE=('Trait-key note (4 Nov 2026 label correction)','The internal trait keys DVT and VTE used in this workbook (column headers, row labels and file names) are retained unchanged so that joins and code are unaffected. They denote the complete three-character Pan-UK Biobank ICD-10 categories I80 (phlebitis and thrombophlebitis) and I82 (other venous embolism and thrombosis), respectively. The key DVT does not denote a deep-vein-specific phenotype and the key VTE is not the FinnGen venous thromboembolism endpoint. No value was recalculated.')
FGNOTE=('Label note (4 Nov 2026)','In this workbook VTE (I9_VTE) is the separately defined FinnGen venous thromboembolism endpoint. It is distinct from the Pan-UK Biobank ICD-10 I82 trait (other venous embolism and thrombosis), and the Pan-UK Biobank ICD-10 I80 trait is not a deep-vein-specific phenotype. No value was recalculated.')
log=[]
def readme(ws,note):
    r=ws.max_row+1;ws.cell(row=r,column=1).value=note[0];ws.cell(row=r,column=2).value=note[1]
    for col in (1,2):
        src=ws.cell(row=r-1,column=col);dst=ws.cell(row=r,column=col)
        if src.has_style: dst._style=src._style
def run(fn,src,note,edits=None):
    wb=openpyxl.load_workbook(src+fn)
    for ws in wb:
        if edits:
            for row in ws.iter_rows():
                for c in row:
                    if isinstance(c.value,str):
                        n=edits(ws.title,c)
                        if n is not None and n!=c.value: log.append((fn,ws.title,c.coordinate,c.value,n));c.value=n
    readme(wb['README'],note);wb.save(N+fn)
def st6(t,c):
    v=c.value
    if t=='trait_effective_N':
        return {'Deep vein thrombosis':'Phlebitis/thrombophlebitis','Venous thromboembolism':'Venous embolism/thrombosis'}.get(v) if c.column==1 else {'DVT':'I80','VTE':'I82'}.get(v) if c.column==2 else None
def st8(t,c):
    v=c.value
    if 'DVT' in v and 'VTE' in v and 'mitral valve' in v: return v.replace('DVT, HF','I80, HF').replace('mitral valve, VTE','mitral valve, I82')
run('ST6_CTAS_effective_N_sensitivity.xlsx',S21,KEYNOTE,st6)
run('ST8_CTAS_obstructive_cluster_sensitivity.xlsx',S21,KEYNOTE,st8)
for fn in ['ST3_CTAS_locus_architecture.xlsx','ST7_CTAS_leave_one_trait_out.xlsx','ST9_CTAS_colocalisation.xlsx','ST12_CTAS_trait_missingness.xlsx']: run(fn,S21,KEYNOTE)
run('ST11_CTAS_FinnGen_support.xlsx',S4,FGNOTE)
csv.writer(open(N+'../SUPP_EDIT_LOG.tsv','w'),delimiter='\t',lineterminator='\n').writerows([('file','sheet','cell','old','new')]+log)
print(len(log));[print(l) for l in log]
