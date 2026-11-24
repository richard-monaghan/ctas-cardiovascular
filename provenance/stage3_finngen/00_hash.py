import hashlib,os,sys,glob,csv
R='<CTAS_WORKSPACE>/'
F=[]
F+=glob.glob(R+'CTAS_FinnGen/data/finngen_gwas/finngen_R12_I9_*.gz')
F+=[R+'CTAS_FinnGen/results/magma_results/%s%s'%(e,s) for e in ['AF','Angina','CAD','Cardiomyopathy','HF','Hypertension','MI','Stroke','TIA','VTE'] for s in ['.genes.out','_annot.genes.annot','.log']]
F+=[R+'CTAS_FinnGen/data/magma_formatted/%s.magma.txt'%e for e in ['AF','Angina','CAD','Cardiomyopathy','HF','Hypertension','MI','Stroke','TIA','VTE']]
for d in ['outputs/final_ctas_finngen_20261021','outputs/final_c_primary_freeze_20261021','outputs/final_ctas_tables_20261021','outputs/table3_final_20261028','panel_v2_locus_lock','target_pipeline/apoe_finngen_forensic_01']:
    for dp,_,fs in os.walk(R+d):
        for f in fs:
            if f!='.DS_Store' and '__pycache__' not in dp: F.append(os.path.join(dp,f))
F+=[R+'CTAS_manuscript_251026.docx',R+'MAGMA_Mac/magma']
F+=glob.glob(R+'ctas_explorer_v2/data/*')+glob.glob(R+'ctas_explorer_v2/*.html')
ref='<REFERENCE_DIR>/reference/'
F+=[ref+x for x in ['NCBI37.3.gene.loc','g1000_eur.bim','g1000_eur.bed','g1000_eur.fam','g1000_eur.synonyms','hg38ToHg19.over.chain']]
out=open(sys.argv[1],'w'); w=csv.writer(out,delimiter='\t'); w.writerow(['path','bytes','sha256'])
for p in sorted(set(F)):
    if not os.path.isfile(p): continue
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<24),b''): h.update(b)
    w.writerow([p,os.path.getsize(p),h.hexdigest()]); out.flush()
