import zipfile,re,os,json,csv
PAT=re.compile(r'\bPAD\b|peripheral\s+art|peripheral\s+vasc|\bI73\b|peripheral\s+vascular|PVD',re.I)
def txt_xml(b):
    s=b.decode('utf8','ignore'); s=re.sub(r'</a:p>|</w:p>|</row>|</si>','\n',s); return re.sub(r'<[^>]+>','',s)
def units(path):
    """yield (location, text) units from any asset file"""
    ext=os.path.splitext(path)[1].lower()
    if ext in('.docx','.pptx','.xlsx'):
        z=zipfile.ZipFile(path)
        for n in z.namelist():
            if n.endswith('.xml') and ('word/document' in n or 'slides/slide' in n or 'notesSlide' in n or 'sharedStrings' in n or 'worksheets/sheet' in n or 'footnotes' in n or 'comments' in n or 'charts/chart' in n or 'header' in n or 'footer' in n):
                yield (n, txt_xml(z.read(n)))
    elif ext in('.txt','.tsv','.csv','.md','.json','.js','.html','.css','.py','.sh','.R'):
        t=open(path,encoding='utf8',errors='ignore').read()
        for i,l in enumerate(t.split('\n'),1): yield ('line %d'%i,l)
def scan(path):
    out=[]
    for loc,t in units(path):
        for m in PAT.finditer(t):
            out.append((loc,m.group(0),t[max(0,m.start()-70):m.end()+70].replace('\n',' ').replace('\t',' ')))
    return out
