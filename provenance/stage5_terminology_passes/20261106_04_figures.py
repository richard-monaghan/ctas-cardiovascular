import shutil, hashlib, csv, re
from pptx import Presentation
N = '<CTAS_WORKSPACE>/outputs/ctas_final_unified_20261106/'
FIG = N + 'figures/'
NEWP = N + 'work/figout/'
log = []
# panels + source data (only the three with trait tick labels)
for n in ['F1a', 'SF3a', 'SF4a']:
    for e in ['png', 'pdf']:
        shutil.copy(NEWP + 'panels/%s.%s' % (n, e), FIG + 'panels/%s.%s' % (n, e))
    shutil.copy(NEWP + 'source_data/panel_%s.tsv' % n, FIG + 'source_data/panel_%s.tsv' % n)
    log.append(('panel', n, 're-plotted with trait tick label PAD -> I73 (figcode figlib.TRAIT_ABBR only)'))

OLD_ENT = 'MV, mitral valve disease; PAD, peripheral artery disease; TIA,'
NEW_ENT = 'MV, mitral valve disease; TIA,'
def legfix(t):
    if OLD_ENT not in t:
        return t
    t = t.replace(OLD_ENT, NEW_ENT)
    assert 'HTN, hypertension; I80,' in t
    return t.replace('HTN, hypertension; I80,', 'HTN, hypertension; I73, other peripheral vascular disease; I80,')

def process(path, swaps):
    p = Presentation(path)
    for si, s in enumerate(p.slides, 1):
        for sh in s.shapes:
            if sh.shape_type == 13 and (si, sh.name) in swaps:
                part = sh.part.related_part(sh._element.blipFill.blip.rEmbed)
                old = hashlib.sha256(part.blob).hexdigest()[:10]
                part._blob = open(NEWP + 'panels/' + swaps[(si, sh.name)] + '.png', 'rb').read()
                log.append((path.split('/')[-1], 's%d %s' % (si, sh.name), 'image swapped %s (old %s)' % (swaps[(si, sh.name)], old)))
            elif sh.has_text_frame:
                for para in sh.text_frame.paragraphs:
                    for run in para.runs:
                        n = legfix(run.text)
                        if n != run.text:
                            log.append((path.split('/')[-1], 's%d %s' % (si, sh.name), 'legend abbreviation list')); run.text = n
        if s.has_notes_slide:
            for para in s.notes_slide.notes_text_frame.paragraphs:
                for run in para.runs:
                    n = legfix(run.text)
                    if n != run.text:
                        log.append((path.split('/')[-1], 's%d notes' % si, 'legend abbreviation list')); run.text = n
    p.save(path)

process(FIG + 'Main Figures/Main Figures.pptx', {(1, 'Picture 2'): 'F1a'})
process(FIG + 'Main Figures/F1.pptx', {(1, 'Picture 2'): 'F1a'})
process(FIG + 'Supplementary Figures/Supplementary Figures.pptx', {(3, 'Picture 2'): 'SF3a', (4, 'Picture 2'): 'SF4a'})
process(FIG + 'Supplementary Figures/SF3.pptx', {(1, 'Picture 2'): 'SF3a'})
process(FIG + 'Supplementary Figures/SF4.pptx', {(1, 'Picture 2'): 'SF4a'})
for lf in ['F1_legend.txt']:
    t = open(FIG + 'legends/' + lf).read()
    n = legfix(t.replace('\n', ' ') if False else t)
    assert n != t
    open(FIG + 'legends/' + lf, 'w').write(n)
    log.append((lf, '', 'legend abbreviation list'))
csv.writer(open(FIG + 'FIGURES_EDIT_LOG.tsv', 'w'), delimiter='\t', lineterminator='\n').writerows([('file', 'location', 'change')] + log)
print(len(log))
