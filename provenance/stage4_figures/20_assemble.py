#!/usr/bin/env python3
"""Assemble A4-portrait editable PPTX figures: panel artwork (600-dpi PNG) top-aligned under editable panel headings;
legend in speaker notes. Individual F1-F5, SF1-SF5 and compiled main / supplementary decks."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from legends import LEG
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from PIL import Image
O = '<CTAS_WORKSPACE>/outputs/final_ctas_figures_20261021/'
W_IN, H_IN, M, ART = 8.2677, 11.6929, 15 / 25.4, 180 / 25.4
HEAD_H, GAP = 0.26, 0.16
LAYOUT = {  # rows of (panel, heading)
 'F1': [[('F1a', 'a  From correlated traits to calibrated CTAS')], [('F1b', 'b  Genome-wide CTAS evidence across 19,080 genes')]],
 'F2': [[('F2a', 'a  Null calibration (behaviour simulation)'), ('F2b', 'b  Fixed-total-signal breadth behaviour')], [('F2c', 'c  Correlation-model and seed sensitivity')]],
 'F3': [[('F3a', 'a  Cardiovascular-family composition of the 91 CTAS loci')], [('F3b', 'b  Trait-family coverage across the nested locus sets')]],
 'F4': [[('F4a', 'a  Colocalisation evidence across the 91 CTAS loci')], [('F4b', 'b  Tissue distribution of the 31 PP4-supported loci')]],
 'F5': [[('F5a', 'a  Prior-study convergence across the nested locus sets')], [('F5b', 'b  FinnGen external support for CTAS genes')]],
 'SF1': [[('SF1a', 'a  Empirical P against the 10⁹-draw null'), ('SF1b', 'b  Benjamini–Hochberg step-up')], [('SF1c', 'c  Most extreme genes and tail treatment')]],
 'SF2': [[('SF2a', 'a  Direct resolution of significant genes'), ('SF2b', 'b  Breadth × strength behaviour')], [('SF2c', 'c  Fixed-total-signal response'), ('SF2d', 'd  Leave-one-trait-out retention')]],
 'SF3': [[('SF3a', 'a  Signed whitened profiles of the 91 locus representatives')], [('SF3b', 'b  Ordination of family composition')]],
 'SF4': [[('SF4a', 'a  Supported cells by trait and tissue'), ('SF4b', 'b  Shared- vs distinct-variant evidence')], [('SF4c', 'c  Matched-SNP count and PP4 (descriptive)'), ('SF4d', 'd  Locus-level maximum PP4')]],
 'SF5': [[('SF5a', 'a  Four-comparator intersections across the 91 CTAS loci')], [('SF5b', 'b  Gene-level FinnGen support by threshold'), ('SF5c', 'c  FinnGen endpoint counts, q<0.05 genes')]]}
def add_fig(prs, fig):
    s = prs.slides.add_slide(prs.slide_layouts[6]); y = M; bottom = 0
    for row in LAYOUT[fig]:
        x = M; rowh = 0
        for name, head in row:
            im = Image.open(O + f'panels/{name}.png'); w_in = im.width / 600; h_in = im.height / 600
            tb = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w_in), Inches(HEAD_H)); tf = tb.text_frame; tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
            p = tf.paragraphs[0]; r = p.add_run(); r.text = head; r.font.name = 'Arial'; r.font.size = Pt(12); r.font.bold = True; r.font.color.rgb = RGBColor(0, 0, 0)
            s.shapes.add_picture(O + f'panels/{name}.png', Inches(x), Inches(y + HEAD_H), Inches(w_in), Inches(h_in))
            x += w_in + (ART - 2 * 3.45 if len(row) > 1 else 0); rowh = max(rowh, HEAD_H + h_in)
        y += rowh + GAP; bottom = y - GAP
    assert bottom <= H_IN - M + 1e-6, (fig, bottom)
    s.notes_slide.notes_text_frame.text = LEG[fig]
    return bottom
def deck(figs, path):
    prs = Presentation(); prs.slide_width = Inches(W_IN); prs.slide_height = Inches(H_IN)
    bots = {f: add_fig(prs, f) for f in figs}; prs.save(path); return bots
for f in LAYOUT:
    sub = 'Main Figures' if f.startswith('F') else 'Supplementary Figures'
    b = deck([f], O + f'{sub}/{f}.pptx'); open(O + f'legends/{f}_legend.txt', 'w').write(LEG[f] + '\n'); print(f, f'artwork bottom {b[f]:.2f} in of {H_IN:.2f}')
deck(['F1', 'F2', 'F3', 'F4', 'F5'], O + 'Main Figures/Main Figures.pptx')
deck(['SF1', 'SF2', 'SF3', 'SF4', 'SF5'], O + 'Supplementary Figures/Supplementary Figures.pptx')
print('compiled decks written')
