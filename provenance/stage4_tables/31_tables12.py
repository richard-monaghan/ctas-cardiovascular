#!/usr/bin/env python3
"""Tables 1 and 2 (A4 portrait, shared table_style). The protected paper-folder files are swapped
(Tables/Table1.pptx holds the method comparison; Tables/Table2.pptx holds the trait panel), so each source
is identified by its CONTENT and written to the correctly named staged file. Protected originals are only read.
Table 1 cell text and footnote are verbatim. Table 2 keeps the comparator rows verbatim; only the CTAS row,
the footnote abbreviations and the legend are updated for the final method."""
import os, sys, textwrap
from pptx import Presentation
from pptx.util import Inches
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import table_style as ts
P = '<CTAS_V2_PAPER_WORKSPACE>/'; O = '<CTAS_WORKSPACE>/outputs/final_ctas_tables_20261021/'

def read(path):
    s = Presentation(path).slides[0]; t = next(sh.table for sh in s.shapes if sh.has_table)
    grid = [[t.cell(r, c).text for c in range(len(t.columns))] for r in range(len(t.rows))]
    foot = next(sh.text_frame.text for sh in s.shapes if sh.has_text_frame and sh.text_frame.text.strip())
    return grid, foot, s.notes_slide.notes_text_frame.text
src = {}
for f in ['Table1.pptx', 'Table2.pptx']:
    g, ft, nt = read(P + 'Tables/' + f)
    kind = 'trait_panel' if g[0][0] == 'Trait' else 'method_comparison' if g[0][0] == 'Method' else None
    assert kind, f; src[kind] = (f, g, ft, nt)
assert src['trait_panel'][3].startswith('Table 1') and src['method_comparison'][3].startswith('Table 2')

def place(slide, text, est_bottom, name):
    f = O + f'work/{name}_footnote_top.txt'
    top = float(open(f).read()) if os.path.exists(f) else est_bottom + ts.FOOTNOTE_GAP_IN
    ts.add_footnote(slide, text, top_in=top)

def new_slide():
    prs = Presentation(); prs.slide_width = Inches(ts.A4_W_IN); prs.slide_height = Inches(ts.A4_H_IN)
    return prs, prs.slides.add_slide(prs.slide_layouts[6])

# ---------------- Table 1: 16-trait panel (verbatim content)
fname, grid, foot, _ = src['trait_panel']
assert len(grid) == 18 and len(grid[0]) == 7 and len(grid) - 2 == 16
prs, slide = new_slide()
t = slide.shapes.add_table(18, 7, Inches(ts.MARGIN_IN), Inches(ts.TABLE_TOP_IN), Inches(ts.ARTWORK_W_IN), Inches(8.5)).table
col_w = [1.51, 0.90, 0.75, 0.75, 0.82, 1.18, 1.18]
for i, w in enumerate(col_w): t.columns[i].width = Inches(w)
for i in range(18): t.rows[i].height = ts.H_GROUP_HEADER if i == 0 else ts.H_SUBHEADER if i == 1 else ts.H_BODY
for r1, c1, r2, c2 in [(0, 2, 0, 4), (0, 5, 0, 6)]: t.cell(r1, c1).merge(t.cell(r2, c2))
for r in range(18):
    for c in range(7):
        cell = t.cell(r, c)
        if r == 0: ts.style_header_cell(cell, grid[r][c], group=True)
        elif r == 1: ts.style_header_cell(cell, grid[r][c], group=False)
        else: ts.style_body_cell(cell, grid[r][c], row_index=r - 2, first_col=(c == 0))
place(slide, foot, ts.TABLE_TOP_IN + 0.32 + 0.30 + 16 * 0.40, 'Table1')
legend1 = """Table 1. Cardiovascular traits included in the CTAS analysis

This table lists the 16 cardiovascular traits that constitute the CTAS trait panel.
All traits are drawn from the Pan-UK Biobank European-ancestry summary statistics.

Column definitions:
- Trait: Full trait name with acronym in parentheses.
- Pan-UKBB phenotype: Exact phenotype identifier (ICD-10 or Phecode).
- Cases, all ancestries: Total case count across all ancestries.
- EUR cases / EUR controls: European-ancestry case and control counts.
- Observed scale, h² (SE) / Liability scale, h² (SE): SNP-based heritability
  estimates with standard errors. Observed-scale reflects raw trait variance;
  liability-scale adjusts for disease prevalence.

Abbreviations: EUR, European ancestry; h², SNP-based heritability;
SE, standard error; Pan-UKBB, Pan-UK Biobank.

Values are unchanged from the signed trait-panel table; nothing was recalculated.
"""
slide.notes_slide.notes_text_frame.text = legend1.strip()
prs.save(O + 'Main Tables/Table1.pptx'); open(O + 'legends/Table1_legend.txt', 'w').write(legend1)
t1grid = grid

# ---------------- Table 2: conceptual method comparison
fname, grid, foot, _ = src['method_comparison']
ctas_r = [i for i, r in enumerate(grid) if r[0] == 'CTAS']; assert ctas_r == [len(grid) - 1]
comparators = [r[0] for r in grid[1:-1]]
grid = [list(r) for r in grid]
grid[-1] = ['CTAS', 'Gene profile across the complete panel',
    'Which genes retain the strongest positive cross-trait evidence after expected panel dependence is discounted?',
    'Trait correlation matrix Σ = 2 sin(πρS/6); ZCA-cor whitening of each gene\'s one-sided trait z-scores; only positive whitened evidence is retained',
    'Gene-level aggregate evidence (CTAS) across the 16-trait panel; empirical P and genome-wide BH q',
    'Ranked dependence-aware gene architecture; no causal-variant, common-effect or causal-gene claim']
old_sigma = 'Σ, 16 × 16 Gaussian-copula correlation matrix of the sixteen one-sided trait z-scores'
assert old_sigma in foot
foot2 = foot.replace(old_sigma, 'Σ, 16 × 16 trait correlation matrix, 2 sin(πρS/6) transform of the Spearman rank correlations ρS of the sixteen trait signals; ZCA-cor, zero-phase component analysis whitening of the correlation matrix (Σ^−1/2); P, empirical P value')
N_ROWS, N_COLS = len(grid), len(grid[0])
prs, slide = new_slide()
t = slide.shapes.add_table(N_ROWS, N_COLS, Inches(ts.MARGIN_IN), Inches(ts.TABLE_TOP_IN), Inches(ts.ARTWORK_W_IN), Inches(6.5)).table
col_w = [0.80, 0.72, 1.10, 1.57, 1.67, 1.22]
for i, w in enumerate(col_w): t.columns[i].width = Inches(w)
for i in range(N_ROWS): t.rows[i].height = ts.H_GROUP_HEADER if i == 0 else ts.H_BODY
for r in range(N_ROWS):
    for c in range(N_COLS):
        cell = t.cell(r, c)
        if r == 0: ts.style_header_cell(cell, grid[r][c], group=True)
        else: ts.style_body_cell(cell, grid[r][c], row_index=r - 1, first_col=(c == 0))
for c in range(N_COLS): ts.set_cell_bg(t.cell(N_ROWS - 1, c), "D9ECF8")
from pptx.util import Pt
for r in range(1, N_ROWS):          # house-style body size (A4_PORTRAIT_HOUSE_STYLE.md section 5: body 9 pt) so the 9-row table and footnote fit on the page
    for c in range(N_COLS):
        for para in t.cell(r, c).text_frame.paragraphs:
            for run in para.runs: run.font.size = Pt(9)
place(slide, foot2, ts.TABLE_TOP_IN + 9.0, 'Table2')
names = ', '.join(comparators[:-1]) + ' and ' + comparators[-1]
legend2 = f"""Table 2. Conceptual comparison of CTAS with established multi-trait genetic methods.

This table compares CTAS with {len(comparators)} established multi-trait genetic approaches ({names}) on the conceptual dimensions that most usefully distinguish them.

Column definitions:
- Method: Name of the method.
- Unit of analysis: The biological unit at which inference is performed.
- Question addressed: The primary scientific question the method answers.
- Role of cross-trait structure: How dependence or shared structure across traits is modelled or used.
- Primary output: The principal statistical product returned to the user.
- Inferential niche: The distinct inferential role the method occupies.

The CTAS row is shaded only to indicate the method introduced in this paper; shading does not imply superiority over the other methods. CTAS identifies aggregate cross-trait gene-level evidence; it does not itself establish shared causal variants, estimate a common effect size or assign causal genes. The trait correlation matrix is Σ = 2 sin(πρS/6), computed from Spearman correlations of the 16 gene-level trait signals, and each gene's one-sided trait z-scores are whitened with Σ^−1/2 (ZCA-cor whitening). Empirical P values are followed by genome-wide Benjamini–Hochberg q-values. Locus-level summaries are collapsed downstream of the gene-level CTAS test.

Abbreviations: SNP, single-nucleotide polymorphism; LD, linkage disequilibrium; GWAS, genome-wide association study; SEM, structural equation modelling; CTAS, Cross-Trait Aggregation Score; Σ, 16 × 16 trait correlation matrix; ZCA-cor, zero-phase component analysis whitening of the correlation matrix; BH, Benjamini–Hochberg.

Comparator-method wording is unchanged from the approved comparison table; only the CTAS row, the abbreviations and this legend were updated for the final method.
"""
slide.notes_slide.notes_text_frame.text = legend2.strip()
prs.save(O + 'Main Tables/Table2.pptx'); open(O + 'legends/Table2_legend.txt', 'w').write(legend2)
# verification echo
for nm, g in [('Table1', t1grid), ('Table2', grid)]:
    s = Presentation(O + f'Main Tables/{nm}.pptx').slides[0]; tt = next(sh.table for sh in s.shapes if sh.has_table)
    print(nm, 'rows', len(tt.rows), 'cells identical to intended grid:', all(tt.cell(r, c).text == g[r][c] for r in range(len(g)) for c in range(len(g[0]))), 'sources:', src['trait_panel'][0] if nm == 'Table1' else src['method_comparison'][0])
print('comparators', len(comparators), comparators)
