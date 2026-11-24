"""
A4 portrait pilot of Table 1 using the shared table style (table_style.py).

Cell text is extracted VERBATIM from the current signed Tables/Table1.pptx
(no scientific drift possible), then re-laid onto an exact A4 portrait page:
15 mm margins, 180 mm native editable table, dark-blue/light-blue headers,
alternating rows, 0.5 pt #CCCCCC borders, Arial, standard footnote.

Output (pilot path; canonical Tables/Table1.pptx is replaced only in Phase B):
  _source_assets/A4_pilots/Table1_A4.pptx
"""
import os
import sys
from pptx import Presentation
from pptx.util import Inches

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import table_style as ts

ROOT = "<CTAS_V2_PAPER_WORKSPACE>"
SRC_PPTX = os.path.join(ROOT, "Tables", "Table1.pptx")
OUT_DIR = os.path.join(ROOT, "_source_assets", "A4_pilots")
OUT_PPTX = os.path.join(OUT_DIR, "Table1_A4.pptx")
LEGEND = os.path.join(ROOT, "_source_assets", "Table1_legend.txt")

# ---- 1. extract the signed content verbatim -------------------------------
src = Presentation(SRC_PPTX)
slide0 = src.slides[0]
tbl = next(sh.table for sh in slide0.shapes if sh.has_table)
footnote = next(sh.text_frame.text for sh in slide0.shapes
                if sh.has_text_frame and sh.text_frame.text.strip())

grid = [[tbl.cell(r, c).text for c in range(7)] for r in range(18)]

# merge structure of the signed table (group-header row)
MERGES = [(0, 2, 0, 4), (0, 5, 0, 6)]  # (r1, c1, r2, c2)

# ---- 2. rebuild on exact A4 portrait --------------------------------------
prs = Presentation()
prs.slide_width = Inches(ts.A4_W_IN)
prs.slide_height = Inches(ts.A4_H_IN)
slide = prs.slides.add_slide(prs.slide_layouts[6])

n_rows, n_cols = 18, 7
shape = slide.shapes.add_table(n_rows, n_cols,
                               Inches(ts.MARGIN_IN), Inches(ts.TABLE_TOP_IN),
                               Inches(ts.ARTWORK_W_IN), Inches(8.5))
table = shape.table

# Table-1-specific column widths (inches, sum = 180 mm): trait column sized so
# names wrap at most to two lines; h2 columns sized so "0.0362 (0.0032)" stays
# on ONE line at 10 pt Arial.
col_w = [1.51, 0.90, 0.75, 0.75, 0.82, 1.18, 1.18]
assert abs(sum(col_w) - ts.ARTWORK_W_IN) < 0.01, sum(col_w)
for i, w in enumerate(col_w):
    table.columns[i].width = Inches(w)

for i in range(n_rows):
    table.rows[i].height = (ts.H_GROUP_HEADER if i == 0
                            else ts.H_SUBHEADER if i == 1 else ts.H_BODY)

for r1, c1, r2, c2 in MERGES:
    table.cell(r1, c1).merge(table.cell(r2, c2))

# style AFTER merges; touching every cell (origins and covered) guarantees
# fills survive the merge
for r in range(n_rows):
    for c in range(n_cols):
        cell = table.cell(r, c)
        if r == 0:
            ts.style_header_cell(cell, grid[r][c], group=True)
        elif r == 1:
            ts.style_header_cell(cell, grid[r][c], group=False)
        else:
            ts.style_body_cell(cell, grid[r][c], row_index=r - 2, first_col=(c == 0))

# ---- 3. standard footnote + speaker notes ----------------------------------
# estimate rendered (grown) table height: PowerPoint expands rows whose text
# wraps; emulate so the footnote never overlaps the table
import textwrap
est_h = 0.0
for r in range(n_rows):
    lines = 1
    size = 10.5 if r == 0 else 9.5 if r == 1 else 10
    for c in range(n_cols):
        if r == 0 and c in (3, 4, 6):
            continue
        w = sum(col_w[c:c + (3 if c == 2 else 2 if c == 5 else 1)]) if r == 0 else col_w[c]
        cpl = max(1, int((w - 0.10) / (0.0069 * size)))
        lines = max(lines, len(textwrap.wrap(grid[r][c], width=cpl, break_long_words=False)))
    base = 0.32 if r == 0 else 0.30 if r == 1 else 0.32
    est_h += max(base, lines * 0.165 + 0.10)
tbl_bottom = ts.TABLE_TOP_IN + est_h
ts.add_footnote(slide, footnote, top_in=tbl_bottom + ts.FOOTNOTE_GAP_IN)

with open(LEGEND) as f:
    slide.notes_slide.notes_text_frame.text = f.read().strip()

os.makedirs(OUT_DIR, exist_ok=True)
prs.save(OUT_PPTX)

# ---- 4. verification echo ---------------------------------------------------
chk = Presentation(OUT_PPTX)
print(f"saved: {OUT_PPTX}")
print(f"slide: {chk.slide_width.inches:.4f} x {chk.slide_height.inches:.4f} in")
t2 = next(sh.table for sh in chk.slides[0].shapes if sh.has_table)
same = all(t2.cell(r, c).text == grid[r][c] for r in range(18) for c in range(7))
print(f"native table: {len(t2.rows)}x{len(t2.columns)}; cell text identical to signed Table 1: {same}")
print(f"notes chars: {len(chk.slides[0].notes_slide.notes_text_frame.text)}")
