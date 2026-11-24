"""
A4 portrait build of Table 2 using the shared table style (table_style.py).

Cell text and the footnote are extracted VERBATIM from the signed landscape
Tables/Table2.pptx (no scientific drift possible), then re-laid onto an exact
A4 portrait page: 15 mm margins, 180 mm native editable table, identical
visual grammar to canonical Table 1 (header #0072B2 white bold, alternating
#FFFFFF/#F5F5F5 body rows, 0.5 pt #CCCCCC borders, Arial, standard footnote).
The signed table's single header row maps to the group-header style; the
CTAS emphasis row keeps its signed #D9ECF8 shading (a content marker the
footnote/legend explicitly reference).

Column widths are content-balanced for Table 2's four long prose columns
(allowed to differ from Table 1; every other constant is shared).

FOOTNOTE PLACEMENT — measured two-pass (2026-08-07): pass 1 placed the
footnote from a wrap model; the ACTUAL grown table bottom was then measured
from a real render of this PPTX (render_slides.py, the OpenAI presentations
runtime renderer) and the footnote reset exactly 6 mm below it via
place_footnote_measured.py — real bottom 9.393 in, footnote top 9.630 in.
slides_test.py: PASS (no overflow). The model-only matplotlib preview route
was rejected; validation evidence is the real render in
_source_assets/Table2_A4_render/slide-1.png.

Signed landscape predecessor preserved in
_source_assets/superseded_pre_A4_portrait_20260806/Table2/.
"""
import os
import sys
import textwrap
from pptx import Presentation
from pptx.util import Inches

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import table_style as ts

ROOT = "<CTAS_V2_PAPER_WORKSPACE>"
SRC_PPTX = os.path.join(ROOT, "_source_assets", "superseded_pre_A4_portrait_20260806",
                        "Table2", "Table2.pptx")
OUT_PPTX = os.path.join(ROOT, "Tables", "Table2.pptx")
LEGEND = os.path.join(ROOT, "_source_assets", "Table2_legend.txt")

# ---- 1. extract the signed content verbatim (from the ARCHIVED landscape
#         original, so extraction is independent of the file being replaced) --
src = Presentation(SRC_PPTX)
slide0 = src.slides[0]
tbl = next(sh.table for sh in slide0.shapes if sh.has_table)
footnote = next(sh.text_frame.text for sh in slide0.shapes
                if sh.has_text_frame and sh.text_frame.text.strip())

N_ROWS, N_COLS = 7, 6
grid = [[tbl.cell(r, c).text for c in range(N_COLS)] for r in range(N_ROWS)]
CTAS_ROW = 6                      # signed emphasis row (fill #D9ECF8)
assert grid[CTAS_ROW][0] == "CTAS"
MERGES = []                       # signed Table 2 has no merged cells

# ---- 2. rebuild on exact A4 portrait --------------------------------------
prs = Presentation()
prs.slide_width = Inches(ts.A4_W_IN)
prs.slide_height = Inches(ts.A4_H_IN)
slide = prs.slides.add_slide(prs.slide_layouts[6])

shape = slide.shapes.add_table(N_ROWS, N_COLS,
                               Inches(ts.MARGIN_IN), Inches(ts.TABLE_TOP_IN),
                               Inches(ts.ARTWORK_W_IN), Inches(6.5))
table = shape.table

# content-balanced column widths (inches, sum = 180 mm): the four prose
# columns are sized proportional to their longest cell so every column wraps
# to a similar line count (keeps the page compact and row growth even);
# Method / Unit columns sized for their short tokens.
col_w = [0.80, 0.72, 1.10, 1.57, 1.67, 1.22]
assert abs(sum(col_w) - ts.ARTWORK_W_IN) < 0.01, sum(col_w)
for i, w in enumerate(col_w):
    table.columns[i].width = Inches(w)

for i in range(N_ROWS):
    table.rows[i].height = (ts.H_GROUP_HEADER if i == 0 else ts.H_BODY)

for r1, c1, r2, c2 in MERGES:
    table.cell(r1, c1).merge(table.cell(r2, c2))

# style AFTER merges; then re-apply the signed CTAS emphasis fill
for r in range(N_ROWS):
    for c in range(N_COLS):
        cell = table.cell(r, c)
        if r == 0:
            ts.style_header_cell(cell, grid[r][c], group=True)
        else:
            ts.style_body_cell(cell, grid[r][c], row_index=r - 1, first_col=(c == 0))
for c in range(N_COLS):
    ts.set_cell_bg(table.cell(CTAS_ROW, c), "D9ECF8")

# ---- 3. footnote (pass-1 model position; the MEASURED two-pass then ran ---
# pass 1 (this script) places the footnote from a textwrap model (0.0049
# in/char, 0.328 in single-line rows + 0.1667 in per extra line + two-line
# safety). Pass 2 (place_footnote_measured.py, run 2026-08-07 against the
# real render_slides.py output) measured the ACTUAL grown table bottom
# (9.393 in) and reset the footnote exactly 6 mm below it (top = 9.630 in).
# If the table content or column widths are ever edited, re-run BOTH passes.
est_h = 0.0
for r in range(N_ROWS):
    lines = 1
    size = 10.5 if r == 0 else 10
    for c in range(N_COLS):
        cpl = max(1, int((col_w[c] - 0.08) / (0.0049 * size)))
        lines = max(lines, len(textwrap.wrap(grid[r][c], width=cpl,
                                             break_long_words=False)))
    est_h += 0.328 + (lines - 1) * 0.1667
est_h += 2 * 0.1667                # two-line safety margin
tbl_bottom = ts.TABLE_TOP_IN + est_h
ts.add_footnote(slide, footnote, top_in=tbl_bottom + ts.FOOTNOTE_GAP_IN)

with open(LEGEND) as f:
    slide.notes_slide.notes_text_frame.text = f.read().strip()

prs.save(OUT_PPTX)

# ---- 4. verification echo ---------------------------------------------------
chk = Presentation(OUT_PPTX)
print(f"saved: {OUT_PPTX}")
print(f"slide: {chk.slide_width.inches:.4f} x {chk.slide_height.inches:.4f} in")
t2 = next(sh.table for sh in chk.slides[0].shapes if sh.has_table)
same = all(t2.cell(r, c).text == grid[r][c] for r in range(N_ROWS) for c in range(N_COLS))
print(f"native table: {len(t2.rows)}x{len(t2.columns)}; cell text identical to signed Table 2: {same}")
box = next(sh for sh in chk.slides[0].shapes
           if sh.has_text_frame and sh.text_frame.text.strip())
print(f"footnote verbatim: {box.text_frame.text == footnote}")
print(f"footnote top (model+headroom): {box.top.inches:.3f} in")
print(f"CTAS row fill: {t2.cell(CTAS_ROW, 0).fill.fore_color.rgb}")
print(f"notes chars: {len(chk.slides[0].notes_slide.notes_text_frame.text)}")
