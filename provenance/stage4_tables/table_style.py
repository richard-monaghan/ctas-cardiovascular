"""
CTAS V2 shared A4 table style — single source of truth for Table 1 and Table 2.
See _source_assets/A4_PORTRAIT_HOUSE_STYLE.md section 5.

Both main tables import these constants and helpers so their visual grammar is
identical: page, margins, colours, fonts, row heights, alternating rows,
borders, alignment, footnote typography and spacing. Column widths are set per
table (content differs); everything else is shared.

Tables remain native editable PowerPoint table objects — never images.
No scientific values, row order, terminology or footnote wording may change;
only mechanical wrapping caused by the narrower A4 measure is allowed.
"""
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.oxml import parse_xml

# ---- page ---------------------------------------------------------------
A4_W_IN = 8.2677          # 210 mm
A4_H_IN = 11.6929         # 297 mm
MARGIN_IN = 15 / 25.4     # 15 mm -> 0.5906 in
ARTWORK_W_IN = 180 / 25.4 # 180 mm -> 7.0866 in

# ---- palette (V2 house colours) ------------------------------------------
BLUE_DARK = "0072B2"      # group-header fill (CTAS data blue)
BLUE_LIGHT = "56B4E9"     # subheader fill (restrained secondary)
ROW_ALT = "F5F5F5"        # alternating body rows
ROW_PLAIN = "FFFFFF"
BORDER = "CCCCCC"
GREY_TEXT = RGBColor(0x66, 0x66, 0x66)
INK = RGBColor(0x00, 0x00, 0x00)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)

# ---- typography -----------------------------------------------------------
FONT = "Arial"
SIZE_GROUP_HEADER = Pt(10.5)
SIZE_SUBHEADER = Pt(9.5)
SIZE_BODY = Pt(10)
SIZE_FOOTNOTE = Pt(8.5)

# ---- row heights (minima; rows grow only for mechanical wrapping) -----------
# Natural publication heights (Richard's 22:37 neutral-canvas directive): no
# artificially enlarged rows; the table ends where the content ends.
H_GROUP_HEADER = Inches(0.32)
H_SUBHEADER = Inches(0.30)
H_BODY = Inches(0.32)

# ---- spacing --------------------------------------------------------------
TABLE_TOP_IN = MARGIN_IN           # 15 mm below top margin
FOOTNOTE_GAP_IN = 6 / 25.4         # 6 mm between table and footnote

_BORDER_XML = (
    '<a:tcBorders xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main">'
    '  <a:top w="6350" cap="flat" cmpd="sng" algn="ctr"><a:solidFill><a:srgbClr val="' + BORDER + '"/></a:solidFill></a:top>'
    '  <a:left w="6350" cap="flat" cmpd="sng" algn="ctr"><a:solidFill><a:srgbClr val="' + BORDER + '"/></a:solidFill></a:left>'
    '  <a:bottom w="6350" cap="flat" cmpd="sng" algn="ctr"><a:solidFill><a:srgbClr val="' + BORDER + '"/></a:solidFill></a:bottom>'
    '  <a:right w="6350" cap="flat" cmpd="sng" algn="ctr"><a:solidFill><a:srgbClr val="' + BORDER + '"/></a:solidFill></a:right>'
    '</a:tcBorders>'
)


def set_cell_text(cell, text, bold=False, font_size=SIZE_BODY, align=PP_ALIGN.LEFT, color=INK):
    cell.text = str(text)
    for paragraph in cell.text_frame.paragraphs:
        paragraph.alignment = align
        for run in paragraph.runs:
            run.font.name = FONT
            run.font.size = font_size
            run.font.bold = bold
            run.font.color.rgb = color
    cell.vertical_anchor = MSO_ANCHOR.MIDDLE
    cell.margin_left = Inches(0.04)
    cell.margin_right = Inches(0.04)
    cell.margin_top = Inches(0.02)
    cell.margin_bottom = Inches(0.02)


def set_cell_bg(cell, hex_color):
    """Solid cell fill via the real python-pptx fill API (a previous version
    injected a raw <a:shd> element, which PowerPoint ignores -> fill.type None).
    Call AFTER any cell merges, on both origin and covered cells."""
    cell.fill.solid()
    cell.fill.fore_color.rgb = RGBColor.from_string(hex_color)


def set_cell_borders(cell):
    tc = cell._tc
    tcPr = tc.get_or_add_tcPr()
    for child in list(tcPr):
        if child.tag.endswith('tcBorders'):
            tcPr.remove(child)
    tcPr.append(parse_xml(_BORDER_XML))


def style_body_cell(cell, text, row_index, first_col=False):
    """Alternating fill + standard alignment for a body cell.
    row_index is 0-based within the body (excluding header rows)."""
    align = PP_ALIGN.LEFT if first_col else PP_ALIGN.CENTER
    set_cell_text(cell, text, bold=False, font_size=SIZE_BODY, align=align, color=INK)
    set_cell_bg(cell, ROW_ALT if row_index % 2 == 1 else ROW_PLAIN)
    set_cell_borders(cell)


def style_header_cell(cell, text, group=True):
    """group=True -> dark-blue group header; False -> light-blue subheader."""
    if group:
        set_cell_text(cell, text, bold=True, font_size=SIZE_GROUP_HEADER,
                      align=PP_ALIGN.CENTER, color=WHITE)
        set_cell_bg(cell, BLUE_DARK)
    else:
        set_cell_text(cell, text, bold=True, font_size=SIZE_SUBHEADER,
                      align=PP_ALIGN.CENTER, color=INK)
        set_cell_bg(cell, BLUE_LIGHT)
    set_cell_borders(cell)


def add_footnote(slide, text, top_in, width_in=ARTWORK_W_IN, left_in=MARGIN_IN):
    """Standard table footnote text box (8 pt grey, left, wrapped)."""
    box = slide.shapes.add_textbox(Inches(left_in), Inches(top_in),
                                   Inches(width_in), Inches(0.8))
    tf = box.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.text = text
    p.font.name = FONT
    p.font.size = SIZE_FOOTNOTE
    p.font.color.rgb = GREY_TEXT
    p.alignment = PP_ALIGN.LEFT
    return box
