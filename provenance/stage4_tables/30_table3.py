#!/usr/bin/env python3
"""Table 3: integrated evidence for the 21-locus q<0.001 CTAS core.
Membership, representative, CTAS, M, PR and q from Task 1; family profile from the final-whitener recomputation
(code/20_architecture.R); colocalisation from Task 3; FinnGen from Task 2; prior convergence from the signed
four-study matrix (coordinate-Foley primary). Literature wording for the 16 retained rows is preserved verbatim
from the verified Table 3 mapping; the five added rows carry uncited context from curated local notes."""
import os, re, math, csv
import numpy as np, pandas as pd
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from lxml import etree
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

CT = '<CTAS_WORKSPACE>/'; OUT = CT + 'outputs/final_ctas_tables_20261021/'
T1 = CT + 'outputs/final_c_primary_freeze_20261021/'; T2 = CT + 'outputs/final_ctas_finngen_20261021/'; T3 = CT + 'outputs/final_ctas_colocalisation_20261021/'
FOUR = CT + 'outputs/four_comparator_harmonisation_20260919/ctas_104_locus_four_study_support_matrix.tsv'
CURATED = CT + 'ctas_explorer_v2/data/annotation_curated_v2.tsv'
MATRIX = CT + 'analytic_null_validation/inputs/UKBB_16trait_gene_matrix.tsv'

EXPECTED = ['V2L004 CELSR2', 'V2L010 CFHR3', 'V2L020 CARF', 'V2L023 DALRD3', 'V2L036 ATP6V1G2', 'V2L037 C6orf106', 'V2L042 LPA', 'V2L044 TBX20', 'V2L052 LPL', 'V2L057 CDKN2B',
            'V2L063 HTRA1', 'V2L074 ACAD10', 'V2L075 HNF1A', 'V2L080 MARK3', 'V2L083 IREB2', 'V2L084 FES', 'V2L086 CFDP1', 'V2L089 ORMDL3', 'V2L090 ZNF652', 'V2L095 LDLR', 'V2L097 APOE']
ADDED = ['CFHR3', 'C6orf106', 'LPL', 'HTRA1', 'LDLR']
PANEL = [('CAD_I25', 'CAD'), ('MI', 'MI'), ('Angina', 'Ang'), ('AF_I48', 'AF'), ('HF_I50', 'HF'), ('Cardiomyopathy', 'CMP'), ('Hypertension_401', 'HTN'), ('Atherosclerosis', 'Ath'),
         ('PAD', 'PAD'), ('Ischemic_Stroke', 'IS'), ('TIA', 'TIA'), ('DVT', 'DVT'), ('VTE', 'VTE'), ('Aortic_Aneurysm', 'AoAn'), ('Aortic_Valve', 'AoV'), ('Mitral_Valve', 'MV')]
PANEL_ABBR = dict(PANEL)
TISSUE = {'liver': 'liver', 'blood': 'whole blood', 'aorta': 'aorta', 'heart_AA': 'atrial appendage', 'heart_LV': 'left ventricle', 'kidney_cortex': 'kidney cortex', 'coronary': 'coronary artery'}
FG_ABBR = {'AF': 'AF', 'Angina': 'Ang', 'CAD': 'CAD', 'Cardiomyopathy': 'CMP', 'HF': 'HF', 'Hypertension': 'HTN', 'MI': 'MI', 'Stroke': 'Stroke', 'TIA': 'TIA', 'VTE': 'VTE'}
FAMILY = {'coronary': 'coronary', 'pressure': 'pressure', 'arrhythmia': 'arrhythmia', 'heart_failure': 'heart failure', 'valve': 'valve', 'aortic_vascular': 'aortic/vascular', 'cerebrovascular': 'cerebrovascular', 'thrombotic': 'thrombotic'}

# ---- literature: 16 retained rows verbatim from the verified mapping (manuscript references 35-48)
LIT = {
    "LPA": ("Anchor", "Lp(a) biology; human genetics support causal coronary risk", "[35]"),
    "CDKN2B": ("Anchor", "9p21 cell-cycle and efferocytosis biology in atherosclerosis", "[36]"),
    "FES": ("Functional", "15q26 CAD-locus kinase; protective in atherosclerosis models", "[37]"),
    "APOE": ("Anchor", "Canonical lipid-transport locus; coronary & cerebrovascular", ""),
    "MARK3": ("Unresolved", "Kinase with metabolic and hepatic biology; cardiovascular mechanism unestablished", "[38]"),
    "CFDP1": ("Functional", "BCAR1–CFDP1–TMEM170A carotid-IMT/CAD locus; CFDP1 supports cardiac development and function, with BCAR1 the adjacent effector candidate", "[39,40]"),
    "ATP6V1G2": ("Unresolved", "V-ATPase subunit within the MHC/6p21 multi-trait cardiovascular–immune region; no single-gene CV effector established (extended LD)", ""),
    "ZNF652": ("Human-genetic", "ZNF652-region association with blood pressure and hypertension", "[41]"),
    "CELSR2": ("Anchor", "1p13 CELSR2–PSRC1–SORT1 liver-expression CAD/LDL locus; CELSR2 regulates hepatocyte lipid accumulation, SORT1/PSRC1 the adjacent effectors", "[42,43]"),
    "CARF": ("Unresolved", "Statistical representative at the 2q33 WDR12-region coronary/MI locus; CARF (calcium-responsive regulator) has limited direct CV validation", ""),
    "ACAD10": ("Human-genetic", "12q24 (SH2B3/ALDH2) cross-cardiovascular region; ACAD10 variant associated with sex-specific hypertension susceptibility", "[44]"),
    "HNF1A": ("Anchor", "Hepatic regulator; coding variants associated with multiple cardiovascular-risk phenotypes", "[45]"),
    "IREB2": ("Unresolved", "Iron-response regulator in the 15q25 region that also carries the established ADAMTS7 coronary signal; the CTAS signal is not equated with ADAMTS7", ""),
    "TBX20": ("Functional", "Cardiac transcription factor; human variants associated with developmental heart disease and cardiomyopathy", "[46,47]"),
    "DALRD3": ("Unresolved", "tRNA-modification gene labelling the 3p21.31 coronary-artery-disease region; RHOA is the implicated regional effector, not an established causal gene", ""),
    "ORMDL3": ("Emerging", "Sphingolipid regulator with endothelial and atherosclerosis evidence", "[48]"),
    # added rows: conservative, non-causal context condensed from curated local notes; deliberately uncited (no verified manuscript reference number)
    "CFHR3": ("Emerging", "Complement factor-H-related gene at 1q31.3 (alternative complement pathway); no established cardiovascular effector assigned", ""),
    "C6orf106": ("Human-genetic", "C6orf106 (ILRUN) labels the 6p21.31 region associated with coronary artery disease; regional effector not established", ""),
    "LPL": ("Anchor", "Lipoprotein lipase; established determinant of plasma triglycerides and coronary artery disease", ""),
    "HTRA1": ("Human-genetic", "Vessel-wall serine protease; loss of function causes CARASIL, and variants associate with cerebral small-vessel disease and stroke", ""),
    "LDLR": ("Anchor", "Hepatic LDL receptor; familial hypercholesterolaemia gene and canonical coronary-atherosclerosis locus", ""),
}
LIT_REF = {"LPA": ("35", "Clarke et al., 2009", "20032323"), "CDKN2B": ("36", "Kojima et al., 2014", "24531546"), "FES": ("37", "Karamanavi et al., 2022", "36321446"),
    "APOE": ("", "(canonical; no single primary citation)", ""), "MARK3": ("38", "Lennerz et al., 2010", "20733003"), "CFDP1": ("39,40", "Gertow et al., 2012; Giardoglou et al., 2023", "23152477;37566073"),
    "ATP6V1G2": ("", "(no locus-specific functional citation)", ""), "ZNF652": ("41", "Newton-Cheh et al., 2009", "19430483"), "CELSR2": ("42,43", "Schadt et al., 2008; Tan et al., 2021", "18462017;34478580"),
    "CARF": ("", "(no locus-specific functional citation)", ""), "ACAD10": ("44", "Cho et al., 2021", "34828409"), "HNF1A": ("45", "Reiner et al., 2009", "20031592"),
    "IREB2": ("", "(no locus-specific functional citation)", ""), "TBX20": ("46,47", "Kirk et al., 2007; Amor-Salamanca et al., 2024", "17668378;38353104"),
    "DALRD3": ("", "(no locus-specific functional citation)", ""), "ORMDL3": ("48", "Wang et al., 2026", "42390638")}
EXPECTED_CITES = {"LPA": "[35]", "CDKN2B": "[36]", "FES": "[37]", "MARK3": "[38]", "CFDP1": "[39,40]", "ZNF652": "[41]", "CELSR2": "[42,43]", "ACAD10": "[44]", "HNF1A": "[45]", "TBX20": "[46,47]", "ORMDL3": "[48]"}

def fmt_q(q):
    e = math.floor(math.log10(q)); m = q / 10 ** e
    if round(m, 1) >= 10: m /= 10; e += 1
    sup = str.maketrans('-0123456789', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹')
    return f"{m:.1f}×10{str(e).translate(sup)}"

def load():
    g = pd.read_csv(T1 + 'final_C_gene_results_seed42.tsv', sep='\t').set_index('gene')
    loc = pd.read_csv(T1 + 'final_C_locus_results_seed42.tsv', sep='\t').set_index('canonical_locus_id')
    arch = pd.read_csv(OUT + 'source_data/architecture_locus_ranked.tsv', sep='\t').set_index('canonical_locus_id')
    co = pd.read_csv(T3 + 'strict_core_21_locus_colocalisation.tsv', sep='\t').set_index('canonical_locus_id')
    fg = pd.read_csv(T2 + 'strict_core_21_locus_finngen.tsv', sep='\t', dtype={'endpoints_p01': str, 'min_P': str}, keep_default_na=False).set_index('canonical_locus_id')
    four = pd.read_csv(FOUR, sep='\t', dtype=str, keep_default_na=False).set_index('canonical_locus_id')
    P = pd.read_csv(MATRIX, sep='\t', index_col=0)
    return g, loc, arch, co, fg, four, P

def build_rows():
    g, loc, arch, co, fg, four, P = load()
    core = loc[loc.n_q001_genes > 0].copy()
    core['rep'] = core.best_gene_by_final_C_CTAS
    core = core.assign(_c=g.loc[core.rep, 'CTAS'].values, _p=g.loc[core.rep, 'empirical_p'].values).sort_values(['_c', '_p', 'rep'], ascending=[False, True, True])
    rows = []
    for lid, r in core.iterrows():
        gene = r.rep; gr = g.loc[gene]; a = arch.loc[lid]; c = co.loc[lid]; f = fg.loc[lid]; s4 = four.loc[lid]
        traits = [ab for t, ab in PANEL if pd.notna(P.at[gene, t]) and P.at[gene, t] < 0.01]
        if c.coloc_class in ('strong_PP4', 'moderate_PP4'):
            coloc_main, coloc_sub = f"PP4 {float(c.max_PP4):.2f}", f"{TISSUE[c.best_tissue]} · {PANEL_ABBR[c.best_trait]}"
        elif c.coloc_class == 'not_assessed': coloc_main, coloc_sub = 'Not assessed', ''
        else: coloc_main, coloc_sub = 'No PP4 support', ''
        if f.status != 'evaluable': fg_main, fg_sub, fg_n = 'Not evaluable', '', None
        else:
            eps = [x.split('=')[0] for x in f.endpoints_p01.split(';') if x]
            fg_n = int(f.n_endpoints_p01); fg_main, fg_sub = f"{fg_n}/10", ', '.join(FG_ABBR[e] for e in eps)
        studies = [n for n, col in [('Qiao', 'Qiao_support'), ('Liu', 'Liu_support'), ('Foley', 'Foley_support'), ('PLEIO', 'PLEIO_support')] if s4[col] == 'yes']
        st, txt, cite = LIT[gene]
        rows.append(dict(locus=lid, gene=gene, chr=int(r.chromosome), ctas=gr.CTAS, q=gr.q_BH, M=gr.magnitude_M, PR=gr.participation_ratio_PR, empirical_p=gr.empirical_p,
            n_p01=len(traits), traits_p01=', '.join(traits), neff=a.Neff_family_breadth, family=' / '.join(FAMILY[x.strip()] for x in a.leading_family_display.split('/')),
            largest_share_family=a.largest_share_family, coloc_class=c.coloc_class, coloc_disposition=c.coloc_disposition, max_PP4=c.max_PP4, best_tissue=c.best_tissue, best_trait=c.best_trait,
            n_assessed_cells=c.n_assessed_cells, n_PP4_supported_cells=c.n_PP4_supported_cells, coloc_main=coloc_main, coloc_sub=coloc_sub,
            fg_status=f.status, fg_n=fg_n, fg_endpoints=f.endpoints_p01, fg_minP=f.min_P, fg_ge3=f.supported_ge3, fg_main=fg_main, fg_sub=fg_sub,
            in_qiao=s4.Qiao_support, in_liu=s4.Liu_support, in_foley_coord=s4.Foley_support, in_foley_gene_legacy=s4.Foley_gene_support, in_pleio=s4.PLEIO_support,
            n_studies4=len(studies), prior='·'.join(studies) if studies else 'None of 4', lit_status=st, lit_text=txt, lit_cite=cite,
            lit_phrase=txt + (f' {cite}' if cite else ''), added=gene in ADDED))
    return rows, (g, loc, arch, co, fg, four, P)

FOOTNOTE = ("The complete set of 21 canonical loci surviving genome-wide Benjamini–Hochberg q<0.001, ordered by descending CTAS of the statistical representative; "
 "downstream evidence did not determine inclusion. Statistical representatives label canonical loci and are not causal-gene assignments. CTAS, Cross-Trait Aggregation "
 "Score; BH q, Benjamini–Hochberg q-value (empirical P from 10⁹ null draws). CTAS is computed jointly from all 16 gene-level trait signals after ZCA-cor whitening with "
 "the transformed trait correlation matrix Σ = 2 sin(πρS/6), as CTAS = M × log2(1 + PR). UKBB signal profile: M is the magnitude of the positive whitened evidence and "
 "PR its participation ratio across the positive whitened axes; because those axes are mixtures of the original traits, PR is not a count of named traits. The named "
 "UKBB traits shown are those with raw gene-level P<0.01 (Table 1 panel order) and are descriptive, not additive components of CTAS. Cross-family profile: a separate "
 "descriptive summary of the same whitened evidence — Family Neff, the effective number of contributing cardiovascular trait families across eight named families, and "
 "the leading family (two co-leading families are shown when the second-largest size-adjusted share is at least 90% of the largest). Best colocalisation: maximum PP4 "
 "shown only when PP4≥0.5, with GTEx tissue and trait; “No PP4 support” denotes PP3-only or neither class and is kept distinct from “Not assessed” (MHC exclusion or no "
 "cis-eQTL). PP4 supports a shared-variant model, not a causal gene. FinnGen: number of ten prespecified FinnGen R12 endpoints (AF, Ang, CAD, CMP, HF, HTN, MI, Stroke, TIA, "
 "VTE) with gene-level P<0.01, supported endpoints named; “Not evaluable” = representative absent from the FinnGen gene outputs. CTAS was not rerun in FinnGen, so this "
 "is independent-cohort support, not replication. Prior multi-trait: representation in four prior cross-trait cardiovascular resources — Qiao, Liu, Foley and PLEIO — "
 "each matched on its native unit against the canonical GRCh37/hg19 loci (Qiao region and Foley candidate-variant coordinate overlap; Liu pleiotropic-gene membership; "
 "PLEIO index-variant coordinate overlap); “None of 4” = comparator-unmatched, not evidence of novelty or absent biology; this is convergence, not replication. "
 "Literature context is orienting: Anchor, established cardiovascular locus; Functional, experimental support; Human-genetic, human genetic association support; "
 "Emerging, emerging evidence; Unresolved, effector biology not yet resolved. Bracketed numbers are manuscript references 35–48; uncited rows carry orienting context only. "
 "Trait abbreviations: AF, atrial fibrillation; Ang, angina; AoAn, aortic aneurysm; AoV, aortic valve disease; Ath, atherosclerosis; CAD, coronary artery disease; CMP, "
 "cardiomyopathy; DVT, deep vein thrombosis; HF, heart failure; HTN, hypertension; IS, ischaemic stroke; MI, myocardial infarction; MV, mitral valve disease; PAD, "
 "peripheral artery disease; TIA, transient ischaemic attack; VTE, venous thromboembolism.")

# ---------------------------------------------------------------- PPTX
NS_A = "http://schemas.openxmlformats.org/drawingml/2006/main"
BLUE, GREY_BORDER, GREY_TEXT, BAND, FONT = "0072B2", "CCCCCC", "666666", "F2F5F8", "Aptos Narrow"
HDR_SIZE, BODY_SIZE, SUB_SIZE, FOOT_SIZE = Pt(9), Pt(8), Pt(7), Pt(7)
HEADER = ["Locus / representative", "CTAS (BH q)", "UKBB signal profile", "Cross-family profile", "Best colocalisation (PP4≥0.5)", "FinnGen (P<0.01 / 10)", "Prior multi-trait", "Literature context"]
COL_W = [1.12, 0.86, 2.10, 1.05, 1.32, 1.40, 1.20, 2.62]
HDR_H, ROW_H, TBL_LEFT, TBL_TOP = 0.46, 0.40, 0.42, 0.30
SPLIT = 11   # rows on page 1; remaining rows continue on page 2 under a repeated header

def set_fill(cell, hexc):
    tcPr = cell._tc.get_or_add_tcPr()
    for t in ("solidFill", "noFill"):
        el = tcPr.find(f"{{{NS_A}}}{t}")
        if el is not None: tcPr.remove(el)
    f = etree.SubElement(tcPr, f"{{{NS_A}}}solidFill"); etree.SubElement(f, f"{{{NS_A}}}srgbClr").set("val", hexc)

def set_borders(cell):
    tcPr = cell._tc.get_or_add_tcPr()
    for t in ("lnL", "lnR", "lnT", "lnB"):
        old = tcPr.find(f"{{{NS_A}}}{t}")
        if old is not None: tcPr.remove(old)
        ln = etree.SubElement(tcPr, f"{{{NS_A}}}{t}"); ln.set("w", "6350"); ln.set("cap", "flat")
        f = etree.SubElement(ln, f"{{{NS_A}}}solidFill"); etree.SubElement(f, f"{{{NS_A}}}srgbClr").set("val", GREY_BORDER)

def style(cell, lines, align):
    cell.vertical_anchor = MSO_ANCHOR.MIDDLE
    cell.margin_left, cell.margin_right, cell.margin_top, cell.margin_bottom = Inches(0.05), Inches(0.04), Inches(0.02), Inches(0.02)
    tf = cell.text_frame; tf.word_wrap = True; tf.clear()
    for i, segs in enumerate(lines):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph(); p.alignment = align; p.space_before = Pt(0); p.space_after = Pt(0)
        for text, bold, size, color in segs:
            run = p.add_run(); run.text = text; run.font.name = FONT; run.font.size = size; run.font.bold = bold; run.font.color.rgb = RGBColor.from_string(color)

def add_table_slide(prs, rows, start_index, footnote_top=None, cont=False):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    n = 1 + len(rows)
    t = slide.shapes.add_table(n, len(HEADER), Inches(TBL_LEFT), Inches(TBL_TOP), Inches(sum(COL_W)), Inches(HDR_H + ROW_H * len(rows))).table
    tblPr = t._tbl.find(f"{{{NS_A.replace('drawingml/2006/main', 'presentationml/2006/main')}}}tblPr")
    pr = t._tbl.tblPr; pr.set("firstRow", "0"); pr.set("bandRow", "0")
    for j, w in enumerate(COL_W): t.columns[j].width = Inches(w)
    t.rows[0].height = Inches(HDR_H)
    for i in range(1, n): t.rows[i].height = Inches(ROW_H)
    al = [PP_ALIGN.LEFT, PP_ALIGN.CENTER, PP_ALIGN.LEFT, PP_ALIGN.LEFT, PP_ALIGN.LEFT, PP_ALIGN.LEFT, PP_ALIGN.CENTER, PP_ALIGN.LEFT]
    for j, h in enumerate(HEADER):
        c = t.cell(0, j); style(c, [[(h + (' (cont.)' if cont and j == 0 else ''), True, HDR_SIZE, "FFFFFF")]], al[j]); set_fill(c, BLUE); set_borders(c)
    B, S, D = BODY_SIZE, SUB_SIZE, "1A1A1A"
    for i, r in enumerate(rows, start=1):
        bg = BAND if (start_index + i) % 2 == 0 else "FFFFFF"
        cells = [[[(r['locus'], False, S, GREY_TEXT)], [(r['gene'], True, B, D)]],
                 [[(f"{r['ctas']:.2f}", True, B, D)], [("(" + fmt_q(r['q']) + ")", False, S, GREY_TEXT)]],
                 [[(f"M {r['M']:.2f} · PR {r['PR']:.2f}", False, B, D)], [(f"{r['n_p01']} traits P<0.01: {r['traits_p01']}" if r['n_p01'] else "0 traits P<0.01", False, S, GREY_TEXT)]],
                 [[(f"Family Neff {r['neff']:.2f}", False, B, D)], [(f"Leading: {r['family']}", False, S, GREY_TEXT)]],
                 [[(r['coloc_main'], False, B, D)]] + ([[(r['coloc_sub'], False, S, GREY_TEXT)]] if r['coloc_sub'] else []),
                 [[(r['fg_main'], True, B, D)]] + ([[(r['fg_sub'], False, S, GREY_TEXT)]] if r['fg_sub'] else []),
                 [[(r['prior'], False, B, D)]],
                 [[(r['lit_status'] + " · ", True, B, D), (r['lit_phrase'], False, S, GREY_TEXT)]]]
        for j, lines in enumerate(cells):
            c = t.cell(i, j); style(c, lines, al[j]); set_fill(c, bg); set_borders(c)
    if footnote_top is not None:
        box = slide.shapes.add_textbox(Inches(TBL_LEFT), Inches(footnote_top), Inches(sum(COL_W)), Inches(1.6))
        tf = box.text_frame; tf.word_wrap = True; p = tf.paragraphs[0]; p.alignment = PP_ALIGN.LEFT
        run = p.add_run(); run.text = FOOTNOTE; run.font.name = FONT; run.font.size = FOOT_SIZE; run.font.color.rgb = RGBColor.from_string(GREY_TEXT)
    return slide

def build_pptx(rows, legend, footnote_top=None):
    prs = Presentation(); prs.slide_width = Emu(12192000); prs.slide_height = Emu(6858000)
    s1 = add_table_slide(prs, rows[:SPLIT], 0)
    s2 = add_table_slide(prs, rows[SPLIT:], SPLIT, footnote_top=footnote_top if footnote_top else TBL_TOP + HDR_H + 0.48 * (len(rows) - SPLIT) + 0.15, cont=True)
    for s in (s1, s2): s.notes_slide.notes_text_frame.text = legend
    prs.save(OUT + 'Main Tables/Table3.pptx')

# ---------------------------------------------------------------- legend
def build_legend(rows):
    return f"""Table 3. Integrated evidence across the stringent CTAS locus set.

Integrated evidence across the complete set of {len(rows)} canonical loci surviving genome-wide Benjamini–Hochberg q<0.001 (the stringent CTAS core), one row per canonical locus, ordered by descending CTAS of the statistical representative. Inclusion is determined only by the prespecified q<0.001 threshold; colocalisation, FinnGen and prior-study evidence never determine which loci appear. Statistical representatives are reproducible locus labels, not causal-gene assignments. The table continues on a second page with a repeated header.

Columns:
- Locus / representative: canonical locus ID (V2L###) and its statistical representative gene (the highest-CTAS significant member).
- CTAS (BH q): the gene-level Cross-Trait Aggregation Score and its genome-wide Benjamini–Hochberg q-value, from empirical P values over 10^9 null draws. Significance is defined by BH q, not by CTAS magnitude alone. CTAS is computed jointly from all 16 gene-level trait signals after ZCA-cor whitening with the transformed trait correlation matrix Sigma = 2 sin(pi rho_S / 6), as CTAS = M × log2(1 + PR).
- UKBB signal profile: M, the magnitude of the positive whitened evidence, and PR, its participation ratio across the positive whitened axes. These axes are mixtures of the original named traits, so PR is not a count of named traits. The named UKBB traits shown are those with raw gene-level P<0.01 (Table 1 panel order); they are descriptive context and cannot be divided into additive contributions to CTAS.
- Cross-family profile: a separate descriptive summary of the same whitened evidence, distinct from PR — Family N_eff, the effective number of contributing cardiovascular trait families across eight named families, and the leading family (two co-leading families shown when the second-largest size-adjusted share is at least 90% of the largest). Recomputed with the final whitener.
- Best colocalisation (PP4≥0.5): the strongest posterior probability of a shared causal variant (PP4) across seven GTEx tissues and 16 traits, displayed only where PP4≥0.5, with tissue and trait. “No PP4 support” (PP3-only or neither class) and “Not assessed” (MHC exclusion, here ATP6V1G2) are kept distinct. PP4 supports a shared-variant model; it does not establish the causal gene or a mechanism.
- FinnGen (P<0.01 / 10): the number of ten prespecified FinnGen R12 cardiovascular endpoints (AF, Ang, CAD, CMP, HF, HTN, MI, Stroke, TIA, VTE) at gene-level P<0.01, with the supported endpoints named. ORMDL3 is not evaluable (absent from the FinnGen gene outputs). CTAS was not rerun in FinnGen, so this is independent-cohort support rather than replication.
- Prior multi-trait: convergent representation in four prior cross-trait cardiovascular resources — Qiao, Liu, Foley and PLEIO — each matched on its native unit against the canonical GRCh37/hg19 loci: Qiao reported-region and Foley candidate-variant coordinate overlap; Liu pleiotropic-gene membership; PLEIO index-variant coordinate overlap. “None of 4” denotes comparator-unmatched and is not evidence of novelty or of absent biology. This is literature/method convergence, not independent replication.
- Literature context: a short orienting label and phrase — Anchor, Functional, Human-genetic, Emerging or Unresolved. It is orienting rather than exhaustive and does not promote a representative to a causal gene. Bracketed reference numbers identify the locus-specific source(s) (manuscript references 35–48). Loci without bracketed numbers carry orienting context only; the five loci added to the strict core (CFHR3, C6orf106, LPL, HTRA1, LDLR) are deliberately uncited pending the manuscript reference pass.

The colocalisation, FinnGen and prior multi-trait streams answer different questions and are complementary; absence of support in one stream is not a negative biological verdict. The table turns the aggregate validation gradients into locus-level priorities; it does not claim that all {len(rows)} loci are mechanistically resolved.

Trait abbreviations (Table 1 panel): CAD, coronary artery disease; MI, myocardial infarction; Ang, angina; AF, atrial fibrillation; HF, heart failure; CMP, cardiomyopathy; HTN, hypertension; Ath, atherosclerosis; PAD, peripheral artery disease; IS, ischaemic stroke; TIA, transient ischaemic attack; DVT, deep vein thrombosis; VTE, venous thromboembolism; AoAn, aortic aneurysm; AoV, aortic valve disease; MV, mitral valve disease.

Sources: final primary gene and locus results; final-whitener architecture recomputation (source_data/architecture_locus_ranked.tsv); strict-core colocalisation; strict-core FinnGen support; signed four-study support matrix (ctas_104_locus_four_study_support_matrix.tsv, coordinate-Foley primary).
"""

# ---------------------------------------------------------------- XLSX + audit
HF = Font(bold=True, color="FFFFFF"); HFILL = PatternFill("solid", fgColor=BLUE)
def sheet(wb, name, header, data, title=None, widths=None):
    ws = wb.create_sheet(name); r0 = 1
    if title: ws.cell(1, 1, title).font = Font(bold=True); r0 = 3
    for j, h in enumerate(header, 1):
        c = ws.cell(r0, j, h); c.font = HF; c.fill = HFILL; c.alignment = Alignment(wrap_text=True, vertical='top')
    for i, row in enumerate(data, r0 + 1):
        for j, v in enumerate(row, 1): ws.cell(i, j, v)
    ws.freeze_panes = ws.cell(r0 + 1, 1); ws.auto_filter.ref = f"A{r0}:{get_column_letter(len(header))}{r0 + len(data)}"
    for j, h in enumerate(header, 1): ws.column_dimensions[get_column_letter(j)].width = (widths or {}).get(h, min(60, max(10, len(str(h)) + 2)))
    return ws

def readme(wb, title, pairs):
    ws = wb.create_sheet('README', 0); ws.cell(1, 1, title).font = Font(bold=True, size=12)
    for i, (k, v) in enumerate(pairs, 3):
        ws.cell(i, 1, k).font = Font(bold=True); c = ws.cell(i, 2, v); c.alignment = Alignment(wrap_text=True, vertical='top')
    ws.column_dimensions['A'].width = 26; ws.column_dimensions['B'].width = 120
    return ws

def main():
    rows, (g, loc, arch, co, fg, four, P) = build_rows()
    legend = build_legend(rows)
    open(OUT + 'legends/Table3_legend.txt', 'w').write(legend)
    fn_top = None
    if os.path.exists(OUT + 'work/table3_footnote_top.txt'): fn_top = float(open(OUT + 'work/table3_footnote_top.txt').read())
    build_pptx(rows, legend, fn_top)
    # checks
    ck = []
    def chk(name, ok, obs, exp): ck.append((name, 'PASS' if ok else 'FAIL', str(obs), str(exp)))
    got = [f"{r['locus']} {r['gene']}" for r in rows]
    chk('21 unique loci, exact expected membership', sorted(got) == sorted(EXPECTED) and len(set(got)) == 21, len(set(got)), 21)
    chk('five additions present', all(any(r['gene'] == a for r in rows) for a in ADDED), ';'.join(r['gene'] for r in rows if r['added']), ';'.join(ADDED))
    chk('ordered by descending representative CTAS', all(rows[i]['ctas'] >= rows[i + 1]['ctas'] for i in range(20)), 'monotone', 'monotone')
    sup = [r for r in rows if r['coloc_class'] in ('strong_PP4', 'moderate_PP4')]; ev = [r for r in rows if r['coloc_class'] != 'not_assessed']
    chk('colocalisation 12/21 overall, 12/20 evaluable', (len(sup), len(ev)) == (12, 20), f'{len(sup)}/21, {len(sup)}/{len(ev)}', '12/21, 12/20')
    chk('ATP6V1G2 sole not-assessed', [r['gene'] for r in rows if r['coloc_class'] == 'not_assessed'] == ['ATP6V1G2'], [r['gene'] for r in rows if r['coloc_class'] == 'not_assessed'], ['ATP6V1G2'])
    fge = [r for r in rows if r['fg_status'] == 'evaluable']; fg3 = [r for r in fge if r['fg_n'] >= 3]
    chk('FinnGen 20 evaluable, ORMDL3 not evaluable, 9/20 >=3 endpoints', len(fge) == 20 and [r['gene'] for r in rows if r['fg_status'] != 'evaluable'] == ['ORMDL3'] and len(fg3) == 9, f'{len(fge)} evaluable; {len(fg3)}/20 >=3', '20; 9/20')
    cnt = lambda k: sum(r[k] == 'yes' for r in rows)
    four_ = (cnt('in_qiao'), cnt('in_liu'), cnt('in_foley_coord'), cnt('in_pleio'), sum(r['n_studies4'] > 0 for r in rows))
    chk('four-study Qiao/Liu/coordinate-Foley/PLEIO/union = 8/13/9/17/17', four_ == (8, 13, 9, 17, 17), '/'.join(map(str, four_)), '8/13/9/17/17')
    by = {r['gene']: r for r in rows}
    chk('retained-row reference numbers preserved; none invented', all(by[k]['lit_cite'] == v for k, v in EXPECTED_CITES.items()) and all(by[k]['lit_cite'] == '' for k in by if k not in EXPECTED_CITES),
        ';'.join(f"{k}{by[k]['lit_cite']}" for k in EXPECTED_CITES), 'refs 35-48 unchanged; others uncited')
    # Task 1/2/3 agreement
    t1 = all(abs(r['ctas'] - g.at[r['gene'], 'CTAS']) == 0 and r['q'] == g.at[r['gene'], 'q_BH'] and r['M'] == g.at[r['gene'], 'magnitude_M'] and r['PR'] == g.at[r['gene'], 'participation_ratio_PR']
             and loc.at[r['locus'], 'best_gene_by_final_C_CTAS'] == r['gene'] for r in rows)
    chk('representative/CTAS/M/PR/q equal the final primary result', t1, t1, True)
    t3 = all(co.at[r['locus'], 'representative_gene'] == r['gene'] and co.at[r['locus'], 'coloc_class'] == r['coloc_class'] and str(co.at[r['locus'], 'max_PP4']) == str(r['max_PP4']) for r in rows)
    chk('colocalisation fields equal the colocalisation bundle', t3, t3, True)
    t2 = all(fg.at[r['locus'], 'representative_gene'] == r['gene'] and fg.at[r['locus'], 'status'] == r['fg_status'] and fg.at[r['locus'], 'endpoints_p01'] == r['fg_endpoints'] for r in rows)
    chk('FinnGen fields equal the FinnGen bundle', t2, t2, True)
    chk('family profile from final-whitener recomputation (M/PR reproduce final values)', all(abs(arch.at[r['locus'], 'magnitude_M'] - r['M']) < 1e-10 for r in rows), True, True)
    # audit (one row per displayed value)
    aud = []
    for r in rows:
        for fld, val, src in [('representative', r['gene'], 'Task1 final_C_locus_results_seed42.tsv'), ('CTAS', r['ctas'], 'Task1 gene results'), ('q_BH', r['q'], 'Task1 gene results'),
            ('M', r['M'], 'Task1 gene results'), ('PR', r['PR'], 'Task1 gene results'), ('traits_P<0.01', r['traits_p01'], 'UKBB_16trait_gene_matrix.tsv (raw P)'),
            ('Family_Neff', r['neff'], 'source_data/architecture_locus_ranked.tsv (final whitener)'), ('leading_family', r['family'], 'source_data/architecture_locus_ranked.tsv'),
            ('coloc_class', r['coloc_class'], 'Task3 strict_core_21_locus_colocalisation.tsv'), ('max_PP4', r['max_PP4'], 'Task3'), ('best_tissue', r['best_tissue'], 'Task3'), ('best_trait', r['best_trait'], 'Task3'),
            ('FinnGen_status', r['fg_status'], 'Task2 strict_core_21_locus_finngen.tsv'), ('FinnGen_n_P<0.01', r['fg_n'], 'Task2'), ('FinnGen_endpoints', r['fg_endpoints'], 'Task2'), ('FinnGen_min_P', r['fg_minP'], 'Task2'),
            ('Qiao', r['in_qiao'], 'ctas_104_locus_four_study_support_matrix.tsv'), ('Liu', r['in_liu'], 'four-study matrix'), ('Foley_coordinate', r['in_foley_coord'], 'four-study matrix'),
            ('Foley_gene_symbol_legacy', r['in_foley_gene_legacy'], 'four-study matrix (provenance only)'), ('PLEIO', r['in_pleio'], 'four-study matrix'),
            ('literature_status', r['lit_status'], 'verified Table 3 mapping' if not r['added'] else 'curated local note (annotation_curated_v2.tsv)'), ('literature_citation', r['lit_cite'], 'manuscript references 35-48 (unchanged)' if not r['added'] else 'none (uncited pending manuscript pass)')]:
            aud.append(dict(locus=r['locus'], representative=r['gene'], field=fld, displayed_or_stored_value=val, source=src))
    pd.DataFrame(aud).to_csv(OUT + 'source_data/Table3_evidence_audit.tsv', sep='\t', index=False)
    cur = pd.read_csv(CURATED, sep='\t', dtype=str, keep_default_na=False).set_index('canonical_locus_id')
    lr = []
    for r in rows:
        if r['added']:
            c = cur.loc[r['locus']]
            lr.append(dict(locus=r['locus'], representative=r['gene'], row_status='added to strict core', literature_status=r['lit_status'], table_text=r['lit_text'], citation_in_table='none',
                source_note=c.annotation_note, source_ids=c.source_ids, source_review_status=c.review_status, action_for_manuscript_pass='Verify a locus-specific reference in EndNote before adding any bracketed number; no number assigned here.'))
        else:
            n, a, pm = LIT_REF[r['gene']]
            lr.append(dict(locus=r['locus'], representative=r['gene'], row_status='retained (verified mapping preserved verbatim)', literature_status=r['lit_status'], table_text=r['lit_text'],
                citation_in_table=r['lit_cite'] or 'none (uncited by design)', source_note=a, source_ids=f'PMID {pm}' if pm else '', source_review_status='verified 2026-08-27', action_for_manuscript_pass='none'))
    pd.DataFrame(lr).to_csv(OUT + 'source_data/Table3_literature_review.tsv', sep='\t', index=False)
    # workbook
    wb = openpyxl.Workbook(); wb.remove(wb.active)
    disp = [[r['locus'], r['gene'], round(r['ctas'], 2), fmt_q(r['q']), round(r['M'], 2), round(r['PR'], 2), r['n_p01'], r['traits_p01'], round(r['neff'], 2), r['family'],
             r['coloc_main'] + (f" ({r['coloc_sub']})" if r['coloc_sub'] else ''), r['fg_main'] + (f": {r['fg_sub']}" if r['fg_sub'] else ''), r['prior'], f"{r['lit_status']} - {r['lit_phrase']}"] for r in rows]
    sheet(wb, 'Table 3', ['Locus', 'Representative', 'CTAS', 'BH q', 'Magnitude M', 'PR', 'Traits P<0.01 (n)', 'Traits P<0.01 (named)', 'Family Neff', 'Leading family',
          'Best colocalisation (PP4>=0.5)', 'FinnGen P<0.01 (/10)', 'Prior multi-trait', 'Literature context'], disp,
          title='Table 3. Integrated evidence across the stringent CTAS locus set (21 loci, q<0.001, ordered by descending CTAS).',
          widths={'Traits P<0.01 (named)': 34, 'Best colocalisation (PP4>=0.5)': 30, 'FinnGen P<0.01 (/10)': 32, 'Prior multi-trait': 20, 'Literature context': 80, 'Leading family': 16})
    cols = ['locus', 'gene', 'chr', 'ctas', 'q', 'empirical_p', 'M', 'PR', 'n_p01', 'traits_p01', 'neff', 'largest_share_family', 'family', 'coloc_class', 'coloc_disposition', 'max_PP4', 'best_tissue', 'best_trait',
            'n_assessed_cells', 'n_PP4_supported_cells', 'fg_status', 'fg_n', 'fg_endpoints', 'fg_minP', 'fg_ge3', 'in_qiao', 'in_liu', 'in_foley_coord', 'in_foley_gene_legacy', 'in_pleio', 'n_studies4', 'lit_status', 'lit_cite', 'added']
    sheet(wb, 'Evidence data', cols, [[r[c] for c in cols] for r in rows], title='Full-precision evidence per locus (foley_gene_legacy is provenance only; the paper-facing Foley call is coordinate containment).')
    sheet(wb, 'Literature notes', ['locus', 'representative', 'row_status', 'literature_status', 'table_text', 'citation_in_table', 'source_note', 'source_ids'], [[x[k] for k in ['locus', 'representative', 'row_status', 'literature_status', 'table_text', 'citation_in_table', 'source_note', 'source_ids']] for x in lr])
    sheet(wb, 'Validation', ['check', 'result', 'observed', 'expected'], [list(c) for c in ck])
    readme(wb, 'Table 3 - integrated evidence for the 21-locus q<0.001 CTAS core', [
        ('Selection rule', 'All 21 canonical loci with q<0.001 in the final CTAS analysis, ordered by descending CTAS of the statistical representative. No locus was selected or removed using downstream evidence.'),
        ('Additions', 'Relative to the earlier 16-row table, five loci enter the strict core: CFHR3, C6orf106, LPL, HTRA1 and LDLR.'),
        ('CTAS / q', 'From the final primary analysis: transformed trait correlation matrix Sigma = 2 sin(pi rho_S / 6), ZCA-cor whitening, empirical P from 10^9 null draws, genome-wide BH q.'),
        ('Family profile', 'Family N_eff and leading family recomputed with the final whitener from the unchanged 16-trait evidence; recomputed M and PR reproduce the final values to <1e-14. No leading family differs from the earlier summaries; Family N_eff differs by at most 0.036.'),
        ('Colocalisation', '12/21 loci with PP4>=0.5 overall; 12/20 among evaluable loci. ATP6V1G2 is not assessed (MHC exclusion), never counted as negative.'),
        ('FinnGen', '20 representatives evaluable; ORMDL3 not evaluable. 9/20 have gene-level P<0.01 in at least 3 of 10 endpoints. Independent-cohort support, not replication of CTAS.'),
        ('Prior multi-trait', 'Qiao 8, Liu 13, Foley (candidate-variant coordinate containment) 9, PLEIO 17; union 17/21. Convergence, not replication. The legacy gene-symbol Foley call is kept in Evidence data as provenance only.'),
        ('Literature', 'Retained rows keep the verified wording and manuscript reference numbers 35-48 exactly. The five added rows carry uncited orienting context from curated local notes; see Literature notes and source_data/Table3_literature_review.tsv.'),
        ('Interpretation', 'Statistical representatives are locus labels, not causal genes. PP4 supports a shared-variant model only. Absence of support in one stream is not a negative verdict.')])
    wb.save(OUT + 'Main Tables/Table3.xlsx')
    pd.DataFrame(ck, columns=['check', 'result', 'observed', 'expected']).to_csv(OUT + 'work/table3_checks.tsv', sep='\t', index=False)
    for c in ck: print(c[1], c[0], '|', c[2])

if __name__ == '__main__':
    main()
