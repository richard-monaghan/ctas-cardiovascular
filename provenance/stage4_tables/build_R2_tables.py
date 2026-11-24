#!/usr/bin/env python3
"""Build the Results Section 2 supplementary workbooks: ST3 (per-gene calibration/
FDR), ST4-5 (simulations), ST15 (score-form robustness). Typed values, frozen
headers, filters, README + validation sheets. Reads locked/real source data only."""
import csv, os
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

CTAS = "<CTAS_WORKSPACE>"
PAPER = "<CTAS_V2_PAPER_WORKSPACE>"
ST = f"{PAPER}/Supplementary Tables"
R2 = f"{PAPER}/_source_assets/R2_review"
NOMINAL_1E6_Q = 19.43226
BLUE, SKY, LIGHT = "0072B2", "56B4E9", "DCEEF8"

def rd(path, delim=","):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter=delim))

def num(x):
    try:
        f = float(x)
        return int(f) if f == int(f) and abs(f) < 1e15 and "." not in str(x) and "e" not in str(x).lower() else f
    except (ValueError, TypeError):
        return x if x not in ("NA", "") else None

def style_header(ws, ncol, height=30):
    hf, hfill = Font(bold=True, color="FFFFFF", size=11), PatternFill("solid", fgColor=BLUE)
    for j in range(1, ncol + 1):
        c = ws.cell(1, j); c.font = hf; c.fill = hfill
        c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
    ws.row_dimensions[1].height = height
    ws.freeze_panes = "A2"
    ws.auto_filter.ref = f"A1:{get_column_letter(ncol)}1"

def write_sheet(ws, headers, rows, widths=None, numeric_cols=None):
    numeric_cols = numeric_cols or set()
    for j, h in enumerate(headers, 1):
        ws.cell(1, j, h)
        ws.column_dimensions[get_column_letter(j)].width = (widths or {}).get(h, max(11, len(h) + 2))
    for i, r in enumerate(rows, 2):
        for j, h in enumerate(headers, 1):
            v = r.get(h)
            if h in numeric_cols:
                v = num(v)
            ws.cell(i, j, v)
    style_header(ws, len(headers))

def readme_sheet(wb, lines, title="README"):
    ws = wb.create_sheet(title, 0)
    for i, t in enumerate(lines, 1):
        c = ws.cell(i, 1, t)
        if i == 1: c.font = Font(bold=True, size=13)
    ws.column_dimensions["A"].width = 118
    return ws

# ======================================================================= ST3
gss = rd(f"{CTAS}/panel_independence_stringency/gene_stringency_scores.csv")
g2l = {r["gene"]: r for r in rd(f"{CTAS}/panel_v2_locus_lock/v2_gene_to_locus.csv")}
master = {r["canonical_locus_id"]: r for r in rd(f"{CTAS}/panel_v2_locus_lock/v2_locus_master.csv")}
# documented per-gene tail-model sensitivity (null_tail_report.md, 5 deepest genes)
dsens = {"LPA": ("1.41e-09", "2.72"), "CDKN2B": ("2.63e-08", "1.49"), "IGF2R": ("5.35e-08", "1.25"),
         "FES": ("9.01e-08", "1.10"), "APOE": ("9.57e-08", "1.08")}

st3_rows = []
for r in gss:
    ctas = float(r["magnitude_score"]); q = float(r["q_BH"]); p = float(r["null_p"])
    dec = float(r["extrapolation_decades_below_band"])
    tier = ("core (q<0.001)" if r["survives_at_q001"]=="TRUE" else "q<0.01" if r["survives_at_q01"]=="TRUE"
            else "q<0.05" if r["survives_at_q05"]=="TRUE" else "ns")
    above = ctas > NOMINAL_1E6_Q
    loc = g2l.get(r["gene"])
    cid = loc["canonical_locus_id"] if loc else ""
    lab = master[cid]["display_locus_label"] if cid else ""
    pw, dl = dsens.get(r["gene"], ("", ""))
    st3_rows.append({
        "gene": r["gene"], "CTAS": ctas, "empirical_p_direct": ("" if above else p),
        "reported_p": p, "tail_extrapolated": ("yes" if dec > 0 else "no"),
        "above_nominal_1e6_quantile": ("yes" if above else "no"),
        "tail_model": ("exponential" if above else "direct"),
        "tail_model_powerlaw_p": pw, "tail_model_sensitivity_dlog10": dl,
        "q_BH": q, "q_tier": tier, "canonical_locus_id": cid, "statistical_locus_label": lab})
st3_rows.sort(key=lambda d: -d["CTAS"])
n_q05 = sum(r["q_tier"] in ("q<0.05","q<0.01","core (q<0.001)") for r in st3_rows)
n_q01 = sum(r["q_tier"] in ("q<0.01","core (q<0.001)") for r in st3_rows)
n_q001 = sum(r["q_tier"]=="core (q<0.001)" for r in st3_rows)
n_above = sum(r["above_nominal_1e6_quantile"]=="yes" for r in st3_rows)
n_extrap = sum(r["tail_extrapolated"]=="yes" for r in st3_rows)
assert (n_q05,n_q01,n_q001)==(250,102,32) and n_above==36 and n_extrap==25, (n_q05,n_q01,n_q001,n_above,n_extrap)

wb = Workbook(); wb.remove(wb.active)
readme_sheet(wb, [
    "ST2 — CTAS null calibration (Results Section 2).",
    "",
    "Source: panel_independence_stringency/gene_stringency_scores.csv (19,080 genes);",
    "canonical loci from panel_v2_locus_lock/ (V2L001-V2L104).",
    "CTAS = magnitude x breadth (Sigma-whitened, positive-part squared, participation-ratio breadth).",
    "reported_p = null_p (empirical vs 1,000,000 Gaussian-copula null draws; exponential tail model where extrapolated).",
    "empirical_p_direct = the direct empirical p, populated where the empirical null resolves the value directly.",
    "",
    "HEADLINE TAIL RESULT: Twenty-five genes required upper-tail extrapolation beyond the directly",
    "resolved empirical null (tail_extrapolated = yes; extrapolation_decades_below_band > 0).",
    "The 'tail_extrapolated_genes' sheet lists exactly these 25 genes.",
    "",
    "The above_nominal_1e6_quantile column (CTAS > 19.43) is a technical score-quantile flag, NOT a",
    "biological result and NOT an extrapolation claim: of the 36 genes above that nominal quantile, 25 required",
    "extrapolation and the other 11 still had directly resolved empirical-null p-values (they remain in the",
    "ordinary gene-level results, ST1, and need no separate claim). See the 'upper_tail_diagnostics' sheet",
    "for the full 36-gene technical audit with the requires_extrapolation field.",
    "",
    "tail_model = 'exponential' where extrapolated (R^2=0.994 vs power-law 0.974); 'direct' otherwise.",
    "tail_model_sensitivity_dlog10 = |log10(exp p) - log10(power-law p)| for the 5 genes with authoritative",
    "  per-model values (null_tail_report.md); LPA 2.72 is the largest. Tail-model choice did not change any q-tier.",
    "q_tier tiers: q<0.05 = 250 genes; q<0.01 = 102; q<0.001 = 32.",
], "README")
ws = wb.create_sheet("gene_results")
write_sheet(ws, ["gene","CTAS","empirical_p_direct","reported_p","tail_extrapolated",
    "above_nominal_1e6_quantile","tail_model","tail_model_powerlaw_p","tail_model_sensitivity_dlog10",
    "q_BH","q_tier","canonical_locus_id","statistical_locus_label"], st3_rows,
    widths={"gene":12,"statistical_locus_label":22,"tail_model":12,"above_nominal_1e6_quantile":15},
    numeric_cols={"CTAS","empirical_p_direct","reported_p","q_BH","tail_model_powerlaw_p","tail_model_sensitivity_dlog10"})
for i in range(2, len(st3_rows)+2):
    ws.cell(i,2).number_format="0.00"; ws.cell(i,3).number_format="0.00E+00"; ws.cell(i,4).number_format="0.00E+00"; ws.cell(i,10).number_format="0.00E+00"
# tail_extrapolated_genes = ONLY the 25 genuinely extrapolated genes
ws2 = wb.create_sheet("tail_extrapolated_genes")
tail_rows = [r for r in st3_rows if r["tail_extrapolated"]=="yes"]
assert len(tail_rows) == 25, len(tail_rows)
write_sheet(ws2, ["gene","CTAS","reported_p","tail_model","tail_model_powerlaw_p",
    "tail_model_sensitivity_dlog10","q_tier","canonical_locus_id"], tail_rows,
    numeric_cols={"CTAS","reported_p","tail_model_powerlaw_p","tail_model_sensitivity_dlog10"})
for i in range(2, len(tail_rows)+2):
    ws2.cell(i,2).number_format="0.00"; ws2.cell(i,3).number_format="0.00E+00"
# upper_tail_diagnostics = all 36 above the nominal 1e-6 score quantile, with an
# explicit requires_extrapolation flag (25 True / 11 False) for technical audit
wsd = wb.create_sheet("upper_tail_diagnostics")
diag_rows = [dict(r, requires_extrapolation=("yes" if r["tail_extrapolated"]=="yes" else "no"))
             for r in st3_rows if r["above_nominal_1e6_quantile"]=="yes"]
assert len(diag_rows) == 36 and sum(r["requires_extrapolation"]=="yes" for r in diag_rows) == 25
write_sheet(wsd, ["gene","CTAS","reported_p","requires_extrapolation","above_nominal_1e6_quantile",
    "tail_model","tail_model_powerlaw_p","tail_model_sensitivity_dlog10","q_tier","canonical_locus_id"], diag_rows,
    numeric_cols={"CTAS","reported_p","tail_model_powerlaw_p","tail_model_sensitivity_dlog10"})
for i in range(2, len(diag_rows)+2):
    wsd.cell(i,2).number_format="0.00"; wsd.cell(i,3).number_format="0.00E+00"
wsb = wb.create_sheet("q_boundary_summary")
write_sheet(wsb, ["q_tier","empirical_p_boundary","neglog10_p_boundary","n_genes","n_canonical_loci"], [
    {"q_tier":"q<0.05","empirical_p_boundary":0.000654999345000655,"neglog10_p_boundary":3.18375913,"n_genes":250,"n_canonical_loci":104},
    {"q_tier":"q<0.01","empirical_p_boundary":5.1999948000052e-05,"neglog10_p_boundary":4.28399709,"n_genes":102,"n_canonical_loci":47},
    {"q_tier":"q<0.001","empirical_p_boundary":9.99999000001e-07,"neglog10_p_boundary":6.00000043,"n_genes":32,"n_canonical_loci":16},
], numeric_cols={"empirical_p_boundary","neglog10_p_boundary","n_genes","n_canonical_loci"})
for i in range(2,5): wsb.cell(i,2).number_format="0.00E+00"
wsv = wb.create_sheet("validation_checks")
write_sheet(wsv, ["check","expected","observed","pass"], [
    {"check":"total genes","expected":19080,"observed":len(st3_rows),"pass":len(st3_rows)==19080},
    {"check":"q<0.05 genes","expected":250,"observed":n_q05,"pass":n_q05==250},
    {"check":"q<0.01 genes","expected":102,"observed":n_q01,"pass":n_q01==102},
    {"check":"q<0.001 genes","expected":32,"observed":n_q001,"pass":n_q001==32},
    {"check":"above_nominal_1e6_quantile","expected":36,"observed":n_above,"pass":n_above==36},
    {"check":"tail_p_extrapolated","expected":25,"observed":n_extrap,"pass":n_extrap==25},
    {"check":"row1 gene","expected":"LPA","observed":st3_rows[0]["gene"],"pass":st3_rows[0]["gene"]=="LPA"},
], numeric_cols={"expected","observed"})
wb.save(f"{ST}/ST2_CTAS_null_calibration.xlsx")
print("ST2 (null calibration) written:", n_q05, n_q01, n_q001, "above", n_above, "extrap", n_extrap)

# ==================================================================== ST4-5
sim1 = rd(f"{R2}/sim1_calibration_typeI.csv"); sim2 = rd(f"{R2}/sim2_fixed_total_signal.csv")
sim3 = rd(f"{R2}/sim3_breadth_by_strength.csv"); prov = rd(f"{R2}/sim_provenance.csv")
sim0 = rd(f"{R2}/sim0_null_reproduction.csv")
wb = Workbook(); wb.remove(wb.active)
readme_sheet(wb, [
    "ST3 — CTAS behaviour simulations (Results Section 2), run on the IMPLEMENTED CTAS statistic.",
    "",
    "CTAS depends on the 16-trait evidence only through the whitened vector z_white = Sigma^(-1/2) z.",
    "Under the Gaussian-copula null z ~ N(0,Sigma), z_white ~ N(0, I_16); controlled evidence is therefore",
    "specified directly in the decorrelated axis space (identity covariance). This is NOT a V1 recurrence sim.",
    "",
    "fixed_total_signal sheet: total non-centrality lambda = sum(mu^2) held constant while spread over k active",
    "  axes (mu = sqrt(lambda/k)). Shows CTAS rewards signal distributed across more decorrelated axes.",
    "breadth_by_strength sheet: grid of k active axes x per-axis strength mu; distinguishes weak-broad,",
    "  strong-narrow and strong-broad evidence.",
    "calibration sheet: nominal vs observed type-I on independent null genes (empirical p vs 1e6 reference null).",
    "selection-probability columns use per-gene NOMINAL null-survival cutoffs (a=0.05/0.01/0.001 -> CTAS",
    "  9.96/11.92/14.25); these are single-test nominal alphas, NOT the genome-wide BH-FDR gene cutoffs.",
    "All seeds, grids, replicate counts and checksums are on the 'parameters_provenance' sheet.",
], "README")
wsc = wb.create_sheet("calibration_typeI")
write_sheet(wsc, list(sim1[0].keys()), sim1, numeric_cols=set(sim1[0].keys())-{"within_MC"})
ws4 = wb.create_sheet("fixed_total_signal")
write_sheet(ws4, list(sim2[0].keys()), sim2, numeric_cols=set(sim2[0].keys()))
ws5 = wb.create_sheet("breadth_by_strength")
write_sheet(ws5, list(sim3[0].keys()), sim3, numeric_cols=set(sim3[0].keys()))
wsp = wb.create_sheet("parameters_provenance")
write_sheet(wsp, list(prov[0].keys()), prov, widths={"value":48})
ws0 = wb.create_sheet("validation_null_reproduction")
write_sheet(ws0, list(sim0[0].keys()), sim0, numeric_cols=set(sim0[0].keys()))
wb.save(f"{ST}/ST3_CTAS_behaviour_simulations.xlsx")
print("ST3 (behaviour simulations) written")

# ===================================================================== ST15
frt = rd(f"{CTAS}/panel_form_robustness/form_robustness_table.csv")
rc = {r["variant"]: r for r in rd(f"{CTAS}/panel_form_robustness/rank_concordance.csv")}
ts = {r["variant"]: r for r in rd(f"{CTAS}/panel_form_robustness/tail_sensitivity.csv")}
hv = {r["variant"]: r for r in rd(f"{CTAS}/panel_form_robustness/headline_verdicts.csv")}
clc = {r["variant"]: r for r in rd(f"{R2}/form_canonical_locus_counts.csv")}
st15 = []
for r in frt:
    v = r["variant"]; c = clc.get(v, {}); h = hv.get(v, {})
    st15.append({
        "form": v, "class": r["class"], "prespecified_sensible": ("yes" if r["class"]=="sensible" else "diagnostic_only"),
        "spearman_vs_reference": rc[v]["spearman_vs_ref"], "jaccard_q05_vs_ref": rc[v]["jaccard_sig_q05"],
        "jaccard_core_q001_vs_ref": rc[v]["jaccard_core_q001"],
        "n_sig_q05": r["n_sig_q05"], "n_sig_q01": r["n_sig_q01"], "n_sig_q001": r["n_sig_q001"],
        "canonical_loci_q05": c.get("q05_canonical_loci",""), "canonical_loci_q01": c.get("q01_canonical_loci",""),
        "canonical_loci_q001": c.get("q001_canonical_loci",""),
        "canonical_map_coverage_q05": (f"{c.get('q05_mappable','')}/{c.get('q05_genes','')}" if c else ""),
        "core_n": r["core_n"], "core_crossfamily_prop": r["core_crossfamily_prop"],
        "noncoronary_branch_OR": h.get("H3_odds_ratio",""), "H1_calibration": h.get("H1",""),
        "H2_core_crossfamily": h.get("H2",""), "H3_noncoronary_branch": h.get("H3",""), "H4_breadth_gradient": h.get("H4",""),
        "tail_better_model": ts.get(v,{}).get("better_model",""),
    })
wb = Workbook(); wb.remove(wb.active)
readme_sheet(wb, [
    "ST4 — CTAS score-form robustness (Results Section 2).",
    "",
    "Compares the implemented CTAS ('reference') against prespecified sensible alternative score forms and",
    "diagnostic-only forms. Source: panel_form_robustness/ (form_robustness_table, rank_concordance,",
    "tail_sensitivity, headline_verdicts). Canonical locus counts recomputed against the LOCKED locus object",
    "(panel_v2_locus_lock/): each form's significant genes mapped through the canonical gene->locus assignment",
    "(coverage shown as mappable/total; the legacy 64-style locus count is NOT used).",
    "",
    "Primary claim: gene RANKING is stable across sensible forms (Spearman vs reference 0.93-0.99).",
    "Caveat: exact set size and membership at a fixed q vary (Jaccard q05 0.35-0.78; core q001 0.06-0.50).",
    "The 250/102/32 membership is NOT form-invariant and must not be presented as such. The principal",
    "biological composition (non-coronary branch, OR ~12-20) holds across sensible forms.",
    "Alternate-form canonical locus counts (canonical_loci_q05/q01/q001) are COVERAGE-LIMITED (mappable/total shown)",
    "and are NOT a headline result. The primary robustness evidence is the gene-level rank correlation",
    "(spearman_vs_reference) and the gene/discovery counts; the reference form anchors at 104/47/16.",
    "H1-H4 = prespecified headline verdicts (calibration, core cross-family, non-coronary branch, breadth gradient).",
], "README")
wsf = wb.create_sheet("form_robustness")
cols15 = list(st15[0].keys())
write_sheet(wsf, cols15, st15, widths={"form":32,"canonical_map_coverage_q05":16},
    numeric_cols={"spearman_vs_reference","jaccard_q05_vs_ref","jaccard_core_q001_vs_ref","n_sig_q05","n_sig_q01",
        "n_sig_q001","canonical_loci_q05","canonical_loci_q01","canonical_loci_q001","core_n","core_crossfamily_prop","noncoronary_branch_OR"})
wsv = wb.create_sheet("validation_checks")
ref = st15[0]
write_sheet(wsv, ["check","expected","observed","pass"], [
    {"check":"reference is first row","expected":"reference","observed":ref["form"],"pass":ref["form"]=="reference"},
    {"check":"reference canonical loci q05","expected":104,"observed":ref["canonical_loci_q05"],"pass":str(ref["canonical_loci_q05"])=="104"},
    {"check":"reference canonical loci q001","expected":16,"observed":ref["canonical_loci_q001"],"pass":str(ref["canonical_loci_q001"])=="16"},
    {"check":"reference spearman","expected":1.0,"observed":ref["spearman_vs_reference"],"pass":float(ref["spearman_vs_reference"])==1.0},
    {"check":"n sensible+diagnostic forms","expected":7,"observed":len(st15),"pass":len(st15)==7},
], numeric_cols={"expected","observed"})
wb.save(f"{ST}/ST4_CTAS_score_form_robustness.xlsx")
print("ST4 (score-form robustness) written: forms", len(st15), "ref loci q05", ref["canonical_loci_q05"])
