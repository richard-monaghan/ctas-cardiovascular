#!/usr/bin/env python3
"""Supplementary-figure panels SF1-SF5 from the accepted bundles (no reanalysis, no new tail model)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figlib import *
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
from matplotlib.patches import Patch
from scipy.stats import spearmanr
g = pd.read_csv(T1 + 'final_C_gene_results_seed42.tsv', sep='\t'); loc = pd.read_csv(T1 + 'final_C_locus_results_seed42.tsv', sep='\t')
GS = 'final_C_gene_results_seed42.tsv'; BLUES = LinearSegmentedColormap.from_list('b', ['#FFFFFF', SKY, BLUE, DARK])
Mnull = 1e9; floor = 1 / (Mnull + 1)
zero = g[g.null_exceedance_count == 0]; tail = g[g.p_tail_extrapolated]
def sci(x):
    e = int(np.floor(np.log10(x))); m = x / 10 ** e
    return f'{m:.1f}' + r'$\times10^{' + str(e) + r'}$'
bounds = {lab: float(g.empirical_p[g.q_BH < t].max()) for lab, t in [('q<0.05', .05), ('q<0.01', .01), ('q<0.001', .001)]}

# ================= SF1a empirical P vs CTAS (all genes)
fig, ax = plt.subplots(figsize=(3.45, 3.0)); fig.subplots_adjust(left=0.17, bottom=0.15, right=0.97, top=0.97)
ns = g[g.q_BH >= .05]; sg = g[g.q_BH < .05]
ax.scatter(ns.CTAS, -np.log10(ns.empirical_p), s=1.5, c=GREY, lw=0, rasterized=True, label='q ≥ 0.05')
ax.scatter(sg.CTAS, -np.log10(sg.empirical_p), s=4, c=VERM, lw=0, label='q < 0.05')
z0 = zero[~zero.p_tail_extrapolated]; ax.scatter(z0.CTAS, -np.log10(z0.empirical_p), s=22, facecolor='none', edgecolor=INK, lw=0.8, label='0 exceedances (at floor)')
ax.scatter(tail.CTAS, -np.log10(tail.empirical_p), s=40, marker='*', c=BLUE, edgecolor=INK, lw=0.4, label='tail estimate')
ax.axhline(-np.log10(floor), color=AXIS, lw=0.6, ls=':'); ax.text(1, -np.log10(floor) + 0.2, 'direct-resolution floor 1/(10$^9$+1)', fontsize=6.5, color=AXIS)
for lab, p in bounds.items(): ax.axhline(-np.log10(p), color='#404040', lw=0.5, ls='--'); ax.text(0.3, -np.log10(p) + 0.12, f'{lab} boundary', fontsize=6, color='#404040')
for r in tail.itertuples(): ax.annotate(r.gene, (r.CTAS, -np.log10(r.empirical_p)), xytext=(-22, -2), textcoords='offset points', fontsize=7, style='italic', ha='right')
ax.set_xlabel('CTAS'); ax.set_ylabel(r'$-\log_{10}$(empirical P)'); ax.legend(frameon=False, fontsize=6.5, loc='lower right', handletextpad=0.3)
num('SF1', 'a', 'null draws', '10^9', 'Task 1 run objects', 'M', '1e9', '1e9')
num('SF1', 'a', 'genes with zero exceedances', len(zero), GS, 'null_exceedance_count == 0', len(zero), 3)
num('SF1', 'a', 'genes with tail estimate', f'{len(tail)} ({", ".join(tail.gene)})', GS, 'p_tail_extrapolated', f'{len(tail)}:{",".join(tail.gene)}', '1:LPA')
src('SF1', 'a', 'empirical P vs CTAS for all 19,080 genes', GS, 'CTAS, empirical_p, null_exceedance_count, p_tail_extrapolated, q_BH', 'replaced (10^6 null-survival panel superseded)')
save(fig, 'SF1a', g[['gene', 'CTAS', 'null_exceedance_count', 'empirical_p', 'p_tail_extrapolated', 'q_BH']])
# SF1b BH step-up
ps = np.sort(g.empirical_p.values); m = len(ps); rk = np.arange(1, m + 1)
fig, ax = plt.subplots(figsize=(3.45, 3.0)); fig.subplots_adjust(left=0.17, bottom=0.15, right=0.97, top=0.97)
ax.loglog(rk, ps, color=BLUE, lw=1.0, label='sorted empirical P')
for (lab, t), c, n_ in zip([('q<0.05', .05), ('q<0.01', .01), ('q<0.001', .001)], [SKY, BLUE, DARK], [223, 99, 38]):
    ax.loglog(rk, t * rk / m, color=c, lw=0.7, ls='--'); k = int((g.q_BH < t).sum())
    ax.scatter([k], [ps[k - 1]], s=14, color=c, zorder=4); ax.annotate(f'{lab}: {k} genes', (k, ps[k - 1]), xytext=(7, -11), textcoords='offset points', fontsize=6.5, color=c, bbox=dict(boxstyle='square,pad=0.1', fc='white', ec='none'))
    num('SF1', 'b', f'BH crossing {lab}', f'{k} genes', GS, 'count q_BH < t', k, n_)
ax.set_xlim(1, 2e4); ax.set_xlabel('Rank (1 = smallest P)'); ax.set_ylabel('Empirical P'); ax.legend(frameon=False, fontsize=6.5, loc='lower right')
ax.text(0.04, 0.96, 'dashed: BH line q·i/m (m = 19,080)', transform=ax.transAxes, fontsize=6.5, va='top')
src('SF1', 'b', 'BH step-up plot of sorted empirical P', GS, 'sorted empirical_p vs rank; lines q*i/m', 'replaced (old tail-model comparison superseded)')
save(fig, 'SF1b', pd.DataFrame(dict(rank=rk, empirical_p=ps)))
# SF1c extreme-gene table
top = g.sort_values('CTAS', ascending=False).head(8)
fig, ax = plt.subplots(figsize=(7.087, 2.35)); ax.axis('off'); fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
rows = [[r.gene, f'{r.CTAS:.2f}', f'{int(r.null_exceedance_count)}', f'{r.empirical_p:.3g}', 'tail estimate' if r.p_tail_extrapolated else ('floor 1/(M+1)' if r.null_exceedance_count == 0 else 'direct count'), f'{r.q_BH:.3g}', r.q_tier] for r in top.itertuples()]
tb = ax.table(cellText=rows, colLabels=['Gene', 'CTAS', 'Null exceedances (of 10$^9$)', 'Empirical P', 'P basis', 'BH q', 'q-set'], loc='upper center', cellLoc='center', colLoc='center',
              colWidths=[0.11, 0.09, 0.2, 0.13, 0.15, 0.12, 0.1])
tb.auto_set_font_size(False); tb.set_fontsize(7.5); tb.scale(1, 1.25)
for (r_, c_), cell in tb.get_celld().items():
    cell.set_edgecolor('#CCCCCC'); cell.set_linewidth(0.5)
    if r_ == 0: cell.set_facecolor(BLUE); cell.get_text().set_color('white'); cell.get_text().set_fontweight('bold')
    elif r_ % 2 == 0: cell.set_facecolor('#F5F5F5')
    if c_ == 0 and r_ > 0: cell.get_text().set_fontstyle('italic')
ax.text(0.5, 0.03, f'{len(zero)} genes had zero exceedances; {len(tail)} ({", ".join(tail.gene)}) received the tail estimate. The q-set boundaries lie at empirical P {sci(bounds["q<0.05"])}, {sci(bounds["q<0.01"])} and {sci(bounds["q<0.001"])},\nfar above the direct-resolution floor 1/(10$^9$ + 1), so the tail treatment changed no q-set boundary.', ha='center', fontsize=7, transform=ax.transAxes)
for lab, p in bounds.items(): num('SF1', 'c', f'empirical-P boundary {lab}', f'{p:.2g}', GS, 'max empirical_p within set', p)
src('SF1', 'c', 'eight highest-CTAS genes: exceedances, P basis, q', GS, 'top 8 by CTAS', 'replaced (old tail-model sensitivity table superseded)')
save(fig, 'SF1c', top[['gene', 'CTAS', 'null_exceedance_count', 'empirical_p', 'p_tail_extrapolated', 'q_BH', 'q_tier']])

# ================= SF2a direct-resolution diagnostic (223 significant genes)
fig, ax = plt.subplots(figsize=(3.45, 3.0)); fig.subplots_adjust(left=0.18, bottom=0.15, right=0.97, top=0.97)
s_ = g[g.q_BH < .05].copy(); s_['c1'] = s_.null_exceedance_count + 1
dr = s_[s_.null_exceedance_count > 0]
ax.scatter(dr.CTAS, dr.c1, s=6, c=BLUE, lw=0, label=f'direct count ≥ 1 ({len(dr)})')
z0 = s_[(s_.null_exceedance_count == 0) & ~s_.p_tail_extrapolated]; ax.scatter(z0.CTAS, z0.c1, s=22, facecolor='none', edgecolor=INK, lw=0.8, label=f'0 exceedances, floor P ({len(z0)})')
tl = s_[s_.p_tail_extrapolated]; ax.scatter(tl.CTAS, tl.c1, s=40, marker='*', c=VERM, edgecolor=INK, lw=0.4, label=f'tail estimate ({len(tl)}: {", ".join(tl.gene)})')
ax.set_yscale('log'); ax.set_xlabel('CTAS (q<0.05 genes)'); ax.set_ylabel('Null exceedances + 1 (of 10$^9$ draws)')
for lab, p in bounds.items():
    y = p * (Mnull + 1); ax.axhline(y, color='#404040', lw=0.5, ls='--'); ax.text(ax.get_xlim()[1] if False else s_.CTAS.max(), y * 1.15, f'{lab} boundary', fontsize=6.5, ha='right', color='#404040')
ax.legend(frameon=False, fontsize=6.5, loc='lower left')
num('SF2', 'a', 'q<0.05 genes directly resolved (>=1 exceedance)', len(dr), GS, 'null_exceedance_count>0 among q<0.05', len(dr), 220)
num('SF2', 'a', 'q<0.05 genes at floor without tail estimate', len(z0), GS, 'count==0 & not tail', len(z0), 2)
num('SF2', 'a', 'q<0.05 genes with tail estimate', len(tl), GS, 'p_tail_extrapolated', len(tl), 1)
src('SF2', 'a', 'null exceedance counts for the 223 q<0.05 genes vs CTAS with q-set boundaries', GS, 'null_exceedance_count, CTAS, p_tail_extrapolated; boundary count = boundary P x (M+1)', 'replaced (10^6 upper-tail diagnostic superseded)')
save(fig, 'SF2a', s_[['gene', 'CTAS', 'null_exceedance_count', 'empirical_p', 'p_tail_extrapolated']])
# SF2b breadth x strength
st4 = SIGNED_ST + 'ST4_CTAS_behaviour_simulations.xlsx'
bbs = pd.read_excel(st4, sheet_name='breadth_by_strength'); fts = pd.read_excel(st4, sheet_name='fixed_total_signal'); rep_b = int(bbs.n_rep.iloc[0]); rep_f = int(fts.n_rep.iloc[0])
pv = bbs.pivot(index='per_axis_mu', columns='k_axes', values='mean_CTAS').sort_index(ascending=False)
fig, ax = plt.subplots(figsize=(3.45, 3.0)); fig.subplots_adjust(left=0.16, bottom=0.15, right=0.79, top=0.90)
im = ax.imshow(pv.values, cmap=BLUES, aspect='auto')
for i in range(pv.shape[0]):
    for j in range(pv.shape[1]): v = pv.values[i, j]; ax.text(j, i, f'{v:.1f}', ha='center', va='center', fontsize=6.5, color='white' if v > pv.values.max() * 0.6 else INK)
ax.set_xticks(range(pv.shape[1])); ax.set_xticklabels(pv.columns); ax.set_yticks(range(pv.shape[0])); ax.set_yticklabels(pv.index)
ax.set_xlabel('Independent whitened axes carrying signal (k)'); ax.set_ylabel('Per-axis strength (μ)'); ax.spines[:].set_visible(False)
cax = fig.add_axes([0.81, 0.15, 0.025, 0.75]); cb = fig.colorbar(im, cax=cax); cb.set_label('Mean CTAS', fontsize=7.5); cb.ax.tick_params(labelsize=7)
ax.set_title(f'{rep_b:,} replicates per cell', fontsize=7, pad=3)
num('SF2', 'b', 'replicates per cell', f'{rep_b:,}', 'ST4 breadth_by_strength', 'n_rep', rep_b, 20000)
src('SF2', 'b', 'mean CTAS by k and mu', 'ST4_CTAS_behaviour_simulations.xlsx (signed)', 'per_axis_mu, k_axes, mean_CTAS', 'retained with terminology updates')
save(fig, 'SF2b', bbs)
# SF2c fixed-total-signal response relative to k = 1
fig, ax = plt.subplots(figsize=(3.45, 3.0)); fig.subplots_adjust(left=0.17, bottom=0.15, right=0.97, top=0.97)
rel = []
for lam, c in zip([20, 40, 80], [SKY, BLUE, DARK]):
    s = fts[fts['lambda'] == lam].sort_values('k_axes'); r = s.mean_CTAS / s.mean_CTAS.iloc[0]; ax.plot(s.k_axes, r, '-o', ms=2.2, lw=1.2, color=c, label=str(lam))
    rel.append(pd.DataFrame(dict(lam=lam, k=s.k_axes, rel=r)))
ax.axhline(1, color=GREY, lw=0.6, ls='--'); ax.set_xticks([1, 4, 8, 12, 16]); ax.set_xlabel('Independent whitened axes carrying the signal (k)'); ax.set_ylabel('Mean CTAS relative to k = 1')
ax.legend(title='total signal λ', frameon=False, fontsize=7, title_fontsize=7, loc='upper left'); ax.text(0.97, 0.20, f'{rep_f:,} replicates per point', transform=ax.transAxes, ha='right', fontsize=7)
num('SF2', 'c', 'replicates per point', f'{rep_f:,}', 'ST4 fixed_total_signal', 'n_rep', rep_f, 30000)
src('SF2', 'c', 'mean CTAS(k)/mean CTAS(k=1)', 'ST4_CTAS_behaviour_simulations.xlsx (signed)', 'mean_CTAS by lambda, k', 'retained with terminology updates')
save(fig, 'SF2c', pd.concat(rel))
# SF2d leave-one-trait-out
st7 = T4 + 'Supplementary Tables/ST7_CTAS_leave_one_trait_out.xlsx'
lr = pd.read_excel(st7, sheet_name='locus_retention_91', header=1); ots = pd.read_excel(st7, sheet_name='omitted_trait_summary', header=1)
cnt = lr.n_runs_retained_q05.value_counts().reindex(range(17), fill_value=0)
fig, ax = plt.subplots(figsize=(3.45, 3.0)); fig.subplots_adjust(left=0.15, bottom=0.15, right=0.97, top=0.76)
ax.bar(cnt.index, cnt.values, color=[DARK if k == 16 else SKY for k in cnt.index], width=0.8)
for k, v in cnt.items():
    if v: ax.text(k, v + 0.6, str(v), ha='center', fontsize=6.5)
ax.set_xlabel('Omitted-trait runs in which the locus is retained (of 16)'); ax.set_ylabel('Final CTAS loci (n = 91)'); ax.set_xticks([0, 4, 8, 12, 16])
mn, lo, hi, s16 = ots.loci_retained_q05_of_91.mean(), ots.loci_retained_q05_of_91.min(), ots.loci_retained_q05_of_91.max(), int(cnt[16])
fig.text(0.15, 0.95, f'Per omitted trait: {lo}–{hi} of 91 loci retained (mean {mn:.1f})\n{s16} loci retained in all 16 runs', fontsize=6.8, va='top')
fig.text(0.15, 0.835, 'Precomputed sensitivity runs: raw-Spearman\nworking matrices, 10$^6$-draw nulls', fontsize=6.3, color=AXIS, va='top')
num('SF2', 'd', 'loci retained in all 16 runs', s16, 'ST7 locus_retention_91', 'n_runs_retained_q05 == 16', s16, 47)
num('SF2', 'd', 'retention range per omitted trait', f'{lo}–{hi}', 'ST7 omitted_trait_summary', 'min/max loci_retained_q05_of_91', f'{lo}-{hi}', '67-89')
num('SF2', 'd', 'mean retained per omitted trait', f'{mn:.1f}', 'ST7 omitted_trait_summary', 'mean loci_retained_q05_of_91', f'{mn:.1f}', '80.1')
src('SF2', 'd', 'distribution of omitted-trait runs retained per final locus', 'outputs/final_ctas_tables_20261021/Supplementary Tables/ST7_CTAS_leave_one_trait_out.xlsx', 'locus_retention_91.n_runs_retained_q05; omitted_trait_summary', 'rebuilt (rebased to final 91 loci)')
save(fig, 'SF2d', lr[['canonical_locus_id', 'representative_gene', 'n_runs_retained_q05']])

# ================= SF3a signed whitened profiles
prof = pd.read_csv(T4 + 'source_data/architecture_trait_profiles.tsv', sep='\t'); A = pd.read_csv(T4 + 'source_data/architecture_locus_ranked.tsv', sep='\t').set_index('canonical_locus_id')
prof['rank'] = prof.canonical_locus_id.map(A.ctas_rank); prof = prof.sort_values('rank')
tr = [t for _, ts in FAMILIES for t in ts]; W = prof[['w_' + t for t in tr]].values
vmax = float(np.ceil(np.nanpercentile(np.abs(W), 99)))
cmap = LinearSegmentedColormap.from_list('d', [VERM, '#F4C7A8', '#FFFFFF', '#A6CEE9', DARK])
fig = plt.figure(figsize=(7.087, 6.0)); ax = fig.add_axes([0.10, 0.04, 0.66, 0.81])
im = ax.imshow(W, aspect='auto', cmap=cmap, norm=TwoSlopeNorm(0, -vmax, vmax), interpolation='nearest')
ax.set_xticks(range(16)); ax.set_xticklabels([TRAIT_ABBR[t] for t in tr], rotation=90, fontsize=7); ax.xaxis.tick_top()
x0 = 0
for f, ts in FAMILIES:
    ax.axvline(x0 - 0.5, color='white', lw=1.2); ax.text(x0 + len(ts) / 2 - 0.5, -6.0 - 2.6 * (FAMILIES.index((f, ts)) % 2), FAM_LABEL[f], ha='center', va='bottom', fontsize=6.3)
    x0 += len(ts)
yt = [0, 9, 19, 29, 39, 49, 59, 69, 79, 89]; ax.set_yticks(yt); ax.set_yticklabels([str(i + 1) for i in yt], fontsize=7); ax.set_ylabel('CTAS rank of locus representative (1 = strongest)')
cax = fig.add_axes([0.80, 0.50, 0.022, 0.30]); cb = fig.colorbar(im, cax=cax, extend='both'); cb.set_label('Signed whitened coordinate', fontsize=7.5); cb.ax.tick_params(labelsize=7)
fig.text(0.79, 0.44, f'{len(prof)} loci; final ZCA-cor\nwhitener. Colour capped\nat ±{vmax:.0f} (99th percentile).\n\nNegative coordinates\nindicate evidence below\nexpectation on a whitened\naxis, not protective effects.', fontsize=6.8, va='top')
num('SF3', 'a', 'loci shown', len(prof), 'architecture_trait_profiles.tsv', 'rows', len(prof), 91)
src('SF3', 'a', 'signed whitened 16-trait profiles, rows by CTAS rank, columns grouped by family', 'outputs/final_ctas_tables_20261021/source_data/architecture_trait_profiles.tsv', 'w_<trait> columns', 'rebuilt (final-whitener recomputation)')
save(fig, 'SF3a', prof)
# SF3b Hellinger ordination
od = pd.read_csv(T4 + 'source_data/architecture_ordination.tsv', sep='\t'); sd = pd.read_csv(T4 + 'source_data/architecture_structure_diagnostics.tsv', sep='\t', keep_default_na=False)
ve = sd[sd.metric == 'PCA_var_explained_hellinger_final'].set_index('key').value.astype(float)
sil = sd[sd.metric == 'silhouette_pam_hellinger_final'].set_index('key').value.astype(float); kbest = sil.idxmax()
stab = float(sd[(sd.metric == 'boot_stability_pam_hellinger_final') & (sd.key == kbest)].value.iloc[0])
verdict = sd[sd.metric == 'verdict_hellinger_final'].key.iloc[0]
od['lead'] = od.canonical_locus_id.map(A.largest_share_family)
FCOL = {'coronary': BLUE, 'pressure': VERM, 'arrhythmia': GREEN, 'heart_failure': PURPLE, 'valve': ORANGE, 'aortic_vascular': SKY, 'cerebrovascular': '#7F7F7F', 'thrombotic': '#000000'}
fig, ax = plt.subplots(figsize=(7.087, 3.5)); fig.subplots_adjust(left=0.08, bottom=0.13, right=0.70, top=0.97)
for f, c in FCOL.items():
    s = od[od.lead == f]
    if len(s): ax.scatter(s.hellinger_PC1, s.hellinger_PC2, s=8 + (s.CTAS - 14) * 6, c=c, alpha=0.85, edgecolor='white', lw=0.4, label=f'{FAM_LABEL[f]} ({len(s)})')
for r in od.nsmallest(8, 'ctas_rank').itertuples(): ax.annotate(r.representative_gene, (r.hellinger_PC1, r.hellinger_PC2), xytext=(4, 3), textcoords='offset points', fontsize=6.5, style='italic')
ax.set_xlabel(f'Hellinger PC1 ({100 * ve["PC1"]:.1f}%)'); ax.set_ylabel(f'Hellinger PC2 ({100 * ve["PC2"]:.1f}%)'); ax.axhline(0, color='#EEEEEE', lw=0.6, zorder=0); ax.axvline(0, color='#EEEEEE', lw=0.6, zorder=0)
ax.legend(title='Leading family (n loci)', frameon=False, fontsize=6.8, title_fontsize=7, loc='upper left', bbox_to_anchor=(1.01, 1.0), markerscale=1.2)
fig.text(0.72, 0.30, f'Best PAM silhouette {sil.max():.3f} (k = {kbest[2:]});\nbootstrap stability {stab:.2f}.\nVerdict: overlapping continuum,\nno stable discrete clusters.\nPoint size: representative CTAS.', fontsize=6.8, va='top')
num('SF3', 'b', 'PC1 variance explained', f'{100 * ve["PC1"]:.1f}%', 'architecture_structure_diagnostics.tsv', 'PCA_var_explained_hellinger_final PC1', f'{100 * ve["PC1"]:.1f}')
num('SF3', 'b', 'PC2 variance explained', f'{100 * ve["PC2"]:.1f}%', 'architecture_structure_diagnostics.tsv', 'PCA_var_explained_hellinger_final PC2', f'{100 * ve["PC2"]:.1f}')
num('SF3', 'b', 'best PAM silhouette', f'{sil.max():.3f}', 'architecture_structure_diagnostics.tsv', 'max silhouette_pam_hellinger_final', f'{sil.max():.3f}')
num('SF3', 'b', 'bootstrap stability at best k', f'{stab:.2f}', 'architecture_structure_diagnostics.tsv', 'boot_stability_pam_hellinger_final', f'{stab:.2f}')
num('SF3', 'b', 'verdict', 'overlapping continuum', 'architecture_structure_diagnostics.tsv', 'verdict_hellinger_final', verdict, 'OVERLAPPING_CONTINUUM')
src('SF3', 'b', 'Hellinger PCA of family shares', 'architecture_ordination.tsv; architecture_structure_diagnostics.tsv', 'hellinger_PC1/PC2; variance explained; silhouette; stability', 'rebuilt (final-whitener recomputation)')
save(fig, 'SF3b', od)

# ================= SF4
cells = pd.read_csv(T3 + 'final_cell_colocalisation.tsv', sep='\t'); ass = cells[cells.disposition == 'assessed'].copy(); sup = pd.read_csv(T3 + 'supported_cell_detail.tsv', sep='\t')
lc = pd.read_csv(T3 + 'final_locus_colocalisation.tsv', sep='\t'); snp = pd.read_csv(T3 + 'matched_snp_diagnostic.tsv', sep='\t').set_index('set').loc['q<0.05']
# a trait x tissue
TR = [t for _, ts in FAMILIES for t in ts]; TI = list(TISSUE)
H = sup.groupby(['trait', 'tissue']).size().unstack(fill_value=0).reindex(index=TR, columns=TI, fill_value=0)
fig, ax = plt.subplots(figsize=(3.45, 3.3)); fig.subplots_adjust(left=0.17, bottom=0.05, right=0.83, top=0.80)
im = ax.imshow(H.values, cmap=BLUES, aspect='auto', vmin=0)
for i in range(H.shape[0]):
    for j in range(H.shape[1]):
        if H.values[i, j]: ax.text(j, i, str(H.values[i, j]), ha='center', va='center', fontsize=6, color='white' if H.values[i, j] > H.values.max() * 0.55 else INK)
ax.set_yticks(range(16)); ax.set_yticklabels([TRAIT_ABBR[t] for t in TR], fontsize=6.8); ax.set_xticks(range(7)); ax.set_xticklabels([TISSUE[t] for t in TI], rotation=50, ha='left', fontsize=6.8); ax.xaxis.tick_top(); ax.spines[:].set_visible(False); ax.tick_params(length=0)
cax = fig.add_axes([0.86, 0.05, 0.03, 0.75]); cb = fig.colorbar(im, cax=cax); cb.set_label(f'Cells with PP4 ≥ 0.5 (total {int(H.values.sum())})', fontsize=7); cb.ax.tick_params(labelsize=6.5)
num('SF4', 'a', 'supported cells', int(H.values.sum()), 'supported_cell_detail.tsv', 'rows (q<0.05, PP4>=0.5)', int(H.values.sum()), 150)
src('SF4', 'a', 'supported cells by trait and tissue', 'supported_cell_detail.tsv', 'count by trait, tissue', 'rebuilt (final mapped cells)')
save(fig, 'SF4a', H.reset_index())
# b PP3 vs PP4
fig, ax = plt.subplots(figsize=(3.45, 3.3)); fig.subplots_adjust(left=0.17, bottom=0.14, right=0.97, top=0.97)
ax.scatter(ass.PP3, ass.PP4, s=2, c=BLUE, alpha=0.15, lw=0, rasterized=True)
ax.axhline(0.5, color=AXIS, lw=0.6, ls='--'); ax.axvline(0.5, color=AXIS, lw=0.6, ls='--'); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.text(0.03, 0.97, 'shared variant\n(PP4 ≥ 0.5)', fontsize=6.5, va='top', color=AXIS); ax.text(0.97, 0.45, 'distinct variants\n(PP3 ≥ 0.5, PP4 < 0.5)', fontsize=6.5, va='top', ha='right', color=AXIS); ax.text(0.03, 0.03, 'neither', fontsize=6.5, color=AXIS)
ax.set_xlabel('PP3 (distinct-variant posterior)'); ax.set_ylabel('PP4 (shared-variant posterior)'); ax.text(0.97, 0.97, f'{len(ass):,} assessed cells', transform=ax.transAxes, ha='right', va='top', fontsize=6.8)
num('SF4', 'b', 'assessed cells', f'{len(ass):,}', 'final_cell_colocalisation.tsv', 'disposition == assessed', len(ass), 9056)
src('SF4', 'b', 'PP3 vs PP4 for assessed cells', 'final_cell_colocalisation.tsv', 'PP3, PP4 where disposition == assessed', 'rebuilt (final mapped cells)')
save(fig, 'SF4b', ass[['canonical_locus_id', 'tissue', 'trait', 'PP3', 'PP4']])
# c matched SNPs vs PP4
rho = spearmanr(np.log10(ass.n_snps), ass.PP4).statistic
fig, ax = plt.subplots(figsize=(3.45, 3.3)); fig.subplots_adjust(left=0.17, bottom=0.14, right=0.97, top=0.97)
ax.scatter(ass.n_snps, ass.PP4, s=2, c=BLUE, alpha=0.15, lw=0, rasterized=True); ax.set_xscale('log')
for y_ in (0.5, 0.8): ax.axhline(y_, color=AXIS, lw=0.6, ls='--')
ax.set_ylim(0, 1); ax.set_xlabel('Matched SNPs per test (log scale)'); ax.set_ylabel('PP4 (shared-variant posterior)')
ax.text(0.03, 0.97, f'Spearman ρ = {rho:.3f} (descriptive)\nN = {len(ass):,} assessed cells\n{int(ass.n_snps.min()):,}–{int(ass.n_snps.max()):,} matched SNPs', transform=ax.transAxes, va='top', fontsize=6.8)
num('SF4', 'c', 'Spearman rho log10(SNPs) vs PP4', f'{rho:.3f}', 'final_cell_colocalisation.tsv; matched_snp_diagnostic.tsv', 'spearman over assessed cells', f'{rho:.3f}', f'{snp.spearman_rho_log10_nsnps_vs_PP4:.3f}')
num('SF4', 'c', 'N assessed cells', f'{len(ass):,}', 'matched_snp_diagnostic.tsv', 'n_assessed_cells', len(ass), int(snp.n_assessed_cells))
num('SF4', 'c', 'matched SNP range', f'{int(ass.n_snps.min()):,}–{int(ass.n_snps.max()):,}', 'matched_snp_diagnostic.tsv', 'min/max', f'{int(ass.n_snps.min())}-{int(ass.n_snps.max())}', f'{int(snp.min_matched_snps)}-{int(snp.max_matched_snps)}')
src('SF4', 'c', 'matched SNPs vs PP4 (assessed cells)', 'final_cell_colocalisation.tsv', 'n_snps, PP4', 'rebuilt (final mapped cells)')
save(fig, 'SF4c', ass[['canonical_locus_id', 'tissue', 'trait', 'n_snps', 'PP4']])
# d ranked max PP4, 87 evaluable loci
evl = lc[lc.evaluable].sort_values('max_PP4', ascending=False).reset_index(drop=True); na = lc[~lc.evaluable]
fig, ax = plt.subplots(figsize=(3.45, 3.3)); fig.subplots_adjust(left=0.17, bottom=0.14, right=0.97, top=0.97)
cc = [DARK if v >= 0.8 else SKY if v >= 0.5 else LGREY for v in evl.max_PP4]; ax.scatter(range(1, len(evl) + 1), evl.max_PP4, s=8, c=cc, lw=0)
for y_ in (0.5, 0.8): ax.axhline(y_, color=AXIS, lw=0.6, ls='--')
for i, r in evl.head(6).iterrows(): ax.annotate(r.representative_gene, (i + 1, r.max_PP4), xytext=(22, 0.985 - 0.052 * i), textcoords='data', fontsize=6.3, style='italic', va='center', arrowprops=dict(arrowstyle='-', color=GREY, lw=0.4))
ax.set_ylim(0, 1.02); ax.set_xlabel(f'Evaluable loci ranked by maximum PP4 (n = {len(evl)})'); ax.set_ylabel('Maximum PP4 per locus')
ax.text(0.04, 0.04, 'Not assessed (not plotted):\n' + '\n'.join(f'{r.canonical_locus_id} {r.representative_gene}' + (' (MHC)' if r.coloc_disposition == 'not_assessed_MHC' else ' (no cis-eQTL)') for r in na.itertuples()), transform=ax.transAxes, ha='left', va='bottom', fontsize=6.3)
num('SF4', 'd', 'evaluable loci plotted', len(evl), 'final_locus_colocalisation.tsv', 'evaluable == TRUE', len(evl), 87)
num('SF4', 'd', 'not-assessed loci listed', len(na), 'final_locus_colocalisation.tsv', 'evaluable == FALSE', len(na), 4)
src('SF4', 'd', 'ranked maximum PP4 per evaluable locus; not-assessed listed', 'final_locus_colocalisation.tsv', 'max_PP4 where evaluable', 'rebuilt (final 91-locus set)')
save(fig, 'SF4d', lc[['canonical_locus_id', 'representative_gene', 'coloc_class', 'evaluable', 'max_PP4']])

# ================= SF5a UpSet of four comparators (91 loci)
fm = pd.read_csv(CT + 'outputs/four_comparator_harmonisation_20260919/ctas_104_locus_four_study_support_matrix.tsv', sep='\t', dtype=str, keep_default_na=False).set_index('canonical_locus_id').loc[loc.canonical_locus_id]
S4 = [('Qiao', 'Qiao_support'), ('Liu', 'Liu_support'), ('Foley (coordinate)', 'Foley_support'), ('PLEIO', 'PLEIO_support')]
combo = fm.apply(lambda r: tuple(int(r[c] == 'yes') for _, c in S4), axis=1); vc = combo.value_counts()
vc = vc.sort_values(ascending=False); assert vc.sum() == 91
fig = plt.figure(figsize=(7.087, 3.2)); axb = fig.add_axes([0.28, 0.42, 0.70, 0.52]); axm = fig.add_axes([0.28, 0.05, 0.70, 0.33]); axs = fig.add_axes([0.03, 0.05, 0.10, 0.33])
xs = np.arange(len(vc))
axb.bar(xs, vc.values, color=[LGREY if sum(k) == 0 else BLUE for k in vc.index], width=0.6)
for x_, v in zip(xs, vc.values): axb.text(x_, v + 0.8, str(v), ha='center', fontsize=7)
axb.set_xlim(-0.6, len(vc) - 0.4); axb.set_xticks([]); axb.set_ylabel('Final CTAS loci'); axb.spines['bottom'].set_visible(False)
axb.text(0.99, 0.95, f'{int(vc[(0, 0, 0, 0)])}/91 loci not represented in any of the four comparators', transform=axb.transAxes, ha='right', va='top', fontsize=7)
for i, (lab, _) in enumerate(S4):
    for x_, k in zip(xs, vc.index): axm.scatter(x_, i, s=22, c=INK if k[i] else '#DDDDDD', zorder=3)
for x_, k in zip(xs, vc.index):
    on = [i for i in range(4) if k[i]]
    if len(on) > 1: axm.plot([x_, x_], [min(on), max(on)], color=INK, lw=1, zorder=2)
axm.set_xlim(-0.6, len(vc) - 0.4); axm.set_ylim(3.5, -0.5); axm.set_yticks(range(4)); axm.set_yticklabels([l for l, _ in S4], fontsize=7.5); axm.set_xticks([]); axm.spines[:].set_visible(False); axm.tick_params(length=0)
tot = [int((fm[c] == 'yes').sum()) for _, c in S4]
axs.barh(range(4), tot, color=BLUE, height=0.5); axs.set_ylim(3.5, -0.5); axs.invert_xaxis(); axs.set_yticks([]); axs.spines['left'].set_visible(False); axs.spines['right'].set_visible(True); axs.spines['top'].set_visible(False)
for i, v in enumerate(tot): axs.text(v + 1, i, str(v), va='center', ha='right', fontsize=7)
axs.set_xlabel('Loci per study', fontsize=7); axs.tick_params(labelsize=6.5)
for (lab, _), v, e in zip(S4, tot, [19, 31, 12, 48]): num('SF5', 'a', f'{lab} loci (q<0.05)', v, 'ctas_104_locus_four_study_support_matrix.tsv', 'support == yes among 91', v, e)
num('SF5', 'a', 'loci in none of four', int(vc[(0, 0, 0, 0)]), 'ctas_104_locus_four_study_support_matrix.tsv', 'all four == no', int(vc[(0, 0, 0, 0)]), 33)
num('SF5', 'a', 'intersection total', int(vc.sum()), 'ctas_104_locus_four_study_support_matrix.tsv', 'sum of intersections', int(vc.sum()), 91)
for k, v in vc.items(): num('SF5', 'a', 'intersection ' + ('+'.join(l for (l, _), b in zip(S4, k) if b) or 'none'), v, 'ctas_104_locus_four_study_support_matrix.tsv', 'exact combination count', int(v))
src('SF5', 'a', 'exact four-comparator intersections over the 91 final loci', 'ctas_104_locus_four_study_support_matrix.tsv', 'Qiao/Liu/Foley(coordinate)/PLEIO support combinations', 'rebuilt (four studies, final 91 loci)')
save(fig, 'SF5a', pd.DataFrame([dict(combination='+'.join(l for (l, _), b in zip(S4, k) if b) or 'none', n_loci=int(v)) for k, v in vc.items()]))
# b gene-level FinnGen at >=1/2/3
gs = pd.read_csv(T2 + 'gene_nested_set_summary.tsv', sep='\t').set_index('set')
fig, ax = plt.subplots(figsize=(3.45, 3.0)); fig.subplots_adjust(left=0.15, bottom=0.15, right=0.98, top=0.83)
items = [('q<0.05', SKY), ('q<0.01', BLUE), ('q<0.001', DARK), ('background_q>=0.05', LGREY)]
for j, (k, c) in enumerate(items):
    for i, thr in enumerate([1, 2, 3]):
        v = 100 * gs.at[k, f'prop_ge{thr}']; ax.bar(i + (j - 1.5) * 0.2, v, 0.19, color=c, label=(k.replace('background_q>=0.05', 'background (q ≥ 0.05)') + f' (n = {int(gs.at[k, "n_evaluable"]):,})') if i == 0 else None)
        if thr == 3: ax.text(i + (j - 1.5) * 0.2, v + 1, f'{v:.0f}', ha='center', fontsize=6.3)
        num('SF5', 'b', f'{k} support >= {thr} endpoints', f'{v:.0f}%' if thr == 3 else '(bar)', 'gene_nested_set_summary.tsv', f'supported_ge{thr}/n_evaluable', f'{int(gs.at[k, f"supported_ge{thr}"])}/{int(gs.at[k, "n_evaluable"])}',
            {3: {'q<0.05': '91/218', 'q<0.01': '50/97', 'q<0.001': '19/37', 'background_q>=0.05': '1843/18326'}[k]}.get(thr))
ax.set_xticks(range(3)); ax.set_xticklabels(['≥1', '≥2', '≥3']); ax.set_xlabel('FinnGen endpoints with gene-level P<0.01'); ax.set_ylabel('CTAS genes supported (%)'); ax.set_ylim(0, 85)
ax.legend(frameon=False, fontsize=6.3, loc='lower left', bbox_to_anchor=(-0.02, 1.0), ncol=2, handlelength=1, columnspacing=0.8)
src('SF5', 'b', 'gene-level FinnGen support at >=1/2/3 endpoints', 'gene_nested_set_summary.tsv', 'prop_ge1/2/3 (evaluable genes)', 'rebuilt (gene unit)')
save(fig, 'SF5b', gs.reset_index())
# c endpoint-count distribution among 218 evaluable q<0.05 genes
gd = pd.read_csv(T2 + 'gene_support_distribution.tsv', sep='\t'); d5 = gd[gd.set == 'q<0.05']
fig, ax = plt.subplots(figsize=(3.45, 3.0)); fig.subplots_adjust(left=0.15, bottom=0.15, right=0.98, top=0.86)
ax.bar(d5.n_endpoints_p01, d5.n_genes, color=[DARK if k >= 3 else SKY for k in d5.n_endpoints_p01], width=0.75)
for r in d5.itertuples():
    if r.n_genes: ax.text(r.n_endpoints_p01, r.n_genes + 0.8, str(r.n_genes), ha='center', fontsize=6.5)
ax.set_xticks(range(11)); ax.set_xlabel('FinnGen endpoints with gene-level P<0.01 (of 10)'); ax.set_ylabel('q<0.05 CTAS genes')
n_ev = int(d5.n_genes.sum()); ge3 = int(d5[d5.n_endpoints_p01 >= 3].n_genes.sum())
ax.set_title(f'{n_ev} evaluable q<0.05 genes (5 unavailable excluded)\n{ge3} with ≥3 endpoints (dark)', fontsize=7, pad=3)
num('SF5', 'c', 'evaluable q<0.05 genes', n_ev, 'gene_support_distribution.tsv', 'sum n_genes (q<0.05)', n_ev, 218); num('SF5', 'c', 'genes with >=3 endpoints', ge3, 'gene_support_distribution.tsv', 'sum n_genes k>=3', ge3, 91)
for r in d5.itertuples(): num('SF5', 'c', f'genes with {r.n_endpoints_p01} endpoints', r.n_genes, 'gene_support_distribution.tsv', 'n_genes', int(r.n_genes))
src('SF5', 'c', 'distribution of FinnGen endpoint counts among evaluable q<0.05 genes', 'gene_support_distribution.tsv', 'set == q<0.05', 'rebuilt (gene unit)')
save(fig, 'SF5c', d5)
dump('supp'); print('supp panels done;', len(LEDGER), 'entries;', sum(r['status'] == 'FAIL' for r in LEDGER), 'FAIL')
