#!/usr/bin/env python3
"""Main-figure panels F1-F5 from the accepted bundles (no reanalysis)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figlib import *
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.colors import LinearSegmentedColormap, Normalize
from scipy.stats import norm

g = pd.read_csv(T1 + 'final_C_gene_results_seed42.tsv', sep='\t'); loc = pd.read_csv(T1 + 'final_C_locus_results_seed42.tsv', sep='\t')
R = pd.read_csv(T1 + 'sigma_C.csv', index_col=0); ev = pd.read_csv(T1 + 'eigenvalues_C.tsv', sep='\t')
GS = 'final_C_gene_results_seed42.tsv'; LS = 'final_C_locus_results_seed42.tsv'
nG = [int((g.q_BH < t).sum()) for t in (.05, .01, .001)]
nL = [int((loc.n_q05_genes > 0).sum()), int((loc.n_q01_genes > 0).sum()), int((loc.n_q001_genes > 0).sum())]
minev = float(ev.eigenvalue.min())
BLUES = LinearSegmentedColormap.from_list('b', ['#FFFFFF', SKY, BLUE, DARK])

# ================= F1a
fig = plt.figure(figsize=(7.087, 4.35))
order = list(R.columns); ab = [TRAIT_ABBR[t] for t in order]
M = R.values.copy(); np.fill_diagonal(M, np.nan)
off = R.values[~np.eye(16, dtype=bool)]
ax = fig.add_axes([0.075, 0.25, 0.36, 0.64])
im = ax.imshow(M, cmap=BLUES, vmin=0, vmax=0.5, interpolation='nearest')
ax.set_xticks(range(16)); ax.set_xticklabels(ab, rotation=90, fontsize=7); ax.set_yticks(range(16)); ax.set_yticklabels(ab, fontsize=7)
for s in ax.spines.values(): s.set_visible(False)
ax.tick_params(length=0)
ax.set_title('Transformed trait correlation matrix R\n' + r'$R = 2\,\sin(\pi\rho_S/6)$, 16 × 16 (diagonal masked)', fontsize=8, pad=4)
cax = fig.add_axes([0.10, 0.13, 0.31, 0.022]); cb = fig.colorbar(im, cax=cax, orientation='horizontal'); cb.ax.tick_params(labelsize=7, length=2)
cb.set_label(f'Off-diagonal correlation (range {off.min():.3f} to {off.max():.3f})', fontsize=7)
num('F1', 'a', 'off-diagonal correlation range', f'{off.min():.3f} to {off.max():.3f}', 'sigma_C.csv', 'min/max off-diagonal', f'{off.min():.3f} to {off.max():.3f}')
# example gene profile (top representative), raw z and whitened
P = pd.read_csv(CT + 'analytic_null_validation/inputs/UKBB_16trait_gene_matrix.tsv', sep='\t', index_col=0)[order]
prof = pd.read_csv(T4 + 'source_data/architecture_trait_profiles.tsv', sep='\t').set_index('representative_gene')
ex = loc.sort_values('best_gene_CTAS', ascending=False).best_gene_by_final_C_CTAS.iloc[0]
z = norm.isf(P.loc[ex].clip(1e-300, 1 - 1e-10).values); z = np.nan_to_num(z, nan=0.0)
w = prof.loc[ex, ['w_' + t for t in order]].values.astype(float)
ax1 = fig.add_axes([0.53, 0.73, 0.44, 0.155]); ax1.bar(range(16), z, color=SKY, width=0.75)
ax1.set_xticks(range(16)); ax1.set_xticklabels(ab, rotation=90, fontsize=6.5); ax1.set_ylabel('z', fontsize=8); ax1.tick_params(labelsize=7)
ax1.set_title(f'Gene-level trait evidence z (example: {ex})', fontsize=8, pad=2); ax1.axhline(0, color=AXIS, lw=0.5)
fig.text(0.75, 0.585, r'ZCA-cor whitening:  $z_{white} = W z$,  $R = V\Lambda V^{\top}$,  $W = V\Lambda^{-1/2}V^{\top}$', ha='center', fontsize=8)
fig.text(0.75, 0.548, f'V: eigenvectors; Λ: eigenvalues (minimum {minev:.2f}; no floor required)', ha='center', fontsize=7, color=AXIS)
num('F1', 'a', 'primary minimum eigenvalue (paper-facing)', f'{minev:.2f}', 'eigenvalues_C.tsv', 'min eigenvalue', f'{minev:.2f}', '0.50')
ax2 = fig.add_axes([0.53, 0.29, 0.44, 0.15])
ax2.bar(range(16), np.where(w > 0, w, 0), color=BLUE, width=0.75, label='kept (positive)'); ax2.bar(range(16), np.where(w < 0, w, 0), color=LGREY, width=0.75, label='set to zero (negative)')
ax2.axhline(0, color=AXIS, lw=0.5); ax2.set_xticks(range(16)); ax2.set_xticklabels([f'{i + 1}' for i in range(16)], fontsize=6.5)
ax2.set_xlabel('Whitened axes (mixtures of traits, not phenotypes)', fontsize=7, labelpad=1); ax2.set_ylabel(r'$z_{white}$', fontsize=8); ax2.tick_params(labelsize=7)
ax2.legend(frameon=False, fontsize=6.5, loc='lower left', bbox_to_anchor=(0.0, 1.0), ncol=2, handlelength=1, borderaxespad=0.1); ax2.set_title('Positive truncation after whitening', fontsize=8, pad=14)
fig.text(0.75, 0.155, r'$w_k = \max(0, z_{white,k})^2$,   $M = \sqrt{\sum_k w_k}$,   $PR = (\sum_k w_k)^2 / \sum_k w_k^2$', ha='center', fontsize=8)
fig.text(0.75, 0.10, r'CTAS $= M \times \log_2(1 + PR)$', ha='center', fontsize=8.5, fontweight='bold', color=BLUE)
fig.text(0.75, 0.035, 'Empirical P against 10$^9$ null draws, whitened and scored identically;\nthen genome-wide Benjamini–Hochberg q', ha='center', fontsize=7.5, linespacing=1.3)
num('F1', 'a', 'primary null draws', '10^9', 'Task 1 run objects (M recorded = 1e9)', 'M', '1e9', '1e9')
src('F1', 'a', 'transformed correlation matrix; example gene z and whitened profile', 'sigma_C.csv; eigenvalues_C.tsv; UKBB_16trait_gene_matrix.tsv; architecture_trait_profiles.tsv', 'off-diagonal R; raw z = Phi^-1(1-P); w_ from final-whitener recomputation', 'rebuilt')
save(fig, 'F1a', pd.DataFrame(M, index=ab, columns=ab).reset_index().rename(columns={'index': 'trait'}))

# ================= F1b genome-wide
mp = pd.read_csv(CT + 'panel_manhattan/manhattan_points_genes.csv', usecols=['gene', 'chrn', 'x'])
d = g.merge(mp, on='gene', how='left'); assert d.x.notna().all() and len(d) == 19080
d['logp'] = -np.log10(d.empirical_p); d['sig'] = d.q_BH < .05
bounds = {t: float(d.empirical_p[d.q_BH < t].max()) for t in (.05, .01, .001)}
fig = plt.figure(figsize=(7.087, 3.35)); ax = fig.add_axes([0.075, 0.13, 0.80, 0.84])
ns = d[~d.sig]
for odd, col in [(1, '#4C78A8'), (0, '#555555')]:
    s = ns[ns.chrn % 2 == odd]; ax.scatter(s.x, s.logp, s=1.2, c=col, lw=0, rasterized=True)
s = d[d.sig]; ax.scatter(s.x, s.logp, s=4, c=VERM, lw=0, zorder=3)
xmax = d.x.max(); gut = xmax * 1.11
for (t, p), lab, n1, n2 in zip(bounds.items(), ['q<0.05', 'q<0.01', 'q<0.001'], nG, nL):
    y = -np.log10(p); ax.axhline(y, color='#404040', lw=0.6, ls='--', xmax=(xmax * 1.005) / gut, zorder=4)
    ax.text(xmax * 1.015, y, f'BH {lab}\n{n1} genes · {n2} loci', va='center', fontsize=6.5, color='#404040')
    num('F1', 'b', f'empirical-P boundary {lab}', f'{p:.3g}', GS, 'max empirical_p among genes with q_BH < t', p)
    num('F1', 'b', f'genes / loci {lab}', f'{n1} genes · {n2} loci', GS + '; ' + LS, 'count q_BH<t; loci with >=1 member q<t', f'{n1}/{n2}', {'q<0.05': '223/91', 'q<0.01': '99/46', 'q<0.001': '38/21'}[lab])
ci = d.groupby('chrn').x.agg(['min', 'max']); mids = (ci['min'] + ci['max']) / 2
ax.set_xticks(mids.values); ax.set_xticklabels(['X' if c == 23 else ('' if c in (19, 21) else str(c)) for c in ci.index], fontsize=7)
ax.set_xlim(0, gut); ax.set_ylim(0, d.logp.max() * 1.16); ax.set_xlabel('Chromosome', fontsize=9); ax.set_ylabel(r'$-\log_{10}$(empirical P)', fontsize=9)
core = loc[loc.n_q001_genes > 0]; lab_genes = core.best_gene_by_final_C_CTAS.tolist()
L = d.set_index('gene').loc[lab_genes].reset_index().sort_values('x')
# greedy label placement: stagger vertically, avoid overlaps in display space
info = ax.text(0.01, 0.98, f'19,080 tested genes; {nG[0]} at BH q<0.05 (vermillion)\nLabels: representatives of the {len(lab_genes)} q<0.001 loci', transform=ax.transAxes, va='top', fontsize=7)
fig.canvas.draw(); placed = [info.get_window_extent(fig.canvas.get_renderer()).expanded(1.05, 1.2)]
for _, r in L.iterrows():
    for dx, dy, ha in [(0, 0.55, 'center'), (0, 1.1, 'center'), (0, 1.65, 'center'), (0, 2.2, 'center'), (0.012, -0.15, 'left'), (-0.012, -0.15, 'right'), (0, -0.55, 'center'), (0, 2.75, 'center'), (0, 3.3, 'center'), (0, -1.1, 'center')]:
        ty = r.logp + dy; tx = r.x + dx * xmax
        t_ = ax.text(tx, ty, r.gene, fontsize=6.5, style='italic', ha=ha, va='bottom', zorder=6)
        bb = t_.get_window_extent(fig.canvas.get_renderer()).expanded(1.08, 1.15)
        axbb = ax.get_window_extent()
        if any(bb.overlaps(p_) for p_ in placed) or bb.y1 > axbb.y1 or bb.x1 > axbb.x0 + (axbb.x1 - axbb.x0) * (xmax / gut): t_.remove(); continue
        placed.append(bb); ax.plot([r.x, tx], [r.logp, ty], color=GREY, lw=0.4, zorder=5); break
    else: raise SystemExit('label placement failed for ' + r.gene)
num('F1', 'b', 'tested genes plotted', '19,080', GS, 'rows', len(d), 19080); num('F1', 'b', 'labelled strict-core representatives', len(lab_genes), LS, 'loci with n_q001_genes>0', len(lab_genes), 21)
src('F1', 'b', '-log10 empirical P for all genes by genomic position', GS + '; panel_manhattan/manhattan_points_genes.csv (coordinates only)', 'empirical_p, q_BH; x = cumulative position; boundaries = max p within each q-set; labels = q<0.001 locus representatives', 'rebuilt')
save(fig, 'F1b', d[['gene', 'chrn', 'x', 'empirical_p', 'q_BH', 'q_tier']])

# ================= F2a calibration (signed simulation)
st4 = SIGNED_ST + 'ST4_CTAS_behaviour_simulations.xlsx'
cal = pd.read_excel(st4, sheet_name='calibration_typeI'); fts = pd.read_excel(st4, sheet_name='fixed_total_signal'); pp = pd.read_excel(st4, sheet_name='parameters_provenance')
ncal = int(cal.n_test.iloc[0]); mref = float(pp.set_index('item').value['M_null']); rep_fts = int(fts.n_rep.iloc[0])
fig, ax = plt.subplots(figsize=(3.45, 2.9)); fig.subplots_adjust(left=0.2, bottom=0.17, right=0.97, top=0.97)
ax.plot([5e-4, 0.3], [5e-4, 0.3], color=GREY, lw=0.7, ls='--')
ax.errorbar(cal.alpha, cal.observed_typeI, yerr=[cal.observed_typeI - cal.ci_lo, cal.ci_hi - cal.observed_typeI], fmt='o', ms=3.5, color=BLUE, ecolor=BLUE, elinewidth=0.8, capsize=1.5)
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('Nominal probability (α)'); ax.set_ylabel('Observed type-I error')
ax.text(0.04, 0.96, f'{ncal:,} null genes\nreference null: {mref:,.0f} draws\nall {len(cal)} levels within MC error', transform=ax.transAxes, va='top', fontsize=7)
num('F2', 'a', 'calibration null genes', f'{ncal:,}', 'ST4 calibration_typeI', 'n_test', ncal, 200000); num('F2', 'a', 'reference null draws', f'{mref:,.0f}', 'ST4 parameters_provenance', 'M_null', mref, 1e6)
num('F2', 'a', 'levels within MC error', f'{int(cal.within_MC.sum())}/{len(cal)}', 'ST4 calibration_typeI', 'within_MC', int(cal.within_MC.sum()), 7)
src('F2', 'a', 'observed vs nominal type-I error', 'ST4_CTAS_behaviour_simulations.xlsx (signed)', 'alpha, observed_typeI, ci_lo, ci_hi', 'retained with terminology updates')
save(fig, 'F2a', cal)
# F2b fixed total signal
fig, ax = plt.subplots(figsize=(3.45, 2.9)); fig.subplots_adjust(left=0.17, bottom=0.17, right=0.97, top=0.97)
for lam, c in zip([20, 40, 80], [SKY, BLUE, DARK]):
    s = fts[fts['lambda'] == lam]; ax.plot(s.k_axes, s.mean_CTAS, '-o', ms=2.2, lw=1.2, color=c, label=f'{lam}')
ax.legend(title='total signal λ', frameon=False, fontsize=7, title_fontsize=7, loc='lower right'); ax.set_xlabel('Independent whitened axes carrying the signal (k)'); ax.set_ylabel('Mean CTAS')
ax.set_xticks([1, 4, 8, 12, 16]); ax.text(0.04, 0.96, f'{rep_fts:,} replicates per point', transform=ax.transAxes, va='top', fontsize=7)
num('F2', 'b', 'replicates per point', f'{rep_fts:,}', 'ST4 fixed_total_signal', 'n_rep', rep_fts, 30000)
src('F2', 'b', 'mean CTAS vs k at fixed total signal', 'ST4_CTAS_behaviour_simulations.xlsx (signed)', 'lambda, k_axes, mean_CTAS', 'retained with terminology updates')
save(fig, 'F2b', fts)
# F2c correlation-model and seed sensitivity
rs = pd.read_csv(T1 + 'raw_spearman_sensitivity.tsv', sep='\t', dtype=str, keep_default_na=False)
rho = float(rs[(rs.section == 'genomewide_correlation') & (rs.metric == 'spearman_CTAS_A_vs_C')].value.iloc[0])
x = rs[(rs.section == 'A_vs_C_same_seed') & (rs.seed == '42')].pivot(index='tier', columns='metric', values='value').astype(float).loc[['q05', 'q01', 'q001']]
seed = pd.read_csv(T1 + 'seed_stability_C.tsv', sep='\t'); symd = int(seed.symmetric_difference.sum())
fig, axs = plt.subplots(1, 2, figsize=(7.087, 2.85)); fig.subplots_adjust(left=0.08, right=0.99, bottom=0.17, top=0.76, wspace=0.25)
labs = ['q<0.05', 'q<0.01', 'q<0.001']; rows = []
for ax, unit in zip(axs, ['genes', 'loci']):
    pr, se, ja = x[f'{unit}_C'].values, x[f'{unit}_A'].values, x[f'{unit}_jaccard'].values; xi = np.arange(3)
    ax.bar(xi - 0.19, pr, 0.36, color=BLUE, label='Primary (transformed correlation)'); ax.bar(xi + 0.19, se, 0.36, color=LGREY, edgecolor=GREY, lw=0.5, label='Raw-Spearman sensitivity')
    for i in range(3):
        ax.text(i - 0.19, pr[i], f'{int(pr[i])}', ha='center', va='bottom', fontsize=7); ax.text(i + 0.19, se[i], f'{int(se[i])}', ha='center', va='bottom', fontsize=7)
        exp = {'genes': [(223, 240), (99, 107), (38, 45)], 'loci': [(91, 98), (46, 50), (21, 26)]}[unit][i]
        num('F2', 'c', f'{unit} {labs[i]} primary / sensitivity', f'{int(pr[i])} / {int(se[i])}', 'raw_spearman_sensitivity.tsv', f'{unit}_C, {unit}_A seed 42', f'{int(pr[i])}/{int(se[i])}', f'{exp[0]}/{exp[1]}')
        num('F2', 'c', f'{unit} {labs[i]} Jaccard', f'{ja[i]:.2f}', 'raw_spearman_sensitivity.tsv', f'{unit}_jaccard', f'{ja[i]:.2f}')
        rows.append(dict(unit=unit, q_set=labs[i], primary=int(pr[i]), sensitivity=int(se[i]), jaccard=ja[i]))
    ax.set_xticks(xi); ax.set_xticklabels([f'{labs[i]}\nJ = {ja[i]:.2f}' for i in range(3)], fontsize=7.5); ax.set_ylabel(f'CTAS {unit}'); ax.set_ylim(0, max(se) * 1.32)
axs[0].legend(frameon=False, fontsize=7, loc='upper right')
fig.text(0.08, 0.94, f'Genome-wide CTAS rank correlation, primary vs raw-Spearman sensitivity: ρ = {rho:.8f}.', fontsize=7.5)
fig.text(0.08, 0.885, 'At every threshold the primary set is a subset of the sensitivity set (J, Jaccard index).', fontsize=7.5)
fig.text(0.08, 0.83, f'Primary membership identical across seeds 42, 1 and 2 (symmetric difference {symd} at every threshold, genes and loci).', fontsize=7.5)
num('F2', 'c', 'rank rho primary vs raw-Spearman', f'{rho:.8f}', 'raw_spearman_sensitivity.tsv', 'spearman (genomewide_correlation)', f'{rho:.8f}', '0.99996171')
num('F2', 'c', 'seed membership symmetric difference', symd, 'seed_stability_C.tsv', 'sum symmetric_difference', symd, 0)
src('F2', 'c', 'primary vs raw-Spearman counts, Jaccard, rank rho; seed stability', 'raw_spearman_sensitivity.tsv; seed_stability_C.tsv', 'seed 42 rows', 'replaced (score-form panel superseded)')
save(fig, 'F2c', pd.DataFrame(rows))

# ================= F3a family composition (91 loci, final whitener)
A = pd.read_csv(T4 + 'source_data/architecture_locus_ranked.tsv', sep='\t').sort_values('ctas_rank')
S = A[['S_' + f for f, _ in FAMILIES]].values
fig = plt.figure(figsize=(7.087, 5.3)); ax = fig.add_axes([0.13, 0.05, 0.61, 0.80]); axt = fig.add_axes([0.105, 0.05, 0.018, 0.80])
im = ax.imshow(S, aspect='auto', cmap=BLUES, vmin=0, vmax=max(0.6, S.max()), interpolation='nearest')
ax.set_xticks(range(8)); ax.set_xticklabels([FAM_LABEL[f] for f, _ in FAMILIES], rotation=40, ha='left', fontsize=7.5); ax.xaxis.tick_top()
yt = [0, 9, 19, 29, 39, 49, 59, 69, 79, 89]; ax.set_yticks([]); ax.tick_params(length=2)
ax.set_ylabel('')
tier = A.locus_q_tier.map({'q<0.001': 2, 'q<0.01': 1, 'q<0.05': 0}).values[:, None]
axt.imshow(tier, aspect='auto', cmap=LinearSegmentedColormap.from_list('t', [LGREY, SKY, DARK]), vmin=0, vmax=2, interpolation='nearest'); axt.set_xticks([]); axt.set_yticks(yt); axt.set_yticklabels([str(i + 1) for i in yt], fontsize=7)
axt.set_ylabel('CTAS rank (1 = strongest representative)', fontsize=8)
cax = fig.add_axes([0.78, 0.55, 0.02, 0.28]); cb = fig.colorbar(im, cax=cax); cb.set_label('Family share $S_f$', fontsize=7.5); cb.ax.tick_params(labelsize=7)
from matplotlib.patches import Patch
fig.legend(handles=[Patch(color=DARK, label='q<0.001'), Patch(color=SKY, label='q<0.01'), Patch(color=LGREY, label='q<0.05')], title='Locus set', loc='lower left', bbox_to_anchor=(0.765, 0.30), frameon=False, fontsize=7, title_fontsize=7)
fig.text(0.765, 0.24, f'{len(A)} loci; shares from the\nfinal ZCA-cor whitened\nprofile of each representative', fontsize=7, va='top')
num('F3', 'a', 'loci shown', len(A), 'architecture_locus_ranked.tsv', 'rows', len(A), 91)
src('F3', 'a', 'eight-family share heatmap ordered by representative CTAS', 'outputs/final_ctas_tables_20261021/source_data/architecture_locus_ranked.tsv', 'S_<family> columns; ctas_rank; locus_q_tier', 'rebuilt (final-whitener recomputation)')
save(fig, 'F3a', A[['canonical_locus_id', 'ctas_rank', 'representative_gene', 'locus_q_tier'] + ['S_' + f for f, _ in FAMILIES]])
# F3b coverage
cs = pd.read_csv(T4 + 'source_data/architecture_raw_family_coverage_summary.tsv', sep='\t')
fig, ax = plt.subplots(figsize=(4.3, 3.2)); fig.subplots_adjust(left=0.15, right=0.72, bottom=0.14, top=0.97)
cats = ['0', '1', '2', '3', '4+']; cols = ['#F7F7F7', '#DCEAF5', SKY, BLUE, DARK]; sets = ['q<0.05', 'q<0.01', 'q<0.001']
for i, sname in enumerate(sets):
    sd_ = cs[cs.locus_set == sname].set_index('family_count'); bottom = 0; n_ = int(sd_.n_loci.iloc[0])
    for c, col in zip(cats, cols):
        v = float(sd_.at[c, 'pct']); k = int(sd_.at[c, 'n'])
        if v > 0:
            ax.bar(i, v, 0.62, bottom=bottom, color=col, edgecolor='white', lw=0.5)
            if v >= 6: ax.text(i, bottom + v / 2, f'{v:.0f}%', ha='center', va='center', fontsize=7, color='white' if c in ('3', '4+') else INK)
        num('F3', 'b', f'{sname} loci with {c} families', f'{v:.0f}%' if v >= 6 else '(not printed)', 'architecture_raw_family_coverage_summary.tsv', f'n/n_loci ({k}/{n_})', f'{v:.4f}')
        bottom += v
    num('F3', 'b', f'{sname} denominator', f'n = {n_}', 'architecture_raw_family_coverage_summary.tsv', 'n_loci', n_, {'q<0.05': 91, 'q<0.01': 46, 'q<0.001': 21}[sname])
ax.set_xticks(range(3)); ax.set_xticklabels([f'{s}\n(n = {int(cs[cs.locus_set == s].n_loci.iloc[0])})' for s in sets], fontsize=7.5); ax.set_ylabel('Loci (%)'); ax.set_ylim(0, 100)
ax.legend(handles=[Patch(facecolor=c_, edgecolor=GREY, lw=0.4, label=l) for c_, l in zip(cols, cats)][::-1], title='Families with\nraw P<0.01', frameon=False, fontsize=7, title_fontsize=7, loc='center left', bbox_to_anchor=(1.02, 0.5))
src('F3', 'b', 'distribution of represented families (raw P<0.01) per locus representative', 'architecture_raw_family_coverage_summary.tsv', 'pct by family_count for each nested set', 'rebuilt')
save(fig, 'F3b', cs)

# ================= F4 colocalisation
cl = pd.read_csv(T3 + 'nested_locus_class_summary.tsv', sep='\t').set_index('set'); r5 = cl.loc['q<0.05']
classes = [('strong_PP4', 'Strong (PP4 ≥ 0.8)', DARK), ('moderate_PP4', 'Moderate (0.5 ≤ PP4 < 0.8)', SKY), ('PP3_only', 'PP3-only', '#7F7F7F'), ('neither', 'Neither', LGREY), ('not_assessed', 'Not assessed', 'white')]
fig, ax = plt.subplots(figsize=(7.087, 1.45)); fig.subplots_adjust(left=0.02, right=0.98, bottom=0.42, top=0.80)
left = 0
for k, lab, col in classes:
    v = int(r5[k]); ax.barh(0, v, left=left, color=col, edgecolor=GREY if col == 'white' else 'white', lw=0.6, height=0.8, label=lab)
    ax.text(left + v / 2, 0, str(v), ha='center', va='center', fontsize=8, color='white' if col in (DARK, '#7F7F7F') else INK)
    num('F4', 'a', f'{k} loci (q<0.05)', v, 'nested_locus_class_summary.tsv', k, v, {'strong_PP4': 17, 'moderate_PP4': 14, 'PP3_only': 25, 'neither': 31, 'not_assessed': 4}[k]); left += v
ax.set_xlim(0, 91); ax.set_yticks([]); ax.set_xticks([0, 25, 50, 75, 91]); ax.spines['left'].set_visible(False); ax.set_xlabel('Final CTAS loci (q<0.05, n = 91)', fontsize=8, labelpad=1)
ax.legend(frameon=False, ncol=5, fontsize=7, loc='upper center', bbox_to_anchor=(0.5, -0.75), handlelength=1.2)
sup, ev_ = int(r5.PP4_supported_loci), int(r5.n_loci_evaluable)
ax.set_title(f'{sup}/91 loci PP4-supported ({sup}/{ev_} evaluable); 4 not assessed, excluded from the evaluable denominator', fontsize=8, pad=3)
num('F4', 'a', 'PP4-supported / all', f'{sup}/91', 'nested_locus_class_summary.tsv', 'PP4_supported_loci/n_loci_total', f'{sup}/91', '31/91'); num('F4', 'a', 'PP4-supported / evaluable', f'{sup}/{ev_}', 'nested_locus_class_summary.tsv', 'PP4_supported_loci/n_loci_evaluable', f'{sup}/{ev_}', '31/87')
src('F4', 'a', 'locus colocalisation classes, 91 loci', 'nested_locus_class_summary.tsv', 'q<0.05 row', 'rebuilt')
save(fig, 'F4a', cl.reset_index())
# F4b tissue distribution of supported loci
sc = pd.read_csv(T3 + 'supported_cell_detail.tsv', sep='\t'); lc = pd.read_csv(T3 + 'final_locus_colocalisation.tsv', sep='\t')
supl = lc[lc.coloc_class.isin(['strong_PP4', 'moderate_PP4'])].sort_values('max_PP4', ascending=False)
tis = list(TISSUE)
agg = sc.groupby(['canonical_locus_id', 'tissue']).agg(n=('PP4', 'size'), mx=('PP4', 'max')).reset_index()
fig = plt.figure(figsize=(7.087, 6.0)); ax = fig.add_axes([0.20, 0.03, 0.56, 0.86])
yl = [f'{r.representative_gene}' for r in supl.itertuples()]; ypos = {lid: i for i, lid in enumerate(supl.canonical_locus_id)}
for i in range(len(supl)): ax.axhline(i, color='#EEEEEE', lw=0.5, zorder=0)
for j in range(len(tis)): ax.axvline(j, color='#EEEEEE', lw=0.5, zorder=0)
cmap = LinearSegmentedColormap.from_list('pp', [SKY, BLUE, DARK]); normp = Normalize(0.5, 1.0)
for r in agg.itertuples():
    ax.scatter(tis.index(r.tissue), ypos[r.canonical_locus_id], s=12 + 14 * r.n, c=[cmap(normp(r.mx))], edgecolor='white', lw=0.4, zorder=3)
ax.set_yticks(range(len(supl))); ax.set_yticklabels(yl, fontsize=7, style='italic'); ax.invert_yaxis(); ax.set_xticks(range(len(tis))); ax.set_xticklabels([TISSUE[t] for t in tis], rotation=40, ha='left', fontsize=7.5); ax.xaxis.tick_top()
ax.set_xlim(-0.6, len(tis) - 0.4); ax.tick_params(length=0)
for s_ in ax.spines.values(): s_.set_visible(False)
cax = fig.add_axes([0.82, 0.62, 0.02, 0.22]); cb = fig.colorbar(plt.cm.ScalarMappable(norm=normp, cmap=cmap), cax=cax); cb.set_label('Max PP4 in tissue', fontsize=7.5); cb.ax.tick_params(labelsize=7)
for nn, yy in zip([1, 3, 6], [0.50, 0.46, 0.42]): fig.add_axes([0.80, yy, 0.03, 0.03]).axis('off'); fig.axes[-1].scatter([0.5], [0.5], s=12 + 14 * nn, c=GREY); fig.text(0.84, yy + 0.012, f'{nn} trait{"s" if nn > 1 else ""} with PP4 ≥ 0.5', fontsize=7, va='center')
fig.text(0.80, 0.37, f'{len(supl)} PP4-supported loci\n{len(sc)} cells with PP4 ≥ 0.5', fontsize=7, va='top')
num('F4', 'b', 'supported loci shown', len(supl), 'final_locus_colocalisation.tsv', 'coloc_class in strong/moderate', len(supl), 31); num('F4', 'b', 'supported cells', len(sc), 'supported_cell_detail.tsv', 'rows', len(sc), 150)
src('F4', 'b', 'tissue x locus: count of traits with PP4>=0.5 and max PP4', 'supported_cell_detail.tsv; final_locus_colocalisation.tsv', 'group by locus, tissue', 'rebuilt (final mapped cells)')
save(fig, 'F4b', agg)

# ================= F5a four-study convergence
fm = pd.read_csv(CT + 'outputs/four_comparator_harmonisation_20260919/ctas_104_locus_four_study_support_matrix.tsv', sep='\t', dtype=str, keep_default_na=False).set_index('canonical_locus_id')
li = loc.set_index('canonical_locus_id'); setsL = {'q<0.05': li.index, 'q<0.01': li.index[li.n_q01_genes > 0], 'q<0.001': li.index[li.n_q001_genes > 0]}
studies = [('Qiao', 'Qiao_support', SKY), ('Liu', 'Liu_support', '#4A9BD1'), ('Foley (coordinate)', 'Foley_support', BLUE), ('PLEIO', 'PLEIO_support', '#1F5F99'), ('Four-study union', None, VERM)]
EXP = {'q<0.05': [19, 31, 12, 48, 58], 'q<0.01': [15, 22, 11, 28, 33], 'q<0.001': [8, 13, 9, 17, 17]}
fig, ax = plt.subplots(figsize=(7.087, 3.0)); fig.subplots_adjust(left=0.08, right=0.99, bottom=0.17, top=0.86)
rows = []
for i, (sname, ids) in enumerate(setsL.items()):
    d_ = fm.loc[ids]; n_ = len(ids)
    for j, (lab, col_, c_) in enumerate(studies):
        k = int((d_[col_] == 'yes').sum()) if col_ else int(((d_.Qiao_support == 'yes') | (d_.Liu_support == 'yes') | (d_.Foley_support == 'yes') | (d_.PLEIO_support == 'yes')).sum())
        xx = i + (j - 2) * 0.16; v = pct(k, n_)
        ax.bar(xx, v, 0.145, color=c_, edgecolor=INK if col_ is None else 'none', lw=0.6, label=lab if i == 0 else None)
        ax.text(xx, v + 1.5, f'{k}', ha='center', fontsize=6.5)
        num('F5', 'a', f'{lab} {sname}', f'{k} ({v:.0f}%)', 'ctas_104_locus_four_study_support_matrix.tsv (subset to final loci)', f'{col_ or "any of four"} == yes / {n_}', k, EXP[sname][j])
        rows.append(dict(q_set=sname, n_loci=n_, study=lab, supported=k, pct=v))
ax.set_xticks(range(3)); ax.set_xticklabels([f'{s} (n = {len(v)} loci)' for s, v in setsL.items()]); ax.set_ylabel('Loci supported (%)'); ax.set_ylim(0, 100)
ax.legend(frameon=False, ncol=5, fontsize=7, loc='upper center', bbox_to_anchor=(0.5, 1.16), handlelength=1.2)
src('F5', 'a', 'per-comparator and union support by nested locus set', 'ctas_104_locus_four_study_support_matrix.tsv', 'Qiao/Liu/Foley(coordinate)/PLEIO support calls; union = any', 'rebuilt (four studies, final sets)')
save(fig, 'F5a', pd.DataFrame(rows))
# F5b FinnGen gene-level
gs = pd.read_csv(T2 + 'gene_nested_set_summary.tsv', sep='\t').set_index('set')
fig, ax = plt.subplots(figsize=(7.087, 2.9)); fig.subplots_adjust(left=0.10, right=0.99, bottom=0.16, top=0.93)
items = [('q<0.05', 'q<0.05 CTAS genes', SKY), ('q<0.01', 'q<0.01 CTAS genes', BLUE), ('q<0.001', 'q<0.001 CTAS genes', DARK), ('background_q>=0.05', 'Background genes (q ≥ 0.05)', LGREY)]
EXPF = {'q<0.05': (91, 218, 42), 'q<0.01': (50, 97, 52), 'q<0.001': (19, 37, 51), 'background_q>=0.05': (1843, 18326, 10)}
for i, (k, lab, c_) in enumerate(items):
    a_, n_ = int(gs.at[k, 'supported_ge3']), int(gs.at[k, 'n_evaluable']); v = pct(a_, n_)
    ax.bar(i, v, 0.6, color=c_); ax.text(i, v + 1.5, f'{v:.0f}%\n{a_:,}/{n_:,}', ha='center', fontsize=7.5)
    num('F5', 'b', f'FinnGen support {k}', f'{v:.0f}% ({a_:,}/{n_:,})', 'gene_nested_set_summary.tsv', 'supported_ge3/n_evaluable', f'{a_}/{n_}={round(v)}', f'{EXPF[k][0]}/{EXPF[k][1]}={EXPF[k][2]}')
ax.set_xticks(range(4)); ax.set_xticklabels([l for _, l, _ in items], fontsize=7.5); ax.set_ylabel('Genes with FinnGen support (%)'); ax.set_ylim(0, 70)
ax.text(0.99, 0.97, 'Support: gene-level P<0.01 in ≥3 of 10 FinnGen endpoints (evaluable genes)\nincreases from q<0.05 to q<0.01, then plateaus at q<0.001', transform=ax.transAxes, ha='right', va='top', fontsize=7)
src('F5', 'b', 'gene-level FinnGen support (>=3 endpoints) for nested CTAS gene sets and background', 'gene_nested_set_summary.tsv', 'supported_ge3 / n_evaluable', 'rebuilt (gene-level headline)')
save(fig, 'F5b', gs.reset_index())
dump('main')
print('main panels done;', len(LEDGER), 'ledger entries;', sum(r['status'] == 'FAIL' for r in LEDGER), 'FAIL')
