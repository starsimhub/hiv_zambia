"""
Plot ZAMPHIA 2016 by age band and sex — 4 panels: prevalence, annual
incidence, ART coverage, viral load suppression. Female (pink) and male
(blue) share axes per panel; model ensemble shown as box-and-whisker
per stratum with the ZAMPHIA point + 95% CI overlaid.

Reads:
  data/zamphia_2016_hiv_by_age_sex.csv (prev, 5-year bands)
  data/hiv_incidence_zamphia_2016.csv (incidence, 3 bands)
  data/hiv_treatment_status_{females,males}_age_rows.csv (ART coverage, 5-year bands)
  data/hiv_vls_conditional_zamphia_2016.csv (VLS, 3 bands)
  results/zam_hiv_calib_stats.df (model ensemble)
Writes:
  figures/zamphia_hiv_age_sex.png

Note: age-stratified ART coverage and VLS are not currently reported per
stratum by stisim's HIV module, so those two panels are data-only until
an analyzer is added.
"""

import re

import numpy as np
import pandas as pd
import pylab as pl
import sciris as sc

from utils import set_font


ZAMPHIA_YEAR = 2016
PREV_BANDS = [(15, 20), (20, 25), (25, 30), (30, 35), (35, 40),
              (40, 45), (45, 50), (50, 55), (55, 60)]
INC_BANDS = [(15, 25), (25, 35), (35, 50)]
INC_MODEL_BINS = {
    (15, 25): [(15, 20), (20, 25)],
    (25, 35): [(25, 30), (30, 35)],
    (35, 50): [(35, 40), (40, 45), (45, 50)],
}
SEX_COLOR = {'f': '#e75480', 'm': '#2b6cb0'}
SEX_LABEL = {'f': 'Female', 'm': 'Male'}
_CI_RE = re.compile(r'\(([\d.]+),\s*([\d.]+)\)')


def _parse_ci(s):
    m = _CI_RE.search(str(s))
    return (float(m.group(1)), float(m.group(2))) if m else (np.nan, np.nan)


def _model_row(df_stats, year):
    idx = int(np.argmin(np.abs(df_stats.index.values - year)))
    return df_stats.iloc[idx]


def _stats_from_row(row, cols):
    """Return a list of `dict(med, q1, q3, whislo, whishi)` per stratum, or None if missing."""
    stats = []
    for col_base in cols:
        need = [(col_base, p) for p in ('10%', '25%', '50%', '75%', '90%')]
        if any(c not in row.index for c in need):
            return None
        stats.append(dict(
            med=row[(col_base, '50%')],
            q1=row[(col_base, '25%')],
            q3=row[(col_base, '75%')],
            whislo=row[(col_base, '10%')],
            whishi=row[(col_base, '90%')],
            fliers=[],
        ))
    return stats


def _model_prev_boxes(row, sex, mult=100.0):
    cols = [f'hiv.prevalence_{sex}_{ab1}_{ab2}' for ab1, ab2 in PREV_BANDS]
    stats = _stats_from_row(row, cols)
    if stats is None:
        return None
    for s in stats:
        for k in ('med', 'q1', 'q3', 'whislo', 'whishi'):
            s[k] *= mult
    return stats


def _model_incidence_boxes(row, sex):
    """Aggregate 5-year new_infections + n_infected across the 3 ZAMPHIA bands.

    Annual incidence % ≈ new_infections_year / (n_alive_year − n_infected_year), all summed
    across the fine 5-year bins that fall inside the coarse ZAMPHIA band.
    n_alive_stratum is derived from n_infected / prevalence per fine bin.
    """
    percentiles = ('10%', '25%', '50%', '75%', '90%')
    box_keys = ('whislo', 'q1', 'med', 'q3', 'whishi')
    stats = []
    for coarse in INC_BANDS:
        vals_by_pct = {}
        for pct in percentiles:
            new_inf = 0.0
            n_sus = 0.0
            for ab1, ab2 in INC_MODEL_BINS[coarse]:
                needed = {
                    'ni': (f'hiv.new_infections_{sex}_{ab1}_{ab2}', pct),
                    'inf': (f'hiv.n_infected_{sex}_{ab1}_{ab2}', pct),
                    'prev': (f'hiv.prevalence_{sex}_{ab1}_{ab2}', pct),
                }
                if any(c not in row.index for c in needed.values()):
                    return None
                inf_i = row[needed['inf']]
                prev_i = row[needed['prev']]
                if prev_i <= 0:
                    return None
                alive_i = inf_i / prev_i
                new_inf += row[needed['ni']]
                n_sus += max(alive_i - inf_i, 0.0)
            vals_by_pct[pct] = (new_inf / max(n_sus, 1.0)) * 100
        stats.append({box_keys[i]: vals_by_pct[p] for i, p in enumerate(percentiles)} | {'fliers': []})
    return stats


def _bxp(ax, stats, positions, color):
    ax.bxp(stats, positions=positions, widths=0.35, showfliers=False, patch_artist=True,
           boxprops=dict(facecolor=color, alpha=0.35, edgecolor=color),
           whiskerprops=dict(color=color),
           capprops=dict(color=color),
           medianprops=dict(color=color, linewidth=1.5))


def _scatter_data(ax, positions, values, lo, hi, color):
    yerr = np.vstack([np.array(values) - np.array(lo), np.array(hi) - np.array(values)])
    ax.errorbar(positions, values, yerr=yerr, fmt='o', color=color, mec='k', capsize=3, zorder=5)


def _positions(n_bands, offset):
    return np.arange(n_bands) + offset


def _load_zamphia_art_stratified():
    frames = []
    for path, sex in [
        ('data/hiv_treatment_status_females_age_rows.csv', 'f'),
        ('data/hiv_treatment_status_males_age_rows.csv', 'm'),
    ]:
        df = pd.read_csv(path)
        df = df.rename(columns={df.columns[0]: 'AgeBin'})
        df = df[df['AgeBin'].str.match(r'^\d+-\d+$', na=False)].copy()
        df['On ART'] = df['On ART'].astype(str).str.replace(r'[()]', '', regex=True)
        df['pct_art'] = pd.to_numeric(df['On ART'], errors='coerce')
        df['sex'] = sex
        frames.append(df[['AgeBin', 'sex', 'pct_art']].dropna())
    return pd.concat(frames, ignore_index=True)


def plot(zamphia_prev, zamphia_inc, zamphia_art, zamphia_vls, df_stats,
         out_path='figures/zamphia_hiv_age_sex.png'):
    set_font(size=14)
    fig, axes = pl.subplots(2, 2, figsize=(15, 10))
    row = _model_row(df_stats, ZAMPHIA_YEAR)

    # --- Panel 1: prevalence, 5-year bands ---
    ax = axes[0, 0]
    prev_labels = [f'{a}-{b-1}' for a, b in PREV_BANDS]
    n = len(PREV_BANDS)
    for sex, offset in [('f', -0.2), ('m', 0.2)]:
        color = SEX_COLOR[sex]
        pos = _positions(n, offset)
        stats = _model_prev_boxes(row, sex)
        if stats is not None:
            _bxp(ax, stats, pos, color)
        sub = zamphia_prev[zamphia_prev.sex == sex].set_index('age_band').reindex(
            [f'{a}_{b-1}' for a, b in PREV_BANDS]
        )
        _scatter_data(ax, pos, sub.prevalence.tolist(),
                      sub.ci_lo.tolist(), sub.ci_hi.tolist(), color)
    ax.set_xticks(np.arange(n)); ax.set_xticklabels(prev_labels, rotation=45, ha='right')
    ax.set_title('HIV prevalence'); ax.set_ylabel('Prevalence (%)'); ax.set_ylim(bottom=0)

    # --- Panel 2: incidence, 3 bands ---
    ax = axes[0, 1]
    inc_labels = [f'{a}-{b-1}' for a, b in INC_BANDS]
    n = len(INC_BANDS)
    sex_key = {'f': ('Females percentage annual incidence', 'Females 95% CI'),
               'm': ('Males percentage annual incidence',   'Males 95% CI')}
    inc_idx = zamphia_inc.set_index('Age')
    for sex, offset in [('f', -0.15), ('m', 0.15)]:
        color = SEX_COLOR[sex]
        pos = _positions(n, offset)
        stats = _model_incidence_boxes(row, sex)
        if stats is not None:
            _bxp(ax, stats, pos, color)
        pct_col, ci_col = sex_key[sex]
        vals = []; lo = []; hi = []
        for a, b in INC_BANDS:
            key = f'{a}-{b-1}'
            vals.append(float(inc_idx.loc[key, pct_col]))
            l, h = _parse_ci(inc_idx.loc[key, ci_col]); lo.append(l); hi.append(h)
        _scatter_data(ax, pos, vals, lo, hi, color)
    ax.set_xticks(np.arange(n)); ax.set_xticklabels(inc_labels)
    ax.set_title('Annual HIV incidence'); ax.set_ylabel('Annual incidence (%)'); ax.set_ylim(bottom=0)

    # --- Panel 3: ART coverage, 5-year bands (data only) ---
    ax = axes[1, 0]
    art_labels = prev_labels
    art_bands = [f'{a}-{b-1}' for a, b in PREV_BANDS]
    n = len(art_bands)
    art_idx = zamphia_art.set_index(['sex', 'AgeBin'])['pct_art']
    for sex, offset in [('f', -0.15), ('m', 0.15)]:
        color = SEX_COLOR[sex]
        pos = _positions(n, offset)
        vals = [art_idx.get((sex, ab), np.nan) for ab in art_bands]
        mask = ~np.isnan(vals)
        ax.scatter(np.array(pos)[mask], np.array(vals)[mask], color=color, s=60,
                   edgecolor='k', label=f'{SEX_LABEL[sex]} (ZAMPHIA)')
    ax.set_xticks(np.arange(n)); ax.set_xticklabels(art_labels, rotation=45, ha='right')
    ax.set_title('ART coverage (% of PLHIV on ART)')
    ax.set_ylabel('% on ART'); ax.set_ylim(0, 100); ax.set_xlabel('Age band')
    ax.legend(frameon=False, loc='lower right')

    # --- Panel 4: viral load suppression conditional on ART, 3 bands (data only) ---
    ax = axes[1, 1]
    vls_bands = [f'{a}-{b-1}' for a, b in INC_BANDS]
    n = len(vls_bands)
    vls_idx = zamphia_vls.set_index(['Gender', 'AgeBin'])['p_vls']
    for sex, offset in [('f', -0.15), ('m', 0.15)]:
        color = SEX_COLOR[sex]
        pos = _positions(n, offset)
        vals = [vls_idx.get((sex, ab), np.nan) * 100 for ab in vls_bands]
        mask = ~np.isnan(vals)
        ax.scatter(np.array(pos)[mask], np.array(vals)[mask], color=color, s=60, edgecolor='k')
    ax.set_xticks(np.arange(n)); ax.set_xticklabels(vls_bands)
    ax.set_title('Viral load suppression (% of ART users VLS)')
    ax.set_ylabel('% VLS'); ax.set_ylim(0, 100); ax.set_xlabel('Age band')

    # Shared legend at figure level
    fig.legend(handles=[
        pl.Rectangle((0,0), 1, 1, facecolor=SEX_COLOR['f'], alpha=0.35, edgecolor=SEX_COLOR['f'], label='Female — model'),
        pl.Rectangle((0,0), 1, 1, facecolor=SEX_COLOR['m'], alpha=0.35, edgecolor=SEX_COLOR['m'], label='Male — model'),
        pl.Line2D([0],[0], marker='o', color='w', markerfacecolor=SEX_COLOR['f'], markeredgecolor='k', label='Female — ZAMPHIA', markersize=8),
        pl.Line2D([0],[0], marker='o', color='w', markerfacecolor=SEX_COLOR['m'], markeredgecolor='k', label='Male — ZAMPHIA', markersize=8),
    ], loc='upper center', ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.02))

    sc.figlayout()
    sc.savefig(out_path, dpi=100)
    return fig


if __name__ == '__main__':
    zamphia_prev = pd.read_csv('data/zamphia_2016_hiv_by_age_sex.csv')
    zamphia_prev = zamphia_prev[zamphia_prev.sex.isin(['f', 'm'])].copy()
    zamphia_inc = pd.read_csv('data/hiv_incidence_zamphia_2016.csv')
    zamphia_art = _load_zamphia_art_stratified()
    zamphia_vls = pd.read_csv('data/hiv_vls_conditional_zamphia_2016.csv')
    df_stats = sc.loadobj('results/zam_hiv_calib_stats.df')
    plot(zamphia_prev, zamphia_inc, zamphia_art, zamphia_vls, df_stats)
    print('Done.')
