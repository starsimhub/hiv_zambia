"""
Plot ZAMPHIA 2016 by age band and sex — prevalence (top row, 5-year
bands) and annual incidence (bottom row, 3 coarser bands as reported).
Model ensemble overlay renders when the corresponding age x sex columns
are present in the calibration stats df; otherwise data-only.

Reads:
  data/zamphia_2016_hiv_by_age_sex.csv
  data/hiv_incidence_zamphia_2016.csv
  results/zam_hiv_calib_stats.df
Writes:
  figures/zamphia_hiv_age_sex.png
"""

import re

import numpy as np
import pandas as pd
import pylab as pl
import sciris as sc

from utils import percentile_pairs, set_font


ZAMPHIA_YEAR = 2016
PREV_BANDS = [(15, 20), (20, 25), (25, 30), (30, 35), (35, 40),
              (40, 45), (45, 50), (50, 55), (55, 60)]
# ZAMPHIA reports incidence at coarser bands; model 5-year bands are aggregated to match.
INC_BANDS = [(15, 25), (25, 35), (35, 50)]
INC_MODEL_BINS = {
    (15, 25): [(15, 20), (20, 25)],
    (25, 35): [(25, 30), (30, 35)],
    (35, 50): [(35, 40), (40, 45), (45, 50)],
}
_CI_RE = re.compile(r'\(([\d.]+),\s*([\d.]+)\)')


def _parse_ci(s):
    m = _CI_RE.search(str(s))
    return (float(m.group(1)), float(m.group(2))) if m else (np.nan, np.nan)


def _model_row(df_stats, year):
    idx = int(np.argmin(np.abs(df_stats.index.values - year)))
    return df_stats.iloc[idx]


def _model_prev(row, sex, stat):
    vals = []
    for ab1, ab2 in PREV_BANDS:
        col = (f'hiv.prevalence_{sex}_{ab1}_{ab2}', stat)
        if col not in row.index:
            return None
        vals.append(row[col] * 100)
    return np.array(vals)


def _model_incidence(row, sex, stat):
    """Aggregate 5-year new_infections + n_infected across the 3 ZAMPHIA bands.

    Returns annual % incidence per ZAMPHIA band, or None if columns absent.
    Approximation: annual rate ≈ new_infections_year_stratum / (n_alive_year_stratum − n_infected_year_stratum).
    Requires `hiv.new_infections_{sex}_{ab1}_{ab2}` and `hiv.n_infected_{sex}_{ab1}_{ab2}` in extras.
    """
    vals = []
    for coarse_band in INC_BANDS:
        new_inf = 0.0
        n_inf = 0.0
        n_alive = 0.0
        for ab1, ab2 in INC_MODEL_BINS[coarse_band]:
            for base, target in [('hiv.new_infections', 'new_inf'),
                                 ('hiv.n_infected', 'n_inf')]:
                col = (f'{base}_{sex}_{ab1}_{ab2}', stat)
                if col not in row.index:
                    return None
                if target == 'new_inf':
                    new_inf += row[col]
                else:
                    n_inf += row[col]
            # n_alive_{sex}_{ab1}_{ab2} is not a standard result; approximate
            # susceptibles from n_infected + PLHIV vs population is not
            # available at that resolution, so use n_infected / prevalence
            # inversion when the prevalence column is present.
            prev_col = (f'hiv.prevalence_{sex}_{ab1}_{ab2}', stat)
            if prev_col in row.index and row[prev_col] > 0:
                n_alive_band = row[(f'hiv.n_infected_{sex}_{ab1}_{ab2}', stat)] / row[prev_col]
                n_alive += n_alive_band
            else:
                return None
        n_sus = max(n_alive - n_inf, 1.0)
        vals.append(new_inf / n_sus * 100)
    return np.array(vals)


def plot(zamphia_prev, zamphia_inc, df_stats, out_path='figures/zamphia_hiv_age_sex.png'):
    set_font(size=16)
    fig, axes = pl.subplots(2, 2, figsize=(14, 9))
    row = _model_row(df_stats, ZAMPHIA_YEAR)

    # --- Row 1: prevalence, 5-year bands ---
    prev_labels = [f'{a}-{b-1}' for a, b in PREV_BANDS]
    x_prev = np.arange(len(PREV_BANDS))
    for ax, sex, label in [(axes[0, 0], 'f', 'Female'), (axes[0, 1], 'm', 'Male')]:
        sub = zamphia_prev[zamphia_prev.sex == sex].set_index('age_band').reindex(
            [f'{a}_{b-1}' for a, b in PREV_BANDS]
        )
        yerr = np.vstack([sub.prevalence - sub.ci_lo, sub.ci_hi - sub.prevalence])
        ax.errorbar(x_prev, sub.prevalence, yerr=yerr, fmt='o', color='k',
                    capsize=3, label='ZAMPHIA 2016')

        med = _model_prev(row, sex, '50%')
        if med is not None:
            line, = ax.plot(x_prev, med, label=f'Model median ({ZAMPHIA_YEAR})')
            for pair, alpha in zip(percentile_pairs, np.linspace(0.15, 0.4, len(percentile_pairs))):
                lo = _model_prev(row, sex, f'{pair[0]:.0%}')
                hi = _model_prev(row, sex, f'{pair[1]:.0%}')
                if lo is not None and hi is not None:
                    ax.fill_between(x_prev, lo, hi, alpha=alpha, facecolor=line.get_color())
        ax.set_title(f'{label} — prevalence')
        ax.set_xticks(x_prev)
        ax.set_xticklabels(prev_labels, rotation=45, ha='right')
        ax.set_ylim(bottom=0)
        ax.legend(frameon=False)
    axes[0, 0].set_ylabel('HIV prevalence (%)')

    # --- Row 2: incidence, 3 bands ---
    inc_labels = [f'{a}-{b-1}' for a, b in INC_BANDS]
    x_inc = np.arange(len(INC_BANDS))
    sex_key = {'f': ('Females percentage annual incidence', 'Females 95% CI'),
               'm': ('Males percentage annual incidence',   'Males 95% CI')}
    inc_indexed = zamphia_inc.set_index('Age')
    for ax, sex, label in [(axes[1, 0], 'f', 'Female'), (axes[1, 1], 'm', 'Male')]:
        pct_col, ci_col = sex_key[sex]
        vals, lo_vals, hi_vals = [], [], []
        for a, b in INC_BANDS:
            key = f'{a}-{b-1}'
            v = float(inc_indexed.loc[key, pct_col])
            lo, hi = _parse_ci(inc_indexed.loc[key, ci_col])
            vals.append(v); lo_vals.append(lo); hi_vals.append(hi)
        vals = np.array(vals); lo_vals = np.array(lo_vals); hi_vals = np.array(hi_vals)
        yerr = np.vstack([vals - lo_vals, hi_vals - vals])
        ax.errorbar(x_inc, vals, yerr=yerr, fmt='o', color='k',
                    capsize=3, label='ZAMPHIA 2016')

        med = _model_incidence(row, sex, '50%')
        if med is not None:
            line, = ax.plot(x_inc, med, label=f'Model median ({ZAMPHIA_YEAR})')
            for pair, alpha in zip(percentile_pairs, np.linspace(0.15, 0.4, len(percentile_pairs))):
                lo = _model_incidence(row, sex, f'{pair[0]:.0%}')
                hi = _model_incidence(row, sex, f'{pair[1]:.0%}')
                if lo is not None and hi is not None:
                    ax.fill_between(x_inc, lo, hi, alpha=alpha, facecolor=line.get_color())
        ax.set_title(f'{label} — annual incidence')
        ax.set_xticks(x_inc)
        ax.set_xticklabels(inc_labels)
        ax.set_ylim(bottom=0)
        ax.legend(frameon=False)
    axes[1, 0].set_ylabel('Annual HIV incidence (%)')
    axes[1, 0].set_xlabel('Age band')
    axes[1, 1].set_xlabel('Age band')

    sc.figlayout()
    sc.savefig(out_path, dpi=100)
    return fig


if __name__ == '__main__':
    zamphia_prev = pd.read_csv('data/zamphia_2016_hiv_by_age_sex.csv')
    zamphia_prev = zamphia_prev[zamphia_prev.sex.isin(['f', 'm'])].copy()
    zamphia_inc = pd.read_csv('data/hiv_incidence_zamphia_2016.csv')
    df_stats = sc.loadobj('results/zam_hiv_calib_stats.df')
    plot(zamphia_prev, zamphia_inc, df_stats)
    print('Done.')
