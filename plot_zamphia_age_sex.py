"""
Plot ZAMPHIA 2016 HIV prevalence by age band and sex, with model
ensemble overlay at 2016 when the age x sex columns are present in the
calibration stats df.

Reads data/zamphia_2016_hiv_by_age_sex.csv and results/zam_hiv_calib_stats.df.
Writes figures/zamphia_hiv_prevalence_age_sex.png.
"""

import sciris as sc
import pandas as pd
import numpy as np
import pylab as pl
from utils import set_font, percentile_pairs


ZAMPHIA_YEAR = 2016
BANDS = [(15, 20), (20, 25), (25, 30), (30, 35), (35, 40),
         (40, 45), (45, 50), (50, 55), (55, 60)]


def _model_row(df_stats, year):
    """Return the row of df_stats closest to `year` (annualised index)."""
    idx = np.argmin(np.abs(df_stats.index.values - year))
    return df_stats.iloc[idx]


def _model_series(row, sex, stat):
    """Return prevalence % across BANDS for given sex and percentile label.

    Returns None if the age x sex columns are absent (pre-refit calibration).
    """
    vals = []
    for ab1, ab2 in BANDS:
        col = (f'hiv.prevalence_{sex}_{ab1}_{ab2}', stat)
        if col not in row.index:
            return None
        vals.append(row[col] * 100)
    return np.array(vals)


def plot(zamphia, df_stats, out_path='figures/zamphia_hiv_prevalence_age_sex.png'):
    set_font(size=18)
    fig, axes = pl.subplots(1, 2, figsize=(14, 5), sharey=True)
    row = _model_row(df_stats, ZAMPHIA_YEAR)
    band_labels = [f'{a}-{b-1}' for a, b in BANDS]
    x = np.arange(len(BANDS))

    for ax, sex, label in [(axes[0], 'f', 'Female'), (axes[1], 'm', 'Male')]:
        sub = zamphia[zamphia.sex == sex].set_index('age_band').reindex(
            [f'{a}_{b-1}' for a, b in BANDS]
        )
        yerr = np.vstack([sub.prevalence - sub.ci_lo, sub.ci_hi - sub.prevalence])
        ax.errorbar(x, sub.prevalence, yerr=yerr, fmt='o', color='k',
                    capsize=3, label='ZAMPHIA 2016')

        med = _model_series(row, sex, '50%')
        if med is not None:
            line, = ax.plot(x, med, label=f'Model median ({ZAMPHIA_YEAR})')
            for pair, alpha in zip(percentile_pairs, np.linspace(0.15, 0.4, len(percentile_pairs))):
                lo = _model_series(row, sex, f'{pair[0]:.0%}')
                hi = _model_series(row, sex, f'{pair[1]:.0%}')
                if lo is not None and hi is not None:
                    ax.fill_between(x, lo, hi, alpha=alpha, facecolor=line.get_color())

        ax.set_title(f'{label}')
        ax.set_xticks(x)
        ax.set_xticklabels(band_labels, rotation=45, ha='right')
        ax.set_xlabel('Age band')
        ax.set_ylim(bottom=0)
        ax.legend(frameon=False)

    axes[0].set_ylabel('HIV prevalence (%)')
    sc.figlayout()
    sc.savefig(out_path, dpi=100)
    return fig


if __name__ == '__main__':
    zamphia = pd.read_csv('data/zamphia_2016_hiv_by_age_sex.csv')
    zamphia = zamphia[zamphia.sex.isin(['f', 'm'])].copy()
    df_stats = sc.loadobj('results/zam_hiv_calib_stats.df')
    plot(zamphia, df_stats)
    print('Done.')
