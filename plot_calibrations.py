"""
Plot the current HIV calibration ensemble against data.

Reads results/zam_hiv_calib_stats.df (per-year describe() with percentile
columns) and results/zam_hiv_par_stats.df (parameter summary). Writes
figures/hiv_calibration_fit.png. Also prints the posterior parameter
table to stdout.
"""

import sciris as sc
from plot_sims import plot_hiv_sims
from utils import percentile_pairs


if __name__ == '__main__':

    df_stats = sc.loadobj('results/zam_hiv_calib_stats.df')
    par_stats = sc.loadobj('results/zam_hiv_par_stats.df')

    plot_hiv_sims(
        df_stats,
        start_year=1985,
        end_year=2030,
        which='multi',
        percentile_pairs=percentile_pairs,
        title='hiv_calibration_fit',
    )

    # Posterior parameter table (mean, 5-95%)
    pars = [p for p in par_stats.columns if p not in ['index', 'mismatch']]
    for p in pars:
        print(f'{p}: {par_stats[p]["mean"]:.3f} ({par_stats[p]["5%"]:.3f}-{par_stats[p]["95%"]:.3f})')

    print('Done.')
