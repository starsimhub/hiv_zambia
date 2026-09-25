"""
Shared calibration driver for the HIV model.

Library used by per-experiment `run.py` scripts under `experiments/`.
Each experiment's `run.py` defines its own `calib_pars` dict and calls
`run_and_save(calib_pars=...)`. This file holds the machinery that
doesn't vary between experiments in the same epoch.
"""

# Additions to handle numpy multithreading
import os
os.environ.update(
    OMP_NUM_THREADS='1',
    OPENBLAS_NUM_THREADS='1',
    NUMEXPR_NUM_THREADS='1',
    MKL_NUM_THREADS='1',
)

import pandas as pd
import sciris as sc
import stisim as sti

from hiv_model import make_sim
from utils import percentiles


def run_and_save(calib_pars, n_trials=1000, n_workers=50, shrink_to=500,
                 raw_path='raw_results/zam_hiv_calib.obj',
                 stats_path='results/zam_hiv_calib_stats.df',
                 par_stats_path='results/zam_hiv_par_stats.df'):
    """Run Optuna calibration for a given calib_pars dict; save results + stats."""
    sim = make_sim(verbose=-1, use_calib=False)
    data = pd.read_csv('data/zambia_hiv_calib.csv')
    extra_results = ['hiv.n_diagnosed', 'hiv.n_on_art', 'n_alive']

    # Age x sex prevalence + counts for ZAMPHIA 2016 comparison (diagnostic only)
    age_bins = [15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 100]
    for sex in ('f', 'm'):
        for ab1, ab2 in zip(age_bins[:-1], age_bins[1:]):
            extra_results.append(f'hiv.prevalence_{sex}_{ab1}_{ab2}')
            extra_results.append(f'hiv.new_infections_{sex}_{ab1}_{ab2}')
            extra_results.append(f'hiv.n_infected_{sex}_{ab1}_{ab2}')

    calib = sti.Calibration(
        calib_pars=calib_pars,
        sim=sim,
        extra_results=extra_results,
        data=data,
        total_trials=n_trials, n_workers=n_workers,
        die=True, reseed=False, storage=None, save_results=True,
    )
    calib.calibrate(load=True)
    print(f'Best pars are {calib.best_pars}')

    print('Shrinking and saving...')
    calib = calib.shrink(n_results=shrink_to)
    sc.saveobj(raw_path, calib)

    print('Making stats...')
    df_stats = calib.resdf.groupby(calib.resdf.time).describe(percentiles=percentiles)
    sc.saveobj(stats_path, df_stats)
    par_stats = calib.df.describe(percentiles=[0.05, 0.95])
    sc.saveobj(par_stats_path, par_stats)

    return sim, calib
