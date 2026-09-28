"""
Experiment 02.05 — re-run 02.04's 12 pars on corrected prev 15-49 targets.

Between 02.04 and 02.05, `data/zambia_hiv_calib.csv`'s
`hiv.prevalence_15_49` column was refreshed with user-supplied UNAIDS
series (see commit `fdfa4e8`). Peak values are higher in mid-1990s /
mid-2000s (15.1% at 1999 vs prior 14%); endpoints lower (2023 = 8.9%
vs prior 9.8%). 2024-2025 targets now populated.

Same 12 calib_pars as 02.04. In-epoch.

Reproduce: `python experiments/02_stratified_treatment/05_corrected_prev_rerun/run.py`
from the repo root.
"""

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))
from run_hiv_calibration import run_and_save  # noqa: E402


# Identical to exp 02.04.
calib_pars = {
    'hiv.beta_m2f':                 dict(low=0.008, high=0.2,  guess=0.012),
    'hiv.rel_death':                dict(low=0.6,   high=1.6,  guess=1.0),
    'hiv.eff_condom':               dict(low=0.5,   high=0.9,  guess=0.75),
    'structuredsexual.prop_f0':     dict(low=0.55,  high=0.9,  guess=0.85),
    'structuredsexual.prop_m0':     dict(low=0.60,  high=0.80, guess=0.70),
    'structuredsexual.prop_f2':     dict(low=0.005, high=0.05, guess=0.01),
    'structuredsexual.prop_m2':     dict(low=0.005, high=0.05, guess=0.02),
    'structuredsexual.f1_conc':     dict(low=0.01,  high=0.2,  guess=0.01),
    'structuredsexual.m1_conc':     dict(low=0.01,  high=0.2,  guess=0.01),
    'structuredsexual.f2_conc':     dict(low=0.05,  high=0.5,  guess=0.1),
    'structuredsexual.m2_conc':     dict(low=0.2,   high=0.8,  guess=0.5),
    'structuredsexual.p_pair_form': dict(low=0.4,   high=0.9,  guess=0.5),
}


if __name__ == '__main__':
    run_and_save(calib_pars=calib_pars, n_trials=1000, n_workers=50)
