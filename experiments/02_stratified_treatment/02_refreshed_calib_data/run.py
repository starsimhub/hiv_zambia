"""
Experiment 02.02 — re-run exp 02.01 pars on refreshed calibration data.

Same calib_pars as `../01_baseline/run.py` (widened beta). What changed
between 02.01 and 02.02:
- `data/zambia_hiv_calib.csv` refreshed on `origin/main` — targets fixed
  after user identified a problem with the previous data.
- `run_hiv_calibration.py` extra_results now include
  `hiv.new_infections_{sex}_{ab1}_{ab2}` and `hiv.n_infected_{sex}_{ab1}_{ab2}`
  so the ZAMPHIA age × sex incidence plot's model overlay populates.

Neither change is a model structural change (kept in-epoch), but the
target refresh limits direct fit-metric comparability with 02.01.

Reproduce: `python experiments/02_stratified_treatment/02_refreshed_calib_data/run.py`
from the repo root.
"""

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))
from run_hiv_calibration import run_and_save  # noqa: E402


# Same as exp 02.01.
calib_pars = {
    'hiv.beta_m2f':                 dict(low=0.008, high=0.2,  guess=0.012),
    'hiv.rel_death':                dict(low=0.6,   high=1.6,  guess=1.0),
    'structuredsexual.prop_f0':     dict(low=0.55,  high=0.9,  guess=0.85),
    'structuredsexual.prop_m0':     dict(low=0.60,  high=0.80, guess=0.70),
    'structuredsexual.f1_conc':     dict(low=0.01,  high=0.2,  guess=0.01),
    'structuredsexual.m1_conc':     dict(low=0.01,  high=0.2,  guess=0.01),
    'structuredsexual.p_pair_form': dict(low=0.4,   high=0.9,  guess=0.5),
}


if __name__ == '__main__':
    run_and_save(calib_pars=calib_pars, n_trials=1000, n_workers=50)
