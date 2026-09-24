"""
Experiment 02 — open hiv.rel_death, widen structuredsexual.prop_m0 upper bound.

Diff vs experiment 01:
- Added `hiv.rel_death` to the search space (low=0.6, high=1.6, guess=1.0)
  so mortality can be pulled up independently — exp 01 HIV deaths undershot 2x.
- Widened `structuredsexual.prop_m0` upper bound from 0.9 to 0.98 —
  exp 01 pegged that boundary in 342/500 top draws.

Reproduce: `python experiments/01_post_dedup/02_open_rel_death_widen_prop_m0/run.py`
from the repo root.
"""

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))
from run_hiv_calibration import run_and_save  # noqa: E402


calib_pars = {
    'hiv.beta_m2f':                 dict(low=0.008, high=0.02, guess=0.012),
    'hiv.rel_death':                dict(low=0.6,   high=1.6,  guess=1.0),
    'structuredsexual.prop_f0':     dict(low=0.55,  high=0.9,  guess=0.85),
    'structuredsexual.prop_m0':     dict(low=0.50,  high=0.98, guess=0.81),
    'structuredsexual.f1_conc':     dict(low=0.01,  high=0.2,  guess=0.01),
    'structuredsexual.m1_conc':     dict(low=0.01,  high=0.2,  guess=0.01),
    'structuredsexual.p_pair_form': dict(low=0.4,   high=0.9,  guess=0.5),
}


if __name__ == '__main__':
    run_and_save(calib_pars=calib_pars, n_trials=1000, n_workers=50)
