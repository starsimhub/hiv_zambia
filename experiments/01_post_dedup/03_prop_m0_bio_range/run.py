"""
Experiment 03 — constrain structuredsexual.prop_m0 to biologically plausible range.

Diff vs experiment 02:
- `structuredsexual.prop_m0` narrowed: [0.50, 0.98] → [0.60, 0.80], guess
  0.81 → 0.70. Exp 02 pegged the widened 0.98 upper (posterior mean
  0.955), which the researcher flags as biologically implausible —
  low-risk males in Zambia are ~60-80% at most. Constrain the prior and
  force the calibration to find other levers to cool transmission (or
  reveal that it cannot).

Reproduce: `python experiments/01_post_dedup/03_prop_m0_bio_range/run.py`
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
    'structuredsexual.prop_m0':     dict(low=0.60,  high=0.80, guess=0.70),
    'structuredsexual.f1_conc':     dict(low=0.01,  high=0.2,  guess=0.01),
    'structuredsexual.m1_conc':     dict(low=0.01,  high=0.2,  guess=0.01),
    'structuredsexual.p_pair_form': dict(low=0.4,   high=0.9,  guess=0.5),
}


if __name__ == '__main__':
    run_and_save(calib_pars=calib_pars, n_trials=1000, n_workers=50)
