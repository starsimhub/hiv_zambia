"""
Experiment 02.08 — bigger VLS swing + condom pattern changes.

Diff vs exp 02.07 (data only; calib_pars identical):

- `data/hiv_vls_conditional_over_time.csv` — larger trajectory swing:
  2000 at ~35% of ZAMPHIA-2016 (was 50%); 2010 at ~70% (was 85%);
  2016 measured; 2020 at 92-96% (was 79-94%); 2025 at 95-97%. Bigger
  DTG-era jump should visibly cool 2020-2025 aggregate FOI.
- `data/condom_use.csv`:
  - Extended with a 2025 column reflecting continued scale-up (0.92-0.97
    across general partnerships).
  - (0,0) stable-stable floor raised from 0.01-0.05 to 0.10-0.20 across
    years condoms existed (family-planning motivation is a real,
    under-modelled use case).
  - Cross-risk partnerships (0,1), (0,2), (1,0), (1,2), (2,0), (2,1)
    raised 2020 from 0.80 to 0.90 and 2015 from 0.70 to 0.80. These
    are the pathways most likely driving the female mid-life prevalence
    residual.

Same 13 pars as 02.07. In-epoch.

Reproduce: `python experiments/02_stratified_treatment/08_bigger_vls_condom_patterns/run.py`
from the repo root. Requires stisim rc1.7.1.
"""

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))
from run_hiv_calibration import run_and_save  # noqa: E402


calib_pars = {
    'hiv.beta_m2f':                     dict(low=0.008, high=0.2,  guess=0.012),
    'hiv.rel_death':                    dict(low=0.6,   high=1.6,  guess=1.0),
    'hiv.eff_condom':                   dict(low=0.5,   high=0.95, guess=0.75),
    'structuredsexual.prop_f0':         dict(low=0.55,  high=0.9,  guess=0.85),
    'structuredsexual.prop_m0':         dict(low=0.60,  high=0.80, guess=0.70),
    'structuredsexual.prop_f2':         dict(low=0.005, high=0.05, guess=0.01),
    'structuredsexual.prop_m2':         dict(low=0.005, high=0.05, guess=0.02),
    'structuredsexual.f1_conc':         dict(low=0.01,  high=0.2,  guess=0.01),
    'structuredsexual.m1_conc':         dict(low=0.01,  high=0.2,  guess=0.01),
    'structuredsexual.f2_conc':         dict(low=0.05,  high=0.5,  guess=0.1),
    'structuredsexual.m2_conc':         dict(low=0.2,   high=0.8,  guess=0.5),
    'structuredsexual.p_pair_form':     dict(low=0.4,   high=0.9,  guess=0.5),
    'structuredsexual.rel_condom_use':  dict(low=0.5,   high=1.5,  guess=1.0),
}


if __name__ == '__main__':
    run_and_save(calib_pars=calib_pars, n_trials=1000, n_workers=50)
