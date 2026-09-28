"""
Experiment 02.07 — time-varying ART efficacy via VLS trajectory; widened
                    eff_condom upper (0.9 → 0.95).

Two changes vs exp 02.06:
- `data/hiv_vls_conditional_over_time.csv` replaces the 2016-only
  file. Adds time-varying rows: 2000 at ~50% of ZAMPHIA-2016 values
  (early cART with poor adherence), 2010 at ~85% (cART established),
  2016 measured, 2020 at 79-94% (DTG rollout), 2025 at 82-95%. sti.ART's
  stratified coverage path interpolates per (age_bin, sex) stratum.
- `hiv.eff_condom` upper widened 0.9 → 0.95. Exp 02.06 pinned it at 0.88
  and the researcher said "could be higher".

13 pars total (same as 02.06). Motivation: infections should cool
throughout 2000-2020 (moderate) and more steeply 2020-2025 (large) —
matching UNAIDS trajectory. Time-varying VLS gives per-agent ART
efficacy that grows over time.

Reproduce: `python experiments/02_stratified_treatment/07_time_varying_vls/run.py`
from the repo root. Requires stisim on rc1.7.1 branch.
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
