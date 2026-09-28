"""
Experiment 02.04 — open five global-FOI levers to attack the uniform ~3-4x
                    incidence overshoot: high-risk fractions, high-risk
                    concurrency, and HIV condom effectiveness.

New calibration parameters (vs 02.03):
- `structuredsexual.prop_f2`  [0.005, 0.05] — highest-risk female fraction
                                             (currently 0.01 fixed).
- `structuredsexual.prop_m2`  [0.005, 0.05] — highest-risk male fraction
                                             (currently 0.02 fixed).
- `structuredsexual.f2_conc`  [0.05, 0.5]   — concurrency, highest-risk F
                                             (currently 0.1 fixed).
- `structuredsexual.m2_conc`  [0.2, 0.8]    — concurrency, highest-risk M
                                             (currently 0.5 fixed).
- `hiv.eff_condom`            [0.5, 0.9]    — condom transmission-blocking
                                             efficacy (currently 0.5;
                                             literature ~0.7-0.85 for HIV).

Everything else unchanged from 02.03. 12 pars total. In-epoch.

Reproduce: `python experiments/02_stratified_treatment/04_open_foi_levers/run.py`
from the repo root.
"""

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))
from run_hiv_calibration import run_and_save  # noqa: E402


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
