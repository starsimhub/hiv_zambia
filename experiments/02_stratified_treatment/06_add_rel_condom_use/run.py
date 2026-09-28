"""
Experiment 02.06 — add structuredsexual.rel_condom_use to the calibration
                    on top of exp 02.05's 12 pars.

`rel_condom_use` was added upstream in stisim 1.7.1 (branch rc1.7.1,
commit 2d30b96): scalar multiplier on condom_data values, product
clipped to [0, 1]. Rationale: `hiv.eff_condom` posterior pinned upper
0.88 in exp 02.05 — the sampler wants more condom cooling than the
old fixed 0.5 can provide. `rel_condom_use` is a coverage lever
orthogonal to `eff_condom`'s efficacy lever.

Prior: [0.5, 1.5], guess 1.0. Lower half lets the sampler EXPAND the
current coverage (e.g. 0.7 → 1.05, clipped to 1.0); upper half also
allowed but likely to be ignored unless the sampler wants MORE
coverage than the data specifies. 13 pars total.

Reproduce: `python experiments/02_stratified_treatment/06_add_rel_condom_use/run.py`
from the repo root. Requires stisim on rc1.7.1 branch.
"""

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))
from run_hiv_calibration import run_and_save  # noqa: E402


calib_pars = {
    'hiv.beta_m2f':                     dict(low=0.008, high=0.2,  guess=0.012),
    'hiv.rel_death':                    dict(low=0.6,   high=1.6,  guess=1.0),
    'hiv.eff_condom':                   dict(low=0.5,   high=0.9,  guess=0.75),
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
