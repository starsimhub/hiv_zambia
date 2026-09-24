"""
Experiment 02.01 — epoch 02 baseline: mostly epoch 01 exp 03 calib_pars
                    on the new stratified-treatment structural setup,
                    with hiv.beta_m2f widened 10x on the upper bound.

The epoch structural changes are documented in `../README.md`:
stratified ART coverage, stratified VLS, ART held as proportion
post-2023, and realistic condom-use scale-up. Combined they cool
transmission substantially — a smoke run at guess-pars already lands
close to data on PLHIV and prev. Widening beta_m2f upper 0.02 -> 0.2
gives the calibration room to explore higher-transmission regimes if
the sampler wants them (the cooling knobs may leave beta relatively
un-constrained now).

Reproduce: `python experiments/02_stratified_treatment/01_baseline/run.py`
from the repo root.
"""

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))
from run_hiv_calibration import run_and_save  # noqa: E402


# From epoch 01 exp 03, with hiv.beta_m2f upper widened 0.02 -> 0.2.
calib_pars = {
    'hiv.beta_m2f':                 dict(low=0.008, high=0.2, guess=0.012),
    'hiv.rel_death':                dict(low=0.6,   high=1.6,  guess=1.0),
    'structuredsexual.prop_f0':     dict(low=0.55,  high=0.9,  guess=0.85),
    'structuredsexual.prop_m0':     dict(low=0.60,  high=0.80, guess=0.70),
    'structuredsexual.f1_conc':     dict(low=0.01,  high=0.2,  guess=0.01),
    'structuredsexual.m1_conc':     dict(low=0.01,  high=0.2,  guess=0.01),
    'structuredsexual.p_pair_form': dict(low=0.4,   high=0.9,  guess=0.5),
}


if __name__ == '__main__':
    run_and_save(calib_pars=calib_pars, n_trials=1000, n_workers=50)
