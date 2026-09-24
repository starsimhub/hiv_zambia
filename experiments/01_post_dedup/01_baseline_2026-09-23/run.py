"""
Experiment 01 — baseline: 6-parameter Optuna calibration on the post-dedup epoch.

Retrospective run.py: original 2026-09-23 execution called the driver
directly (at commit 625cd89). This file captures the exact calib_pars
used so the run is reproducible from this folder alone.

Reproduce: `python experiments/01_post_dedup/01_baseline_2026-09-23/run.py`
from the repo root.
"""

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))
from run_hiv_calibration import run_and_save  # noqa: E402


calib_pars = {
    'hiv.beta_m2f':                 dict(low=0.008, high=0.02, guess=0.012),
    'structuredsexual.prop_f0':     dict(low=0.55,  high=0.9,  guess=0.85),
    'structuredsexual.prop_m0':     dict(low=0.50,  high=0.9,  guess=0.81),
    'structuredsexual.f1_conc':     dict(low=0.01,  high=0.2,  guess=0.01),
    'structuredsexual.m1_conc':     dict(low=0.01,  high=0.2,  guess=0.01),
    'structuredsexual.p_pair_form': dict(low=0.4,   high=0.9,  guess=0.5),
}


if __name__ == '__main__':
    run_and_save(calib_pars=calib_pars, n_trials=1000, n_workers=50)
