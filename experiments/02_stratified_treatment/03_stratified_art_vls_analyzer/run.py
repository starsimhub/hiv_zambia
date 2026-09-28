"""
Experiment 02.03 — add stratified ART/VLS analyzer, same calib_pars as 02.02.

New `analyzers.HIVArtVlsStrat` records per-timestep n_infected / n_on_art /
n_vls by 5-year age band × sex. Extras list is expanded to include those
columns so the ZAMPHIA age × sex plot's ART coverage and VLS panels
finally get model overlays.

The analyzer is a downstream diagnostic — it does not affect the sim
dynamics or the calibration loss. `calib_pars` identical to 02.02.

Reproduce: `python experiments/02_stratified_treatment/03_stratified_art_vls_analyzer/run.py`
from the repo root.
"""

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))
from run_hiv_calibration import run_and_save  # noqa: E402


# Identical to exp 02.02.
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
