"""
CROI2027 recency abstract scenarios — MVP.

Three scenarios × N draws, paired seeds. Sim runs 1985-2050; analysis window 2027-2040:
  baseline   — no PN
  pn_all     — single-hop PN triggered by every new diagnosis
  pn_rtri    — single-hop PN triggered only when RTRI flags the index recent
                (sens 25% on infections < 1 yr, spec 96%)

Outputs written to results/recency_scens_2027_2040.csv: per-draw, per-scenario
cumulative 2027-2040 new infections, new diagnoses, PN attended-contacts.
"""

import numpy as np
import pandas as pd
import starsim as ss

from hiv_model import make_sim
from run_pn_scens import make_pn_pars

DRAWS_CSV = 'experiments/02_stratified_treatment/08_bigger_vls_condom_patterns/outputs/draws_used.csv'
STOP = 2051  # sim endpoint (inclusive-of-2050); analysis window narrower
PN_START = 2026  # when PN interventions turn on
ANALYSIS_START = 2027
ANALYSIS_END = 2041  # exclusive => reports 2027-2040 inclusive
PN_PARS_ENHANCED = dict(pnc=0.5, pnp=0.1, pac=0.5, pap=0.2)  # "high" from run_pn_scens

# RTRI operating characteristics (matches recency_fermi.py defaults)
RTRI_SENS = 0.25
RTRI_SPEC = 0.96


def make_rtri_eligibility(sens=RTRI_SENS, spec=RTRI_SPEC, recent_window_days=365):
    """Newly-diagnosed AND flagged recent by RTRI (sens if truly recent, 1-spec otherwise)."""
    def just_diagnosed_rtri(sim):
        hiv = sim.diseases.hiv
        new_dx = (hiv.ti_diagnosed == hiv.ti).uids
        if len(new_dx) == 0:
            return new_dx
        dt_days = (sim.ti - hiv.ti_infected[new_dx]) * sim.dt * 365
        truly_recent = dt_days < recent_window_days
        p_flag = np.where(truly_recent, sens, 1 - spec)
        flagged = np.random.random(len(new_dx)) < p_flag
        return new_dx[flagged]
    return just_diagnosed_rtri


def build_one(draw_row, scenario, seed):
    """One sim for one draw × scenario."""
    calib_pars = draw_row.drop(['index', 'mismatch']).to_dict()
    pn_pars = None
    if scenario in ('pn_all', 'pn_rtri'):
        pn_pars = make_pn_pars(**PN_PARS_ENHANCED)

    sim = make_sim(seed=seed, stop=STOP, calib_pars=calib_pars, pn_pars=pn_pars, verbose=-1)

    # Swap RTRI eligibility onto the PN intervention post-hoc
    if scenario == 'pn_rtri':
        sim.interventions['notify_partners'].eligibility = make_rtri_eligibility()

    sim.scenario = scenario
    sim.draw_idx = int(draw_row['index'])
    return sim


def _cum(result, y_start=ANALYSIS_START, y_end=ANALYSIS_END):
    df = result.to_df(resample='year', use_years=True)
    yrs = df['timevec'].dt.year
    return float(df.loc[(yrs >= y_start) & (yrs < y_end), 'value'].sum())


def extract(sim):
    """Cumulative 2027-2040 outcomes for one sim."""
    hiv = sim.results['hiv']
    out = dict(
        scenario=sim.scenario,
        draw_idx=sim.draw_idx,
        cum_new_infections=_cum(hiv['new_infections']),
        cum_new_diagnoses=_cum(hiv['new_diagnoses']),
    )
    if 'notify_partners' in sim.results:
        out['cum_attended'] = _cum(sim.results['notify_partners']['new_attended'])
    else:
        out['cum_attended'] = 0.0
    return out


def main(scenarios=('baseline', 'pn_all', 'pn_rtri'), n_draws=None, n_cpus=50):
    draws = pd.read_csv(DRAWS_CSV)
    if n_draws is not None:
        draws = draws.head(n_draws)
    print(f'Building {len(scenarios)} scenarios × {len(draws)} draws = {len(scenarios)*len(draws)} sims')

    sims = []
    for scen in scenarios:
        for i, row in draws.iterrows():
            seed = int(row['index'])  # paired seed = draw index
            sims.append(build_one(row, scen, seed))

    print(f'Running {len(sims)} sims on {n_cpus} workers…')
    sims = ss.parallel(sims, n_cpus=n_cpus).sims

    rows = [extract(s) for s in sims]
    out = pd.DataFrame(rows)
    out.to_csv('results/recency_scens_2027_2040.csv', index=False)
    print(f'Wrote results/recency_scens_2027_2040.csv ({len(out)} rows)')
    return out


if __name__ == '__main__':
    main()
