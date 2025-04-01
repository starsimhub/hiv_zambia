"""
Create HIV model and interventions
"""

# %% Imports and settings
import numpy as np
import starsim as ss
import sciris as sc
import pandas as pd
import stisim as sti
from interventions import make_hiv_intvs

def make_hiv():
    """ Make HIV arguments for sim"""
    hiv = sti.HIV(
        beta_m2f=0.035,
        eff_condom=0.95,
        init_prev_data=pd.read_csv('data/init_prev_hiv.csv'),
        rel_init_prev=1.,
    )
    return hiv


def make_sim(seed=1, n_agents=None, dt=1/12, start=1990, stop=2030, debug=False, verbose=1/12, analyzers=None):

    total_pop = {1970: 5.203e6, 1980: 7.05e6, 1985: 8.691e6, 1990: 9980999, 2000: 11.83e6}[start]
    if n_agents is None: n_agents = [int(10e3), int(5e2)][debug]
    if dt is None: dt = [1/12, 1][debug]

    ####################################################################################################################
    # Demographic modules
    ####################################################################################################################
    fertility_data = pd.read_csv(f'data/asfr.csv')
    pregnancy = ss.Pregnancy(unit='month', fertility_rate=fertility_data)
    death_data = pd.read_csv(f'data/deaths.csv')
    death = ss.Deaths(death_rate=death_data, rate_units=1)

    ####################################################################################################################
    # People and networks
    ####################################################################################################################
    ppl = ss.People(n_agents, age_data=pd.read_csv(f'data/age_dist_{start}.csv', index_col='age')['value'])
    sexual = sti.FastStructuredSexual(
        prop_f0=0.79,
        prop_m0=0.83,
        f1_conc=0.16,
        m1_conc=0.11,
        p_pair_form=0.58,
        condom_data=pd.read_csv(f'data/condom_use.csv'),
    )
    maternal = ss.MaternalNet(unit='month')

    ####################################################################################################################
    # Diseases and interventions
    ####################################################################################################################
    hiv = make_hiv()
    diseases = [hiv]
    intvs = make_hiv_intvs()

    sim = ss.Sim(
        dt=dt,
        rand_seed=seed,
        total_pop=total_pop,
        start=start,
        stop=stop,
        people=ppl,
        diseases=diseases,
        networks=[sexual, maternal],
        demographics=[pregnancy, death],
        interventions=intvs,
        analyzers=analyzers,
        verbose=verbose,
    )

    return sim


if __name__ == '__main__':

    # SETTINGS
    debug = False
    seed = 1
    do_save = True
    do_run = True

    sim = make_sim(start=1990, stop=2041)
    sim.run()
    df = sim.to_df(resample='year', use_years=True, sep='.')  # Use dots to separate columns
    if do_save: sc.saveobj(f'results/zambia_sim.df', df)

