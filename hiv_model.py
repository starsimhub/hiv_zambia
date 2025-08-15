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


def make_sim(n_agents=10e3, start=1990, stop=2030, verbose=1/12, analyzers=None):

    ####################################################################################################################
    # Diseases and interventions
    ####################################################################################################################
    hiv = make_hiv()
    diseases = [hiv]
    intvs = make_hiv_intvs()

    sim = sti.Sim(
        n_agents=n_agents,
        datafolder='data/',
        location='zambia',
        start=start,
        stop=stop,
        diseases=diseases,
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
    do_plot = True

    sim = make_sim(start=1990, stop=2041)
    sim.run()
    df = sim.to_df(resample='year', use_years=True, sep='.')  # Use dots to separate columns
    if do_save: sc.saveobj(f'results/zambia_sim.df', df)

    if do_plot:
        from plot_sims import plot_hiv_sims
        plot_hiv_sims(df)
