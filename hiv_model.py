"""
Create HIV model and interventions
"""

# %% Imports and settings
import sciris as sc
import pandas as pd
import stisim as sti
import starsim as ss
from interventions import make_hiv_intvs
ss.options.warnings = 'error'


def make_sim(verbose=1/12, analyzers=None):

    nw_pars = dict(
        prop_f0=0.75,
        prop_m0=0.6,
        f1_conc=0.1,
        m1_conc=0.2,
        p_pair_form=0.5,
        condom_data=pd.read_csv(f'data/condom_use.csv'),
    )

    hiv = sti.HIV(
        beta_m2f=0.15,
        eff_condom=0.9,
        dist_ti_init_infected = ss.uniform(low=-24, high=0),
        init_prev_data=pd.read_csv('data/init_prev_hiv.csv'),
        rel_init_prev=20.,
    )

    intvs = make_hiv_intvs()

    sim = sti.Sim(
        n_agents=10e3,
        start=1990,
        stop=2030,
        demographics='zambia',
        diseases=hiv,
        nw_pars=nw_pars,
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

    sim = make_sim()
    sim.run()

    # Export to dataframe and plot
    df = sim.to_df(resample='year', use_years=True, sep='_')  # Use dots to separate columns
    df.index = df['timevec']

    from plot_sims import plot_hiv_sims
    plot_hiv_sims(df, start_year=1990, title='hiv_plots')

    if do_save: sc.saveobj(f'results/zambia_sim.df', df)

