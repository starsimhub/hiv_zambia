"""
Create HIV model and interventions
"""

# %% Imports and settings
import sciris as sc
import pandas as pd
import stisim as sti
import starsim as ss
from interventions import make_hiv_intvs
# ss.options.warnings = 'error'


def make_sim(verbose=1/12, analyzers=None):

    nw = sti.StructuredSexual(
        prop_f0=0.79,
        prop_m0=0.83,
        f1_conc=0.16,
        m1_conc=0.11,
        p_pair_form=0.58,
        condom_data=pd.read_csv(f'data/condom_use.csv'),
    )

    hiv = sti.HIV(
        beta_m2f=0.035,
        eff_condom=0.95,
        init_prev_data=pd.read_csv('data/init_prev_hiv.csv'),
        rel_init_prev=1.,
    )

    intvs = make_hiv_intvs()

    sim = sti.Sim(
        n_agents=10e3,
        start=1990,
        stop=2030,
        datafolder='data/',
        demographics='zambia',
        diseases=hiv,
        networks=[nw, ss.MaternalNet()],
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

