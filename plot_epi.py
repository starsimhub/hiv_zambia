"""
Plot prevalence by age and infections by sex work.

If updates are needed:
    1. Make required updates to the model (hiv_model.py) or data folder
    2. Run run_calibration to generate the files 'results/zam_sti_calib_stats.obj'
    3. Run run_plot_data.py to generate the epi result files:
            epi_df = 'results/epi_df.df'
            sw_df = 'results/sw_df.df'
    4. Run this script to generate the figure

"""

# Import packages
import sciris as sc
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as pl
from matplotlib.gridspec import GridSpec
from utils import set_font


# %% Plotting functions
def plot_infections_by_sw(df, disease='hiv', ax=None, start_year=2010, end_year=2030):
    set_font(size=20)
    width = 0.6

    groups = {'fsw': 'FSW', 'client': 'Client', 'non_fsw': 'F', 'non_client': 'M'}
    colors = ['#d46e9c', '#2f734a', '#d46e9c', '#2f734a']
    alphas = [0.9, 0.9, 0.3, 0.3]

    si = sc.findfirst(df.index, start_year)
    ei = sc.findfirst(df.index, end_year)

    # New infections acquired by sex and sex work status
    bottom = np.zeros(2)
    x = np.array([0.5, 1.5])
    g = 0
    for group, glabel in groups.items():
        vals = np.array([
            df[f'new_infections_{group}_{disease}'][si:ei].mean(),
            df[f'new_transmissions_{group}_{disease}'][si:ei].mean(),
        ])
        p = ax.bar(x, vals, width, label=glabel, bottom=bottom, color=colors[g], alpha=alphas[g])
        ax.bar_label(p, labels=[glabel, glabel], label_type='center')
        bottom += vals
        g += 1

    ax.set_title(disease.upper()+f" infections\n{start_year}-{end_year} average")
    sc.SIticks(ax, axis='y')
    ax.set_ylim(bottom=0)
    ax.set_xticks(x, ['Acquired', 'Transmitted'])
    ax.set_xlabel('')
    ax.set_ylabel('')

    total_trans_hiv = sw_df[f'new_transmissions_fsw_hiv'][si:ei].mean()+sw_df[f'new_transmissions_client_hiv'][10:30].mean()+sw_df[f'new_transmissions_non_fsw_hiv'][10:30].mean()+sw_df[f'new_transmissions_non_client_hiv'][10:30].mean()
    print(f'HIV SW share: {(sw_df[f"new_transmissions_fsw_hiv"][si:ei].mean()+sw_df[f"new_transmissions_client_hiv"][10:30].mean())/total_trans_hiv}')

    # ax.set_ylim(bottom=0)
    return ax


# %% Run as a script
if __name__ == '__main__':

    show = False
    epi_df = sc.loadobj(f'results/epi_df.df')
    sw_df = sc.loadobj(f'results/sw_df.df')

    # Initialize plot
    set_font(size=20)
    fig, axes = pl.subplots(1, 3, figsize=(15, 6))
    color = sc.vectocolor([.5, .8, 1])[1]
    axes = axes.ravel()
    scolors = ['#ee7989', '#4682b4']

    ax = axes[0]
    ax = plot_infections_by_sw(sw_df, ax=ax)

    ax = axes[1]
    thisdf = epi_df.loc[(epi_df.age != '0-15') & (epi_df.age != '65+')].copy()
    # sns.barplot(data=thisdf, x="age", y="new_infections", hue="sex", ax=ax, palette=scolors)
    thisdf['prevalence'] *= 100
    sns.barplot(data=thisdf, x="age", y="prevalence", hue="sex", ax=ax, palette=scolors)
    ax.set_ylabel('%')
    ax.set_title('HIV prevalence\n2025')
    ax.set_xlabel('')
    ax.legend(frameon=False, prop={'size': 16})
    ax.set_ylim(bottom=0)

    ax = axes[2]
    sns.barplot(data=thisdf, x="age", y="new_infections", hue="sex", ax=ax, palette=scolors)
    ax.set_title('Mean new HIV infections/year\n2010-2030')
    ax.set_ylabel('')
    ax.set_xlabel('')
    ax.legend(frameon=False, prop={'size': 16})
    ax.set_ylim(bottom=0)

    sc.figlayout()
    sc.savefig("figures/hiv_epi.png", dpi=100)

    print('Done.')
