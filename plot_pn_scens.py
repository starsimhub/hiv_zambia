"""
Plot partner notification scenarios
Generate files by running run_pn_scens.py on VMs -- takes about 10min
"""

# %% Imports and settings
import pandas as pd
import sciris as sc
import matplotlib.pyplot as pl
import seaborn as sns
import utils as ut


def plot_scens(df, show=False, savefig=True):

    # Plot settings
    ut.set_font(size=20)
    fig, axes = pl.subplots(1, 2, figsize=(15, 10))
    legendfont = 14

    clist = sc.vectocolor(4)[::-1]
    colors = {'Base':'k', 'PN - low':clist[0], 'PN - med':clist[1], 'PN - high':clist[2]}
    disease = 'hiv'

    res_list = ['new_infections', 'prevalence']
    for rn, res in enumerate(res_list):
        ax = axes[rn]
        for scenario in ['Base', 'PN - low', 'PN - med', 'PN - high']:
            idx = pd.IndexSlice
            socdf = df.loc[idx[2020:2040, scenario], f'{disease}.{res}']['50%']
            socy = socdf.rolling(3, min_periods=1).mean()
            if res == 'prevalence':
                socy = socy * 100
            socdf = socdf.reset_index()
            socx = socdf.timevec.values
            ax.plot(socx, socy, label=scenario, color=colors[scenario])

        if rn == 0: ax.legend(loc='upper left', frameon=False, prop={'size': legendfont})
        ax.set_ylim(bottom=0)
        ax.set_title(f'{disease.upper()} {res.replace("new_", "")}')
        sc.SIticks(ax)

    fig.tight_layout()

    pl.savefig(f"figures/pn_scens.png", dpi=100)
    if show:
        pl.show()
    return


if __name__ == '__main__':

    show = False

    df = sc.loadobj('results/pn_scens.df')
    plot_scens(df)

    print('Done!')





