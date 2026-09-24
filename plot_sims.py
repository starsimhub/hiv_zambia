# %% Imports and settings
import sciris as sc
import pylab as pl
import numpy as np
import pandas as pd
from utils import set_font, get_y

location = 'zambia'


def _detect_model_sep(df):
    """Return '.' if model columns use dot separator, else '_'."""
    top = df.columns.get_level_values(0) if df.columns.nlevels > 1 else df.columns
    return '.' if any(str(c).startswith('hiv.') for c in top) else '_'


def plot_hiv_sims(df, start_year=2000, end_year=2025, which='single', percentile_pairs=[[.1, .99]], title='hiv_plots'):
    """ Create quantile or individual plots of HIV epi dynamics """
    set_font(size=20)
    fig, axes = pl.subplots(2, 3, figsize=(18, 7))
    axes = axes.ravel()
    alphas = np.linspace(0.2, 0.5, len(percentile_pairs))

    hiv_data = pd.read_csv(f'data/{location}_hiv_data.csv')
    hiv_data = hiv_data.loc[(hiv_data.year >= start_year) & (hiv_data.year <= end_year)]
    dfplot = df.loc[(df.index >= start_year) & (df.index <= end_year)]

    sep = _detect_model_sep(dfplot)
    def m(name):  # model column name for a given epi field
        return name if sep == '_' else name.replace('hiv_', 'hiv.', 1)

    pn = 0
    x = dfplot.index

    # Population size
    ax = axes[pn]
    resname = 'n_alive'
    ax.scatter(hiv_data.year, hiv_data[resname], color='k', label='Data')
    y = get_y(dfplot, which, resname)
    line, = ax.plot(x, y, label='Modeled')
    if which == 'multi':
        for idx, percentile_pair in enumerate(percentile_pairs):
            yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
            yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
            ax.fill_between(x, yl, yu, alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('Population size')
    ax.legend(frameon=False)
    sc.SIticks(ax)
    ax.set_ylim(bottom=0)
    pn += 1

    # PLHIV
    ax = axes[pn]
    data_col, mod_col = 'hiv_n_infected', m('hiv_n_infected')
    ax.scatter(hiv_data.year, hiv_data[data_col], label='Data', color='k')
    y = get_y(dfplot, which, mod_col)
    line, = ax.plot(x, y, label='PLHIV')
    if which == 'multi':
        for idx, percentile_pair in enumerate(percentile_pairs):
            yl = dfplot[(mod_col, f"{percentile_pair[0]:.0%}")]
            yu = dfplot[(mod_col, f"{percentile_pair[1]:.0%}")]
            ax.fill_between(x, yl, yu, alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('PLHIV')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    # HIV prevalence
    ax = axes[pn]
    data_col, mod_col = 'hiv_prevalence_15_49', m('hiv_prevalence_15_49')
    ax.scatter(hiv_data.year, hiv_data[data_col] * 100, label='Data', color='k')
    x = dfplot.index
    y = get_y(dfplot, which, mod_col)
    line, = ax.plot(x, y * 100, label='Prevalence')
    if which == 'multi':
        for idx, percentile_pair in enumerate(percentile_pairs):
            yl = dfplot[(mod_col, f"{percentile_pair[0]:.0%}")]
            yu = dfplot[(mod_col, f"{percentile_pair[1]:.0%}")]
            ax.fill_between(x, yl * 100, yu * 100, alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('HIV prevalence 15-49 (%)')
    ax.set_ylim(bottom=0)
    pn += 1

    # Infections
    ax = axes[pn]
    data_col, mod_col = 'hiv_new_infections', m('hiv_new_infections')
    ax.scatter(hiv_data.year, hiv_data[data_col], label='UNAIDS', color='k')
    x = dfplot.index
    y = get_y(dfplot, which, mod_col)
    line, = ax.plot(x, y, label='HIV infections')
    if which == 'multi':
        for idx, percentile_pair in enumerate(percentile_pairs):
            yl = dfplot[(mod_col, f"{percentile_pair[0]:.0%}")]
            yu = dfplot[(mod_col, f"{percentile_pair[1]:.0%}")]
            ax.fill_between(x, yl, yu, alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('New HIV infections/yr')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    # HIV deaths
    ax = axes[pn]
    data_col, mod_col = 'hiv_new_deaths', m('hiv_new_deaths')
    ax.scatter(hiv_data.year, hiv_data[data_col], label='UNAIDS', color='k')
    x = dfplot.index
    y = get_y(dfplot, which, mod_col)
    line, = ax.plot(x, y, label='HIV deaths')
    if which == 'multi':
        for idx, percentile_pair in enumerate(percentile_pairs):
            yl = dfplot[(mod_col, f"{percentile_pair[0]:.0%}")]
            yu = dfplot[(mod_col, f"{percentile_pair[1]:.0%}")]
            ax.fill_between(x[:-1], yl[:-1], yu[:-1], alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('HIV-related deaths')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    # On ART
    ax = axes[pn]
    data_col, mod_col = 'hiv_n_on_art', m('hiv_n_on_art')
    if data_col in hiv_data.columns:
        ax.scatter(hiv_data.year, hiv_data[data_col], color='k', label='Data')
    y = get_y(dfplot, which, mod_col)
    line, = ax.plot(x, y, label='On ART')
    if which == 'multi':
        for idx, percentile_pair in enumerate(percentile_pairs):
            yl = dfplot[(mod_col, f"{percentile_pair[0]:.0%}")]
            yu = dfplot[(mod_col, f"{percentile_pair[1]:.0%}")]
            ax.fill_between(x, yl, yu, alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('On ART')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    sc.figlayout()
    sc.savefig(f"figures/{title}.png", dpi=100)

    return fig
