"""
Plot the recency-test decision analysis (recency_fermi.py).
Fast -- pure arithmetic, no simulation needed.
"""

# %% Imports and settings
import numpy as np
import sciris as sc
import matplotlib.pyplot as pl
import utils as ut
import recency_fermi as rf


def plot_fermi(pars=None, show=False, savefig=True):

    ut.set_font(size=14)
    fig, axes = pl.subplots(1, 3, figsize=(20, 6.5), gridspec_kw=dict(width_ratios=[1, 1.3, 1.3]))
    pars = rf.default_pars(**sc.mergedicts(pars))
    out = rf.evaluate(pars)

    # A. Cost per additional diagnosis by strategy
    ax = axes[0]
    strats = sc.objdict(targeted='With RTRI', budget_matched='Budget-matched\n(no RTRI)', full='Full coverage\n(no RTRI)')
    colors = ['#c0392b', '#2c7fb8', '#7fcdbb']
    vals = [out[s].cost_per_dx for s in strats]
    bars = ax.bar(list(strats.values()), vals, color=colors)
    for bar, s in zip(bars, strats):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(), f'${out[s].cost_per_dx:,.0f}\n{out[s].n_dx:.1f} dx',
                ha='center', va='bottom', fontsize=11)
    ax.tick_params(axis='x', labelsize=11)
    ax.set_ylabel('Cost per additional diagnosis ($)')
    ax.set_ylim(top=max(vals) * 1.2)
    ax.set_title(f'A. Targeting costs {out.penalty:.1f}× more per diagnosis')
    sc.boxoff(ax)

    # B. Required enhanced-PN cost to justify the RTRI, across sensitivity x RTRI cost
    ax = axes[1]
    sens = np.linspace(0.05, 0.95, 91)
    cost = np.linspace(0.5, 25, 99)
    req = rf.sweep(pars, x='sens', xvals=sens, y='cost_rtri', yvals=cost, outcome='required_pn_cost')
    im = ax.pcolormesh(sens * 100, cost, np.log10(req.values), cmap='viridis_r', shading='auto')
    cs = ax.contour(sens * 100, cost, req.values, levels=[pars.cost_pn, 100, 445], colors=['w', 'w', 'w'],
                    linestyles=['-', '--', ':'])
    ax.clabel(cs, fmt=lambda v: f'${v:.0f}', fontsize=10)
    cb = fig.colorbar(im, ax=ax)
    ticks = [10, 30, 100, 300, 1000, 3000]
    cb.set_ticks(np.log10(ticks)); cb.set_ticklabels([f'${t:,}' for t in ticks])
    cb.set_label('Enhanced-PN cost needed to justify RTRI ($/index)')
    ax.plot(pars.sens * 100, pars.cost_rtri, 'r*', ms=16, label='Assumption')
    ax.set_xlabel(f'Sensitivity for infection < 1 year (%), specificity {pars.spec:.0%}')
    ax.set_ylabel('RTRI cost ($ per test)')
    ax.set_title(f'B. Targeting pays only below the ${pars.cost_pn:.0f} contour')
    ax.legend(loc='center right', frameon=False, fontsize=11, labelcolor='w')

    # C. Penalty across published assay parameter sets
    ax = axes[2]
    scens = rf.assay_scenarios(pars).sort_values('penalty')
    y = np.arange(len(scens))
    cols = ['#c0392b' if 'assumption' in l else ('#e67e22' if 'ART' in l else '#2c7fb8') for l in scens.index]
    ax.hlines(y, 1, scens.penalty, color='0.8')
    ax.scatter(scens.penalty, y, c=cols, s=60, zorder=3)
    ax.axvline(1, color='k', lw=1)
    ax.text(1.05, len(scens) - 0.5, 'break-even', fontsize=10, va='top')
    ax.set_yticks(y)
    ax.set_yticklabels([l.replace('_', ' ') for l in scens.index], fontsize=10)
    ax.set_xlim(left=0)
    ax.set_xlabel('Cost per diagnosis, targeted ÷ untargeted')
    ax.set_title('C. Targeting loses at every published parameter set')
    sc.boxoff(ax)

    fig.tight_layout()
    if savefig:
        pl.savefig('figures/recency_fermi.png', dpi=150)
    if show:
        pl.show()
    return fig


if __name__ == '__main__':
    plot_fermi()
    print('Done!')
