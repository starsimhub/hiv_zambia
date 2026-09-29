"""
Decision-analytic (Fermi) model: does recency-test triggering earn its cost?

Companion to the agent-based scenarios. Takes a cohort of people newly
diagnosed with HIV and compares three ways of spending money on enhanced
partner notification (PN):

    targeted       -- recency test (RTRI) for everyone; enhanced PN only
                      for those flagged recent
    budget_matched -- no RTRI; the same total budget spent on enhanced PN
                      for as many indexes as it covers
    full           -- no RTRI; enhanced PN for the whole cohort

Every function is plain numpy arithmetic, so any parameter can be passed as
an array and the outputs broadcast -- that is how the sweeps work.

Usage:
    python recency_fermi.py          # validate, print summary, write results/recency_fermi_*.csv

    import recency_fermi as rf
    out = rf.evaluate(sens=0.44, spec=0.979)
    out.penalty                      # cost/dx targeted vs budget-matched

Provenance: the defaults reproduce the spreadsheet
`recency_testing_PN_model.xlsx` (see validate()). Prevalence, sensitivity,
specificity and the cost figures came from Robyn; contact positivity (40% /
17%), contacts per enhanced PN (1.0) and the undiagnosed share (15%) are
structural placeholders, not literature values. Contacts per PN and the
undiagnosed share cancel in the targeted:untargeted ratio; the contact
positivity pair does not, and drives the conclusion.
"""

# %% Imports
import numpy as np
import pandas as pd
import sciris as sc

ASSAY_FILE = 'data/recency_assays.csv'


# %% Parameters

def default_pars(**kwargs):
    """ Default parameters; override any by keyword """
    pars = sc.objdict(
        # Cohort and epidemiology
        n_cohort   = 1000,   # New HIV diagnoses modelled; only scales counts and dollars
        prevalence = 0.10,   # Adult HIV prevalence -- context only, not used in any calculation
        p_recent   = 0.10,   # Share of new diagnoses truly infected < 1 year (the prior)

        # Recency assay
        sens       = 0.25,   # P(flagged recent | infected < 1 yr)
        spec       = 0.96,   # P(flagged long-term | infected > 1 yr), among ART-naive
        cost_rtri  = 10.0,   # $ per person tested, all-in

        # Undisclosed prior ART (off by default). Treated people are often misread as recent,
        # which degrades effective specificity unless viral load is added to the algorithm.
        p_prior_art = 0.0,   # Share of apparent new diagnoses who are actually ART-experienced
        frr_art     = 0.535, # False-recent rate among ART-experienced (CEPHIA visual read)

        # Partner notification
        cost_pn    = 20.0,   # $ per index for enhanced PN
        contacts   = 1.0,    # Additional contacts reached per enhanced PN (cancels in ratio)

        # Contact-level outcomes
        pos_recent = 0.40,   # Contact HIV positivity if index recently infected
        pos_old    = 0.17,   # Contact HIV positivity if index not recently infected
        p_undx     = 0.15,   # Share of HIV+ contacts not already diagnosed (cancels in ratio)
    )
    unknown = set(kwargs) - set(pars)
    if unknown:
        errormsg = f'Unknown parameter(s): {sc.strjoin(unknown)}. Valid: {sc.strjoin(pars.keys())}'
        raise KeyError(errormsg)
    pars.update(kwargs)
    return pars


# %% Core model

def contact_positivity(q, pars):
    """ Contact positivity for a group of indexes of whom a fraction q is truly recent """
    return q * pars.pos_recent + (1 - q) * pars.pos_old


def effective_spec(pars):
    """ Specificity after mixing in undisclosed ART-experienced people """
    frr = (1 - pars.p_prior_art) * (1 - pars.spec) + pars.p_prior_art * pars.frr_art
    return 1 - frr


def test_performance(pars):
    """ 2x2 table and derived performance of the recency test """
    n, p, se = pars.n_cohort, pars.p_recent, pars.sens
    sp = effective_spec(pars)
    out = sc.objdict()
    out.tp = n * p * se             # Truly recent, flagged
    out.fn = n * p * (1 - se)       # Truly recent, missed
    out.fp = n * (1 - p) * (1 - sp) # Not recent, flagged
    out.tn = n * (1 - p) * sp       # Not recent, not flagged
    out.flagged    = out.tp + out.fp
    out.flag_rate  = out.flagged / n
    out.ppv        = out.tp / out.flagged
    out.p_rec_neg  = out.fn / (out.fn + out.tn)
    with np.errstate(divide='ignore', invalid='ignore'):
        out.lr_pos = np.divide(se, 1 - sp)  # inf for a perfectly specific test
        out.lr_neg = np.divide(1 - se, sp)
    out.spec_eff   = sp
    return out


def _strategy(n_pn, q, cost_assay, pars):
    """ Outcomes of giving enhanced PN to n_pn indexes, of whom fraction q is truly recent """
    s = sc.objdict()
    s.cost_assay   = cost_assay
    s.n_pn         = n_pn
    s.cost_pn      = n_pn * pars.cost_pn
    s.cost_total   = s.cost_assay + s.cost_pn
    s.n_contacts   = n_pn * pars.contacts
    s.q_recent     = q
    s.pos          = contact_positivity(q, pars)
    s.n_hiv        = s.n_contacts * s.pos
    s.n_dx         = s.n_hiv * pars.p_undx
    with np.errstate(divide='ignore', invalid='ignore'):
        s.cost_per_dx = s.cost_total / s.n_dx
    return s


def evaluate(pars=None, **kwargs):
    """
    Run the whole model. Returns an objdict with the test performance, the three
    strategies, and the headline breakeven/penalty quantities.

    Args:
        pars: dict of parameters (from default_pars()); keyword args override it
    """
    pars = default_pars(**sc.mergedicts(pars, kwargs))
    tp = test_performance(pars)
    n = pars.n_cohort

    tgt = _strategy(tp.flagged, tp.ppv, n * pars.cost_rtri, pars)
    bm  = _strategy(tgt.cost_total / pars.cost_pn, pars.p_recent, 0 * n, pars)
    full = _strategy(n + 0 * tp.flagged, pars.p_recent, 0 * n, pars)

    out = sc.objdict(pars=pars, test=tp, targeted=tgt, budget_matched=bm, full=full)

    # Decomposition of the penalty
    out.reach_lost   = bm.n_pn / tgt.n_pn             # Indexes served, untargeted vs targeted
    out.enrichment   = tgt.pos / bm.pos               # E: positivity bought by targeting
    out.penalty      = tgt.cost_per_dx / bm.cost_per_dx
    out.share_doing  = tgt.cost_pn / tgt.cost_total   # Share of targeted budget that reaches PN

    # Breakeven, closed form: targeting wins iff cost_rtri < f * cost_pn * (E - 1)
    f, E = tp.flag_rate, out.enrichment
    out.breakeven_rtri_cost = f * pars.cost_pn * (E - 1)
    with np.errstate(divide='ignore', invalid='ignore'):
        out.required_pn_cost = pars.cost_rtri / (f * (E - 1))
    out.targeting_wins = pars.cost_rtri < out.breakeven_rtri_cost

    # Decision window: how different are the test-positive and test-negative groups?
    pos_neg = contact_positivity(tp.p_rec_neg, pars)
    out.neg_vs_pos = pos_neg / tgt.pos
    out.window     = (tgt.pos - pos_neg) / bm.pos
    return out


# %% Tables

def comparison_table(out):
    """ The three-column strategy comparison, as in the spreadsheet """
    rows = sc.objdict(
        cost_assay  = 'RTRI assay ($)',
        n_pn        = 'Indexes receiving enhanced PN',
        cost_pn     = 'Enhanced PN ($)',
        cost_total  = 'Total cost ($)',
        n_contacts  = 'Extra contacts reached',
        q_recent    = 'Truly recent among indexes',
        pos         = 'Weighted contact positivity',
        n_hiv       = 'HIV+ contacts identified',
        n_dx        = 'Extra new diagnoses',
        cost_per_dx = 'Cost per new diagnosis ($)',
    )
    cols = sc.objdict(targeted='With RTRI', budget_matched='Without RTRI (budget-matched)', full='Without RTRI (full coverage)')
    df = pd.DataFrame({label: [float(out[k][r]) for r in rows] for k, label in cols.items()}, index=list(rows.values()))
    return df


def load_assays(fn=ASSAY_FILE):
    """
    Published assay characteristics, converted to model sensitivity/specificity.
    Sensitivity for infection < 1 yr is approximated as MDRI/365 (the expected share of a
    uniformly distributed under-1-year cohort classified recent) unless reported directly.
    The assay developers reject 'sensitivity' as an assay property; this is only a bridge
    into the decision model. Uganda FRR is 0/200 observed; 0.5% is used as a conservative value.
    """
    df = pd.read_csv(fn)
    df['sens'] = np.where(df.sens_reported.notna(), df.sens_reported, df.mdri_days / 365)
    df['spec'] = 1 - df.frr
    return df


def assay_scenarios(pars=None, fn=ASSAY_FILE):
    """ Re-run the model at every published assay parameter set, plus the original assumption """
    pars = default_pars(**sc.mergedicts(pars))
    assays = load_assays(fn)
    rows = [dict(label='assumption', sens=pars.sens, spec=pars.spec)]
    rows += assays[['label', 'sens', 'spec']].to_dict('records')
    rows.append(dict(label='asante_cdc_linear + 15% undisclosed ART', sens=160/365, spec=0.979, p_prior_art=0.15))

    records = []
    for row in rows:
        kw = {k: v for k, v in row.items() if k != 'label'}
        o = evaluate(pars, **kw)
        records.append(dict(
            label      = row['label'],
            sens       = o.pars.sens,
            spec_eff   = float(o.test.spec_eff),
            flag_rate  = float(o.test.flag_rate),
            ppv        = float(o.test.ppv),
            p_rec_neg  = float(o.test.p_rec_neg),
            enrichment = float(o.enrichment),
            breakeven_rtri_cost = float(o.breakeven_rtri_cost),
            required_pn_cost    = float(o.required_pn_cost),
            cost_per_dx_targeted  = float(o.targeted.cost_per_dx),
            cost_per_dx_untargeted = float(o.budget_matched.cost_per_dx),
            penalty    = float(o.penalty),
        ))
    return pd.DataFrame(records).set_index('label')


def sweep(pars=None, x='sens', xvals=None, y='cost_rtri', yvals=None, outcome='required_pn_cost'):
    """
    Two-way sweep of any outcome attribute of evaluate() (dotted paths allowed, e.g. 'targeted.cost_per_dx').
    Returns a DataFrame indexed by y, columns x.
    """
    xvals = np.array([0.10, 0.25, 0.40, 0.55, 0.70, 0.85, 0.95]) if xvals is None else np.asarray(xvals)
    yvals = np.array([1, 2, 5, 10, 15, 20, 25]) if yvals is None else np.asarray(yvals)
    X, Y = np.meshgrid(xvals, yvals)
    o = evaluate(pars, **{x: X, y: Y})
    Z = o
    for key in outcome.split('.'):
        Z = Z[key]
    return pd.DataFrame(Z, index=pd.Index(yvals, name=y), columns=pd.Index(xvals, name=x))


# %% Validation against the spreadsheet

def validate(verbose=True):
    """ Check the defaults reproduce recency_testing_PN_model.xlsx """
    o = evaluate()
    checks = [
        ('flagged',                o.test.flagged,               61),
        ('ppv',                    o.test.ppv,                   0.409836),
        ('p_rec_neg',              o.test.p_rec_neg,             0.079872),
        ('lr_pos',                 o.test.lr_pos,                6.25),
        ('lr_neg',                 o.test.lr_neg,                0.78125),
        ('targeted.cost_total',    o.targeted.cost_total,        11220),
        ('budget_matched.n_pn',    o.budget_matched.n_pn,        561),
        ('targeted.n_dx',          o.targeted.n_dx,              2.418),
        ('budget_matched.n_dx',    o.budget_matched.n_dx,        16.24095),
        ('full.n_dx',              o.full.n_dx,                  28.95),
        ('targeted.cost_per_dx',   o.targeted.cost_per_dx,       4640.1985),
        ('budget_matched.cost_per_dx', o.budget_matched.cost_per_dx, 690.8463),
        ('penalty',                o.penalty,                    6.716687),
        ('enrichment',             o.enrichment,                 1.369235),
        ('breakeven_rtri_cost',    o.breakeven_rtri_cost,        0.450466),
        ('required_pn_cost',       o.required_pn_cost,           443.9844),
        ('window',                 o.window,                     0.393221),
    ]
    for name, got, want in checks:
        assert np.isclose(got, want, rtol=1e-5), f'{name}: got {got}, expected {want}'
    if verbose:
        print(f'✓ All {len(checks)} checks match the spreadsheet')
    return True


# %% Run as a script

def summarize(out):
    """ Print headline numbers """
    t, p = out.test, out.pars
    sc.heading('Recency test performance')
    print(f'Flagged: {t.flagged:.0f} of {p.n_cohort} ({t.flag_rate:.1%}); PPV {t.ppv:.1%}; '
          f'P(recent | negative) {t.p_rec_neg:.1%} vs prior {p.p_recent:.0%}; LR+ {t.lr_pos:.2f}, LR- {t.lr_neg:.2f}')
    sc.heading('Strategy comparison')
    print(comparison_table(out).round(2).to_string())
    sc.heading('Breakeven')
    print(f'Penalty: {out.penalty:.2f}x  (reach lost {out.reach_lost:.1f}x ÷ enrichment {out.enrichment:.2f}x)')
    print(f'Breakeven RTRI cost: ${out.breakeven_rtri_cost:.2f} (actual ${p.cost_rtri:.0f})')
    print(f'Required PN cost to justify RTRI: ${out.required_pn_cost:.0f} (actual ${p.cost_pn:.0f})')
    print(f'Test-negative group is {out.neg_vs_pos:.0%} as productive as test-positive; decision window {out.window:.0%}')
    return


if __name__ == '__main__':

    validate()
    out = evaluate()
    summarize(out)

    scens = assay_scenarios()
    sc.heading('Published assay parameter sets')
    print(scens.round(3).to_string())

    comparison_table(out).to_csv('results/recency_fermi_strategies.csv')
    scens.to_csv('results/recency_fermi_assay_scenarios.csv')
    sweep().to_csv('results/recency_fermi_sweep_required_pn_cost.csv')
    print('\nWrote results/recency_fermi_*.csv')
