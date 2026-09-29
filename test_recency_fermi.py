"""
Tests for recency_fermi.py. Run with: pytest test_recency_fermi.py
"""
import numpy as np
import pytest
import recency_fermi as rf


def test_reproduces_spreadsheet():
    assert rf.validate(verbose=False)


def test_budget_matching():
    o = rf.evaluate()
    assert np.isclose(o.targeted.cost_total, o.budget_matched.cost_total)


def test_nuisance_parameters_cancel():
    """ Contacts per PN and the undiagnosed share scale every strategy equally """
    base = rf.evaluate().penalty
    for kw in [dict(contacts=2.5), dict(p_undx=0.5), dict(n_cohort=37)]:
        assert np.isclose(rf.evaluate(**kw).penalty, base)


def test_closed_form_breakeven():
    """ At the breakeven RTRI cost, targeted and budget-matched cost per diagnosis are equal """
    be = rf.evaluate().breakeven_rtri_cost
    assert np.isclose(rf.evaluate(cost_rtri=be).penalty, 1.0)
    assert rf.evaluate(cost_rtri=0.9 * be).targeting_wins
    assert not rf.evaluate(cost_rtri=1.1 * be).targeting_wins


def test_required_pn_cost_consistent():
    o = rf.evaluate()
    assert np.isclose(rf.evaluate(cost_pn=o.required_pn_cost).penalty, 1.0)


def test_perfect_test():
    t = rf.evaluate(sens=1.0, spec=1.0).test
    assert np.isclose(t.ppv, 1.0) and np.isclose(t.p_rec_neg, 0.0)


def test_undisclosed_art_lowers_spec():
    assert rf.evaluate(p_prior_art=0.15).test.spec_eff < rf.evaluate().test.spec_eff


def test_sweep_matches_pointwise():
    df = rf.sweep(xvals=[0.25, 0.5], yvals=[10, 20])
    assert np.isclose(df.loc[10, 0.25], rf.evaluate().required_pn_cost)
    assert np.isclose(df.loc[20, 0.5], rf.evaluate(sens=0.5, cost_rtri=20).required_pn_cost)


def test_unknown_par_raises():
    with pytest.raises(KeyError):
        rf.evaluate(sensitivity=0.3)


def test_assay_table_loads():
    df = rf.load_assays()
    assert df.sens.between(0, 1).all() and df.spec.between(0, 1).all()
