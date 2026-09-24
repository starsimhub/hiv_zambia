# Experiment 02 — open `hiv.rel_death`, widen `prop_m0` upper bound

**Commit:** *TBD after run*
**Date run:** *TBD*
**Trials / workers:** 1000 / 50
**Sustainability:** *TBD*

## Diff vs experiment 01

- Added `hiv.rel_death` to search space: `low=0.6, high=1.6, guess=1.0`.
- Widened `structuredsexual.prop_m0` upper bound: `0.9` → `0.98`.
- Added `hiv.prevalence_{sex}_{ab1}_{ab2}` for all 5-year bins from 15
  upward to `extra_results` so the ZAMPHIA plot's model overlay populates.
  Diagnostic only — not a calibration target.

## `calib_pars`

```python
{
    'hiv.beta_m2f':                 dict(low=0.008, high=0.02, guess=0.012),
    'hiv.rel_death':                dict(low=0.6,   high=1.6,  guess=1.0),   # NEW
    'structuredsexual.prop_f0':     dict(low=0.55,  high=0.9,  guess=0.85),
    'structuredsexual.prop_m0':     dict(low=0.50,  high=0.98, guess=0.81),  # widened
    'structuredsexual.f1_conc':     dict(low=0.01,  high=0.2,  guess=0.01),
    'structuredsexual.m1_conc':     dict(low=0.01,  high=0.2,  guess=0.01),
    'structuredsexual.p_pair_form': dict(low=0.4,   high=0.9,  guess=0.5),
}
```

## Posterior parameters (mean, 5–95%)

*TBD after run.*

## Fit at 2023

*TBD after run.*

## What we learned

*TBD after run. Questions this experiment should answer:*

- Does opening `rel_death` bring the HIV death rate up to ~17k/yr?
- Does the widened `prop_m0` prior express as an even higher upper draw,
  or does the posterior land inside [0.5, 0.98]? If it still pegs at
  0.98, the pegging is not a prior-width problem.
- Does bringing mortality up reduce PLHIV overshoot and thereby cool
  incidence toward the 23k/yr target?

## Next-step candidate

*TBD after run.*

## Figures

*TBD after run — will embed:*

```markdown
![6-panel time series fit](figures/hiv_calibration_fit.png)
![ZAMPHIA 2016 vs model (age × sex)](figures/zamphia_hiv_prevalence_age_sex.png)
```
