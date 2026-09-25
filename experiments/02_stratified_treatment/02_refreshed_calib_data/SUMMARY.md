# Experiment 02.02 — re-run of exp 02.01 pars on refreshed calibration data

**Date run:** *TBD*.
**Commit:** *TBD*.
**Trials / workers:** 1000 / 50.

## Question

Exp 02.01 hit HIV deaths on data for the first time in the project via
the epoch 02 structural rework. Between 02.01 and 02.02, the calibration
target CSV `data/zambia_hiv_calib.csv` was refreshed on origin (previous
targets had a problem, per user). Extras were also expanded to include
`hiv.new_infections_{sex}_{ab1}_{ab2}` and
`hiv.n_infected_{sex}_{ab1}_{ab2}` so the ZAMPHIA age × sex incidence
plot's model overlay renders. Neither change is a model structural
change — kept in-epoch, but direct fit-metric comparability with 02.01
is limited by the target refresh.

## `calib_pars`

Identical to `../01_baseline/run.py`.

## Result

*TBD after run.*

## Fit at 2023 (median, 10-90%)

*TBD after run. Focus vs 02.01:*
- Do the same posteriors hit the refreshed targets, or does the sampler
  shift?
- Does the age × sex incidence overlay reveal a comparable F/M asymmetry
  to what the prevalence plot showed?

## Posterior parameters (mean, 5-95%)

*TBD after run.*

## Observations

*TBD.*

## Next-step candidate

*TBD.*

## Figures

*TBD after run — will embed:*

```markdown
![6-panel time series fit](figures/hiv_calibration_fit.png)
![ZAMPHIA 2016 age × sex — prevalence + incidence](figures/zamphia_hiv_age_sex.png)
```
