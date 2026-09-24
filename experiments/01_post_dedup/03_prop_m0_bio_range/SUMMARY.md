# Experiment 03 — constrain `structuredsexual.prop_m0` to biologically plausible range

**Date run:** *TBD*.
**Commit:** *TBD*.
**Trials / workers:** 1000 / 50.

## Question

Exp 02 posterior mean for `structuredsexual.prop_m0` was 0.955 with the
5-95% at 0.92-0.98, still pegging the widened 0.98 upper bound. The
researcher flags this as biologically implausible — the low-risk share
of males in Zambia should sit around 60-80%, not 95%+. If we constrain
the prior to [0.60, 0.80], can the calibration find other knobs to cool
transmission, or does the fit collapse?

## Diff vs exp 02

- `structuredsexual.prop_m0`: `[0.50, 0.98]` → `[0.60, 0.80]`, guess
  0.81 → 0.70. All other parameters unchanged.

## `calib_pars`

See `run.py`.

## Result

*TBD after run.*

## Fit at 2023 (median, 10-90%)

*TBD after run.*

## Posterior parameters (mean, 5-95%)

*TBD after run.*

## Observations

*TBD. Questions this experiment should answer:*

- Where does `prop_m0` land inside [0.60, 0.80]? If it pegs the new
  upper (0.80), the biological ceiling is still tight against the fit.
- Do other parameters shift to compensate? Candidates: `beta_m2f` drops
  toward the low end of its prior, `prop_f0` shifts, `p_pair_form`
  drops.
- Does aggregate mismatch worsen? Ensemble mismatch in exp 02 was
  23.8-30.5 (tight). A meaningful degradation here signals that the
  next lever (VL coverage, condom effectiveness, differential f2m
  transmission) is genuinely needed.

## Next-step candidate

*TBD after run.*

## Figures

*TBD after run — will embed:*

```markdown
![6-panel time series fit](figures/hiv_calibration_fit.png)
![ZAMPHIA 2016 vs model (age × sex)](figures/zamphia_hiv_prevalence_age_sex.png)
```
