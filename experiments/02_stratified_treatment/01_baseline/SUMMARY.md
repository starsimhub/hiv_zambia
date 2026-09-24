# Experiment 02.01 — baseline: epoch 01 exp 03 pars on stratified-treatment setup

**Date run:** *TBD*.
**Commit:** *TBD*.
**Trials / workers:** 1000 / 50.

## Question

Epoch 01 closed with a structural fit ceiling: female mid-life
prevalence overshoot, HIV deaths stuck ~2× low, `structuredsexual.prop_m0`
pegging its upper prior. Epoch 02 makes three model-side changes without
touching priors ([see epoch README](../README.md)): stratified ART
coverage, stratified VLS, and ART held as proportion post-2023. This
experiment isolates the effect of that structural rework on the same
7-parameter search — does the treatment-cascade wiring alone move the
fit?

## `calib_pars`

From `experiments/01_post_dedup/03_prop_m0_bio_range/run.py`, with
`hiv.beta_m2f` upper widened 10x (0.02 → 0.2). Reasoning: the four
structural changes cool transmission substantially at the guess pars,
so beta may want a wider search range now. See `run.py`.

## Result

*TBD after run.*

## Fit at 2023 (median, 10-90%)

*TBD. Comparison vs epoch 01 exp 03 will focus on:*

- **HIV deaths**: does the improved mortality dynamics from stratified
  ART change the death count meaningfully?
- **Prev 15-49 / PLHIV**: does the growing (rather than flat) ART pool
  cool onward transmission enough to bring aggregate prev toward data?
- **Female age × sex overshoot**: does the female-skewed ART share
  (~70% of ART users) selectively cool female prevalence?
- **`hiv.rel_death` posterior**: still push toward the upper bound, or
  does the mortality problem move?
- **`prop_m0`**: still peg 0.80, or does it relax?

## Posterior parameters (mean, 5-95%)

*TBD after run.*

## Observations

*TBD.*

## Next-step candidate

*TBD. Depends on what changes and what doesn't.*

## Figures

*TBD after run — will embed:*

```markdown
![6-panel time series fit](figures/hiv_calibration_fit.png)
![ZAMPHIA 2016 vs model (age × sex)](figures/zamphia_hiv_prevalence_age_sex.png)
```
