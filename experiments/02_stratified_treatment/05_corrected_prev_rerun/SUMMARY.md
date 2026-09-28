# Experiment 02.05 — 02.04 pars re-run on corrected prev 15-49 targets

**Date run:** *TBD*.
**Commit:** *TBD*.
**Trials / workers:** 1000 / 50.

## Question

Exp 02.04 halved mismatch (12.7 vs 24.0) with the new FOI levers, but
its fit was against the pre-refresh `hiv.prevalence_15_49` column. User
supplied a corrected UNAIDS series (commit `fdfa4e8`): peak higher in
mid-1990s / mid-2000s (15.1% at 1999 vs prior 14%), endpoints lower
(2023 = 8.9% vs prior 9.8%). Same 12 calib_pars as 02.04. Does the
02.04 posterior structure hold on the corrected targets, and does
`hiv.eff_condom` still pin upper?

## `calib_pars`

Identical to `../04_open_foi_levers/run.py`.

## Result

*TBD after run.*

## Fit at 2023 (median, 10-90%)

*TBD.*

## Posterior parameters (mean, 5-95%)

*TBD.*

## Observations

*TBD. Focus:*
- Does `eff_condom` still pin upper on corrected targets, or does the
  higher-peak prev pull it back?
- Does `beta_m2f` stay outside old [0.008, 0.02] range?
- Does the female mid-life residual persist as the isolated problem,
  or does the peak-lifting change the age × sex pattern?

## Next-step candidate

*TBD. Likely one of:*
- Add `structuredsexual.rel_condom_use` upstream in stisim and reopen
  the condom cooling question with a second lever.
- Attack the female mid-life residual with a differential parameter
  or age-mixing knob.
- Declare epoch 02 complete and move to scenarios.

## Figures

*TBD.*
