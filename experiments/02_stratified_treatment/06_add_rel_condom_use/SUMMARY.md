# Experiment 02.06 — add rel_condom_use (13-par calibration; stisim 1.7.1 rc1.7.1)

**Date run:** *TBD*.
**Commit:** *TBD*.
**stisim commit:** `2d30b96` on rc1.7.1 (local branch; not pushed).
**Trials / workers:** 1000 / 50.

## Question

Exp 02.05 pinned `hiv.eff_condom` at 0.88 against its 0.9 upper — the
sampler wants more condom cooling than the biology-consistent efficacy
range allows. `rel_condom_use` (new stisim par, scalar multiplier on
condom_data) is an orthogonal lever: efficacy per condom act vs
coverage/uptake. Does opening it pull the prev/incidence overshoot
further toward data, and where does the posterior land?

## `calib_pars`

Same 12 pars as exp 02.05 plus `structuredsexual.rel_condom_use` in
[0.5, 1.5], guess 1.0.

## Result

*TBD after run.*

## Fit at 2023 (median, 10-90%)

*TBD.*

## Posterior parameters (mean, 5-95%)

*TBD.*

## Observations

*TBD. Focus:*
- Does `rel_condom_use` pin the upper (`>1.0` — sampler wants more
  coverage than data), or stay near 1.0?
- Does `eff_condom` relax off the upper now that coverage can scale?
- Does the prev overshoot cool toward 8.9% target?
- Does the extinction tail widen further, or stay similar?

## Next-step candidate

*TBD.*

## Figures

*TBD.*
