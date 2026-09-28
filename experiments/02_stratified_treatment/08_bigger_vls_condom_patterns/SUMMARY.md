# Experiment 02.08 — bigger VLS swing + condom pattern changes

**Date run:** *TBD*.
**Commit:** *TBD*.
**Trials / workers:** 1000 / 50.

## Question

Exp 02.07's gentle VLS swing (50% → 94% across 2000-2020) shifted
posteriors but didn't cool 2020-2025 infections. This experiment tries
a bigger swing (~35% → ~95%) and extends condom_use.csv with 2025 +
raises specific partnership types. Does the aggregate FOI in 2020-2025
finally track the UNAIDS ~26 k target?

## Diff vs exp 02.07 (data only)

- `hiv_vls_conditional_over_time.csv` — trajectory scaled 0.35 (2000) /
  0.70 (2010) / 1.0 (2016) / boosted (2020: 92-96%, 2025: 95-97%).
- `condom_use.csv`:
  - Extended with 2025 column.
  - (0,0) floor raised 0.05 → 0.15-0.20 across data years.
  - Cross-risk partnerships raised 0.80 → 0.90 at 2020, 0.70 → 0.80 at 2015.

Same 13 pars as 02.07.

## Result

*TBD.*

## Fit at 2023 (median, 10-90%)

*TBD. Focus:*
- Does 2020-2023 new infections drop from 55 k (02.07) toward 26 k data?
- Does prev 15-49 at 2023 drop from 12.7% toward 8.9% data?
- Does the female mid-life prev overshoot cool at all?

## Posterior parameters (mean, 5-95%)

*TBD. Focus:*
- Does `hiv.eff_condom` relax further (VLS + more condom coverage doing more)?
- Does `rel_condom_use` come off the upper (data now higher, less need to scale)?
- Does `beta_m2f` drop further toward biological ~0.02?

## Observations

*TBD.*

## Next-step candidate

*TBD. Depending on how the fit lands — either declare epoch 02 done or
try one more lever (time-varying nonsupp_art_efficacy?).*

## Figures

*TBD.*
