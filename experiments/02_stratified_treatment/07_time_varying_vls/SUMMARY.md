# Experiment 02.07 — time-varying VLS trajectory + widened eff_condom upper

**Date run:** *TBD*.
**Commit:** *TBD*.
**Trials / workers:** 1000 / 50.

## Question

Exp 02.06 hit best mismatch 10.8 but new infections in 2020-2025
declined too slowly relative to UNAIDS (model 50 k vs data 26 k at
2023). The researcher's diagnosis: model's ART efficacy is time-invariant
at ZAMPHIA-2016 rates, missing the DTG-era improvement in VLS + the
poor early-cART efficacy in the 2000s. Making per-agent ART efficacy
grow with time — via a time-varying `vls_coverage` trajectory —
should cool 2020-2025 infections without needing new calibration
parameters.

Also: `hiv.eff_condom` posterior pinned at 0.88 upper in 02.06;
widening prior to 0.95 gives it further room.

## Diff vs exp 02.06

- `data/hiv_vls_conditional_zamphia_2016.csv` → `hiv_vls_conditional_over_time.csv`.
  Trajectory scaled from ZAMPHIA 2016 anchor:
  - 2000: ~50% of ZAMPHIA (early cART, poor adherence)
  - 2010: ~85% (cART established, before DTG)
  - 2016: ZAMPHIA-measured (77-92% by stratum)
  - 2020: 79-94% (DTG rollout)
  - 2025: 82-95% (target VLS)
- `hiv.eff_condom` prior upper widened 0.9 → 0.95.
- 13 calib pars, same list as 02.06.

## Result

*TBD after run.*

## Fit at 2023 (median, 10-90%)

*TBD. Focus:*
- Does 2020-2025 infection curve steepen toward the ~26 k target?
- Does prev 15-49 at 2023 drop from 13.3% toward the 8.9% target?

## Posterior parameters (mean, 5-95%)

*TBD. Focus:*
- Does `eff_condom` still pin — now at 0.95?
- Does `rel_condom_use` stay at 1.26 or drop (since VLS now doing more)?
- Does `beta_m2f` relax from 0.042 back toward biological range?

## Observations

*TBD.*

## Next-step candidate

*TBD.*

## Figures

*TBD.*
