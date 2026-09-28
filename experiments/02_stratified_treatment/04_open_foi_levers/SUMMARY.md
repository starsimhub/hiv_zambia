# Experiment 02.04 — open five global-FOI levers to attack uniform ~3-4× incidence overshoot

**Date run:** *TBD*.
**Commit:** *TBD*.
**Trials / workers:** 1000 / 50.

## Question

Exps 02.01-02.03 showed model over-transmits by ~3-4× uniformly across
all F/M age × sex strata. The parameter set to date (7 pars, mostly
network structure at lower-risk tiers) cannot pull FOI down enough.
This experiment opens levers on:

- The **highest-risk fractions** (`prop_f2`, `prop_m2`) — currently
  fixed at 0.01 / 0.02. If the model has too much high-risk activity
  driving concurrency-heavy transmission, tuning these down cools.
- **High-risk concurrency** (`f2_conc`, `m2_conc`) — currently fixed at
  0.1 / 0.5. Analogous to f1/m1 already open, but in the high-risk
  tail.
- **HIV condom effectiveness** (`hiv.eff_condom`) — currently 0.5. HIV
  literature is 0.7-0.85; 0.5 likely understates condom cooling given
  the realistic scale-up we wired (0.7-0.95 coverage). Not a network
  knob but a biological calibration correction the user flagged as
  the same class of decision.

Twelve pars total.

## `calib_pars`

See `run.py`.

## Result

*TBD after run.*

## Fit at 2023 (median, 10-90%)

*TBD.*

## Posterior parameters (mean, 5-95%)

*TBD.*

## Observations

*TBD. Focus:*
- Does the aggregate FOI cool enough to bring PLHIV/prev toward data?
- Where do the new posteriors land — does `eff_condom` pin near the
  upper (biology-consistent) or lower (0.5 was fine)?
- Do the high-risk fractions want to shrink (fewer concurrent
  partnerships) or widen (already needed)?
- Does age × sex incidence overshoot remain uniform, or does it
  become more sex/age-asymmetric now that the cooling levers can
  work?

## Next-step candidate

*TBD.*

## Figures

*TBD after run.*
