# Experiment 02.03 — add stratified ART/VLS analyzer, same calib_pars as 02.02

**Date run:** *TBD*.
**Commit:** *TBD*.
**Trials / workers:** 1000 / 50.

## Question

02.02's ZAMPHIA age × sex plot had data-only panels for ART coverage
and VLS because stisim doesn't stratify those quantities. Adding a
downstream `HIVArtVlsStrat` analyzer emits per-stratum `n_infected`,
`n_on_art`, and `n_vls` counts. Same 7-parameter `calib_pars` as 02.02.
Does the model reproduce ZAMPHIA's age × sex ART coverage and VLS
patterns given the stratified inputs we've already wired?

## `calib_pars`

Identical to `../02_refreshed_calib_data/run.py`.

## Result

*TBD after run.*

## Fit at 2023 (median, 10-90%)

*TBD.*

## Posterior parameters (mean, 5-95%)

*TBD.*

## Observations

*TBD. Focus:*
- Does model ART coverage box-plot cover the ZAMPHIA points at each
  5-year band × sex? Youngest bands (15-19) are the most likely
  mismatch — data ART coverage there is only 26.5% F, hard to reach
  via the diagnosis-then-treat pipeline.
- Does model VLS cover the ZAMPHIA points (75-92% conditional on ART)?
  Values come from the `vls_coverage` DataFrame + `p_effective_art` at
  ART start; anything off would indicate our stratified VLS wiring
  isn't landing correctly.

## Next-step candidate

*TBD. Depending on findings, likely a network-side cooling lever
investigation to attack the uniform ~3-4x FOI overshoot.*

## Figures

*TBD after run.*
