# Epoch 02 — stratified ART & VLS, projected as proportion

**Status:** open (2026-09-24).

## What's distinct from epoch 01

Four structural changes, all landed together — three on the treatment
cascade, one on condom use:

1. **ART coverage stratified by age × sex.** Replaces the aggregate
   `n_art` time series with a long-form `(Year, AgeBin, Gender, p_art)`
   DataFrame. Age × sex ratios come from ZAMPHIA 2016 `% On ART` by
   5-year band × sex; aggregate proportion is derived per year from
   `n_art.csv` (n_art / PLHIV from `zambia_hiv_data.csv`) and scaled by
   the stratum ratios. See `art_data.build_art_coverage`.
2. **VLS stratified by age × sex.** Passes ZAMPHIA 2016 conditional-on-
   ART viral load suppression by 3-band × sex (`data/hiv_vls_conditional_zamphia_2016.csv`)
   as `sti.ART(vls_coverage=…)`. Prior epochs used the module default
   (all ART users virally suppressed). See `art_data.build_vls_coverage`.
3. **ART held as proportion post-2023 (not count).** Prior epoch's dual-
   column trick (`p_art = 0.97` for `year >= 2024`) was a no-op because
   `n_art.csv` had no post-2023 rows. New behavior: uniform 0.95 across
   all strata for 2024-2030 (matches UNAIDS 95-95-95 country target).
   ART count now grows with PLHIV over the projection horizon rather
   than staying flat at 1.27 M.
4. **Realistic condom use scale-up.** Replaced `data/condom_use.csv`
   with substantially higher values across all partnerships: general
   partnerships ramp from 0 pre-1990 to 0.7-0.95 by 2015-2020 (was
   flat 0.01-0.1); stable-stable partnerships reach 0.9-0.95; fsw-client
   reaches 0.95 by 2015 (was 0.5). Prior values had condom use flat and
   low, understating a major transmission-cooling lever.

## Comparability caveats vs epoch 01

- Aggregate ART count trajectory differs: epoch 01 hit 1.27 M in 2023
  then flat. Epoch 02 hits 1.27 M in 2023 (matches data) then grows
  to ~1.60 M by 2030 (diagnosis-limited toward 95% coverage of
  growing PLHIV).
- VLS was implicitly ~100% in epoch 01; now age × sex differentiated
  around 85-92%, which slightly relaxes onward-transmission suppression.
- Any epoch-01 posterior with `hiv.rel_death` interpreted against
  epoch-02 mortality dynamics is nonsense — different ART cascade.

## Known open issues at start of epoch

- ART aggregate coverage plateaus around 72% in the model by 2030
  despite the 95% target — sti.ART is diagnosis-limited (can only
  initiate agents already diagnosed). Whether this is the right
  epidemiologic behavior for Zambia is a scenario-relevant question,
  not a calibration blocker.
- Data inconsistency: `n_art.csv[2023] / zambia_hiv_data.csv[2023]` =
  98% aggregate p_art. Either UNAIDS n_art is optimistic or whole-pop
  PLHIV in the calib CSV is underestimated. Documented, not resolved.
- Dedup anchor (`base_year=1990`) still contaminated. Same as epoch 01.

## Experiments

- `01_baseline/` — epoch 01 exp 03's calib_pars (with `hiv.beta_m2f`
  upper widened 10×) re-run on the new structural setup. First
  experiment on data to hit HIV deaths (19 k median vs 17 k data);
  PLHIV/prev overshot because the widened beta let Optuna find hot
  regimes. `prop_m0` still pegs; female mid-life overshoot persists.
- `02_refreshed_calib_data/` — same 7-par `calib_pars` as 02.01, but
  `data/zambia_hiv_calib.csv` was refreshed on origin between the two
  runs (previous targets had a problem, per researcher). Also adds
  age × sex incidence + PLHIV extras so the ZAMPHIA incidence plot's
  model overlay populates. Direct comparability with 02.01 on the fit
  table is limited. Uniform ~3-4× incidence overshoot across all
  strata identified.
- `03_stratified_art_vls_analyzer/` — same `calib_pars` as 02.02; adds
  a downstream `HIVArtVlsStrat` analyzer so the ZAMPHIA age × sex ART
  coverage and VLS panels get model overlays. Analyzer is
  diagnostic-only; doesn't affect fit dynamics.
