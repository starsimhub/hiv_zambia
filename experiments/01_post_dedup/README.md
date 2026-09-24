# Epoch 01 — post-dedup

## What's distinct

- All-cause deaths in `data/zambia_deaths.csv` passed through
  `stisim.data.dedup_deaths(base_year=1990, end_year=2030)`. The raw file
  is preserved as `data/zambia_deaths_all_cause.csv`.
- `sti.HIV` initialised with `age_bins=[0, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 100]`
  (ZAMPHIA-aligned).
- Calibration harness uses `sti.default_build_fn` with dot-notation
  parameter routing (`hiv.beta_m2f`, `structuredsexual.prop_f0`, …).
- Calibration target CSV `data/zambia_hiv_calib.csv` uses dot notation for
  column names (`hiv.prevalence_15_49`, `hiv.n_infected`, …).

## Comparability caveats vs prior state

The pre-dedup and post-dedup runs are not comparable — the HIV death
observation is on a different scale and shape after dedup. Any pre-dedup
run should be treated as historical and not benchmarked against numbers
in this epoch.

## Known open issues at start of epoch

- Dedup anchor is `base_year=1990`, but Zambia HIV prevalence is already
  ~9% by 1990. Peak AIDS-share therefore lands at ~59% vs the expected
  70–85%. A pre-1990 UN WPP row would fix this but has not been sourced.
  Deferred; may need a new epoch when addressed.

## Experiments

- `01_baseline_2026-09-23/` — 6-parameter baseline. First end-to-end
  Optuna run of the epoch. Fit not acceptable.
- `02_open_rel_death_widen_prop_m0/` — adds `hiv.rel_death` and widens
  `structuredsexual.prop_m0` upper bound to 0.98.
