# Epoch 01 — post-dedup

**Status:** closed 2026-09-24. Parameter-tweak budget exhausted; residual
misfit is structural. Next epoch opens once the missing data pieces
(VL coverage, and any sex-asymmetric transmission evidence) are in place.

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
- `02_open_rel_death_widen_prop_m0/` (2026-09-24) — adds `hiv.rel_death`
  and widens `structuredsexual.prop_m0` to 0.98. Partial improvement
  (PLHIV/prev overshoot came down); HIV deaths unchanged despite
  `rel_death` posterior near upper bound; `prop_m0` still pegs. Age × sex
  diagnostic reveals systematic female overshoot at 25-49 that's the
  likely driver of the residual aggregate misfit.
- `03_prop_m0_bio_range/` (2026-09-24) — constrained
  `structuredsexual.prop_m0` to biologically plausible [0.60, 0.80].
  Fit barely moved; `prop_m0` still pegs (now at 0.80); compensating
  parameter shifts physically ambiguous; female mid-life overshoot
  unchanged. Confirms the epoch has hit a structural fit ceiling.

## What we learned across the epoch

1. **Structural fit ceiling around mismatch ~24-31.** All three
   experiments plateau at similar mismatch despite widely different
   parameter configurations. Parameter-only tweaks within this epoch
   cannot break the ceiling.

2. **Female mid-life prevalence overshoot is the largest residual
   signal.** Model overshoots ZAMPHIA at ages 25-49 for women (~35% at
   30-39 vs ~22% data); male fit is clean across all bands. The
   current parameter set is symmetric in transmission direction and
   cannot express the F/M asymmetry the data requires. Highest-leverage
   next-epoch move: open `hiv.beta_f2m` differentially.

3. **`structuredsexual.prop_m0` pegs upper bounds robustly.** Pegged
   at 0.9 (exp 01), 0.98 (exp 02), and 0.8 (exp 03). The model uses
   this knob as its only handle on cooling male-side aggregate
   transmission. Adding an orthogonal transmission-cooling knob (VL
   coverage uncertainty, condom effectiveness) is needed to relieve
   the load.

4. **HIV deaths stuck at 6.8k/yr across all three experiments.**
   `hiv.rel_death` posterior pushes 1.3-1.4 (near upper bounds)
   without moving the death count. Death rate per PLHIV is stable
   ~0.4%. Root cause not identified — could be age-at-death mismatch
   (ART scale-up preventing deaths at ages the data expects them),
   or observation-window drop. Warrants a separate investigation
   experiment before a next-epoch structural change.

## What's blocking the next epoch

- **VL coverage data for Zambia** (researcher sourcing 2026-09-24) to
  parameterise `sti.ART.vls_coverage` with uncertainty.
- **Any evidence for sex-asymmetric HIV transmission in Zambia** to
  motivate opening `hiv.beta_f2m` differentially rather than
  symmetrically.
- **Dedup anchor upgrade** (pre-1990 UN WPP row) still deferred; may
  or may not warrant epoch separation.
