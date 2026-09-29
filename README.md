# hiv_zambia

Agent-based model of HIV transmission in Zambia, built on [STIsim](https://github.com/starsimhub/stisim). The analytical question is whether **partner notification** — operationally hard and potentially harmful — is worth doing at all in a generalized HIV epidemic, which partner types (current vs prior) deliver the highest yield, and whether recency-flagged targeting earns its cost. Companion to [`sti_notification`](https://github.com/starsimhub/sti_notification), where the partner-notification framing originates for bacterial STIs.

## Status

**Phase 1 (calibration) — closed 2026-09-28.** Baseline: `experiments/02_stratified_treatment/08_bigger_vls_condom_patterns/`. 13-parameter Optuna posterior on stisim `rc1.7.1`, 1000 trials shrunk to the top 100 by mismatch (committed at `outputs/draws_used.csv`, best mismatch 11.4). See the epoch history in [`experiments/02_stratified_treatment/README.md`](experiments/02_stratified_treatment/README.md) for the sequence of levers that got there — stratified ART + VLS coverage, time-varying viral load suppression, and realistic condom scale-up.

**Phase 2 (interventions) — up next.** Add `RecencyTest` (100% uptake among just-diagnosed, ~25% sensitivity on infections < 1 year, ~96% specificity) and extend the existing `PartnerNotification` intervention with an all-new-diagnoses vs recency-flagged-only trigger.

**Phase 3 (scenarios).** Factorial across PN triggering × partner-type reach × recency filter, with contacts-notified-per-infection-averted and false-recent-notification burden as the headline outcomes.

## Repo layout

- `hiv_model.py` — sim builder (HIV + StructuredSexual + PriorPartners + MaternalNet, Zambia demographics)
- `interventions.py` — testing arms + ART + PrEP + single-hop `PartnerNotification`
- `analyzers.py` — `HIVArtVlsStrat` age × sex diagnostic overlays
- `art_data.py` — stratified ART / VLS coverage builders from ZAMPHIA 2016 + national totals
- `run_hiv_calibration.py` — Optuna calibration driver (single canonical entry point; per-experiment `run.py` files just define `calib_pars` and call `run_and_save`)
- `run_pn_scens.py` — scenario runner (Phase 3)
- `recency_fermi.py`, `plot_recency_fermi.py`, `test_recency_fermi.py` — companion closed-form cost analysis of recency-triggered PN (see below)
- `plot_calibrations.py`, `plot_zamphia_age_sex.py` — reproducible calibration figures
- `data/` — Zambia demography, national HIV surveillance, ZAMPHIA 2016 age × sex prevalence / incidence / ART / VLS
- `experiments/` — calibration provenance in a two-level `NN_epoch/NN_experiment/` scheme where an **epoch** is a structural model change (new data preprocessing, new module) and an **experiment** is a parameter change within an epoch. Each experiment folder holds `SUMMARY.md` (commit hash + `calib_pars` + fit table + notes) and `figures/`. Raw outputs stay gitignored; the recoverability contract is commit hash + `calib_pars`.

## Environment

- stisim `1.7.1` on branch [`rc1.7.1`](https://github.com/starsimhub/stisim/tree/rc1.7.1), editable install at `/home/robyn/stisim/` (adds `structuredsexual.rel_condom_use` — a scalar multiplier on the condom-use matrix, not upstreamed as of 2026-09-28)
- starsim `3.6.1`
- Python 3.11 (conda env `starsim`)

## Reproducing the baseline

```bash
# From repo root, with the `starsim` env active:
python experiments/02_stratified_treatment/08_bigger_vls_condom_patterns/run.py
```

This runs 1000 Optuna TPE trials, shrinks to the top 500 by mismatch, and writes `raw_results/zam_hiv_calib.obj`. Then:

```bash
python plot_calibrations.py           # 6-panel time-series fit
python plot_zamphia_age_sex.py        # 2 × 2 ZAMPHIA 2016 overlay (prev + incidence + ART + VLS)
```

Runtime: ~2–3 h on 50 workers.

## Companion decision analysis

`recency_fermi.py` is a closed-form (Fermi) cost model that asks whether recency-test triggering of enhanced partner notification earns its cost, independent of the transmission model. It compares, for a cohort of new diagnoses: RTRI for everyone with enhanced PN for the RTRI-recent; the same budget spent on untargeted enhanced PN; and enhanced PN for everyone. Targeting wins only if `cost_rtri < flag_rate × cost_pn × (enrichment − 1)`.

```bash
python recency_fermi.py        # validate against the original spreadsheet, print summary, write results/recency_fermi_*.csv
python plot_recency_fermi.py   # figures/recency_fermi.png
pytest test_recency_fermi.py
```

All functions broadcast over numpy arrays, so `rf.evaluate(sens=np.linspace(0, 1, 101))` or `rf.sweep(...)` gives sensitivity analyses directly. `data/recency_assays.csv` holds published MDRI / false-recent rates for current recency assays, used as a robustness check (sensitivity approximated as MDRI/365). Contact positivity (40% / 17%) is a structural placeholder and the input that drives the result; the ABM's recency-stratified PN yield is the natural replacement once Phase 2 lands.

## Related projects

- **[sti_notification](https://github.com/starsimhub/sti_notification)** — parent project; PN "minimum meaningful benefit" framing originates there for bacterial STIs
- **[hivsim_zim](https://github.com/starsimhub/hivsim_zim)** — sibling HIV-only Zimbabwe validation model
- **[STIsim](https://github.com/starsimhub/stisim)** — underlying simulation framework
