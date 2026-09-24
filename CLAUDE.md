# CLAUDE.md

Agent-based HIV transmission model for Zambia. Companion analysis to `sti_notification`'s partner-notification framing, adapted from bacterial STIs to HIV: given that partner notification is operationally hard and potentially harmful, what is the minimum meaningful benefit, and does recency-flagged targeting earn its cost?

See `README.md` for orientation and `ANALYSIS_PLAN.md` for scope, intake, and current state.

## Repo layout

Root-level Python for a small sprint-scale project.

- `hiv_model.py` — model builder
- `interventions.py` — HIV testing + ART + PrEP + single-hop `PartnerNotification` (recency test to be added)
- `run_hiv_calibration.py` — Optuna calibration entry point (single canonical driver; not duplicated per experiment)
- `run_pn_scens.py` — scenario runner
- `plot_*.py` — figure scripts
- `utils.py` — helpers
- `data/` — Zambia demography + national HIV surveillance
- `docs/` — abstract, ZAMPHIA 2016 final report
- `experiments/` — calibration provenance (see `experiments/README.md`).
  Two-level: `NN_epoch/NN_experiment/` where an **epoch** is a structural
  model change (new data preprocessing, new module, etc.) and an
  **experiment** is a parameter change within an epoch. Each experiment
  folder holds `SUMMARY.md` (commit hash + calib_pars + fit table + notes)
  and `figures/` copies. Raw outputs (`raw_results/*.obj`) stay gitignored;
  recoverability contract is commit hash + calib_pars.

Results go to `results/` (gitignored bulk; committable summaries only).

## State of play

**Phase 1 (in progress).** Calibration modernization on current stisim (1.7.0 / starsim 3.6.1). Completed: `sti.Prep` kwarg fix (`eff_prep` → `prep_eff`); calibration data columns migrated to dot notation (`hiv.prevalence_15_49` etc.); custom `make_sim_pars` deleted in favor of stisim's `default_build_fn` with dot-notation calib-par routing (`hiv.beta_m2f`, `structuredsexual.prop_f0`, …); `sti.HIV` configured with ZAMPHIA-aligned 5-year `age_bins`; ZAMPHIA 2016 age×sex prevalence extracted to `data/zamphia_2016_hiv_by_age_sex.csv` (held out as validation, not a calibration target). The `run_hiv_calibration.py` script runs 1000 Optuna TPE trials and shrinks to the top 500 draws.

**Phase 2.** Add `RecencyTest` intervention: 100% uptake among just-diagnosed, sensitivity 25% on infections < 1 year, specificity 96%. Extend `PartnerNotification` with a triggering-source variant (all-new-diagnoses vs recency-flagged-only).

**Phase 3.** Scenario factorial across PN triggering × partner-type reach × recency filter, with burden analyzers (contacts notified per infection averted, false-recent notifications, missed true-recent indices).

Multi-hop tracing is explicitly out of scope — the analytical question is minimum meaningful benefit, not tracing-depth optimization.

## Intake

**Model.** `sti.HIV` + `sti.StructuredSexual` (FSW-segmented) + `sti.PriorPartners(dur_recall=3mo)` + `MaternalNet`, `demographics='zambia'`, 10k agents, 1985 start. Baseline services: FSW / general-population / CD4 < 200 HIV testing arms with historical scale-up curves; ART targeting 97 % by 2024 (from `data/n_art.csv`); PrEP scaling to 80 % among FSW by 2025.

**Question.** Is partner notification worth doing at all in a generalized HIV epidemic? Which partner types (current vs prior; possibly by relationship type) deliver the highest yield? Does recency filtering — despite ~75 % missed-true-recent and ~4 % false-recent burden — improve the efficacy-to-harm ratio?

**Data.** National HIV surveillance 1990–2023 in `data/zambia_hiv_calib.csv`. ZAMPHIA 2016 final report PDF in `docs/` for age × sex prevalence targets.

**Constraints.** Solo development this week (Robyn). IDM compute if needed. End-of-week: calibrated ensemble ready for scenarios.

## Environment

- stisim 1.7.0 (editable at `/home/robyn/stisim/`)
- starsim 3.6.1
- Python at `/home/robyn/miniconda/bin/python`; no conda env activation required

## Conventions

- Any modification to the editable `stisim` install is committed, pushed, and PR'd immediately per `stisim:editable-dep-hygiene`.
- Downstream vs upstream decisions run through `stisim:extending-stisim` — real bugs go upstream; opt-in project knobs stay downstream.
- Comment discipline in shared library code per `stisim:comment-hygiene` — no project-history in stisim source.

## Related projects

- **[sti_notification](https://github.com/starsimhub/sti_notification)** — parent project; the PN "minimum meaningful benefit" framing originates there for bacterial STIs.
- **[hivsim_zim](https://github.com/starsimhub/hivsim_zim)** — sibling HIV-only slim model for Zimbabwe validation.
