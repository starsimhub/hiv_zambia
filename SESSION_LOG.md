## 2026-09-23

**Done this session:**
- Project intake: reframed the recency-testing analysis away from the abstract's multi-hop tracing question toward a minimum-meaningful-PN + role-of-recency narrative (adapting the `sti_notification` framing). Scaffolding landed: `CLAUDE.md`, `ANALYSIS_PLAN.md` (living brief with status markers), agent memory at `~/.claude/projects/-home-robyn-hiv-zambia/memory/`. ZAMPHIA 2016 held out as validation, not calibration.
- Modernization pass on the dormant Oct-2024 code: `sti.Prep(eff_prep=)` → `prep_eff=`; calibration CSV columns migrated to dot notation (`hiv.prevalence_15_49`, etc.) to match current `sti.Calibration` lookup convention; custom `make_sim_pars` deleted in favor of `sti.default_build_fn`; calibration parameter keys renamed to dot notation (`hiv.beta_m2f`, `structuredsexual.prop_f0`, …); `sti.HIV` given ZAMPHIA-aligned `age_bins` so per-band results are produced for validation without re-running.
- ZAMPHIA 2016 Table C.2 extracted from the report PDF into `data/zamphia_2016_hiv_by_age_sex.csv` (weighted prevalence + SE + 95 % CI + unweighted N per 5-year band × sex).
- First 1000-trial Optuna calibration on the modernized path ran end-to-end. Fit revealed a transmission-too-high signature. Diagnosed the underlying issue as the mortality double-count via `stisim:model-primer`'s `calibration-knobs.md` reference — UN WPP all-cause deaths include AIDS; `sti.HIV` also kills agents; double-counted → PLHIV pool inflated → prevalence and incidence stuck too high.
- Applied `stisim.data.dedup_deaths` (with `base_year=1990` — the file has no 1985 anchor) and reran. Preserved raw as `data/zambia_deaths_all_cause.csv`. Fit still not great: dedup traded a death overshoot for a death undershoot; PLHIV overshoots by ~500 k; incidence 3× too high. Enough to pause and think tomorrow.

**State of the analysis:**
- Working: model builds and runs cleanly on stisim 1.7.0 / starsim 3.6.1; calibration harness runs end-to-end in ~5 min on 50 workers; ART cascade fits well; ZAMPHIA validation data ready.
- Not yet: fit to PLHIV, incidence, prevalence, HIV deaths is off (see numbers in ANALYSIS_PLAN.md). Phase 2 (recency intervention + PN triggering variants) not started.
- Ruled out: multi-hop tracing (dropped from abstract framing); ZAMPHIA as a calibration target (held out for validation).

**Next steps:**
1. Improve the dedup anchor. Peak AIDS-share came out at 59 % (expected 70-85 %) because the 1990 anchor is already contaminated. Source pre-1990 UN WPP data (1985 or earlier) so `dedup_deaths` can interpolate from a clean pre-epidemic baseline.
2. If step 1 alone doesn't fix the ~2× HIV-death undershoot, add `sti.HIV`'s `rel_death` (and possibly `rel_death_f`, `art_death_age`) to the calibration search space — see `stisim:model-primer`'s `calibration-knobs.md` for the full mortality-knob catalog.
3. Investigate why incidence stays 3× too high after mortality correction. `structuredsexual.prop_m0` hit its upper bound (0.9) in 342/500 top draws — either widen the prior (up to ~0.98) or check whether `sti.ART`'s default `vls_coverage` is doing enough transmission reduction in the ART-scale-up era.
4. Only once the historical fit is credible: proceed to Phase 2 (`RecencyTest` intervention + `PartnerNotification` triggering-source variant per the ANALYSIS_PLAN.md scope).

**Unresolved:**
- Round-2 intake questions (recency spec details, partner-type stratification axis, harm-metric set) still deferred — properly belong to Phase 2, not blocking Phase 1.
- Skill-discovery gap: `stisim:analysis-intake` and `stisim:analysis-selector` fired at project start but did not surface the mortality-dedup requirement, which lives in `stisim:model-writer`'s and `stisim:model-primer`'s references. Memory saved to both `sti_notification` and `hiv-zambia` agent memory (`feedback_mortality_dedup_always_fires.md`) to auto-fire this check regardless of intake path.
- Whether `hivsim_zim` (sibling repo) applied the same dedup — worth checking before its calibration is reused.

**Workspace state:**
- Branch: `main`, up to date with `origin/main` before this session's commit.
- Uncommitted: 8 modified + scaffolding + extracted data + preserved raw + figures + pre-dedup snapshot. About to commit.
- Output structure: `raw_results/` gitignored (bulk `.obj` files); `results/` tracked (minimal `.df` summary tables needed for replication).
