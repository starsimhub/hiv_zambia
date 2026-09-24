# Analysis plan — hiv_zambia

Living brief for the HIV recency + partner-notification analysis in Zambia. Fields carry an explicit status: **Decided** / **Provisional** / **Needs evidence** / **Needs researcher decision**. Update in place as the project evolves.

Origin: `docs/Abstract_recency_NY1.docx` (Yamamoto, Stuart, Platais, Bershteyn). Scope has since been reframed away from multi-hop tracing toward a minimum-meaningful-benefit narrative adapted from `sti_notification`.

---

## Scientific frame

**Provisional research question** *(Provisional — supersedes the abstract's multi-hop framing)*.
In a generalized HIV epidemic like Zambia, given that partner notification is operationally hard and can cause harm (notification burden, false-recent notifications, opportunity cost from missed true-recents), what is the minimum meaningful benefit of PN, which partner types deliver the highest yield, and does recency-assay triggering improve the efficacy-to-harm ratio over PN triggered on all new diagnoses?

**Decision / scientific purpose** *(Provisional)*.
Inform whether and how recency-guided partner notification should feature in a national HIV response, and identify the minimum operational commitment that produces meaningful population impact.

**Setting** *(Decided)*.
Zambia.

**Population** *(Decided)*.
Whole simulated population, ages 0–100 (analysis focused on sexually-active adults). FSW segment resolved via `sti.StructuredSexual`.

**Intervention** *(Decided in principle, details deferred to Phase 2)*.
1. A new `RecencyTest` intervention: 100 % uptake among just-diagnosed (eligibility `ti_diagnosed == sim.ti`); Bernoulli(0.25) flags recent for true infections < 1 year old; Bernoulli(0.04) flags recent for older infections (false-recent).
2. A triggering-source variant on the existing single-hop `PartnerNotification`: switches between "all new diagnoses" and "recency-flagged only" as the index-case pool.

**Comparator** *(Provisional)*.
Baseline standard-of-care (existing testing / ART / PrEP with no PN) vs. PN-on-all-new-diagnoses vs. PN-on-recency-flagged-only, across current-only / prior-only / current+prior partner-network axes.

**Outcomes** *(Provisional)*.
- **Primary:** HIV infections averted 2026–2050.
- **Burden:** contacts notified per infection averted; false-recent notifications per infection averted; count of true-recent index cases missed by the recency test per year.
- **Secondary:** HIV incidence trajectory; year first reaching annual incidence < 0.1 per 100 person-years (per abstract, retained but no longer headline).

**Time horizon** *(Decided)*.
1985 model start (calibration); scenarios 2026–2050.

**Scenarios** *(Provisional — factorial not yet fully specified)*.
Axis A (PN triggering): none / all-new-diagnoses / recency-flagged-only.
Axis B (partner-type reach): current only / prior only / current + prior.
Axis C (relationship-type slicing: stable / casual / commercial): *Needs researcher decision* — add only if axes A × B motivate it.

---

## Methodological assessment

*(From `stisim:analysis-selector`, 2026-09-23.)*

- **Analytical objective:** dynamic-mechanistic estimation of population-level HIV impact of a partner-notification intervention under alternative triggering rules.
- **Dynamic transmission required:** yes — the outcomes depend on infections averted through interrupted onward transmission chains, not just on directly-tested contacts.
- **Individual-level representation potentially valuable:** yes — per-individual `ti_infected` powers the recency flag, per-individual partnership history powers PN reach, and combinatorial targeting (recency-flagged × partner-type) has no compartmental analogue at useful resolution.
- **HIVsim / STIsim suitability:** plausibly strong. HIV disease, structured sexual network with FSW segmentation, prior-partners network with configurable recall, HIV testing / ART / PrEP, and existing single-hop `PartnerNotification` are all in place. Recency intervention is a new class to add downstream.
- **Evidence gaps:** ZAMPHIA age × sex prevalence tables need extraction from the report PDF into a machine-readable table for calibration wiring.
- **Multi-method plan:** not required — single-method dynamic simulation across calibration and scenario phases.

---

## Evidence

**Known evidence / data.**
- National HIV surveillance 1990–2023 in `data/zambia_hiv_calib.csv` (whole-pop `n_alive`, `hiv.prevalence_15_49`, `hiv.n_infected`, `hiv.new_infections`, `hiv.new_deaths`) — the calibration target set.
- ZAMPHIA 2016 age × sex HIV prevalence in `data/zamphia_2016_hiv_by_age_sex.csv` — extracted from Table C.2 of the final report (weighted estimate + SE + 95% CI + unweighted N per 5-year band, sexes and totals). Held out as validation.
- STIsim capability audit (2026-09-23): `sti.HIV.ti_infected` is a FloatArr, set at infection, cleared only on death — recency arithmetic `(sim.ti - ti_infected) < ss.months(12)` works cleanly. `sti.PriorPartners` at [stisim/networks/layered_networks.py:95](/home/robyn/stisim/stisim/networks/layered_networks.py#L95) with configurable `dur_recall` (default `ss.years(1)`); construction verified working with `dur_recall=ss.years(0.25)`.

**Evidence still needed.**
- (Phase 2) Literature reference to justify the 25 % sensitivity / 96 % specificity operational assumptions for the recency assay.

**Questions blocked on research.**
- None currently blocking.

**Important uncertainties.**
- Recency test sensitivity on truly-recent infections (~25 %) is a headline assumption. False-recent rate (~4 %) is also a headline assumption. Both should be sensitivity-checked in scenarios.
- PrEP eligibility default in current `sti.Prep` is FSW-only, which differs from the abstract's implicit "population PrEP" — closer to reality but worth flagging in the writeup.

---

## Constraints

**Deliverables** *(Provisional)*.
- Calibrated 500-draw ensemble by end-of-week 2026-09-27.
- Recency + PN scenario factorial results in the following two weeks.
- Manuscript-style writeup adapting the `sti_notification` PN narrative to HIV.

**Timeline** *(Provisional beyond this week)*.
- Week of 2026-09-23: calibration modernization + refresh (Phase 1).
- Subsequent weeks: intervention build (Phase 2), scenarios (Phase 3), writeup.
- No externally-imposed hard deadline recorded.

**Stakeholders** *(Provisional)*.
Co-authors: Nao Yamamoto, Ingrida Platais, Anna Bershteyn (NYU). Consume figures + writeups; not modifying the repo directly at this stage.

**Compute environment** *(Decided)*.
Local development on Robyn's laptop. IDM Azure VMs (see `calib:idm-azure`) available if the 500-draw calibration or the scenario factorial exceeds laptop-scale.

**Collaborators** *(Provisional — see stakeholders)*.
Solo development this week; co-authors read-only for now.

---

## Working style

**Researcher STIsim experience** *(Decided)*.
Robyn is the primary developer of the STIsim package.

**Desired guidance level** *(Provisional — implicit)*.
Don't explain STIsim internals; surface only consequential design decisions; make conventional stisim/starsim choices without asking. Any modification to the editable `stisim` package itself gets flagged in advance per `stisim:extending-stisim`.

**Memory strategy** *(Decided)*.
- Curated project memory: `CLAUDE.md` (project instructions) + `ANALYSIS_PLAN.md` (this file, living) in the repo, version-controlled.
- Agent memory: `/home/robyn/.claude/projects/-home-robyn-hiv-zambia/memory/` for user preferences and cross-session recall.
- Calibration provenance: `experiments/` (introduced 2026-09-24 when the first return-loop iteration was scoped). Two-level epoch/experiment scheme documented in `experiments/README.md`.

**Repository / git strategy** *(Decided for this week)*.
Work on `main` directly. Commit frequently. Feature branch + PR workflow revisited if co-authors start touching code.

---

## Open items

**Open design decisions.**
- Recency intervention: sensitivity operationalization — flat Bernoulli(0.25) on `age_at_diagnosis < 1yr` vs. LAg-like MDRI ramp. *Recommendation: flat.*
- Partner-type stratification: axis A × B first; axis C (relationship type) only if motivated. *Recommendation: defer C.*
- Harm-metric list beyond notifications-per-infection-averted + false-recents + missed-true-recents. *Recommendation: this set is enough for the "is PN worth it" narrative; add DALYs only if the manuscript needs a QALY headline.*

**Modernization decisions taken (2026-09-23).**
- `sti.Prep` `eff_prep` → `prep_eff` kwarg rename (API drift fix). Default PrEP eligibility is now FSW-only (per current `sti.Prep`), which is closer to Zambia programme reality than "80 % population-wide".
- Custom `make_sim_pars` builder deleted; calibration and post-calibration paths now use stisim's `sti.default_build_fn` with dot-notation parameter routing.
- Calibration data columns migrated from underscore prefix (`hiv_prevalence_15_49`) to dot notation (`hiv.prevalence_15_49`) to match the stisim `sti.Calibration` result-lookup convention.
- `sti.HIV` given ZAMPHIA-aligned `age_bins=[0, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 100]` so per-band results are produced for post-calibration validation without re-running.
- Deaths file passed through `stisim.data.dedup_deaths(base_year=1990, end_year=2030)`. Raw preserved as `data/zambia_deaths_all_cause.csv`. Peak AIDS-share landed at 59 % vs the expected 70-85 %; the 1990 anchor is contaminated (Zambia HIV prev already ~9 % by 1990) and the anchor should be replaced with a pre-1990 UN WPP row when convenient.
- Output structure: `raw_results/` gitignored for bulk `.obj` files (~1.5 MB each); `results/` tracked for minimal `.df` summary tables needed by co-authors to replicate figures.

## Calibration state

Per-experiment fits live in `experiments/NN_epoch/NN_experiment/SUMMARY.md`.
Current epoch: `experiments/01_post_dedup/`. Latest completed experiment:
`01_baseline_2026-09-23` (fit not acceptable). In flight:
`02_open_rel_death_widen_prop_m0`.

**Superseded framings.**
- **Multi-hop tracing (0 / 1 / 2 / 3 hops).** Original abstract framing. Dropped 2026-09-23 in favor of the minimum-meaningful-benefit narrative. The abstract PDF at `docs/Abstract_recency_NY1.docx` still carries the superseded framing and will need revision at manuscript stage.
- **Calibration to 5 DHS waves 2002–2024.** Original abstract framing (co-author draft). Dropped 2026-09-23: no such 5-wave DHS series with fielded HIV testing exists at that cadence for Zambia; ZAMPHIA 2016 is the primary age × sex target instead. The abstract will need revision.
