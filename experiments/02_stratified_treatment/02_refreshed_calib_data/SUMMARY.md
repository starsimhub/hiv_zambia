# Experiment 02.02 — re-run of exp 02.01 pars on refreshed calibration data

**Date run:** 2026-09-28 (rerun after data fix at `ec3782d`; earlier close
was against my old scaffold's calib CSV, not origin's refresh, because
of a `git rebase --theirs` mix-up).
**Commit:** `ec3782d`.
**Trials / workers:** 1000 / 50. **Ensemble size after shrink:** 500 draws.
**Sustainability:** 0/500 extinct at 2030 (min prev 15-49 = 4.1% at 2030).
**Mismatch:** min=25.28, mean=31.08, max=35.10.

## Question

In-epoch: same `calib_pars` as 02.01, but on `origin/main`'s refreshed
`data/zambia_hiv_calib.csv` (previous targets had a problem per
researcher — mostly a 2020-2023 trajectory adjustment plus higher
epidemic-peak deaths in early 2000s) and with age × sex extras
(`hiv.new_infections_{sex}_{ab1}_{ab2}`, `hiv.n_infected_{sex}_{ab1}_{ab2}`)
so the ZAMPHIA age × sex incidence overlay populates.

## Result

**Fit is measurably worse than exp 02.01 on the refreshed (lower)
targets — model over-transmits by 3-4× on incidence at every stratum
of the ZAMPHIA age × sex overlay.** Best mismatch 25.3 (vs 22.6 on the
old-data 02.01). `prop_m0` re-approaches its upper bound (0.77) and
`rel_death` returns to ~1.4 — the sampler is grasping for cooling. The
new age × sex incidence overlay confirms the transmission overshoot is
**uniform across strata**, not sex-asymmetric. Prevalence overshoot is
an accumulation artifact of uniformly too-high incidence, not an F/M
transmission bias.

![6-panel time series fit](figures/hiv_calibration_fit.png)

![ZAMPHIA 2016 age × sex — 4-panel: prevalence, incidence, ART coverage (data-only), VLS (data-only)](figures/zamphia_hiv_age_sex.png)

## Fit at 2023 (median, 10-90%)

| Metric | Data | Exp 02.01 (old data) | Exp 02.02 (correct data) | Change |
|---|---|---|---|---|
| Population | 20.3 M | 20.3 M | 20.3 M [18.7 – 20.7] | ✓ |
| PLHIV | 1.30 M | 1.80 M | 1.86 M [1.15 – 3.31] | slightly worse |
| New infections/yr | 23 k | 79 k | 77 k [34 – 174] | ~same |
| HIV deaths/yr | 17 k | 19 k | 22 k [14 – 33] | hotter |
| On ART | 1.27 M | 1.33 M | 1.36 M [0.86 – 2.37] | +3% |
| Prev 15-49 | 9.8 % | 14.5 % | 14.7 % [7.9 – 31.9] | ~same |

Note: 02.01 targets and 02.02 targets differ — 02.02's are the correct
refreshed data. Direct comparison of fit tables across the boundary is
approximate; refreshed data has ~lower incidence targets and higher
epidemic-peak deaths in early 2000s.

## Posterior parameters (mean, 5-95%)

| Parameter | Prior | Exp 02.01 mean | Exp 02.02 mean | Note |
|---|---|---|---|---|
| `hiv.beta_m2f` | [0.008, 0.20] | 0.0164 | 0.0168 (0.0104-0.0256) | still inside old [0.008, 0.02] |
| `hiv.rel_death` | [0.6, 1.6] | 1.236 | 1.409 (0.959-1.554) | up; back near upper |
| `structuredsexual.prop_f0` | [0.55, 0.9] | 0.766 | 0.755 (0.691-0.885) | ~same |
| `structuredsexual.prop_m0` | [0.60, 0.80] | 0.778 | 0.768 (0.638-0.799) | re-approaches upper |
| `structuredsexual.f1_conc` | [0.01, 0.2] | 0.073 | 0.079 | ~same |
| `structuredsexual.m1_conc` | [0.01, 0.2] | 0.103 | 0.077 | dropped |
| `structuredsexual.p_pair_form` | [0.4, 0.9] | 0.534 | 0.616 (0.438-0.731) | up |

## Observations

1. **Age × sex incidence overlay confirms uniform over-transmission.**
   Model boxes cover annual incidence 1-3% at F 15-24 (vs ZAMPHIA 1.07%),
   1-3% at F 25-34 (vs 1.16%), 0.3-1% at F 35-49 (vs 1.06%). Male
   incidence boxes 0-1.2% at M 15-24 (vs 0.08%), 0.4-2.6% at M 25-34
   (vs 0.25%), 0.3-2% at M 35-49 (vs 0.87%). The 3-4× over-transmission
   is roughly uniform across strata — model doesn't get F/M ratios
   badly wrong at each band, but has too much force of infection
   globally.
2. **Female mid-life prevalence overshoot has the same character.**
   Model medians 33-42% for F 25-49 vs data 20-30%. Male medians match
   data more closely. This asymmetry in the PREVALENCE overshoot,
   despite roughly-uniform incidence overshoot, is consistent with
   accumulation: higher-female incidence per ZAMPHIA + longer years-
   with-HIV = more overshoot on prevalence for females at mid-life.
3. **`rel_death` posterior back near upper prior** (1.41, 90% CI
   0.96-1.55). Deaths at 22 k vs 17 k target — mildly hot. Sampler is
   pulling `rel_death` UP anyway, which suggests it's coupled with
   other pars (probably `p_pair_form` up → more transmission → more
   PLHIV → need higher `rel_death` to hit death count in absolute
   terms with a smaller ratio).
4. **ART coverage and VLS panels are data-only for now.** stisim's
   `sti.HIV` reports only aggregate `p_on_art` and doesn't emit age × sex
   n_on_art or VLS per stratum. Model overlay would require adding an
   analyzer or extending `sti.HIV`'s stratified results. Both panels
   show the ZAMPHIA data — F ART coverage 27-72% climbing with age;
   male 15-19 suppressed; VLS 77-92% conditional on ART.
5. **Widened `beta_m2f` still bought nothing** (three experiments now).
   Recommend narrowing back to [0.008, 0.025].

## Next-step candidates

Ordered by expected leverage:

1. **Add age × sex ART coverage + VLS analyzer to the model** so the
   two bottom panels of the ZAMPHIA plot get model overlays. Small
   downstream work; doesn't affect calibration itself.
2. **Cool the aggregate force of infection.** The 3-4× uniform
   incidence overshoot is the residual problem. Candidates:
   `eff_condom` (currently 0.5 — is it too low?), the diagnosis-limited
   ART cascade capping at 72% (holding back cooling), condom-use
   scale-up curve (already realistic; probably not).
3. **Narrow `beta_m2f` back to [0.008, 0.025]** — cleanup, in-epoch.
