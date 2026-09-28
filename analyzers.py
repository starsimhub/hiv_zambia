"""
Downstream analyzers for the Zambia HIV model.

`HIVArtVlsStrat` records the age × sex breakdown of PLHIV, on-ART, and
virally-suppressed agents at each timestep — quantities that stisim's
`sti.HIV` does not stratify by default. Feeds the ART coverage and VLS
panels of `plot_zamphia_age_sex.py`.
"""

import numpy as np
import starsim as ss


class HIVArtVlsStrat(ss.Analyzer):
    """Stratified ART coverage + VLS by age band × sex, on the HIV module.

    Emits per-timestep counts under `sim.results.hiv_art_vls_strat`:
    - `n_infected_{sex}_{ab1}_{ab2}` — PLHIV in stratum
    - `n_on_art_{sex}_{ab1}_{ab2}`   — on ART in stratum
    - `n_vls_{sex}_{ab1}_{ab2}`      — on effective ART (virally
      suppressed) in stratum

    Proportions (% on ART among PLHIV; % VLS among on-ART) are computed
    in the plot from these three counts to keep the analyzer
    denominator-agnostic.
    """

    def __init__(self, age_bins=None, sex_keys=('f', 'm'), **kwargs):
        super().__init__(**kwargs)
        self.age_bins = age_bins or [15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 100]
        self.sex_keys = sex_keys

    def init_results(self):
        super().init_results()
        results = []
        for sex in self.sex_keys:
            for ab1, ab2 in zip(self.age_bins[:-1], self.age_bins[1:]):
                for base in ('n_infected', 'n_on_art', 'n_vls'):
                    results.append(ss.Result(f'{base}_{sex}_{ab1}_{ab2}', dtype=int))
        self.define_results(*results)

    def step(self):
        sim = self.sim
        ppl = sim.people
        hiv = sim.diseases.hiv
        ti = self.ti

        age = ppl.age
        for sex in self.sex_keys:
            sex_mask = ppl.female if sex == 'f' else ppl.male
            for ab1, ab2 in zip(self.age_bins[:-1], self.age_bins[1:]):
                stratum = sex_mask & (age >= ab1) & (age < ab2)
                self.results[f'n_infected_{sex}_{ab1}_{ab2}'][ti] = int((stratum & hiv.infected).sum())
                self.results[f'n_on_art_{sex}_{ab1}_{ab2}'][ti] = int((stratum & hiv.on_art).sum())
                self.results[f'n_vls_{sex}_{ab1}_{ab2}'][ti] = int((stratum & hiv.on_effective_art).sum())
