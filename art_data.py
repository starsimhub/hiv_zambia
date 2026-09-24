"""
Build age/sex-stratified ART coverage + VLS DataFrames for sti.ART().

Combines:
- Aggregate ART count time series (data/n_art.csv)
- Aggregate PLHIV time series (data/zambia_hiv_data.csv, hiv_n_infected column)
- ZAMPHIA 2016 %ART by 5-year age band x sex
  (data/hiv_treatment_status_females_age_rows.csv,
  data/hiv_treatment_status_males_age_rows.csv)
- ZAMPHIA 2016 VLS conditional on ART by 3-band x sex
  (data/hiv_vls_conditional_zamphia_2016.csv)

Method (ART coverage): aggregate proportion `p_art_agg(year) = n_art / PLHIV`
scaled by stratum-specific ratio `pct_art_stratum / pct_art_zamphia_aggregate`.
Post-2023 held at `p_art_projected` (UNAIDS 95-95-95 default 0.95).
"""

import pandas as pd


def _load_pct_art_by_stratum():
    """Load ZAMPHIA 2016 %ART by age x sex, long-form."""
    frames = []
    for path, sex in [
        ('data/hiv_treatment_status_females_age_rows.csv', 'f'),
        ('data/hiv_treatment_status_males_age_rows.csv', 'm'),
    ]:
        df = pd.read_csv(path)
        first_col = df.columns[0]
        df = df.rename(columns={first_col: 'AgeBin'})
        df = df[df['AgeBin'].str.match(r'^\d+-\d+$', na=False)].copy()
        # Strip parenthesized flags used for marginally-reliable estimates in ZAMPHIA
        df['On ART'] = df['On ART'].astype(str).str.replace(r'[()]', '', regex=True)
        df['Number'] = df['Number'].astype(str).str.replace(r'[()]', '', regex=True)
        df['pct_art'] = pd.to_numeric(df['On ART'], errors='coerce') / 100
        df['n_plhiv'] = pd.to_numeric(df['Number'], errors='coerce')
        df['Gender'] = sex
        frames.append(df[['AgeBin', 'Gender', 'pct_art', 'n_plhiv']].dropna())
    return pd.concat(frames, ignore_index=True)


def build_art_coverage(sim_start=1985, sim_end=2030, p_art_projected=0.95):
    """Stratified p_art coverage DataFrame for sti.ART(coverage=...)."""
    strat = _load_pct_art_by_stratum()

    n_art = pd.read_csv('data/n_art.csv').set_index('year')['n_art']
    plhiv = pd.read_csv('data/zambia_hiv_data.csv').set_index('year')['hiv_n_infected']
    yrs = sorted(set(n_art.index) & set(plhiv.index))
    agg_p_art = pd.Series({y: (n_art[y] / plhiv[y]) if plhiv[y] > 0 else 0.0 for y in yrs})

    # ZAMPHIA aggregate %ART, weighted by stratum PLHIV n
    agg_zamphia = (strat['pct_art'] * strat['n_plhiv']).sum() / strat['n_plhiv'].sum()
    strat['ratio'] = strat['pct_art'] / agg_zamphia

    rows = []
    data_end = int(agg_p_art.index.max())
    data_start = int(agg_p_art.index.min())
    for year in range(sim_start, sim_end + 1):
        if year > data_end:
            # UNAIDS 95-95-95 aspiration: uniform aggregate proportion across strata
            for _, r in strat.iterrows():
                rows.append({'Year': year, 'AgeBin': r['AgeBin'], 'Gender': r['Gender'],
                             'p_art': p_art_projected})
        else:
            agg_year = agg_p_art[year] if year >= data_start else 0.0
            for _, r in strat.iterrows():
                rows.append({'Year': year, 'AgeBin': r['AgeBin'], 'Gender': r['Gender'],
                             'p_art': min(r['ratio'] * agg_year, 1.0)})
    return pd.DataFrame(rows)


def build_vls_coverage():
    """Stratified VLS (conditional on ART) DataFrame for sti.ART(vls_coverage=...)."""
    return pd.read_csv('data/hiv_vls_conditional_zamphia_2016.csv')
