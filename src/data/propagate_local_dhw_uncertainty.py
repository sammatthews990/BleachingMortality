'''Propagate local DHW calibration uncertainty into operational predictions.

The deterministic fully local-calibrated INLA prediction remains the central
estimate. Correction uncertainty is propagated with a local response slope
estimated from the current-event-only INLA contrast, then combined with
programme-specific residual draws from other held-out events. This is a
transparent first-order approximation for the operational model; the formal
production fit should pool repeated correction-layer imputations.
'''

from pathlib import Path

import numpy as np
import pandas as pd


PREDICTION_INPUT = Path('output/inla_spatiotemporal/cv_predictions.csv')
CORRECTION_INPUT = Path(
    'data/processed/noaa_dhw_correction_layer_local_first_validation.csv'
)
CONTEXT_INPUT = Path(
    'output/prediction_outliers/formal_ensemble_residuals.csv'
)
RRN_INPUT = Path('data/processed/rrn_pressure_reef_year.csv')
CYCLONE_INPUT = Path('data/processed/bom_cyclone_reef_year.csv')
COTS_INPUT = Path('data/gbrPredsAdj_20262408.csv')
DISEASE_INPUT = Path('data/processed/noaa_disease_risk_validation.csv')
EXPLANATION_INPUT = Path('data/curated/extreme_miss_explanations.csv')
OUTPUT_DIR = Path('output/local_calibrated_operational')

BASELINE = 'persistent_rw1_cyclone_cots_partial_pool'
LOCAL = (
    'persistent_rw1_local_dhw_decomposed_hazards_freshwater_partial_pool'
)
CURRENT_ONLY = (
    'persistent_rw1_operational_update_decomposed_hazards_freshwater_'
    'partial_pool'
)
SIMULATIONS = 4000
RANDOM_SEED = 20260828


def prediction_wide() -> pd.DataFrame:
    predictions = pd.read_csv(PREDICTION_INPUT)
    predictions = predictions[
        predictions['scheme'].eq('leave_one_event_out')
        & predictions['candidate'].isin([BASELINE, LOCAL, CURRENT_ONLY])
    ].copy()
    keys = [
        'programme_key', 'source_observation_id', 'ReefID', 'ReefName',
        'event_year', 'observed_mortality', 'observed_occurrence', 'fold',
    ]
    wide = predictions.pivot(
        index=keys, columns='candidate', values='predicted_mortality'
    ).reset_index()
    wide = wide.rename(columns={
        BASELINE: 'predicted_noaa',
        LOCAL: 'predicted_local_calibrated',
        CURRENT_ONLY: 'predicted_current_event_only',
    })
    required = {
        'predicted_noaa', 'predicted_local_calibrated',
        'predicted_current_event_only',
    }
    if not required.issubset(wide.columns):
        raise ValueError('Required INLA prediction candidates are missing')
    local_exposure = predictions[
        predictions['candidate'].eq(LOCAL)
    ][keys + [
        'ann_maxdhw_original', 'ann_maxdhw',
        'cyc_interval_maxHrs4mw', 'cyc_interval_peak_year',
        'tc_interval_min_distance_km', 'tc_interval_peak_name',
        'tc_interval_wind_distance_index',
        'cot_interval_idw_max', 'cots_outbreak_probability',
        'wqc_freqcc12', 'wqc_prior10_percentile', 'wqc_10yr_sum',
        'log_coastal_rain30',
    ]].rename(columns={
        'ann_maxdhw_original': 'noaa_dhw',
        'ann_maxdhw': 'corrected_dhw',
    })
    return wide.merge(local_exposure, on=keys, how='left', validate='one_to_one')


def add_correction_response(rows: pd.DataFrame) -> pd.DataFrame:
    layer = pd.read_csv(CORRECTION_INPUT)[[
        'ReefID', 'event_year', 'local_first_correction', 'correction_sd',
        'correction_source', 'correction_uncertainty_method',
        'nearest_local_logger_km', 'effective_local_loggers',
        'local_loggers_used',
    ]]
    rows = rows.merge(
        layer, on=['ReefID', 'event_year'], how='left', validate='many_to_one'
    )
    if rows['correction_sd'].isna().any():
        raise ValueError('Some mortality rows lack a correction uncertainty')

    usable = rows['local_first_correction'].abs().ge(0.25)
    rows['raw_mortality_per_dhw'] = np.where(
        usable,
        (
            rows['predicted_current_event_only'] - rows['predicted_noaa']
        ) / rows['local_first_correction'],
        np.nan,
    )
    rows['raw_mortality_per_dhw'] = rows[
        'raw_mortality_per_dhw'
    ].clip(lower=0)
    positive = rows['raw_mortality_per_dhw'].dropna()
    global_median = float(positive.median())

    rows['mortality_per_dhw'] = np.nan
    for programme, index in rows.groupby('programme_key').groups.items():
        values = rows.loc[index, 'raw_mortality_per_dhw'].dropna()
        if len(values) < 10:
            median = global_median
            upper = float(positive.quantile(.95))
        else:
            median = float(values.median())
            upper = float(values.quantile(.95))
        rows.loc[index, 'mortality_per_dhw'] = (
            rows.loc[index, 'raw_mortality_per_dhw']
            .clip(upper=upper)
            .fillna(median)
        )
    rows['dhw_uncertainty_mortality_sd'] = (
        rows['mortality_per_dhw'] * rows['correction_sd']
    )
    return rows


def simulate_intervals(rows: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(RANDOM_SEED)
    output = []
    for _, row in rows.iterrows():
        calibration = rows[
            rows['programme_key'].eq(row['programme_key'])
            & ~rows['event_year'].eq(row['event_year'])
        ]
        residuals = (
            calibration['observed_mortality']
            - calibration['predicted_local_calibrated']
        ).dropna().to_numpy()
        if len(residuals) < 20:
            residuals = (
                rows.loc[
                    ~rows['event_year'].eq(row['event_year']),
                    'observed_mortality',
                ]
                - rows.loc[
                    ~rows['event_year'].eq(row['event_year']),
                    'predicted_local_calibrated',
                ]
            ).dropna().to_numpy()
        correction_error = rng.normal(
            0, row['correction_sd'], SIMULATIONS
        ) * row['mortality_per_dhw']
        dhw_only = np.clip(
            row['predicted_local_calibrated'] + correction_error,
            0, 1,
        )
        event_radii = []
        for _, event_rows in calibration.groupby('event_year'):
            event_residuals = np.sort(np.abs(
                event_rows['observed_mortality']
                - event_rows['predicted_local_calibrated']
            ))
            event_rank = min(
                int(np.ceil((len(event_residuals) + 1) * .90)),
                len(event_residuals),
            )
            event_radii.append(float(event_residuals[event_rank - 1]))
        residual_radius = max(event_radii)
        dhw_radius = 1.644854 * row['dhw_uncertainty_mortality_sd']
        total_radius = float(np.sqrt(
            residual_radius ** 2 + dhw_radius ** 2
        ))
        output.append({
            'prediction_q05': max(
                0, row['predicted_local_calibrated'] - total_radius
            ),
            'prediction_q50': row['predicted_local_calibrated'],
            'prediction_q95': min(
                1, row['predicted_local_calibrated'] + total_radius
            ),
            'dhw_only_q05': float(np.quantile(dhw_only, .05)),
            'dhw_only_q95': float(np.quantile(dhw_only, .95)),
            'residual_conformal_radius90': residual_radius,
            'combined_radius90': total_radius,
            'residual_calibration_n': len(residuals),
        })
    intervals = pd.DataFrame(output, index=rows.index)
    result = pd.concat([rows, intervals], axis=1)
    result['interval_width90'] = (
        result['prediction_q95'] - result['prediction_q05']
    )
    result['dhw_interval_width90'] = (
        result['dhw_only_q95'] - result['dhw_only_q05']
    )
    result['covered90'] = (
        result['observed_mortality'].ge(result['prediction_q05'])
        & result['observed_mortality'].le(result['prediction_q95'])
    )
    result['residual'] = (
        result['observed_mortality']
        - result['predicted_local_calibrated']
    )
    result['outside_upper90'] = result['observed_mortality'].gt(
        result['prediction_q95']
    )
    result['outside_lower90'] = result['observed_mortality'].lt(
        result['prediction_q05']
    )
    return result


def calibration_summary(rows: pd.DataFrame) -> pd.DataFrame:
    grouping = ['programme_key', 'event_year', 'correction_source']
    detail = rows.groupby(grouping, as_index=False).agg(
        observations=('source_observation_id', 'size'),
        empirical_coverage90=('covered90', 'mean'),
        mean_interval_width90=('interval_width90', 'mean'),
        mean_dhw_interval_width90=('dhw_interval_width90', 'mean'),
        rmse=('residual', lambda x: float(np.sqrt(np.mean(x ** 2)))),
        bias_observed_minus_predicted=('residual', 'mean'),
    )
    overall = pd.DataFrame([{
        'programme_key': 'all',
        'event_year': 0,
        'correction_source': 'all',
        'observations': len(rows),
        'empirical_coverage90': rows['covered90'].mean(),
        'mean_interval_width90': rows['interval_width90'].mean(),
        'mean_dhw_interval_width90': rows['dhw_interval_width90'].mean(),
        'rmse': float(np.sqrt(np.mean(rows['residual'] ** 2))),
        'bias_observed_minus_predicted': rows['residual'].mean(),
    }])
    return pd.concat([overall, detail], ignore_index=True)


def add_ecological_context(reef_events: pd.DataFrame) -> pd.DataFrame:
    context = pd.read_csv(CONTEXT_INPUT)
    context = context.groupby(
        ['programme_key', 'ReefID', 'event_year'], as_index=False
    ).agg(
        prop_acropora_pre=('prop_acropora_pre', 'median'),
        observed_pre_cover=('observed_pre_cover', 'median'),
        era5_coastal_rain_max_30day=(
            'era5_coastal_rain_max_30day', 'median'
        ),
        era5_wind_calm_fraction=('era5_wind_calm_fraction', 'median'),
        freshwater_risk30_percentile=(
            'freshwater_risk30_percentile', 'median'
        ),
        disturbance_type=('DISTURBANCE_TYPE', 'first'),
        survey_description=('description', 'first'),
    )
    rrn = pd.read_csv(RRN_INPUT).rename(columns={'LABEL_ID': 'ReefID'})
    cyclone = pd.read_csv(CYCLONE_INPUT)[[
        'ReefID', 'event_year', 'tc_min_distance_km', 'tc_nearest_name',
        'tc_nearest_max_wind_ms', 'tc_storms_within300km',
    ]]
    cots = pd.read_csv(COTS_INPUT).rename(columns={
        'year': 'event_year', 'reefName': 'ReefName',
        'outbrProb': 'cots_outbreak_probability',
    })[['ReefName', 'event_year', 'cots_outbreak_probability']]
    disease = pd.read_csv(DISEASE_INPUT)[[
        'programme_key', 'ReefID', 'year', 'disease_risk_max',
        'disease_risk_days_ge1',
    ]].rename(columns={'year': 'event_year'})
    registry = pd.read_csv(EXPLANATION_INPUT).rename(columns={
        'reef_id': 'ReefID', 'reef_name': 'registry_reef_name',
    })
    result = (
        reef_events
        .merge(context, on=['programme_key', 'ReefID', 'event_year'],
               how='left', validate='one_to_one')
        .merge(rrn, on=['ReefID', 'event_year'], how='left',
               validate='many_to_one')
        .merge(cyclone, on=['ReefID', 'event_year'], how='left',
               validate='many_to_one')
        .merge(cots, on=['ReefName', 'event_year'], how='left',
               validate='many_to_one')
        .merge(disease, on=['programme_key', 'ReefID', 'event_year'],
               how='left', validate='one_to_one')
        .merge(registry, on=['ReefID', 'event_year'], how='left',
               validate='many_to_one')
    )
    result['priority_class'] = np.select(
        [
            result['outside_upper90'] & result['primary_explanation'].notna(),
            result['outside_upper90'],
            result['residual'].ge(.20),
            result['outside_lower90'],
        ],
        [
            'outside_90_explanation_registered',
            'outside_90_unexplained',
            'large_underprediction_within_interval',
            'outside_90_overprediction',
        ],
        default='lower_priority',
    )
    return result.sort_values(
        ['outside_upper90', 'residual'], ascending=[False, False]
    )


def reef_event_table(rows: pd.DataFrame) -> pd.DataFrame:
    grouped = rows.groupby(
        ['programme_key', 'ReefID', 'ReefName', 'event_year'], as_index=False
    ).agg(
        observations=('source_observation_id', 'size'),
        observed_mortality=('observed_mortality', 'mean'),
        predicted_noaa=('predicted_noaa', 'mean'),
        predicted_mortality=('predicted_local_calibrated', 'mean'),
        prediction_q05=('prediction_q05', 'mean'),
        prediction_q95=('prediction_q95', 'mean'),
        local_first_correction=('local_first_correction', 'first'),
        correction_sd=('correction_sd', 'first'),
        correction_source=('correction_source', 'first'),
        nearest_local_logger_km=('nearest_local_logger_km', 'first'),
        noaa_dhw=('noaa_dhw', 'first'),
        corrected_dhw=('corrected_dhw', 'first'),
        cyc_interval_maxHrs4mw=('cyc_interval_maxHrs4mw', 'first'),
        cyc_interval_peak_year=('cyc_interval_peak_year', 'first'),
        tc_interval_min_distance_km=(
            'tc_interval_min_distance_km', 'first'
        ),
        tc_interval_peak_name=('tc_interval_peak_name', 'first'),
        tc_interval_wind_distance_index=(
            'tc_interval_wind_distance_index', 'first'
        ),
        cot_interval_idw_max=('cot_interval_idw_max', 'first'),
        wqc_current=('wqc_freqcc12', 'first'),
        wqc_decadal_percentile=('wqc_prior10_percentile', 'first'),
        model_wqc_10yr_sum=('wqc_10yr_sum', 'first'),
        log_coastal_rain30=('log_coastal_rain30', 'first'),
    )
    grouped['residual'] = (
        grouped['observed_mortality'] - grouped['predicted_mortality']
    )
    grouped['prediction_gain_from_local_dhw'] = (
        grouped['predicted_mortality'] - grouped['predicted_noaa']
    )
    grouped['absolute_error'] = grouped['residual'].abs()
    grouped['outside_upper90'] = grouped['observed_mortality'].gt(
        grouped['prediction_q95']
    )
    grouped['outside_lower90'] = grouped['observed_mortality'].lt(
        grouped['prediction_q05']
    )
    return add_ecological_context(grouped)


def residual_feature_screen(rows: pd.DataFrame) -> pd.DataFrame:
    candidates = {
        'Local DHW correction uncertainty': 'correction_sd',
        'Pre-event Acropora proportion': 'prop_acropora_pre',
        'Pre-event coral cover': 'observed_pre_cover',
        'Coastal 30-day rainfall': 'era5_coastal_rain_max_30day',
        'Calm-wind fraction': 'era5_wind_calm_fraction',
        'Freshwater-risk percentile': 'freshwater_risk30_percentile',
        'WQC decadal percentile': 'wqc_prior10_percentile',
        'WQC continuous 10-year sum': 'model_wqc_10yr_sum',
        'Interval RRN cyclone wave hours': 'cyc_interval_maxHrs4mw',
        'Interval cyclone track distance': 'tc_interval_min_distance_km',
        'Interval cyclone wind-distance index': (
            'tc_interval_wind_distance_index'
        ),
        'Cyclone track distance': 'tc_min_distance_km',
        'Cyclone maximum wind': 'tc_nearest_max_wind_ms',
        'COTS outbreak probability': 'cots_outbreak_probability',
        'Disease risk maximum': 'disease_risk_max',
    }
    records = []
    severe = rows[rows['observed_mortality'].ge(.30)]
    for label, column in candidates.items():
        complete = severe[['residual', column]].dropna()
        if len(complete) < 8 or complete[column].nunique() < 3:
            continue
        records.append({
            'feature': label,
            'column': column,
            'n': len(complete),
            'spearman_with_residual': complete[column].corr(
                complete['residual'], method='spearman'
            ),
            'mean_residual_upper_quartile': complete.loc[
                complete[column].ge(complete[column].quantile(.75)),
                'residual',
            ].mean(),
            'mean_residual_other': complete.loc[
                complete[column].lt(complete[column].quantile(.75)),
                'residual',
            ].mean(),
        })
    result = pd.DataFrame(records)
    result['upper_quartile_residual_difference'] = (
        result['mean_residual_upper_quartile']
        - result['mean_residual_other']
    )
    return result.sort_values(
        'spearman_with_residual', key=lambda x: x.abs(), ascending=False
    )


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = simulate_intervals(add_correction_response(prediction_wide()))
    reef_events = reef_event_table(rows)
    calibration = calibration_summary(rows)
    features = residual_feature_screen(reef_events)

    rows.to_csv(OUTPUT_DIR / 'uncertainty_predictions.csv', index=False)
    calibration.to_csv(
        OUTPUT_DIR / 'uncertainty_calibration.csv', index=False
    )
    reef_events.to_csv(
        OUTPUT_DIR / 'largest_misses_local_calibrated.csv', index=False
    )
    reef_events.head(30).to_csv(
        OUTPUT_DIR / 'priority_misses_local_calibrated.csv', index=False
    )
    features.to_csv(
        OUTPUT_DIR / 'residual_feature_priorities.csv', index=False
    )
    pd.DataFrame([{
        'simulation_draws_per_observation': SIMULATIONS,
        'random_seed': RANDOM_SEED,
        'central_model': LOCAL,
        'dhw_sensitivity_contrast': CURRENT_ONLY,
        'residual_calibration': 'same_programme_other_events',
        'interval': (
            'other_event_conformal_radius_plus_delta_method_dhw_variance'
        ),
    }]).to_csv(OUTPUT_DIR / 'uncertainty_configuration.csv', index=False)

    print(calibration.head(1).to_string(index=False))
    print('\nLargest remaining underpredictions')
    print(reef_events[[
        'programme_key', 'ReefName', 'event_year', 'observed_mortality',
        'predicted_mortality', 'prediction_q95', 'residual',
        'outside_upper90', 'priority_class',
    ]].head(15).to_string(index=False))
    print('\nResidual feature screen')
    print(features.head(10).to_string(index=False))


if __name__ == '__main__':
    main()
