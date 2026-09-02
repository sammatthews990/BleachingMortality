'''Build a local-first, ENSO-aware correction to NOAA DHW.

Current-event reef loggers take precedence. Automated-station records are used
only when reef loggers are unavailable at that site/event. Historical evidence
is an ENSO-matched fallback and is never added to a supported local correction.
'''

from pathlib import Path
import re

import numpy as np
import pandas as pd

from build_event_dhw_correction_layer import haversine_matrix


SUMMARY_INPUT = Path('data/processed/aims_logger_noaa_dhw_validation.csv')
DAILY_INPUT = Path('data/processed/aims_logger_noaa_dhw_daily.csv.gz')
REEF_GRID_INPUT = Path('data/processed/environment_screen_grid_validation.csv')
ENSO_INPUT = Path('data/processed/enso_event_context.csv')
OUTPUT_DIR = Path('output/dhw_correction_local_first')
LAYER_OUTPUT = Path(
    'data/processed/noaa_dhw_correction_layer_local_first_validation.csv'
)

MIN_PEAK_DAYS = 70
LOCAL_RADII_KM = (10.0, 25.0, 50.0)
SELECTED_LOCAL_RADIUS_KM = 25.0
HISTORICAL_RADIUS_KM = 100.0
MAX_LOCAL_SOURCES = 4
DIRECT_MATCH_RADIUS_KM = 5.0
MIN_CONSISTENT_SIGN_FRACTION = 0.75
MIN_PHASE_MATCHED_EVENTS = 2


def peak_qualified_series() -> pd.DataFrame:
    '''Require coverage at the actual event DHW peak, not an earlier peak.'''
    summary = pd.read_csv(SUMMARY_INPUT)
    daily = pd.read_csv(DAILY_INPUT)
    keys = ['event_year', 'source', 'site', 'series', 'subsite']
    for frame in (summary, daily):
        frame['subsite'] = frame['subsite'].fillna('')
    valid = daily.dropna(subset=['logger_dhw'])
    peaks = valid.loc[
        valid.groupby(keys)['logger_dhw'].idxmax(),
        keys + ['logger_days_in_84d', 'logger_dhw', 'date'],
    ].rename(columns={
        'logger_days_in_84d': 'true_peak_days',
        'logger_dhw': 'unconstrained_peak_dhw',
        'date': 'unconstrained_peak_date',
    })
    result = summary.merge(peaks, on=keys, how='left', validate='one_to_one')
    result['peak_qualified'] = (
        result['adequate_dhw_coverage']
        & result['true_peak_days'].ge(MIN_PEAK_DAYS)
    )
    return result


def aggregate_site_events(series: pd.DataFrame) -> pd.DataFrame:
    '''Create one site-event discrepancy under a source and habitat hierarchy.'''
    records = []
    for (event_year, site), rows in series.groupby(['event_year', 'site']):
        reef = rows[rows['source'] == 'temperature_logger']
        source_rows = reef if not reef.empty else rows[
            rows['source'] == 'automated_weather'
        ]
        shallow = source_rows[
            source_rows['depth_class'].eq('shallow_0_5m')
            | source_rows['habitat_position'].eq('flat')
        ]
        used = shallow if not shallow.empty else source_rows
        if used.empty:
            continue
        records.append({
            'event_year': int(event_year),
            'site': site,
            'lat': float(used['lat'].median()),
            'lon': float(used['lon'].median()),
            'discrepancy': float(
                used['dhw_discrepancy_logger_minus_noaa'].median()
            ),
            'discrepancy_sd': float(
                used['dhw_discrepancy_logger_minus_noaa'].std()
            ),
            'logger_series': len(used),
            'median_depth_m': float(used['depth_m'].median()),
            'measurement_source': (
                'reef_logger' if not reef.empty else 'automated_fallback'
            ),
            'position_rule': (
                'shallow_or_flat_preferred'
                if not shallow.empty else 'all_available'
            ),
        })
    return pd.DataFrame(records)


def local_predict(
    sources: pd.DataFrame,
    targets: pd.DataFrame,
    value: str,
    radius_km: float,
    max_sources: int = MAX_LOCAL_SOURCES,
) -> pd.DataFrame:
    '''Predict only from sources inside a hard radius; return missing otherwise.'''
    n = len(targets)
    empty = pd.DataFrame({
        'prediction': np.full(n, np.nan),
        'nearest_km': np.full(n, np.nan),
        'effective_sources': np.zeros(n),
        'local_sd': np.full(n, np.nan),
        'sources_used': np.zeros(n, dtype=int),
    })
    if sources.empty:
        return empty
    distance = haversine_matrix(targets, sources)
    weights = np.exp(-distance / radius_km)
    weights[distance > radius_km] = 0
    if len(sources) > max_sources:
        distant = np.argsort(distance, axis=1)[:, max_sources:]
        weights[np.arange(n)[:, None], distant] = 0
    weight_sum = weights.sum(axis=1)
    values = sources[value].to_numpy(dtype=float)
    mean = np.divide(
        weights @ values, weight_sum,
        out=np.full(n, np.nan), where=weight_sum > 0,
    )
    square_sum = (weights ** 2).sum(axis=1)
    effective = np.divide(
        weight_sum ** 2, square_sum,
        out=np.zeros(n), where=square_sum > 0,
    )
    variance = np.divide(
        (weights * (values[None, :] - mean[:, None]) ** 2).sum(axis=1),
        weight_sum,
        out=np.full(n, np.nan), where=weight_sum > 0,
    )
    nearest = distance.min(axis=1)
    return pd.DataFrame({
        'prediction': mean,
        'nearest_km': nearest,
        'effective_sources': effective,
        'local_sd': np.sqrt(variance),
        'sources_used': (weights > 0).sum(axis=1),
    })


def normalise_site_name(value: str) -> str:
    '''Normalise survey and logger labels for conservative same-reef matching.'''
    value = re.sub(r'\([^)]*\)', ' ', str(value).lower())
    value = re.sub(r'\b(reef|island)\b', ' ', value)
    return re.sub(r'[^a-z0-9]+', ' ', value).strip()


def direct_site_predict(
    sources: pd.DataFrame, targets: pd.DataFrame, value: str
) -> pd.DataFrame:
    '''Return a direct site-event measurement before any spatial averaging.'''
    rows = []
    source_keys = sources.assign(
        site_key=sources['site'].map(normalise_site_name)
    )
    for _, target in targets.iterrows():
        target_key = normalise_site_name(target['ReefName'])
        matched = source_keys[source_keys['site_key'] == target_key]
        if matched.empty:
            rows.append({
                'prediction': np.nan,
                'distance_km': np.nan,
                'site': '',
                'measurement_sd': np.nan,
            })
            continue
        distances = haversine_matrix(
            target.to_frame().T, matched
        )[0]
        nearest_index = int(np.argmin(distances))
        nearest = matched.iloc[nearest_index]
        if distances[nearest_index] > DIRECT_MATCH_RADIUS_KM:
            rows.append({
                'prediction': np.nan,
                'distance_km': np.nan,
                'site': '',
                'measurement_sd': np.nan,
            })
            continue
        rows.append({
            'prediction': float(nearest[value]),
            'distance_km': float(distances[nearest_index]),
            'site': nearest['site'],
            'measurement_sd': float(nearest['discrepancy_sd']),
        })
    return pd.DataFrame(rows)


def enso_historical_effects(
    history: pd.DataFrame, target_phase: str
) -> pd.DataFrame:
    '''Use phase-matched sites only when their within-phase sign is consistent.'''
    matched = history[history['enso_phase'] == target_phase]
    if matched.empty:
        return pd.DataFrame(columns=[
            'site', 'lat', 'lon', 'effect', 'events', 'sign_fraction',
        ])
    effects = matched.groupby('site', as_index=False).agg(
        lat=('lat', 'median'),
        lon=('lon', 'median'),
        events=('event_year', 'nunique'),
        effect=('discrepancy', 'median'),
        positive_fraction=('discrepancy', lambda x: (x > 0).mean()),
    )
    effects['sign_fraction'] = np.maximum(
        effects['positive_fraction'], 1 - effects['positive_fraction']
    )
    return effects[
        effects['events'].ge(MIN_PHASE_MATCHED_EVENTS)
        & effects['sign_fraction'].ge(MIN_CONSISTENT_SIGN_FRACTION)
    ].copy()


def make_layer(
    site_events: pd.DataFrame,
    reef_grid: pd.DataFrame,
    local_radius_km: float,
    uncertainty: dict[str, float],
) -> pd.DataFrame:
    '''Use local current measurements first; ENSO history only as fallback.'''
    parts = []
    for event_year, targets in reef_grid.groupby('event_year', sort=True):
        phase = targets['enso_phase'].iloc[0]
        current = site_events[site_events['event_year'] == event_year]
        history = site_events[site_events['event_year'] < event_year]
        historical_sources = enso_historical_effects(history, phase)
        local = local_predict(
            current, targets, 'discrepancy', local_radius_km
        )
        direct = direct_site_predict(current, targets, 'discrepancy')
        historical = local_predict(
            historical_sources, targets, 'effect', HISTORICAL_RADIUS_KM
        )
        direct_supported = direct['prediction'].notna().to_numpy()
        local_supported = local['prediction'].notna().to_numpy()
        historical_supported = historical['prediction'].notna().to_numpy()
        correction = np.where(
            direct_supported,
            direct['prediction'].to_numpy(),
            np.where(
                local_supported,
                local['prediction'].to_numpy(),
                np.where(
                    historical_supported,
                    historical['prediction'].to_numpy(),
                    0,
                ),
            ),
        )
        source = np.where(
            direct_supported,
            'current_direct',
            np.where(
                local_supported,
                'current_local',
                np.where(
                    historical_supported, 'enso_historical_fallback', 'none'
                ),
            ),
        )
        result = targets.copy()
        result['direct_site_correction'] = direct['prediction'].to_numpy()
        result['direct_logger_site'] = direct['site'].to_numpy()
        result['local_event_correction'] = np.where(
            direct_supported,
            direct['prediction'].to_numpy(),
            local['prediction'].to_numpy(),
        )
        result['enso_historical_correction'] = (
            historical['prediction'].fillna(0).to_numpy()
        )
        result['local_first_correction'] = correction
        result['correction_source'] = source
        direct_measurement_sd = direct['measurement_sd'].to_numpy()
        direct_sd = np.where(
            np.isfinite(direct_measurement_sd),
            np.maximum(
                direct_measurement_sd,
                uncertainty['direct_measurement_sd_floor'],
            ),
            uncertainty['direct_measurement_sd_floor'],
        )
        interpolated_sd = np.sqrt(
            np.nan_to_num(local['local_sd'].to_numpy(), nan=0) ** 2
            + uncertainty['local_interpolation_sd_floor'] ** 2
        )
        result['correction_sd'] = np.where(
            direct_supported,
            direct_sd,
            np.where(
                local_supported,
                interpolated_sd,
                uncertainty['unsupported_correction_sd'],
            ),
        )
        result['correction_uncertainty_method'] = np.where(
            direct_supported,
            'replicate_or_pooled_direct_floor',
            np.where(
                local_supported,
                'source_spread_plus_heldout_floor',
                'unsupported_zero_mean_heldout_scale',
            ),
        )
        result['nearest_local_logger_km'] = np.where(
            direct_supported,
            direct['distance_km'].to_numpy(),
            local['nearest_km'].to_numpy(),
        )
        result['effective_local_loggers'] = np.where(
            direct_supported, 1, local['effective_sources'].to_numpy()
        )
        result['local_loggers_used'] = np.where(
            direct_supported, 1, local['sources_used'].to_numpy()
        )
        result['local_discrepancy_sd'] = np.where(
            direct_supported, 0, local['local_sd'].to_numpy()
        )
        result['nearest_historical_logger_km'] = (
            historical['nearest_km'].to_numpy()
        )
        result['effective_historical_loggers'] = (
            historical['effective_sources'].to_numpy()
        )
        result['historical_loggers_used'] = (
            historical['sources_used'].to_numpy()
        )
        result['enso_history_phase'] = phase
        result['current_logger_sites'] = current['site'].nunique()
        result['historical_phase_sites'] = len(historical_sources)
        parts.append(result)
    return pd.concat(parts, ignore_index=True)


def validation_predictions(
    site_events: pd.DataFrame,
    local_radius_km: float,
) -> pd.DataFrame:
    '''Leave each site out while preserving the local-first hierarchy.'''
    rows = []
    for index, target in site_events.iterrows():
        same_event = site_events[
            (site_events['event_year'] == target['event_year'])
            & (site_events['site'] != target['site'])
        ]
        local = local_predict(
            same_event, target.to_frame().T,
            'discrepancy', local_radius_km,
        ).iloc[0]
        history = site_events[
            site_events['event_year'] < target['event_year']
        ]
        historical_sources = enso_historical_effects(
            history, target['enso_phase']
        )
        historical = local_predict(
            historical_sources, target.to_frame().T,
            'effect', HISTORICAL_RADIUS_KM,
        ).iloc[0]
        if np.isfinite(local['prediction']):
            prediction = local['prediction']
            source = 'current_local'
        elif np.isfinite(historical['prediction']):
            prediction = historical['prediction']
            source = 'enso_historical_fallback'
        else:
            prediction = 0.0
            source = 'none'
        row = target.to_dict()
        row.update({
            'prediction': prediction,
            'correction_source': source,
            'nearest_local_logger_km': local['nearest_km'],
            'effective_local_loggers': local['effective_sources'],
            'local_loggers_used': local['sources_used'],
            'local_discrepancy_sd': local['local_sd'],
            'historical_prediction': (
                historical['prediction']
                if np.isfinite(historical['prediction']) else 0.0
            ),
            'local_radius_km': local_radius_km,
        })
        rows.append(row)
    return pd.DataFrame(rows)


def estimate_uncertainty_calibration(
    site_events: pd.DataFrame, validation: pd.DataFrame
) -> dict[str, float]:
    '''Estimate correction-error scales without using mortality outcomes.'''
    selected = validation[
        validation['local_radius_km'].eq(SELECTED_LOCAL_RADIUS_KM)
    ].copy()
    selected['error'] = selected['prediction'] - selected['discrepancy']
    local = selected[selected['correction_source'].eq('current_local')]
    unsupported = selected[selected['correction_source'].eq('none')]

    replicate_sd = site_events['discrepancy_sd'].dropna()
    replicate_sd = replicate_sd[replicate_sd > 0]
    direct_floor = float(replicate_sd.median())

    local_error_variance = float(np.mean(local['error'] ** 2))
    local_spread_variance = float(np.mean(
        np.nan_to_num(local['local_discrepancy_sd'], nan=0) ** 2
    ))
    interpolation_floor = float(np.sqrt(max(
        local_error_variance - local_spread_variance, 0
    )))
    unsupported_sd = float(np.sqrt(np.mean(unsupported['error'] ** 2)))
    return {
        'direct_measurement_sd_floor': direct_floor,
        'local_interpolation_sd_floor': interpolation_floor,
        'unsupported_correction_sd': unsupported_sd,
        'local_heldout_rmse': float(np.sqrt(local_error_variance)),
        'unsupported_heldout_rmse': unsupported_sd,
    }


def score_predictions(rows: pd.DataFrame) -> dict[str, float]:
    error = rows['prediction'] - rows['discrepancy']
    extreme = rows['discrepancy'].abs() >= 3
    return {
        'local_radius_km': float(rows['local_radius_km'].iloc[0]),
        'n': len(rows),
        'rmse': float(np.sqrt(np.mean(error ** 2))),
        'mae': float(np.mean(np.abs(error))),
        'bias': float(np.mean(error)),
        'extreme_n': int(extreme.sum()),
        'extreme_rmse': float(np.sqrt(np.mean(error[extreme] ** 2))),
        'local_support_fraction': float(
            (rows['correction_source'] == 'current_local').mean()
        ),
        'fallback_fraction': float(
            (rows['correction_source'] == 'enso_historical_fallback').mean()
        ),
        'no_correction_fraction': float(
            (rows['correction_source'] == 'none').mean()
        ),
    }


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    series_audit = peak_qualified_series()
    qualified = series_audit[series_audit['peak_qualified']].copy()
    site_events = aggregate_site_events(qualified)
    enso = pd.read_csv(ENSO_INPUT)[
        [
            'event_year', 'enso_phase', 'enso_category',
            'roni_bleaching_summer_mean',
        ]
    ]
    site_events = site_events.merge(
        enso, on='event_year', how='left', validate='many_to_one'
    )
    reef_grid = (
        pd.read_csv(REEF_GRID_INPUT)
        .rename(columns={'year': 'event_year'})
        .groupby(
            ['ReefID', 'ReefName', 'event_year'], as_index=False
        )
        .agg(lon=('lon', 'median'), lat=('lat', 'median'))
        .merge(enso, on='event_year', how='left', validate='many_to_one')
    )

    validations = []
    metrics = []
    for radius in LOCAL_RADII_KM:
        prediction = validation_predictions(site_events, radius)
        validations.append(prediction)
        metrics.append(score_predictions(prediction))
    validation = pd.concat(validations, ignore_index=True)
    metrics = pd.DataFrame(metrics)
    uncertainty = estimate_uncertainty_calibration(site_events, validation)
    layer = make_layer(
        site_events, reef_grid, SELECTED_LOCAL_RADIUS_KM, uncertainty
    )

    enso_summary = site_events.groupby(
        ['event_year', 'enso_phase', 'enso_category'], as_index=False
    ).agg(
        site_events=('site', 'size'),
        sites=('site', 'nunique'),
        median_discrepancy=('discrepancy', 'median'),
        mean_discrepancy=('discrepancy', 'mean'),
        positive_fraction=('discrepancy', lambda x: (x > 0).mean()),
        q10=('discrepancy', lambda x: x.quantile(.10)),
        q90=('discrepancy', lambda x: x.quantile(.90)),
    )
    source_summary = qualified.groupby(
        ['event_year', 'source'], as_index=False
    ).agg(
        qualified_series=('series', 'size'),
        sites=('site', 'nunique'),
        median_discrepancy=(
            'dhw_discrepancy_logger_minus_noaa', 'median'
        ),
    )

    series_audit.to_csv(
        OUTPUT_DIR / 'series_peak_coverage_audit.csv', index=False
    )
    site_events.to_csv(
        OUTPUT_DIR / 'site_event_discrepancies.csv', index=False
    )
    validation.to_csv(
        OUTPUT_DIR / 'leave_site_out_predictions.csv', index=False
    )
    metrics.to_csv(
        OUTPUT_DIR / 'local_radius_metrics.csv', index=False
    )
    enso_summary.to_csv(
        OUTPUT_DIR / 'enso_event_discrepancy_summary.csv', index=False
    )
    source_summary.to_csv(
        OUTPUT_DIR / 'qualified_source_summary.csv', index=False
    )
    pd.DataFrame([uncertainty]).to_csv(
        OUTPUT_DIR / 'uncertainty_calibration.csv', index=False
    )
    layer.to_csv(LAYER_OUTPUT, index=False)
    pd.DataFrame([{
        'minimum_days_at_actual_peak': MIN_PEAK_DAYS,
        'selected_local_radius_km': SELECTED_LOCAL_RADIUS_KM,
        'maximum_local_sources': MAX_LOCAL_SOURCES,
        'direct_match_radius_km': DIRECT_MATCH_RADIUS_KM,
        'historical_radius_km': HISTORICAL_RADIUS_KM,
        'minimum_consistent_sign_fraction':
            MIN_CONSISTENT_SIGN_FRACTION,
        'minimum_phase_matched_events': MIN_PHASE_MATCHED_EVENTS,
        'qualified_series_events': len(qualified),
        'site_events': len(site_events),
        'sites': site_events['site'].nunique(),
        **uncertainty,
    }]).to_csv(OUTPUT_DIR / 'configuration.csv', index=False)

    print('Peak-qualified source coverage')
    print(source_summary.to_string(index=False))
    print('\nLocal-radius validation')
    print(metrics.to_string(index=False))
    print('\nENSO event discrepancy summary')
    print(enso_summary.to_string(index=False))


if __name__ == '__main__':
    main()
