'''Estimate persistent and within-event spatial corrections to NOAA DHW.

``rolling_prior`` uses only logger events before the target event.
``within_event`` adds an interpolation of current-event logger residuals.
Site effects are partially pooled across years and shrink toward zero where
the nearest logger is distant. Spatial ranges are selected by blocked tests.
'''

from pathlib import Path
import numpy as np
import pandas as pd

LOGGER_INPUT = Path('data/processed/aims_logger_noaa_dhw_validation.csv')
REEF_GRID_INPUT = Path('data/processed/environment_screen_grid_validation.csv')
OUTPUT_DIR = Path('output/dhw_correction')
LAYER_OUTPUT = Path('data/processed/noaa_dhw_correction_layer_validation.csv')
FUTURE_OUTPUT = Path('data/processed/noaa_dhw_correction_layer_future.csv')
RANGES_KM = (25.0, 50.0, 100.0, 200.0)
MAX_NEARBY_SOURCES = 8


def haversine_matrix(targets, sources):
    '''Return target-by-source great-circle distances in kilometres.'''
    tlat = np.radians(targets['lat'].to_numpy(dtype=float))[:, None]
    tlon = np.radians(targets['lon'].to_numpy(dtype=float))[:, None]
    slat = np.radians(sources['lat'].to_numpy(dtype=float))[None, :]
    slon = np.radians(sources['lon'].to_numpy(dtype=float))[None, :]
    dlat, dlon = slat - tlat, slon - tlon
    a = np.sin(dlat / 2) ** 2 + np.cos(tlat) * np.cos(slat) * np.sin(dlon / 2) ** 2
    return 6371.0088 * 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))


def aggregate_site_events(logger):
    '''Prefer shallow records within site-event; otherwise use all depths.'''
    records = []
    for (event_year, site), rows in logger.groupby(['event_year', 'site']):
        shallow = rows[
            (rows['depth_class'] == 'shallow_0_5m')
            | (rows['habitat_position'] == 'flat')
        ]
        used = shallow if not shallow.empty else rows
        records.append({
            'event_year': int(event_year), 'site': site,
            'lat': float(used['lat'].median()),
            'lon': float(used['lon'].median()),
            'discrepancy': float(used['dhw_discrepancy_logger_minus_noaa'].median()),
            'logger_series': len(used),
            'median_depth_m': float(used['depth_m'].median()),
            'position_rule': (
                'shallow_preferred' if not shallow.empty else 'all_available'
            ),
        })
    return pd.DataFrame(records)


def depth_diagnostics(logger):
    '''Create paired shallow/deep contrasts within the same site and event.'''
    known = logger[logger['depth_class'].isin(['shallow_0_5m', 'deep_gt5m'])]
    grouped = known.groupby(
        ['event_year', 'site', 'depth_class'], as_index=False
    ).agg(
        discrepancy=('dhw_discrepancy_logger_minus_noaa', 'median'),
        depth_m=('depth_m', 'median'), series=('series', 'nunique'),
    )
    discrepancy = grouped.pivot(
        index=['event_year', 'site'], columns='depth_class', values='discrepancy'
    )
    depth = grouped.pivot(
        index=['event_year', 'site'], columns='depth_class', values='depth_m'
    )
    pairs = discrepancy.join(depth, lsuffix='_dhw', rsuffix='_depth').reset_index()
    required = [
        'shallow_0_5m_dhw', 'deep_gt5m_dhw',
        'shallow_0_5m_depth', 'deep_gt5m_depth',
    ]
    pairs = pairs.dropna(subset=required).copy()
    pairs['shallow_minus_deep_dhw'] = (
        pairs['shallow_0_5m_dhw'] - pairs['deep_gt5m_dhw']
    )
    pairs['dhw_discrepancy_change_per_m'] = (
        pairs['deep_gt5m_dhw'] - pairs['shallow_0_5m_dhw']
    ) / (
        pairs['deep_gt5m_depth'] - pairs['shallow_0_5m_depth']
    )
    return grouped, pairs


def persistent_site_effects(rows):
    '''Partially pool repeated site discrepancies toward the training median.'''
    if rows.empty:
        return pd.DataFrame(columns=['site', 'lat', 'lon', 'effect', 'events'])
    centre = float(rows['discrepancy'].median())
    sites = rows.groupby('site', as_index=False).agg(
        lat=('lat', 'median'), lon=('lon', 'median'),
        events=('event_year', 'nunique'), site_median=('discrepancy', 'median'),
        site_iqr=('discrepancy', lambda x: x.quantile(.75) - x.quantile(.25)),
    )
    shrink = sites['events'] / (sites['events'] + 1.0)
    sites['effect'] = shrink * sites['site_median'] + (1 - shrink) * centre
    return sites


def spatial_predict(sources, targets, value, range_km):
    '''Interpolate values and expose local support diagnostics.'''
    if sources.empty:
        return pd.DataFrame({
            'prediction': np.zeros(len(targets)),
            'nearest_km': np.full(len(targets), np.nan),
            'effective_sources': np.zeros(len(targets)),
            'local_sd': np.full(len(targets), np.nan),
        })
    distance = haversine_matrix(targets, sources)
    weights = np.exp(-distance / range_km)
    weights[distance > 3 * range_km] = 0
    if len(sources) > MAX_NEARBY_SOURCES:
        distant = np.argsort(distance, axis=1)[:, MAX_NEARBY_SOURCES:]
        weights[np.arange(len(targets))[:, None], distant] = 0
    weight_sum = weights.sum(axis=1)
    values = sources[value].to_numpy()
    mean = np.divide(weights @ values, weight_sum, out=np.zeros(len(targets)), where=weight_sum > 0)
    nearest = distance.min(axis=1)
    prediction = mean * np.exp(-((nearest / (2 * range_km)) ** 2))
    square_sum = (weights ** 2).sum(axis=1)
    effective = np.divide(weight_sum ** 2, square_sum, out=np.zeros(len(targets)), where=square_sum > 0)
    variance = np.divide(
        (weights * (values[None, :] - mean[:, None]) ** 2).sum(axis=1),
        weight_sum, out=np.full(len(targets), np.nan), where=weight_sum > 0,
    )
    return pd.DataFrame({
        'prediction': prediction, 'nearest_km': nearest,
        'effective_sources': effective, 'local_sd': np.sqrt(variance),
    })


def score_predictions(rows, scheme, range_km):
    error = rows['prediction'] - rows['discrepancy']
    denominator = ((rows['discrepancy'] - rows['discrepancy'].mean()) ** 2).sum()
    return {
        'scheme': scheme, 'range_km': range_km, 'n': len(rows),
        'rmse': float(np.sqrt(np.mean(error ** 2))),
        'mae': float(np.mean(np.abs(error))), 'bias': float(np.mean(error)),
        'predictive_r2': float(1 - (error ** 2).sum() / denominator),
        'zero_correction_rmse': float(np.sqrt(np.mean(rows['discrepancy'] ** 2))),
    }


def persistent_cv(site_events, range_km):
    predictions = []
    for event_year in sorted(site_events['event_year'].unique()):
        assessment = site_events[site_events['event_year'] == event_year].copy()
        effects = persistent_site_effects(site_events[site_events['event_year'] != event_year])
        assessment['prediction'] = spatial_predict(effects, assessment, 'effect', range_km)['prediction'].to_numpy()
        assessment['fold'], assessment['scheme'] = str(event_year), 'leave_event_out'
        predictions.append(assessment)
    for site in sorted(site_events['site'].unique()):
        assessment = site_events[site_events['site'] == site].copy()
        effects = persistent_site_effects(site_events[site_events['site'] != site])
        assessment['prediction'] = spatial_predict(effects, assessment, 'effect', range_km)['prediction'].to_numpy()
        assessment['fold'], assessment['scheme'] = site, 'leave_site_out'
        predictions.append(assessment)
    result = pd.concat(predictions, ignore_index=True)
    metrics = pd.DataFrame([
        score_predictions(rows, scheme, range_km)
        for scheme, rows in result.groupby('scheme')
    ])
    return result, metrics


def event_updated_cv(site_events, persistent_range, update_range):
    '''Predict each site using other events plus other sites in its event.'''
    predictions = []
    for event_year in sorted(site_events['event_year'].unique()):
        current = site_events[site_events['event_year'] == event_year].copy()
        effects = persistent_site_effects(site_events[site_events['event_year'] != event_year])
        current['persistent_prediction'] = spatial_predict(
            effects, current, 'effect', persistent_range
        )['prediction'].to_numpy()
        current['event_residual'] = current['discrepancy'] - current['persistent_prediction']
        for index, target in current.iterrows():
            update_sources = current.drop(index=index).copy()
            event_offset = float(update_sources['event_residual'].median())
            update_sources['event_anomaly'] = (
                update_sources['event_residual'] - event_offset
            )
            local_update = spatial_predict(
                update_sources, target.to_frame().T,
                'event_anomaly', update_range,
            ).iloc[0]
            row = target.to_dict()
            row.update({
                'prediction': target['persistent_prediction'] +
                    event_offset + local_update['prediction'],
                'event_offset': event_offset,
                'local_event_update': local_update['prediction'],
                'event_update': event_offset + local_update['prediction'],
                'nearest_current_logger_km': local_update['nearest_km'],
                'effective_current_loggers': local_update['effective_sources'],
                'scheme': 'within_event_leave_site_out', 'fold': target['site'],
            })
            predictions.append(row)
    return pd.DataFrame(predictions)


def make_layer(site_events, reef_grid, persistent_range, update_range):
    parts = []
    for event_year, targets in reef_grid.groupby('event_year', sort=True):
        history = site_events[site_events['event_year'] < event_year]
        current = site_events[site_events['event_year'] == event_year].copy()
        effects = persistent_site_effects(history)
        persistent = spatial_predict(effects, targets, 'effect', persistent_range)
        rolling = persistent['prediction'].to_numpy()
        if current.empty:
            update = pd.DataFrame({
                'prediction': np.zeros(len(targets)), 'nearest_km': np.nan,
                'effective_sources': np.zeros(len(targets)), 'local_sd': np.nan,
            })
        else:
            current_base = spatial_predict(effects, current, 'effect', persistent_range)['prediction'].to_numpy()
            current['event_residual'] = current['discrepancy'] - current_base
            event_offset = float(current['event_residual'].median())
            current['event_anomaly'] = current['event_residual'] - event_offset
            update = spatial_predict(current, targets, 'event_anomaly', update_range)
        result = targets.copy()
        result['rolling_prior_correction'] = rolling
        result['event_offset'] = event_offset if not current.empty else 0
        result['local_event_update'] = update['prediction'].to_numpy()
        result['within_event_update'] = (
            result['event_offset'] + result['local_event_update']
        )
        result['within_event_correction'] = result['rolling_prior_correction'] + result['within_event_update']
        result['nearest_historical_logger_km'] = persistent['nearest_km'].to_numpy()
        result['effective_historical_loggers'] = persistent['effective_sources'].to_numpy()
        result['nearest_current_logger_km'] = np.asarray(update['nearest_km'])
        result['effective_current_loggers'] = np.asarray(update['effective_sources'])
        result['historical_events_available'] = history['event_year'].nunique()
        result['current_logger_sites'] = current['site'].nunique()
        parts.append(result)
    return pd.concat(parts, ignore_index=True)


def make_future_layer(site_events, reefs, persistent_range):
    effects = persistent_site_effects(site_events)
    prediction = spatial_predict(effects, reefs, 'effect', persistent_range)
    result = reefs.copy()
    result['persistent_correction'] = prediction['prediction'].to_numpy()
    result['nearest_logger_km'] = prediction['nearest_km'].to_numpy()
    result['effective_loggers'] = prediction['effective_sources'].to_numpy()
    result['local_sd'] = prediction['local_sd'].to_numpy()
    result['events_available'] = site_events['event_year'].nunique()
    return result


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    logger = pd.read_csv(LOGGER_INPUT)
    logger = logger[logger['adequate_dhw_coverage']].copy()
    depth_series, depth_pairs = depth_diagnostics(logger)
    site_events = aggregate_site_events(logger)
    screen_predictions, screen_metrics = [], []
    for range_km in RANGES_KM:
        predictions, metrics = persistent_cv(site_events, range_km)
        predictions['range_km'] = range_km
        screen_predictions.append(predictions)
        screen_metrics.append(metrics)
    persistent_metrics = pd.concat(screen_metrics, ignore_index=True)
    selected_persistent = float(
        persistent_metrics[persistent_metrics['scheme'] == 'leave_site_out']
        .sort_values('rmse').iloc[0]['range_km']
    )
    update_predictions, update_metrics = [], []
    for range_km in RANGES_KM:
        predictions = event_updated_cv(site_events, selected_persistent, range_km)
        predictions['range_km'] = range_km
        update_predictions.append(predictions)
        update_metrics.append(score_predictions(predictions, 'within_event_leave_site_out', range_km))
    update_metrics = pd.DataFrame(update_metrics)
    selected_update = float(update_metrics.sort_values('rmse').iloc[0]['range_km'])

    reef_grid = (
        pd.read_csv(REEF_GRID_INPUT).rename(columns={'year': 'event_year'})
        .groupby(['ReefID', 'ReefName', 'event_year'], as_index=False)
        .agg(lon=('lon', 'median'), lat=('lat', 'median'))
    )
    layer = make_layer(site_events, reef_grid, selected_persistent, selected_update)
    future_reefs = reef_grid.groupby(['ReefID', 'ReefName'], as_index=False).agg(
        lon=('lon', 'median'), lat=('lat', 'median')
    )
    future = make_future_layer(site_events, future_reefs, selected_persistent)
    persistence = site_events.groupby('site', as_index=False).agg(
        lat=('lat', 'median'), lon=('lon', 'median'), events=('event_year', 'nunique'),
        median_discrepancy=('discrepancy', 'median'),
        positive_fraction=('discrepancy', lambda x: (x > 0).mean()),
        discrepancy_iqr=('discrepancy', lambda x: x.quantile(.75) - x.quantile(.25)),
    )

    pd.concat(screen_predictions, ignore_index=True).to_csv(OUTPUT_DIR / 'persistent_cv_predictions.csv', index=False)
    pd.concat(update_predictions, ignore_index=True).to_csv(OUTPUT_DIR / 'within_event_cv_predictions.csv', index=False)
    pd.concat([persistent_metrics, update_metrics], ignore_index=True).to_csv(OUTPUT_DIR / 'range_screen_metrics.csv', index=False)
    persistence.to_csv(OUTPUT_DIR / 'site_persistence.csv', index=False)
    site_events.to_csv(OUTPUT_DIR / 'site_event_discrepancies.csv', index=False)
    depth_series.to_csv(OUTPUT_DIR / 'depth_class_discrepancies.csv', index=False)
    depth_pairs.to_csv(OUTPUT_DIR / 'paired_depth_diagnostics.csv', index=False)
    pd.DataFrame([{
        'paired_site_events': len(depth_pairs),
        'median_shallow_minus_deep_dhw': depth_pairs['shallow_minus_deep_dhw'].median(),
        'median_dhw_discrepancy_change_per_m': depth_pairs['dhw_discrepancy_change_per_m'].median(),
    }]).to_csv(OUTPUT_DIR / 'depth_effect_summary.csv', index=False)
    layer.to_csv(LAYER_OUTPUT, index=False)
    future.to_csv(FUTURE_OUTPUT, index=False)
    pd.DataFrame([{
        'selected_persistent_range_km': selected_persistent,
        'selected_within_event_range_km': selected_update,
        'site_events': len(site_events), 'sites': site_events['site'].nunique(),
        'events': site_events['event_year'].nunique(),
    }]).to_csv(OUTPUT_DIR / 'selected_configuration.csv', index=False)
    print('Persistent range screen')
    print(persistent_metrics.to_string(index=False))
    print('\nWithin-event range screen')
    print(update_metrics.to_string(index=False))
    print(f'\nSelected ranges: persistent={selected_persistent:.0f} km; within-event={selected_update:.0f} km')


if __name__ == '__main__':
    main()
