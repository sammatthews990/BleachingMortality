'''Fetch and harmonise the pilot Queensland WMIP river-discharge layer.

This acquisition is intentionally outside the production pipeline. It writes
immutable raw response snapshots, station-day discharge, catchment-event
summaries, a preliminary distance-routed reef-event table, and an audit. The
routed layer is a freshwater-source sensitivity, never observed salinity.
'''

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd


BASE_URL = 'https://water-monitoring.information.qld.gov.au/cgi/webservice.exe'
DEFAULT_CROSSWALK = Path('config/freshwater_gauge_crosswalk.csv')
DEFAULT_REEF_GRID = Path('data/processed/freshwater_predictor_grid_2016_2025.csv')
DEFAULT_CHLA = Path('data/processed/sst_chla_features_validation.csv')
OUTPUT_DIR = Path('data/processed')
CACHE_DIR = Path('data/cache/wmip_discharge')
EVENT_YEARS = (2016, 2017, 2020, 2022, 2024)
UNAVAILABLE_LABELS = {'data not yet available'}


def wmip_url(payload: dict) -> str:
    '''Return the provider's JSON-in-query-string request URL.'''
    encoded = urllib.parse.quote(
        json.dumps(payload, separators=(',', ':')), safe=''
    )
    return f'{BASE_URL}?{encoded}'


def trace_payload(station_ids, start_time, end_time):
    return {
        'function': 'get_ts_traces',
        'version': '2',
        'params': {
            'site_list': ','.join(station_ids),
            'datasource': 'ATQ',
            'varfrom': '140.00',
            'varto': '140.00',
            'start_time': start_time,
            'end_time': end_time,
            'data_type': 'mean',
            'interval': 'day',
            'multiplier': '1',
            'report_time': 'end',
        },
    }


def fetch_json(payload, timeout=180):
    request = urllib.request.Request(
        wmip_url(payload), headers={'User-Agent': 'GBR-freshwater-pilot/1.0'}
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.load(response)


def quality_status(label):
    '''Map native quality text without discarding the native code or label.'''
    normalised = str(label or '').strip().lower()
    if normalised in UNAVAILABLE_LABELS or 'not available' in normalised:
        return 'unavailable'
    if 'provisional' in normalised or 'unverified' in normalised:
        return 'provisional'
    if normalised in {'', 'good', 'approved', 'validated'}:
        return 'validated'
    return 'flagged'


def parse_wmip_traces(response, snapshot_id, retrieved_at_utc, event_year):
    '''Convert a WMIP daily-mean response to the gauge-day contract.'''
    if int(response.get('error_num', 0)) != 0:
        raise RuntimeError(
            f"WMIP error {response.get('error_num')}: "
            f"{response.get('error_msg', 'unspecified')}"
        )
    traces = response.get('return', {}).get('traces', [])
    rows = []
    for series in traces:
        if int(series.get('error_num', 0)) != 0:
            raise RuntimeError(
                f"WMIP trace error for {series.get('site')}: "
                f"{series.get('error_num')}"
            )
        site = str(series.get('site', '')).strip()
        details = series.get('site_details', {})
        quality_codes = {
            str(key): value
            for key, value in series.get('quality_codes', {}).items()
        }
        units = series.get('varto_details', {}).get('units')
        if units and str(units).lower() not in {'cumecs', 'm3/s', 'm^3/s'}:
            raise ValueError(f'Unexpected discharge units for {site}: {units}')
        timezone_hours = float(details.get('timezone', 10.0))
        for point in series.get('trace', []):
            period_end = pd.to_datetime(str(point.get('t')), format='%Y%m%d%H%M%S')
            # report_time=end labels a daily mean at the start of the next day.
            date_local = (period_end - pd.Timedelta(days=1)).date()
            code = str(point.get('q', ''))
            label = quality_codes.get(code, '')
            status = quality_status(label)
            try:
                native_value = float(point.get('v'))
            except (TypeError, ValueError):
                native_value = math.nan
            value = native_value
            if status == 'unavailable' or not math.isfinite(native_value):
                value = math.nan
            if math.isfinite(value) and value < 0:
                raise ValueError(f'Negative discharge for {site} on {date_local}')
            rows.append({
                'snapshot_id': snapshot_id,
                'event_year': int(event_year),
                'station_id': site,
                'station_name': details.get('name'),
                'station_lon': pd.to_numeric(details.get('longitude'), errors='coerce'),
                'station_lat': pd.to_numeric(details.get('latitude'), errors='coerce'),
                'date_local': date_local.isoformat(),
                'period_end_local': period_end.isoformat(),
                'timezone_offset_hours': timezone_hours,
                'discharge_native_value': native_value,
                'discharge_mean_m3_s': value,
                'discharge_volume_ml_day': value * 86.4,
                'quality_code': code,
                'quality_label': label,
                'value_status': status,
                'datasource': 'ATQ',
                'varfrom': '140.00',
                'varto': '140.00',
                'retrieved_at_utc': retrieved_at_utc,
                'source_url': BASE_URL,
                'licence': response.get('licence'),
            })
    result = pd.DataFrame(rows)
    if result.empty:
        raise RuntimeError(f'WMIP returned no traces for event {event_year}')
    keys = ['snapshot_id', 'station_id', 'date_local']
    if result.duplicated(keys).any():
        raise RuntimeError('WMIP response has duplicate snapshot/station/day rows')
    return result


def event_window(event_year, cutoff):
    start = pd.Timestamp(event_year - 1, 12, 1)
    month, day = (3, 31) if cutoff == 'march' else (4, 30)
    end = pd.Timestamp(event_year, month, day)
    return start, end


def summarise_catchment_events(gauge_day, crosswalk):
    '''Aggregate audited gauges without applying model-time transformations.'''
    linked = gauge_day.merge(
        crosswalk,
        on=['station_id'],
        how='left',
        validate='many_to_one',
        suffixes=('', '_crosswalk'),
    )
    if linked['source_id'].isna().any():
        missing = sorted(linked.loc[linked['source_id'].isna(), 'station_id'].unique())
        raise RuntimeError(f'WMIP stations missing from crosswalk: {missing}')
    linked['date_local'] = pd.to_datetime(linked['date_local'])
    outputs = []
    for event_year in sorted(linked['event_year'].unique()):
        event_rows = linked.loc[linked['event_year'].eq(event_year)].copy()
        for cutoff in ('march', 'april'):
            start, end = event_window(int(event_year), cutoff)
            selected = event_rows.loc[event_rows['date_local'].between(start, end)]
            expected_days = int((end - start).days + 1)
            for source_id, group in selected.groupby('source_id', sort=True):
                group = group.sort_values('date_local')
                weights = group['station_weight'].astype(float)
                discharge = group['discharge_mean_m3_s'].astype(float) * weights
                discharge_by_day = discharge.groupby(group['date_local']).sum(min_count=1)
                rolling7 = discharge_by_day.rolling(7, min_periods=7).mean()
                rolling30 = discharge_by_day.rolling(30, min_periods=30).mean()
                first = group.iloc[0]
                observed = int(discharge_by_day.notna().sum())
                unavailable = int(group['value_status'].eq('unavailable').sum())
                flagged = int(group['value_status'].isin(['flagged', 'provisional']).sum())
                issue_time = (
                    end + pd.Timedelta(days=1, hours=-10)
                ).tz_localize('UTC')
                retrieved_time = pd.to_datetime(first['retrieved_at_utc'], utc=True)
                product_mode = (
                    'initial_forecast' if retrieved_time <= issue_time
                    else 'environmental_hindcast'
                )
                outputs.append({
                    'snapshot_id': first['snapshot_id'],
                    'source_id': source_id,
                    'basin_id': first['basin_id'],
                    'catchment_name': first['catchment_name'],
                    'event_year': int(event_year),
                    'product_mode': product_mode,
                    'cutoff_name': cutoff,
                    'window_start_local': start.date().isoformat(),
                    'window_end_local': end.date().isoformat(),
                    'issue_time_utc': issue_time.isoformat(),
                    'retrieved_after_issue': bool(retrieved_time > issue_time),
                    'station_ids': '|'.join(sorted(group['station_id'].unique())),
                    'expected_days': expected_days,
                    'observed_days': observed,
                    'unavailable_rows': unavailable,
                    'flagged_rows': flagged,
                    'coverage_fraction': observed / expected_days,
                    'latest_observation_local': (
                        discharge_by_day.dropna().index.max().date().isoformat()
                        if observed else None
                    ),
                    'discharge_total_ml': float((discharge_by_day * 86.4).sum(min_count=1)),
                    'discharge_max_1day_m3_s': float(discharge_by_day.max()),
                    'discharge_max_7day_mean_m3_s': float(rolling7.max()),
                    'discharge_max_30day_mean_m3_s': float(rolling30.max()),
                    'routing_lon': float(first['routing_lon']),
                    'routing_lat': float(first['routing_lat']),
                    'kernel_scale_km': float(first['kernel_scale_km']),
                    'max_distance_km': float(first['max_distance_km']),
                    'pilot_role': first['pilot_role'],
                    'selection_status': first['selection_status'],
                    'regulation_flag': first['regulation_flag'],
                    'rating_audit_status': first['rating_audit_status'],
                    'use_for_proxy': observed / expected_days >= 0.8,
                })
    return pd.DataFrame(outputs)


def haversine_km(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = map(
        np.radians, [lon1, lat1, lon2, lat2]
    )
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    return 6371.0088 * 2 * np.arcsin(np.minimum(1, np.sqrt(a)))


def build_reef_events(catchment_events, reef_grid, chla=None):
    '''Route source summaries with a declared preliminary distance kernel.'''
    reef_columns = ['ReefID', 'ReefName', 'reef_longitude', 'reef_latitude', 'event_year']
    reefs = reef_grid[reef_columns + [
        'era5_coastal_rain_dec_mar_max30', 'era5_coastal_rain_q1_total',
        'era5_coastal_rain_q1_max7', 'imos_kd490_q90', 'imos_secchi_p10',
        'ereefs_salinity_min', 'ereefs_freshwater_exposure30',
    ]].drop_duplicates(reef_columns)
    if chla is not None:
        chla_small = chla[['ReefID', 'year', 'chla_wetseason_median', 'chla_wetseason_n']].copy()
        chla_small = chla_small.groupby(['ReefID', 'year'], as_index=False).agg({
            'chla_wetseason_median': 'median', 'chla_wetseason_n': 'max'
        })
        reefs = reefs.merge(
            chla_small, left_on=['ReefID', 'event_year'], right_on=['ReefID', 'year'],
            how='left', validate='many_to_one'
        ).drop(columns='year')
    else:
        reefs['chla_wetseason_median'] = np.nan
        reefs['chla_wetseason_n'] = np.nan

    outputs = []
    for (event_year, cutoff), summaries in catchment_events.groupby(
        ['event_year', 'cutoff_name'], sort=True
    ):
        event_reefs = reefs.loc[reefs['event_year'].eq(event_year)].copy()
        contributions = []
        for source in summaries.itertuples(index=False):
            distance = haversine_km(
                event_reefs['reef_longitude'].to_numpy(),
                event_reefs['reef_latitude'].to_numpy(),
                source.routing_lon,
                source.routing_lat,
            )
            weight = np.exp(-distance / source.kernel_scale_km)
            weight[distance > source.max_distance_km] = 0
            contributions.append(pd.DataFrame({
                'ReefID': event_reefs['ReefID'].to_numpy(),
                'source_id': source.source_id,
                'distance_km': distance,
                'connection_weight': weight,
                'weighted_total_ml': weight * source.discharge_total_ml,
                'weighted_max7_m3_s': weight * source.discharge_max_7day_mean_m3_s,
                'source_usable': bool(source.use_for_proxy),
            }))
        contribution = pd.concat(contributions, ignore_index=True)
        usable = contribution.loc[contribution['source_usable']]
        aggregate = usable.groupby('ReefID', as_index=False).agg(
            wmip_routed_discharge_total_ml=('weighted_total_ml', 'sum'),
            wmip_routed_discharge_max7_m3_s=('weighted_max7_m3_s', 'sum'),
            wmip_connection_weight_sum=('connection_weight', 'sum'),
            wmip_nearest_source_distance_km=('distance_km', 'min'),
        )
        dominant = contribution.sort_values(
            ['ReefID', 'connection_weight'], ascending=[True, False]
        ).drop_duplicates('ReefID')[['ReefID', 'source_id']].rename(
            columns={'source_id': 'wmip_dominant_source_id'}
        )
        result = event_reefs.merge(aggregate, on='ReefID', how='left').merge(
            dominant, on='ReefID', how='left'
        )
        result['cutoff_name'] = cutoff
        result['product_mode'] = summaries['product_mode'].iloc[0]
        result['issue_time_utc'] = summaries['issue_time_utc'].iloc[0]
        result['wmip_routing_method'] = 'provisional exponential distance from configured outlet'
        result['freshwater_risk_percentile'] = np.nan
        result['salinity_observed'] = False
        result['use_in_primary_mortality_model'] = False
        outputs.append(result)
    return pd.concat(outputs, ignore_index=True)


def write_audit(gauge_day, catchment_events, crosswalk, snapshots, output_path):
    station_audit = []
    for station_id, group in gauge_day.groupby('station_id', sort=True):
        expected = sum(
            (event_window(int(year), 'april')[1] -
             event_window(int(year), 'april')[0]).days + 1
            for year in sorted(group['event_year'].unique())
        )
        station_audit.append({
            'station_id': station_id,
            'station_name': group['station_name'].dropna().iloc[0],
            'first_date': group['date_local'].min(),
            'last_date': group['date_local'].max(),
            'rows': int(len(group)),
            'finite_discharge_rows': int(group['discharge_mean_m3_s'].notna().sum()),
            'unavailable_rows': int(group['value_status'].eq('unavailable').sum()),
            'flagged_or_provisional_rows': int(group['value_status'].isin(['flagged', 'provisional']).sum()),
            'expected_event_window_days': expected,
        })
    audit = {
        'contract_version': 1,
        'created_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
        'scope': 'pilot only; no mortality refit or salinity claim',
        'mackay_control_correction': (
            'Mackay Reef 16-015 is near 16S; Daintree/Bloomfield are pilot '
            'sources. Pioneer is a separate southern API control.'
        ),
        'snapshots': snapshots,
        'crosswalk_rows': int(len(crosswalk)),
        'station_audit': station_audit,
        'catchment_event_rows': int(len(catchment_events)),
        'catchment_event_coverage': catchment_events[[
            'source_id', 'event_year', 'cutoff_name', 'coverage_fraction',
            'flagged_rows', 'regulation_flag', 'rating_audit_status',
            'use_for_proxy'
        ]].to_dict('records'),
        'promotion_blockers': [
            'routing outlet coordinates are provisional',
            'rating tables and high-flow rating limits have not been independently reviewed',
            'regulated and ungauged catchment fractions require hydrological review',
            'distance routing has not been validated against plume transport',
            'freshwater proxy must pass event- and catchment-held-out low-salinity validation',
        ],
    }
    output_path.write_text(json.dumps(audit, indent=2), encoding='utf-8')


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--crosswalk', type=Path, default=DEFAULT_CROSSWALK)
    parser.add_argument('--reef-grid', type=Path, default=DEFAULT_REEF_GRID)
    parser.add_argument('--chla', type=Path, default=DEFAULT_CHLA)
    parser.add_argument('--event-years', nargs='+', type=int, default=list(EVENT_YEARS))
    return parser.parse_args()


def main():
    args = parse_args()
    crosswalk = pd.read_csv(args.crosswalk, dtype={'station_id': str})
    required = {
        'source_id', 'basin_id', 'catchment_name', 'station_id',
        'station_weight', 'routing_lon', 'routing_lat', 'kernel_scale_km',
        'max_distance_km', 'pilot_role', 'selection_status',
        'regulation_flag', 'rating_audit_status'
    }
    missing = required.difference(crosswalk.columns)
    if missing:
        raise ValueError(f'Crosswalk missing columns: {sorted(missing)}')
    if crosswalk['station_id'].duplicated().any():
        raise ValueError('Each pilot station must map to one source')

    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    retrieved_at = dt.datetime.now(dt.timezone.utc).replace(microsecond=0)
    retrieved_at_text = retrieved_at.isoformat()
    gauge_outputs = []
    snapshots = []
    station_ids = sorted(crosswalk['station_id'].astype(str).tolist())
    for event_year in args.event_years:
        start = f'{event_year - 1}1201000000'
        end = f'{event_year}0501000000'
        payload = trace_payload(station_ids, start, end)
        print(f'Fetching WMIP discharge for event {event_year}', flush=True)
        response = fetch_json(payload)
        content = json.dumps(response, sort_keys=True, separators=(',', ':')).encode('utf-8')
        digest = hashlib.sha256(content).hexdigest()
        snapshot_id = f"{retrieved_at.strftime('%Y%m%dT%H%M%SZ')}_{event_year}_{digest[:12]}"
        snapshot_path = CACHE_DIR / f'{snapshot_id}.json'
        if snapshot_path.exists():
            if snapshot_path.read_bytes() != content:
                raise RuntimeError(f'Immutable snapshot collision: {snapshot_path}')
        else:
            snapshot_path.write_bytes(content)
        parsed = parse_wmip_traces(
            response, snapshot_id, retrieved_at_text, event_year
        )
        window_start, window_end = event_window(event_year, 'april')
        parsed_dates = pd.to_datetime(parsed['date_local'])
        parsed = parsed.loc[
            parsed_dates.between(window_start, window_end)
        ].copy()
        gauge_outputs.append(parsed)
        snapshots.append({
            'snapshot_id': snapshot_id,
            'event_year': event_year,
            'sha256': digest,
            'path': str(snapshot_path),
            'request': payload,
        })

    gauge_day = pd.concat(gauge_outputs, ignore_index=True)
    catchment_events = summarise_catchment_events(gauge_day, crosswalk)
    gauge_metadata = crosswalk[[
        'station_id', 'source_id', 'basin_id', 'catchment_name'
    ]].rename(columns={'source_id': 'catchment_id'})
    gauge_day = gauge_day.merge(
        gauge_metadata, on='station_id', how='left', validate='many_to_one'
    )
    gauge_day['coordinate_datum'] = 'GDA94'

    reef_grid = pd.read_csv(args.reef_grid, low_memory=False)
    chla = pd.read_csv(args.chla) if args.chla.exists() else None
    reef_events = build_reef_events(catchment_events, reef_grid, chla)

    gauge_path = OUTPUT_DIR / 'freshwater_gauge_day.csv'
    catchment_path = OUTPUT_DIR / 'freshwater_catchment_event.csv'
    reef_path = OUTPUT_DIR / 'freshwater_reef_event.csv'
    audit_path = OUTPUT_DIR / 'freshwater_ingestion_audit.json'
    gauge_day.to_csv(gauge_path, index=False)
    catchment_events.to_csv(catchment_path, index=False)
    reef_events.to_csv(reef_path, index=False)
    write_audit(gauge_day, catchment_events, crosswalk, snapshots, audit_path)
    print(f'Wrote {len(gauge_day):,} gauge-days to {gauge_path}')
    print(f'Wrote {len(catchment_events):,} catchment-events to {catchment_path}')
    print(f'Wrote {len(reef_events):,} reef-events to {reef_path}')
    print(f'Wrote audit to {audit_path}')


if __name__ == '__main__':
    main()
