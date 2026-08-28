'''Build a GBR-wide AIMS logger versus NOAA CoralTemp DHW validation.

The public AIMS daily API is cached by bleaching event. Logger DHW is then
reconstructed with the collocated NOAA maximum monthly mean (MMM) and the NOAA
84-day DHW definition. Event summaries require at least 70 observed days in
the 84-day window ending at the logger-DHW peak; incomplete records remain in
the daily output but are not treated as reliable DHW comparisons.
'''

from __future__ import annotations

import json
import os
import time
from pathlib import Path
from urllib.parse import urlencode

import numpy as np
import pandas as pd
import requests
import xarray as xr


API = 'https://api.aims.gov.au/data-v2.0/10.25845/5b4eb0f9bb848/data/daily'
SST_URL = 's3://gbr-dms-data-public/noaa-crw-chs-sst/data.zarr'
MMM_URL = 's3://gbr-dms-data-public/noaa-crw-climatology/data.zarr'
DHW_URL = 's3://gbr-dms-data-public/noaa-crw-chs-dhw/data.zarr'
DEFAULT_EVENT_YEARS = (2016, 2017, 2020, 2022, 2024)
RAW_DIR = Path('data/raw/aims_temperature/events')
DAILY_OUTPUT = Path('data/processed/aims_logger_noaa_dhw_daily.csv.gz')
SUMMARY_OUTPUT = Path('data/processed/aims_logger_noaa_dhw_validation.csv')
MANIFEST_OUTPUT = Path('data/processed/aims_logger_noaa_dhw_manifest.json')
DEPLOYMENT_INPUT = Path(
    'data/raw/aims_temperature/deployments/temp-logger-deployments.csv'
)
AUTOMATED_INPUT = Path('data/processed/aims_temperature_daily_with_depth.csv.gz')
PAGE_SIZE = 10_000
MIN_PEAK_WINDOW_DAYS = 70
GBR_BOUNDS = {'min_lon': 142.0, 'max_lon': 154.0, 'min_lat': -25.0, 'max_lat': -10.0}


def event_dates(event_year: int) -> tuple[pd.Timestamp, pd.Timestamp]:
    return pd.Timestamp(event_year - 1, 11, 1), pd.Timestamp(event_year, 5, 31)


def fetch_event(event_year: int, refresh: bool = False) -> pd.DataFrame:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    path = RAW_DIR / f'aims_temperature_daily_event_{event_year}.csv.gz'
    if path.exists() and not refresh:
        return pd.read_csv(path, parse_dates=['date'])

    start, end = event_dates(event_year)
    url = API + '?' + urlencode({
        'from_date': start.date().isoformat(),
        'thru_date': end.date().isoformat(),
        'size': PAGE_SIZE,
    })
    pages: list[pd.DataFrame] = []
    seen: set[str] = set()
    session = requests.Session()
    while url and url not in seen:
        seen.add(url)
        for attempt in range(5):
            try:
                response = session.get(url, timeout=180)
                response.raise_for_status()
                payload = response.json()
                break
            except (requests.RequestException, json.JSONDecodeError):
                if attempt == 4:
                    raise
                time.sleep(2 ** attempt)
        rows = payload.get('results') or []
        if not rows:
            break
        pages.append(pd.DataFrame(rows))
        url = (payload.get('links') or {}).get('next')

    if not pages:
        raise RuntimeError(f'AIMS returned no daily rows for event {event_year}')
    data = pd.concat(pages, ignore_index=True)
    data['date'] = pd.to_datetime(data['time']).dt.normalize()
    data['logger_temperature'] = data['qc_val'].fillna(data['cal_val'])
    keep = [
        'date', 'site', 'series', 'subsite', 'lat', 'lon', 'deployment_id',
        'logger_temperature', 'qc_count',
    ]
    data = data[keep].dropna(subset=['date', 'series', 'lat', 'lon'])
    data.to_csv(path, index=False, compression='gzip')
    return data


def daily_series_median(rows: pd.DataFrame, event_year: int) -> pd.DataFrame:
    daily = (
        rows.groupby(
            ['date', 'source', 'site', 'series', 'subsite'], as_index=False
        )
        .agg(
            lat=('lat', 'median'),
            lon=('lon', 'median'),
            depth_m=('depth_m', 'median'),
            logger_temperature=('logger_temperature', 'median'),
            logger_observations=('qc_count', 'sum'),
            logger_deployments=('deployment_id', 'nunique'),
        )
    )
    daily['event_year'] = event_year
    return daily


def extract_noaa(daily: pd.DataFrame) -> pd.DataFrame:
    print('Opening NOAA CoralTemp and climatology stores', flush=True)
    sst_store = xr.open_zarr(SST_URL, storage_options={'anon': True})
    mmm_store = xr.open_zarr(MMM_URL, storage_options={'anon': True})
    dhw_store = xr.open_zarr(DHW_URL, storage_options={'anon': True})

    parts: list[pd.DataFrame] = []
    for event_year, event_rows in daily.groupby('event_year', sort=True):
        locations = (
            event_rows.groupby(['series_key'], as_index=False)
            .agg(lat=('lat', 'median'), lon=('lon', 'median'))
            .sort_values('series_key')
            .reset_index(drop=True)
        )
        points = xr.DataArray(np.arange(len(locations)), dims='point')
        lats = xr.DataArray(locations['lat'].to_numpy(), dims='point')
        lons = xr.DataArray(locations['lon'].to_numpy(), dims='point')
        start, end = event_dates(int(event_year))

        sst = (
            sst_store['analysed_sst']
            .sel(time=slice(start, end))
            .sel(lat=lats, lon=lons, method='nearest')
            .assign_coords(point=points)
            .compute()
        )
        dhw = (
            dhw_store['degree_heating_week']
            .sel(time=slice(start, end))
            .sel(lat=lats, lon=lons, method='nearest')
            .assign_coords(point=points)
            .compute()
        )
        mmm = (
            mmm_store['sst_clim_mmm']
            .isel(time=0)
            .sel(lat=lats, lon=lons, method='nearest')
            .assign_coords(point=points)
            .compute()
        )

        noaa = xr.Dataset({
            'noaa_sst': sst,
            'noaa_official_dhw': dhw,
        }).to_dataframe().reset_index()
        noaa['date'] = pd.to_datetime(noaa['time']).dt.normalize()
        noaa = noaa.merge(
            locations.reset_index(names='point'), on='point',
            how='left', validate='many_to_one'
        )
        noaa = noaa.merge(
            pd.DataFrame({
                'point': np.arange(len(locations)),
                'noaa_mmm': np.asarray(mmm),
            }), on='point', how='left', validate='many_to_one'
        )
        noaa['event_year'] = int(event_year)
        parts.append(noaa[[
            'event_year', 'series_key', 'date', 'noaa_sst',
            'noaa_official_dhw', 'noaa_mmm',
        ]])
        print(
            f'NOAA event {event_year}: {len(locations)} logger series',
            flush=True,
        )
    return pd.concat(parts, ignore_index=True)


def reconstruct_dhw(rows: pd.DataFrame) -> pd.DataFrame:
    start, end = event_dates(int(rows['event_year'].iloc[0]))
    calendar = pd.DataFrame({'date': pd.date_range(start, end, freq='D')})
    fixed = rows.iloc[0][
        [
            'event_year', 'source', 'series_key', 'series', 'site', 'subsite',
            'lat', 'lon', 'depth_m', 'noaa_mmm',
        ]
    ].to_dict()
    full = calendar.merge(
        rows.drop(columns=[
            'event_year', 'source', 'series_key', 'site', 'subsite',
            'lat', 'lon', 'depth_m'
        ]),
        on='date', how='left', validate='one_to_one'
    )
    for name, value in fixed.items():
        full[name] = value
    full['matched'] = (
        full['logger_temperature'].notna() & full['noaa_sst'].notna()
    )
    full['logger_hotspot'] = full['logger_temperature'] - full['noaa_mmm']
    full['logger_qualifying_hotspot'] = np.where(
        full['logger_hotspot'] >= 1, full['logger_hotspot'], 0
    )
    full.loc[full['logger_temperature'].isna(), 'logger_qualifying_hotspot'] = np.nan
    full['logger_dhw'] = (
        full['logger_qualifying_hotspot']
        .rolling(84, min_periods=1).sum() / 7
    )
    full['logger_days_in_84d'] = (
        full['logger_temperature'].notna().astype(int)
        .rolling(84, min_periods=1).sum()
    )
    full['logger_minus_noaa_sst'] = (
        full['logger_temperature'] - full['noaa_sst']
    )
    return full


def summarise_series(rows: pd.DataFrame) -> dict[str, object]:
    observed = rows.dropna(subset=['logger_temperature', 'noaa_sst'])
    eligible_peaks = rows[
        rows['logger_days_in_84d'] >= MIN_PEAK_WINDOW_DAYS
    ]
    adequate = not eligible_peaks.empty
    peak_pool = eligible_peaks if adequate else rows
    logger_peak = peak_pool.loc[peak_pool['logger_dhw'].idxmax()]
    official_peak = rows.loc[rows['noaa_official_dhw'].idxmax()]
    return {
        'event_year': int(rows['event_year'].iloc[0]),
        'source': rows['source'].iloc[0],
        'site': rows['site'].iloc[0],
        'series': rows['series'].iloc[0],
        'subsite': rows['subsite'].iloc[0],
        'lat': rows['lat'].iloc[0],
        'lon': rows['lon'].iloc[0],
        'depth_m': rows['depth_m'].iloc[0],
        'depth_class': (
            'unknown' if not np.isfinite(rows['depth_m'].iloc[0])
            else 'shallow_0_5m' if rows['depth_m'].iloc[0] <= 5
            else 'deep_gt5m'
        ),
        'habitat_position': (
            'flat' if 'FL' in rows['series'].iloc[0]
            else 'slope' if 'SL' in rows['series'].iloc[0]
            else 'other'
        ),
        'matched_days': len(observed),
        'logger_days_at_peak_84d': int(logger_peak['logger_days_in_84d']),
        'adequate_dhw_coverage': bool(adequate),
        'noaa_mmm': rows['noaa_mmm'].iloc[0],
        'logger_max_dhw': logger_peak['logger_dhw'],
        'logger_max_dhw_date': logger_peak['date'].date(),
        'noaa_official_max_dhw': official_peak['noaa_official_dhw'],
        'noaa_official_max_dhw_date': official_peak['date'].date(),
        'dhw_discrepancy_logger_minus_noaa': (
            logger_peak['logger_dhw'] - official_peak['noaa_official_dhw']
            if adequate else np.nan
        ),
        'mean_logger_minus_noaa_sst': observed['logger_minus_noaa_sst'].mean(),
        'median_logger_minus_noaa_sst': observed['logger_minus_noaa_sst'].median(),
    }


def main() -> None:
    refresh = os.environ.get('REFRESH_AIMS_TEMPERATURE') == '1'
    event_years = tuple(
        int(value.strip())
        for value in os.environ.get(
            'AIMS_EVENT_YEARS',
            ','.join(str(year) for year in DEFAULT_EVENT_YEARS),
        ).split(',')
        if value.strip()
    )
    raw = []
    raw_counts = {}
    gbr_counts = {}
    deployments = pd.read_csv(
        DEPLOYMENT_INPUT,
        usecols=['deployment_id', 'depth'],
        dtype={'deployment_id': str},
    ).rename(columns={'depth': 'depth_m'})
    deployments['depth_m'] = pd.to_numeric(
        deployments['depth_m'], errors='coerce'
    )
    for event_year in event_years:
        event = fetch_event(event_year, refresh=refresh)
        raw_counts[str(event_year)] = len(event)
        event = event[
            event['lon'].between(GBR_BOUNDS['min_lon'], GBR_BOUNDS['max_lon'])
            & event['lat'].between(GBR_BOUNDS['min_lat'], GBR_BOUNDS['max_lat'])
        ].copy()
        event['source'] = 'temperature_logger'
        event['deployment_id'] = event['deployment_id'].astype(str)
        event = event.merge(
            deployments, on='deployment_id', how='left',
            validate='many_to_one'
        )
        gbr_counts[str(event_year)] = len(event)
        raw.append(daily_series_median(event, event_year))
        print(
            f'AIMS GBR event {event_year}: {len(event):,} daily rows',
            flush=True,
        )
    daily = pd.concat(raw, ignore_index=True)
    if AUTOMATED_INPUT.exists():
        automated = pd.read_csv(AUTOMATED_INPUT, parse_dates=['date'])
        automated = automated[
            (automated['source'] == 'automated_weather')
            & automated['event_year'].isin(event_years)
        ].copy()
        if not automated.empty:
            automated = automated.rename(columns={
                'temperature_c': 'logger_temperature',
                'daily_observations': 'logger_observations',
            })
            automated['logger_deployments'] = 1
            daily = pd.concat([
                daily,
                automated[[
                    'date', 'source', 'site', 'series', 'subsite',
                    'lat', 'lon', 'depth_m', 'logger_temperature',
                    'logger_observations', 'logger_deployments', 'event_year',
                ]],
            ], ignore_index=True)
    # Automated products can reuse a display series name across deployments
    # and occasionally contribute multiple rows for one calendar date. Use a
    # stable source/site/series/subsite identity, then enforce one daily value
    # before NOAA colocation and rolling-DHW reconstruction.
    daily['subsite'] = daily['subsite'].fillna('')
    daily['series_key'] = (
        daily['source'].astype(str) + '::'
        + daily['site'].astype(str) + '::'
        + daily['series'].astype(str) + '::'
        + daily['subsite'].astype(str)
    )
    daily = (
        daily.groupby(
            [
                'date', 'event_year', 'series_key', 'source',
                'site', 'series', 'subsite',
            ],
            as_index=False,
        )
        .agg(
            lat=('lat', 'median'),
            lon=('lon', 'median'),
            depth_m=('depth_m', 'median'),
            logger_temperature=('logger_temperature', 'median'),
            logger_observations=('logger_observations', 'sum'),
            logger_deployments=('logger_deployments', 'sum'),
        )
    )
    noaa = extract_noaa(daily)
    joined = daily.merge(
        noaa, on=['event_year', 'series_key', 'date'],
        how='left', validate='many_to_one'
    )
    complete = pd.concat([
        reconstruct_dhw(rows)
        for _, rows in joined.groupby(['event_year', 'series_key'], sort=True)
    ], ignore_index=True)
    summary = pd.DataFrame([
        summarise_series(rows)
        for _, rows in complete.groupby(['event_year', 'series_key'], sort=True)
    ])

    DAILY_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    complete.to_csv(DAILY_OUTPUT, index=False, compression='gzip')
    summary.to_csv(SUMMARY_OUTPUT, index=False)
    manifest = {
        'aims_api': API,
        'noaa_sst': SST_URL,
        'noaa_mmm': MMM_URL,
        'noaa_dhw': DHW_URL,
        'event_years': list(event_years),
        'raw_row_counts': raw_counts,
        'gbr_bounds': GBR_BOUNDS,
        'gbr_row_counts': gbr_counts,
        'minimum_observed_days_in_peak_84d_window': MIN_PEAK_WINDOW_DAYS,
        'series_events': len(summary),
        'adequate_series_events': int(summary['adequate_dhw_coverage'].sum()),
    }
    MANIFEST_OUTPUT.write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(summary.groupby('event_year').agg(
        series_events=('series', 'size'),
        adequate=('adequate_dhw_coverage', 'sum'),
        median_dhw_discrepancy=('dhw_discrepancy_logger_minus_noaa', 'median'),
    ).to_string(), flush=True)


if __name__ == '__main__':
    main()
