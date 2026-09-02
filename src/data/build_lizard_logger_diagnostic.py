'''Compare AIMS in-water Lizard temperature loggers with NOAA CoralTemp.

The AIMS daily endpoint is public and the two input JSON files are intentionally
kept as immutable raw data. Multiple concurrent logger deployments are reduced
to a daily subsite median before comparison with the same nearest 0.05-degree
NOAA cell used by the environmental feature pipeline.
'''

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


SST_URL = 's3://gbr-dms-data-public/noaa-crw-chs-sst/data.zarr'
MMM_URL = 's3://gbr-dms-data-public/noaa-crw-climatology/data.zarr'
DHW_URL = 's3://gbr-dms-data-public/noaa-crw-chs-dhw/data.zarr'
INPUTS = [
    Path('data/raw/aims_temperature/lizard_LIZSL1_daily_2023_2024.json'),
    Path('data/raw/aims_temperature/lizard_LIZFL1_daily_2023_2024.json'),
    Path('data/raw/aims_temperature/north_direction_FL1_daily_2023_2024.json'),
    Path('data/raw/aims_temperature/north_direction_SL1_daily_2023_2024.json'),
]
OUTPUT = Path('data/processed/lizard_cluster_logger_noaa_daily_2024.csv')
SUMMARY = Path('data/processed/lizard_cluster_logger_noaa_summary_2024.csv')
MONTHLY = Path('data/processed/lizard_cluster_logger_noaa_monthly_2024.csv')
DHW_SUMMARY = Path(
    'data/processed/lizard_cluster_logger_dhw_uplift_2024.csv'
)
BLEACHING_START = pd.Timestamp('2023-11-01')
BLEACHING_END = pd.Timestamp('2024-04-30')


def read_aims_daily(path: Path) -> pd.DataFrame:
    with path.open(encoding='utf-8') as stream:
        payload = json.load(stream)
    rows = pd.DataFrame(payload['results'])
    rows['date'] = pd.to_datetime(rows['time'])
    rows['logger_temperature'] = rows['qc_val'].fillna(rows['cal_val'])
    return rows[[
        'date', 'site', 'series', 'subsite', 'lat', 'lon', 'deployment_id',
        'logger_temperature', 'qc_min', 'qc_max', 'qc_count',
    ]]


def main() -> None:
    missing = [str(path) for path in INPUTS if not path.exists()]
    if missing:
        raise FileNotFoundError(f'Missing AIMS logger inputs: {missing}')

    logger = pd.concat(
        [read_aims_daily(path) for path in INPUTS], ignore_index=True
    )
    daily = (
        logger.groupby(
            ['date', 'site', 'series', 'subsite', 'lat', 'lon'], as_index=False
        )
        .agg(
            logger_temperature=('logger_temperature', 'median'),
            logger_daily_min=('qc_min', 'min'),
            logger_daily_max=('qc_max', 'max'),
            logger_deployments=('deployment_id', 'nunique'),
            logger_observations=('qc_count', 'sum'),
        )
        .sort_values(['series', 'date'])
    )

    print('Opening NOAA CoralTemp DMS store', flush=True)
    ds = xr.open_zarr(SST_URL, storage_options={'anon': True})
    climatology = xr.open_zarr(MMM_URL, storage_options={'anon': True})
    dhw_store = xr.open_zarr(DHW_URL, storage_options={'anon': True})
    extracted = []
    for series, rows in daily.groupby('series', sort=True):
        lat = float(rows['lat'].iloc[0])
        lon = float(rows['lon'].iloc[0])
        noaa = (
            ds['analysed_sst']
            .sel(time=slice(rows['date'].min(), rows['date'].max()))
            .sel(lat=lat, lon=lon, method='nearest')
            .compute()
            .to_dataframe(name='noaa_sst')
            .reset_index()
        )
        noaa['date'] = pd.to_datetime(noaa['time']).dt.normalize()
        noaa['series'] = series
        official_dhw = (
            dhw_store['degree_heating_week']
            .sel(time=slice(rows['date'].min(), rows['date'].max()))
            .sel(lat=lat, lon=lon, method='nearest')
            .compute()
            .to_dataframe(name='noaa_official_dhw')
            .reset_index()
        )
        official_dhw['date'] = pd.to_datetime(
            official_dhw['time']
        ).dt.normalize()
        noaa = noaa.merge(
            official_dhw[['date', 'noaa_official_dhw']],
            on='date', how='left', validate='one_to_one'
        )
        mmm = float(
            climatology['sst_clim_mmm']
            .isel(time=0)
            .sel(lat=lat, lon=lon, method='nearest')
            .compute()
        )
        noaa['noaa_mmm'] = mmm
        extracted.append(noaa[[
            'date', 'series', 'noaa_sst', 'noaa_official_dhw',
            'noaa_mmm', 'lat', 'lon'
        ]].rename(columns={
            'lat': 'noaa_grid_lat', 'lon': 'noaa_grid_lon'
        }))

    comparison = daily.merge(
        pd.concat(extracted, ignore_index=True),
        on=['date', 'series'], how='left', validate='one_to_one'
    )
    comparison['logger_minus_noaa'] = (
        comparison['logger_temperature'] - comparison['noaa_sst']
    )
    comparison['logger_7day'] = comparison.groupby('series')[
        'logger_temperature'
    ].transform(lambda values: values.rolling(7, min_periods=5).mean())
    comparison['noaa_7day'] = comparison.groupby('series')[
        'noaa_sst'
    ].transform(lambda values: values.rolling(7, min_periods=5).mean())
    comparison['delta_7day'] = (
        comparison['logger_7day'] - comparison['noaa_7day']
    )
    comparison['logger_hotspot'] = (
        comparison['logger_temperature'] - comparison['noaa_mmm']
    )
    comparison['noaa_hotspot_reconstructed'] = (
        comparison['noaa_sst'] - comparison['noaa_mmm']
    )
    comparison['logger_qualifying_hotspot'] = np.where(
        comparison['logger_hotspot'] >= 1,
        comparison['logger_hotspot'], 0
    )
    comparison['noaa_qualifying_hotspot'] = np.where(
        comparison['noaa_hotspot_reconstructed'] >= 1,
        comparison['noaa_hotspot_reconstructed'], 0
    )
    comparison['logger_dhw_reconstructed'] = comparison.groupby('series')[
        'logger_qualifying_hotspot'
    ].transform(lambda values: values.rolling(84, min_periods=1).sum() / 7)
    comparison['noaa_dhw_reconstructed'] = comparison.groupby('series')[
        'noaa_qualifying_hotspot'
    ].transform(lambda values: values.rolling(84, min_periods=1).sum() / 7)

    summary_rows = []
    for series, rows in comparison.groupby('series', sort=True):
        valid = rows.dropna(subset=['logger_temperature', 'noaa_sst']).copy()
        peak_cut = valid['noaa_sst'].quantile(0.75)
        peak = valid[valid['noaa_sst'] >= peak_cut]
        summary_rows.append({
            'series': series,
            'site': valid['site'].iloc[0],
            'subsite': valid['subsite'].iloc[0],
            'logger_depth_context': {
                'LIZFL1': '5 m reef flat',
                'LIZSL1': '10 m reef slope',
                'NDIRECTIONFL1': '2 m reef flat',
                'NDIRECTIONSL1': '9 m reef slope',
            }.get(series, 'see AIMS deployment index'),
            'first_date': valid['date'].min().date(),
            'last_date': valid['date'].max().date(),
            'matched_days': len(valid),
            'mean_logger_temperature': valid['logger_temperature'].mean(),
            'mean_noaa_sst': valid['noaa_sst'].mean(),
            'mean_logger_minus_noaa': valid['logger_minus_noaa'].mean(),
            'median_logger_minus_noaa': valid['logger_minus_noaa'].median(),
            'q10_logger_minus_noaa': valid['logger_minus_noaa'].quantile(0.10),
            'q90_logger_minus_noaa': valid['logger_minus_noaa'].quantile(0.90),
            'mean_delta_on_warmest_noaa_quartile': peak['logger_minus_noaa'].mean(),
            'max_7day_logger_minus_noaa': valid['delta_7day'].max(),
            # This is an offset integral, not DHW: no MMM threshold is applied.
            'positive_offset_degree_weeks': np.maximum(
                valid['logger_minus_noaa'], 0
            ).sum() / 7,
        })

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    comparison.to_csv(OUTPUT, index=False)
    pd.DataFrame(summary_rows).to_csv(SUMMARY, index=False)
    monthly = (
        comparison.dropna(subset=['logger_temperature', 'noaa_sst'])
        .assign(month=lambda rows: rows['date'].dt.to_period('M').astype(str))
        .groupby(['site', 'series', 'month'], as_index=False)
        .agg(
            matched_days=('logger_minus_noaa', 'size'),
            mean_logger_minus_noaa=('logger_minus_noaa', 'mean'),
            signed_offset_degree_weeks=(
                'logger_minus_noaa', lambda values: values.sum() / 7
            ),
            positive_offset_degree_weeks=(
                'logger_minus_noaa',
                lambda values: np.maximum(values, 0).sum() / 7,
            ),
            max_7day_logger_minus_noaa=('delta_7day', 'max'),
        )
    )
    monthly.to_csv(MONTHLY, index=False)
    event = comparison[
        comparison['date'].between(BLEACHING_START, BLEACHING_END)
    ].copy()
    dhw_rows = []
    for series, rows in event.groupby('series', sort=True):
        logger_peak = rows.loc[rows['logger_dhw_reconstructed'].idxmax()]
        noaa_peak = rows.loc[rows['noaa_dhw_reconstructed'].idxmax()]
        official_peak = rows.loc[rows['noaa_official_dhw'].idxmax()]
        dhw_rows.append({
            'site': rows['site'].iloc[0],
            'series': series,
            'subsite': rows['subsite'].iloc[0],
            'bleaching_start': max(BLEACHING_START, rows['date'].min()).date(),
            'bleaching_end': min(BLEACHING_END, rows['date'].max()).date(),
            'matched_days': len(rows),
            'noaa_mmm': rows['noaa_mmm'].iloc[0],
            'logger_max_dhw': logger_peak['logger_dhw_reconstructed'],
            'logger_max_dhw_date': logger_peak['date'].date(),
            'noaa_reconstructed_max_dhw': noaa_peak['noaa_dhw_reconstructed'],
            'noaa_reconstructed_max_dhw_date': noaa_peak['date'].date(),
            'noaa_official_max_dhw': official_peak['noaa_official_dhw'],
            'noaa_official_max_dhw_date': official_peak['date'].date(),
            'dhw_uplift_vs_reconstructed_noaa': (
                logger_peak['logger_dhw_reconstructed'] -
                noaa_peak['noaa_dhw_reconstructed']
            ),
            'dhw_uplift_vs_official_noaa': (
                logger_peak['logger_dhw_reconstructed'] -
                official_peak['noaa_official_dhw']
            ),
        })
    pd.DataFrame(dhw_rows).to_csv(DHW_SUMMARY, index=False)
    print(pd.DataFrame(summary_rows).to_string(index=False), flush=True)


if __name__ == '__main__':
    main()
