'''Extract leakage-safe NOAA thermal intensity-duration-frequency metrics.

The fixed 1986-2012 reference predates every mortality-validation event. For
each reef and summer it contains the maximum mean NOAA CoralTemp SST over 1,
3, 7, 14, 28, 56 and 84 consecutive days. Event return periods use a smoothed
empirical exceedance probability at each duration. This is an IDF-style
screening product, not a fitted copula or an extrapolated GEV return level.
'''

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


SST_URL = 's3://gbr-dms-data-public/noaa-crw-chs-sst/data.zarr'
DURATIONS = (1, 3, 7, 14, 28, 56, 84)
BASELINE_EVENT_YEARS = tuple(range(1986, 2013))


def summer_block_maxima(
    dates, values, event_year: int, durations=DURATIONS
) -> dict[int, float]:
    '''Maximum rolling mean SST for each duration in a Nov-Apr summer.'''
    series = pd.Series(
        np.asarray(values, dtype=float),
        index=pd.DatetimeIndex(dates).normalize(), dtype=float,
    ).groupby(level=0).mean()
    summer = series.reindex(pd.date_range(
        f'{event_year - 1}-11-01', f'{event_year}-04-30', freq='D'
    ))
    maxima = {}
    for duration in durations:
        rolling = summer.rolling(duration, min_periods=duration).mean()
        maxima[duration] = (
            float(rolling.max()) if rolling.notna().any() else np.nan
        )
    return maxima


def empirical_return_period(reference, value: float) -> float:
    '''Smoothed empirical return period in years for an upper-tail value.'''
    reference = np.asarray(reference, dtype=float)
    reference = reference[np.isfinite(reference)]
    if not np.isfinite(value) or len(reference) == 0:
        return np.nan
    exceedances = int(np.sum(reference >= value))
    return float((len(reference) + 1) / (exceedances + 1))


def summarise_return_curve(reference: dict[int, np.ndarray], event: dict[int, float]):
    '''Create duration-specific and compressed empirical IDF fields.'''
    return_periods = {
        duration: empirical_return_period(reference[duration], event[duration])
        for duration in DURATIONS
    }
    result = {}
    for duration in DURATIONS:
        result[f'noaa_idf_maxmean_sst_{duration}d_c'] = event[duration]
        result[f'noaa_idf_return_period_{duration}d_years'] = return_periods[duration]
    finite = {d: r for d, r in return_periods.items() if np.isfinite(r)}
    result['noaa_idf_max_return_period_years'] = np.nan
    result['noaa_idf_mean_log_return_period'] = np.nan
    result['noaa_idf_critical_duration_days'] = np.nan
    result['noaa_idf_persistence_logrp_difference'] = np.nan
    if finite:
        critical = max(finite, key=finite.get)
        result['noaa_idf_max_return_period_years'] = finite[critical]
        result['noaa_idf_mean_log_return_period'] = float(np.mean(
            np.log(list(finite.values()))
        ))
        result['noaa_idf_critical_duration_days'] = float(critical)
        short = [np.log(finite[d]) for d in (1, 3, 7) if d in finite]
        long = [np.log(finite[d]) for d in (28, 56, 84) if d in finite]
        if short and long:
            result['noaa_idf_persistence_logrp_difference'] = float(
                np.mean(long) - np.mean(short)
            )
    return result


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--input', default='data/processed/environment_screen_grid_validation.csv'
    )
    parser.add_argument(
        '--output', default='data/processed/noaa_idf_metrics_validation.csv'
    )
    parser.add_argument('--force', action='store_true')
    return parser.parse_args()


def main():
    args = parse_args()
    output_path = Path(args.output)
    if output_path.exists() and not args.force:
        existing = pd.read_csv(output_path)
        required = {
            'noaa_idf_max_return_period_years',
            'noaa_idf_mean_log_return_period',
            'noaa_idf_persistence_logrp_difference',
        }
        if required.issubset(existing.columns):
            print(f'Using cached {output_path}', flush=True)
            return

    rows = pd.read_csv(args.input)
    required = {'ReefID', 'year', 'lon', 'lat'}
    missing = required - set(rows.columns)
    if missing:
        raise ValueError(f'Input lacks columns: {sorted(missing)}')
    if rows.duplicated(['ReefID', 'year']).any():
        raise ValueError('Input has duplicate ReefID-year keys')
    sites = rows.groupby('ReefID', as_index=False).agg(
        lon=('lon', 'median'), lat=('lat', 'median')
    ).sort_values('ReefID').reset_index(drop=True)
    lats = xr.DataArray(sites['lat'].to_numpy(), dims='site')
    lons = xr.DataArray(sites['lon'].to_numpy(), dims='site')

    print('Opening NOAA CoralTemp daily SST', flush=True)
    dataset = xr.open_zarr(SST_URL, storage_options={'anon': True})
    block = dataset['analysed_sst'].sel(
        time=slice('1985-11-01', f'{int(rows.year.max())}-04-30')
    ).sel(lat=lats, lon=lons, method='nearest').compute()
    dataset.close()
    dates = pd.DatetimeIndex(block.time.values)
    values = np.asarray(block.values, dtype=float)
    if np.nanmedian(values) > 100:
        values -= 273.15

    site_summaries = {}
    all_event_years = sorted(set(BASELINE_EVENT_YEARS) | set(rows.year.astype(int)))
    for site_index, reef_id in enumerate(sites.ReefID):
        maxima = {
            year: summer_block_maxima(dates, values[:, site_index], year)
            for year in all_event_years
        }
        reference = {
            duration: np.array([
                maxima[year][duration] for year in BASELINE_EVENT_YEARS
            ]) for duration in DURATIONS
        }
        site_summaries[reef_id] = {
            year: summarise_return_curve(reference, maxima[year])
            for year in rows.loc[rows.ReefID == reef_id, 'year'].astype(int)
        }

    records = []
    for row in rows.itertuples(index=False):
        record = row._asdict()
        record.update(site_summaries[row.ReefID][int(row.year)])
        records.append(record)
    output = pd.DataFrame(records).sort_values(['year', 'ReefID'])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(output_path, index=False)
    manifest = {
        'input': args.input,
        'output': args.output,
        'rows': int(len(output)),
        'sst_store': SST_URL,
        'source': 'NOAA Coral Reef Watch CoralTemp v3.1 daily 5-km nighttime SST',
        'baseline_event_years': [min(BASELINE_EVENT_YEARS), max(BASELINE_EVENT_YEARS)],
        'durations_days': list(DURATIONS),
        'return_period': '(N + 1) / (number of baseline annual maxima >= event + 1)',
        'frequency_model': 'smoothed empirical; no GEV extrapolation or copula',
    }
    output_path.with_suffix('.manifest.json').write_text(
        json.dumps(manifest, indent=2), encoding='utf-8'
    )
    print(f'Wrote {len(output)} rows to {output_path}', flush=True)


if __name__ == '__main__':
    main()
