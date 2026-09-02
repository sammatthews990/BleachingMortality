'''Extract NOAA coral-disease outbreak-risk features for model reef-events.

The product uses negative values as masks when the preceding winter does not
meet its applicability criterion. Those values are not disease risk. This
script records applicability separately and assigns zero event risk when the
entire wet-season window is masked.
'''

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


DISEASE_URL = 's3://gbr-dms-data-public/noaa-crw-cdr-hdw-wdw/data.zarr'
SOURCE = 'NOAA CRW Daily 5-km Coral Disease Outbreak Risk v2.0'


def summarise_risk(values: np.ndarray) -> dict[str, np.ndarray]:
    values = np.asarray(values, dtype=float)
    if values.ndim == 1:
        values = values[:, None]
    valid = np.isfinite(values) & (values >= 0) & (values <= 15)
    risk = np.where(valid, values, np.nan)
    valid_days = valid.sum(axis=0)
    positive = np.where(valid, np.maximum(values, 0), 0)

    maximum = np.zeros(values.shape[1], dtype=float)
    q90 = np.zeros(values.shape[1], dtype=float)
    mean_positive = np.zeros(values.shape[1], dtype=float)
    for site in range(values.shape[1]):
        site_values = risk[:, site]
        finite = site_values[np.isfinite(site_values)]
        if len(finite):
            maximum[site] = finite.max()
            q90[site] = np.quantile(finite, 0.90)
        positive_values = finite[finite > 0]
        if len(positive_values):
            mean_positive[site] = positive_values.mean()

    return {
        'disease_risk_applicable_days': valid_days.astype(int),
        'disease_risk_applicable': (valid_days > 0).astype(int),
        'disease_risk_max': maximum,
        'disease_risk_q90': q90,
        'disease_risk_mean_positive': mean_positive,
        'disease_risk_days_positive': (positive > 0).sum(axis=0).astype(int),
        'disease_risk_days_ge1': (positive >= 1).sum(axis=0).astype(int),
        'disease_risk_burden': positive.sum(axis=0) / 7,
    }


def extract_year(dataset, rows: pd.DataFrame, year: int) -> pd.DataFrame:
    rows = rows.reset_index(drop=True).copy()
    lats = xr.DataArray(rows['lat'].to_numpy(), dims='site')
    lons = xr.DataArray(rows['lon'].to_numpy(), dims='site')
    values = (
        dataset['current_summer_outbreak_risk']
        .sel(time=slice(f'{year - 1}-11-01', f'{year}-04-30'))
        .sel(lat=lats, lon=lons, method='nearest')
        .compute()
        .values
    )
    result = rows.copy()
    for name, feature in summarise_risk(values).items():
        result[name] = feature
    result['disease_risk_source'] = SOURCE
    return result


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input', default='data/processed/environment_screen_grid_validation.csv'
    )
    parser.add_argument(
        '--output', default='data/processed/noaa_disease_risk_validation.csv'
    )
    parser.add_argument('--years', nargs='*', type=int)
    parser.add_argument('--force', action='store_true')
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = pd.read_csv(args.input)
    required = {'year', 'lon', 'lat'}
    missing = required - set(rows.columns)
    if missing:
        raise ValueError(f'Input lacks columns: {sorted(missing)}')
    rows['year'] = pd.to_numeric(rows['year'], errors='raise').astype(int)
    if args.years:
        rows = rows[rows['year'].isin(args.years)].copy()
    if rows.empty:
        raise ValueError('No reef-year rows selected')

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    part_dir = output.parent / f'{output.stem}_parts'
    part_dir.mkdir(parents=True, exist_ok=True)
    dataset = xr.open_zarr(DISEASE_URL, storage_options={'anon': True})

    parts = []
    for year in sorted(rows['year'].unique()):
        part = part_dir / f'{output.stem}_{year}.csv'
        if part.exists() and not args.force:
            result = pd.read_csv(part)
            print(f'Loaded cached disease risk for {year}', flush=True)
        else:
            result = extract_year(dataset, rows[rows['year'] == year], year)
            result.to_csv(part, index=False)
            print(f'Extracted disease risk for {year}', flush=True)
        parts.append(result)

    combined = pd.concat(parts, ignore_index=True)
    combined.to_csv(output, index=False)
    print(
        combined.groupby('year').agg(
            reef_years=('year', 'size'),
            applicable=('disease_risk_applicable', 'sum'),
            maximum_risk=('disease_risk_max', 'max'),
            median_risk=('disease_risk_max', 'median'),
        ).to_string(),
        flush=True,
    )


if __name__ == '__main__':
    main()
