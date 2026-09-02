'''Extract daily SST-distribution and median chlorophyll features from GBR DMS.

Primary features are aligned to ecological event windows:

* ``sst_summer_*``: 1 November of the previous year through 30 April.
* ``sst_prevyear_*``: the complete calendar year before the response.
* ``chla_wetseason_*``: 1 December through 30 April.
* ``chla_q1_*``: 1 January through 31 March, matching existing cloud/Secchi.

The input is a reef-year CSV with ``year``, ``lon`` and ``lat``. ``ReefID`` or
``LABEL_ID`` is retained as the feature key. NOAA CoralTemp is gap-free and is
sampled at the nearest 0.05-degree cell. IMOS OC3 chlorophyll uses the nearest
valid water cell within 5 km when a reef centroid is masked. Extraction is
cached one year at a time so an interrupted full-GBR run can resume.
'''

from __future__ import annotations

import argparse
import json
import math
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from scipy.stats import kurtosis, skew


SST_URL = 's3://gbr-dms-data-public/noaa-crw-chs-sst/data.zarr'
CHLA_PRODUCTS = {
    'oc3': (
        's3://gbr-dms-data-public/imos-srs-aqua-oc-chla-4km/data.zarr',
        'chl_oc3',
    ),
    'oci': (
        's3://gbr-dms-data-public/imos-srs-aqua-oc-chla-oci/data.zarr',
        'chl_oci',
    ),
}
SST_SOURCE = 'NOAA Coral Reef Watch CoralTemp v3.1 daily 5-km nighttime SST'
CHLA_SOURCE = {
    'oc3': 'IMOS SRS MODIS-Aqua daily chlorophyll-a, OC3 algorithm',
    'oci': 'IMOS SRS MODIS-Aqua daily chlorophyll-a, OCI algorithm',
}


def distribution_summary(values: np.ndarray, prefix: str) -> dict[str, np.ndarray]:
    '''Summarise time x site values with finite-sample shape estimates.'''
    values = np.asarray(values, dtype=float)
    if values.ndim == 1:
        values = values[:, None]
    valid = np.isfinite(values)
    n = valid.sum(axis=0)
    with np.errstate(invalid='ignore'):
        median = np.nanmedian(values, axis=0)
        mean = np.nanmean(values, axis=0)
        sd = np.nanstd(values, axis=0, ddof=1)
        skewness = skew(values, axis=0, bias=False, nan_policy='omit')
        excess_kurtosis = kurtosis(
            values, axis=0, fisher=True, bias=False, nan_policy='omit'
        )
    # Higher moments are unstable with very few observations. SST normally has
    # complete daily coverage; this gate principally protects synthetic tests
    # and any future non-gap-filled product.
    skewness[n < 8] = np.nan
    excess_kurtosis[n < 8] = np.nan
    return {
        f'{prefix}_n': n.astype(int),
        f'{prefix}_median': median,
        f'{prefix}_mean': mean,
        f'{prefix}_sd': sd,
        f'{prefix}_skewness': skewness,
        f'{prefix}_excess_kurtosis': excess_kurtosis,
    }


def chlorophyll_summary(values: np.ndarray, prefix: str) -> dict[str, np.ndarray]:
    '''Summarise positive time x site chlorophyll observations.'''
    values = np.asarray(values, dtype=float)
    if values.ndim == 1:
        values = values[:, None]
    values[values <= 0] = np.nan
    valid = np.isfinite(values)
    n = valid.sum(axis=0)
    with warnings.catch_warnings(), np.errstate(invalid='ignore'):
        warnings.simplefilter('ignore', category=RuntimeWarning)
        median = np.nanmedian(values, axis=0)
        q10 = np.nanquantile(values, 0.10, axis=0)
        q90 = np.nanquantile(values, 0.90, axis=0)
    median[n == 0] = np.nan
    q10[n == 0] = np.nan
    q90[n == 0] = np.nan
    return {
        f'{prefix}_n': n.astype(int),
        f'{prefix}_median': median,
        f'{prefix}_q10': q10,
        f'{prefix}_q90': q90,
    }


def haversine_km(lat0: float, lon0: float, lats, lons):
    lat0_r = np.radians(lat0)
    lats_r = np.radians(lats)
    dlat = np.radians(np.asarray(lats) - lat0)
    dlon = np.radians(np.asarray(lons) - lon0)
    a = np.sin(dlat / 2) ** 2 + np.cos(lat0_r) * np.cos(lats_r) * np.sin(dlon / 2) ** 2
    return 6371.0 * 2 * np.arcsin(np.sqrt(a))


def nearest_valid_pixel(
    valid_grid: np.ndarray,
    grid_lats: np.ndarray,
    grid_lons: np.ndarray,
    reef_lat: float,
    reef_lon: float,
    radius_km: float = 5.0,
) -> tuple[int, int, float] | None:
    '''Find the nearest valid grid cell within a conservative local patch.'''
    lat_i = int(np.abs(grid_lats - reef_lat).argmin())
    lon_i = int(np.abs(grid_lons - reef_lon).argmin())
    if valid_grid[lat_i, lon_i]:
        distance = float(haversine_km(
            reef_lat, reef_lon, grid_lats[lat_i], grid_lons[lon_i]
        ))
        return lat_i, lon_i, distance

    lat_step = max(float(np.median(np.abs(np.diff(grid_lats)))), 1e-6)
    lon_step = max(float(np.median(np.abs(np.diff(grid_lons)))), 1e-6)
    lat_cells = int(math.ceil(radius_km / (111.0 * lat_step))) + 1
    lon_scale = max(111.0 * math.cos(math.radians(reef_lat)), 1e-6)
    lon_cells = int(math.ceil(radius_km / (lon_scale * lon_step))) + 1
    lat_slice = slice(max(0, lat_i - lat_cells), min(len(grid_lats), lat_i + lat_cells + 1))
    lon_slice = slice(max(0, lon_i - lon_cells), min(len(grid_lons), lon_i + lon_cells + 1))
    local_lats, local_lons = np.meshgrid(
        grid_lats[lat_slice], grid_lons[lon_slice], indexing='ij'
    )
    distances = haversine_km(reef_lat, reef_lon, local_lats, local_lons)
    selected = valid_grid[lat_slice, lon_slice] & (distances <= radius_km)
    if not selected.any():
        return None
    flat = np.where(selected, distances, np.inf).argmin()
    local_i, local_j = np.unravel_index(flat, selected.shape)
    return (
        int((lat_slice.start or 0) + local_i),
        int((lon_slice.start or 0) + local_j),
        float(distances[local_i, local_j]),
    )


def extract_sst_window(ds_sst, rows: pd.DataFrame, start: str, end: str, prefix: str):
    lats = xr.DataArray(rows['lat'].to_numpy(), dims='site')
    lons = xr.DataArray(rows['lon'].to_numpy(), dims='site')
    values = (
        ds_sst['analysed_sst']
        .sel(time=slice(start, end))
        .sel(lat=lats, lon=lons, method='nearest')
        .compute()
        .values
    )
    return distribution_summary(values, prefix)


def extract_chla_windows(
    ds_chla,
    variable: str,
    rows: pd.DataFrame,
    windows: dict[str, tuple[str, str]],
    radius_km: float = 5.0,
):
    '''Extract centroid series, with nearest-water fallback shared by windows.'''
    padding = radius_km / 111.0 + 0.02
    lat_min = float(rows['lat'].min() - padding)
    lat_max = float(rows['lat'].max() + padding)
    lon_min = float(rows['lon'].min() - padding)
    lon_max = float(rows['lon'].max() + padding)
    lat_values = ds_chla.latitude.values
    lat_slice = slice(lat_max, lat_min) if lat_values[0] > lat_values[-1] else slice(lat_min, lat_max)

    window_arrays = {}
    valid_any = None
    for prefix, (start, end) in windows.items():
        array = ds_chla[variable].sel(
            time=slice(start, end), latitude=lat_slice,
            longitude=slice(lon_min, lon_max),
        ).compute()
        window_arrays[prefix] = array
        current_valid = np.isfinite(array.values).any(axis=0)
        valid_any = current_valid if valid_any is None else (valid_any | current_valid)

    sample = next(iter(window_arrays.values()))
    grid_lats = sample.latitude.values
    grid_lons = sample.longitude.values
    selected_pixels = [
        nearest_valid_pixel(valid_any, grid_lats, grid_lons, row.lat, row.lon, radius_km)
        for row in rows.itertuples(index=False)
    ]
    output = {
        'chla_match_method': np.array([
            'nearest_valid_water' if pixel is not None else 'missing_no_valid_water'
            for pixel in selected_pixels
        ], dtype=object),
        'chla_grid_distance_km': np.array([
            pixel[2] if pixel is not None else np.nan for pixel in selected_pixels
        ]),
    }
    for prefix, array in window_arrays.items():
        series = np.full((array.sizes['time'], len(rows)), np.nan)
        for site, pixel in enumerate(selected_pixels):
            if pixel is not None:
                series[:, site] = array.values[:, pixel[0], pixel[1]]
        output.update(chlorophyll_summary(series, prefix))
    return output


def extract_year(ds_sst, ds_chla, variable: str, rows: pd.DataFrame, year: int):
    rows = rows.reset_index(drop=True).copy()
    previous = year - 1
    result = rows.copy()
    summaries = {}
    summaries.update(extract_sst_window(
        ds_sst, rows, f'{previous}-11-01', f'{year}-04-30', 'sst_summer'
    ))
    summaries.update(extract_sst_window(
        ds_sst, rows, f'{previous}-01-01', f'{previous}-12-31', 'sst_prevyear'
    ))
    summaries.update(extract_chla_windows(
        ds_chla, variable, rows,
        {
            'chla_wetseason': (f'{previous}-12-01', f'{year}-04-30'),
            'chla_q1': (f'{year}-01-01', f'{year}-03-31'),
        },
    ))
    for name, values in summaries.items():
        result[name] = values
    return result


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input', default='data/processed/environment_screen_grid_validation.csv'
    )
    parser.add_argument(
        '--output', default='data/processed/sst_chla_features_validation.csv'
    )
    parser.add_argument('--product', choices=sorted(CHLA_PRODUCTS), default='oc3')
    parser.add_argument('--years', nargs='*', type=int)
    parser.add_argument('--force', action='store_true')
    return parser.parse_args()


def main():
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)
    rows = pd.read_csv(input_path)
    required = {'year', 'lon', 'lat'}
    missing = required - set(rows.columns)
    if missing:
        raise ValueError(f'Input lacks columns: {sorted(missing)}')
    rows['year'] = pd.to_numeric(rows['year'], errors='raise').astype(int)
    if args.years:
        rows = rows[rows['year'].isin(args.years)].copy()
    if rows.empty:
        raise ValueError('No reef-year rows selected')

    output_path.parent.mkdir(parents=True, exist_ok=True)
    part_dir = output_path.parent / f'{output_path.stem}_parts'
    part_dir.mkdir(parents=True, exist_ok=True)

    print('Opening NOAA CoralTemp and IMOS chlorophyll DMS stores', flush=True)
    ds_sst = xr.open_zarr(SST_URL, storage_options={'anon': True})
    chla_url, chla_variable = CHLA_PRODUCTS[args.product]
    ds_chla = xr.open_zarr(chla_url, storage_options={'anon': True})

    parts = []
    for year in sorted(rows['year'].unique()):
        part_path = part_dir / f'{output_path.stem}_{year}.csv'
        if part_path.exists() and not args.force:
            print(f'Loading cached {year}: {part_path}', flush=True)
            part = pd.read_csv(part_path)
        else:
            selected = rows[rows['year'] == year].copy()
            print(f'Extracting {year}: {len(selected)} reef-years', flush=True)
            part = extract_year(ds_sst, ds_chla, chla_variable, selected, int(year))
            part['sst_source'] = SST_SOURCE
            part['chla_source'] = CHLA_SOURCE[args.product]
            part['chla_product'] = args.product
            part.to_csv(part_path, index=False)
        parts.append(part)

    combined = pd.concat(parts, ignore_index=True)
    combined.to_csv(output_path, index=False)
    manifest = {
        'input': str(input_path),
        'output': str(output_path),
        'rows': int(len(combined)),
        'years': sorted(map(int, combined['year'].unique())),
        'sst_store': SST_URL,
        'chlorophyll_store': chla_url,
        'chlorophyll_variable': chla_variable,
        'chlorophyll_product': args.product,
        'sst_windows': {
            'summer': 'previous-year-11-01/current-year-04-30',
            'previous_year': 'previous-year-01-01/previous-year-12-31',
        },
        'chlorophyll_windows': {
            'wetseason': 'previous-year-12-01/current-year-04-30',
            'q1': 'current-year-01-01/current-year-03-31',
        },
        'chlorophyll_fallback_radius_km': 5.0,
    }
    output_path.with_suffix('.manifest.json').write_text(
        json.dumps(manifest, indent=2), encoding='utf-8'
    )
    print(f'Wrote {len(combined)} rows to {output_path}', flush=True)


if __name__ == '__main__':
    main()
