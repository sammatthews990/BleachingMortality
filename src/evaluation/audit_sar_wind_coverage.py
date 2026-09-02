'''Audit Q1 GBR coverage of the IMOS Sentinel-1 gridded SAR wind product.

Daily files are nominal 0.01-degree grids, but only pixels intersected by an
available Sentinel swath contain observations. This script samples the GBR at
0.1-degree spacing to determine whether the product is dense enough for calm-
duration metrics. It is a coverage audit, not the final reef extraction.
'''

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import re
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import requests


THREDDS = 'https://thredds.aodn.org.au/thredds'
CATALOG_ROOT = (
    f'{THREDDS}/catalog/IMOS/SRS/Surface-Waves/SAR_Wind/'
    'DELAYED/GRIDDED/SENTINEL-1'
)
DODS_ROOT = (
    f'{THREDDS}/dodsC/IMOS/SRS/Surface-Waves/SAR_Wind/'
    'DELAYED/GRIDDED/SENTINEL-1'
)
YEAR = 2024
MONTHS = (1, 2, 3)

# The source grid starts at -48.995, 95.005 with 0.01-degree spacing.
# Indices below cover approximately 24--9 S and 142--154 E. Stride 10 gives a
# 0.1-degree audit grid while preserving the satellite swath geometry.
LAT_INDICES = np.arange(2500, 4001, 10)
LON_INDICES = np.arange(4700, 5901, 10)
QUERY = (
    '?WSPD_CAL[0:1:0]'
    f'[{LAT_INDICES[0]}:10:{LAT_INDICES[-1]}]'
    f'[{LON_INDICES[0]}:10:{LON_INDICES[-1]}]'
)
EXPECTED_VALUES = len(LAT_INDICES) * len(LON_INDICES)
ROW_PATTERN = re.compile(r'^\[0\]\[\d+\],\s*(.*)$')


def catalogue_files(year, month):
    '''Return NetCDF URL paths from one official monthly THREDDS catalogue.'''
    url = f'{CATALOG_ROOT}/{year}/{month:02d}/catalog.xml'
    response = requests.get(url, timeout=60)
    response.raise_for_status()
    root = ET.fromstring(response.content)
    namespace = {'t': 'http://www.unidata.ucar.edu/namespaces/thredds/InvCatalog/v1.0'}
    files = []
    for dataset in root.findall('.//t:dataset[@urlPath]', namespace):
        path = dataset.attrib['urlPath']
        if path.endswith('.nc'):
            files.append(path)
    return sorted(files)


def read_daily_grid(path):
    '''Read the constrained calibrated-speed grid and return unpacked m/s.'''
    url = f'{THREDDS}/dodsC/{path}.ascii{QUERY}'
    last_error = None
    for attempt in range(4):
        try:
            response = requests.get(url, timeout=90)
            response.raise_for_status()
            values = []
            for line in response.text.splitlines():
                match = ROW_PATTERN.match(line)
                if match:
                    values.extend(int(value.strip()) for value in match.group(1).split(','))
            if len(values) != EXPECTED_VALUES:
                raise ValueError(
                    f'Expected {EXPECTED_VALUES} values, received {len(values)}'
                )
            packed = np.asarray(values, dtype=np.int16)
            unpacked = packed.astype(float) * 0.001 + 30.0
            unpacked[packed == -32768] = np.nan
            return unpacked
        except (requests.RequestException, ValueError) as error:
            last_error = error
            if attempt == 3:
                raise
    raise last_error


def main():
    paths = []
    for month in MONTHS:
        paths.extend(catalogue_files(YEAR, month))
    if len(paths) not in (90, 91):
        raise RuntimeError(f'Expected a full Q1 daily catalogue, found {len(paths)} files')

    daily_values = {}
    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = {executor.submit(read_daily_grid, path): path for path in paths}
        completed = 0
        for future in as_completed(futures):
            path = futures[future]
            date_match = re.search(r'_M_(\d{8})_', path)
            if date_match is None:
                raise ValueError(f'No date in file path: {path}')
            date = pd.to_datetime(date_match.group(1), format='%Y%m%d')
            daily_values[date] = future.result()
            completed += 1
            if completed % 15 == 0 or completed == len(paths):
                print(f'  {completed}/{len(paths)} daily grids audited')

    dates = sorted(daily_values)
    matrix = np.vstack([daily_values[date] for date in dates])
    valid = np.isfinite(matrix)
    daily_coverage = pd.DataFrame({
        'date': dates,
        'audit_grid_cells': matrix.shape[1],
        'valid_grid_cells': valid.sum(axis=1),
        'valid_grid_fraction': valid.mean(axis=1),
        'median_wind_speed': np.nanmedian(matrix, axis=1),
        'cells_below_3ms': np.sum(matrix < 3, axis=1),
        'cells_below_5ms': np.sum(matrix < 5, axis=1),
    })

    lats = -48.995 + LAT_INDICES * 0.01
    lons = 95.005 + LON_INDICES * 0.01
    grid_lats, grid_lons = np.meshgrid(lats, lons, indexing='ij')
    grid_summary = pd.DataFrame({
        'lat': grid_lats.ravel(),
        'lon': grid_lons.ravel(),
        'q1_days': len(dates),
        'observed_days': valid.sum(axis=0),
        'observed_day_fraction': valid.mean(axis=0),
        'wind_mean_observed_swaths': np.nanmean(matrix, axis=0),
        'wind_q10_observed_swaths': np.nanquantile(matrix, 0.10, axis=0),
        'observed_fraction_below_3ms': np.nanmean(matrix < 3, axis=0),
        'observed_fraction_below_5ms': np.nanmean(matrix < 5, axis=0),
    })
    grid_summary.loc[grid_summary['observed_days'] == 0, [
        'wind_mean_observed_swaths', 'wind_q10_observed_swaths',
        'observed_fraction_below_3ms', 'observed_fraction_below_5ms'
    ]] = np.nan

    output_dir = Path('output/wind_audit')
    output_dir.mkdir(parents=True, exist_ok=True)
    daily_coverage.to_csv(
        output_dir / 'sar_wind_2024_q1_daily_coverage.csv', index=False
    )
    grid_summary.to_csv(
        output_dir / 'sar_wind_2024_q1_grid_coverage.csv', index=False
    )

    registry = pd.read_csv('data/processed/cheung_recreated_gbr_full.csv')
    registry = registry.loc[
        registry['year'].eq(YEAR),
        ['LABEL_ID', 'LOC_NAME_S', 'lon', 'lat']
    ].copy()
    nearest_lat = np.clip(
        np.rint((registry['lat'].to_numpy() - lats[0]) / 0.1).astype(int),
        0, len(lats) - 1
    )
    nearest_lon = np.clip(
        np.rint((registry['lon'].to_numpy() - lons[0]) / 0.1).astype(int),
        0, len(lons) - 1
    )
    flat_index = nearest_lat * len(lons) + nearest_lon
    reef_coverage = pd.concat(
        [
            registry.reset_index(drop=True),
            grid_summary.iloc[flat_index].reset_index(drop=True).add_prefix(
                'audit_'
            ),
        ],
        axis=1,
    )
    reef_coverage.to_csv(
        output_dir / 'sar_wind_2024_q1_reef_coverage.csv', index=False
    )
    print('\nDaily GBR coverage fraction:')
    print(daily_coverage['valid_grid_fraction'].describe())
    print('\nObserved days per 0.1-degree grid cell:')
    print(grid_summary['observed_days'].describe())
    print(
        'Cells with at least 10 observed Q1 days:',
        int((grid_summary['observed_days'] >= 10).sum()),
        'of', len(grid_summary)
    )
    print('\nObserved days at registry features (nearest audit cell):')
    print(reef_coverage['audit_observed_days'].describe())
    print(
        'Reef features with at least 10 observed Q1 days:',
        int((reef_coverage['audit_observed_days'] >= 10).sum()),
        'of', len(reef_coverage)
    )


if __name__ == '__main__':
    main()
