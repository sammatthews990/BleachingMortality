'''Extract compact ERA5 rainfall, wind, and air-temperature event summaries.

The preferred research route is the Copernicus Climate Data Store (CDS).
Because CDS retrieval requires a user account, accepted licence, and API key,
the executable default uses Open-Meteo's documented ERA5-only mirror. It does
not use Open-Meteo's changing 'best match' product.

For each reef, the script samples the nearest 0.25-degree ERA5 cell and the
nearest 0.25-degree land-cell centre. The latter is labelled a coastal-rainfall
proxy: it is useful for testing but is not catchment runoff or river discharge.
'''

from __future__ import annotations

import argparse
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from shapely import contains_xy


EVENT_YEARS = (2016, 2017, 2020, 2022, 2024)
ERA5_STEP = 0.25
OPEN_METEO_URL = 'https://archive-api.open-meteo.com/v1/archive'
NATURAL_EARTH_URL = (
    'https://naturalearth.s3.amazonaws.com/10m_physical/ne_10m_land.zip'
)
NATURAL_EARTH_FILE = Path(
    'data/raw/natural_earth/ne_10m_land.zip'
)
ENVIRONMENT_FILE = Path('data/processed/cheung_recreated_gbr_full.csv')
OUTPUT_FILE = Path('data/processed/era5_weather_reef_year.csv')
PART_DIR = Path('data/processed/era5_weather_parts')
SOURCE = (
    'ECMWF ERA5 0.25-degree hourly reanalysis accessed through the '
    'Open-Meteo Historical Weather API with models=era5'
)
SOURCE_DOI = '10.24381/cds.adbb2d47'


def nearest_grid(values, step=ERA5_STEP):
    '''Round coordinates to the nearest regular grid without banker's ties.'''
    values = np.asarray(values, dtype=float)
    return np.floor(values / step + 0.5) * step


def longest_true_run(values):
    '''Length of the longest consecutive True run.'''
    longest = 0
    current = 0
    for value in np.asarray(values, dtype=bool):
        current = current + 1 if value else 0
        longest = max(longest, current)
    return int(longest)


def summarise_hourly_weather(hourly):
    '''Summarise one December--March hourly ERA5 time series.'''
    frame = hourly.copy()
    frame['time'] = pd.to_datetime(frame['time'], utc=True)
    frame = frame.set_index('time').sort_index()
    event_year = int(frame.index.max().year)
    q1 = frame.loc[
        (frame.index >= pd.Timestamp(f'{event_year}-01-01', tz='UTC'))
        & (frame.index <= pd.Timestamp(f'{event_year}-03-31 23:59', tz='UTC'))
    ]
    december = frame.loc[frame.index.year == event_year - 1]

    daily = frame.resample('1D').agg(
        precipitation=('precipitation', 'sum'),
        wind_speed_10m=('wind_speed_10m', 'mean'),
        temperature_2m_mean=('temperature_2m', 'mean'),
        temperature_2m_max=('temperature_2m', 'max'),
    )
    q1_daily = daily.loc[daily.index.year == event_year]
    rolling_7 = q1_daily['precipitation'].rolling(7, min_periods=1).sum()
    rolling_30 = daily['precipitation'].rolling(30, min_periods=1).sum()

    valid_wind = q1_daily['wind_speed_10m'].dropna()
    summary = {
        'period_start': frame.index.min().date().isoformat(),
        'period_end': frame.index.max().date().isoformat(),
        'n_hours': int(frame['wind_speed_10m'].notna().sum()),
        'q1_n_days': int(q1_daily['wind_speed_10m'].notna().sum()),
        'rain_december_total': float(december['precipitation'].sum()),
        'rain_q1_total': float(q1['precipitation'].sum()),
        'rain_q1_max_1day': float(q1_daily['precipitation'].max()),
        'rain_q1_max_7day': float(rolling_7.max()),
        'rain_dec_mar_max_30day': float(rolling_30.max()),
        'rain_q1_days_ge_20mm': int((q1_daily['precipitation'] >= 20).sum()),
        'wind_q1_mean': float(valid_wind.mean()),
        'wind_q1_q10': float(valid_wind.quantile(0.10)),
        'wind_q1_min': float(valid_wind.min()),
        'wind_q1_fraction_below_3': float((valid_wind < 3).mean()),
        'wind_q1_fraction_below_5': float((valid_wind < 5).mean()),
        'wind_q1_longest_spell_below_3': longest_true_run(valid_wind < 3),
        'wind_q1_longest_spell_below_5': longest_true_run(valid_wind < 5),
        'temperature_q1_mean': float(q1['temperature_2m'].mean()),
        'temperature_q1_max': float(q1['temperature_2m'].max()),
    }
    return summary


def unit_sphere_coordinates(latitude, longitude):
    lat = np.deg2rad(np.asarray(latitude, dtype=float))
    lon = np.deg2rad(np.asarray(longitude, dtype=float))
    return np.column_stack(
        [np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)]
    )


def chord_to_km(chord_distance):
    angle = 2 * np.arcsin(np.minimum(np.asarray(chord_distance) / 2, 1))
    return angle * 6371.0088


def build_land_grid(reef_grid):
    '''Return ERA5 land-cell centres and nearest-land mapping for reef cells.'''
    if not NATURAL_EARTH_FILE.exists():
        NATURAL_EARTH_FILE.parent.mkdir(parents=True, exist_ok=True)
        print(f'Downloading Natural Earth land mask: {NATURAL_EARTH_URL}')
        urllib.request.urlretrieve(NATURAL_EARTH_URL, NATURAL_EARTH_FILE)
    land = gpd.read_file(f'zip://{NATURAL_EARTH_FILE.resolve()}')
    land_geometry = land.geometry.union_all()

    latitude = np.arange(-27.0, -7.75, ERA5_STEP)
    longitude = np.arange(139.0, 155.25, ERA5_STEP)
    lon_grid, lat_grid = np.meshgrid(longitude, latitude)
    on_land = contains_xy(land_geometry, lon_grid.ravel(), lat_grid.ravel())
    land_grid = pd.DataFrame(
        {
            'grid_lat': lat_grid.ravel()[on_land],
            'grid_lon': lon_grid.ravel()[on_land],
        }
    )
    if land_grid.empty:
        raise RuntimeError('Natural Earth mask produced no land cells near GBR')

    tree = cKDTree(
        unit_sphere_coordinates(land_grid['grid_lat'], land_grid['grid_lon'])
    )
    distance, index = tree.query(
        unit_sphere_coordinates(reef_grid['grid_lat'], reef_grid['grid_lon'])
    )
    mapping = reef_grid.copy()
    mapping['coastal_grid_lat'] = land_grid.iloc[index]['grid_lat'].to_numpy()
    mapping['coastal_grid_lon'] = land_grid.iloc[index]['grid_lon'].to_numpy()
    mapping['coastal_grid_distance_km'] = chord_to_km(distance)
    return mapping


def fetch_open_meteo_batch(coordinates, event_year, attempts=5):
    '''Fetch a batch of exact ERA5 grid-cell time series.'''
    start_date = f'{event_year - 1}-12-01'
    end_date = f'{event_year}-03-31'
    params = {
        'latitude': ','.join(f'{lat:.2f}' for lat, _ in coordinates),
        'longitude': ','.join(f'{lon:.2f}' for _, lon in coordinates),
        'start_date': start_date,
        'end_date': end_date,
        'hourly': 'precipitation,wind_speed_10m,temperature_2m',
        'models': 'era5',
        'cell_selection': 'nearest',
        'wind_speed_unit': 'ms',
        'precipitation_unit': 'mm',
        'timezone': 'GMT',
    }
    url = OPEN_METEO_URL + '?' + urllib.parse.urlencode(params)
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(url, timeout=180) as response:
                result = json.load(response)
            if isinstance(result, dict):
                result = [result]
            if len(result) != len(coordinates):
                raise RuntimeError(
                    f'ERA5 API returned {len(result)} locations for '
                    f'{len(coordinates)} requests'
                )
            return result
        except (urllib.error.URLError, TimeoutError) as error:
            if attempt == attempts - 1:
                raise
            wait_seconds = 2 ** attempt
            print(f'  API retry after {error}: {wait_seconds}s')
            time.sleep(wait_seconds)
    raise RuntimeError('ERA5 request exhausted retries')


def extract_grid_year(grid, event_year, batch_size=20):
    records = []
    coordinates = list(zip(grid['grid_lat'], grid['grid_lon']))
    for start in range(0, len(coordinates), batch_size):
        batch = coordinates[start:start + batch_size]
        print(
            f'  {event_year} ERA5 cells {start + 1}--'
            f'{min(start + len(batch), len(coordinates))}/{len(coordinates)}',
            flush=True,
        )
        responses = fetch_open_meteo_batch(batch, event_year)
        for requested, response in zip(batch, responses):
            hourly = pd.DataFrame(response['hourly'])
            summary = summarise_hourly_weather(hourly)
            summary.update(
                {
                    'year': event_year,
                    'grid_lat': requested[0],
                    'grid_lon': requested[1],
                    'returned_lat': response['latitude'],
                    'returned_lon': response['longitude'],
                    'source': SOURCE,
                    'source_doi': SOURCE_DOI,
                    'access_route': 'open-meteo-era5-mirror',
                }
            )
            records.append(summary)
    return pd.DataFrame(records)


def prefix_summary(summary, prefix):
    identity = {'year', 'grid_lat', 'grid_lon'}
    return summary.rename(
        columns={column: f'{prefix}{column}' for column in summary if column not in identity}
    )


def main(batch_size=20):
    environment = pd.read_csv(ENVIRONMENT_FILE)
    reef_columns = ['LABEL_ID', 'LOC_NAME_S', 'lon', 'lat']
    reefs = environment[reef_columns].drop_duplicates().reset_index(drop=True)
    if len(reefs) != 7063:
        raise ValueError(f'Expected 7,063 reef features, found {len(reefs)}')
    reefs['grid_lat'] = nearest_grid(reefs['lat'])
    reefs['grid_lon'] = nearest_grid(reefs['lon'])
    reef_grid = reefs[['grid_lat', 'grid_lon']].drop_duplicates().reset_index(drop=True)
    mapping = build_land_grid(reef_grid)
    reefs = reefs.merge(mapping, on=['grid_lat', 'grid_lon'], validate='many_to_one')

    all_grid = pd.concat(
        [
            reef_grid,
            mapping[['coastal_grid_lat', 'coastal_grid_lon']].rename(
                columns={
                    'coastal_grid_lat': 'grid_lat',
                    'coastal_grid_lon': 'grid_lon',
                }
            ),
        ],
        ignore_index=True,
    ).drop_duplicates().sort_values(['grid_lat', 'grid_lon']).reset_index(drop=True)
    print(
        f'{len(reefs)} reef features, {len(reef_grid)} reef ERA5 cells, '
        f'{len(all_grid)} total reef/coastal ERA5 cells'
    )

    PART_DIR.mkdir(parents=True, exist_ok=True)
    grid_years = []
    for year in EVENT_YEARS:
        part_path = PART_DIR / f'era5_grid_{year}.csv'
        if part_path.exists():
            print(f'Loading cached ERA5 summaries: {part_path}')
            grid_years.append(pd.read_csv(part_path))
            continue
        extracted = extract_grid_year(all_grid, year, batch_size=batch_size)
        extracted.to_csv(part_path, index=False)
        grid_years.append(extracted)
    grid_summary = pd.concat(grid_years, ignore_index=True)

    outputs = []
    for year in EVENT_YEARS:
        year_summary = grid_summary.loc[grid_summary['year'].eq(year)].copy()
        reef_summary = prefix_summary(year_summary, 'era5_reef_')
        coastal_summary = prefix_summary(year_summary, 'era5_coastal_').rename(
            columns={
                'grid_lat': 'coastal_grid_lat',
                'grid_lon': 'coastal_grid_lon',
            }
        )
        year_output = reefs.merge(
            reef_summary, on=['grid_lat', 'grid_lon'], validate='many_to_one'
        ).merge(
            coastal_summary,
            on=['coastal_grid_lat', 'coastal_grid_lon'],
            validate='many_to_one',
        )
        year_output['year'] = year
        outputs.append(year_output)
    output = pd.concat(outputs, ignore_index=True)
    if output.duplicated(reef_columns + ['year']).any():
        raise RuntimeError('Duplicate reef-feature-year rows in ERA5 output')
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(OUTPUT_FILE, index=False)
    coverage = output.groupby('year').agg(
        rows=('LABEL_ID', 'size'),
        q1_days=('era5_reef_q1_n_days', 'min'),
        reef_rain=('era5_reef_rain_q1_total', 'mean'),
        coastal_rain=('era5_coastal_rain_q1_total', 'mean'),
        wind=('era5_reef_wind_q1_mean', 'mean'),
    )
    print(coverage)
    print(f'Wrote {len(output)} rows to {OUTPUT_FILE}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--batch-size', type=int, default=20)
    arguments = parser.parse_args()
    main(batch_size=arguments.batch_size)
