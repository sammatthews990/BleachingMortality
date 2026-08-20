'''Extract consistent Q1 reef wind metrics from IMOS Metop-B scatterometry.

The user-linked Sentinel-1 SAR product has high spatial resolution but sparse
swath coverage and starts after the 2016 event. Metop-B spans every event in
this analysis. Good-quality calibrated 10-m neutral wind observations within
30 km of each reef are aggregated to daily medians before Q1 summaries are
calculated, so days with multiple swath pixels do not receive extra weight.
'''

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import requests


FILE_ROOT = (
    'https://thredds.aodn.org.au/thredds/fileServer/IMOS/SRS/'
    'Surface-Waves/Wind-Scatterometry-DM00/METOP-B'
)
RAW_DIR = Path('data/raw/imos_metop_b_wind')
OUTPUT = Path('data/processed/imos_metop_b_wind_reef_year.csv')
ENVIRONMENT = Path('data/processed/cheung_recreated_gbr_full.csv')
EVENT_YEARS = (2016, 2017, 2020, 2022, 2024)
SEARCH_RADIUS_KM = 30.0
SOURCE = (
    'IMOS Metop-B calibrated equivalent-neutral 10-m scatterometer wind; '
    'good QC; daily median within 30 km; Q1'
)


def tile_details(lat, lon):
    '''Return the THREDDS folder and file name for one reef coordinate.'''
    latitude_cell = int(np.ceil(-lat))
    longitude_cell = int(np.floor(lon))
    return tile_details_from_cells(latitude_cell, longitude_cell)


def tile_details_from_cells(latitude_cell, longitude_cell):
    '''Return archive names for integer south-latitude/east-longitude cells.'''
    folder_latitude = (latitude_cell // 20 + 1) * 20
    folder_longitude = longitude_cell // 20 * 20
    folder = f'{folder_latitude:03d}S_{folder_longitude:03d}E'
    filename = (
        'IMOS_SRS-Surface-Waves_M_Wind-METOP-B_FV02_'
        f'{latitude_cell:03d}S-{longitude_cell:03d}E-DM00.nc'
    )
    return folder, filename


def haversine_km(lat, lon, observation_lat, observation_lon):
    '''Great-circle distance from one reef to arrays of observations.'''
    reef_latitude = np.radians(lat)
    observation_latitude = np.radians(observation_lat)
    latitude_difference = observation_latitude - reef_latitude
    longitude_difference = np.radians(observation_lon - lon)
    a = (
        np.sin(latitude_difference / 2) ** 2
        + np.cos(reef_latitude) * np.cos(observation_latitude)
        * np.sin(longitude_difference / 2) ** 2
    )
    return 6371.0 * 2 * np.arcsin(np.sqrt(a))


def longest_true_run(values):
    '''Length of the longest consecutive True run; missing dates are False.'''
    longest = 0
    current = 0
    for value in values:
        if value:
            current += 1
            longest = max(longest, current)
        else:
            current = 0
    return longest


def download_tile(folder, filename):
    '''Download one small official tile unless a non-empty cache exists.'''
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    destination = RAW_DIR / filename
    if destination.exists() and destination.stat().st_size > 0:
        return destination
    url = f'{FILE_ROOT}/{folder}/{filename}'
    partial = destination.with_suffix('.nc.partial')
    last_error = None
    for attempt in range(4):
        try:
            response = requests.get(url, timeout=120)
            if response.status_code == 404:
                if partial.exists():
                    partial.unlink()
                return None
            response.raise_for_status()
            partial.write_bytes(response.content)
            partial.replace(destination)
            return destination
        except requests.RequestException as error:
            last_error = error
            if partial.exists():
                partial.unlink()
            if attempt == 3:
                raise
    raise last_error


def summarise_daily_wind(daily, requested_dates):
    '''Return coverage and calm-condition summaries from a daily wind series.'''
    daily = daily.reindex(requested_dates)
    observed = daily.dropna()
    n_days = len(observed)
    summary = {
        'wind_n_days': n_days,
        'wind_day_coverage': n_days / len(requested_dates),
        'wind_mean': np.nan,
        'wind_q10': np.nan,
        'wind_min': np.nan,
        'wind_days_below_3': np.nan,
        'wind_fraction_below_3': np.nan,
        'wind_longest_spell_below_3': np.nan,
        'wind_days_below_5': np.nan,
        'wind_fraction_below_5': np.nan,
        'wind_longest_spell_below_5': np.nan,
    }
    if n_days == 0:
        return summary
    summary.update({
        'wind_mean': float(observed.mean()),
        'wind_q10': float(observed.quantile(0.10)),
        'wind_min': float(observed.min()),
        'wind_days_below_3': int((observed < 3).sum()),
        'wind_fraction_below_3': float((observed < 3).mean()),
        'wind_longest_spell_below_3': longest_true_run(
            daily.lt(3).fillna(False).to_numpy()
        ),
        'wind_days_below_5': int((observed < 5).sum()),
        'wind_fraction_below_5': float((observed < 5).mean()),
        'wind_longest_spell_below_5': longest_true_run(
            daily.lt(5).fillna(False).to_numpy()
        ),
    })
    return summary


def main():
    environment = pd.read_csv(ENVIRONMENT)
    reefs = environment.loc[
        environment['year'].eq(2024),
        ['LABEL_ID', 'LOC_NAME_S', 'lon', 'lat']
    ].copy()
    if len(reefs) != 7063:
        raise ValueError(f'Expected 7,063 registry features, found {len(reefs)}')
    reefs['wind_latitude_cell'] = np.ceil(-reefs['lat']).astype(int)
    reefs['wind_longitude_cell'] = np.floor(reefs['lon']).astype(int)
    home_cells = reefs[
        ['wind_latitude_cell', 'wind_longitude_cell']
    ].drop_duplicates()
    tile_records = []
    for home in home_cells.itertuples(index=False):
        for latitude_cell in range(
            home.wind_latitude_cell - 1, home.wind_latitude_cell + 2
        ):
            for longitude_cell in range(
                home.wind_longitude_cell - 1, home.wind_longitude_cell + 2
            ):
                folder, filename = tile_details_from_cells(
                    latitude_cell, longitude_cell
                )
                tile_records.append({
                    'latitude_cell': latitude_cell,
                    'longitude_cell': longitude_cell,
                    'folder': folder,
                    'filename': filename,
                })
    tiles = pd.DataFrame(tile_records).drop_duplicates(
        ['latitude_cell', 'longitude_cell']
    )
    print(f'Downloading or checking {len(tiles)} neighbouring Metop-B tiles')
    tile_paths = {}
    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = {
            executor.submit(
                download_tile, row.folder, row.filename
            ): (row.latitude_cell, row.longitude_cell)
            for row in tiles.itertuples(index=False)
        }
        completed = 0
        for future in as_completed(futures):
            tile_key = futures[future]
            tile_paths[tile_key] = future.result()
            completed += 1
            if completed % 10 == 0 or completed == len(futures):
                print(f'  {completed}/{len(futures)} tiles ready')

    available_tiles = {
        key: path for key, path in tile_paths.items() if path is not None
    }
    print(
        f'{len(available_tiles)} of {len(tile_paths)} neighbouring tiles '
        'contain archived ocean observations'
    )

    tile_event_observations = {}
    for tile_number, (tile_key, tile_path) in enumerate(
        available_tiles.items(), 1
    ):
        with h5py.File(tile_path, 'r') as source:
            time_values = source['TIME'][:]
            times = pd.Timestamp('1985-01-01') + pd.to_timedelta(
                time_values, unit='D'
            )
            latitudes = source['LATITUDE'][:]
            longitudes = source['LONGITUDE'][:]
            packed_speed = source['WSPD_CAL'][:]
            quality = source['WSPD_CAL_quality_control'][:]
            scale = float(source['WSPD_CAL'].attrs['scale_factor'][0])
        for year in EVENT_YEARS:
            requested_dates = pd.date_range(
                f'{year}-01-01', f'{year}-03-31', freq='D'
            )
            selected_times = (
                (times >= requested_dates[0])
                & (times < requested_dates[-1] + pd.Timedelta(days=1))
            )
            selected_rows = np.where(selected_times)[0]
            observation_time = np.repeat(
                times[selected_rows].values, latitudes.shape[1]
            )
            observation_lat = latitudes[selected_rows].reshape(-1)
            observation_lon = longitudes[selected_rows].reshape(-1)
            observation_packed = packed_speed[selected_rows].reshape(-1)
            observation_quality = quality[selected_rows].reshape(-1)
            good = (
                (observation_quality == 1)
                & (observation_packed != -32768)
                & np.isfinite(observation_lat)
                & np.isfinite(observation_lon)
            )
            observation_time = pd.to_datetime(observation_time[good])
            observation_lat = observation_lat[good]
            observation_lon = observation_lon[good]
            observation_speed = observation_packed[good].astype(float) * scale
            tile_event_observations[(*tile_key, year)] = (
                observation_time.to_numpy(),
                observation_lat,
                observation_lon,
                observation_speed,
            )
        if tile_number % 20 == 0 or tile_number == len(available_tiles):
            print(
                f'  {tile_number}/{len(available_tiles)} tile files loaded'
            )

    output_rows = []
    grouped_reefs = reefs.groupby(
        ['wind_latitude_cell', 'wind_longitude_cell'], sort=False
    )
    for group_number, (home_cell, tile_reefs) in enumerate(grouped_reefs, 1):
        neighbour_cells = [
            (latitude_cell, longitude_cell)
            for latitude_cell in range(home_cell[0] - 1, home_cell[0] + 2)
            for longitude_cell in range(home_cell[1] - 1, home_cell[1] + 2)
        ]
        home_folder, home_filename = tile_details_from_cells(*home_cell)
        for year in EVENT_YEARS:
            requested_dates = pd.date_range(
                f'{year}-01-01', f'{year}-03-31', freq='D'
            )
            observation_chunks = [
                tile_event_observations[(*cell, year)]
                for cell in neighbour_cells
                if (*cell, year) in tile_event_observations
            ]
            if observation_chunks:
                observation_time = pd.to_datetime(np.concatenate([
                    chunk[0] for chunk in observation_chunks
                ]))
                observation_lat = np.concatenate([
                    chunk[1] for chunk in observation_chunks
                ])
                observation_lon = np.concatenate([
                    chunk[2] for chunk in observation_chunks
                ])
                observation_speed = np.concatenate([
                    chunk[3] for chunk in observation_chunks
                ])
            else:
                observation_time = pd.DatetimeIndex([])
                observation_lat = np.array([], dtype=float)
                observation_lon = np.array([], dtype=float)
                observation_speed = np.array([], dtype=float)
            for reef in tile_reefs.itertuples(index=False):
                distances = haversine_km(
                    reef.lat, reef.lon, observation_lat, observation_lon
                )
                within_radius = distances <= SEARCH_RADIUS_KM
                if within_radius.any():
                    daily = pd.Series(
                        observation_speed[within_radius],
                        index=observation_time[within_radius].normalize()
                    ).groupby(level=0).median()
                    nearest_distance = float(distances[within_radius].min())
                    n_observations = int(within_radius.sum())
                else:
                    daily = pd.Series(dtype=float)
                    nearest_distance = np.nan
                    n_observations = 0
                summary = summarise_daily_wind(daily, requested_dates)
                output_rows.append({
                    'LABEL_ID': reef.LABEL_ID,
                    'LOC_NAME_S': reef.LOC_NAME_S,
                    'lon': reef.lon,
                    'lat': reef.lat,
                    'year': year,
                    'wind_source': SOURCE,
                    'wind_search_radius_km': SEARCH_RADIUS_KM,
                    'wind_nearest_observation_km': nearest_distance,
                    'wind_n_observations': n_observations,
                    'wind_requested_days': len(requested_dates),
                    'wind_home_tile_file': home_filename,
                    'wind_home_tile_available': home_cell in available_tiles,
                    'wind_neighbour_tiles_available': len(observation_chunks),
                    **summary,
                })
        if group_number % 10 == 0 or group_number == len(grouped_reefs):
            print(
                f'  {group_number}/{len(grouped_reefs)} reef groups summarised'
            )

    output = pd.DataFrame(output_rows)
    if len(output) != len(reefs) * len(EVENT_YEARS):
        raise RuntimeError('Unexpected number of reef-year wind rows')
    if output.duplicated(
        ['LABEL_ID', 'LOC_NAME_S', 'lon', 'lat', 'year']
    ).any():
        raise RuntimeError('Duplicate reef-feature-year wind rows')
    output['wind_coverage_adequate'] = output['wind_day_coverage'] >= 0.70
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(OUTPUT, index=False)

    coverage = output.groupby('year').agg(
        reefs=('LABEL_ID', 'size'),
        median_days=('wind_n_days', 'median'),
        q10_days=('wind_n_days', lambda values: values.quantile(0.10)),
        adequate=('wind_coverage_adequate', 'sum'),
        median_wind=('wind_mean', 'median'),
        median_calm_fraction=('wind_fraction_below_3', 'median'),
    )
    print('\n', coverage)
    print(f'Wrote {len(output)} rows to {OUTPUT}')


if __name__ == '__main__':
    main()
