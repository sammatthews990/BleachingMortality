'''Extract continuous Q1 eReefs wind metrics at all GBR registry features.

The GBR4 v2 daily store contains the atmospheric wind forcing used by eReefs.
It provides complete Q1 coverage for historical events but stops on 17 January
2024. The partial 2024 values are retained and never labelled full-quarter
doldrum exposure.
'''

from pathlib import Path

import numpy as np
import pandas as pd
import s3fs
import zarr


STORE = 'gbr-dms-data-public/aims-ereefs-agg-hydrodynamic-4km-daily/data.zarr'
EVENT_YEARS = (2016, 2017, 2020, 2022, 2024)
MAPPING_FILE = Path('data/processed/ereefs_salinity_imos_kd490_reef_year.csv')
OUTPUT_FILE = Path('data/processed/ereefs_wind_reef_year.csv')
PART_DIR = Path('data/processed/ereefs_wind_parts')
SOURCE = (
    'AIMS eReefs GBR4 hydrodynamic v2 daily mean_wspeed atmospheric forcing; Q1'
)


def summarise_wind_matrix(daily_wind):
    '''Summarise rows of days by columns of unique eReefs grid cells.'''
    wind = np.asarray(daily_wind, dtype=float)
    valid = np.isfinite(wind)
    n_days = valid.sum(axis=0)
    has_data = n_days > 0
    summary = {'ereefs_wind_n_days': n_days}
    for name in ('mean', 'q10', 'min'):
        summary[f'ereefs_wind_{name}'] = np.full(n_days.shape, np.nan)
    observed = wind[:, has_data]
    summary['ereefs_wind_mean'][has_data] = np.nanmean(observed, axis=0)
    summary['ereefs_wind_q10'][has_data] = np.nanquantile(
        observed, 0.10, axis=0
    )
    summary['ereefs_wind_min'][has_data] = np.nanmin(observed, axis=0)

    for threshold in (3.0, 5.0):
        label = int(threshold)
        below = valid & (wind < threshold)
        days_below = below.sum(axis=0)
        fraction = np.full(n_days.shape, np.nan)
        fraction[has_data] = days_below[has_data] / n_days[has_data]
        current_run = np.zeros(n_days.shape, dtype=int)
        longest_run = np.zeros(n_days.shape, dtype=int)
        for day_below in below:
            current_run = np.where(day_below, current_run + 1, 0)
            longest_run = np.maximum(longest_run, current_run)
        longest_run = longest_run.astype(float)
        longest_run[~has_data] = np.nan
        summary[f'ereefs_wind_days_below_{label}'] = days_below.astype(float)
        summary[f'ereefs_wind_days_below_{label}'][~has_data] = np.nan
        summary[f'ereefs_wind_fraction_below_{label}'] = fraction
        summary[f'ereefs_wind_longest_spell_below_{label}'] = longest_run
    return summary


def main():
    mapping = pd.read_csv(MAPPING_FILE)
    reef_columns = ['LABEL_ID', 'LOC_NAME_S', 'lon', 'lat']
    reefs = mapping.loc[
        mapping['year'].eq(2024),
        reef_columns + [
            'ereefs_latitude', 'ereefs_longitude',
            'ereefs_lat_index', 'ereefs_lon_index'
        ]
    ].copy()
    if len(reefs) != 7063:
        raise ValueError(f'Expected 7,063 mapped features, found {len(reefs)}')

    unique_cells = reefs[[
        'ereefs_lat_index', 'ereefs_lon_index',
        'ereefs_latitude', 'ereefs_longitude'
    ]].drop_duplicates().reset_index(drop=True)
    cell_key = list(zip(
        unique_cells['ereefs_lat_index'], unique_cells['ereefs_lon_index']
    ))
    cell_lookup = {key: index for index, key in enumerate(cell_key)}
    reef_cell_index = np.asarray([
        cell_lookup[key] for key in zip(
            reefs['ereefs_lat_index'], reefs['ereefs_lon_index']
        )
    ])

    PART_DIR.mkdir(parents=True, exist_ok=True)
    fs = s3fs.S3FileSystem(anon=True)
    store = s3fs.S3Map(STORE, s3=fs)
    ereefs = zarr.open_group(store, mode='r')
    times = pd.to_datetime(
        ereefs['time'][:], unit='D', origin='1990-01-01'
    )
    row_indices = unique_cells['ereefs_lat_index'].to_numpy(dtype=int)
    column_indices = unique_cells['ereefs_lon_index'].to_numpy(dtype=int)

    outputs = []
    for year in EVENT_YEARS:
        part_path = PART_DIR / f'ereefs_wind_{year}.csv'
        if part_path.exists():
            print(f'Loading cached eReefs wind for {year}: {part_path}')
            outputs.append(pd.read_csv(part_path))
            continue

        requested_start = pd.Timestamp(f'{year}-01-01')
        requested_end = pd.Timestamp(f'{year}-03-31')
        available = (
            (times >= requested_start)
            & (times <= requested_end)
        )
        time_indices = np.where(available)[0]
        if len(time_indices) == 0:
            raise RuntimeError(f'No eReefs wind days available for {year}')
        print(
            f'Extracting {year}: {times[time_indices[0]].date()} to '
            f'{times[time_indices[-1]].date()} ({len(time_indices)} days)'
        )
        daily = np.empty((len(time_indices), len(unique_cells)), dtype=float)
        for day_number, time_index in enumerate(time_indices, 1):
            surface = np.asarray(
                ereefs['mean_wspeed'][time_index, :, :], dtype=float
            )
            daily[day_number - 1, :] = surface[
                row_indices, column_indices
            ]
            if day_number % 15 == 0 or day_number == len(time_indices):
                print(f'  {day_number}/{len(time_indices)} days')

        summaries = summarise_wind_matrix(daily)
        year_output = reefs[reef_columns].copy()
        year_output['year'] = year
        year_output['ereefs_latitude'] = reefs['ereefs_latitude']
        year_output['ereefs_longitude'] = reefs['ereefs_longitude']
        year_output['ereefs_lat_index'] = reefs['ereefs_lat_index']
        year_output['ereefs_lon_index'] = reefs['ereefs_lon_index']
        year_output['ereefs_wind_source'] = SOURCE
        year_output['ereefs_wind_period_start'] = times[
            time_indices[0]
        ].date()
        year_output['ereefs_wind_period_end'] = times[
            time_indices[-1]
        ].date()
        year_output['ereefs_wind_requested_days'] = (
            requested_end - requested_start
        ).days + 1
        year_output['ereefs_wind_available_days'] = len(time_indices)
        year_output['ereefs_wind_complete_q1'] = (
            times[time_indices[-1]] >= requested_end
        )
        for name, values in summaries.items():
            year_output[name] = values[reef_cell_index]
        year_output.to_csv(part_path, index=False)
        outputs.append(year_output)

    output = pd.concat(outputs, ignore_index=True)
    if output.duplicated(reef_columns + ['year']).any():
        raise RuntimeError('Duplicate eReefs wind feature-year rows')
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(OUTPUT_FILE, index=False)
    coverage = output.groupby('year').agg(
        rows=('LABEL_ID', 'size'),
        days=('ereefs_wind_available_days', 'first'),
        complete_q1=('ereefs_wind_complete_q1', 'first'),
        mean_wind=('ereefs_wind_mean', 'mean'),
        median_calm_fraction=(
            'ereefs_wind_fraction_below_3', 'median'
        ),
    )
    print('\n', coverage)
    print(f'Wrote {len(output)} rows to {OUTPUT_FILE}')


if __name__ == '__main__':
    main()
