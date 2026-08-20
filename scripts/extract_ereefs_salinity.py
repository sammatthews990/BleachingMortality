'''Extract Q1 surface salinity and freshwater exposure from eReefs GBR4 v2.

The source ends on 17 January 2024. Values for 2024 are therefore retained as
partial-period observations and are never labelled complete Q1 exposure.
Freshwater exposure follows the AIMS definition: daily salinity deficit below
26 or 28 PSU, accumulated in PSU-days.
'''

from pathlib import Path

import numpy as np
import pandas as pd
import s3fs
import zarr
from scipy.spatial import cKDTree


SOURCE = (
    'AIMS eReefs GBR4 hydrodynamic v2 daily aggregation; surface k=16 (-0.5 m)'
)
STORE = 'gbr-dms-data-public/aims-ereefs-agg-hydrodynamic-4km-daily/data.zarr'
EVENT_YEARS = (2016, 2017, 2020, 2022, 2024)


def summarise_salinity(salinity_daily):
    '''Return daily salinity summaries for columns representing reef cells.'''
    salinity = np.asarray(salinity_daily, dtype=float)
    valid = np.isfinite(salinity)
    n_days = valid.sum(axis=0)
    valid_columns = n_days > 0
    result = {'salinity_n_days': n_days}
    for name in ('salinity_mean', 'salinity_q10', 'salinity_min'):
        result[name] = np.full(n_days.shape, np.nan, dtype=float)
    observed = salinity[:, valid_columns]
    result['salinity_mean'][valid_columns] = np.nanmean(observed, axis=0)
    result['salinity_q10'][valid_columns] = np.nanquantile(
        observed, 0.10, axis=0
    )
    result['salinity_min'][valid_columns] = np.nanmin(observed, axis=0)
    for threshold in (26.0, 28.0):
        deficit = np.where(valid, np.maximum(threshold - salinity, 0), np.nan)
        exposure = np.nansum(deficit, axis=0)
        exposure[n_days == 0] = np.nan
        result[f'freshwater_exposure_{int(threshold)}'] = exposure
        below = np.where(valid, salinity < threshold, False).sum(axis=0)
        result[f'days_below_{int(threshold)}'] = below.astype(float)
        result[f'days_below_{int(threshold)}'][n_days == 0] = np.nan
    return result


def main():
    environment_path = Path('data/processed/cheung_recreated_gbr_full.csv')
    output_path = Path(
        'data/processed/ereefs_salinity_imos_kd490_reef_year.csv'
    )
    part_dir = Path('data/processed/ereefs_salinity_parts')
    part_dir.mkdir(parents=True, exist_ok=True)
    environment = pd.read_csv(environment_path)
    required = {
        'year', 'LABEL_ID', 'LOC_NAME_S', 'lon', 'lat',
        'secc3m', 'k490_q90', 'secc3m_p10',
    }
    missing = required.difference(environment.columns)
    if missing:
        raise ValueError(f'Environmental table lacks columns: {sorted(missing)}')

    reef_columns = ['LABEL_ID', 'LOC_NAME_S', 'lon', 'lat']
    reefs = environment[reef_columns].drop_duplicates().reset_index(drop=True)
    if len(reefs) != 7063:
        raise ValueError(f'Expected 7063 reef features, found {len(reefs)}')

    fs = s3fs.S3FileSystem(anon=True)
    store = s3fs.S3Map(STORE, s3=fs)
    ereefs = zarr.open_group(store, mode='r')
    times = pd.to_datetime(
        ereefs['time'][:], unit='D', origin='1990-01-01'
    )
    grid_lats = ereefs['latitude'][:]
    grid_lons = ereefs['longitude'][:]

    surface_sample = ereefs['salt'][0, 16, :, :]
    valid_water = np.isfinite(surface_sample)
    water_lat, water_lon = np.meshgrid(grid_lats, grid_lons, indexing='ij')
    valid_coords = np.column_stack(
        [water_lat[valid_water], water_lon[valid_water]]
    )
    tree = cKDTree(valid_coords)
    _, nearest = tree.query(reefs[['lat', 'lon']].to_numpy())
    water_rows, water_cols = np.where(valid_water)
    reef_rows = water_rows[nearest]
    reef_cols = water_cols[nearest]

    outputs = []
    for year in EVENT_YEARS:
        part_path = part_dir / f'ereefs_salinity_{year}.csv'
        if part_path.exists():
            print(f'Loading cached salinity for {year}: {part_path}', flush=True)
            outputs.append(pd.read_csv(part_path))
            continue
        requested_start = pd.Timestamp(f'{year}-01-01')
        requested_end = pd.Timestamp(f'{year}-03-31')
        time_index = np.where(
            (times >= requested_start) & (times <= requested_end)
        )[0]
        if len(time_index) == 0:
            raise RuntimeError(f'No eReefs salinity observations for {year}')

        print(
            f'Extracting {year}: {times[time_index[0]].date()} to '
            f'{times[time_index[-1]].date()} ({len(time_index)} days)',
            flush=True,
        )
        frames = []
        for position, time_index_value in enumerate(time_index, start=1):
            surface = ereefs['salt'][time_index_value, 16, :, :]
            frames.append(surface[reef_rows, reef_cols])
            if position % 15 == 0 or position == len(time_index):
                print(f'  {position}/{len(time_index)} days', flush=True)
        daily = np.asarray(frames)
        summaries = summarise_salinity(daily)

        year_output = reefs.copy()
        year_output['year'] = year
        year_output['ereefs_latitude'] = grid_lats[reef_rows]
        year_output['ereefs_longitude'] = grid_lons[reef_cols]
        year_output['ereefs_lat_index'] = reef_rows
        year_output['ereefs_lon_index'] = reef_cols
        year_output['salinity_source'] = SOURCE
        year_output['salinity_period_start'] = times[time_index[0]].date()
        year_output['salinity_period_end'] = times[time_index[-1]].date()
        year_output['salinity_requested_days'] = (
            requested_end - requested_start
        ).days + 1
        year_output['salinity_available_days'] = len(time_index)
        year_output['salinity_complete_q1'] = (
            times[time_index[-1]] >= requested_end
        )
        year_output['salinity_forcing_warning'] = year == 2022
        for name, values in summaries.items():
            year_output[name] = values
        year_output.to_csv(part_path, index=False)
        outputs.append(year_output)

    salinity = pd.concat(outputs, ignore_index=True)
    # Some registry labels and names refer to multiple spatial features. Use
    # rounded coordinates to preserve those features without relying on the
    # final decimal places changed by CSV round-tripping.
    join_keys = [
        'LABEL_ID', 'LOC_NAME_S', '_lon_key', '_lat_key', 'year'
    ]
    salinity['_lon_key'] = salinity['lon'].round(6)
    salinity['_lat_key'] = salinity['lat'].round(6)
    optical = environment.assign(
        _lon_key=environment['lon'].round(6),
        _lat_key=environment['lat'].round(6),
    )[join_keys + ['secc3m', 'k490_q90', 'secc3m_p10']]
    combined = salinity.merge(
        optical,
        on=join_keys,
        how='left',
        validate='one_to_one',
        indicator='_optical_join',
    )
    if not combined['_optical_join'].eq('both').all():
        raise RuntimeError('At least one salinity row did not match an IMOS row')
    combined.drop(
        columns=['_optical_join', '_lon_key', '_lat_key'], inplace=True
    )
    if combined.duplicated(reef_columns + ['year']).any():
        raise RuntimeError('Duplicate reef-feature-year rows after optical join')

    output_path.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output_path, index=False)
    coverage = combined.groupby('year').agg(
        rows=('LABEL_ID', 'size'),
        salinity_days=('salinity_available_days', 'first'),
        complete_q1=('salinity_complete_q1', 'first'),
        mean_salinity=('salinity_mean', 'mean'),
        reefs_exposed_28=('freshwater_exposure_28', lambda x: int((x > 0).sum())),
    )
    print(coverage)
    print(f'Wrote {len(combined)} rows to {output_path}')


if __name__ == '__main__':
    main()
