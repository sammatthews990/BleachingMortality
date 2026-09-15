'''Extract IMOS night-SST, MHW and paired day/night metrics at reef-events.

The source-control feature reconstructs NOAA-style Degree Heating Weeks from
the IMOS 0.02-degree night-only SST and the SSTAARS 1992-2016 climatology:

* HotSpot = night SST - SSTAARS maximum monthly mean (MMM);
* retain the full HotSpot only when HotSpot is at least 1 degree C;
* sum over 84 calendar days and divide by seven;
* take the maximum during the event calendar year, matching ``ann_maxdhw``.

Because the IMOS daily composite is not gap-filled, short gaps of at most two
days are linearly interpolated. Remaining missing days are coverage-normalised
only when at least 70 percent of an 84-day window is available. Both the raw
observed and coverage-normalised estimates are retained for audit.

The same event window also yields Hobday MHW summaries and a multiscale
intensity-duration profile. A separate percentile-based heat dose follows Li
Shing Hiung et al. (2026): days above the local SSTAARS 90th percentile gate
the full anomaly above the daily climatological mean, accumulated within
prespecified 10-, 11- and 12-week windows centred on the local climatological
summer peak. Daytime SST is used only for an explicitly labelled observed
day-to-following-night sensitivity.
'''

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import xarray as xr


NIGHT_SST_URL = (
    's3://gbr-dms-data-public/imos-srs-sst-ms-1day-night/data.zarr'
)
DAY_SST_URL = (
    's3://gbr-dms-data-public/imos-srs-sst-ms-1day-day/data.zarr'
)
SSTAARS_NCSS_URL = (
    'https://thredds.aodn.org.au/thredds/ncss/grid/'
    'CSIRO/Climatology/SSTAARS/2017/SSTAARS.nc'
)
SSTAARS_VARIABLES = (
    'tm', 'ta_real', 'ta_imag', 't2a_real', 't2a_imag',
    't3a_real', 't3a_imag', 't4a_real', 't4a_imag',
)
SSTAARS_PERCENTILE_VARIABLES = ('TEMP_50th_perc', 'TEMP_90th_perc')
SST_OFFSET_C = 0.17
ROLLING_DAYS = 84
INTENSITY_DURATION_DAYS = (1, 3, 7, 14, 28, 56, 84)
PERCENTILE_HEAT_DOSE_WEEKS = (10, 11, 12)


def haversine_km(lat0: float, lon0: float, lats, lons):
    '''Great-circle distance from one point to arrays of points.'''
    lat0_r = np.radians(lat0)
    lats_r = np.radians(np.asarray(lats))
    dlat = np.radians(np.asarray(lats) - lat0)
    dlon = np.radians(np.asarray(lons) - lon0)
    a = (
        np.sin(dlat / 2) ** 2
        + np.cos(lat0_r) * np.cos(lats_r) * np.sin(dlon / 2) ** 2
    )
    return 6371.0 * 2 * np.arcsin(np.sqrt(a))


def reconstruct_sstaars_daily(
    coefficients: dict[str, np.ndarray], day_of_year
) -> np.ndarray:
    '''Reconstruct the SSTAARS daily climatology from harmonic coefficients.'''
    day = np.asarray(day_of_year, dtype=float)
    scalar_day = day.ndim == 0
    day = np.atleast_1d(day)
    tm = np.atleast_1d(np.asarray(coefficients['tm'], dtype=float))
    theta = 2 * np.pi * day[:, None] / 365.25
    result = np.broadcast_to(tm, (len(day), len(tm))).astype(float).copy()
    for harmonic in range(1, 5):
        stem = 'ta' if harmonic == 1 else f't{harmonic}a'
        real = np.atleast_1d(np.asarray(
            coefficients[f'{stem}_real'], dtype=float
        ))
        imag = np.atleast_1d(np.asarray(
            coefficients[f'{stem}_imag'], dtype=float
        ))
        # The source lineage uses the conjugate of exp(i * theta), yielding
        # real*cos(theta) + imag*sin(theta). This was checked against the
        # published daily SSTAARS field before implementation.
        result += (
            real[None, :] * np.cos(harmonic * theta)
            + imag[None, :] * np.sin(harmonic * theta)
        )
    if scalar_day:
        return result[0]
    return result


def sstaars_mmm(coefficients: dict[str, np.ndarray]) -> np.ndarray:
    '''Maximum monthly mean derived from the same SSTAARS climatology.'''
    dates = pd.date_range('2005-01-01', '2005-12-31', freq='D')
    daily = reconstruct_sstaars_daily(coefficients, dates.dayofyear.to_numpy())
    monthly = np.vstack([
        np.nanmean(daily[dates.month == month], axis=0)
        for month in range(1, 13)
    ])
    return np.nanmax(monthly, axis=0)


def fill_short_gaps(values, max_gap: int = 2) -> np.ndarray:
    '''Linearly fill only internally bounded missing runs no longer than max_gap.'''
    output = np.asarray(values, dtype=float).copy()
    missing = ~np.isfinite(output)
    index = 0
    while index < len(output):
        if not missing[index]:
            index += 1
            continue
        start = index
        while index < len(output) and missing[index]:
            index += 1
        end = index
        length = end - start
        if (
            length <= max_gap
            and start > 0
            and end < len(output)
            and np.isfinite(output[start - 1])
            and np.isfinite(output[end])
        ):
            output[start:end] = np.linspace(
                output[start - 1], output[end], length + 2
            )[1:-1]
    return output


def longest_true_run(flags) -> int:
    '''Length of the longest consecutive True run.'''
    longest = 0
    current = 0
    for flag in np.asarray(flags, dtype=bool):
        current = current + 1 if flag else 0
        longest = max(longest, current)
    return int(longest)


def interpolate_monthly_climatology(monthly_values, dates) -> np.ndarray:
    '''Periodically interpolate twelve mid-month values to daily dates.'''
    values = np.asarray(monthly_values, dtype=float)
    if values.shape[0] != 12:
        raise ValueError('Expected twelve monthly climatology values')
    dates = pd.DatetimeIndex(dates)
    centres = pd.DatetimeIndex([
        pd.Timestamp(2005, month, 1)
        + (pd.Timestamp(2005 if month < 12 else 2006, month % 12 + 1, 1)
           - pd.Timestamp(2005, month, 1)) / 2
        for month in range(1, 13)
    ]).dayofyear.to_numpy(dtype=float)
    x = np.concatenate(([centres[-1] - 365], centres, [centres[0] + 365]))
    flat = values.reshape(12, -1)
    extended = np.vstack((flat[-1], flat, flat[0]))
    target = dates.dayofyear.to_numpy(dtype=float)
    target -= ((dates.is_leap_year) & (dates.month > 2)).astype(float)
    output = np.column_stack([
        np.interp(target, x, extended[:, column])
        for column in range(extended.shape[1])
    ])
    return output.reshape((len(dates),) + values.shape[1:])


def summarise_marine_heatwaves(
    sst_c, climatology_c, median50_c, threshold90_c, minimum_duration: int = 5,
    join_gaps: int = 2,
) -> dict[str, float]:
    '''Summarise standard Hobday events from a complete daily summer series.'''
    sst = np.asarray(sst_c, dtype=float)
    climatology = np.asarray(climatology_c, dtype=float)
    median = np.asarray(median50_c, dtype=float)
    threshold = np.asarray(threshold90_c, dtype=float)
    finite = (
        np.isfinite(sst) & np.isfinite(climatology)
        & np.isfinite(median) & np.isfinite(threshold)
    )
    exceedance = finite & (sst > threshold)
    runs = []
    start = None
    for index, hot in enumerate(np.append(exceedance, False)):
        if hot and start is None:
            start = index
        elif not hot and start is not None:
            if index - start >= minimum_duration:
                runs.append([start, index - 1])
            start = None
    joined = []
    for run in runs:
        if joined and run[0] - joined[-1][1] - 1 <= join_gaps:
            joined[-1][1] = run[1]
        else:
            joined.append(run)
    result = {
        'imos_mhw_event_count': float(len(joined)),
        'imos_mhw_max_duration_days': 0.0,
        'imos_mhw_max_intensity_c': 0.0,
        'imos_mhw_cumulative_intensity_c_days': 0.0,
        'imos_mhw_category2_days': 0.0,
        'imos_mhw_max_category': 0.0,
    }
    if not joined:
        return result
    anomaly = sst - climatology
    scale = threshold - median
    category = np.zeros(len(sst), dtype=float)
    eligible = finite & (sst > threshold) & (scale > 0)
    category[eligible] = np.floor(
        1 + (sst[eligible] - threshold[eligible]) / scale[eligible]
    )
    category = np.clip(category, 0, 4)
    result['imos_mhw_max_duration_days'] = float(max(
        end - start + 1 for start, end in joined
    ))
    result['imos_mhw_max_intensity_c'] = float(max(
        np.nanmax(anomaly[start:end + 1]) for start, end in joined
    ))
    result['imos_mhw_cumulative_intensity_c_days'] = float(sum(
        np.nansum(anomaly[start:end + 1]) for start, end in joined
    ))
    event_mask = np.zeros(len(sst), dtype=bool)
    for start, end in joined:
        event_mask[start:end + 1] = True
    result['imos_mhw_category2_days'] = float(
        np.sum(event_mask & (category >= 2))
    )
    result['imos_mhw_max_category'] = float(np.max(category[event_mask]))
    return result


def summarise_percentile_heat_dose(
    dates,
    sst_c,
    climatology_c,
    threshold90_c,
    event_year: int,
    minimum_coverage: float = 0.70,
    raw_sst_c=None,
) -> dict[str, object]:
    '''Accumulate percentile-gated climatological anomaly in warm windows.

    Unlike a Hobday event summary, this metric has no minimum run duration.
    Every observed day above the daily 90th percentile contributes its full
    SST anomaly above the daily climatological mean. Windows are centred on
    the local climatological peak between November and April.
    '''
    dates = pd.DatetimeIndex(dates).normalize()
    sst = np.asarray(sst_c, dtype=float)
    climatology = np.asarray(climatology_c, dtype=float)
    threshold = np.asarray(threshold90_c, dtype=float)
    raw_sst = sst if raw_sst_c is None else np.asarray(raw_sst_c, dtype=float)
    if not (
        len(dates) == len(sst) == len(climatology)
        == len(threshold) == len(raw_sst)
    ):
        raise ValueError('Percentile heat-dose inputs must have equal lengths')

    result: dict[str, object] = {
        'imos_pbd_peak_climatology_date': '',
    }
    for weeks in PERCENTILE_HEAT_DOSE_WEEKS:
        result.update({
            f'imos_pbd{weeks}_c_weeks': np.nan,
            f'imos_pbd{weeks}_c_weeks_observed': np.nan,
            f'imos_pbd{weeks}_coverage': np.nan,
            f'imos_pbd{weeks}_raw_coverage': np.nan,
            f'imos_pbd{weeks}_exceedance_days': np.nan,
        })

    warm_season = (
        (dates >= pd.Timestamp(f'{event_year - 1}-11-01'))
        & (dates <= pd.Timestamp(f'{event_year}-04-30'))
        & np.isfinite(climatology)
    )
    if not warm_season.any():
        return result
    warm_indices = np.flatnonzero(warm_season)
    peak_index = warm_indices[np.nanargmax(climatology[warm_season])]
    peak_date = dates[peak_index]
    result['imos_pbd_peak_climatology_date'] = peak_date.strftime('%Y-%m-%d')

    finite = (
        np.isfinite(sst) & np.isfinite(climatology) & np.isfinite(threshold)
    )
    anomaly = sst - climatology
    gated_anomaly = np.where(
        finite,
        np.where(sst > threshold, np.maximum(anomaly, 0.0), 0.0),
        np.nan,
    )
    for weeks in PERCENTILE_HEAT_DOSE_WEEKS:
        days = weeks * 7
        start = peak_date - pd.Timedelta(days=(days - 1) // 2)
        window = (dates >= start) & (dates < start + pd.Timedelta(days=days))
        if int(window.sum()) != days:
            continue
        valid_days = int(finite[window].sum())
        coverage = valid_days / days
        raw_finite = (
            np.isfinite(raw_sst[window])
            & np.isfinite(climatology[window])
            & np.isfinite(threshold[window])
        )
        result[f'imos_pbd{weeks}_coverage'] = float(coverage)
        result[f'imos_pbd{weeks}_raw_coverage'] = float(
            raw_finite.sum() / days
        )
        if coverage < minimum_coverage:
            continue
        observed_dose = float(np.nansum(gated_anomaly[window]) / 7)
        result.update({
            f'imos_pbd{weeks}_c_weeks': observed_dose / coverage,
            f'imos_pbd{weeks}_c_weeks_observed': observed_dose,
            f'imos_pbd{weeks}_exceedance_days': float(np.sum(
                finite[window] & (sst[window] > threshold[window])
            )),
        })
    return result


def summarise_intensity_duration_profile(
    hotspot_c, minimum_coverage: float = 0.70,
) -> dict[str, float]:
    '''Compress multiscale maximum mean HotSpot into level and persistence.'''
    hotspot = np.asarray(hotspot_c, dtype=float)
    maxima = []
    result = {}
    for duration in INTENSITY_DURATION_DAYS:
        series = pd.Series(hotspot)
        rolling = series.rolling(
            duration, min_periods=math.ceil(duration * minimum_coverage)
        ).mean()
        if duration > 1:
            rolling.iloc[:duration - 1] = np.nan
        maximum = float(rolling.max()) if rolling.notna().any() else np.nan
        result[f'imos_id_maxmean_hotspot_{duration}d_c'] = maximum
        maxima.append(maximum)
    valid = np.isfinite(maxima)
    result['imos_id_hotspot_level_c'] = np.nan
    result['imos_id_hotspot_persistence_slope_c_per_log_day'] = np.nan
    if valid.sum() >= 4:
        log_duration = np.log(np.asarray(INTENSITY_DURATION_DAYS)[valid])
        values = np.asarray(maxima)[valid]
        result['imos_id_hotspot_level_c'] = float(
            np.trapezoid(values, log_duration) /
            (log_duration[-1] - log_duration[0])
        )
        result['imos_id_hotspot_persistence_slope_c_per_log_day'] = float(
            np.polyfit(log_duration, values, 1)[0]
        )
    return result


def summarise_site_year(
    dates,
    night_sst_c,
    mmm_c: float,
    event_year: int,
    climatology_coefficients: dict[str, float] | None = None,
    monthly_p50_c=None,
    monthly_p90_c=None,
    minimum_coverage: float = 0.70,
    max_gap_days: int = 2,
) -> dict[str, object]:
    '''Calculate matched DHW, hot-spell, peak and hot-night metrics.'''
    series = pd.Series(
        np.asarray(night_sst_c, dtype=float),
        index=pd.DatetimeIndex(dates).normalize(),
    ).groupby(level=0).mean()
    full_index = pd.date_range(
        f'{event_year - 1}-08-01', f'{event_year}-12-31', freq='D'
    )
    original = series.reindex(full_index).to_numpy()
    filled = fill_short_gaps(original, max_gap=max_gap_days)
    finite = np.isfinite(filled)
    hot_spot = filled - float(mmm_c)
    dhw_contribution = np.where(
        finite,
        np.where(hot_spot >= 1.0, hot_spot, 0.0),
        np.nan,
    )

    contribution = pd.Series(dhw_contribution, index=full_index)
    valid = pd.Series(finite.astype(int), index=full_index)
    original_valid = pd.Series(
        np.isfinite(original).astype(int), index=full_index
    )
    rolling_sum = contribution.rolling(
        ROLLING_DAYS, min_periods=1
    ).sum()
    rolling_n = valid.rolling(ROLLING_DAYS).sum()
    rolling_original_n = original_valid.rolling(ROLLING_DAYS).sum()
    enough = rolling_n >= math.ceil(ROLLING_DAYS * minimum_coverage)
    full_window = pd.Series(1, index=full_index).rolling(ROLLING_DAYS).sum()
    enough &= full_window == ROLLING_DAYS
    dhw_observed = (rolling_sum / 7).where(enough)
    dhw_scaled = (
        rolling_sum * ROLLING_DAYS / rolling_n / 7
    ).where(enough)
    calendar = (full_index.year == event_year)
    calendar_scaled = dhw_scaled[calendar]
    calendar_observed = dhw_observed[calendar]

    result: dict[str, object] = {
        'imos_dhw1_84_max': np.nan,
        'imos_dhw1_84_max_observed': np.nan,
        'imos_dhw1_84_peak_date': '',
        'imos_dhw1_84_peak_coverage': np.nan,
        'imos_dhw1_84_peak_raw_coverage': np.nan,
        'imos_dhw1_84_valid_windows': int(calendar_scaled.notna().sum()),
    }
    if calendar_scaled.notna().any():
        peak_date = calendar_scaled.idxmax()
        result.update({
            'imos_dhw1_84_max': float(calendar_scaled.loc[peak_date]),
            'imos_dhw1_84_max_observed': float(
                calendar_observed.loc[peak_date]
            ),
            'imos_dhw1_84_peak_date': peak_date.strftime('%Y-%m-%d'),
            'imos_dhw1_84_peak_coverage': float(
                rolling_n.loc[peak_date] / ROLLING_DAYS
            ),
            'imos_dhw1_84_peak_raw_coverage': float(
                rolling_original_n.loc[peak_date] / ROLLING_DAYS
            ),
        })

    daily_climatology = None
    daily_p50 = None
    daily_p90 = None
    if (
        climatology_coefficients is not None
        and monthly_p50_c is not None
        and monthly_p90_c is not None
    ):
        daily_climatology = reconstruct_sstaars_daily(
            climatology_coefficients,
            full_index.dayofyear.to_numpy(),
        ).reshape(-1)
        daily_p50 = interpolate_monthly_climatology(
            monthly_p50_c, full_index
        ).reshape(-1)
        daily_p90 = interpolate_monthly_climatology(
            monthly_p90_c, full_index
        ).reshape(-1)
        result.update(summarise_percentile_heat_dose(
            full_index,
            filled,
            daily_climatology,
            daily_p90,
            event_year,
            minimum_coverage=minimum_coverage,
            raw_sst_c=original,
        ))
    else:
        result.update(summarise_percentile_heat_dose(
            [], [], [], [], event_year,
            minimum_coverage=minimum_coverage,
        ))

    event_mask = (
        (full_index >= pd.Timestamp(f'{event_year - 1}-11-01'))
        & (full_index <= pd.Timestamp(f'{event_year}-04-30'))
    )
    event_sst = filled[event_mask]
    event_hot_spot = event_sst - float(mmm_c)
    event_finite = np.isfinite(event_sst)
    event_raw_finite = np.isfinite(original[event_mask])
    event_days = int(event_mask.sum())
    event_coverage = float(event_finite.sum() / event_days)
    result.update({
        'imos_event_n_days': event_days,
        'imos_event_n_valid': int(event_finite.sum()),
        'imos_event_n_raw_valid': int(event_raw_finite.sum()),
        'imos_event_coverage': event_coverage,
        'imos_event_raw_coverage': float(event_raw_finite.sum() / event_days),
        'imos_hotspell_mmm1_max_days': np.nan,
        'imos_hotspot_3d_max': np.nan,
        'imos_extreme_dhd_mmm2': np.nan,
        'imos_hot_nights_mmm1_count': np.nan,
        'imos_hot_night_fraction': np.nan,
        'imos_mhw_event_count': np.nan,
        'imos_mhw_max_duration_days': np.nan,
        'imos_mhw_max_intensity_c': np.nan,
        'imos_mhw_cumulative_intensity_c_days': np.nan,
        'imos_mhw_category2_days': np.nan,
        'imos_mhw_max_category': np.nan,
        'imos_id_hotspot_level_c': np.nan,
        'imos_id_hotspot_persistence_slope_c_per_log_day': np.nan,
    })
    result.update({
        f'imos_id_maxmean_hotspot_{duration}d_c': np.nan
        for duration in INTENSITY_DURATION_DAYS
    })
    if event_coverage >= minimum_coverage:
        hot_nights = event_finite & (event_hot_spot >= 1.0)
        three_day_peak = pd.Series(event_hot_spot).rolling(
            3, min_periods=2
        ).mean().max()
        extreme = np.where(
            event_finite, np.maximum(event_hot_spot - 2.0, 0), np.nan
        )
        result.update({
            'imos_hotspell_mmm1_max_days': longest_true_run(hot_nights),
            'imos_hotspot_3d_max': float(three_day_peak),
            'imos_extreme_dhd_mmm2': float(
                np.nansum(extreme) * event_days / event_finite.sum()
            ),
            'imos_hot_nights_mmm1_count': int(hot_nights.sum()),
            'imos_hot_night_fraction': float(
                hot_nights.sum() / event_finite.sum()
            ),
        })
        result.update(summarise_intensity_duration_profile(
            event_hot_spot, minimum_coverage=minimum_coverage
        ))
        if daily_climatology is not None:
            # AusTemp categories use the published SSTAARS median and 90th
            # percentile; Hobday event duration additionally applies the
            # five-day and two-day joining rules here.
            result.update(summarise_marine_heatwaves(
                event_sst,
                daily_climatology[event_mask],
                daily_p50[event_mask],
                daily_p90[event_mask],
            ))
    return result


def download_sstaars_subset(
    rows: pd.DataFrame, path: Path, padding_degrees: float = 0.10
) -> None:
    '''Download a small GBR coefficient subset from the official AODN NCSS.'''
    parameters = [
        ('var', variable)
        for variable in SSTAARS_VARIABLES + SSTAARS_PERCENTILE_VARIABLES
    ]
    parameters.extend([
        ('north', str(float(rows['lat'].max() + padding_degrees))),
        ('south', str(float(rows['lat'].min() - padding_degrees))),
        ('west', str(float(rows['lon'].min() - padding_degrees))),
        ('east', str(float(rows['lon'].max() + padding_degrees))),
        ('time_start', '2005-01-01T00:00:00Z'),
        ('time_end', '2005-12-31T23:59:59Z'),
        ('horizStride', '1'),
        ('addLatLon', 'true'),
        ('accept', 'netcdf3'),
    ])
    response = requests.get(SSTAARS_NCSS_URL, params=parameters, timeout=300)
    response.raise_for_status()
    if not response.content.startswith(b'CDF'):
        raise ValueError('SSTAARS subset response was not NetCDF3')
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix('.tmp.nc')
    temporary.write_bytes(response.content)
    temporary.replace(path)


def nearest_sstaars_pixels(
    rows: pd.DataFrame, coefficients: xr.Dataset, radius_km: float = 5.0
) -> pd.DataFrame:
    '''Choose one stable, nearest finite-water climatology pixel per row.'''
    grid_lats = coefficients['LATITUDE'].values
    grid_lons = coefficients['LONGITUDE'].values
    valid = np.isfinite(coefficients['tm'].values)
    matches = []
    for row in rows.itertuples(index=False):
        lat_i = int(np.abs(grid_lats - row.lat).argmin())
        lon_i = int(np.abs(grid_lons - row.lon).argmin())
        lat_step = float(np.median(np.abs(np.diff(grid_lats))))
        lon_step = float(np.median(np.abs(np.diff(grid_lons))))
        lat_cells = int(math.ceil(radius_km / (111 * lat_step))) + 1
        lon_cells = int(math.ceil(
            radius_km
            / (111 * max(math.cos(math.radians(row.lat)), 0.1) * lon_step)
        )) + 1
        lat_slice = slice(
            max(0, lat_i - lat_cells), min(len(grid_lats), lat_i + lat_cells + 1)
        )
        lon_slice = slice(
            max(0, lon_i - lon_cells), min(len(grid_lons), lon_i + lon_cells + 1)
        )
        local_lats, local_lons = np.meshgrid(
            grid_lats[lat_slice], grid_lons[lon_slice], indexing='ij'
        )
        distance = haversine_km(row.lat, row.lon, local_lats, local_lons)
        eligible = valid[lat_slice, lon_slice] & (distance <= radius_km)
        if not eligible.any():
            matches.append((np.nan, np.nan, np.nan, -1, -1))
            continue
        flat = np.where(eligible, distance, np.inf).argmin()
        local_i, local_j = np.unravel_index(flat, eligible.shape)
        i = int((lat_slice.start or 0) + local_i)
        j = int((lon_slice.start or 0) + local_j)
        matches.append((
            float(grid_lats[i]), float(grid_lons[j]),
            float(distance[local_i, local_j]), i, j,
        ))
    result = pd.DataFrame(matches, columns=[
        'imos_grid_lat', 'imos_grid_lon', 'imos_grid_distance_km',
        '_clim_lat_i', '_clim_lon_i',
    ])
    mmm = np.full(len(rows), np.nan)
    good = result['_clim_lat_i'] >= 0
    if good.any():
        i = result.loc[good, '_clim_lat_i'].astype(int).to_numpy()
        j = result.loc[good, '_clim_lon_i'].astype(int).to_numpy()
        selected = {
            variable: coefficients[variable].values[i, j]
            for variable in SSTAARS_VARIABLES
        }
        mmm[good.to_numpy()] = sstaars_mmm(selected)
    result['sstaars_mmm_c'] = mmm
    for variable in SSTAARS_VARIABLES:
        values = np.full(len(rows), np.nan)
        if good.any():
            values[good.to_numpy()] = coefficients[variable].values[i, j]
        result[f'_clim_{variable}'] = values
    for variable in SSTAARS_PERCENTILE_VARIABLES:
        for month_index in range(12):
            values = np.full(len(rows), np.nan)
            if good.any():
                values[good.to_numpy()] = coefficients[variable].values[
                    month_index, i, j
                ]
            result[f'_clim_{variable}_{month_index + 1:02d}'] = values
    return result


def load_corrected_sst_block(
    dataset: xr.Dataset, year: int, lats: xr.DataArray, lons: xr.DataArray,
    eligible_neighbor: np.ndarray, quality_min: int,
) -> tuple[pd.DatetimeIndex, np.ndarray]:
    '''Load one year-range of corrected SST for site-neighbour coordinates.'''
    archive_time = pd.DatetimeIndex(dataset['time'].values)
    start = pd.Timestamp(f'{year - 1}-07-31')
    end = pd.Timestamp(f'{year}-12-31T23:59:59')
    time_positions = np.flatnonzero(
        (archive_time >= start) & (archive_time <= end)
    )
    if len(time_positions) == 0:
        raise ValueError(f'No SST observations for {year}')
    block = dataset[[
        'sea_surface_temperature', 'sses_bias', 'quality_level'
    ]].isel(time=time_positions).sel(
        lat=lats, lon=lons, method='nearest'
    ).compute()
    local_dates = (
        pd.DatetimeIndex(block['time'].values)
        .tz_localize('UTC')
        .tz_convert('Australia/Brisbane')
        .tz_localize(None)
        .normalize()
    )
    corrected = (
        block['sea_surface_temperature'].values - 273.15
        + SST_OFFSET_C - block['sses_bias'].values
    )
    corrected[block['quality_level'].values < quality_min] = np.nan
    corrected[~np.isfinite(corrected)] = np.nan
    corrected = np.where(eligible_neighbor[None, :, :], corrected, np.nan)
    return local_dates, corrected


def choose_daily_neighbor(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    '''Prefer the stable centre pixel, then nearest ordered valid neighbour.'''
    has_value = np.isfinite(values)
    chosen = np.argmax(has_value, axis=1)
    has_any = has_value.any(axis=1)
    daily = np.full(values.shape[0], np.nan)
    daily[has_any] = values[np.flatnonzero(has_any), chosen[has_any]]
    return daily, chosen


def summarise_day_night_relief(
    day_dates, day_sst_c, night_dates, night_sst_c, mmm_c: float,
    event_year: int, minimum_hot_pairs: int = 5,
    minimum_pair_coverage: float = 0.50,
) -> dict[str, float]:
    '''Pair daytime SST with the following local night on observed hot days.'''
    day = pd.Series(
        np.asarray(day_sst_c, dtype=float),
        index=pd.DatetimeIndex(day_dates).normalize(), dtype=float,
    ) \
        .groupby(level=0).mean()
    night = pd.Series(
        np.asarray(night_sst_c, dtype=float),
        index=pd.DatetimeIndex(night_dates).normalize(), dtype=float,
    ).groupby(level=0).mean()
    dates = pd.date_range(
        f'{event_year - 1}-11-01', f'{event_year}-04-29', freq='D'
    )
    day_values = day.reindex(dates).to_numpy()
    following_night = night.reindex(dates + pd.Timedelta(days=1)).to_numpy()
    hot_day = np.isfinite(day_values) & (day_values >= float(mmm_c) + 1)
    paired = hot_day & np.isfinite(following_night)
    hot_count = int(hot_day.sum())
    paired_count = int(paired.sum())
    pair_coverage = paired_count / hot_count if hot_count else np.nan
    result = {
        'imos_hot_day_observed_count': float(hot_count),
        'imos_hot_day_night_pair_count': float(paired_count),
        'imos_hot_day_night_pair_coverage': float(pair_coverage),
        'imos_day_to_following_night_drop_c': np.nan,
        'imos_night_relief_fraction': np.nan,
        'imos_unrelieved_hot_pair_fraction': np.nan,
    }
    if (
        paired_count >= minimum_hot_pairs
        and np.isfinite(pair_coverage)
        and pair_coverage >= minimum_pair_coverage
    ):
        relieved = following_night[paired] < float(mmm_c) + 1
        result.update({
            'imos_day_to_following_night_drop_c': float(np.median(
                day_values[paired] - following_night[paired]
            )),
            'imos_night_relief_fraction': float(np.mean(relieved)),
            'imos_unrelieved_hot_pair_fraction': float(np.mean(~relieved)),
        })
    return result


def extract_year(
    night: xr.Dataset,
    day: xr.Dataset,
    rows: pd.DataFrame,
    matches: pd.DataFrame,
    year: int,
    quality_min: int,
    minimum_coverage: float,
    max_gap_days: int,
    daily_fallback_radius_km: float,
    maximum_match_radius_km: float,
) -> pd.DataFrame:
    '''Extract one event-year block and calculate site summaries.'''
    valid_match = matches['imos_grid_lat'].notna().to_numpy()
    result = pd.concat([
        rows.reset_index(drop=True),
        matches.loc[:, ~matches.columns.str.startswith('_')].reset_index(drop=True),
    ], axis=1)
    summaries = [None] * len(rows)
    if valid_match.any():
        selected = matches.loc[valid_match].reset_index(drop=True)
        lat_step = float(np.median(np.abs(np.diff(night['lat'].values))))
        lon_step = float(np.median(np.abs(np.diff(night['lon'].values))))
        offsets = [(0, 0)] + [
            (di, dj)
            for di in range(-2, 3)
            for dj in range(-2, 3)
            if (di, dj) != (0, 0)
        ]
        candidate_lats = np.column_stack([
            selected['imos_grid_lat'].to_numpy() + di * lat_step
            for di, _ in offsets
        ])
        candidate_lons = np.column_stack([
            selected['imos_grid_lon'].to_numpy() + dj * lon_step
            for _, dj in offsets
        ])
        reef_lats = rows.loc[valid_match, 'lat'].to_numpy()[:, None]
        reef_lons = rows.loc[valid_match, 'lon'].to_numpy()[:, None]
        centre_distance = haversine_km(
            selected['imos_grid_lat'].to_numpy()[:, None],
            selected['imos_grid_lon'].to_numpy()[:, None],
            candidate_lats,
            candidate_lons,
        )
        reef_candidate_distance = haversine_km(
            reef_lats, reef_lons, candidate_lats, candidate_lons
        )
        eligible_neighbor = (
            (centre_distance <= daily_fallback_radius_km)
            & (reef_candidate_distance <= maximum_match_radius_km)
        )
        eligible_neighbor[:, 0] = True
        lats = xr.DataArray(candidate_lats, dims=('site', 'neighbor'))
        lons = xr.DataArray(candidate_lons, dims=('site', 'neighbor'))
        # The public archive contains a small number of overlapping product
        # records and is therefore not strictly monotonic. Select by explicit
        # integer positions instead of label-based slicing; duplicate local
        # dates are collapsed in ``summarise_site_year``.
        local_dates, corrected_c = load_corrected_sst_block(
            night, year, lats, lons, eligible_neighbor, quality_min
        )
        day_dates, day_corrected_c = load_corrected_sst_block(
            day, year, lats, lons, eligible_neighbor, quality_min
        )
        valid_indices = np.flatnonzero(valid_match)
        for site, row_index in enumerate(valid_indices):
            site_values = corrected_c[:, site, :]
            daily_values, chosen_neighbor = choose_daily_neighbor(site_values)
            has_any = np.isfinite(daily_values)
            day_values, _ = choose_daily_neighbor(
                day_corrected_c[:, site, :]
            )
            selected_distance = haversine_km(
                float(rows.loc[row_index, 'lat']),
                float(rows.loc[row_index, 'lon']),
                candidate_lats[site], candidate_lons[site],
            )
            used_distance = selected_distance[chosen_neighbor[has_any]]
            result.loc[row_index, 'imos_daily_spatial_fallback_fraction'] = (
                float(np.mean(chosen_neighbor[has_any] != 0))
                if has_any.any() else np.nan
            )
            result.loc[row_index, 'imos_daily_match_distance_median_km'] = (
                float(np.median(used_distance)) if len(used_distance) else np.nan
            )
            result.loc[row_index, 'imos_daily_match_distance_max_km'] = (
                float(np.max(used_distance)) if len(used_distance) else np.nan
            )
            summaries[row_index] = summarise_site_year(
                local_dates,
                daily_values,
                result.loc[row_index, 'sstaars_mmm_c'],
                year,
                climatology_coefficients={
                    variable: matches.loc[row_index, f'_clim_{variable}']
                    for variable in SSTAARS_VARIABLES
                },
                monthly_p50_c=np.array([
                    matches.loc[
                        row_index, f'_clim_TEMP_50th_perc_{month:02d}'
                    ] for month in range(1, 13)
                ]),
                monthly_p90_c=np.array([
                    matches.loc[
                        row_index, f'_clim_TEMP_90th_perc_{month:02d}'
                    ] for month in range(1, 13)
                ]),
                minimum_coverage=minimum_coverage,
                max_gap_days=max_gap_days,
            )
            summaries[row_index].update(summarise_day_night_relief(
                day_dates, day_values, local_dates, daily_values,
                result.loc[row_index, 'sstaars_mmm_c'], year,
            ))
    empty_summary = summarise_site_year(
        [], [], np.nan, year,
        minimum_coverage=minimum_coverage,
        max_gap_days=max_gap_days,
    )
    empty_summary.update(summarise_day_night_relief(
        [], [], [], [], np.nan, year
    ))
    summaries = [summary or empty_summary for summary in summaries]
    result = pd.concat([result, pd.DataFrame(summaries)], axis=1)
    result['imos_quality_min'] = quality_min
    result['imos_sst_offset_c'] = SST_OFFSET_C
    result['imos_max_gap_days'] = max_gap_days
    result['imos_minimum_coverage'] = minimum_coverage
    result['imos_daily_fallback_radius_km'] = daily_fallback_radius_km
    result['imos_maximum_match_radius_km'] = maximum_match_radius_km
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--input', default='data/processed/environment_screen_grid_validation.csv'
    )
    parser.add_argument(
        '--output', default='data/processed/imos_thermal_metrics_validation.csv'
    )
    parser.add_argument(
        '--climatology-cache',
        default='data/cache/imos_thermal/sstaars_gbr_coefficients.nc',
    )
    parser.add_argument('--years', nargs='*', type=int)
    parser.add_argument('--quality-min', type=int, default=3, choices=(3, 4, 5))
    parser.add_argument('--minimum-coverage', type=float, default=0.70)
    parser.add_argument('--max-gap-days', type=int, default=2)
    parser.add_argument('--water-radius-km', type=float, default=5.0)
    parser.add_argument('--daily-fallback-radius-km', type=float, default=3.0)
    parser.add_argument('--force', action='store_true')
    parser.add_argument('--force-climatology', action='store_true')
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = pd.read_csv(args.input)
    required = {'ReefID', 'year', 'lon', 'lat'}
    missing = required - set(rows.columns)
    if missing:
        raise ValueError(f'Input lacks columns: {sorted(missing)}')
    rows['year'] = pd.to_numeric(rows['year'], errors='raise').astype(int)
    if args.years:
        rows = rows[rows['year'].isin(args.years)].copy()
    rows = rows.reset_index(drop=True)
    if rows.empty:
        raise ValueError('No reef-event rows selected')
    if rows.duplicated(['ReefID', 'year']).any():
        raise ValueError('Input has duplicate ReefID-year keys')
    if not 0 < args.minimum_coverage <= 1:
        raise ValueError('--minimum-coverage must be in (0, 1]')

    climatology_path = Path(args.climatology_cache)
    climatology_valid = False
    if climatology_path.exists() and not args.force_climatology:
        with xr.open_dataset(
            climatology_path, engine='scipy', decode_times=False
        ) as cached:
            climatology_valid = all(
                variable in cached
                for variable in SSTAARS_VARIABLES + SSTAARS_PERCENTILE_VARIABLES
            ) and cached.sizes.get('MONTH_OF_YEAR') == 12
    if not climatology_valid:
        print('Downloading SSTAARS harmonic coefficient subset', flush=True)
        download_sstaars_subset(rows, climatology_path)
    with xr.open_dataset(
        climatology_path, engine='scipy', decode_times=False
    ) as opened:
        climatology = opened.load()
    matches = nearest_sstaars_pixels(
        rows, climatology, radius_km=args.water_radius_km
    )

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    part_dir = output_path.parent / f'{output_path.stem}_parts'
    part_dir.mkdir(parents=True, exist_ok=True)
    night = None
    day = None

    parts = []
    for year in sorted(rows['year'].unique()):
        part_path = part_dir / f'{output_path.stem}_{year}.csv'
        use_cache = part_path.exists() and not args.force
        if use_cache:
            part = pd.read_csv(part_path)
            expanded_fields = {
                'imos_mhw_max_duration_days',
                'imos_mhw_max_intensity_c',
                'imos_mhw_cumulative_intensity_c_days',
                'imos_mhw_category2_days',
                'imos_pbd12_c_weeks',
                'imos_pbd12_coverage',
                'imos_day_to_following_night_drop_c',
                'imos_night_relief_fraction',
                'imos_id_hotspot_level_c',
                'imos_id_hotspot_persistence_slope_c_per_log_day',
            }
            cache_ok = (
                len(part) == int((rows['year'] == year).sum())
                and expanded_fields.issubset(part.columns)
                and (part['imos_quality_min'] == args.quality_min).all()
                and np.allclose(
                    part['imos_minimum_coverage'], args.minimum_coverage
                )
                and (part['imos_max_gap_days'] == args.max_gap_days).all()
                and np.allclose(
                    part['imos_daily_fallback_radius_km'],
                    args.daily_fallback_radius_km,
                )
                and np.allclose(
                    part['imos_maximum_match_radius_km'], args.water_radius_km
                )
            )
            use_cache = bool(cache_ok)
        if use_cache:
            print(f'Loading cached {year}: {part_path}', flush=True)
        else:
            if night is None:
                print(
                    'Opening IMOS 0.02-degree daily day/night SST', flush=True
                )
                night = xr.open_zarr(
                    NIGHT_SST_URL, storage_options={'anon': True}
                )
                day = xr.open_zarr(
                    DAY_SST_URL, storage_options={'anon': True}
                )
            selected_index = rows.index[rows['year'] == year].to_numpy()
            print(
                f'Extracting {year}: {len(selected_index)} reef-events',
                flush=True,
            )
            part = extract_year(
                night,
                day,
                rows.loc[selected_index].reset_index(drop=True),
                matches.loc[selected_index].reset_index(drop=True),
                int(year),
                args.quality_min,
                args.minimum_coverage,
                args.max_gap_days,
                args.daily_fallback_radius_km,
                args.water_radius_km,
            )
            part.to_csv(part_path, index=False)
        parts.append(part)

    if night is not None:
        night.close()
    if day is not None:
        day.close()

    combined = pd.concat(parts, ignore_index=True).sort_values(
        ['year', 'ReefID']
    )
    combined.to_csv(output_path, index=False)
    manifest = {
        'input': args.input,
        'output': args.output,
        'rows': int(len(combined)),
        'years': sorted(map(int, combined['year'].unique())),
        'night_sst_store': NIGHT_SST_URL,
        'day_sst_store': DAY_SST_URL,
        'climatology_service': SSTAARS_NCSS_URL,
        'climatology_cache': str(climatology_path),
        'climatology': 'SSTAARS 1992-2016 night SST harmonic fit',
        'quality_min': args.quality_min,
        'bias_correction': 'sst_celsius + 0.17 - sses_bias',
        'minimum_coverage': args.minimum_coverage,
        'max_gap_days': args.max_gap_days,
        'water_radius_km': args.water_radius_km,
        'daily_fallback_radius_km': args.daily_fallback_radius_km,
        'maximum_daily_match_radius_km': args.water_radius_km,
        'daily_spatial_fallback': (
            'use the centre 0.02-degree pixel when valid; otherwise the '
            'nearest quality-controlled pixel within the configured radius'
        ),
        'dhw_definition': (
            'calendar-year maximum of 84-day sum of full HotSpot when '
            'night SST - SSTAARS MMM >= 1 C, divided by 7; remaining '
            'missingness coverage-normalised'
        ),
        'event_metric_window': 'previous-year-11-01/current-year-04-30',
        'mhw_definition': (
            'night SST above periodically interpolated SSTAARS 90th '
            'percentile for at least 5 days; qualifying events separated '
            'by at most 2 days are joined; categories use SSTAARS p50/p90'
        ),
        'percentile_heat_dose_definition': (
            'Li Shing Hiung et al. (2026) alpha=1 analogue: each night-SST '
            'day above the periodically interpolated SSTAARS 90th percentile '
            'contributes its full positive anomaly above the SSTAARS daily '
            'climatological mean; accumulated in prespecified 10-, 11- and '
            '12-week windows centred on the local November-April climatology '
            'peak, with no minimum event duration; primary screen uses 12 weeks'
        ),
        'percentile_heat_dose_reference': 'doi:10.1029/2025GL119516',
        'intensity_duration_profile': (
            'maximum mean SSTAARS-MMM HotSpot at 1, 3, 7, 14, 28, 56 '
            'and 84 days, compressed to log-duration level and slope; '
            'no frequency or return period is inferred from IMOS'
        ),
        'cooling_definition': (
            'quality-controlled daytime skin SST paired by local GBR date '
            'to the following nighttime skin SST; summaries require at '
            'least 5 observed hot-day pairs and 50 percent pair coverage'
        ),
    }
    output_path.with_suffix('.manifest.json').write_text(
        json.dumps(manifest, indent=2), encoding='utf-8'
    )
    print(f'Wrote {len(combined)} rows to {output_path}', flush=True)


if __name__ == '__main__':
    main()
