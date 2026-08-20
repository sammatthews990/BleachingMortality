"""Extract Q1 Secchi-depth proxies from the IMOS MODIS-Aqua Kd490 product.

Reef centroids often fall on pixels masked as land or shallow reef. Values are
taken from the centroid grid cell when it contains observations and, only when
it does not, from the smallest surrounding-water neighbourhood with data (2,
3, then 5 km). Two complementary Q1 summaries are returned:

* ``secc3m`` = 1.7 / mean(K_490), representing background water clarity.
* ``secc3m_p10`` = 1.7 / Q90(K_490), representing acute low-clarity events.

The second metric is deliberately based on the same IMOS product in every
event year. It is an optical plume proxy, not a direct salinity measurement.
"""

import numpy as np
import pandas as pd
import xarray as xr


SECCHI_SOURCE = (
    "IMOS SRS MODIS-Aqua K_490 v202302; Secchi = 1.7 / Q1 mean(K_490)"
)


def haversine_km(lat0, lon0, lats, lons):
    """Great-circle distance from one point to arrays of points."""
    lat0_rad = np.radians(lat0)
    lats_rad = np.radians(lats)
    dlat = np.radians(lats - lat0)
    dlon = np.radians(lons - lon0)
    a = (
        np.sin(dlat / 2.0) ** 2
        + np.cos(lat0_rad) * np.cos(lats_rad) * np.sin(dlon / 2.0) ** 2
    )
    return 6371.0 * 2.0 * np.arcsin(np.sqrt(a))


def extract_q1_secchi(ds_k490, reefs, year, search_radii_km=(2.0, 3.0, 5.0)):
    """Return Secchi estimates and extraction diagnostics for a reef table.

    ``reefs`` must contain ``lon`` and ``lat``. Results preserve its row order.
    Longitude groups separated by more than one degree are read separately so
    focused repairs do not download large empty areas of the source grid.
    """
    missing_columns = {"lon", "lat"} - set(reefs.columns)
    if missing_columns:
        raise ValueError(f"Reef table lacks columns: {sorted(missing_columns)}")

    work = reefs[["lon", "lat"]].copy().reset_index(drop=True)
    work["source_row"] = np.arange(len(work))
    result = work.copy()
    result['k490_q90'] = np.nan
    result['secc3m_p10'] = np.nan
    result["secc3m"] = np.nan
    result["secc_source"] = SECCHI_SOURCE
    result["secc_match_method"] = "missing_no_valid_water"
    result["secc_radius_km"] = np.nan
    result["secc_n_observations"] = 0
    result["secc_n_grid_cells"] = 0
    result["secc_centroid_grid_distance_km"] = np.nan
    result["secc_nearest_valid_distance_km"] = np.nan

    if len(work) == 0:
        return result.drop(columns="source_row")

    ordered = work.sort_values("lon").copy()
    ordered["spatial_group"] = ordered["lon"].diff().gt(1.0).cumsum()
    max_radius_km = max(search_radii_km)
    padding_deg = max_radius_km / 111.0 + 0.02

    for _, group in ordered.groupby("spatial_group", sort=False):
        lat_hi = float(group["lat"].max() + padding_deg)
        lat_lo = float(group["lat"].min() - padding_deg)
        lon_lo = float(group["lon"].min() - padding_deg)
        lon_hi = float(group["lon"].max() + padding_deg)

        latitude_slice = (
            slice(lat_hi, lat_lo)
            if ds_k490.latitude.values[0] > ds_k490.latitude.values[-1]
            else slice(lat_lo, lat_hi)
        )
        window = ds_k490["K_490"].sel(
            time=slice(f"{year}-01-01", f"{year}-03-31"),
            latitude=latitude_slice,
            longitude=slice(lon_lo, lon_hi),
        )
        if window.sizes.get("time", 0) == 0:
            continue

        # Sums and counts weight all valid observations equally even when cloud
        # cover gives adjacent water pixels different sample sizes.
        positive = window.where(window > 0)
        stats = xr.Dataset(
            {
                "k490_sum": positive.sum("time", skipna=True),
                "n_obs": positive.count("time"),
            }
        ).compute()
        q90_grid = positive.quantile(
            0.90, dim='time', skipna=True
        ).compute()

        grid_lats = stats.latitude.values
        grid_lons = stats.longitude.values
        k490_sum = stats.k490_sum.values
        n_obs = stats.n_obs.values
        k490_q90 = q90_grid.values

        for reef in group.itertuples(index=False):
            row = int(reef.source_row)
            lat_idx = int(np.abs(grid_lats - reef.lat).argmin())
            lon_idx = int(np.abs(grid_lons - reef.lon).argmin())
            centroid_distance = float(
                haversine_km(
                    reef.lat,
                    reef.lon,
                    grid_lats[lat_idx],
                    grid_lons[lon_idx],
                )
            )
            result.loc[row, "secc_centroid_grid_distance_km"] = centroid_distance

            point_n = int(n_obs[lat_idx, lon_idx])
            if point_n > 0:
                point_mean = k490_sum[lat_idx, lon_idx] / point_n
                point_q90 = k490_q90[lat_idx, lon_idx]
                result.loc[row, 'k490_q90'] = point_q90
                result.loc[row, 'secc3m_p10'] = 1.7 / point_q90
                result.loc[row, "secc3m"] = 1.7 / point_mean
                result.loc[row, "secc_match_method"] = "centroid_grid_cell"
                result.loc[row, "secc_radius_km"] = 0.0
                result.loc[row, "secc_n_observations"] = point_n
                result.loc[row, "secc_n_grid_cells"] = 1
                result.loc[row, "secc_nearest_valid_distance_km"] = centroid_distance
                continue

            # Five kilometres is at most about six 0.01-degree grid cells in
            # either direction over the GBR. Work only on this local patch.
            lat_step = max(float(np.median(np.abs(np.diff(grid_lats)))), 1e-6)
            lon_step = max(float(np.median(np.abs(np.diff(grid_lons)))), 1e-6)
            lat_cells = int(np.ceil(max_radius_km / (111.0 * lat_step))) + 1
            lon_km_per_degree = 111.0 * np.cos(np.radians(reef.lat))
            lon_cells = int(
                np.ceil(max_radius_km / (lon_km_per_degree * lon_step))
            ) + 1
            lat_slice = slice(
                max(0, lat_idx - lat_cells),
                min(len(grid_lats), lat_idx + lat_cells + 1),
            )
            lon_slice = slice(
                max(0, lon_idx - lon_cells),
                min(len(grid_lons), lon_idx + lon_cells + 1),
            )

            local_lats, local_lons = np.meshgrid(
                grid_lats[lat_slice], grid_lons[lon_slice], indexing="ij"
            )
            distances = haversine_km(reef.lat, reef.lon, local_lats, local_lons)
            local_n = n_obs[lat_slice, lon_slice]
            local_sum = k490_sum[lat_slice, lon_slice]
            local_q90 = k490_q90[lat_slice, lon_slice]
            valid_water = local_n > 0
            if valid_water.any():
                result.loc[row, "secc_nearest_valid_distance_km"] = float(
                    distances[valid_water].min()
                )

            for radius_km in search_radii_km:
                selected = valid_water & (distances <= radius_km)
                selected_n = int(local_n[selected].sum())
                if selected_n == 0:
                    continue
                spatial_temporal_mean = float(local_sum[selected].sum() / selected_n)
                spatial_q90 = float(
                    np.average(local_q90[selected], weights=local_n[selected])
                )
                result.loc[row, 'k490_q90'] = spatial_q90
                result.loc[row, 'secc3m_p10'] = 1.7 / spatial_q90
                result.loc[row, "secc3m"] = 1.7 / spatial_temporal_mean
                result.loc[row, "secc_match_method"] = "water_neighbourhood"
                result.loc[row, "secc_radius_km"] = radius_km
                result.loc[row, "secc_n_observations"] = selected_n
                result.loc[row, "secc_n_grid_cells"] = int(selected.sum())
                break

    return (
        result.sort_values("source_row")
        .drop(columns="source_row")
        .reset_index(drop=True)
    )
