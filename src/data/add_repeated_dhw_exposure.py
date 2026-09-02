"""Add leakage-safe repeated-DHW exposure metrics to the model table.

This local enrichment avoids re-downloading the large remote environmental
stores when the annual reef DHW archive is already present. It updates the
canonical CSV in place, validates it, and writes a compact reef-year audit.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from validate_environmental_data import validate_environmental_file


ENVIRONMENT = Path("data/processed/cheung_recreated_gbr_full.csv")
AUDIT = Path("data/processed/repeated_dhw_exposure_audit.csv")
HISTORY_START = 1985
HISTORY_END = 2023


def extract_annual_history(features: pd.DataFrame) -> tuple[list[int], np.ndarray]:
    """Extract the same nearest-cell NOAA CRW annual maxima as the main build."""
    dataset = xr.open_zarr(
        "s3://gbr-dms-data-public/noaa-crw-chs-dhw/data.zarr",
        storage_options={"anon": True},
    )
    reef_lats = xr.DataArray(features["lat"].to_numpy(), dims="reef")
    reef_lons = xr.DataArray(features["lon"].to_numpy(), dims="reef")
    years = list(range(HISTORY_START, HISTORY_END + 1))
    annual = []
    for year in years:
        values = dataset["degree_heating_week"].sel(
            time=slice(f"{year}-01-01", f"{year}-12-31")
        ).sel(lat=reef_lats, lon=reef_lons, method="nearest")
        annual.append(values.max(dim="time").compute().values)
        print(f"Extracted annual maximum DHW for {year}", flush=True)
    return years, np.asarray(annual, dtype=float)


def main() -> None:
    environment = pd.read_csv(ENVIRONMENT)
    feature_keys = ["LABEL_ID", "lon", "lat"]
    features = environment[feature_keys].drop_duplicates().reset_index(drop=True)
    years, annual = extract_annual_history(features)
    years_array = np.asarray(years)

    metric_frames = []
    for target_year in sorted(environment["year"].unique()):
        prior = annual[years_array < target_year, :]
        prior_years = years_array[years_array < target_year]
        fixed = prior[prior_years >= 2016, :]
        rolling = prior[prior_years >= target_year - 8, :]
        current = features.copy()
        current["year"] = int(target_year)
        current["dhw_events_since2016_n6"] = np.sum(fixed > 6, axis=0)
        current["dhw_events_prior8_n6"] = np.sum(rolling > 6, axis=0)
        current["dhw_history_years_since2016"] = np.sum(
            np.isfinite(fixed), axis=0
        )
        event = np.isfinite(prior) & (prior > 6)
        years_since = np.full(len(features), np.nan)
        no_prior = np.ones(len(features))
        for reef_index in range(len(features)):
            event_indices = np.where(event[:, reef_index])[0]
            if len(event_indices) > 0:
                years_since[reef_index] = (
                    target_year - prior_years[event_indices[-1]]
                )
                no_prior[reef_index] = 0
        current["dhw_years_since_last_n6"] = years_since
        current["dhw_years_since_last_n6_capped8"] = np.where(
            np.isfinite(years_since), np.minimum(years_since, 8), 8
        )
        current["dhw_no_prior_n6"] = no_prior
        metric_frames.append(current)
    metrics = pd.concat(metric_frames, ignore_index=True)

    new_columns = [
        "dhw_events_since2016_n6",
        "dhw_events_prior8_n6",
        "dhw_history_years_since2016",
        "dhw_years_since_last_n6",
        "dhw_years_since_last_n6_capped8",
        "dhw_no_prior_n6",
    ]
    environment = environment.drop(columns=new_columns, errors="ignore").merge(
        metrics, on=[*feature_keys, "year"], how="left", validate="one_to_one"
    )
    complete_columns = [
        "dhw_events_since2016_n6", "dhw_events_prior8_n6",
        "dhw_history_years_since2016", "dhw_no_prior_n6",
    ]
    if environment[complete_columns].isna().any().any():
        raise ValueError("Repeated-exposure enrichment produced missing values")

    # The current event must never contribute to a history predictor.
    check_2016 = environment.loc[environment["year"] == 2016]
    if not (check_2016["dhw_events_since2016_n6"] == 0).all():
        raise AssertionError("2016 history includes the response year")

    environment.to_csv(ENVIRONMENT, index=False)
    validate_environmental_file(
        ENVIRONMENT, ENVIRONMENT.with_suffix(".manifest.json")
    )

    audit = environment[[
        "LABEL_ID", "LOC_NAME_S", "year", "ann_maxdhw", "yrsince6",
        *new_columns,
    ]].copy()
    audit.to_csv(AUDIT, index=False)
    print(
        f"Added repeated exposure to {len(environment):,} reef-years; "
        f"audit: {AUDIT}"
    )


if __name__ == "__main__":
    main()
