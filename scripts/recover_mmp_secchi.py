"""Recover missing MMP Secchi estimates from nearby valid IMOS water pixels.

This focused repair only replaces missing MMP reef-years in the current
processed table. It writes an audit table containing the source radius and
observation count. Full future rebuilds use the same extraction rule in
``fetch_dms_environmental_data.py``.
"""

from pathlib import Path

import pandas as pd
import xarray as xr

from secchi_extraction import SECCHI_SOURCE, extract_q1_secchi
from validate_environmental_data import validate_environmental_file


ENVIRONMENT_FILE = Path("data/processed/cheung_recreated_gbr_full.csv")
MMP_ROWS_FILE = Path("data/processed/validation_rows_mmp.csv")
AUDIT_FILE = Path("data/processed/mmp_secchi_recovery_audit.csv")
MANIFEST_FILE = Path("data/processed/cheung_recreated_gbr_full.manifest.json")


environment = pd.read_csv(ENVIRONMENT_FILE)
mmp_rows = pd.read_csv(MMP_ROWS_FILE)

targets = (
    mmp_rows.loc[
        mmp_rows["secc3m"].isna(),
        ["event_year", "ReefID", "ReefName", "lon", "lat"],
    ]
    .drop_duplicates()
    .sort_values(["event_year", "ReefID"])
    .reset_index(drop=True)
)
if targets.empty:
    print("No missing MMP Secchi values require recovery.")
    raise SystemExit(0)

ds_k490 = xr.open_zarr(
    "s3://gbr-dms-data-public/imos-srs-aqua-oc-k490/data.zarr",
    storage_options={"anon": True},
)

audit_parts = []
for year, year_targets in targets.groupby("event_year", sort=True):
    estimates = extract_q1_secchi(
        ds_k490,
        year_targets[["lon", "lat"]],
        int(year),
    )
    audit_parts.append(
        pd.concat(
            [
                year_targets.reset_index(drop=True),
                estimates.drop(columns=["lon", "lat"]),
            ],
            axis=1,
        )
    )

audit = pd.concat(audit_parts, ignore_index=True)
if audit["secc3m"].isna().any():
    failed = audit.loc[audit["secc3m"].isna(), ["event_year", "ReefID"]]
    failed_labels = ", ".join(
        f"{row.ReefID}/{row.event_year}" for row in failed.itertuples()
    )
    raise RuntimeError(f"No valid IMOS water pixels within 5 km for: {failed_labels}")

metadata_columns = {
    "secc_source": None,
    "secc_match_method": None,
    "secc_radius_km": float("nan"),
    "secc_n_observations": float("nan"),
    "secc_n_grid_cells": float("nan"),
    "secc_centroid_grid_distance_km": float("nan"),
    "secc_nearest_valid_distance_km": float("nan"),
}
for column, default in metadata_columns.items():
    if column not in environment.columns:
        environment[column] = default

existing = environment["secc3m"].notna()
environment.loc[
    existing & environment["secc_source"].isna(), "secc_source"
] = SECCHI_SOURCE
environment.loc[
    existing & environment["secc_match_method"].isna(), "secc_match_method"
] = "centroid_grid_cell_legacy"
environment.loc[
    existing & environment["secc_radius_km"].isna(), "secc_radius_km"
] = 0.0
environment.loc[
    existing & environment["secc_n_grid_cells"].isna(), "secc_n_grid_cells"
] = 1.0

patched_rows = 0
for recovered in audit.itertuples(index=False):
    selected = (
        (environment["year"] == int(recovered.event_year))
        & (environment["LABEL_ID"] == recovered.ReefID)
        & environment["secc3m"].isna()
    )
    if selected.sum() == 0:
        raise RuntimeError(
            f"No missing environment row matched "
            f"{recovered.ReefID}/{recovered.event_year}"
        )
    environment.loc[selected, "secc3m"] = recovered.secc3m
    for column in metadata_columns:
        environment.loc[selected, column] = getattr(recovered, column)
    patched_rows += int(selected.sum())

environment.to_csv(ENVIRONMENT_FILE, index=False)
audit.to_csv(AUDIT_FILE, index=False)
validate_environmental_file(ENVIRONMENT_FILE, MANIFEST_FILE)

print(
    f"Recovered {len(audit)} MMP reef-year estimates and patched "
    f"{patched_rows} environment rows."
)
print(
    audit[
        [
            "event_year",
            "ReefID",
            "secc3m",
            "secc_radius_km",
            "secc_n_observations",
            "secc_nearest_valid_distance_km",
        ]
    ].to_string(index=False)
)
