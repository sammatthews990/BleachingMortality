"""Validate the GBR reef-year environmental modelling table and write a manifest.

This is intentionally a small boundary check, not a workflow framework. It
defines the schema and scientific sanity checks required before the table is
consumed by the exploratory and modelling notebooks.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


CONTRACT_VERSION = "1.0"
EXPECTED_YEARS = [2016, 2017, 2020, 2022, 2024]
KEY_COLUMNS = ["LABEL_ID", "lon", "lat", "year"]
REQUIRED_COLUMNS = {
    "year",
    "LABEL_ID",
    "LOC_NAME_S",
    "lon",
    "lat",
    "ann_maxdhw",
    "histmDHW6",
    "yrsince6",
    "histmDHW4",
    "yrsince4",
    "ann_maxsst",
    "winyear_mean",
    "winyear_sd",
    "mcur_90",
    "dist_to_er_km",
    "secc3m",
    "cloudp_90",
    "cloud_n_days",
    "cloud_n_asc",
    "cloud_n_des",
    "cloud_n_paired_days",
    "cloud_n_preliminary_files",
    "cloud_platform",
    "cloud_product_version",
    "cloud_aggregation",
    "cloud_spatial_match",
    'dhw10_mean',
    'dhw10_max',
    'dhw10_load4',
    'dhw10_n4',
    'dhw10_n6',
    'dhw_novelty10',
    'k490_q90',
    'secc3m_p10',
}

# Bounds are broad sanity limits, not ecological filtering thresholds.
VALID_RANGES = {
    "lon": (140.0, 155.0),
    "lat": (-26.0, -8.0),
    "ann_maxdhw": (0.0, 50.0),
    "histmDHW6": (0.0, 50.0),
    "histmDHW4": (0.0, 50.0),
    "yrsince6": (0.0, 100.0),
    "yrsince4": (0.0, 100.0),
    "ann_maxsst": (15.0, 40.0),
    "winyear_mean": (15.0, 35.0),
    "winyear_sd": (0.0, 10.0),
    "mcur_90": (0.0, 5.0),
    "dist_to_er_km": (0.0, 500.0),
    "secc3m": (0.0, 100.0),
    "cloudp_90": (0.0, 1.0),
    "cloud_n_days": (60.0, 92.0),
    "cloud_n_asc": (60.0, 92.0),
    "cloud_n_des": (60.0, 92.0),
    "cloud_n_paired_days": (60.0, 92.0),
    "cloud_n_preliminary_files": (0.0, 184.0),
    'dhw10_mean': (0.0, 50.0),
    'dhw10_max': (0.0, 50.0),
    'dhw10_load4': (0.0, 500.0),
    'dhw10_n4': (0.0, 10.0),
    'dhw10_n6': (0.0, 10.0),
    'dhw_novelty10': (-50.0, 50.0),
    'k490_q90': (0.0, 20.0),
    'secc3m_p10': (0.0, 100.0),
}

# Known remote-sensing edge gaps are tolerated but quantified in the manifest.
MAX_MISSING_FRACTION = {
    'dhw10_mean': 0.01,
    'dhw10_max': 0.01,
    'dhw10_load4': 0.01,
    'dhw10_n4': 0.01,
    'dhw10_n6': 0.01,
    'dhw_novelty10': 0.01,
    'k490_q90': 0.10,
    'secc3m_p10': 0.10,
    "ann_maxdhw": 0.01,
    "ann_maxsst": 0.01,
    "winyear_mean": 0.01,
    "winyear_sd": 0.01,
    "secc3m": 0.10,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_environmental_dataframe(frame: pd.DataFrame) -> dict[str, object]:
    """Return a serialisable contract report with all errors collected."""
    errors: list[str] = []
    warnings: list[str] = []
    missing_columns = sorted(REQUIRED_COLUMNS.difference(frame.columns))
    if missing_columns:
        errors.append(f"Missing required columns: {missing_columns}")
        return {
            "contract_version": CONTRACT_VERSION,
            "valid": False,
            "errors": errors,
            "warnings": warnings,
            "rows": int(len(frame)),
            "columns": list(frame.columns),
        }

    data = frame.copy()
    data["LABEL_ID"] = data["LABEL_ID"].astype(str).str.strip().str.upper()
    data["year"] = pd.to_numeric(data["year"], errors="coerce")
    for coordinate in ("lon", "lat"):
        data[coordinate] = pd.to_numeric(data[coordinate], errors="coerce").round(6)

    observed_years = sorted(int(value) for value in data["year"].dropna().unique())
    if observed_years != EXPECTED_YEARS:
        errors.append(f"Expected years {EXPECTED_YEARS}, found {observed_years}")
    if data[KEY_COLUMNS].isna().any().any():
        errors.append("Reef-feature-year keys contain missing values")
    duplicate_keys = int(data.duplicated(KEY_COLUMNS).sum())
    if duplicate_keys:
        errors.append(f"Found {duplicate_keys} duplicate reef-feature-year keys")

    feature_columns = ["LABEL_ID", "lon", "lat"]
    feature_year_counts = data.groupby(feature_columns, dropna=False)["year"].nunique()
    incomplete_features = int((feature_year_counts != len(EXPECTED_YEARS)).sum())
    if incomplete_features:
        errors.append(f"Found {incomplete_features} reef features without all expected years")

    missingness: dict[str, dict[str, float | int]] = {}
    for column in sorted(REQUIRED_COLUMNS):
        count = int(data[column].isna().sum())
        fraction = float(count / len(data)) if len(data) else 1.0
        missingness[column] = {"count": count, "fraction": fraction}
        allowed = MAX_MISSING_FRACTION.get(column, 0.0)
        if fraction > allowed:
            errors.append(
                f"{column} missing fraction {fraction:.4f} exceeds allowed {allowed:.4f}"
            )
        elif count:
            warnings.append(f"{column} has {count} missing values ({fraction:.2%})")

    ranges: dict[str, dict[str, float | int | None]] = {}
    for column, (lower, upper) in VALID_RANGES.items():
        numeric = pd.to_numeric(data[column], errors="coerce")
        valid = numeric.dropna()
        below = int((valid < lower).sum())
        above = int((valid > upper).sum())
        ranges[column] = {
            "min": float(valid.min()) if len(valid) else None,
            "max": float(valid.max()) if len(valid) else None,
            "below_range": below,
            "above_range": above,
        }
        if below or above:
            errors.append(
                f"{column} has {below} values below {lower} and {above} above {upper}"
            )

    if not data["cloud_platform"].dropna().eq("NOAA-18").all():
        errors.append("cloud_platform contains values other than NOAA-18")
    if not data["cloud_product_version"].dropna().eq("v06r00").all():
        errors.append("cloud_product_version contains values other than v06r00")

    year_summary = {}
    for year, group in data.groupby("year"):
        year_summary[str(int(year))] = {
            "rows": int(len(group)),
            "features": int(group[feature_columns].drop_duplicates().shape[0]),
            "cloud_n_days": sorted(int(value) for value in group["cloud_n_days"].dropna().unique()),
            "cloud_n_preliminary_files": sorted(
                int(value) for value in group["cloud_n_preliminary_files"].dropna().unique()
            ),
        }

    return {
        "contract_version": CONTRACT_VERSION,
        "valid": not errors,
        "errors": errors,
        "warnings": warnings,
        "rows": int(len(data)),
        "reef_features": int(data[feature_columns].drop_duplicates().shape[0]),
        "years": observed_years,
        "duplicate_keys": duplicate_keys,
        "incomplete_features": incomplete_features,
        "columns": list(data.columns),
        "missingness": missingness,
        "ranges": ranges,
        "year_summary": year_summary,
    }


def validate_environmental_file(input_path: Path, manifest_path: Path) -> dict[str, object]:
    data = pd.read_csv(input_path)
    report = validate_environmental_dataframe(data)
    report.update(
        {
            "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "input": str(input_path),
            "input_sha256": sha256(input_path),
        }
    )
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    if not report["valid"]:
        raise ValueError("Environmental data contract failed:\n- " + "\n- ".join(report["errors"]))
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/processed/cheung_recreated_gbr_full.csv"),
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("data/processed/cheung_recreated_gbr_full.manifest.json"),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    report = validate_environmental_file(args.input, args.manifest)
    print(
        f"Contract passed: {report['rows']} rows, {report['reef_features']} reef features, "
        f"years {report['years']}. Manifest: {args.manifest}"
    )
    for warning in report["warnings"]:
        print(f"Warning: {warning}")


if __name__ == "__main__":
    main()
