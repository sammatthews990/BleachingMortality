"""Reconstruct the Cheung et al. (2025) JFM PATMOS-x cloud predictor.

The paper defines ``cloudp_90`` as mean January--March cloud fraction from
NOAA-18 AVHRR PATMOS-x v6.0 at 0.1 degree resolution. This script streams only
the required HDF5 chunks from NOAA's public S3 archive, calculates several
explicit pass aggregations, samples the nearest PATMOS-x cell to each reef,
and writes a provenance-rich reef-year table.

The archived Cheung values may be supplied with ``--benchmark`` for validation,
but they are never used to construct or fill the predictor.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import datetime as dt
import hashlib
import json
import re
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
from xml.etree import ElementTree

import numpy as np
import pandas as pd
import requests
import s3fs


PRODUCT_VERSION = "v06r00"
PLATFORM = "NOAA-18"
VARIABLE = "cloud_fraction"
PAPER_YEARS = (2016, 2017, 2020)
DEFAULT_YEARS = (2016, 2017, 2020, 2022, 2024)
S3_LIST_URL = "https://noaa-cdr-patmosx-radiances-and-clouds-pds.s3.amazonaws.com/"
S3_BUCKET = "noaa-cdr-patmosx-radiances-and-clouds-pds"
REEF_REGISTRY_URL = "s3://gbr-dms-data-public/gbrmpa-complete-gbr-features/data.parquet"
PRODUCT_DOI = "https://doi.org/10.7289/V5X9287S"
PAPER_DOI = "https://doi.org/10.1111/geb.70105"
FILENAME_RE = re.compile(
    r"patmosx_(?P<version>v\d{2}r\d{2})(?P<preliminary>-preliminary)?_"
    r"(?P<platform>[^_]+)_(?P<orbit>asc|des)_d(?P<date>\d{8})_c(?P<created>\d{8})\.nc$"
)


@dataclass(frozen=True)
class PatmosFile:
    s3_key: str
    name: str
    date: dt.date
    orbit: str
    created: dt.date
    preliminary: bool


def request_with_retry(
    url: str,
    *,
    params: dict[str, object] | None = None,
    timeout: int = 90,
    attempts: int = 5,
) -> requests.Response:
    """GET a public NOAA resource with bounded exponential backoff."""
    last_error: Exception | None = None
    for attempt in range(attempts):
        try:
            response = requests.get(url, params=params, timeout=timeout)
            response.raise_for_status()
            return response
        except requests.RequestException as exc:
            last_error = exc
            if attempt + 1 < attempts:
                time.sleep(2**attempt)
    raise RuntimeError(f"NOAA request failed after {attempts} attempts: {url}") from last_error


def discover_files(year: int, platform: str = PLATFORM) -> list[PatmosFile]:
    """Discover and deterministically deduplicate v6 JFM files for one year."""
    candidates: list[PatmosFile] = []
    namespace = {"s3": "http://s3.amazonaws.com/doc/2006-03-01/"}
    # NOAA's S3 inventory is authoritative for v6; the legacy THREDDS catalog
    # advertises v5.3 even where the corresponding v6 files are queryable.
    prefixes = (
        f"data/{year}/patmosx_{PRODUCT_VERSION}_{platform}_",
        f"data/{year}/patmosx_{PRODUCT_VERSION}-preliminary_{platform}_",
    )
    for prefix in prefixes:
        continuation: str | None = None
        while True:
            params: dict[str, object] = {"list-type": "2", "prefix": prefix, "max-keys": 1000}
            if continuation:
                params["continuation-token"] = continuation
            response = request_with_retry(S3_LIST_URL, params=params, timeout=180)
            root = ElementTree.fromstring(response.content)
            for key_node in root.findall("s3:Contents/s3:Key", namespace):
                key = key_node.text or ""
                name = key.rsplit("/", 1)[-1]
                match = FILENAME_RE.fullmatch(name)
                if not match:
                    continue
                date = dt.datetime.strptime(match.group("date"), "%Y%m%d").date()
                if date.year != year or date.month not in (1, 2, 3):
                    continue
                candidates.append(
                    PatmosFile(
                        s3_key=key,
                        name=name,
                        date=date,
                        orbit=match.group("orbit"),
                        created=dt.datetime.strptime(match.group("created"), "%Y%m%d").date(),
                        preliminary=bool(match.group("preliminary")),
                    )
                )
            if root.findtext("s3:IsTruncated", default="false", namespaces=namespace) != "true":
                break
            continuation = root.findtext("s3:NextContinuationToken", namespaces=namespace)

    # Prefer a final record over preliminary, then the newest creation date.
    chosen: dict[tuple[dt.date, str], PatmosFile] = {}
    for item in candidates:
        key = (item.date, item.orbit)
        incumbent = chosen.get(key)
        rank = (not item.preliminary, item.created)
        if incumbent is None or rank > (not incumbent.preliminary, incumbent.created):
            chosen[key] = item

    return sorted(chosen.values(), key=lambda item: (item.date, item.orbit))


def decode_hdf5(raw: np.ndarray, dataset: object) -> np.ndarray:
    """Apply NetCDF fill, scale, and offset attributes to a raw HDF5 array."""
    values = np.asarray(raw, dtype=np.float32)
    fill = np.asarray(dataset.attrs.get("_FillValue", np.nan)).item()
    values[values == fill] = np.nan
    scale = np.asarray(dataset.attrs.get("scale_factor", 1.0)).item()
    offset = np.asarray(dataset.attrs.get("add_offset", 0.0)).item()
    return values * scale + offset


def stream_subset(
    item: PatmosFile,
    bbox: tuple[float, float, float, float],
    cache_dir: Path,
) -> Path:
    """Stream the required v6 HDF5 chunks and cache a small decoded GBR grid."""
    try:
        import h5py
    except ImportError as exc:
        raise RuntimeError(
            "PATMOS-x v6 streaming requires h5py. Install dependencies with "
            "`python -m pip install -r scripts/requirements-patmosx.txt`."
        ) from exc

    destination = cache_dir / str(item.date.year) / item.name.replace(".nc", ".npz")
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size > 1_000:
        return destination

    west, south, east, north = bbox
    temporary = destination.with_suffix(destination.suffix + ".part")
    last_error: Exception | None = None
    for attempt in range(4):
        try:
            filesystem = s3fs.S3FileSystem(anon=True, default_block_size=1024 * 1024)
            remote_path = f"{S3_BUCKET}/{item.s3_key}"
            with filesystem.open(
                remote_path,
                "rb",
                block_size=1024 * 1024,
                cache_type="readahead",
            ) as remote:
                with h5py.File(remote, "r") as source:
                    latitude = decode_hdf5(source["latitude"][:], source["latitude"])
                    longitude = decode_hdf5(source["longitude"][:], source["longitude"])
                    lat_index = np.flatnonzero((latitude >= south) & (latitude <= north))
                    lon_index = np.flatnonzero((longitude >= west) & (longitude <= east))
                    if len(lat_index) == 0 or len(lon_index) == 0:
                        raise ValueError(f"Bounding box {bbox} does not intersect the PATMOS-x grid")
                    lat_slice = slice(int(lat_index[0]), int(lat_index[-1]) + 1)
                    lon_slice = slice(int(lon_index[0]), int(lon_index[-1]) + 1)
                    cloud = source[VARIABLE]
                    values = decode_hdf5(cloud[0, lat_slice, lon_slice], cloud)
                    latitude = latitude[lat_slice]
                    longitude = longitude[lon_slice]
            with temporary.open("wb") as handle:
                np.savez_compressed(
                    handle,
                    cloud_fraction=values.astype(np.float32),
                    latitude=latitude.astype(np.float64),
                    longitude=longitude.astype(np.float64),
                )
            temporary.replace(destination)
            return destination
        except Exception as exc:
            last_error = exc
            temporary.unlink(missing_ok=True)
            if attempt + 1 < 4:
                time.sleep(2**attempt)
    raise RuntimeError(f"Failed to stream {item.s3_key}") from last_error


def read_grid(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read a decoded cloud grid and its coordinates."""
    with np.load(path) as dataset:
        values = np.asarray(dataset[VARIABLE], dtype=np.float32)
        latitude = np.asarray(dataset["latitude"], dtype=np.float64)
        longitude = np.asarray(dataset["longitude"], dtype=np.float64)
    if values.ndim != 2:
        raise ValueError(f"Expected a two-dimensional cloud grid in {path}, got {values.shape}")
    return values, latitude, longitude


def nanmean_stack(grids: Iterable[np.ndarray]) -> np.ndarray:
    """Average grids without emitting all-missing slice warnings."""
    stack = np.stack(list(grids)).astype(np.float32, copy=False)
    valid = np.isfinite(stack)
    count = valid.sum(axis=0)
    total = np.where(valid, stack, 0.0).sum(axis=0, dtype=np.float64)
    return np.divide(
        total,
        count,
        out=np.full(total.shape, np.nan, dtype=np.float64),
        where=count > 0,
    ).astype(np.float32)


def aggregate_year(
    files: list[PatmosFile],
    downloaded: dict[str, Path],
) -> tuple[dict[str, np.ndarray], np.ndarray, np.ndarray, dict[str, object]]:
    """Calculate pass-specific, pass-pooled, and equal-day JFM means."""
    arrays: dict[tuple[dt.date, str], np.ndarray] = {}
    latitude: np.ndarray | None = None
    longitude: np.ndarray | None = None
    preliminary_count = 0

    for item in files:
        grid, item_latitude, item_longitude = read_grid(downloaded[item.name])
        if latitude is None:
            latitude, longitude = item_latitude, item_longitude
        elif not (np.array_equal(latitude, item_latitude) and np.array_equal(longitude, item_longitude)):
            raise ValueError("PATMOS-x subset grids are not aligned")
        arrays[(item.date, item.orbit)] = grid
        preliminary_count += int(item.preliminary)

    if latitude is None or longitude is None:
        raise ValueError("No PATMOS-x arrays were available for aggregation")

    asc = [grid for (date, orbit), grid in arrays.items() if orbit == "asc"]
    des = [grid for (date, orbit), grid in arrays.items() if orbit == "des"]
    daily = []
    dates = sorted({date for date, _ in arrays})
    paired_days = 0
    for date in dates:
        day_grids = [arrays[(date, orbit)] for orbit in ("asc", "des") if (date, orbit) in arrays]
        paired_days += int(len(day_grids) == 2)
        daily.append(nanmean_stack(day_grids))

    grids = {
        "cloud_fraction_asc_jfm": nanmean_stack(asc),
        "cloud_fraction_des_jfm": nanmean_stack(des),
        # Each available overpass has equal weight.
        "cloud_fraction_passmean_jfm": nanmean_stack(arrays.values()),
        # Each calendar day has equal weight after within-day pass averaging.
        "cloud_fraction_dailymean_jfm": nanmean_stack(daily),
    }
    metadata = {
        "n_files": len(files),
        "n_asc_files": len(asc),
        "n_des_files": len(des),
        "n_days": len(dates),
        "n_paired_days": paired_days,
        "n_preliminary_files": preliminary_count,
        "first_date": min(dates).isoformat(),
        "last_date": max(dates).isoformat(),
    }
    return grids, latitude, longitude, metadata


def load_reefs(path: Path | None) -> pd.DataFrame:
    """Load a reef registry and normalize its coordinate column names."""
    if path is None:
        reefs = pd.read_parquet(REEF_REGISTRY_URL, storage_options={"anon": True})
        if "LEVEL_1" in reefs.columns:
            reefs = reefs.loc[reefs["LEVEL_1"] == "Reef"].copy()
        source_label = REEF_REGISTRY_URL
    else:
        reefs = pd.read_csv(path)
        source_label = str(path)
    aliases = {
        "lon": ("lon", "longitude", "X_COORD", "X"),
        "lat": ("lat", "latitude", "Y_COORD", "Y"),
    }
    rename: dict[str, str] = {}
    for target, candidates in aliases.items():
        source = next((name for name in candidates if name in reefs.columns), None)
        if source is None:
            raise ValueError(f"{source_label} needs one of these {target} columns: {candidates}")
        rename[source] = target
    reefs = reefs.rename(columns=rename)
    required = {"LABEL_ID", "lon", "lat"}
    missing = required.difference(reefs.columns)
    if missing:
        raise ValueError(f"Missing reef columns: {sorted(missing)}")
    keep = ["LABEL_ID"] + (["LOC_NAME_S"] if "LOC_NAME_S" in reefs.columns else []) + ["lon", "lat"]
    reefs = reefs[keep].copy()
    reefs["LABEL_ID"] = reefs["LABEL_ID"].astype(str).str.strip()
    for column in ("lon", "lat"):
        reefs[column] = pd.to_numeric(reefs[column], errors="raise")
    reefs = reefs.drop_duplicates(["LABEL_ID", "lon", "lat"]).copy()
    if reefs.duplicated(["LABEL_ID", "lon", "lat"]).any():
        raise ValueError("Reef registry has duplicate label-coordinate records")
    return reefs


def nearest_indices(values: np.ndarray, targets: np.ndarray) -> np.ndarray:
    return np.abs(values[:, None] - targets[None, :]).argmin(axis=0)


def sample_reefs(
    reefs: pd.DataFrame,
    year: int,
    grids: dict[str, np.ndarray],
    latitude: np.ndarray,
    longitude: np.ndarray,
    metadata: dict[str, object],
) -> pd.DataFrame:
    """Sample nearest 0.1 degree cells and attach source/coverage fields."""
    lat_index = nearest_indices(latitude, reefs["lat"].to_numpy())
    lon_index = nearest_indices(longitude, reefs["lon"].to_numpy())
    # The archived paper values reveal the extraction convention used by its
    # raster workflow: nearest latitude and the first longitude cell centre at
    # or east of the reef. Retain conventional nearest-neighbour values too.
    paper_lon_index = np.clip(
        np.searchsorted(longitude, reefs["lon"].to_numpy(), side="left"),
        0,
        len(longitude) - 1,
    )
    result = reefs.copy()
    result.insert(1, "year", year)
    result["patmos_lat"] = latitude[lat_index]
    result["patmos_lon"] = longitude[lon_index]
    result["patmos_lon_paper_cell"] = longitude[paper_lon_index]
    for name, grid in grids.items():
        result[name] = grid[lat_index, lon_index]
        result[f"{name}_paper_cell"] = grid[lat_index, paper_lon_index]
    # Archived-value validation identifies the ascending pass as the paper's
    # operational definition (2017 r=0.9955; daily two-pass mean r=0.9619).
    result["cloudp_90"] = result["cloud_fraction_asc_jfm_paper_cell"]
    result["cloud_n_days"] = metadata["n_days"]
    result["cloud_n_asc"] = metadata["n_asc_files"]
    result["cloud_n_des"] = metadata["n_des_files"]
    result["cloud_n_paired_days"] = metadata["n_paired_days"]
    result["cloud_n_preliminary_files"] = metadata["n_preliminary_files"]
    result["cloud_platform"] = PLATFORM
    result["cloud_product_version"] = PRODUCT_VERSION
    result["cloud_aggregation"] = "mean(Jan-Mar NOAA-18 ascending pass)"
    result["cloud_spatial_match"] = "paper convention: nearest latitude, east-containing longitude cell"
    return result


def validation_metrics(recreated: pd.DataFrame, benchmark_path: Path) -> pd.DataFrame:
    """Compare candidate aggregations with archived values; never use them as inputs."""
    benchmark = pd.read_csv(benchmark_path, usecols=["LABEL_ID", "year", "cloudp_90"])
    benchmark["LABEL_ID"] = benchmark["LABEL_ID"].astype(str).str.strip()
    benchmark = benchmark.rename(columns={"cloudp_90": "cloudp_90_archived"})
    merged = recreated.merge(benchmark, on=["LABEL_ID", "year"], how="inner")
    records = []
    candidates = [
        "cloud_fraction_asc_jfm",
        "cloud_fraction_des_jfm",
        "cloud_fraction_passmean_jfm",
        "cloud_fraction_dailymean_jfm",
    ]
    candidates += [f"{candidate}_paper_cell" for candidate in candidates]
    for (year,), group in merged.groupby(["year"]):
        for candidate in candidates:
            valid = group[[candidate, "cloudp_90_archived"]].dropna()
            delta = valid[candidate] - valid["cloudp_90_archived"]
            records.append(
                {
                    "year": int(year),
                    "candidate": candidate,
                    "n": len(valid),
                    "correlation": valid[candidate].corr(valid["cloudp_90_archived"]),
                    "bias": delta.mean(),
                    "mae": delta.abs().mean(),
                    "rmse": np.sqrt(np.mean(np.square(delta))),
                }
            )
    return pd.DataFrame.from_records(records)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reefs",
        type=Path,
        default=None,
        help="Optional reef CSV; defaults to the public GBR DMS reef-feature registry.",
    )
    parser.add_argument("--years", nargs="+", type=int, default=list(DEFAULT_YEARS))
    parser.add_argument("--output", type=Path, default=Path("data/processed/patmosx_cloud_reef_year.csv"))
    parser.add_argument("--metadata", type=Path, default=Path("data/processed/patmosx_cloud_metadata.json"))
    parser.add_argument("--cache", type=Path, default=Path("data/cache/patmosx_s3_gbr"))
    parser.add_argument(
        "--benchmark",
        type=Path,
        default=None,
        help="Optional archived cloud table used only for validation, never construction.",
    )
    parser.add_argument("--validation-output", type=Path, default=Path("data/processed/patmosx_cloud_validation.csv"))
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--min-days", type=int, default=60)
    parser.add_argument("--bbox", nargs=4, type=float, metavar=("WEST", "SOUTH", "EAST", "NORTH"), default=(142.0, -25.0, 155.0, -9.0))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    reefs = load_reefs(args.reefs)
    frames = []
    run_metadata: dict[str, object] = {
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "paper": PAPER_DOI,
        "product": PRODUCT_DOI,
        "archive": S3_LIST_URL,
        "product_version": PRODUCT_VERSION,
        "platform": PLATFORM,
        "variable": VARIABLE,
        "temporal_window": "January 1 through March 31, inclusive",
        "primary_aggregation": "equal-weight mean of available NOAA-18 ascending passes",
        "spatial_resolution_degrees": 0.1,
        "spatial_sampling": "nearest cell centre",
        "bbox_west_south_east_north": list(args.bbox),
        "reef_source": str(args.reefs) if args.reefs is not None else REEF_REGISTRY_URL,
        "years": {},
    }

    for year in args.years:
        print(f"Discovering {year} {PLATFORM} {PRODUCT_VERSION} files...")
        files = discover_files(year)
        unique_days = {item.date for item in files}
        if len(unique_days) < args.min_days:
            raise RuntimeError(
                f"{year} has only {len(unique_days)} discoverable JFM days; "
                f"minimum is {args.min_days}. No fallback will be substituted."
            )
        print(f"Downloading/caching {len(files)} GBR subsets across {len(unique_days)} days...")
        downloaded: dict[str, Path] = {}
        # h5py serializes threads through its global library lock. Independent
        # processes provide real parallel range reads and isolate failed files.
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as executor:
            future_map = {
                executor.submit(stream_subset, item, tuple(args.bbox), args.cache): item
                for item in files
            }
            for completed, future in enumerate(concurrent.futures.as_completed(future_map), start=1):
                item = future_map[future]
                downloaded[item.name] = future.result()
                if completed % 25 == 0 or completed == len(files):
                    print(f"  {year}: {completed}/{len(files)} subsets ready")

        grids, latitude, longitude, year_metadata = aggregate_year(files, downloaded)
        frames.append(sample_reefs(reefs, year, grids, latitude, longitude, year_metadata))
        run_metadata["years"][str(year)] = year_metadata

    result = pd.concat(frames, ignore_index=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)
    run_metadata["output"] = str(args.output)
    run_metadata["output_rows"] = len(result)
    run_metadata["output_sha256"] = sha256(args.output)

    if args.benchmark is not None and args.benchmark.exists():
        validation = validation_metrics(result, args.benchmark)
        validation.to_csv(args.validation_output, index=False)
        run_metadata["validation_output"] = str(args.validation_output)
        print("\nValidation against archived Cheung values (validation only):")
        print(validation.to_string(index=False, float_format=lambda value: f"{value:.6f}"))
    else:
        run_metadata["validation_output"] = None
        print(f"Benchmark not found at {args.benchmark}; reconstruction was still completed.")

    args.metadata.write_text(json.dumps(run_metadata, indent=2), encoding="utf-8")
    print(f"\nWrote {args.output} and {args.metadata}")


if __name__ == "__main__":
    main()
