"""Download and summarise ENSO indices for GBR bleaching-event summers."""

from __future__ import annotations

import argparse
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd


RONI_URL = "https://www.cpc.ncep.noaa.gov/data/indices/RONI.ascii.txt"
SOI_URL = "https://www.bom.gov.au/clim_data/IDCKGSH000/soi_monthly.txt"
EVENT_YEARS = (2016, 2017, 2020, 2022, 2024)
SUMMER_SEASONS = ("DJF", "JFM", "FMA")


def enso_category(value: float) -> tuple[str, str, str]:
    """Return phase, strength, and a readable category from continuous RONI."""
    magnitude = abs(value)
    if magnitude < 0.5:
        return "Neutral", "Neutral", "Neutral"
    phase = "El Nino" if value > 0 else "La Nina"
    if magnitude < 1.0:
        strength = "Weak"
    elif magnitude < 1.5:
        strength = "Moderate"
    else:
        strength = "Strong"
    return phase, strength, f"{strength} {phase}"


def soi_phase(value: float) -> str:
    """Classify mean Troup SOI using BoM's sustained +/-7 thresholds."""
    if value < -7:
        return "El Nino-like"
    if value > 7:
        return "La Nina-like"
    return "Neutral/mixed"


def agreement_label(roni_phase: str, atmospheric_phase: str) -> str:
    expected = {
        "El Nino": "El Nino-like",
        "La Nina": "La Nina-like",
        "Neutral": "Neutral/mixed",
    }
    return "coupled" if expected[roni_phase] == atmospheric_phase else "mixed"


def download(url: str, destination: Path, refresh: bool) -> None:
    if destination.exists() and not refresh:
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = Request(url, headers={"User-Agent": "GBR-mortality-pipeline/1.0"})
    with urlopen(request, timeout=90) as response:
        payload = response.read()
    if len(payload) < 100:
        raise RuntimeError(f"Unexpectedly short response from {url}")
    destination.write_bytes(payload)


def read_roni(path: Path) -> pd.DataFrame:
    data = pd.read_csv(path, sep=r"\s+")
    data.columns = [column.lower() for column in data.columns]
    data = data.rename(columns={"yr": "year", "anom": "roni"})
    required = {"seas", "year", "roni"}
    if not required.issubset(data.columns):
        raise ValueError(f"RONI file lacks columns: {sorted(required)}")
    return data


def read_soi(path: Path) -> pd.DataFrame:
    data = pd.read_csv(path, names=["year_month", "soi"])
    data["date"] = pd.to_datetime(
        data["year_month"].astype(str), format="%Y%m", errors="coerce"
    )
    data["soi"] = pd.to_numeric(data["soi"], errors="coerce")
    return data.dropna(subset=["date", "soi"])


def summarise_events(
    roni: pd.DataFrame, soi: pd.DataFrame, event_years: tuple[int, ...]
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for year in event_years:
        summer = roni.loc[
            (roni["year"] == year) & roni["seas"].isin(SUMMER_SEASONS)
        ].set_index("seas")["roni"]
        if set(summer.index) != set(SUMMER_SEASONS):
            raise ValueError(f"Incomplete RONI summer seasons for {year}")

        start = pd.Timestamp(year=year - 1, month=12, day=1)
        end = pd.Timestamp(year=year, month=3, day=31)
        soi_summer = soi.loc[soi["date"].between(start, end)]
        if len(soi_summer) != 4:
            raise ValueError(f"Incomplete December-March SOI months for {year}")

        roni_mean = float(summer.loc[list(SUMMER_SEASONS)].mean())
        phase, strength, category = enso_category(roni_mean)
        soi_mean = float(soi_summer["soi"].mean())
        atmosphere = soi_phase(soi_mean)
        rows.append(
            {
                "event_year": year,
                "roni_djf": float(summer["DJF"]),
                "roni_jfm": float(summer["JFM"]),
                "roni_fma": float(summer["FMA"]),
                "roni_bleaching_summer_mean": roni_mean,
                "roni_bleaching_summer_max_abs": float(
                    summer.iloc[np.argmax(np.abs(summer.to_numpy()))]
                ),
                "enso_phase": phase,
                "enso_strength": strength,
                "enso_category": category,
                "soi_dec_mar_mean": soi_mean,
                "soi_dec_mar_min": float(soi_summer["soi"].min()),
                "soi_dec_mar_max": float(soi_summer["soi"].max()),
                "soi_phase": atmosphere,
                "enso_ocean_atmosphere_state": agreement_label(phase, atmosphere),
                "primary_model_candidate": "roni_bleaching_summer_mean",
                "categorical_use": "interpretation_only",
                "roni_source": RONI_URL,
                "soi_source": SOI_URL,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--raw-dir", type=Path, default=Path("data/raw/enso"))
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/processed/enso_event_context.csv"),
    )
    args = parser.parse_args()

    roni_path = args.raw_dir / "noaa_roni.ascii.txt"
    soi_path = args.raw_dir / "bom_soi_monthly.txt"
    download(RONI_URL, roni_path, args.refresh)
    download(SOI_URL, soi_path, args.refresh)
    event_context = summarise_events(
        read_roni(roni_path), read_soi(soi_path), EVENT_YEARS
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    event_context.to_csv(args.output, index=False)
    print(event_context.to_string(index=False))
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
