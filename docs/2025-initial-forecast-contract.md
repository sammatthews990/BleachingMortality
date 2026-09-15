# Locked 2025 initial-forecast reconstruction

## Purpose and interpretation

The 2025 product is a retrospective, locked-input reconstruction used as a
forward validation event. It was not archived or issued before the 2025
outcomes were observed, so it must not be described as a prospective forecast.
The outcome is excluded from fitting, and no aerial or RHIS observations enter
the initial prediction.

This product is deliberately separate from the canonical development table and
from the optional within-event update. It does not change the selected model,
`config/model_terms.yml`, or the initial-forecast environmental contract.

## Outcome boundary

`src/data/harmonise_mortality_outcomes.R` assigns 2025 bleaching mortality
using the existing May cutoff: January-April 2026 surveys map to event 2025,
while assessments beginning in May 2025 are already at least two months after
the final aerial observation. The current 2025 assessment contains 63 rows at
47 reefs: 36 Manta, 5 LTMP benthic and 22 MMP rows.

## Locked inputs

`src/data/build_2025_initial_forecast_rows.R` reconstructs predictor rows using:

- RRN DHW history through 2025, with pre-2025 history used for recurrence and
  prior-exposure summaries;
- ERA5 weather from `data/processed/era5_weather_reef_year_2025.csv`;
- PATMOS-x cloud from `data/processed/patmosx_cloud_reef_year_2025.csv`;
- pre-event coral state and the selected raw RRN plus dated Manta COTS state;
- reef-specific 2016-2024 climatology for Kd490, Secchi and eReefs current where
  operational 2025 products are unavailable.

The Kd490, Secchi and current substitutions are explicit imputation, not 2025
observations. The row audit records their source and missingness. Three target
reefs without exact RRN identifiers use the nearest RRN reef for DHW and are
flagged in the audit.

## Fit and validation boundary

`src/evaluation/build_locked_2025_initial_forecast.R` refits the selected
cause-aware architecture using outcomes no later than 2024, predicts the 2025
rows, and writes:

- `output/initial_forecast_2025/locked_predictions.csv`;
- `output/initial_forecast_2025/locked_metrics.csv`;
- `output/initial_forecast_2025/locked_manifest.csv`;
- `output/initial_forecast_2025/locked_thermal_fit.rds`.

The manifest must show a maximum training event of 2024 and both aerial and
RHIS flags as false. The resulting prediction is then used as the unchanged
offset in `src/evaluation/test_spatial_early_bleaching_update.R`. Only that
second, separately labelled product may consume early aerial or RHIS data.

The locked 2025 initial reconstruction has RMSE 0.0609, MAE 0.0489 and
occurrence Brier score 0.2301 over the 63 assessment rows. These metrics test
the architecture under the stated imputation contract; they do not validate
the unavailable optical/current inputs as measured 2025 conditions.
