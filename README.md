# BleachingMortality

## Reconstructing PATMOS-x cloud cover

Install the isolated cloud-pipeline dependencies and reconstruct all five
bleaching-event years before assembling the environmental table:

```powershell
python -m pip install -r scripts/requirements-patmosx.txt
python scripts/fetch_patmosx_cloud.py
python scripts/fetch_dms_environmental_data.py
python scripts/validate_environmental_data.py
```

The cloud step streams only GBR chunks from NOAA's public PATMOS-x v6 archive,
caches them under `data/cache/patmosx_s3_gbr`, and writes the reef-year predictor
and provenance metadata under `data/processed`. It does not use the archived
Cheung cloud values as inputs or fill later years from a historical climatology.

See `reports/research_pipeline_implementation_plan.md` for the ordered pipeline
improvements and cloud validation summary.

The final validation command checks the modelling-table schema, spatially
unique reef-feature-year keys, expected event years, broad physical ranges,
and explicit missing-data limits. It writes a checksum and QA summary to
`data/processed/cheung_recreated_gbr_full.manifest.json`.

## Building harmonised mortality outcomes

After rendering the exploratory analysis, or from its saved workspace, rebuild
the three programme-specific outcome tables with:

```powershell
Rscript scripts/harmonise_mortality_outcomes.R
Rscript scripts/build_validation_splits.R
Rscript scripts/fit_candidate_models.R
quarto render reports/model_validation_report.qmd
```

This writes separate manta tow, LTMP benthic, and MMP inshore CSV/RDS tables
plus `data/processed/mortality_outcomes_manifest.csv`. Surveys follow the
May-to-April event window: May--December maps to the current bleaching year and
January--April to the preceding year. Baselines must be earlier than the
assessment and no more than 24 months old. Boundary outcomes are retained
without nudging.

The validation step joins the reconstructed predictors and saves fixed
leave-one-event-out, reef-blocked, and region-blocked assignments. The 2024
event is assessed as an explicit held-out event fold and then included in the
final production fit. See
`reports/validation_protocol.md` before fitting or comparing models.

The model step compares four prespecified two-part benchmark models under the
blocked assignments. These provide interpretable baselines; the formal learners
use all attainable predictors and are not limited by the benchmark selection.
The final command renders the focused blocked-validation report with candidate,
event-transfer, 2024, and observed-versus-predicted diagnostics.

## Fitting the formal prediction models

After the benchmark outputs exist, fit and compare the two formal frameworks:

```powershell
Rscript scripts/fit_formal_brt.R
Rscript scripts/fit_formal_brms.R
Rscript scripts/summarise_formal_models.R
Rscript scripts/build_formal_ensemble.R
quarto render reports/formal_model_report.qmd
```

Both learners use the fixed event-, reef-, and region-blocked assignments. BRT
tuning occurs inside each outer analysis fold and compares Gaussian with
logit-Gaussian bounded positive magnitude. BRMS uses boundary-aware inflated
beta likelihoods and fold-local preprocessing. Both learners include ten-year
cumulative heat load, signed thermal novelty, and an acute IMOS Kd490 optical
extreme in addition to the established predictors. The scripts save matched
out-of-fold predictions plus full-data production models trained on every
modern event, including 2024. The ensemble step retains a BRMS--BRT framework
blend and an accuracy-focused benchmark--BRMS--BRT stack; region-blocked stress
tests do not determine production weights.

Prediction modes are explicit. Leave-one-event-out is a future-event forecast
and excludes event/region random effects. Reef-blocked validation represents
mapping a known event to unsurveyed reefs and includes known event/region
effects; region-blocked validation retains only the known event effect.

The current validation does not justify forcing one complex model onto every
programme. The rendered report shows the selected prediction route, 2022/2024
component calibration under both prediction modes, variable influence,
production weights, and convergence diagnostics.

## Testing the high-DHW ecological response

The shape-aware screen is separate from the formal benchmark and can be rerun
with:

```powershell
Rscript scripts/fit_shape_aware_brms.R
Rscript scripts/fit_shape_aware_brt.R
Rscript scripts/compare_shape_aware_models.R
quarto render reports/shape_aware_model_comparison.qmd
```

It removes the collinear absolute Acropora-cover predictor, adds a DHW by
proportional-Acropora response, and tests directional ecological priors and a
monotone BRT. The comparison produces a programme-specific routed model and
automatically uses BRMS rather than BRT beyond each programme's observed DHW
support. See `reports/shape_aware_model_comparison.html` for held-out event,
severe-tail, 2024, dose-response, and sampling-diagnostic results.

## Testing joint programme and freshwater-event effects

The joint observation and compact ERA5 weather sensitivity can be rerun with:

```powershell
Rscript scripts/fit_joint_compound_brms.R
Rscript scripts/fit_joint_compound_brt.R
$env:JOINT_BRT_VARIANT='core'
Rscript scripts/fit_joint_compound_brt.R
Rscript scripts/compare_joint_compound_models.R
quarto render reports/joint_compound_stress_assessment.qmd
```

The comparison keeps future-event and 2024 reef-blocked reconstruction targets
separate. Rainfall/calm weather is retained only as a 2024 BRT sensitivity; it
does not replace the current routed production models or the missing direct
2024 salinity record.

## Testing RRN water colour and cyclone pressure

The supplied RRN workbook can be extracted and evaluated with:

```powershell
Rscript scripts/extract_rrn_pressure_data.R
Rscript scripts/evaluate_rrn_pressure_proxies.R
$env:JOINT_BRT_VARIANT='rrn_wq'
Rscript scripts/fit_joint_compound_brt.R
$env:JOINT_BRT_VARIANT='rrn_cyclone'
Rscript scripts/fit_joint_compound_brt.R
$env:JOINT_BRT_VARIANT='rrn'
Rscript scripts/fit_joint_compound_brt.R
$env:JOINT_BRT_VARIANT='compound_rrn'
Rscript scripts/fit_joint_compound_brt.R
Rscript scripts/compare_rrn_pressure_models.R
$env:RRN_BRMS_VARIANT='weather_rrn'
Rscript scripts/fit_rrn_pressure_brms.R
$env:RRN_BRMS_VARIANT='weather_rrn_cots'
Rscript scripts/fit_rrn_pressure_brms.R
Rscript scripts/compare_rrn_brms_models.R
quarto render reports/rrn_pressure_assessment.qmd
```

Fit the shared latent reef-event mortality model and compare it with the
pooled BRMS and BRT benchmarks:

```powershell
Rscript scripts/fit_latent_mortality_brms.R
Rscript scripts/compare_latent_mortality_models.R
```

Build and benchmark annual coral-cover transitions containing both gains and
losses:

```powershell
Rscript scripts/extract_rrn_pressure_data.R
Rscript scripts/build_annual_coral_transitions.R
Rscript scripts/fit_annual_coral_change_benchmarks.R
quarto render reports/latent_and_annual_cover_assessment.qmd
```

Fit the decomposed annual gain/loss candidate while retaining the bleaching
modifiers and restricted relative-mortality response:

```powershell
Rscript scripts/build_annual_coral_transitions.R
Rscript scripts/fit_decomposed_annual_change.R
Rscript scripts/summarise_decomposed_annual_change.R
quarto render reports/decomposed_annual_change_assessment.qmd
```

Austral summer `202324` is aligned to the 2024 mortality event. Water colour
is evaluated as a delayed freshwater/plume proxy; cyclone wave-hours remain a
separate mechanical-disturbance predictor. The report compares both against
complete eReefs salinity years, IMOS Kd490, ERA5 rainfall, and matched blocked
mortality predictions. Neither variable currently replaces the routed
production models.

The RRN extraction also creates a within-reef WQC percentile using only the
preceding ten event years and retains observed and IDW-modelled COTS pressure.
The raw WQC/cyclone and relative-WQC/COTS candidates are evaluated separately;
COTS is not automatically promoted when its blocked prediction is worse.

## Auditing severe 2024 manta misses

The five targeted reefs can be audited against every RRN pressure and metric,
then checked for sensitivity to the preceding Acropora estimate, with:

```powershell
Rscript scripts/audit_rrn_target_reefs.R
Rscript scripts/assess_manta_acropora_counterfactual.R
quarto render reports/rrn_target_reef_diagnostics.qmd
```

The audit preserves the original reef-blocked predictions. It distinguishes
observed from interpolated Acropora estimates, avoids using RRN metrics that
include summer 2024--25 to explain 2024 mortality, and reports the remaining
miss after an Acropora counterfactual rather than treating composition as a
complete explanation.
