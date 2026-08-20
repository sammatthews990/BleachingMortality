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
