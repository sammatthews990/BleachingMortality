# Canonical methods

This document records the stable scientific contract. Exact predictor terms and selected model IDs are machine-readable in `config/model_terms.yml` and `config/model_registry.yml`.

## Estimand and response

The target is relative coral mortality plausibly associated with a bleaching event at a monitored reef. Programme-specific observations from LTMP, manta tow and MMP are retained rather than treated as interchangeable measurements.

The event window is May to April. Surveys in May to December map to the current calendar event; surveys in January to April map to the preceding event. This intentionally allows at least two to three months after bleaching before mortality assessment and reduces mixing of successive events. Baselines precede assessment and are no more than 24 months old.

Boundary outcomes are retained. The primary model is two-part:

1. a Bernoulli component models whether observed mortality is greater than zero;
2. a beta component models the magnitude among positive observations.

Their product is expected relative mortality. A binomial response using reconstructed counts is retained as a likelihood sensitivity, not assumed to be literally observed colony survival.

The standalone raw-versus-logger-adjusted DHW diagnostic instead uses a
Bernoulli endpoint for whether relative mortality is greater than zero. It
contains only an intercept and DHW slope, excludes explicitly attributed COTS,
cyclone and flood rows, and reports event, sector, reef and programme transfer.
It does not model positive mortality magnitude or participate in selection.

## Disturbance attribution

The thermal/freshwater component is trained on the restricted bleaching-mortality rows. Rows explicitly attributed to COTS or cyclones do not train the thermal response. Separate COTS and cyclone occurrence/magnitude hazards use the longer cause-labelled annual transition record and are combined with the thermal component on the mortality scale:

`1 - (1 - thermal) * (1 - COTS) * (1 - cyclone)`.

This retains non-thermal hazards for operational prediction without forcing background loss at low DHW to flatten the thermal dose-response. Records labelled as multiple disturbances are retained; cause labels are evidence, not assumed perfect truth.

## Predictor families

The operational terms include:

- local-first annual maximum DHW with hinges at 4 and 8 DHW;
- pre-event coral cover and preceding-year Acropora proportion;
- leakage-safe prior heat load, novelty, event count and recovery interval;
- IMOS-derived water clarity and chlorophyll, PATMOS-x cloud and current speed;
- SST distribution summaries and ERA5 coastal rainfall;
- RRN coloured-water history, with prior-only history separated from current exposure;
- cyclone wave/proximity exposure and prospective COTS pressure;
- prespecified interactions between thermal stress and ecological/environmental modifiers.

eReefs salinity is used for years with coverage and for proxy calibration. It is not substituted for missing 2024 measurements. IMOS Kd490/Secchi and chlorophyll are preferred to eReefs optical outputs because the satellite products are the more direct observations. ENSO indices remain an event-state sensitivity and are not selected operational predictors.

## Selected statistical structure

The selected model is a cause-aware INLA composite with programme-specific observation layers, shared ecological effects, a persistent spatial field and event-level structure. It uses the prospective COTS state defined in the registry. Event identity is explanatory only and is excluded from future-event prediction.

The independent direct BRT uses the full attainable operational predictor set and nested tuning. It is retained for within-event spatial interpolation and nonlinear effect diagnostics. BRMS inflated-beta and BRT formal models remain framework comparators.

## Validation and selection

Leave-one-event-out validation is primary because the deployment target is a new bleaching event. Reef-blocked validation evaluates mapping within a known event, and region-blocked validation is a stress test. The 2024 event remains in cross-validation and in the final full-data fit; it is not permanently withheld because it is the broadest and most severe observed event.

The 2025 outcome is a separate forward assessment. Its initial prediction is reconstructed with training outcomes ending in 2024 and without aerial or RHIS inputs. The later aerial/RHIS nowcast is separately labelled and must improve both the blocked historical tests and locked 2025 metrics before promotion.

Primary predictive evidence is overall and event-specific RMSE, predictive R-squared, severe-event RMSE, false-extreme rate and occurrence calibration. Information criteria and coefficient interpretation supplement but do not override held-out prediction.

See `docs/validation.md` for fold construction and prediction modes.
