# Current model status

Last reconciled with registry version 7 on 2026-09-11.

## Selected initial forecast

`config/model_registry.yml` selects `operational_rrn_raw_plus_manta_state`, a cause-aware INLA composite trained across LTMP, manta and MMP observations.

Its registered leave-one-event-out performance is:

| Metric | Value |
|---|---:|
| RMSE | 0.1485 |
| Predictive R-squared | 0.2998 |
| Severe-event RMSE | 0.3267 |
| False-extreme rate | 0.0458 |

The model combines the selected thermal/freshwater INLA component with a logistic cyclone hazard and raw RRN event pressure plus prospective Manta COTS state. Its fit artefact is registered at `output/cause_aware_competing_hazards/full_thermal_component.rds`; all generated artefacts remain local and must be recreated or obtained with the data bundle.

## Machine-learning role

The standalone direct BRT improves reef-blocked interpolation but transfers poorly to an unseen event: leave-one-event-out RMSE is 0.1799 versus 0.1485 for the selected INLA composite, and severe RMSE is 0.4650 versus 0.3267. It therefore remains a within-event learner, diagnostic and uncertainty candidate. The residual BRT is not in the selected ensemble.

## BRMS/BRT formal comparison

The formal inflated-beta BRMS and BRT pipeline remains in `src/models/` and has its own `formal-comparators` pipeline profile. It is retained because boundary-aware beta, binomial sensitivity and machine-learning comparisons are scientifically important. It is not the registry-selected production route after the later spatial/cause-aware work.

## Known gaps

- Severe Gannett Cay 2020 mortality remains underpredicted even after stronger COTS pressure.
- Penrith Reef 2017 remains underpredicted after cyclone adjustment.
- Mackay 2024 remains unresolved by available freshwater and cyclone proxies.
- Jasper/Kirrily cyclone exposure and its uncertainty need replacement with improved track-intensity data.
- The reliability-weighted, event-centred aerial/RHIS anomaly has the expected residual direction but negligible skill: sector-excluded RMSE improves 0.00030, block-excluded RMSE worsens 0.00006, and the locked 2025 both-design guard fails. It remains a separate condition layer, off by default.
- In a separate bleaching-compatible DHW-only occurrence diagnostic, logger adjustment improves pooled event-held-out Brier score but worsens reef-, sector- and programme-held-out transfer and changes event slopes inconsistently; it does not alter selection.
- The 2025 result is a retrospective locked-input reconstruction, not an archived issued forecast; 2025 outcomes and aerial/RHIS signals are excluded from its fit.
- ENSO plus current rainfall/WQC helps 2022 but worsens 2020 and suppresses 2024; it is not promoted.
- Operational prediction still needs full uncertainty propagation across DHW correction, INLA components, cyclone activation, COTS pressure and any ensemble.

The reviewed implementation priorities are maintained in `docs/model-next-steps.md`; transient work plans belong in `.ai/plans/`.
