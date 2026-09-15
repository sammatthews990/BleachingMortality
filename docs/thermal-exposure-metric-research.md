# Thermal exposure metrics beyond Degree Heating Weeks

## Decision summary

Implementation update (15 September 2026): the expanded IMOS/SSTAARS
extraction, NOAA empirical intensity-duration-frequency screen and hourly-logger
cooling validation are complete.
IMOS-equivalent DHW did not deliver a stable cross-programme validation gain.
Standard Hobday duration/intensity gave modest but validation-dependent gains;
the NOAA frequency block improved some within-event spatial tests but worsened
event transfer, while the satellite day/night contrast failed both mortality
transfer and logger construct validation. No thermal field was promoted; see
`analysis/investigations/imos_thermal_metric_assessment.qmd`.

The current model does **not** contain an explicit measure of marine-heatwave
duration, peak thermal intensity, or night-time cooling. It contains
local-first annual maximum Degree Heating Weeks (DHW), DHW spline hinges at 4
and 8, summer SST skewness and excess kurtosis, and indirect cooling modifiers
(cloud and current speed). DHW combines intensity and duration, but it discards
the temporal arrangement of hot days: a continuous hot spell and several
interrupted spells can produce the same DHW.

The completed IMOS screen tested three mechanisms alongside the current DHW
without replacing the selected model:

1. **Continuous exposure:** longest run of night-time SST at or above the
   IMOS/SSTAARS maximum monthly mean (MMM) + 1 degrees C.
2. **Acute intensity:** maximum 3-day mean HotSpot and, as an alternative,
   degree-heating-days above MMM + 2 degrees C.
3. **Lack of thermal relief:** longest run of unrelieved hot nights. A paired
   day-to-following-night cooling measure should be a sensitivity until it is
   validated against loggers.

The IMOS AusTemp product remains the preferred higher-coverage ingestion target. It already
backfills 2012 onward using bias-corrected, night-time, 0.02-degree MultiSensor
L3S SST and the SSTAARS climatology; it supplies daily SST, SST anomalies,
Hobday marine-heatwave categories, Degree Heating Days (DHD), and contributing
day counts.[^1] All five modelled bleaching events (2016, 2017, 2020, 2022 and
2024) are covered. Custom extraction from the three supplied Zarr stores is
still needed for coral-threshold spell metrics and the experimental day/night
cooling comparison.

A full-data BRT is reasonable as **exploratory triage**, but its split-based
importance cannot establish that a metric generalises or is more important
than a correlated DHW term. The defensible sequence is full-data visual and
importance screening, then grouped candidate ablations on the fixed
leave-one-event-out folds, then a small prespecified formal INLA comparison.
The selected model remains unchanged until that sequence is complete.

## What is already represented

The canonical predictor contract in `config/model_terms.yml` and
`docs/methods.md` includes the following thermal information.

| Existing term | What it represents | What it does not identify |
|---|---|---|
| `ann_maxdhw_z` and DHW hinges | Maximum accumulated HotSpot dose and nonlinear response at 4 and 8 DHW | Whether exposure was continuous; the event's peak HotSpot; night-time recovery |
| `dhw10_load4_z`, `dhw_novelty10_z`, prior-event count and recovery interval | Prior heat history, novelty and recurrence | Shape of the current heatwave |
| `sst_summer_skewness_z`, `sst_summer_excess_kurtosis_z` | Shape of the daily summer SST distribution | An interpretable number of hot days, hot-spell duration, or peak intensity |
| `cloudp_90_z`, `mcur_90_z` | Indirect potential for shading, mixing and cooling | Realised reef temperature decline overnight |

The existing SST-shape experiment found plausible within-event information but
no stable gain across programmes and validation schemes. That result is a
reason to use controlled mechanism-level replacements and ablations, not to
add a large collection of correlated SST summaries.

## What the literature supports

### DHW remains a strong baseline, but its algorithm is not uniquely optimal

NOAA Coral Reef Watch DHW sums daily positive HotSpots over the prior 84 days,
but only when SST is at least MMM + 1 degrees C, and divides the sum by seven.
This is a heat-dose metric: it deliberately combines magnitude and duration.[^2]
It does not retain the sequence of qualifying days or distinguish a sharp peak
from a lower, longer exposure with the same integrated dose.

Two large-scale algorithm studies show that the conventional 1-degree cutoff
and 12-week window should be treated as a validated baseline rather than a
biologically immutable definition. DeCarlo found better global bleaching-event
detection after removing the 1-degree accumulation cutoff and shortening the
window to nine weeks, as part of a model that also used regional thresholds and
ENSO state.[^3] Lachs and colleagues tested 234 configurations against 37,871
bleaching observations; performance peaked with thresholds at or below MMM and
4-8 week accumulation windows.[^4] These studies concern bleaching occurrence,
not reef-level post-bleaching mortality, so the present project's held-out
event tests remain decisive.

An IMOS-native `dhw0_56d_max` (maximum 56-day positive anomaly sum divided by
seven, with no +1 cutoff) is therefore a useful **replacement sensitivity**.
It should not enter the same unconstrained candidate set as several nearly
identical accumulation windows.

### Marine-heatwave duration and intensity are standard, separable metrics

Hobday and colleagues define a marine heatwave as at least five days above a
seasonally varying 90th-percentile threshold. Events separated by gaps of two
days or less are joined. Their primary descriptors include duration, maximum
and mean intensity, cumulative intensity, and rates of onset and decline.[^5]
The later category system scales temperature exceedance relative to the local
distance between the climatological median and 90th percentile.[^6]

These metrics directly answer the duration and peak-intensity questions, but
they are not coral bleaching thresholds. A 90th-percentile MHW can occur below
MMM + 1 degrees C, while an exceptionally warm tropical event may be more
biologically meaningful when expressed relative to MMM. Standard Hobday
metrics should therefore be kept as a distinct candidate block:

- maximum MHW duration;
- cumulative MHW intensity;
- maximum MHW category and days at category 2 or higher.

McClanahan and colleagues provide particularly relevant coral evidence. In a
coordinated 2016 Indo-Pacific bleaching study, combinations of warm-spell peak,
warm-spell duration, cool-spell duration and temperature bimodality explained
about half of bleaching variation, while maximum DHW alone explained 9% in
that dataset.[^7] Warm and cool spells were based on site-specific 90th and
10th percentiles, and spells within five days were merged. This is strong
motivation for the proposed screen, but it is a single-event, bleaching-severity
analysis with geographic interactions; it cannot establish transfer to GBR
mortality events.

### Acute extremes can be biologically distinct from accumulated dose

Severe, rapid heatwaves can cause immediate tissue mortality and rapid reef
decay rather than the slower bleaching pathway normally represented by DHW.[^8]
Experimental work also shows that temperature at observation and degree-heating
hours above a high threshold can dominate physiological decline under acute
exposure, with heating rate adding information.[^9] Absolute experimental
thresholds (for example 34 degrees C) are species- and site-specific and should
not be transferred directly to the GBR-wide model.

For daily satellite data, a maximum single observation is too sensitive to
retrieval error. The preferred acute metric is the maximum **3-day mean**
HotSpot. A secondary tail-dose metric, `extreme_dhd_mmm2`, integrates only the
excess above MMM + 2 degrees C and asks whether unusually intense heat adds
information after total DHW is known.

### Intensity-duration-frequency curves add rarity, not just another dose

The proposed intensity-duration-frequency (IDF, or IFD) analogy is useful, but
the linked Scientific Data paper is a rainfall-engineering study rather than a
thermal-ecology study.[^17] It extracts annual maxima at several aggregation
durations, fits single-gauge and regionally pooled extreme-value models, and
estimates intensities associated with specified return periods. The marine
translation has precedent: Gregory and colleagues counted how frequently SST-
anomaly events exceeded combinations of intensity and duration thresholds and
displayed the result as marine-heatwave IDF surfaces.[^18]

Mazdiyasni and colleagues provide the closer heatwave analogue.[^19] They take
the hottest mean temperature over every duration from one to ten consecutive
days in each year, then use fitted marginal distributions and a bivariate
copula to estimate joint intensity-duration exceedance or non-exceedance
probabilities. They deliberately use daily mean rather than daily maximum air
temperature to retain the effect of absent night cooling. For reefs, that
argument supports daily logger means; it does not make the IMOS `dan` foundation
SST an arithmetic daily mean or turn two satellite overpass composites into a
true diel average.

This would be distinct from the present predictors:

| Construct | Duration information | Frequency/rarity information |
|---|---|---|
| DHW | One 84-day thresholded accumulation window | None for the current event; separate project history terms count prior DHW events |
| Hobday block | Realised duration and intensity of threshold-defined events | Event count, but no local return period for the realised joint severity |
| Thermal IDF | Maximum intensity across several fixed durations | Joint historical exceedance frequency or return period at each duration |

For a thermal prototype, define daily anomaly `A(t) = S_n(t) - C(t)` and, for
durations such as 1, 3, 7, 14, 28, 56 and 84 days, calculate each summer's
maximum rolling mean anomaly. A long historical series could then estimate the
return level or empirical exceedance frequency at every duration. A mortality
event could be compressed to at most two or three candidate fields: the mean
log return period across durations, the duration at which rarity is greatest,
and a long-minus-short return-period contrast distinguishing short sharp
exposure from persistent moderate exposure.

It should not yet be fitted from the supplied IMOS archive. The record contains
only about 13 complete summers, whereas the rainfall paper required at least 15
annual maxima for regional analysis and 30 for a single-site fit. Reef-by-reef
GEV estimates would be unstable, and stationary return periods would be
misleading under rapid ocean warming. A regionally pooled, non-stationary model
would also need to account for strong spatial dependence between neighbouring
reef pixels.

A conservative NOAA-native prototype is now implemented. It uses CoralTemp
daily nighttime SST, a fixed pre-validation reference of 1986-2012 summers,
and block maxima at 1, 3, 7, 14, 28, 56 and 84 days. Smoothed empirical return
periods avoid unsupported tail extrapolation; the BRT receives only mean log
return period across durations and the long-minus-short persistence contrast.
The maximum return period is retained for audit but saturates at the 28-year
empirical ceiling for most recent events.

The NOAA IDF block improves reef-blocked RMSE for LTMP by 0.0008 and MMP by
0.0154, but worsens leave-one-event-out RMSE for LTMP by 0.0045 and manta by
0.0070; MMP event-held-out change is negligible. It is therefore potentially
useful for the retained within-event BRT, not for the selected initial forecast.
A copula or GEV sensitivity could distinguish events beyond the empirical
ceiling, but it must quantify tail and non-stationarity uncertainty and must
not be selected merely because extrapolation separates the severe recent
events more strongly.

### Cooling and high-frequency variability are plausible modifiers

Safaie and colleagues compared 20 in-situ environmental variables and seven
satellite metrics across 81 bleaching observations. Mean daily temperature
range over the preceding 30 days was the most influential predictor; greater
daily range was associated with lower bleaching severity.[^10] This supports a
logger-derived daily-temperature-range metric. It does not prove that cooling
within a particular night caused lower mortality: daily range can also encode
long-term acclimatisation, habitat, tides, flow and depth.

Kumagai and colleagues explicitly defined Degree Cooling Weeks as the 84-day
sum of positive `MMM - SST` divided by seven. DCW was not strongly correlated
with DHW and was an important additional bleaching predictor in their Japanese
analysis.[^11] McClanahan's cool-spell duration result points in the same
direction.[^7] A GBR study also showed that sub-bleaching thermal exposure can
prime corals and reduce later cellular mortality and symbiont loss, illustrating
that temperature trajectory and respite can matter independently of total
dose.[^12]

Direct evidence specific to night cooling is promising but limited. A 2016
study of the Bonaparte Archipelago found no widespread mortality despite
regional thermal stress comparable with affected locations and proposed the
greater magnitude of night-time cooling as a cause of the difference.[^13]
Physical and ecological syntheses also identify wind-driven night cooling and
cellular recovery as plausible processes.[^14] These are sufficient grounds
for a targeted test, not for assuming a universal protective coefficient.

## IMOS product assessment

The supplied stores are 0.02 degrees (roughly 2 km), not 2 degrees. Direct
metadata inspection shows that they begin in 2012 and cover all current event
years.[^16]

| Source | Temperature represented | Recommended role |
|---|---|---|
| `imos-srs-sst-ms-1day-night` | Night-time skin SST | Primary raw source for SSTAARS-compatible threshold and hot-night metrics |
| `imos-srs-sst-ms-1day-day` | Daytime skin SST | Paired day/night cooling sensitivity only |
| `imos-srs-sst-ms-1day-dan` | Foundation SST in the absence of diurnal variation | Sensitivity for integrated heat dose; **not** an arithmetic daily mean and not a night-cooling measure |
| IMOS AusTemp | Bias-corrected night SST, anomaly, MHW/MCS category and summer DHD/DHD counts, backfilled from 2012 | First source for standard IMOS/SSTAARS MHW and DHD metrics |
| SSTAARS | Daily 1992-2016 night-SST climatology and daily percentiles | Sole baseline for the IMOS-native metrics; do not reuse NOAA MMM |

IMOS describes L3S day-only and night-only products as skin SST and the
day+night product as foundation SST. The multi-sensor composites are
non-interpolated, so clouds create missing values. IMOS recommends filtering
on quality and subtracting the sensor-specific bias estimate.[^15] AusTemp
implements a documented correction and quality mask: conversion to Celsius,
`+0.17 - sses_bias`, `quality_level > 2`, and removal of the stated poor-quality
flag.[^1]

AusTemp already uses SSTAARS daily 10th, 50th and 90th percentiles for Hobday
categories. Its DHD is the cumulative positive anomaly relative to the daily
SSTAARS climatology, resets on 1 December, and uses a zero-degree threshold.[^1]
That DHD must be labelled as such; it is not the same as NOAA's MMM+1,
84-day DHW.

The 14-day AusTemp mosaic is appropriate for visual continuity but should not
be primary for spell duration. Carrying the most recent valid observation
forward can manufacture long hot spells. Prefer the raw night product or the
GeoPolar AusTemp source, which adds Himawari observations and has about 20%
more coverage.[^1]

There is also a processing-version risk. IMOS documents a transition between
reprocessed `fv02` historical files and `fv01` operational files after 2022.[^15]
Because 2024 is both the severest event and on the new side of that boundary,
the extraction audit must test for a 2022-2023 discontinuity in bias, missingness,
platform mix and variance before attributing greater acute intensity to ecology.

## Proposed metric definitions

Let `S_n(t)` be the quality-controlled, bias-corrected night SST at a reef,
`C(t)` the SSTAARS daily climatological mean, and `MMM` the maximum of the 12
monthly means derived from that same climatology. Define coral HotSpot
`H(t) = S_n(t) - MMM`. Use the fixed project summer window from 1 November of
the preceding year through 30 April of the event year; do not anchor windows to
survey date.

### Priority block A: continuity and acute intensity

| Field | Definition | Hypothesis |
|---|---|---|
| `imos_hotspell_mmm1_max_days` | Longest continuous run with `H(t) >= 1` | Mortality is higher when the same DHW occurs without breaks |
| `imos_hotspot_3d_max` | Maximum rolling 3-day mean of `H(t)` | Short, very hot periods add damage beyond accumulated DHW |
| `imos_extreme_dhd_mmm2` | Sum of `max(H(t) - 2, 0)` in degree C days | The tail above MMM + 2 carries additional mortality information |

The first formal candidate should use the first two fields only. The tail dose
is an alternative to the 3-day peak if their correlation is high.

### Priority block B: standard marine-heatwave structure

Derive these from AusTemp's SSTAARS-based categories or its underlying
percentiles, applying the Hobday minimum five-day duration and two-day gap
joining rule:

| Field | Definition |
|---|---|
| `imos_mhw_max_duration_days` | Longest standard MHW in the fixed summer window |
| `imos_mhw_max_intensity_c` | Maximum SST anomaly above daily climatology within a standard MHW |
| `imos_mhw_cumulative_intensity_c_days` | Sum of daily anomalies within standard MHW days |
| `imos_mhw_category2_days` | Number of days at Hobday category 2 or higher |

This block is a candidate replacement for the current SST skewness and
kurtosis terms. It should not initially replace coral-threshold DHW.

### Priority block C: night-time thermal relief

The most defensible satellite measure is the duration of hot nights because it
uses the climatology-matched night product directly:

| Field | Definition | Status |
|---|---|---|
| `imos_unrelieved_hot_nights` | Longest run of nights with `H(t) >= 1` in the 28 or 56 days around the predictor-defined thermal peak | Primary satellite proxy |
| `imos_night_relief_fraction` | Fraction of hot daytime observations followed by a night below MMM + 1 | Sensitivity |
| `imos_day_to_night_drop_c` | Median bias-corrected day SST minus following-night SST on hot days | Sensitivity requiring logger validation |

Day and night L3S values are multi-sensor, multi-swath composites rather than
paired measurements of the same water parcel. Their retrieval algorithms and
observation times differ. Consequently `day_to_night_drop_c` is a spatial
surface-temperature contrast, not a literal nocturnal cooling rate. Pair by
local GBR date using observation time plus `sst_dtime`, and assess sensitivity
to quality level and platform.

True night cooling should be estimated from high-frequency in-situ data:

- mean daily temperature range over the 30 days around peak stress
  (`logger_dtr30_c`), matching Safaie et al.;
- sunset-to-pre-dawn temperature decline (`logger_nocturnal_drop30_c`);
- fraction of night-time hours below MMM + 1 and longest run of hot nights;
- degree-heating-hours above MMM + 2 for acute local exposure.

The repository currently contains broad daily AIMS temperature coverage for
the five events (`data/processed/aims_temperature_depth_coverage.csv`) and a
much smaller hourly water-quality logger summary through 2023. The main AIMS
temperature pipeline currently requests daily values, so sub-daily records
must be acquired or recovered before DTR or cooling rate can be calculated at
scale. Logger metrics should remain a calibration/subset analysis unless a
fold-safe mapping model proves they can be predicted GBR-wide.

## Extraction and quality contract

1. **Climatology:** use SSTAARS throughout the IMOS analysis. Derive MMM from
   the SSTAARS daily mean and use the published SSTAARS percentiles for Hobday
   metrics. Never mix IMOS SST with NOAA MMM in the primary comparison.
2. **Window:** use 1 November-30 April, keyed to `event_year`; preserve the
   project's May cutoff and established event-year mapping. Peak-centred
   windows may be defined from temperature alone, never from mortality survey
   timing.
3. **Spatial match:** choose one stable water pixel per reef, record its
   distance, and add a local 3-by-3-water-pixel median as a sensitivity. Do not
   select a different nearest valid pixel each cloudy day.
4. **Quality:** reproduce the AusTemp correction/mask for comparable raw
   extraction. Report observed-day fraction, longest missing gap, SSES
   uncertainty and platform/product version.
5. **Missing days:** do not count a missing day as cool or hot. Do not use the
   14-day last-observation mosaic for duration. Pre-specify a conservative
   coverage rule and report left/right-censored spells; test at most short-gap
   interpolation as a sensitivity.
6. **Day/night pairing:** convert observation time to the GBR local date and
   pair a daytime observation with the following night. Keep unpaired values
   missing.
7. **Fold safety:** learn imputation, scaling, any satellite-to-logger
   calibration and correlation selection within each analysis fold.
8. **Provenance:** write reef-event fields to `data/processed/`, generated
   diagnostics to a new `output/imos_thermal_metrics/` directory, and record
   source URLs, access date, code version, correction, quality threshold,
   climatology and temporal coverage in a manifest.

## Model-testing design

### Completed expanded screen (14 September 2026)

The implemented extractor now includes the matched nighttime DHW, coral-
threshold continuity and acute metrics, formal Hobday events, and the paired
daytime-to-following-night sensitivity. Hobday events require at least five
days above the periodically interpolated SSTAARS 90th percentile; qualifying
events separated by at most two days are joined. Categories use the official
SSTAARS median/90th-percentile formula. The paired cooling fields use only
observed quality-controlled values and require at least five hot-day pairs and
50% pair coverage.

Across 304 mortality-validation reef-events, 283 have IMOS DHW, 211 meet the
70% MHW coverage gate, and 260 have an estimable paired cooling contrast. The
MHW duration/intensity block improved LTMP in both full-row validation schemes
and improved manta and MMP leave-one-event-out RMSE, but those latter gains did
not reproduce under reef blocking. On matched complete cases, LTMP's spatial
gain disappeared, manta retained an event-held-out gain but was approximately
neutral spatially, and MMP event-held-out fitting was underpowered.

The cooling block is rejected in its present form. Its apparent MMP
leave-one-event-out gain reversed on observed complete cases, and it worsened
LTMP and manta mortality-scale error. More importantly, at 59 AIMS hourly
logger site-events its Spearman correlation with peak-window afternoon-to-
following-predawn cooling was -0.18 overall and -0.13 after event centring.
The corresponding correlations with 30-day logger DTR were only 0.14 and 0.10.
This supports the original warning that polar-orbiting day/night composites do
not directly measure nocturnal reef cooling.

The next defensible thermal step is an official higher-coverage GeoPolar/
AusTemp extraction and a formal, small MHW-only assessment. Broader subdaily
logger acquisition is required before revisiting nighttime cooling as a
reef-wide predictor.

### Stage 1: data and construct validation

Before using mortality outcomes:

- compare IMOS/SSTAARS MMM+1 84-day dose with raw NOAA DHW and local-first DHW;
- compare IMOS night SST with the existing daily AIMS logger series by depth,
  event and habitat;
- validate the satellite day-night proxy against any available sub-daily
  logger data;
- plot coverage and each metric by event, region and programme;
- audit the 2022-2023 product-version boundary;
- calculate rank correlations and environmental-envelope overlap, with
  particular attention to DHW, summer skewness/kurtosis, cloud, calm wind and
  current speed.

The primary comparison must distinguish whether IMOS provides a better
temperature estimate from whether a new temporal metric provides new biology.
This requires a matched `IMOS_DHW1_84` reconstruction as a source-control:

- NOAA/local-first DHW versus matched IMOS DHW tests the data product;
- matched IMOS DHW versus IMOS duration/peak/cooling tests metric structure.

### Stage 2: exploratory BRT triage

Fit a full-data, no-event-identity BRT only to the cause-restricted
bleaching-compatible rows used by the thermal component. Preserve programme
and reef-event weighting. Prefer the existing two-part occurrence/positive
magnitude architecture over a single Gaussian mortality fit because a metric
may affect whether loss occurs differently from how large positive loss is.

Run these prespecified blocks rather than one unstructured feature dump:

| Candidate | Thermal fields |
|---|---|
| `thermal_current` | Current local-first DHW + current SST shape terms |
| `thermal_imos_source_control` | Matched IMOS DHW1-84 + current nonthermal terms |
| `thermal_duration_peak` | Current DHW + longest MMM+1 spell + 3-day peak |
| `thermal_mhw_structure` | Current DHW + standard MHW duration/intensity; replace SST skewness/kurtosis |
| `thermal_night_relief` | Current DHW + unrelieved hot nights; paired cooling only as a sensitivity |
| `thermal_imos_compact` | Best one field from each stable mechanism block |

Use clustered bootstrap importance, held-out permutation importance and
candidate-block ablation in addition to GBM relative influence. Split-based
importance can divide credit arbitrarily among correlated temperature fields.
Inspect accumulated-local-effect or partial-dependence curves and two-way
surfaces with DHW, but treat them as explanatory.

### Stage 3: blocked predictive screen

The first actual selection test should reuse the fixed folds and report:

- leave-one-event-out RMSE and predictive R-squared (primary);
- event-specific errors, especially 2020/2022 low-mortality controls and 2024;
- occurrence Brier score and calibration;
- positive-magnitude error;
- severe-event RMSE and false-extreme rate;
- reef- and region-blocked results as supporting evidence;
- programme-specific performance and environmental novelty.

Tune the BRT only within each outer analysis set. Candidate definitions must be
fixed before looking at the held-out event. Rank a **block** by out-of-event
gain and stability, not by its full-data variable-importance position.

### Stage 4: formal assessment and promotion boundary

Only one or two stable fields should progress to INLA. Compare:

1. selected thermal core;
2. source-control replacement using matched IMOS DHW;
3. selected core plus the shortlisted metric(s);
4. a compact IMOS replacement if the source-control and temporal metrics both
   outperform.

Fit occurrence and positive-magnitude effects separately, retain the existing
cause restriction, and do not use event identity operationally. Examine
whether a duration or cooling term belongs as a main effect, a DHW interaction,
or only an event-risk sensitivity. Information criteria are secondary to the
matched held-out metrics.

No production ID, selected model, response definition, event window or
registry metric should change until a candidate passes the existing promotion
criteria. When the experiment is completed, add one record to
`config/experiments.yml` with the decision and primary evidence.

## Recommended first implementation slice

The smallest informative implementation is:

1. ingest AusTemp/night SST at validation reef locations for 2016, 2017, 2020,
   2022 and 2024;
2. construct matched IMOS DHW1-84, `imos_hotspell_mmm1_max_days`,
   `imos_hotspot_3d_max`, and `imos_unrelieved_hot_nights`;
3. retain AusTemp `MHW_category` and DHD summaries as secondary fields;
4. run data-source validation and the full-data two-part BRT triage;
5. immediately follow with the existing leave-one-event-out and reef-blocked
   candidate ablations.

The day/night difference and logger DTR work should be a parallel calibration
substudy, not a prerequisite for the first duration/peak screen. This order
answers the strongest questions with daily night SST while avoiding the claim
that polar-orbiting daily composites directly observe nocturnal cooling.

Suggested new paths, if implementation is approved:

- `src/data/extract_imos_thermal_metrics.py`
- `tests/test_imos_thermal_metrics.py`
- `src/evaluation/screen_imos_thermal_metrics.R`
- `analysis/investigations/imos_thermal_metric_assessment.qmd`
- `output/imos_thermal_metrics/` for generated evidence

Remote acquisition should remain outside the default production profile until
the feature contract and coverage checks pass. The investigation can then be
added as a separate pipeline profile without changing the selected model.

## Sources

[^1]: Li, L. (2026). [*IMOS AusTemp Marine Heatwave Data Product, version 1.0*](https://content.aodn.org.au/Documents/IMOS/Data_product/AusTemp-Heatwave_v1.0.pdf). See also the [AODN dataset record](https://researchdata.edu.au/imos-austemp-heat-australian-coast/4014866).
[^2]: NOAA Coral Reef Watch. [*Daily Global 5 km Satellite Coral Bleaching Heat Stress Monitoring Product Suite Methodology*](https://coralreefwatch.noaa.gov/product/5km/methodology.php).
[^3]: DeCarlo, T. M. (2020). [*Treating coral bleaching as weather: a framework to validate and optimize prediction skill*](https://doi.org/10.7717/peerj.9449). *PeerJ*, 8, e9449.
[^4]: Lachs, L., et al. (2021). [*Fine-Tuning Heat Stress Algorithms to Optimise Global Predictions of Mass Coral Bleaching*](https://doi.org/10.3390/rs13142677). *Remote Sensing*, 13, 2677.
[^5]: Hobday, A. J., et al. (2016). [*A hierarchical approach to defining marine heatwaves*](https://doi.org/10.1016/j.pocean.2015.12.014). *Progress in Oceanography*, 141, 227-238.
[^6]: Hobday, A. J., et al. (2018). [*Categorizing and Naming Marine Heatwaves*](https://doi.org/10.5670/oceanog.2018.205). *Oceanography*, 31(2), 162-173.
[^7]: McClanahan, T. R., et al. (2019). [*Temperature patterns and mechanisms influencing coral bleaching during the 2016 El Nino*](https://doi.org/10.1038/s41558-019-0576-8). *Nature Climate Change*, 9, 845-851.
[^8]: Fordyce, A. J., Ainsworth, T. D., Heron, S. F., and Leggat, W. (2019). [*Marine Heatwave Hotspots in Coral Reef Environments*](https://doi.org/10.3389/fmars.2019.00498). *Frontiers in Marine Science*, 6, 498; Leggat, W., et al. (2019). [*Rapid Coral Decay Is Associated with Marine Heatwave Mortality Events on Reefs*](https://www.sciencedirect.com/science/article/pii/S0960982219308048). *Current Biology*, 29, 2723-2730.e4.
[^9]: Evensen, N. R., et al. (2023). [*The roles of heating rate, intensity, and duration on the response of corals and their endosymbiotic algae to thermal stress*](https://doi.org/10.1016/j.jembe.2023.151930). *Journal of Experimental Marine Biology and Ecology*, 567, 151930.
[^10]: Safaie, A., et al. (2018). [*High frequency temperature variability reduces the risk of coral bleaching*](https://doi.org/10.1038/s41467-018-04074-2). *Nature Communications*, 9, 1671.
[^11]: Kumagai, N. H., Yamano, H., and Committee Sango-Map-Project (2018). [*High-resolution modeling of thermal thresholds and environmental influences on coral bleaching for local and regional reef management*](https://doi.org/10.7717/peerj.4382). *PeerJ*, 6, e4382.
[^12]: Ainsworth, T. D., et al. (2016). [*Climate change disables coral bleaching protection on the Great Barrier Reef*](https://doi.org/10.1126/science.aac7125). *Science*, 352, 338-342.
[^13]: Richards, Z. T., et al. (2019). [*A tropical Australian refuge for photosymbiotic benthic fauna*](https://doi.org/10.1007/s00338-019-01809-5). *Coral Reefs*, 38, 669-676.
[^14]: Fordyce et al. (2019), source [^8].
[^15]: IMOS. [*Satellite Remote Sensing Sea Surface Temperature data and product guidance*](https://imos.org.au/srs-sst-data); Govekar, P. D., Griffin, C., and Beggs, H. (2022). [*Multi-Sensor Sea Surface Temperature Products from the Australian Bureau of Meteorology*](https://doi.org/10.3390/rs14153785). *Remote Sensing*, 14, 3785.
[^16]: [IMOS satellite SST STAC item](https://stac.reefdata.io/browser/collections/imos-satellite-remote-sensing/items/imos-srs-sst-ms-1day-dan?.language=en-AU); public Zarr stores: `s3://gbr-dms-data-public/imos-srs-sst-ms-1day-dan/data.zarr`, `s3://gbr-dms-data-public/imos-srs-sst-ms-1day-night/data.zarr`, and `s3://gbr-dms-data-public/imos-srs-sst-ms-1day-day/data.zarr`.
[^17]: Green, A. C., Guerreiro, S. B., and Fowler, H. J. (2026). [*Global Intensity-Duration-Frequency curves based on observed sub-daily rainfall (GSDR-IDF)*](https://doi.org/10.1038/s41597-026-06858-4). *Scientific Data*, 13, 455.
[^18]: Gregory, J. M., et al. (2022). [*An increase in marine heatwaves without significant changes in surface ocean temperature variability*](https://doi.org/10.1038/s41467-022-34934-x). *Nature Communications*, 13, 7396.
[^19]: Mazdiyasni, O., Sadegh, M., Chiang, F., and AghaKouchak, A. (2019). [*Heat wave Intensity Duration Frequency Curve: A Multivariate Approach for Hazard and Attribution Analysis*](https://doi.org/10.1038/s41598-019-50643-w). *Scientific Reports*, 9, 14117.

SSTAARS climatology reference: Wijffels, S. E., et al. (2018).
[*A fine spatial-scale sea surface temperature atlas of the Australian regional
seas (SSTAARS)*](https://doi.org/10.1016/j.jmarsys.2018.07.005). *Journal of
Marine Systems*, 187, 156-196.
