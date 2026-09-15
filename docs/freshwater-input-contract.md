# Operational freshwater-source input contract

Status: proposed ingestion contract; no mortality model refit or predictor
promotion is authorised by this document.

## Decision

Use the Queensland Water Monitoring Information Portal (WMIP) as the primary
operational source of river discharge. WMIP logs and transmits monitoring data
hourly and exposes station metadata, telemetry, archive data and rated discharge
through a documented JSON API. Recent observations are unverified and can later
be revised.

The Bureau of Meteorology Water Data Online service is a useful daily archive
and independent cross-check, but it is not the primary near-real-time feed. The
Bureau states that most Water Data Online data arrive daily and that the portal
is not real time.

This addition creates a **freshwater-source/plume-risk layer**, not a salinity
map. Discharge, ERA5 catchment rainfall, IMOS Kd490 and IMOS chlorophyll remain
separate measured or modelled channels. A value in PSU may be emitted only by a
calibration model that passes held-out low-salinity validation against AIMS
observations and trustworthy eReefs years.

Authoritative sources:

- [Queensland WMIP API guide](https://water-monitoring.information.qld.gov.au/wini/documents/WMIP_API_2025.pdf)
- [Queensland water monitoring description](https://www.qld.gov.au/environment/water/quality/monitoring)
- [Queensland hydrology archive and station network](https://www.data.qld.gov.au/dataset/hydrology-archive-data)
- [Queensland GBR river-basin boundaries](https://www.data.qld.gov.au/dataset/great-barrier-reef-catchment-and-river-basins)
- [Bureau Water Data Online FAQ and SOS2 service](https://www.bom.gov.au/waterdata/wiski-web-public/faq.htm)
- [Bureau Geofabric](https://www.bom.gov.au/water/geofabric/index.shtml)

## Source and acquisition contract

The machine-readable constants are in `config/freshwater_sources.yml`.

- Query WMIP `get_ts_traces` with datasource `ATQ`, source variable `140.00`
  and destination variable `140.00`. This returns stored discharge where
  available and fills gaps by rating telemetry-stage data. Units are cubic
  metres per second.
- Request daily means on the Queensland local calendar. Preserve the station's
  reported timezone and never infer a UTC day boundary.
- Retain the raw API response or a content hash, retrieval time, request
  parameters, licence, station metadata, quality code and quality label.
- Treat API error codes, non-numeric values, and quality labels such as `Data
  not yet available` as missing. In particular, never interpret a numerical
  zero carrying an unavailable quality code as zero discharge.
- Retain provisional/fair/poor values with flags. Do not silently replace an
  operational snapshot when WMIP later validates or revises it; write a new
  snapshot and retain its `retrieved_at_utc`.
- Limit normal polling to no more than hourly at approximately hh:30, as the
  provider recommends. Historical backfills should be batched by basin and
  date range.

The implemented pilot uses the reviewed-but-provisional crosswalk in
`config/freshwater_gauge_crosswalk.csv`: Pioneer as a southern API control;
Daintree and Bloomfield as the relevant sources for Mackay Reef and nearby
northern reefs; and Mulgrave and Russell as additional northern sources. Mackay
Reef is registry feature `16-015` near 16 degrees south. It is not near Mackay
city, so Pioneer discharge must not be labelled its freshwater exposure.

The five-event backfill wrote 3,815 gauge-days. All five selected gauges had
complete requested December--April daily coverage, although many days retained
native fair or poor quality flags. Earlier feasibility queries also showed that
Dumbleton Weir Headwater (`125013A`) can return numerical zeros labelled `Data
not yet available`; the extractor converts such values to missing. Gauge
selection is suitable for a pilot only: high-flow rating limits, regulation,
ungauged area and provisional routing origins remain promotion blockers.

## Harmonised tables

All tables use snake case, UTF-8 CSV, ISO-8601 dates and explicit units. Native
fields are retained in the cached source response.

### `data/processed/freshwater_gauge_day.csv`

One row per station and Queensland local date:

| Field | Contract |
|---|---|
| `station_id`, `station_name` | Stable WMIP identity and native name |
| `basin_id`, `catchment_id` | Provider basin and audited coastal catchment identity |
| `station_lon`, `station_lat`, `coordinate_datum` | Native station position and datum |
| `date_local`, `timezone_offset_hours` | Observation day and source timezone |
| `discharge_mean_m3_s` | Daily mean rated discharge; missing when unavailable |
| `discharge_volume_ml_day` | `discharge_mean_m3_s * 86.4`; derived, not a second observation |
| `quality_code`, `quality_label` | Native WMIP quality information |
| `value_status` | `validated`, `provisional`, `flagged`, or `unavailable` |
| `datasource`, `varfrom`, `varto` | Must retain `ATQ`, `140.00`, `140.00` |
| `retrieved_at_utc`, `snapshot_id` | Availability cutoff and immutable snapshot identity |
| `source_url`, `licence` | Provenance; WMIP responses currently report CC BY 4.0 |

Rows are unique on `(snapshot_id, station_id, date_local)`. Negative discharge
is invalid. Duplicate observations, timezone ambiguity and unit changes fail
the ingestion check.

### `data/processed/freshwater_catchment_event.csv`

One row per catchment, event year, issue time and product mode. It contains only
summaries based on gauge days observable by `issue_time_utc`:

- source station IDs and combination rule;
- window start/end, latest observation time and latency hours;
- expected, observed, unavailable and flagged day counts;
- raw total discharge volume, maximum daily discharge, maximum 7-day mean and
  maximum 30-day mean;
- counts above prespecified historical flow thresholds, with the baseline
  period and support count explicit;
- matching catchment-rainfall summaries and their source;
- regulated-flow, rating-range, ungauged-area and coverage warnings.

Raw summaries are the ingestion product. Log transforms, standardisation,
percentiles and weights must be estimated inside each training fold. Do not
select a favourable gauge or flow window using the held-out event.

### `data/processed/freshwater_reef_event.csv`

One row per `(ReefID, event_year, issue_time_utc, product_mode)`, linked through
an audited river-source-to-reef crosswalk. It retains:

- contributing catchment and gauge IDs;
- static connection method, source-to-reef distance, direction/bearing,
  connection weight and crosswalk version;
- the catchment discharge and rainfall summaries above;
- contemporaneous Kd490 and chlorophyll values, anomaly baselines, observation
  counts, window ends and latency;
- `freshwater_risk_percentile` and uncertainty bounds only after calibration;
- `salinity_observed = false` for every gauge-derived row and no field named
  `salinity_min_psu` unless it comes from a direct salinity observation.

The initial map should display routed discharge/rainfall forcing, Kd490 and
chlorophyll as separate panels. A combined plume-risk percentile is a
sensitivity until it beats rainfall-only, optics-only and spatial-prior
benchmarks under event- and catchment-held-out validation.

## Spatial linkage

Use Queensland GBR river-basin boundaries to define coastal source basins and
Bureau Geofabric to identify the gauge's upstream network and represented
catchment. Start with a curated downstream gauge per major coastal catchment,
plus explicitly weighted tributary gauges only where no defensible downstream
series exists.

Do not use nearest river mouth or nearest gauge alone as the exposure map.
Create a versioned crosswalk with source-to-reef distance and bearing. The
smallest pilot may use a prespecified distance-decay kernel modified by observed
surface-current direction; its output is a connectivity sensitivity. A later
historical eReefs source-release experiment may replace those weights, but its
weights must be frozen before evaluating a held-out event.

## Forecast-time products

`initial_forecast` is the environmental-only product. Its snapshot includes
only gauge, rainfall and satellite observations available by the declared issue
time; it contains no aerial or RHIS observations. An operational map can be
reissued as environmental data accumulate, but every issue has its own cutoff
and immutable snapshot.

`within_event_update` starts from a named initial/environmental snapshot and may
add aerial or RHIS observations available by its later cutoff. If later river
or satellite data are also incorporated, both the environmental window end and
biological-observation cutoff must remain explicit. It is never scored as the
initial forecast.

For retrospective comparisons, produce matched 31 March and 30 April snapshots
where source latency permits. Do not use later validated revisions to represent
what was available at the original operational issue time unless the product is
clearly labelled a hindcast.
The current 2016--2024 WMIP backfill is therefore labelled
`environmental_hindcast`: it uses only environmental observations, but its
retrieval occurred after the historical issue times and it is not scored as an
operational initial forecast or a within-event update.

## Smallest necessary ingestion change

Before any refit:

1. Add one dependency-light Python extractor,
   `src/data/fetch_wmip_discharge.py`, that reads a reviewed station/catchment
   crosswalk, fetches immutable WMIP snapshots, and writes the gauge-day and
   catchment-event tables plus a coverage/provenance audit.
2. Pilot Pioneer as a southern acquisition control and Daintree/Bloomfield plus
   Mulgrave/Russell for Mackay Reef and the northern 2024 region. Audit gauges,
   ratings, regulation, missingness and historical coverage back to 2016 before
   expanding GBR-wide.
3. Add a left join of the audited catchment-event table to
   `src/data/build_aims_freshwater_calibration.R`. Do not alter the canonical
   environmental table, `config/model_terms.yml`, selected model registry, or
   default production pipeline yet.
4. Create component maps and validate the freshwater proxy against AIMS and
   trustworthy eReefs years. Only after that decision should the new table enter
   the environmental contract and matched model refits.

These steps are now implemented as a pilot. Adding routed discharge to rainfall
and Kd490 reduced event-mean average precision from 0.383 to 0.373, while
balanced Brier score improved from 0.314 to 0.300 across 655 AIMS rows. The
mixed result does not support promotion. The output remains smaller than an
unvalidated salinity surface: one acquisition boundary, one static crosswalk,
one audit and one optional calibration-grid join.
