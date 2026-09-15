# Build a separate 2025 assessment table for a leakage-safe initial-forecast
# reconstruction. This never appends 2025 to the canonical development table.

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(sf)
})

root <- normalizePath('.', winslash = '/', mustWork = TRUE)
processed <- file.path(root, 'data', 'processed')
programme_keys <- c('manta', 'ltmp', 'mmp')
outcome_files <- setNames(file.path(
    processed, paste0('mortality_outcomes_', programme_keys, '.rds')
), programme_keys)
environment_file <- file.path(processed, 'cheung_recreated_gbr_full.csv')
rrn_file <- file.path(processed, 'rrn_pressure_reef_year.csv')
cloud_file <- file.path(processed, 'patmosx_cloud_reef_year_2025.csv')
weather_file <- file.path(processed, 'era5_weather_reef_year_2025.csv')
workspace_file <- file.path(processed, '01_exploratory_workspace.RData')
required <- c(outcome_files, environment_file, rrn_file, cloud_file,
              weather_file, workspace_file)
if (any(!file.exists(required))) {
    stop('Missing 2025 initial-forecast inputs: ',
         paste(required[!file.exists(required)], collapse = ', '))
}

normalise_id <- function(x) toupper(trimws(as.character(x)))
mean_or_na <- function(x) if (all(!is.finite(x))) NA_real_ else mean(x, na.rm = TRUE)
median_or_na <- function(x) if (all(!is.finite(x))) NA_real_ else median(x, na.rm = TRUE)

outcomes <- bind_rows(lapply(programme_keys, function(programme_key_value) {
    readRDS(outcome_files[[programme_key_value]]) |>
        filter(event_year == 2025L) |>
        mutate(programme_key = .env$programme_key_value)
})) |>
    mutate(
        ReefID = normalise_id(ReefID),
        reef_name_key = normalise_id(ReefName),
        legacy_survey_max_dhw = MaxDHW.mean
    )
if (nrow(outcomes) == 0L) stop('No eligible 2025 mortality outcomes')

environment <- read_csv(environment_file, show_col_types = FALSE) |>
    mutate(
        ReefID = normalise_id(LABEL_ID),
        reef_name_key = normalise_id(LOC_NAME_S),
        event_year = as.integer(year)
    )
reef_locations_name <- environment |>
    group_by(reef_name_key) |>
    summarise(lon = mean_or_na(lon), lat = mean_or_na(lat), .groups = 'drop')
reef_locations_id <- environment |>
    group_by(ReefID) |>
    summarise(lon_id = mean_or_na(lon), lat_id = mean_or_na(lat),
              .groups = 'drop')

target_locations <- outcomes |>
    distinct(ReefID, ReefName, reef_name_key) |>
    left_join(reef_locations_name, by = 'reef_name_key',
              relationship = 'many-to-one') |>
    left_join(reef_locations_id, by = 'ReefID', relationship = 'many-to-one') |>
    mutate(lon = coalesce(lon, lon_id), lat = coalesce(lat, lat_id)) |>
    select(-lon_id, -lat_id)
if (any(!is.finite(target_locations$lon) | !is.finite(target_locations$lat))) {
    stop('Some 2025 target reefs lack canonical coordinates')
}

rrn <- read_csv(rrn_file, show_col_types = FALSE) |>
    transmute(
        source_ReefID = normalise_id(LABEL_ID),
        event_year = as.integer(event_year),
        ann_maxdhw = as.numeric(sst_maxdhw)
    )
rrn_ids <- rrn |> filter(event_year == 2025L, is.finite(ann_maxdhw)) |>
    distinct(source_ReefID) |>
    left_join(
        reef_locations_id |> rename(source_ReefID = ReefID,
                                    source_lon = lon_id, source_lat = lat_id),
        by = 'source_ReefID', relationship = 'one-to-one'
    ) |>
    filter(is.finite(source_lon), is.finite(source_lat))

haversine_km <- function(lon1, lat1, lon2, lat2) {
    rad <- pi / 180
    dlon <- (lon2 - lon1) * rad
    dlat <- (lat2 - lat1) * rad
    a <- sin(dlat / 2)^2 + cos(lat1 * rad) * cos(lat2 * rad) * sin(dlon / 2)^2
    6371.0088 * 2 * atan2(sqrt(a), sqrt(pmax(1 - a, 0)))
}

source_map <- bind_rows(lapply(seq_len(nrow(target_locations)), function(i) {
    target <- target_locations[i, ]
    exact <- rrn_ids |> filter(source_ReefID == target$ReefID)
    if (nrow(exact)) {
        return(tibble(
            ReefID = target$ReefID, source_ReefID = target$ReefID,
            dhw_match_method = 'exact_reef_id', dhw_source_distance_km = 0
        ))
    }
    distance <- haversine_km(
        target$lon, target$lat, rrn_ids$source_lon, rrn_ids$source_lat
    )
    nearest <- which.min(distance)
    tibble(
        ReefID = target$ReefID,
        source_ReefID = rrn_ids$source_ReefID[[nearest]],
        dhw_match_method = 'nearest_rrn_reef',
        dhw_source_distance_km = distance[[nearest]]
    )
}))

summarise_dhw <- function(source_id) {
    history <- rrn |>
        filter(source_ReefID == source_id, event_year <= 2025L) |>
        arrange(event_year)
    current <- history$ann_maxdhw[history$event_year == 2025L]
    if (length(current) != 1L || !is.finite(current)) {
        stop('Incomplete 2025 RRN DHW series for ', source_id)
    }
    previous <- history |> filter(event_year < 2025L, is.finite(ann_maxdhw))
    recent10 <- previous |> filter(event_year >= 2015L)
    repeated <- previous |> filter(event_year >= 2016L)
    rolling8 <- previous |> filter(event_year >= 2017L)
    last_over <- function(threshold) {
        hit <- previous |> filter(ann_maxdhw >= threshold) |> slice_tail(n = 1)
        if (!nrow(hit)) return(c(value = 0, years_since = 40))
        c(value = hit$ann_maxdhw[[1]], years_since = 2025L - hit$event_year[[1]])
    }
    over4 <- last_over(4)
    over6 <- last_over(6)
    last_n6 <- previous |> filter(ann_maxdhw > 6) |> slice_tail(n = 1)
    years_since_n6 <- if (nrow(last_n6)) 2025L - last_n6$event_year[[1]] else NA_integer_
    tibble(
        ann_maxdhw = current,
        histmDHW6 = over6[['value']], yrsince6 = over6[['years_since']],
        histmDHW4 = over4[['value']], yrsince4 = over4[['years_since']],
        dhw10_load4 = sum(pmax(recent10$ann_maxdhw - 4, 0)),
        dhw_novelty10 = current - max(recent10$ann_maxdhw),
        dhw_events_since2016_n6 = sum(repeated$ann_maxdhw > 6),
        dhw_events_prior8_n6 = sum(rolling8$ann_maxdhw > 6),
        dhw_history_years_since2016 = nrow(repeated),
        dhw_years_since_last_n6 = years_since_n6,
        dhw_years_since_last_n6_capped8 = if_else(
            is.finite(years_since_n6), pmin(years_since_n6, 8), 8
        ),
        dhw_no_prior_n6 = as.numeric(!nrow(last_n6))
    )
}
dhw_features <- bind_rows(lapply(unique(source_map$source_ReefID), function(id) {
    summarise_dhw(id) |> mutate(source_ReefID = id, .before = 1)
}))

climatology <- environment |>
    group_by(ReefID) |>
    summarise(
        mcur_90 = median_or_na(mcur_90),
        dist_to_er_km = median_or_na(dist_to_er_km),
        secc3m = median_or_na(secc3m),
        secc3m_p10 = median_or_na(secc3m_p10),
        .groups = 'drop'
    )
cloud <- read_csv(cloud_file, show_col_types = FALSE) |>
    transmute(
        ReefID = normalise_id(LABEL_ID),
        cloudp_90 = as.numeric(cloudp_90),
        cloud_product_version = as.character(cloud_product_version),
        cloud_aggregation = as.character(cloud_aggregation)
    ) |>
    group_by(ReefID) |>
    summarise(
        cloudp_90 = mean_or_na(cloudp_90),
        cloud_product_version = paste(sort(unique(cloud_product_version)),
                                      collapse = ';'),
        cloud_aggregation = paste(sort(unique(cloud_aggregation)), collapse = ';'),
        .groups = 'drop'
    )

features <- target_locations |>
    left_join(source_map, by = 'ReefID', relationship = 'one-to-one') |>
    left_join(dhw_features, by = 'source_ReefID', relationship = 'many-to-one') |>
    left_join(climatology, by = 'ReefID', relationship = 'one-to-one') |>
    left_join(cloud, by = 'ReefID', relationship = 'one-to-one')
if (any(!complete.cases(features[, c(
    'ann_maxdhw', 'dhw10_load4', 'dhw_novelty10',
    'dhw_events_since2016_n6', 'dhw_years_since_last_n6_capped8',
    'mcur_90', 'secc3m_p10', 'cloudp_90'
)]))) stop('Incomplete required 2025 environmental features')

rows <- outcomes |>
    left_join(features, by = c('ReefID', 'ReefName', 'reef_name_key'),
              relationship = 'many-to-one') |>
    mutate(
        Region = case_when(
            SECTOR %in% c('CG', 'CL', 'PB') ~ 'Northern GBR',
            SECTOR %in% c('CA', 'CU', 'IN', 'TO', 'WH') ~ 'Central GBR',
            SECTOR %in% c('CB', 'PO', 'SW') ~ 'Southern GBR',
            TRUE ~ NA_character_
        ),
        region_block = Region,
        ann_maxsst = NA_real_, winyear_mean = NA_real_, winyear_sd = NA_real_,
        environmental_feature_count = 1L,
        environment_match_method = 'locked_2025_dedicated_snapshot',
        analysis_set = 'prospective_2025_assessment', temporal_fold = 2025L,
        reef_fold = NA_integer_
    )
if (any(is.na(rows$Region))) stop('2025 assessment rows lack a region block')

# Rebuild pre-event Acropora composition from observations no later than each
# row's baseline report year. This mirrors the development-table rule while
# keeping the 2025 assessment separate.
load(workspace_file)
reef_reference <- read_csv(file.path(root, 'data', 'AIMS-Reef_Reference.csv'),
                           show_col_types = FALSE)
reference_exact <- reef_reference |>
    transmute(raw_key = normalise_id(AIMS_REEF_NAME), ReefName = ReefName) |>
    add_count(raw_key, name = 'matches') |> filter(matches == 1L) |>
    select(-matches)
reference_clean <- reef_reference |>
    transmute(
        clean_key = gsub(' REEF(S)?| ISLAND| IS', '', normalise_id(AIMS_REEF_NAME)),
        ReefName = ReefName
    ) |>
    add_count(clean_key, name = 'matches') |> filter(matches == 1L) |>
    select(-matches)
acropora_depth <- df.AIMS.full |>
    filter(data_type == 'photo-transect', domain_category == 'reef',
           purpose == 'COMPOSITION', variable == 'HARD CORAL') |>
    mutate(
        raw_key = normalise_id(domain_name),
        clean_key = gsub(' REEF(S)?| ISLAND| IS', '', normalise_id(domain_name))
    ) |>
    left_join(reference_exact |> rename(ReefName_exact = ReefName), by = 'raw_key',
              relationship = 'many-to-one') |>
    left_join(reference_clean |> rename(ReefName_clean = ReefName), by = 'clean_key',
              relationship = 'many-to-one') |>
    mutate(ReefName = coalesce(ReefName_exact, ReefName_clean)) |>
    filter(!is.na(ReefName)) |>
    group_by(ReefName, report_year, depth) |>
    summarise(
        total = sum(mean, na.rm = TRUE),
        acropora = sum(mean[grepl('Acropora', reefpage_category,
                                  ignore.case = TRUE)], na.rm = TRUE),
        .groups = 'drop'
    ) |>
    mutate(prop_acropora_pre = if_else(total > 0, pmin(pmax(acropora / total, 0), 1),
                                       NA_real_)) |>
    filter(is.finite(prop_acropora_pre))
acropora_reef <- acropora_depth |>
    group_by(ReefName, report_year) |>
    summarise(prop_acropora_reef = mean(prop_acropora_pre), .groups = 'drop')

rows <- rows |>
    left_join(
        acropora_depth |> select(ReefName, report_year, depth, prop_acropora_pre),
        by = c('ReefName', 'baseline_report_year' = 'report_year', 'depth')
    ) |>
    left_join(
        acropora_reef,
        by = c('ReefName', 'baseline_report_year' = 'report_year')
    ) |>
    mutate(
        prop_acropora_pre = coalesce(prop_acropora_pre, prop_acropora_reef),
        acropora_source = case_when(
            is.finite(prop_acropora_pre) ~ 'observed_baseline_year',
            TRUE ~ NA_character_
        ),
        acropora_latest_source_year = if_else(
            is.finite(prop_acropora_pre), baseline_report_year, NA_integer_
        ),
        acropora_source_reefs = if_else(is.finite(prop_acropora_pre), 1L, NA_integer_)
    ) |>
    select(-prop_acropora_reef)

missing_acropora <- which(!is.finite(rows$prop_acropora_pre))
for (index in missing_acropora) {
    available <- acropora_reef |>
        filter(ReefName == rows$ReefName[[index]],
               report_year <= rows$baseline_report_year[[index]]) |>
        arrange(desc(report_year)) |>
        slice_head(n = 1)
    if (nrow(available)) {
        rows$prop_acropora_pre[[index]] <- available$prop_acropora_reef[[1]]
        rows$acropora_source[[index]] <- 'latest_past_same_reef'
        rows$acropora_latest_source_year[[index]] <- available$report_year[[1]]
        rows$acropora_source_reefs[[index]] <- 1L
    }
}
if (any(!is.finite(rows$prop_acropora_pre))) {
    training <- bind_rows(lapply(programme_keys, function(programme_key_value) {
        readRDS(file.path(processed, paste0('validation_rows_', programme_key_value, '.rds')))
    }))
    fallback <- median(training$prop_acropora_pre, na.rm = TRUE)
    rows <- rows |>
        mutate(
            acropora_source = if_else(
                is.finite(prop_acropora_pre), acropora_source,
                'pre2025_training_median'
            ),
            prop_acropora_pre = coalesce(prop_acropora_pre, fallback),
            acropora_source_reefs = coalesce(acropora_source_reefs, 0L)
        )
}
rows <- rows |>
    mutate(
        acropora_cover_pre = prop_acropora_pre * observed_pre_cover,
        acropora_interpolated = acropora_source != 'observed_baseline_year'
    ) |>
    select(-reef_name_key, -source_ReefID)

if (any(rows$acropora_latest_source_year > rows$baseline_report_year,
        na.rm = TRUE)) stop('Future Acropora data entered the 2025 forecast rows')
if (anyDuplicated(paste(rows$programme_key, rows$source_observation_id))) {
    stop('2025 assessment row IDs are not unique')
}

for (programme_key_value in programme_keys) {
    saveRDS(
        rows |> filter(programme_key == .env$programme_key_value),
        file.path(processed, paste0('initial_forecast_rows_2025_', programme_key_value, '.rds'))
    )
}

canonical_validation <- environment |>
    select(ReefID, event_year, canonical_dhw = ann_maxdhw) |>
    inner_join(rrn |> rename(ReefID = source_ReefID, rrn_dhw = ann_maxdhw),
               by = c('ReefID', 'event_year')) |>
    filter(is.finite(canonical_dhw), is.finite(rrn_dhw)) |>
    group_by(event_year) |>
    summarise(
        n = n(), correlation = cor(canonical_dhw, rrn_dhw),
        rmse = sqrt(mean((canonical_dhw - rrn_dhw)^2)),
        bias = mean(rrn_dhw - canonical_dhw), .groups = 'drop'
    )
write_csv(canonical_validation,
          file.path(processed, 'initial_forecast_2025_dhw_validation.csv'))

audit <- tibble(
    contract_version = 1L,
    event_year = 2025L,
    product = 'retrospectively reconstructed locked initial forecast',
    training_outcomes_end = 2024L,
    aerial_or_rhis_in_initial_forecast = FALSE,
    assessment_rows = nrow(rows),
    target_reefs = n_distinct(rows$ReefID),
    exact_rrn_dhw_reefs = n_distinct(rows$ReefID[rows$dhw_match_method == 'exact_reef_id']),
    nearest_rrn_dhw_reefs = n_distinct(rows$ReefID[rows$dhw_match_method == 'nearest_rrn_reef']),
    maximum_dhw_source_distance_km = max(rows$dhw_source_distance_km),
    kd490_2025_available = FALSE,
    kd490_handling = 'reef-specific median of canonical 2016-2024 event values',
    ereefs_2025_available = FALSE,
    current_handling = 'reef-specific median of canonical 2016-2024 event values',
    chla_2025_available = FALSE,
    chla_handling = 'pre-2025 training median during model preprocessing',
    sst_shape_2025_available = FALSE,
    sst_shape_handling = 'pre-2025 training median during model preprocessing',
    cloud_handling = 'NOAA PATMOS-x 2025 January-March snapshot',
    weather_handling = 'ERA5 2024-12-01 through 2025-03-31 snapshot',
    outcome_use = 'scoring only; assessment responses are NA in model fitting',
    environment_md5 = as.character(tools::md5sum(environment_file)),
    rrn_md5 = as.character(tools::md5sum(rrn_file)),
    cloud_md5 = as.character(tools::md5sum(cloud_file)),
    weather_md5 = as.character(tools::md5sum(weather_file))
)
write_csv(audit, file.path(processed, 'initial_forecast_2025_input_audit.csv'))

message('Locked 2025 assessment rows: ', nrow(rows),
        ' across ', n_distinct(rows$ReefID), ' reefs')
